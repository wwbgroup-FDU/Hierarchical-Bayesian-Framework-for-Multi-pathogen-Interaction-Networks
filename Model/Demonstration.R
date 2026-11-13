
###############################################################################
# Illustrative case study: IAV -> RSV
# Five evidence configurations:
# M1: VCI
# M2: LOG
# M3: VCI + LOG
# M4: VCI + LOG + OR
#
# By Suyi
###############################################################################

# 0. Packages -----------------------------------------------------------------

library(readxl)
library(dplyr)
library(tidyr)
library(stringr)
library(tibble)
library(ggplot2)
library(rstan)

rstan_options(auto_write = TRUE)
options(mc.cores = 4)
set.seed(2026)


# 1. Read and classify data ---------------------------------------------------

file_path <- "file_path"

df_raw <- read_excel(file_path, sheet = 2)
classify_layer <- function(index, unit) {
  if (index == "OR") {
    return("OR")
  } else if (!is.na(unit) && tolower(unit) == "log10") {
    return("LOG")
  } else {
    return("VCI")
  }
}

df <- df_raw %>%
  mutate(
    layer = mapply(classify_layer, `OR/Index`, 单位)
  )


# 2. Identify virus pairs and directions -------------------------------------

parse_pair <- function(pair_str) {
  parts <- unlist(str_split(pair_str, "-"))
  if (length(parts) == 2) {
    return(parts)
  } else {
    return(c(NA, NA))
  }
}

pair_info <- t(sapply(df$`Virus-pair`, parse_pair))
colnames(pair_info) <- c("pair1", "pair2")

df$from <- ifelse(!is.na(df$From), df$From, pair_info[, 1])
df$to   <- ifelse(!is.na(df$To), df$To, pair_info[, 2])

df$direction <- ifelse(
  is.na(df$From) & is.na(df$To),
  "undirected",
  "directed"
)


# 3. Calculate standard errors -----------------------------------------------

compute_se_from_ci <- function(lower, upper) {
  if (is.na(lower) || is.na(upper) ||
      lower <= 0 || upper <= 0 || lower >= upper) {
    return(NA_real_)
  }
  
  (log(upper) - log(lower)) / (2 * 1.96)
}

se_input <- if ("SE" %in% names(df)) {
  as.numeric(df$SE)
} else {
  rep(NA_real_, nrow(df))
}

df$se <- mapply(
  function(layer, se_raw, lower, upper) {
    if (layer == "OR") {
      if (!is.na(lower) && !is.na(upper)) {
        return(compute_se_from_ci(lower, upper))
      } else {
        return(se_raw)
      }
    } else if (layer == "LOG") {
      return(1.5)
    } else if (layer == "VCI") {
      return(1.0)
    } else {
      return(NA_real_)
    }
  },
  df$layer,
  se_input,
  df$CI_lower,
  df$CI_upper
)

# 4. Calculate sigma and transform values ------------------------------------

df <- df %>%
  mutate(
    sigma = ifelse(
      !is.na(IF) & IF > 0 & !is.na(se),
      se / sqrt(log(IF + 1)),
      NA_real_
    ),
    log_value = case_when(
      layer == "OR"  ~ log(Value),
      layer == "LOG" ~ Value,
      layer == "VCI" ~ log(Value),
      TRUE ~ NA_real_
    )
  )

# Standardize LOG layer only

log_range <- df %>%
  filter(layer == "LOG") %>%
  summarise(
    min_neg = min(log_value[log_value < 0], na.rm = TRUE),
    max_pos = max(log_value[log_value > 0], na.rm = TRUE)
  )
if (!is.finite(log_range$min_neg) ||
    !is.finite(log_range$max_pos)) {
  stop("LOG layer must contain both negative and positive values.")
}

cleaned_data <- df %>%
  mutate(
    log_value_scaled = case_when(
      layer != "LOG" ~ log_value,
      log_value < 0 ~ (log_value / abs(log_range$min_neg)) * 2,
      log_value > 0 ~ (log_value / log_range$max_pos) * 0.3,
      log_value == 0 ~ 0,
      TRUE ~ NA_real_
    )
  ) %>%
  select(
    from, to, direction, layer,
    value = Value,
    log_value,
    log_value_scaled,
    se, sigma, IF,
    source = Num
  )

#Check input data
if (anyNA(cleaned_data$log_value_scaled) ||
    anyNA(cleaned_data$sigma) ||
    any(cleaned_data$sigma <= 0)) {
  stop("Missing or invalid transformed values / sigma. Check source data.")
}


# 5. Generate virus and pair IDs ---------------------------------------------

excluded_node <- "IV"

virus_nodes <- sort(
  unique(c(cleaned_data$from, cleaned_data$to))
)
virus_nodes <- setdiff(virus_nodes, excluded_node)
pair_id <- data.frame(
  virus = virus_nodes,
  id = seq_along(virus_nodes)
)

# Retain pairwise observations among individual virus nodes

network_data <- cleaned_data %>%
  filter(!(from == excluded_node | to == excluded_node))

# Assign node IDs

network_data <- network_data %>%
  left_join(pair_id, by = c("from" = "virus")) %>%
  rename(from_id = id) %>%
  left_join(pair_id, by = c("to" = "virus")) %>%
  rename(to_id = id)

# Assign pair keys

network_data <- network_data %>%
  mutate(
    pair_key = ifelse(
      direction == "undirected",
      paste0(
        pmin(from_id, to_id), "_",
        pmax(from_id, to_id)
      ),
      paste0(from_id, "_", to_id)
    )
  )
pair_lookup <- data.frame(
  pair_key = unique(network_data$pair_key)
) %>%
  mutate(pair_index = row_number())
network_data <- network_data %>%
  left_join(pair_lookup, by = "pair_key") %>%
  filter(!is.na(from_id), !is.na(to_id))
N <- max(pair_lookup$pair_index)



# 6. Identify IAV -> RSV ------------------------------------------------------

if (!all(c("IAV", "RSV") %in% pair_id$virus)) {
  stop("IAV or RSV is missing from the virus node list.")
}

IAV_id <- pair_id$id[pair_id$virus == "IAV"]
RSV_id <- pair_id$id[pair_id$virus == "RSV"]

target_key <- paste0(IAV_id, "_", RSV_id)
target_info <- pair_lookup %>%
  filter(pair_key == target_key)
if (nrow(target_info) != 1) {
  stop(
    "No directed IAV -> RSV pair found. ",
    "Check the original From/To information."
  )
}

target_index <- target_info$pair_index

cat("\nTarget pair: IAV -> RSV\n")
cat("Pair index:", target_index, "\n")

# Inspect evidence availability

target_evidence <- network_data %>%
  filter(pair_index == target_index) %>%
  arrange(layer)
print(
  target_evidence %>%
    select(
      from, to, direction, layer,
      value, log_value_scaled,
      se, sigma, IF, source
    )
)
cat("\nEvidence counts:\n")
print(table(target_evidence$layer))


# 7. Construct layer-specific Stan inputs ------------------------------------

or_layer <- network_data %>%
  filter(layer == "OR")
log_layer <- network_data %>%
  filter(layer == "LOG")
vic_layer <- network_data %>%
  filter(layer == "VCI")
K_or  <- nrow(or_layer)
K_log <- nrow(log_layer)
K_vic <- nrow(vic_layer)
if (K_or == 0 || K_log == 0 || K_vic == 0) {
  stop("The full dataset must contain OR, LOG and VCI observations.")
}

stan_data <- list(
  N = N,
  
  K_or = K_or,
  or_pair = as.integer(or_layer$pair_index),
  theta_or = as.numeric(or_layer$log_value_scaled),
  sigma_or = as.numeric(or_layer$sigma),
  
  K_log = K_log,
  log_pair = as.integer(log_layer$pair_index),
  theta_log = as.numeric(log_layer$log_value_scaled),
  sigma_log = as.numeric(log_layer$sigma),
  
  K_vic = K_vic,
  vic_pair = as.integer(vic_layer$pair_index),
  theta_vic = as.numeric(vic_layer$log_value_scaled),
  sigma_vic = as.numeric(vic_layer$sigma),
  
  target_index = as.integer(target_index)
)


# 8. Stan model ---------------------------------------------------------------

model_code <- "

data {

  int<lower=1> N;

  int<lower=0> K_or;
  int<lower=1,upper=N> or_pair[K_or];
  vector[K_or] theta_or;
  vector<lower=0>[K_or] sigma_or;

  int<lower=0> K_log;
  int<lower=1,upper=N> log_pair[K_log];
  vector[K_log] theta_log;
  vector<lower=0>[K_log] sigma_log;

  int<lower=0> K_vic;
  int<lower=1,upper=N> vic_pair[K_vic];
  vector[K_vic] theta_vic;
  vector<lower=0>[K_vic] sigma_vic;

  int<lower=1,upper=N> target_index;

  int<lower=0,upper=1> use_log;
  int<lower=0,upper=1> use_vic;
  int<lower=0,upper=1> use_or;
}

parameters {

  real w_log;
  real w_vic;

  real intercept;

  real<lower=0> tau;

  vector[N] mu_raw;
}

transformed parameters {

  vector[N] mu_or;

  for (i in 1:N) {

    real mean_log = 0;
    real weight_log = 0;

    real mean_vic = 0;
    real weight_vic = 0;

    if (use_log == 1) {

      for (k in 1:K_log) {

        if (log_pair[k] == i) {

          mean_log += theta_log[k] / square(sigma_log[k]);

          weight_log += 1 / square(sigma_log[k]);
        }
      }
    }

    if (use_vic == 1) {

      for (k in 1:K_vic) {

        if (vic_pair[k] == i) {

          mean_vic += theta_vic[k] / square(sigma_vic[k]);

          weight_vic += 1 / square(sigma_vic[k]);
        }
      }
    }

    real prior_mean = intercept;

    if (weight_log > 0) {
      prior_mean += w_log * (mean_log / weight_log);
    }

    if (weight_vic > 0) {
      prior_mean += w_vic * (mean_vic / weight_vic);
    }

    mu_or[i] = prior_mean + mu_raw[i] * tau;
  }
}

model {

  // Hyperpriors
  w_log ~ normal(1, 1);
  w_vic ~ normal(1, 1);
  intercept ~ normal(0, 2);
  tau ~ normal(0, 1);
  mu_raw ~ normal(0, 1);

  // Population-level OR likelihood
  if (use_or == 1) {

    for (k in 1:K_or) {

      theta_or[k] ~ normal(
        mu_or[or_pair[k]],
        sigma_or[k]
      );
    }
  }
}
"



# 9. Fit evidence configurations ----------------------------------------

fit_model <- function(
    use_log,
    use_vic,
    use_or,
    model_name,
    seed_number
) {
  
  current_data <- stan_data
  
  current_data$use_log <- use_log
  current_data$use_vic <- use_vic
  current_data$use_or <- use_or
  
  cat("\nFitting:", model_name, "\n")
  stan(
    model_code = model_code,
    data = current_data,
    iter = 4000,
    warmup = 2000,
    chains = 4,
    cores = 4,
    seed = seed_number,
    control = list(
      adapt_delta = 0.99,
      max_treedepth = 15
    ),
    pars = c(
      "mu_or",
      "w_log",
      "w_vic",
      "intercept",
      "tau",
      "mu_raw"
    ),
    include = TRUE
  )
}

# M1: VCI only

fit_VCI <- fit_model(
  use_log = 0,
  use_vic = 1,
  use_or = 0,
  model_name = "VCI only",
  seed_number = 202601
)


# M2: LOG only

fit_LOG <- fit_model(
  use_log = 1,
  use_vic = 0,
  use_or = 0,
  model_name = "LOG only",
  seed_number = 202602
)

# M3: VCI + LOG

fit_VCI_LOG <- fit_model(
  use_log = 1,
  use_vic = 1,
  use_or = 0,
  model_name = "VCI + LOG",
  seed_number = 202603
)

# M4: Full model

fit_full <- fit_model(
  use_log = 1,
  use_vic = 1,
  use_or = 1,
  model_name = "Full model",
  seed_number = 202604
)

# 10. Convergence diagnostics -------------------------------------------------

fits <- list(
  VCI = fit_VCI,
  LOG = fit_LOG,
  VCI_LOG = fit_VCI_LOG,
  Full = fit_full
)
diagnostics <- bind_rows(
  lapply(names(fits), function(nm) {
    
    sm <- summary(fits[[nm]])$summary
    
    data.frame(
      Model = nm,
      Max_Rhat = max(sm[, "Rhat"], na.rm = TRUE),
      Min_bulk_ESS = min(sm[, "n_eff"], na.rm = TRUE)
    )
  })
)
print(diagnostics)

# 11. Extract IAV -> RSV posterior samples -----------------------------------

extract_target <- function(fit, model_name) {
  
  draws <- rstan::extract(
    fit,
    pars = c("mu_or", "intercept")
  )
  mu <- draws$mu_or[, target_index]
  
  intercept <- draws$intercept
  
  delta <- mu - intercept
  
  data.frame(
    Model = model_name,
    mu = mu,
    delta = delta
  )
}

draws_all <- bind_rows(
  
  extract_target(fit_VCI, "VCI"),
  
  extract_target(fit_LOG, "LOG"),
  
  extract_target(fit_VCI_LOG, "VCI + LOG"),
  
  extract_target(fit_full, "Full model")
)

# 12. Posterior / prior summary ----------------------------------------------

summary_results <- draws_all %>%
  group_by(Model) %>%
  summarise(
    
    # Absolute VII estimate
    
    mu_median = median(mu),
    
    mu_lower = quantile(mu, 0.025),
    
    mu_upper = quantile(mu, 0.975),
    
    # Relative to population intercept
    
    delta_median = median(delta),
    
    delta_lower = quantile(delta, 0.025),
    
    delta_upper = quantile(delta, 0.975),
    
    # Directional probability
    
    P_positive = mean(delta > 0),
    
    P_negative = mean(delta < 0),
    
    .groups = "drop"
  )
print(summary_results)
write.csv(
  summary_results,
  "IAV_RSV_four_models_summary.csv",
  row.names = FALSE
)


# 13. Save all model outputs --------------------------------------------------

save(
  fit_VCI,
  fit_LOG,
  fit_VCI_LOG,
  fit_VCI_LOG_OR,
  fit_full,
  draws_all,
  summary_results,
  diagnostics,
  file = "IAV_RSV_demonstration_results.RData"
)

cat("\nAll models completed.\n")
cat("Target: IAV -> RSV\n")
