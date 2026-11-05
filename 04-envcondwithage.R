##### -------------------------------------------------------------------- #####
##### -------------------------------------------------------------------- #####
##### -------------------------------------------------------------------- #####
##### -------------------------------------------------------------------- #####
##### -------------------------------------------------------------------- #####

# Annual mean age analysis with environmental variables

# -----------------------------
# 4.1 Summarize hake data into annual mean age and combine with environmental variables
# -----------------------------

# summarise annual mean age
annual_age <- lwa %>%
  group_by(Year) %>% 
  summarise(mean_age = mean(Age_yrs_numeric, na.rm = TRUE),
            median_age = median(Age_yrs_numeric, na.rm = TRUE),
            age_sd = sd(Age_yrs_numeric, na.rm = TRUE),
            age_q25 = quantile(Age_yrs_numeric, probs = 0.25, na.rm = TRUE),
            age_q75 = quantile(Age_yrs_numeric, probs = 0.75, na.rm = TRUE),
            age_iqr = IQR(Age_yrs_numeric, na.rm = TRUE),
            n_fish = sum(!is.na(Age_yrs_numeric)),
            n_sources = n_distinct(Source[!is.na(Source)]),
            .groups = "drop") %>%
  ungroup() %>%
  mutate(Year = as.numeric(as.character(Year)))

# summarise seasonal environmental variables
env_seasonal_summarise <- env_long %>%
  group_by(year, covariate) %>%
  summarise(jantomar = mean(value[month %in% 1:3], na.rm = TRUE),
            aprtojune = mean(value[month %in% 4:6], na.rm = TRUE),
            julytosept = mean(value[month %in% 7:9], na.rm = TRUE),
            octtodec = mean(value[month %in% 10:12], na.rm = TRUE),
            .groups = "drop") %>%
  pivot_longer(cols = c(jantomar, aprtojune, julytosept, octtodec),
               names_to = "season",
               values_to = "value") %>%
  mutate(covariate = paste(covariate, season, sep = "_")) %>%
  mutate(yearlag1 = year + 1) %>%
  select(-season) %>%
  pivot_wider(names_from = covariate, values_from = value)

# join environmental data as an AR(1) (1 year lag) to annual mean age data
annual_age_env <- annual_age %>%
  left_join(env_seasonal_summarise %>% select(-year), 
            by = c("Year" = "yearlag1")) %>%
  arrange(Year) 

# -----------------------------
# 4.2 Plot mean sampled age through time and env variables over time
# -----------------------------

# plot env variables over time
p_env_variables_time <- annual_age_env %>%
  pivot_longer(cols = starts_with(c("pdo", "npgo", "oni", "bi", "npi", "mei")),
               names_to = "covariate",
               values_to = "value") %>%
  mutate(index = toupper(sub("_.*", "", covariate)),
         window = sub(".*_", "", covariate)) %>%
  filter(window == "aprtojune") %>%
  ggplot(aes(x = Year, y = value)) +
    geom_point() +
    geom_smooth(method = "lm", se = TRUE, color = "blue") +
    facet_wrap(~ index, scales = "free") +
    labs(title = "Environmental variables over time (Apr-Jun)",
         x = NULL,
         y = "Standardized value") +
    theme_classic()

p_env_variables_time

# save plot
ggsave(p_env_variables_time, filename = "figures/env_variables_over_time_spring_plot.png", width = 8, height = 5, units = "in", dpi = 300)


# mean age over time
p_annual_mean_age <- ggplot(annual_age_env, aes(x = Year, y = mean_age)) +
  geom_line(linewidth = 0.5) +
  geom_point(aes(size = n_fish), alpha = 0.5) +
  geom_smooth(method = "lm", se = TRUE, color = "blue") +
  scale_size_continuous(trans = "sqrt") +
  labs(title = "Annual mean age of sampled fish",
       x = NULL,
       y = "Mean age (years)",
       size = "Fish sampled") +
  theme_classic()

p_annual_mean_age

# save plot
ggsave(p_annual_mean_age, filename = "figures/annual_mean_age_plot.png", width = 8, height = 5, units = "in", dpi = 300)

# -----------------------------
# 4.3 Plot mean age correlations with environmental variables
# -----------------------------

p_annual_mean_age_env_correlations <- annual_age_env %>%
  pivot_longer(cols = starts_with(c("pdo", "npgo", "oni", "bi", "npi", "mei")),
               names_to = "covariate",
               values_to = "value") %>%
  mutate(index = toupper(sub("_.*", "", covariate)),
         window = sub(".*_", "", covariate)) %>%
  ggplot(aes(x = value, y = mean_age)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "lm", se = TRUE, color = "blue") +
  facet_wrap(index ~ window, scales = "free", ncol = 4) +
  labs(title = "Annual mean age correlations with environmental variables",
       x = "Environmental variable (standardized)",
       y = "Mean age (years)") +
  theme_classic()

p_annual_mean_age_env_correlations

# save figure
ggsave(p_annual_mean_age_env_correlations, filename = "figures/")


# ggpairs plot with one season at a time
annual_age_env %>%
  pivot_longer(cols = starts_with(c("pdo", "npgo", "oni", "bi", "npi", "mei")),
               names_to = "covariate",
               values_to = "value") %>%
  mutate(index = toupper(sub("_.*", "", covariate)),
         window = sub(".*_", "", covariate)) %>%
  filter(window == "aprtojun") %>%
  ggplot() +


# -----------------------------
# 4.4 Analysis of mean age and environmental variables
# -----------------------------
























# -----------------------------
# 4.1 Define the environmental variables
# -----------------------------

# list of environmental variables to use
env_variables <- c(
  "pdo_lagprev1month",
  "npgo_lagprev1month",
  "oni_lagprev1month",
  "npi_lagprev1month",
  "mei_lagprev1month",
  "pdo_lagprev1to3months",
  "npgo_lagprev1to3months",
  "oni_lagprev1to3months",
  "npi_lagprev1to3months",
  "mei_lagprev1to3months",
  "pdo_lagprev1to6months",
  "npgo_lagprev1to6months",
  "oni_lagprev1to6months",
  "npi_lagprev1to6months",
  "mei_lagprev1to6months"
)

# check if there are missing environmental variables in the data to us
missing_env_variables <- setdiff(
  env_variables,
  names(lwa_env)
)

# should be zero
missing_env_variables


# -----------------------------
# 4.2 Prepare the individual-level source data
# -----------------------------

age_source_data <- lwa_env %>%
  filter(
    !is.na(Age_yrs),
    !is.na(Year),
    !is.na(Month)
  ) %>%
  mutate(
    Year = as.integer(Year),
    Month = as.integer(Month),

    date = make_date(
      year = Year,
      month = Month,
      day = 1
    ),

    year_f = factor(Year),
    month_f = factor(
      Month,
      levels = 1:12,
      labels = month.abb
    ),

    Source = factor(Source),
    Sex = factor(Sex)
  ) %>%
  arrange(date)


age_source_data %>%
  mutate(Age_yrs = as.numeric(as.character(Age_yrs))) %>%
  summarise(
    n_fish = n(),
    minimum_age = min(Age_yrs),
    median_age = median(Age_yrs),
    mean_age = mean(Age_yrs),
    maximum_age = max(Age_yrs),
    n_years = n_distinct(Year),
    n_months = n_distinct(date)
  )


# -----------------------------
# 4.3 Check whether environmental values are unique within each month
# -----------------------------

env_uniqueness_check <- age_source_data %>%
  group_by(Year, Month) %>%
  summarise(
    across(
      all_of(env_variables),
      \(x) n_distinct(x[!is.na(x)]),
      .names = "n_{.col}"
    ),
    .groups = "drop"
  )

# each year month combination should have just one unique env variable value. max value should be 1.
env_uniqueness_check %>%
  summarise(
    across(
      starts_with("n_"),
      \(x) max(x, na.rm = TRUE)
    )
  ) %>% as.data.frame()

# -----------------------------
# 4.4 Create monthly mean-age dataset
# -----------------------------

# summarize age by month
monthly_mean_age <- age_source_data %>%
  mutate(Age_yrs = as.numeric(as.character(Age_yrs))) %>%
  group_by(
    Year,
    Month,
    date,
    year_f,
    month_f
  ) %>%
  summarise(
    median_age = median(
      Age_yrs,
      na.rm = TRUE
    ),

    mean_age = mean(
      Age_yrs,
      na.rm = TRUE
    ),

    age_sd = sd(
      Age_yrs,
      na.rm = TRUE
    ),

    age_q25 = quantile(
      Age_yrs,
      probs = 0.25,
      na.rm = TRUE
    ),

    age_q75 = quantile(
      Age_yrs,
      probs = 0.75,
      na.rm = TRUE
    ),

    age_iqr = IQR(
      Age_yrs,
      na.rm = TRUE
    ),

    n_fish = sum(
      !is.na(Age_yrs)
    ),

    n_sources = n_distinct(
      Source[!is.na(Source)]
    ),

    prop_female = mean(
      Sex == "F",
      na.rm = TRUE
    ),

    across(
      all_of(env_variables),
      \(x) {
        values <- unique(x[!is.na(x)])

        if (length(values) == 0L) {
          NA_real_
        } else if (length(values) == 1L) {
          values
        } else {
          stop(
            "More than one environmental value occurred within a year-month."
          )
        }
      }
    ),

    .groups = "drop"
  ) %>%
  arrange(date)

# inspect monthly sample sizes
summary(monthly_mean_age$n_fish)

# determine minimum sample size before modelling
minimum_monthly_sample <- 20

# remove months with less than the minimum_monthly_sample of age estimates
monthly_mean_age <- monthly_mean_age %>%
  filter(n_fish >= minimum_monthly_sample)


# -----------------------------
# 4.5 Standardize the environmental predictors
# -----------------------------

# scale variables and log number of fish
monthly_mean_age <- monthly_mean_age %>%
  mutate(
    across(
      all_of(env_variables),
      \(x) as.numeric(scale(x)),
      .names = "{.col}_z"
    ),
    log_n_fish = log(n_fish),
    time_index = row_number()
  )

# view data
monthly_mean_age %>%
  select(
    date,
    mean_age,
    n_fish,
    all_of(paste0(env_variables, "_z"))
  ) %>%
  glimpse()


# -----------------------------
# 4.6 Plot median and mean sampled age through time
# -----------------------------

# median line plot
ggplot(
  monthly_mean_age,
  aes(
    x = date,
    y = median_age
  )
) +
  geom_line(
    linewidth = 0.5
  ) +
  geom_point(
    aes(size = n_fish),
    alpha = 0.5
  ) +
  scale_size_continuous(
    trans = "sqrt"
  ) +
  labs(
    title = "Monthly median age of sampled fish",
    x = NULL,
    y = "Median age (years)",
    size = "Fish sampled"
  ) +
  theme_classic()

# mean line plot
ggplot(
  monthly_mean_age,
  aes(
    x = date,
    y = mean_age
  )
) +
  geom_line(
    linewidth = 0.5
  ) +
  geom_point(
    aes(size = n_fish),
    alpha = 0.5
  ) +
  scale_size_continuous(
    trans = "sqrt"
  ) +
  labs(
    title = "Monthly mean age of sampled fish",
    x = NULL,
    y = "Mean age (years)",
    size = "Fish sampled"
  ) +
  theme_classic()

# mean box plot
ggplot(
  monthly_mean_age,
  aes(
    x = month_f,
    y = mean_age
  )
) +
  geom_boxplot(
    outlier.alpha = 0.3
  ) +
  labs(
    title = "Seasonal variation in sampled age composition",
    x = "Month",
    y = "Mean age (years)"
  ) +
  theme_classic()

# mean age accross sample size
ggplot(
  monthly_mean_age,
  aes(
    x = log_n_fish,
    y = mean_age
  )
) +
  geom_point(
    alpha = 0.5
  ) +
  geom_smooth(
    method = "lm",
    se = TRUE
  ) +
  labs(
    x = "Log monthly sample size",
    y = "Mean age (years)"
  ) +
  theme_classic()

# mean age by source
ggplot(
  monthly_mean_age,
  aes(
    x = n_sources,
    y = mean_age
  )
) +
  geom_point(
    alpha = 0.5
  ) +
  geom_smooth(
    method = "lm",
    se = TRUE
  ) +
  labs(
    x = "Number of sources represented",
    y = "Mean age (years)"
  ) +
  theme_classic()


# -----------------------------
# 4.7 Fit a monthly null model
# -----------------------------

# monthly mean age
monthly_mean_age_null_model <- brm(
  formula = bf(
    mean_age ~
      month_f +
      s(time_index, k = 10) +
      log_n_fish +
      n_sources +
      (1 | year_f)
  ),

  data = monthly_mean_age,

  family = student(),

  prior = c(
    prior(
      normal(0, 1),
      class = "b"
    ),

    prior(
      normal(8, 3),
      class = "Intercept"
    ),

    prior(
      exponential(1),
      class = "sd"
    ),

    prior(
      exponential(1),
      class = "sigma"
    ),

    prior(
      gamma(2, 0.1),
      class = "nu"
    )
  ),

  chains = 4,
  cores = 4,
  iter = 4000,
  warmup = 2000,

  seed = 123,
  backend = "cmdstanr",

  control = list(
    adapt_delta = 0.97,
    max_treedepth = 13
  ),

  save_pars = save_pars(
    all = TRUE
  ),

  file = "stock-specific/hake/model-outputs/monthly_mean_age_null_model",
  file_refit = "on_change"
)


# --------------------------------------------------
# 4.8 Basic model summaries
# --------------------------------------------------

summary(monthly_mean_age_null_model)

# Useful compact fixed-effect tables
fixef(monthly_mean_age_null_model)

# Random-effect SDs and other variance parameters
VarCorr(monthly_mean_age_null_model)


# --------------------------------------------------
# 4.9 model output diagnostics
# --------------------------------------------------

library(brms)
library(bayesplot)
library(posterior)
library(dplyr)
library(ggplot2)
library(tidyr)

check_nuts <- function(model) {

  np <- nuts_params(model)

  cat("\nDIVERGENCES\n")

  print(
    np %>%
      filter(Parameter == "divergent__") %>%
      count(Value)
  )

  cat("\nTREE DEPTH\n")

  print(
    np %>%
      filter(Parameter == "treedepth__") %>%
      count(Value)
  )

  cat("\nLEAPFROG STEPS\n")

  print(
    np %>%
      filter(Parameter == "n_leapfrog__") %>%
      summarise(
        min = min(Value),
        median = median(Value),
        mean = mean(Value),
        max = max(Value)
      )
  )

  invisible(np)
}

# diagnostics - want zero divergences, and minimal iterations at max tree depth
nuts_mean <- check_nuts(
  monthly_mean_age_null_model
)

# check how many iterations were at max tree depth
nuts_mean %>%
  filter(
    Parameter == "treedepth__",
    Value == 13
  )

# proportion of iterations at max tree depth
treedepth_rate <- function(model, max_depth = 13) {

  nuts_params(model) %>%
    filter(Parameter == "treedepth__") %>%
    summarise(
      n = n(),
      n_at_max = sum(Value >= max_depth),
      proportion_at_max = mean(Value >= max_depth)
    )
}

treedepth_rate(
  monthly_mean_age_null_model
)

# chheck E-BFMI - plots compare distribution
bayesplot::mcmc_nuts_energy(
  nuts_params(monthly_mean_age_null_model)
)

# traceplots
plot(monthly_mean_age_null_model)

draws_mean <- as_draws_df(
  monthly_mean_age_null_model
)

mcmc_trace(
  draws_mean,
  pars = c(
    "b_Intercept",
    "sigma",
    "nu"
  )
)

# -----------------------------
# 4.10 posterior predictive checks
# -----------------------------

pp_check(
  monthly_mean_age_null_model,
  ndraws = 100,
  type = "dens_overlay"
)

pp_check(
  monthly_mean_age_null_model,
  type = "scatter_avg",
  ndraws = 100
)

pp_check(
  monthly_mean_age_null_model,
  type = "stat",
  stat = "max",
  ndraws = 100
)

pp_check(
  monthly_mean_age_null_model,
  type = "stat",
  stat = "sd",
  ndraws = 100
)

# -----------------------------
# 4.11 Smooth temporal trend
# -----------------------------

conditional_smooths(
  monthly_mean_age_null_model
)

plot(
  conditional_smooths(
    monthly_mean_age_null_model
  )
)

smooth_mean <- conditional_smooths(
  monthly_mean_age_null_model
)

smooth_mean

# inspect month effects
conditional_effects(
  monthly_mean_age_null_model,
  effects = "month_f"
)

# sample size effects
plot(
  conditional_effects(
    monthly_mean_age_null_model,
    effects = "log_n_fish"
  )
)

plot(
  conditional_effects(
    monthly_mean_age_null_model,
    effects = "n_sources"
  )
)


# -----------------------------
# 4.12 fitted versus observed value
# -----------------------------
pred_mean <- fitted(
  monthly_mean_age_null_model,
  summary = TRUE
)

mean_check <- monthly_mean_age %>%
  mutate(
    fitted = pred_mean[, "Estimate"],
    fitted_lower = pred_mean[, "Q2.5"],
    fitted_upper = pred_mean[, "Q97.5"]
  )

ggplot(
  mean_check,
  aes(
    x = fitted,
    y = mean_age
  )
) +
  geom_point(alpha = 0.5) +
  geom_abline(
    slope = 1,
    intercept = 0,
    linetype = "dashed"
  ) +
  labs(
    x = "Fitted mean age",
    y = "Observed mean age"
  ) +
  theme_classic()


# -----------------------------
# 4.13 residual structure
# -----------------------------

mean_check <- mean_check %>%
  mutate(
    residual = mean_age - fitted
  )

ggplot(
  mean_check,
  aes(
    x = time_index,
    y = residual
  )
) +
  geom_point(alpha = 0.5) +
  geom_smooth(
    method = "gam",
    formula = y ~ s(x)
  ) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed"
  ) +
  labs(
    x = "Time",
    y = "Observed - fitted mean age"
  ) +
  theme_classic()

# residuals by month
ggplot(
  mean_check,
  aes(
    x = month_f,
    y = residual
  )
) +
  geom_boxplot() +
  geom_hline(
    yintercept = 0,
    linetype = "dashed"
  ) +
  theme_classic()























##### -------------------------------------------------------------------- #####
##### -------------------------------------------------------------------- #####
##### -------------------------------------------------------------------- #####
##### -------------------------------------------------------------------- #####
##### -------------------------------------------------------------------- #####

# Fit an annual median age age model for each env variable

# -----------------------------
# 5.1 Prep data
# -----------------------------


# Will need to do this for seasonal covariates (as opposed to annual below)
# annual_env_variables <- c(
#   "pdo_m_lag1_annual_z",
#   "npgo_m_lag1_annual_z",
#   "oni_m_lag1_annual_z",
#   "npi_m_lag1_annual_z",
#
#   "pdo_m_prev3_annual_z",
#   "npgo_m_prev3_annual_z",
#   "oni_m_prev3_annual_z",
#   "npi_m_prev3_annual_z",
#
#   "pdo_m_prev6_annual_z",
#   "npgo_m_prev6_annual_z",
#   "oni_m_prev6_annual_z",
#   "npi_m_prev6_annual_z"
# )
#
# setdiff(
#   annual_env_variables,
#   names(annual_median_age)
# )


# function to summarize covariates annually - currently doing the three month season Apr-Jun
format_covariate <- function(dat, value.col, name) {
  dat %>%
    as_tibble()  %>%
    filter(month %in% c(4,5,6)) %>% # only spring months
    group_by(year) %>%
    summarize(value = mean(.data[[value.col]])) %>%
    mutate(covariate = name)
}

env_vars <- c("pdo", "npgo", "oni", "bi", "npi", "mei")


# Get environmental covariates
# annual for prior year and current year
env_pdo <- format_covariate(pdo, "anomaly", "pdo")
env_npgo <- format_covariate(npgo, "anomaly", "npgo")
env_oni <- format_covariate(oni, "anomaly", "oni")
env_mei <- format_covariate(mei, "anomaly", "mei")
env_npi <- format_covariate(npi_monthly, "value", "npi")

env_bi <- pacea::bi %>%
  mutate(value = anomaly,
         covariate = "bi")  %>%
  select(year, value, covariate)
# env_npi <- pacea::npi_annual %>% # annual npi covariate
#   mutate(value = anomaly,
#          covariate = "npi")  %>%
#   select(year, value, covariate)

env_long <- bind_rows(env_pdo, env_npgo, env_oni, env_bi, env_npi, env_mei)

env_wide <- env_long %>%
  pivot_wider(names_from = covariate, values_from = value)


# filter only canada or us data
filter_index <- unique(lwa_cond$Source)
filter_index_can <- c("CAN_acoustic", "CAN_shoreside", "CAN_JV","CAN_polish", "Acoustic Canada")
filter_index_us <- c("US_foreign", "US_JV", "US_shore", "US_acoustic", "US_atsea", "Acoustic U.S.")

lwa_cond_can <- lwa_cond %>%
  filter(Source %in% filter_index_can)

lwa_cond_us <- lwa_cond %>%
  filter(Source %in% filter_index_us)


## select which option
selected <- "all"
# selected <- "us"
# selected <- "can"

if(selected == "all"){
  lwa_env <- lwa_cond

} else {
  if(selected == "us"){
    lwa_env <- lwa_cond_us

  } else {
    if(selected == "can"){
      lwa_env <- lwa_cond_can

    } else {
      rm(lwa_env)
    }
  }
}

# join env data with condition data
lwa_env <- lwa_env %>%
  mutate(
    year = as.integer(as.character(year)),
    month = as.integer(as.character(month))
  ) %>%
  left_join(env_wide, by = c("year")) %>%
  mutate(
    year_f = factor(year)
  )

# Scale environmental covariates
lwa_env <- lwa_env %>%
  mutate(across(all_of(env_vars), scale))

# copy to original
lwa_env_original <- lwa_env

annual_env_column_names <- c(
  pdo = "pdo",
  npgo = "npgo",
  oni = "oni",
  bi = "bi",
  npi = "npi",
  mei = "mei"
)

# should be nothing
setdiff(
  unname(annual_env_column_names),
  names(lwa_env)
)

# create new dataframe of annual env variables - will be similar to env_wide dataframe
annual_environment <- lwa_env %>%
  transmute(
    Year = as.integer(as.character(Year)),

    pdo_annual_same =
      .data[[annual_env_column_names[["pdo"]]]],

    npgo_annual_same =
      .data[[annual_env_column_names[["npgo"]]]],

    oni_annual_same =
      .data[[annual_env_column_names[["oni"]]]],

    bi_annual_same =
      .data[[annual_env_column_names[["bi"]]]],

    npi_annual_same =
      .data[[annual_env_column_names[["npi"]]]],

    mei_annual_same =
      .data[[annual_env_column_names[["mei"]]]]
  ) %>%
  group_by(Year) %>%
  summarise(
    across(
      ends_with("_same"),
      \(x) {
        values <- unique(x[!is.na(x)])

        if (length(values) == 0L) {
          return(NA_real_)
        }

        if (length(values) > 1L) {
          stop(
            paste0(
              "More than one annual environmental value was found in ",
              "Year ",
              first(Year),
              "."
            ),
            call. = FALSE
          )
        }

        values
      }
    ),

    .groups = "drop"
  ) %>%
  arrange(Year)

# create lag of 1 year for env variables
annual_environment_prior <- annual_environment %>%
  transmute(
    Year = Year + 1L,

    pdo_annual_prior = pdo_annual_same,
    npgo_annual_prior = npgo_annual_same,
    oni_annual_prior = oni_annual_same,
    bi_annual_prior = bi_annual_same,
    npi_annual_prior = npi_annual_same,
    mei_annual_prior = mei_annual_same
  )

# join lagged and current anv variables
annual_environment_paired <- annual_environment %>%
  left_join(
    annual_environment_prior,
    by = "Year"
  ) %>%
  arrange(Year)

# names of env variables
annual_env_variables <- c(
  "pdo_annual_same",
  "pdo_annual_prior",
  "npgo_annual_same",
  "npgo_annual_prior",
  "oni_annual_same",
  "oni_annual_prior",
  "bi_annual_same",
  "bi_annual_prior",
  "npi_annual_same",
  "npi_annual_prior",
  "mei_annual_same",
  "mei_annual_prior"
)

# bothe same year and prior year env variables
env_variables <- annual_env_variables

# new lwa_env dataframe with env variables (prior and current year)
lwa_env <- lwa_env %>%
  mutate(
    Year = as.integer(as.character(Year))
  ) %>%
  select(
    -any_of(annual_env_variables)
  ) %>%
  left_join(
    annual_environment_paired,
    by = "Year"
  )

# scale all env variables
lwa_env <- lwa_env %>%
  mutate(
    across(
      all_of(env_variables),
      \(x) as.numeric(scale(x)),
      .names = "{.col}_z"
    )
  )

# summarise mean ages by year
# maybe should exclude certain months (e.g. winter months)
annual_mean_age <- lwa_env %>%
  filter(Month %in% c(4:10)) %>%
  group_by(Year) %>%
  mutate(Age_yrs = as.numeric(as.character(Age_yrs))) %>%
  summarise(
    median_age = median(
      Age_yrs,
      na.rm = TRUE
    ),

    mean_age = mean(
      Age_yrs,
      na.rm = TRUE
    ),

    age_sd = sd(
      Age_yrs,
      na.rm = TRUE
    ),

    age_iqr = IQR(
      Age_yrs,
      na.rm = TRUE
    ),

    n_fish = sum(
      !is.na(Age_yrs)
    ),

    n_months = n_distinct(
      Month[!is.na(Month)]
    ),

    n_sources = n_distinct(
      Source[!is.na(Source)]
    ),

    across(
      all_of(env_variables),
      \(x) {
        values <- unique(x[!is.na(x)])

        if (length(values) == 0L) {
          NA_real_
        } else if (length(values) == 1L) {
          values
        } else {
          stop(
            paste0(
              "More than one annual value occurred in Year ",
              first(Year),
              "."
            ),
            call. = FALSE
          )
        }
      }
    ),

    .groups = "drop"
  ) %>%
  arrange(Year) %>%
  mutate(
    across(
      all_of(env_variables),
      \(x) as.numeric(scale(x)),
      .names = "{.col}_z"
    ),

    year_z = as.numeric(
      scale(Year)
    ),

    log_n_fish = log(n_fish)
  )



# -----------------------------
# 5.2 prepare dataset
# -----------------------------

annual_age_env <- annual_mean_age %>%
  mutate(
    Year = as.numeric(Year),

    year_z = as.numeric(
      scale(Year)
    ),

    log_n_fish = log(n_fish)
  ) %>%
  arrange(Year)


# correlation plots
annual_mean_age %>%
  select(mean_age, pdo_annual_same, npgo_annual_same, oni_annual_same, bi_annual_same, npi_annual_same, mei_annual_same) %>%
  ggpairs()

annual_mean_age %>%
  select(mean_age, pdo_annual_prior, npgo_annual_prior, oni_annual_prior, bi_annual_prior, npi_annual_prior, mei_annual_prior) %>%
  ggpairs()

# -----------------------------
# 5.3 function to fit annual environmental model
# -----------------------------

fit_annual_age_env_model <- function(
    data,
    env_variable,
    output_directory = "stock-specific/hake/output-models/annual_age_env",
    seed = 123
) {

  required_variables <- c(
    "mean_age",
    "year_z",
    "log_n_fish",
    "n_months",
    "n_sources",
    env_variable
  )

  missing_variables <- setdiff(
    required_variables,
    names(data)
  )

  if (length(missing_variables) > 0L) {
    stop(
      "Missing variables: ",
      paste(missing_variables, collapse = ", "),
      call. = FALSE
    )
  }

  model_data <- data %>%
    filter(
      if_all(
        all_of(required_variables),
        ~ !is.na(.x)
      )
    )

  model_formula <- bf(
    as.formula(
      paste(
        # "median_age ~",
        "mean_age ~",

        env_variable,
        "+ year_z",
        "+ log_n_fish",
        "+ n_months",
        "+ n_sources"
      )
    )
  )

  dir.create(
    output_directory,
    recursive = TRUE,
    showWarnings = FALSE
  )

  brm(
    formula = model_formula,

    data = model_data,

    family = student(),

    prior = c(
      prior(
        normal(0, 0.5),
        class = "b"
      ),

      prior(
        normal(8, 3),
        class = "Intercept"
      ),

      prior(
        exponential(1),
        class = "sigma"
      ),

      prior(
        gamma(2, 0.1),
        class = "nu"
      )
    ),

    chains = 4,
    cores = 4,
    iter = 4000,
    warmup = 2000,

    backend = "cmdstanr",
    seed = seed,

    control = list(
      adapt_delta = 0.97,
      max_treedepth = 13
    ),

    save_pars = save_pars(
      all = TRUE
    ),

    file = file.path(
      output_directory,
      paste0(
        "annual_mean_age_",

        env_variable
      )
    ),

    file_refit = "on_change",

    refresh = 500
  )
}


# -----------------------------
# 5.4 fit each env variable separately
# -----------------------------

# run model
annual_age_env_models <- annual_env_variables %>%
  set_names() %>%
  map(
    \(env_variable) {
      fit_annual_age_env_model(
        data = annual_age_env,
        env_variable = env_variable
      )
    }
  )


# check model diagnostics
annual_age_env_diagnostics <- imap_dfr(
  annual_age_env_models,
  \(model, variable) {

    # Posterior parameter draws
    draws <- posterior::as_draws_array(model)

    # Parameter-level convergence diagnostics
    diagnostics <- posterior::summarise_draws(
      draws,
      "rhat",
      "ess_bulk",
      "ess_tail"
    )

    # NUTS sampler diagnostics
    nuts <- brms::nuts_params(model)

    tibble(
      covariate = variable,

      maximum_rhat = max(
        diagnostics$rhat,
        na.rm = TRUE
      ),

      minimum_bulk_ess = min(
        diagnostics$ess_bulk,
        na.rm = TRUE
      ),

      minimum_tail_ess = min(
        diagnostics$ess_tail,
        na.rm = TRUE
      ),

      divergences = nuts %>%
        filter(Parameter == "divergent__") %>%
        summarise(
          n = sum(Value)
        ) %>%
        pull(n)
    )
  }
)

annual_age_env_diagnostics

# -----------------------------
# 5.5 extract possterior from each model
# -----------------------------

# function to extract across all models
extract_annual_env_draws <- function(
    model,
    env_variable
) {

  coefficient_name <- paste0(
    "b_",
    env_variable
  )

  model_draws <- as_draws_df(model)

  if (!coefficient_name %in% names(model_draws)) {
    stop(
      "Coefficient not found: ",
      coefficient_name,
      call. = FALSE
    )
  }

  tibble(
    covariate = env_variable,
    beta = model_draws[[coefficient_name]]
  )
}

# run function
annual_age_env_draws <- imap_dfr(
  annual_age_env_models,
  \(model, variable) {
    extract_annual_env_draws(
      model = model,
      env_variable = variable
    )
  }
)

# add labels for plotting
annual_age_env_draws <- annual_age_env_draws %>%
  mutate(
    index = case_when(
      grepl("^pdo", covariate)  ~ "PDO",
      grepl("^npgo", covariate) ~ "NPGO",
      grepl("^oni", covariate)  ~ "ONI",
      grepl("^bi", covariate)  ~ "BI",
      grepl("^npi", covariate)  ~ "NPI",
      grepl("^mei", covariate)  ~ "MEI",
      TRUE                      ~ covariate
    ),

    exposure_window = case_when(
      grepl("_same", covariate) ~
        "Same year",

      grepl("_prior", covariate) ~
        "Previous year",

      TRUE ~ "Other"
    ),

    index = factor(
      index,
      levels = c(
        "PDO",
        "NPGO",
        "ONI",
        "BI",
        "NPI",
        "MEI"
      )
    ),

    exposure_window = factor(
      exposure_window,
      levels = c(
        "Same year",
        "Previous year"
      )
    )
  )

# -----------------------------
# 5.6 summarise posterior coefficints
# -----------------------------

annual_age_env_summary <- annual_age_env_draws %>%
  group_by(
    covariate,
    index,
    exposure_window
  ) %>%
  summarise(
    median_beta = median(beta),

    mean_beta = mean(beta),

    lower_beta = quantile(
      beta,
      0.025
    ),

    upper_beta = quantile(
      beta,
      0.975
    ),

    probability_positive = mean(
      beta > 0
    ),

    probability_negative = mean(
      beta < 0
    ),

    probability_abs_gt_0_25 = mean(
      abs(beta) > 0.25
    ),

    probability_abs_gt_0_50 = mean(
      abs(beta) > 0.50
    ),

    .groups = "drop"
  )

# table of posterior coefficient estimates for each env variables
annual_age_env_summary

# -----------------------------
# 5.7 plot all posterior coefficient distributions
# -----------------------------

# plot and save in object
annual_age_env_posterior_plot <- ggplot(
  annual_age_env_draws,
  aes(
    x = beta,
    y = exposure_window
  )
) +
  stat_halfeye(
    .width = c(
      0.50,
      0.80,
      0.95
    ),
    point_interval = "mean_qi",
    normalize = "groups"
  ) +
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    linewidth = 0.6
  ) +
  facet_wrap(
    vars(index),
    ncol = 2,
    scales = "free_x"
  ) +
  labs(
    title = "Environmental associations with annual mean age",
    subtitle = paste(
      "Each environmental variable was assessed",
      "in a separate model"
    ),
    x = "Change in annual mean age per 1-SD increase",
    y = "Environmental exposure window"
  ) +
  theme_classic(base_size = 12) +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(
      face = "bold"
    )
  )

annual_age_env_posterior_plot

# -----------------------------
# 5.8 forest plot of all posteriors
# -----------------------------

annual_age_env_forest_data <- annual_age_env_summary %>%
  mutate(
    covariate_label = paste(
      index,
      exposure_window,
      sep = ": "
    ),

    covariate_label = factor(
      covariate_label,
      levels = rev(
        unique(covariate_label)
      )
    )
  )


annual_age_env_forest_plot <- ggplot(
  annual_age_env_forest_data,
  aes(
    x = mean_beta,
    y = covariate_label
  )
) +
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    linewidth = 0.6
  ) +
  geom_errorbarh(
    aes(
      xmin = lower_beta,
      xmax = upper_beta
    ),
    height = 0
  ) +
  geom_point(
    size = 2.5
  ) +
  labs(
    title = "Posterior environmental effects on annual mean age",
    subtitle = "Points are means; lines are 95% credible intervals",
    x = "Change in annual mean age (years)",
    y = NULL
  ) +
  theme_classic(base_size = 12)

annual_age_env_forest_plot


# -----------------------------
# 5.9 environmental effects plot - data organizing
# -----------------------------

get_model_data <- function(
    model,
    env_variable,
    source_data = annual_age_env
) {

  required_variables <- c(
    "Year",
    "mean_age",
    "n_fish",
    "year_z",
    "log_n_fish",
    "n_months",
    "n_sources",
    env_variable
  )

  missing_variables <- setdiff(
    required_variables,
    names(source_data)
  )

  if (length(missing_variables) > 0L) {
    stop(
      "Missing variables in source_data: ",
      paste(missing_variables, collapse = ", "),
      call. = FALSE
    )
  }

  source_data %>%
    filter(
      if_all(
        all_of(required_variables),
        ~ !is.na(.x)
      )
    ) %>%
    arrange(Year)
}

get_annual_env_effect <- function(
    model,
    env_variable,
    source_data = annual_age_env,
    prediction_length = 100
) {

  model_data <- get_model_data(
    model = model,
    env_variable = env_variable,
    source_data = source_data
  )

  prediction_data <- tibble(
    env_value = seq(
      min(
        model_data[[env_variable]],
        na.rm = TRUE
      ),
      max(
        model_data[[env_variable]],
        na.rm = TRUE
      ),
      length.out = prediction_length
    ),

    year_z = 0,

    log_n_fish = median(
      model_data$log_n_fish,
      na.rm = TRUE
    ),

    n_months = median(
      model_data$n_months,
      na.rm = TRUE
    ),

    n_sources = median(
      model_data$n_sources,
      na.rm = TRUE
    )
  )

  prediction_data[[env_variable]] <-
    prediction_data$env_value

  fitted_summary <- fitted(
    model,
    newdata = prediction_data,
    re_formula = NA,
    probs = c(0.025, 0.975)
  ) %>%
    as_tibble()

  bind_cols(
    prediction_data,
    fitted_summary
  ) %>%
    mutate(
      covariate = env_variable
    )
}

# -----------------------------
# 5.10 function to plot one env relationship
# -----------------------------

plot_annual_age_env_effect <- function(
    model,
    env_variable,
    source_data = annual_age_env
) {

  model_data <- get_model_data(
    model = model,
    env_variable = env_variable,
    source_data = source_data
  )

  effect_data <- get_annual_env_effect(
    model = model,
    env_variable = env_variable,
    source_data = source_data
  )

  index_label <- case_when(
    grepl("^pdo", env_variable)  ~ "PDO",
    grepl("^npgo", env_variable) ~ "NPGO",
    grepl("^oni", env_variable)  ~ "ONI",
    grepl("^npi", env_variable)  ~ "NPI",
    TRUE                         ~ env_variable
  )

  timing_label <- case_when(
    grepl("_same", env_variable)  ~ "Same year",
    grepl("_prior", env_variable) ~ "Previous year",
    TRUE                          ~ ""
  )

  ggplot(
    effect_data,
    aes(
      x = env_value,
      y = Estimate
    )
  ) +
    geom_ribbon(
      aes(
        ymin = Q2.5,
        ymax = Q97.5
      ),
      alpha = 0.20,
      linewidth = 0
    ) +
    geom_line(
      linewidth = 0.9
    ) +
    geom_point(
      data = model_data,
      aes(
        x = .data[[env_variable]],
        y = mean_age,
        size = n_fish
      ),
      inherit.aes = FALSE,
      alpha = 0.50
    ) +
    geom_rug(
      data = model_data,
      aes(
        x = .data[[env_variable]]
      ),
      inherit.aes = FALSE,
      sides = "b",
      alpha = 0.4
    ) +
    scale_size_continuous(
      trans = "sqrt"
    ) +
    labs(
      title = paste(
        "Annual mean age and",
        index_label
      ),
      subtitle = timing_label,
      x = index_label,
      y = "Predicted annual mean age (years)",
      size = "Fish sampled"
    ) +
    theme_classic(base_size = 12)
}


# -----------------------------
# 5.11 effects plot for each env variable
# -----------------------------

annual_age_env_effect_plots <- imap(
  annual_age_env_models,
  \(model, variable) {
    plot_annual_age_env_effect(
      model = model,
      env_variable = variable,
      source_data = annual_age_env
    )
  }
)

# display all variables in separate plots
walk(
  annual_age_env_effect_plots,
  print
)

# all variables in one plot
annual_age_env_effect_combined <-
  wrap_plots(
    annual_age_env_effect_plots,
    ncol = 2,
    guides = "collect"
  ) &
  theme(
    legend.position = "bottom"
  )

# show single plot
annual_age_env_effect_combined


# -----------------------------
# 5.12 posterior predictions and model evaluation
# -----------------------------

get_annual_fitted_values <- function(
    model,
    env_variable
) {

  model_data <- get_model_data(model = model,
                               env_variable = env_variable)

  fitted_values <- fitted(
    model,
    newdata = model_data,
    re_formula = NA,
    probs = c(
      0.025,
      0.975
    )
  ) %>%
    as_tibble()

  bind_cols(
    model_data %>%
      select(
        Year,
        mean_age,
        n_fish
      ),
    fitted_values
  ) %>%
    mutate(
      covariate = env_variable
    )
}

# table of posterior values with error
annual_age_env_fitted <- imap_dfr(
  annual_age_env_models,
  \(model, variable) {
    get_annual_fitted_values(
      model = model,
      env_variable = variable
    )
  }
)

# calculate residuals
annual_age_env_fitted <- annual_age_env_fitted %>%
  mutate(
    residual = mean_age - Estimate
  )

# plots of observed vs fitted through time
ggplot(
  annual_age_env_fitted,
  aes(x = Year)
) +
  geom_ribbon(
    aes(
      ymin = Q2.5,
      ymax = Q97.5
    ),
    alpha = 0.2
  ) +
  geom_line(
    aes(y = Estimate),
    linewidth = 0.8
  ) +
  geom_point(
    aes(y = mean_age),
    size = 1.5
  ) +
  facet_wrap(
    ~ covariate,
    ncol = 2
  ) +
  labs(
    x = "Year",
    y = "Mean age",
    title = "Observed and fitted annual mean age"
  ) +
  theme_classic()

# plot residuals through time
ggplot(
  annual_age_env_fitted,
  aes(
    x = Year,
    y = residual
  )
) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed"
  ) +
  geom_point() +
  geom_smooth(
    method = "loess",
    se = FALSE
  ) +
  facet_wrap(
    ~ covariate,
    ncol = 2
  ) +
  labs(
    x = "Year",
    y = "Observed - fitted mean age",
    title = "Residual patterns through time"
  ) +
  theme_classic()


# observed vs fitted
ggplot(
  annual_age_env_fitted,
  aes(
    x = Estimate,
    y = mean_age
  )
) +
  geom_point() +
  geom_abline(
    slope = 1,
    intercept = 0,
    linetype = "dashed"
  ) +
  facet_wrap(
    ~ covariate,
    ncol = 2
  ) +
  labs(
    x = "Fitted mean age",
    y = "Observed mean age",
    title = "Observed versus fitted annual mean age"
  ) +
  theme_classic()

# rmse and mean absolute error
annual_age_env_fit_summary <- annual_age_env_fitted %>%
  mutate(
    residual = mean_age - Estimate
  ) %>%
  group_by(covariate) %>%
  summarise(
    n = n(),

    RMSE = sqrt(
      mean(residual^2, na.rm = TRUE)
    ),

    MAE = mean(
      abs(residual),
      na.rm = TRUE
    ),

    bias = mean(
      residual,
      na.rm = TRUE
    ),

    .groups = "drop"
  ) %>%
  arrange(RMSE)

annual_age_env_fit_summary

# -----------------------------
# 5.13 function to plot observed vs fitted ages over time
# -----------------------------


plot_annual_fitted <- function(
    model,
    env_variable
) {

  fitted_data <- get_annual_fitted_values(
    model = model,
    env_variable = env_variable
  )

  plot_title <- gsub(
    "_annual_z$",
    "",
    env_variable
  )

  ggplot(
    fitted_data,
    aes(x = Year)
  ) +
    geom_ribbon(
      aes(
        ymin = Q2.5,
        ymax = Q97.5
      ),
      alpha = 0.20
    ) +
    geom_line(
      aes(
        y = Estimate
      ),
      linewidth = 0.9
    ) +
    geom_point(
      aes(
        y = mean_age,
        size = n_fish
      ),
      alpha = 0.65
    ) +
    scale_size_continuous(
      trans = "sqrt"
    ) +
    labs(
      title = paste(
        "Observed and fitted annual mean age:",
        plot_title
      ),
      subtitle = paste(
        "Line is the posterior expected value;",
        "ribbon is the 95% credible interval"
      ),
      x = "Year",
      y = "Annual mean age (years)",
      size = "Fish sampled"
    ) +
    theme_classic(base_size = 12)
}

annual_age_env_fitted_plots <- imap(
  annual_age_env_models,
  \(model, variable) {
    plot_annual_fitted(
      model = model,
      env_variable = variable
    )
  }
)

# only need one plot to see, but worth it to view other if needed
annual_age_env_fitted_plots[[1]]

# -----------------------------
# 5.14 posterior predictive intervals
# -----------------------------

# posterior predictive intervals
get_annual_predictive_summary <- function(
    model,
    env_variable,
    ndraws = 2000
) {

  model_data <- get_model_data(
    model = model,
    env_variable = env_variable
  )

  # Posterior predictive draws
  predictive_draws <- posterior_predict(
    model,
    newdata = model_data,
    re_formula = NA,
    ndraws = ndraws
  )

  # Summarize posterior predictions
  predictive_summary <- apply(
    predictive_draws,
    MARGIN = 2,
    FUN = function(x) {

      c(
        mean = mean(x),

        lower_50 = unname(quantile(x, 0.25)),
        upper_50 = unname(quantile(x, 0.75)),

        lower_95 = unname(quantile(x, 0.025)),
        upper_95 = unname(quantile(x, 0.975))
      )
    }
  ) %>%
    t() %>%
    as_tibble()

  bind_cols(
    model_data %>%
      select(
        Year,
        mean_age,
        n_fish
      ),
    predictive_summary
  ) %>%
    mutate(
      covariate = env_variable
    )
}

set.seed(123)

annual_age_env_predictive <- imap_dfr(
  annual_age_env_models,
  \(model, variable) {
    get_annual_predictive_summary(
      model = model,
      env_variable = variable,
      ndraws = 2000
    )
  }
)

annual_age_env_predictive

# posterior predictions with observations
plot_annual_predictive <- function(
    predictive_data,
    env_variable
) {

  plot_title <- gsub(
    "_annual_z$",
    "",
    env_variable
  )

  ggplot(
    predictive_data,
    aes(x = Year)
  ) +
    geom_ribbon(
      aes(
        ymin = lower_95,
        ymax = upper_95
      ),
      alpha = 0.15
    ) +
    geom_ribbon(
      aes(
        ymin = lower_50,
        ymax = upper_50
      ),
      alpha = 0.30
    ) +
    geom_line(
      aes(y = mean),
      linewidth = 0.9
    ) +
    geom_point(
      aes(
        y = mean_age,
        size = n_fish
      ),
      alpha = 0.70
    ) +
    scale_size_continuous(
      trans = "sqrt"
    ) +
    labs(
      title = paste(
        "Posterior predictions:",
        plot_title
      ),
      subtitle = paste(
        "Dark ribbon is the 50% predictive interval;",
        "light ribbon is the 95% predictive interval"
      ),
      x = "Year",
      y = "Annual mean age (years)",
      size = "Fish sampled"
    ) +
    theme_classic(base_size = 12)
}

annual_age_env_predictive_plots <- annual_age_env_predictive %>%
  split(.$covariate) %>%
  imap(
    \(data, variable) {
      plot_annual_predictive(
        predictive_data = data,
        env_variable = variable
      )
    }
  )

walk(
  annual_age_env_predictive_plots,
  print
)

# combine plots
annual_age_env_predictive_combined <- wrap_plots(
  annual_age_env_predictive_plots,
  ncol = 3
)

annual_age_env_predictive_combined


# -----------------------------
# 5.15 compare each model to a matched null model
# -----------------------------

# function to create null model

fit_annual_age_null_model <- function(
    data,
    env_variable,
    output_directory = "models/annual_age_env_null",
    seed = 123
) {

  required_variables <- c(
    "mean_age",
    "year_z",
    "log_n_fish",
    "n_months",
    "n_sources",
    env_variable
  )

  model_data <- data %>%
    filter(
      if_all(
        all_of(required_variables),
        ~ !is.na(.x)
      )
    )

  dir.create(
    output_directory,
    recursive = TRUE,
    showWarnings = FALSE
  )

  brm(
    formula = bf(
      median_age ~
        year_z +
        log_n_fish +
        n_months +
        n_sources
    ),

    data = model_data,

    family = student(),

    prior = c(
      prior(
        normal(0, 0.5),
        class = "b"
      ),

      prior(
        normal(8, 3),
        class = "Intercept"
      ),

      prior(
        exponential(1),
        class = "sigma"
      ),

      prior(
        gamma(2, 0.1),
        class = "nu"
      )
    ),

    chains = 4,
    cores = 4,
    iter = 4000,
    warmup = 2000,

    backend = "cmdstanr",
    seed = seed,

    control = list(
      adapt_delta = 0.97,
      max_treedepth = 13
    ),

    save_pars = save_pars(
      all = TRUE
    ),

    file = file.path(
      output_directory,
      paste0(
        "null_",
        env_variable
      )
    ),

    file_refit = "on_change"
  )
}

annual_age_null_models <- annual_env_variables %>%
  set_names() %>%
  map(
    \(env_variable) {
      fit_annual_age_null_model(
        data = annual_age_env,
        env_variable = env_variable
      )
    }
  )


# add LOO to all models
annual_age_env_models <- annual_age_env_models %>%
  map(
    \(model) {
      add_criterion(
        model,
        criterion = "loo",
        moment_match = FALSE
      )
    }
  )

loo_results <- map(
  annual_age_env_models,
  \(model) model$criteria$loo
)

walk(loo_results, print)


# pareto-k diagnostics
loo::pareto_k_table(
  loo_results[["pdo_annual_same"]]
)

imap(
  loo_results,
  \(loo_result, model_name) {
    cat("\nModel:", model_name, "\n")
    print(loo::pareto_k_table(loo_result))
  }
)





# -----------------------------
# 6.3 Model annual condition relate to annual age
# -----------------------------

# Simple association
m_age_condition <- brm(
  mean_age ~ condition_z,
  data = annual_joint,
  family = student(),
  backend = "cmdstanr",
  chains = 4,
  iter = 3000,
  warmup = 1500,
  seed = 22
)

# Adjust for linear time trend
m_age_condition_year <- brm(
  mean_age ~ condition_z + year_z,
  data = annual_joint,
  family = student(),
  backend = "cmdstanr",
  chains = 4,
  iter = 3000,
  warmup = 1500,
  seed = 22
)

fixef(m_age_condition)
fixef(m_age_condition_year)


m_age_condition_smooth <- brm(
  mean_age ~ condition_z + s(year, k = 8),
  data = annual_joint,
  family = student(),
  backend = "cmdstanr",
  chains = 4,
  iter = 3000,
  warmup = 1500,
  seed = 22
)


