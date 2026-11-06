#####
#
# Hake body condition with environmental ocean indices
#
#


##### -------------------------------------------------------------------- #####
##### -------------------------------------------------------------------- #####
##### -------------------------------------------------------------------- #####
##### -------------------------------------------------------------------- #####
##### -------------------------------------------------------------------- #####


# Using posteriors of body condition for bayesian analysis with env covariates

# -----------------------------
# 3.1 Add environmental covariates from pacea
# -----------------------------

# list of pacea indices
pacea_indices

# function for formatting covariates from pacea
format_covariate_month <- function(dat, value.col, name) {
  dat %>%
    as_tibble()  %>%
    group_by(year, month) %>%
    summarize(value = mean(.data[[value.col]])) %>%
    mutate(covariate = name)
}

# Get environmental covariates
env_pdo <- format_covariate_month(pdo, "anomaly", "pdo")
env_npgo <- format_covariate_month(npgo, "anomaly", "npgo")
env_oni <- format_covariate_month(oni, "anomaly", "oni")
env_npi <- format_covariate_month(npi_monthly , "value", "npi")
env_mei <- format_covariate_month(mei , "anomaly", "mei")

# combine variables
env_long <- bind_rows(env_pdo, env_npgo, env_oni, env_npi, env_mei)

##
# seasonal covariates of previous 1, 3, and 6 months - combined into one data set

# make data wide
env_wide <- env_long %>%
  pivot_wider(names_from = covariate, values_from = value)

# estimate prev 1 variables first
# include lag of 1 month for each variable to match with body condition data
lag_months <- 1

# use lag 1 variable names to then find the average across 3 and 6 months
lag1_vars <- c(
  "pdo_lagprev1month",
  "npgo_lagprev1month",
  "oni_lagprev1month",
  "npi_lagprev1month",
  "mei_lagprev1month"
)

# formatting data and making 3 and 6 month rolling averages
env_wide_rolling <- env_long %>%
  mutate(covariate = paste(covariate, "lagprev1month", sep = "_"),
        #  month_orig = month,  # remove this for now as it confuses things
         month = month + lag_months) %>%
  mutate(year = case_when(month > 12 ~ year + 1,
                          .default = year),
         month = case_when(month > 12 ~ (month - 12),
                           .default = month)) %>%
  pivot_wider(names_from = covariate, values_from = value) %>%
  ungroup() %>%

  # Create a proper monthly date
  mutate(
    month_date = make_date(
      year = year,
      month = month,
      day = 1
    )
  ) %>%
  arrange(month_date) %>%

  # Insert missing calendar months, if any
  complete(
    month_date = seq(
      min(month_date, na.rm = TRUE),
      max(month_date, na.rm = TRUE),
      by = "month"
    )
  ) %>%
  mutate(
    year = lubridate::year(month_date),
    month = lubridate::month(month_date)
  ) %>%
  arrange(month_date) %>%

  # Calculate rolling means from the lag-1 variables
  mutate(
    across(
      all_of(lag1_vars),
      list(
        lagprev1to3months = ~ slide_dbl(
          .x,
          mean,
          na.rm = FALSE,
          .before = 2,
          .complete = TRUE
        ),
        lagprev1to6months = ~ slide_dbl(
          .x,
          mean,
          na.rm = FALSE,
          .before = 5,
          .complete = TRUE
        )
      ),
      .names = "{sub('_lagprev1month$', '', .col)}_{.fn}"
    )
  )

names(env_wide_rolling)

# -----------------------------
# 3.2 Correlations with environmental covariates
# -----------------------------

# join with body condition data
dat_env <- cond_monthly %>%
  left_join(env_wide_rolling, by = c("year", "month"))

# Names of dat_env environmental data columns:
env_vars <- c( 
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

# correlations of environmental variables with mean body condition
cor_results <- map_dfr(env_vars, function(v) {
  tmp <- dat_env %>%
    select(mean_log_Kn, all_of(v)) %>%
    drop_na()

  ct <- cor.test(tmp$mean_log_Kn, tmp[[v]], method = "pearson")

  tibble(
    covariate = v,
    n = nrow(tmp),
    correlation = unname(ct$estimate),
    p_value = ct$p.value,
    conf_low = ct$conf.int[1],
    conf_high = ct$conf.int[2]
  )
}) %>%
  mutate(p_adj = p.adjust(p_value, method = "BH"))

# print correlation results
cor_results

# correlation plot of mean condition for previous 3 months only
ggpairs(dat_env %>% select(mean_Kn, contains("lagprev1to3months"))) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

# save figure
ggsave(
  "figures/condition_env_correlations.png",
  width = 8,
  height = 6
)


# -----------------------------
# 3.3 Mixed effects models for environmental effects
# -----------------------------

# filter only canada or us data
filter_index <- unique(lwa_cond$Source)
filter_index_can <- c("CAN_acoustic", "CAN_shoreside", "CAN_JV","CAN_polish", "Acoustic Canada")
filter_index_us <- c("US_foreign", "US_JV", "US_shore", "US_acoustic", "US_atsea", "Acoustic U.S.")

# selected data
selected <- "all"
# selected <- "us"
# selected <- "can"

if(selected == "all"){
  lwa_env <- lwa_cond %>%
    filter(Source %in% filter_index)

} else {
  if(selected == "us"){
    lwa_env <- lwa_cond %>%
      filter(Source %in% filter_index_us)

  } else {
    if(selected == "can"){
      lwa_env <- lwa_cond %>%
        filter(Source %in% filter_index_can)

    } else {
      rm(lwa_env)
    }
  }
}

# Fish-level model: body condition as response.
# Environmental covariates joined to each fish by Year-Month.
lwa_env <- lwa_env %>%
  mutate(
    year = as.integer(as.character(year)),
    month = as.integer(as.character(month))
  ) %>%
  left_join(env_wide_rolling, by = c("year", "month")) %>%
  mutate(
    year_f = factor(year),
    month_f = factor(month)
  )

# Scale environmental covariates
lwa_env <- lwa_env %>%
  mutate(across(all_of(env_vars), scale))

# Fit one model per environmental covariate - potentially need to change formula
env_models <- map(env_vars, function(v) {
  form <- as.formula(
    paste0(
      "log_Kn ~ ", v,
      " + Sex + ",
      "(1 | Source) + (1 | year_f) + (1 | month_f)"
    )
  )

  glmmTMB(
    form,
    data = lwa_env,
    family = gaussian(),
    control = glmmTMBControl(
      optCtrl = list(iter.max = 1e4, eval.max = 1e4)
    )
  )
})

# naming each object in list with variable name
names(env_models) <- env_vars

# extract model results output - coefficient estimate with error and p value
env_model_results <- map_dfr(names(env_models), function(v) {
  tidy(env_models[[v]], effects = "fixed") %>%
    filter(term == v) %>%
    mutate(covariate = v)
}) %>%
  mutate(p_adj = p.adjust(p.value, method = "BH"))

env_model_results


# -----------------------------
# 3.4 Plot environmental relationships
# -----------------------------

# filter out Source
cond_monthly <- lwa_cond %>%
  filter(Source %in% filter_index) %>%
  group_by(year, month) %>%
  summarise(
    mean_Kn = mean(Kn, na.rm = TRUE),
    median_Kn = median(Kn, na.rm = TRUE),
    mean_log_Kn = mean(log_Kn, na.rm = TRUE),
    sd_log_Kn = sd(log_Kn, na.rm = TRUE),
    se_log_Kn = sd_log_Kn / sqrt(n()),
    n_fish = n(),
    .groups = "drop"
  )

dat_env <- cond_monthly %>%
  left_join(env_wide_rolling, by = c("year", "month"))

# plotting function for environmental covariates vs mean body condition
plot_env_relationship <- function(v) {
  ggplot(dat_env, aes(x = .data[[v]], y = mean_Kn)) +
    geom_hline(yintercept = 1, linetype = "dashed") +
    geom_point(aes(size = n_fish), alpha = 0.6) +
    scale_size_continuous(range = c(1, 3)) +
    geom_smooth(method = "lm", se = TRUE) +
    labs(
      x = v,
      y = "Monthly mean Kn",
      size = "n fish"
    #   title = paste("Body condition vs", v, "-", selected, "data")
    ) +
    theme_bw()
}

# create list of plots for each environmental covariate
env_plots <- map(env_vars, plot_env_relationship)

# combine plots into one figure
wrap_plots(env_plots)

# save figure
ggsave(
  "figures/glmm_condition_env_relationships.png",
  width = 12,
  height = 8
)

rm(env_plots)


# -----------------------------
# 3.5 Model-based predicted plots
# -----------------------------

pred_lw <- ggpredict(
  m_lw_selected,
  terms = c("log_length [all]", "Sex")
) %>%
  as_tibble() %>%
  mutate(
    Length_cm = exp(x),
    Weight_kg = exp(predicted),
    conf.low = exp(conf.low),
    conf.high = exp(conf.high)
  )

# model fitted lines for male and female length-weight relationship
ggplot() +
  geom_point(
    data = lwa,
    aes(Length_cm, Weight_kg),
    alpha = 0.03,
    size = 0.4
  ) +
  geom_ribbon(
    data = pred_lw,
    aes(Length_cm, ymin = conf.low, ymax = conf.high, fill = group),
    alpha = 0.25
  ) +
  geom_line(
    data = pred_lw,
    aes(Length_cm, Weight_kg, colour = group),
    linewidth = 1
  ) +
  scale_x_log10() +
  scale_y_log10() +
  labs(
    x = "Length (cm)",
    y = "Weight (kg)",
    colour = "Sex",
    fill = "Sex",
    title = "Mixed-effects length-weight model predictions"
  ) +
  theme_bw()

# save figure
ggsave(
  "figures/glmm_length_weight_model_predictions.png",
  width = 8,
  height = 6
)






##### -------------------------------------------------------------------- #####
##### -------------------------------------------------------------------- #####
##### -------------------------------------------------------------------- #####
##### -------------------------------------------------------------------- #####
##### -------------------------------------------------------------------- #####

# Bayesian analysis of hake body condition with environmental covariates
# Using posteriors of body condition for bayesian analysis with env covariates     

# -----------------------------
# 3.11 Join monthly environmental covariates
# -----------------------------

# select which subset of data to use for analysis: all, us, or can
selected <- "all"
# selected <- "us"
# selected <- "can"

# condition monthly draws from previous section 02-lengthweightanalysis.R
if(selected == "all"){
  cond_monthly_draws_used <- cond_monthly_draws

} else {
  if(selected == "us"){
    cond_monthly_draws_used <- cond_monthly_draws_us

  } else {
    if(selected == "can"){
      cond_monthly_draws_used <- cond_monthly_draws_can

    } else {
      rm(cond_monthly_draws_used)
    }
  }
}

# combine condition draws with environmental covariates
cond_env_draws <- cond_monthly_draws_used %>%
  left_join(
    env_wide_rolling,
    by = c("year", "month")
  ) %>%
  arrange(.draw_index, year, month)

glimpse(cond_env_draws)

# List env variables to use in analysis
env_vars <- c(
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

# scale environmental variables
env_scaling <- env_wide_rolling %>%
  ungroup() %>%
  summarise(
    across(
      all_of(env_vars),
      list(
        mean = ~ mean(.x, na.rm = TRUE),
        sd = ~ sd(.x, na.rm = TRUE)
      )
    )
  )

# function to scale environmental variables
scale_env_variable <- function(x, variable, scaling_table) {

  variable_mean <- scaling_table[[paste0(variable, "_mean")]]

  variable_sd <- scaling_table[[paste0(variable, "_sd")]]

  as.numeric(
    (x - variable_mean) / variable_sd
  )
}

# scale environmental variables in the cond_env_draws data frame
for (v in env_vars) {
  cond_env_draws[[paste0(v, "_z")]] <-
    scale_env_variable(
      x = cond_env_draws[[v]],
      variable = v,
      scaling_table = env_scaling
    )
}

# -----------------------------
# 3.12 Create posterior-imputed datasets
# -----------------------------

set.seed(22)

# number of posterior draws to pull from condition model output
n_condition_datasets <- 50

# create list of posterior draws of condition that are matched to environmental variables
available_condition_draws <- cond_env_draws %>%
  distinct(.draw_index, .draw)

# select random sample of posterior draws to use for analysis
selected_condition_draws <- available_condition_draws %>%
  slice_sample(
    n = min(
      n_condition_datasets,
      nrow(available_condition_draws)
    )
  )

# create list of 50 dataset (from posterior draws) that are matched to environmental variables
environment_data_list <- selected_condition_draws$.draw_index %>%
  map(function(draw_number) {

    cond_env_draws %>%
      filter(.draw_index == draw_number) %>%
      mutate(
        month_f = factor(
          month,
          levels = 1:12,
          labels = month.abb
        ),

        year_f = factor(year),

        # Consecutive monthly time index
        time = year * 12L + month
      ) %>%
      arrange(time) %>%
      select(
        mean_log_Kn,
        year,
        year_f,
        month,
        month_f,
        time,
        all_of(env_vars),
        ends_with("_z")
      ) %>%
      drop_na(
        mean_log_Kn
      )
  })

length(environment_data_list)

# -----------------------------
# 3.13 General function for environmental models
# -----------------------------

# function for running all models with each environmental variable
fit_environment_model <- function(
    variable,
    data_list,
    model_file = NULL,
    seed = 22
) {

  variable_z <- paste0(variable, "_z")

  formula_text <- paste0(
    "mean_log_Kn ~ ",
    variable_z,
    " + month_f + (1 | year_f)"
  )

  model_formula <- bf(
    as.formula(formula_text)
  )

  brm_multiple(
    formula = model_formula,
    data = map(
      data_list,
      ~ .x %>%
        drop_na(
          mean_log_Kn,
          all_of(variable_z)
        )
    ),
    family = gaussian(),
    prior = priors_env,

    backend = "cmdstanr",
    chains = 4,
    iter = 1000,
    warmup = 500,
    seed = seed,

    control = list(
      adapt_delta = 0.95,
      max_treedepth = 12
    ),

    combine = TRUE,

    file = model_file,
    file_refit = "always"
    # file_refit = if (is.null(model_file)) NULL else "on_change"
  )
}

# set priors
priors_env <- c(
  prior(
    normal(0, 0.1),
    class = "Intercept"
  ),
  prior(
    normal(0, 0.1),
    class = "b"
  ),
  prior(
    student_t(3, 0, 0.1),
    class = "sd"
  ),
  prior(
    student_t(3, 0, 0.1),
    class = "sigma"
  )
)

# -----------------------------
# 3.14 Fit models for each environmental variable
# -----------------------------

# fit models for each environmental variable using the function above
env_models_bayes <- list(
  pdo_lagprev1month = fit_environment_model(
    variable = "pdo_lagprev1month",
    data_list = environment_data_list,
    model_file = paste0("model-outputs/m_env_pdo_lagprev1month_", selected)
  ),

  npgo_lagprev1month = fit_environment_model(
    variable = "npgo_lagprev1month",
    data_list = environment_data_list,
    model_file = paste0("model-outputs/m_env_npgo_lagprev1month_", selected)
  ),

  oni_lagprev1month = fit_environment_model(
    variable = "oni_lagprev1month",
    data_list = environment_data_list,
    model_file = paste0("model-outputs/m_env_oni_lagprev1month_", selected)
  ),

  npi_lagprev1month = fit_environment_model(
    variable = "npi_lagprev1month",
    data_list = environment_data_list,
    model_file = paste0("model-outputs/m_env_npi_lagprev1month_", selected)
  ),
  
  mei_lagprev1month = fit_environment_model(
    variable = "mei_lagprev1month",
    data_list = environment_data_list,
    model_file = paste0("model-outputs/m_env_mei_lagprev1month_", selected)
  ),

  pdo_lagprev1to3months = fit_environment_model(
    variable = "pdo_lagprev1to3months",
    data_list = environment_data_list,
    model_file = paste0("model-outputs/m_env_pdo_lagprev1to3months_", selected)
  ),

  npgo_lagprev1to3months = fit_environment_model(
    variable = "npgo_lagprev1to3months",
    data_list = environment_data_list,
    model_file = paste0("model-outputs/m_env_npgo_lagprev1to3months_", selected)
  ),

  oni_lagprev1to3months = fit_environment_model(
    variable = "oni_lagprev1to3months",
    data_list = environment_data_list,
    model_file = paste0("model-outputs/m_env_oni_lagprev1to3months_", selected)
  ),

  npi_lagprev1to3months = fit_environment_model(
    variable = "npi_lagprev1to3months",
    data_list = environment_data_list,
    model_file = paste0("model-outputs/m_env_npi_lagprev1to3months_", selected)
  ),

  mei_lagprev1to3months = fit_environment_model(
    variable = "mei_lagprev1to3months",
    data_list = environment_data_list,
    model_file = paste0("model-outputs/m_env_mei_lagprev1to3months_", selected)
  ),

  pdo_lagprev1to6months = fit_environment_model(
    variable = "pdo_lagprev1to6months",
    data_list = environment_data_list,
    model_file = paste0("model-outputs/m_env_pdo_lagprev1to6months_", selected)
  ),

  npgo_lagprev1to6months = fit_environment_model(
    variable = "npgo_lagprev1to6months",
    data_list = environment_data_list,
    model_file = paste0("model-outputs/m_env_npgo_lagprev1to6months_", selected)
  ),

  oni_lagprev1to6months = fit_environment_model(
    variable = "oni_lagprev1to6months",
    data_list = environment_data_list,
    model_file = paste0("model-outputs/m_env_oni_lagprev1to6months_", selected)
  ),

  npi_lagprev1to6months = fit_environment_model(
    variable = "npi_lagprev1to6months",
    data_list = environment_data_list,
    model_file = paste0("model-outputs/m_env_npi_lagprev1to6months_", selected)
  ),

  mei_lagprev1to6months = fit_environment_model(
    variable = "mei_lagprev1to6months",
    data_list = environment_data_list,
    model_file = paste0("model-outputs/m_env_mei_lagprev1to6months_", selected)
  )
)

# summarize

# slow try below before deleting
# env_effect_summary <- imap_dfr(
#   env_models_bayes,
#   function(model, variable) {
#
#     coefficient_name <- paste0(
#       "b_",
#       variable,
#       "_z"
#     )
#
#     model_draws <- as_draws_df(model)
#
#     beta <- model_draws[[coefficient_name]]
#
#     tibble(
#       covariate = variable,
#
#       median_beta = median(beta),
#       lower_beta = quantile(beta, 0.025),
#       upper_beta = quantile(beta, 0.975),
#
#       median_Kn_ratio = median(exp(beta)),
#       lower_Kn_ratio = quantile(exp(beta), 0.025),
#       upper_Kn_ratio = quantile(exp(beta), 0.975),
#
#       probability_positive = mean(beta > 0),
#       probability_negative = mean(beta < 0)
#     )
#   }
# )

# faster try below 
# summarize environmental effects from bayesian models
env_effect_summary <- imap_dfr(
  env_models_bayes,
  function(model, variable) {

    coefficient_name <- paste0(
      "b_",
      variable,
      "_z"
    )

    beta <- as_draws_matrix(
      model,
      variable = coefficient_name
    )[, coefficient_name]

    beta_exp <- exp(beta)

    tibble(
      covariate = variable,

      median_beta = median(beta),
      lower_beta = quantile(beta, 0.025),
      upper_beta = quantile(beta, 0.975),

      median_Kn_ratio = median(beta_exp),
      lower_Kn_ratio = quantile(beta_exp, 0.025),
      upper_Kn_ratio = quantile(beta_exp, 0.975),

      probability_positive = mean(beta > 0),
      probability_negative = mean(beta < 0)
    )
  }
)

env_effect_summary


# -----------------------------
# 3.15 Save posterior outputs
# -----------------------------

dir.create(
  "model-outputs/posterior-outputs",
  showWarnings = FALSE,
  recursive = TRUE
)

saveRDS(
  cond_monthly_draws,
  file = paste0("model-outputs/posterior-outputs/condition_monthly_posterior_draws_", selected, ".rds")
)

saveRDS(
  cond_monthly_summary,
  file = paste0("model-outputs/posterior-outputs/condition_monthly_posterior_summary_", selected, ".rds")
)

saveRDS(
  env_scaling,
  file = paste0("model-outputs/posterior-outputs/environment_covariate_scaling_", selected, ".rds")
)

write.csv(
  cond_monthly_summary,
  file = paste0("model-outputs/posterior-outputs/condition_monthly_posterior_summary_", selected, ".csv"),
  row.names = FALSE
)


# -----------------------------
# 3.16 Plot posterior outputs
# -----------------------------

# plot individual env variable posterior
plot_env_posterior <- function(model, variable) {

  coefficient_name <- paste0(
    "b_",
    variable,
    "_z"
  )

  beta <- posterior::as_draws_matrix(
    model,
    variable = coefficient_name
  )[, coefficient_name]

  beta <- as.numeric(beta)

  posterior_data <- tibble::tibble(
    beta = beta
  )

  variable_label <- variable |>
    gsub("_m_prev1$", ": previous month", x = _) |>
    gsub("_m_prev3$", ": previous 3 months", x = _) |>
    gsub("_m_prev6$", ": previous 6 months", x = _) |>
    toupper()

  ggplot2::ggplot(
    posterior_data,
    ggplot2::aes(
      x = beta,
      y = variable_label
    ) +
    ggdist::stat_halfeye(
      .width = c(0.50, 0.80, 0.95),
      point_interval = ggdist::median_qi,
      normalize = "groups"
    ) +
    ggplot2::geom_vline(
      xintercept = 0,
      linetype = "dashed",
      linewidth = 0.6
    ) +
    ggplot2::labs(
      title = paste(
        "Posterior effect of",
        variable_label
      ),
      subtitle = paste(
        "Coefficient:",
        coefficient_name
      ),
      x = expression(
        beta ~ "(change in log condition per 1 SD)"
      ),
      y = NULL
    ) +
    ggplot2::theme_classic(
      base_size = 12
    )
}

# plot example for PDO previous month
plot_env_posterior(
  model = env_models_bayes[["pdo_lagprev1month"]],
  variable = "pdo_lagprev1month"
)

# function to extract posterior draws for a given environmental model and variable
extract_env_posterior <- function(model, variable) {

  coefficient_name <- paste0(
    "b_",
    variable,
    "_z"
  )

  draws <- posterior::as_draws_matrix(
    model,
    variable = coefficient_name
  )

  if (ncol(draws) != 1L) {
    stop(
      paste0(
        "Expected coefficient '",
        coefficient_name,
        "' in model '",
        variable,
        "', but it was not found uniquely."
      ),
      call. = FALSE
    )
  }

  tibble::tibble(
    model = variable,
    coefficient = coefficient_name,
    beta = as.numeric(draws[, 1])
  )
}

# extract posterior draws for all environmental models
env_posterior_draws <- imap_dfr(
  env_models_bayes,
  \(model, variable) {
    extract_env_posterior(
      model = model,
      variable = variable
    )
  }
)

# add index and exposure window information to the posterior draws
env_posterior_draws <- env_posterior_draws %>%
  mutate(
    index = case_when(
      grepl("^pdo", model)  ~ "PDO",
      grepl("^npgo", model) ~ "NPGO",
      grepl("^oni", model)  ~ "ONI",
      grepl("^npi", model)  ~ "NPI",
      grepl("^mei", model)  ~ "MEI",
      TRUE                  ~ model
    ),

    exposure_window = case_when(
      grepl("_lagprev1month", model) ~ "Previous month",
      grepl("_lagprev1to3", model) ~ "Previous 3 months",
      grepl("_lagprev1to6", model) ~ "Previous 6 months",
      TRUE                    ~ "Other"
    ),

    exposure_window = factor(
      exposure_window,
      levels = c(
        "Previous month",
        "Previous 3 months",
        "Previous 6 months"
      )
    ),

    index = factor(
      index,
      levels = c("PDO", "NPGO", "ONI", "NPI", "MEI")
    )
  )

#####
# plot combined posterior effects for all environmental variables
env_posterior_combined_plot <- ggplot(
  env_posterior_draws,
  aes(
    x = beta,
    y = exposure_window
  )
) +
  stat_halfeye(
    .width = c(0.50, 0.80, 0.95),
    point_interval = median_qi,
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
    title = "Posterior environmental effects on body condition",
    subtitle = paste(
      "Points show posterior medians;",
      "intervals show 50%, 80%, and 95% credible intervals"
    ),
    x = expression(beta ~ "(change in log condition per 1 SD)"),
    y = "Environmental exposure window"
  ) +
  theme_classic(base_size = 12) +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(face = "bold"),
    axis.text.y = element_text(size = 10)
  )

env_posterior_combined_plot

#####
#  plot percent change in body condition
env_posterior_draws <- env_posterior_draws %>%
  mutate(
    percent_change = 100 * (exp(beta) - 1)
  )

env_percent_change_plot <- ggplot(
  env_posterior_draws,
  aes(
    x = percent_change,
    y = exposure_window
  )
) +
  stat_halfeye(
    .width = c(0.50, 0.80, 0.95),
    point_interval = median_qi,
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
    title = "Estimated changes in body condition",
    subtitle = "Change associated with a 1-SD increase in each index",
    x = "Change in body condition (%)",
    y = "Environmental exposure window"
  ) +
  theme_classic(base_size = 12) +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(face = "bold"),
    axis.text.y = element_text(size = 10)
  )

env_percent_change_plot


# funtion to plot linear relationship
plot_env_linear <- function(
    variable,
    data = dat_env,
    posterior_draws = env_posterior_draws,
    response = "mean_log_Kn",
    n_grid = 100,
    scaling_table = env_scaling
) {

  # Extract existing posterior slopes
  variable_draws <- posterior_draws %>%
    filter(model == variable)

  beta_z <- variable_draws %>%
    pull(beta)

  if (length(beta_z) == 0) {
    stop("No posterior draws found for ", variable)
  }

  variable_info <- variable_draws %>%
    distinct(index, exposure_window)

  if (nrow(variable_info) != 1L) {
    stop("Expected one index and exposure window for ", variable)
  }

  index_label <- as.character(variable_info$index[[1]])
  exposure_window <- as.character(variable_info$exposure_window[[1]])

  window_subtitle <- dplyr::case_when(
    exposure_window == "Previous month" ~ "Previous month (1-month lag)",
    exposure_window == "Previous 3 months" ~
      "Previous 1-3 months (3-month mean)",
    exposure_window == "Previous 6 months" ~
      "Previous 1-6 months (6-month mean)",
    TRUE ~ exposure_window
  )

  # Observed data
  observed <- data %>%
    transmute(
      x = .data[[variable]],
      y = .data[[response]]
    ) %>%
    filter(
      is.finite(x),
      is.finite(y)
    )

  # IMPORTANT:
  # These must match the scaling used to fit the model.
  mean_x <- mean(observed$x)
  sd_x <- scaling_table[[paste0(variable, "_sd")]]

  # Convert standardized posterior slopes to original units
  if (length(sd_x) != 1L || !is.finite(sd_x) || sd_x <= 0) {
    stop("Invalid model scaling standard deviation for ", variable)
  }

  beta_original <- beta_z / sd_x

  # Anchor the regression lines at the observed means
  mean_y <- mean(observed$y)

  # Grid of environmental values
  x_grid <- seq(
    min(observed$x),
    max(observed$x),
    length.out = n_grid
  )

  # Compute posterior expected relationships
  # Each row corresponds to an environmental value
  prediction_summary <- lapply(
    x_grid,
    function(x) {

      predicted_y <- mean_y +
        beta_original * (x - mean_x)

      tibble(
        x = x,
        median = median(predicted_y),
        lower = quantile(predicted_y, 0.025),
        upper = quantile(predicted_y, 0.975)
      )
    }
  ) %>%
    bind_rows()

  # Plot
  ggplot() +
    geom_point(
      data = observed,
      aes(x = x, y = y),
      alpha = 0.4,
      size = 1.5
    ) +
    geom_ribbon(
      data = prediction_summary,
      aes(
        x = x,
        ymin = lower,
        ymax = upper
      ),
      alpha = 0.2
    ) +
    geom_line(
      data = prediction_summary,
      aes(
        x = x,
        y = median
      ),
      linewidth = 0.9
    ) +
    labs(
    #   title = index_label,
      subtitle = window_subtitle,
      x = index_label,
      y = "Mean log body condition"
    ) +
    theme_classic(base_size = 12)
}

# create plots for each environmental variable posterior linear relationship
env_linear_variables <- env_posterior_draws %>%
  distinct(model, index, exposure_window) %>%
  arrange(index, exposure_window) %>%
  pull(model)

env_linear_plots <- map(
  env_linear_variables,
  plot_env_linear
)

# combine plots into one figure
env_linear_combined <- wrap_plots(
  env_linear_plots,
  ncol = 3,
  byrow = TRUE
)

# these are the plots of the linear relationships between each environmental covariate and body condition, with posterior uncertainty shown as ribbons around the median predicted line. The points represent the observed data.
# environmental data are on original scale, and the predicted lines are based on the posterior draws of the slopes from the Bayesian models, anchored at the mean of the observed data.
env_linear_combined


##### -------------------------------------------------------------------- #####
##### -------------------------------------------------------------------- #####
##### -------------------------------------------------------------------- #####
##### -------------------------------------------------------------------- #####
##### -------------------------------------------------------------------- #####

# Seasonal body condition responses to environmental covariates

# -----------------------------
# 3.21 Join seasonal condition and environmental covariates
# -----------------------------

# function to summarise environmental covariates over a specified window of months
summarise_environment_window <- function(data, months, label) {
  data %>%
    filter(month %in% (months + 1)) %>%
    group_by(year) %>%
    summarise(
      across(
        matches("^(pdo|npgo|oni|npi|mei)_lagprev1month$"),
        ~ if (sum(!is.na(.x)) == length(months)) {
          mean(.x, na.rm = TRUE)
        } else {
          NA_real_
        }
      ),
      .groups = "drop"
    ) %>%
    rename_with(
      ~ str_replace(.x, "_lagprev1month$", paste0("_", label)),
      matches("^(pdo|npgo|oni|npi|mei)_lagprev1month$")
    )
}

# summarise environmental covariates for each seasonal window
summarise_environment <- full_join(
  summarise_environment_window(env_wide_rolling, 1:3, "jantomar"),  
  summarise_environment_window(env_wide_rolling, 4:6, "aprtojun"),
  by = "year") %>% 
  full_join(
    summarise_environment_window(env_wide_rolling, 5:7, "maytojul"),
    by = "year") %>%
  full_join(
    summarise_environment_window(env_wide_rolling, 8:10, "augtooct"),
    by = "year")

# summarise seasonal body condition for spring (april-june) and late summer/early fall (august-october) periods
summarise_condition <- cond_monthly %>%
  group_by(year) %>%
  summarise(
    mean_log_Kn_aprtojun = mean(mean_log_Kn[month %in% 4:6], na.rm = TRUE),
    mean_log_Kn_augtooct = mean(mean_log_Kn[month %in% 8:10], na.rm = TRUE),
    .groups = "drop"
  )

summarise_condition <- cond_monthly %>%
  filter(month %in% 4:6) %>%
  group_by(year) %>%
  summarise(
    mean_log_Kn = mean(mean_log_Kn, na.rm = TRUE),
    n_fish = sum(n_fish),
    .groups = "drop"
  ) %>% 
  mutate(window = "aprtojun") %>%
  bind_rows(
    cond_monthly %>%
      filter(month %in% 8:10) %>%
      group_by(year) %>%
      summarise(
        mean_log_Kn = mean(mean_log_Kn, na.rm = TRUE),
        n_fish = sum(n_fish),
        .groups = "drop"
      ) %>%
    mutate(window = "augtooct")
  )

# join seasonal condition and environmental covariates
seasonal_joint <- summarise_condition %>%
  left_join(
    summarise_environment,
    by = c("year")
  ) %>%
  arrange(year) %>%
  mutate(   
    year_z = as.numeric(scale(year)),
    condition_z = as.numeric(scale(mean_log_Kn))
  )

glimpse(seasonal_joint)


# -----------------------------
# 3.22 Plot seasonal body condition with environmental covariates
# -----------------------------

# convert to long format for plotting
seasonal_joint_long <- seasonal_joint %>%
  filter(window == "aprtojun") %>%
  select(-n_fish, -window, -year_z, -condition_z) %>%
  pivot_longer(
    cols = c(-year, -mean_log_Kn),
    names_to = c("index", "window"),
    names_sep = "_",      
    values_to = "value")  %>%
  mutate(
    window = recode(
      window,
      jantomar = "Jan to Mar",
      aprtojun = "Apr to Jun",
      maytojul = "May to Jul",
      augtooct = "Aug to Oct"
    ),
    window = factor(
      window,
      levels = c("Jan to Mar", "Apr to Jun", "May to Jul", "Aug to Oct")
    )
  ) %>%
  group_by(index, window) %>%
  mutate(value_z = as.numeric(scale(value))) %>%
  ungroup()

# plot seasonal environmental covariates over time
ggplot(seasonal_joint_long, aes(year, value_z, colour = index)) +
  facet_grid(index ~ window) +
  geom_line(alpha = 0.6) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.5) +
  geom_point(size = 1.5) +
  geom_smooth(
    method = "loess",
    se = FALSE,
    linewidth = 1,
  ) +
  labs(
    x = "Year",
    y = "Standardized value",
    colour = NULL,
    title = "Temporal patterns in seasonal environmental covariates"
  ) +
  theme_classic()

# plot spring body condition and environmental covariates
spring_join_long <- seasonal_joint %>%
  filter(window == "aprtojun") %>%
  select(-n_fish, -window, -year_z, -condition_z, -contains("maytojul"), -contains("augtooct")) %>%
  pivot_longer(
    cols = c(-year, -mean_log_Kn),
    names_to = c("index", "window"),
    names_sep = "_",      
    values_to = "value") %>%
  mutate(
    window = recode(
      window,
      jantomar = "Jan to Mar",
      aprtojun = "Apr to Jun",
    ),
    window = factor(
      window,
      levels = c("Jan to Mar", "Apr to Jun")
    )
  ) %>%
  group_by(index, window) %>%
  mutate(value_z = as.numeric(scale(value))) %>%
  ungroup()

ggplot(spring_join_long, aes(value_z, mean_log_Kn, colour = index)) +
  facet_grid(index ~ window) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.5) +
  geom_point(size = 1.5) +
  geom_smooth(
    method = "lm",
    se = FALSE,
    linewidth = 1,
  ) +
  labs(
    x = "Standardized value",
    y = "Standardized body condition",
    colour = NULL,
    title = "Temporal patterns in spring body condition and environmental covariates"
  ) +
  theme_classic()


# plot fall body condition and environmental covariates over time
fall_join_long <- seasonal_joint %>%
  filter(window == "augtooct") %>%
  select(-n_fish, -window, -year_z, -condition_z, -contains("jantomar"), -contains("aprtojun")) %>%
  pivot_longer(
    cols = c(-year, -mean_log_Kn),
    names_to = c("index", "window"),
    names_sep = "_",      
    values_to = "value") %>%
  mutate(
    window = recode(
      window,
      maytojul = "May to Jul",
      augtooct = "Aug to Oct"
    ),
    window = factor(
      window,
      levels = c("May to Jul", "Aug to Oct")
    )
  ) %>%
  group_by(index, window) %>%
  mutate(value_z = as.numeric(scale(value))) %>%
  ungroup()

ggplot(fall_join_long, aes(value_z, mean_log_Kn, colour = index)) +
  facet_grid(index ~ window) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.5) +
  geom_point(size = 1.5) +
  geom_smooth(
    method = "lm",
    se = FALSE,
    linewidth = 1,
  ) +
  labs(
    x = "Standardized value",
    y = "Standardized body condition",
    colour = NULL,
    title = "Temporal patterns in fall body condition and environmental covariates"
  ) +
  theme_classic()

# -----------------------------
# 3.23 Environment effects on body condition
# -----------------------------

# Select environmental indices
spring_env_vars <- c(
  "pdo_jantomar",
  "npgo_jantomar",
  "oni_jantomar",
  "npi_jantomar",
  "mei_jantomar",
  "pdo_aprtojun",
  "npgo_aprtojun",
  "oni_aprtojun",
  "npi_aprtojun",
  "mei_aprtojun"
)

# Fit one model per environmental variable
spring_condition_models <- set_names(
  map(spring_env_vars, function(v) {

    dat <- seasonal_joint %>%
      filter(window == "aprtojun") %>%
      select(year, year_z, mean_log_Kn, all_of(v)) %>%
      drop_na() %>%
      mutate(env_z = as.numeric(scale(.data[[v]])))

    brm(
      mean_log_Kn ~ env_z + s(year, k = 8),
      # mean_log_Kn ~ env_z,
      data = dat,
      family = student(),
      backend = "cmdstanr",
      chains = 4,
      iter = 3000,
      warmup = 1500,
      seed = 22
    )

  }),
  spring_env_vars
)

# Select environmental indices
fall_env_vars <- c(
  "pdo_maytojul",
  "npgo_maytojul",
  "oni_maytojul",
  "npi_maytojul",
  "mei_maytojul",
  "pdo_augtooct",
  "npgo_augtooct",
  "oni_augtooct",
  "npi_augtooct",
  "mei_augtooct"
)

# Fit one model per environmental variable
fall_condition_models <- set_names(
  map(fall_env_vars, function(v) {

    dat <- seasonal_joint %>%
      filter(window == "augtooct") %>%
      select(year, year_z, mean_log_Kn, all_of(v)) %>%
      drop_na() %>%
      mutate(env_z = as.numeric(scale(.data[[v]])))

    brm(
      mean_log_Kn ~ env_z + s(year, k = 8),
      # mean_log_Kn ~ env_z,
      data = dat,
      family = student(),
      backend = "cmdstanr",
      chains = 4,
      iter = 3000,
      warmup = 1500,
      seed = 22
    )

  }),
  fall_env_vars
)

# combine spring and fall results
full_condition_models <- c(
    spring_condition_models,
    fall_condition_models)

# -----------------------------
# 3.24 Summarize posterior effects of environmental covariates on seasonal body condition
# -----------------------------

# extract posterior draws for all environmental models
full_condition_draws <- purrr::imap_dfr(
  full_condition_models,
  function(model, variable) {

    draws <- posterior::as_draws_matrix(
      model,
      variable = "b_env_z"
    )

    tibble(
      covariate = variable,
      beta = as.numeric(draws[, 1])
    )
  }
) %>%
  mutate(
    index = toupper(sub("_.*", "", covariate)),
    window = sub(".*_", "", covariate),
    season = case_when(
      window %in% c("jantomar", "aprtojun") ~ "Spring",
      window %in% c("maytojul", "augtooct") ~ "Fall",
      TRUE ~ NA_character_
    )
  )

# Summarize posterior effects
full_condition_posterior_summary <- full_condition_draws %>%
  group_by(covariate, index, window, season) %>%
  summarise(
    median_beta = median(beta),
    lower_95 = quantile(beta, 0.025),
    upper_95 = quantile(beta, 0.975),
    prob_positive = mean(beta > 0),
    prob_negative = mean(beta < 0),

    # Convert log-condition effects to percentage change
    median_pct = 100 * (exp(median_beta) - 1),
    lower_pct = 100 * (exp(lower_95) - 1),
    upper_pct = 100 * (exp(upper_95) - 1),

    .groups = "drop"
  )

full_condition_posterior_summary

# -----------------------------
# 3.25 Plot posterior effects of environmental covariates on seasonal body condition
# -----------------------------

full_condition_draws <- full_condition_draws %>%
  mutate(
    pct_change = 100 * (exp(beta) - 1),
    index = factor(
      index,
      levels = c("PDO", "NPGO", "ONI", "NPI", "MEI")
    ),
    window = factor(
      window,
      levels = c("augtooct", "maytojul", "aprtojun", "jantomar")
    ),
    season = case_when(
      window %in% c("jantomar", "aprtojun") ~ "Spring",
      window %in% c("maytojul", "augtooct") ~ "Fall",
      TRUE ~ NA_character_
    ),
    season = factor(
      season,
      levels = c("Spring", "Fall")
    )
  )

p_full_condition_posterior <- ggplot(
  full_condition_draws,
  aes(
    x = pct_change,
    y = window,
    fill = season
  )
) + 
  scale_y_discrete(
    labels = c(
      jantomar = "Jan to Mar",
      aprtojun = "Apr to Jun",
      maytojul = "May to Jul",
      augtooct = "Aug to Oct"
    )
  ) +
  facet_wrap(
    ~ index,
    ncol = 2,
    scales = "free_x"
  ) +
  stat_halfeye(
    .width = c(0.50, 0.80, 0.95),
    point_interval = median_qi,
    normalize = "groups",
    slab_colour = "black",
    slab_linewidth = 0.3
  ) +
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    linewidth = 0.6
  ) +
  labs(
    # title = "Environmental effects on Seasonal body condition",
    # subtitle = "Adjusted for temporal trends",
    x = "Change in body condition (%) per 1-SD increase in index",
    y = NULL
  ) +
  theme_classic(base_size = 12)

p_full_condition_posterior

# save plot
ggsave(
  filename = "figures/seasonal_condition_environment_effects.png",
  plot = p_full_condition_posterior,
  width = 10,
  height = 8,
  dpi = 300
)


# -----------------------------
# 3.26 Observed and fitted relationships between seasonal body condition and environmental covariates
# -----------------------------

get_full_condition_relationship <- function(
    model,
    variable,
    data = seasonal_joint,
    season = c("aprtojun", "augtooct"),
    ndraws = 1000
) {

  # Original data used for this predictor
  dat <- data %>%
    filter(window %in% season) %>%
    select(
      year,
      mean_log_Kn,
      all_of(variable)
    ) %>%
    drop_na() %>%
    mutate(
      env_z = as.numeric(
        scale(.data[[variable]])
      )
    )

  # Environmental grid
  env_grid <- seq(
    min(dat$env_z),
    max(dat$env_z),
    length.out = 100
  )

  # Hold year constant to isolate environmental effect
  newdata <- tibble(
    env_z = env_grid,
    year = median(dat$year)
  )

  # Posterior expected values
  epred <- posterior_epred(
    model,
    newdata = newdata,
    ndraws = ndraws,
    re_formula = NA
  )

  # Summarize expected values
  pred_summary <- tibble(
    env_z = env_grid,
    median = apply(epred, 2, median),
    lower_95 = apply(
      epred, 2,
      \(x) unname(quantile(x, 0.025))
    ),
    upper_95 = apply(
      epred, 2,
      \(x) unname(quantile(x, 0.975))
    )
  )

  # Return both observed data and predictions
  list(
    observed = dat %>%
      mutate(covariate = variable,
             season = case_when(
               season == "aprtojun" ~ "Spring",
               season == "augtooct" ~ "Fall",
               TRUE ~ NA_character_
             ),
             index = toupper(sub("_.*", "", variable)),
             window = sub(".*_", "", variable)),

    predicted = pred_summary %>%
      mutate(covariate = variable,
             season = case_when(
               season == "aprtojun" ~ "Spring",
               season == "augtooct" ~ "Fall",
               TRUE ~ NA_character_
             ),
             index = toupper(sub("_.*", "", variable)),
             window = sub(".*_", "", variable)))
}

# Run for all environmental variables
spring_condition_relationships <- imap(
  spring_condition_models,
  ~ get_full_condition_relationship(
    model = .x,
    variable = .y,
    season = "aprtojun"
  )
)

length(spring_condition_relationships)
str(spring_condition_relationships[[1]])


fall_condition_relationships <- imap(
  fall_condition_models,
  ~ get_full_condition_relationship(
    model = .x,
    variable = .y,
    season = "augtooct"
  )
)

# combine spring and fall relationships
full_condition_relationships <- c(
  spring_condition_relationships,
  fall_condition_relationships
)

# Combine observations
full_condition_observed <- map_dfr(
  full_condition_relationships,
  "observed"
)

# Combine predictions
full_condition_predicted <- map_dfr(
  full_condition_relationships,
  "predicted"
)

# order the levels of the window factor for consistent plotting
window_levels <- c("jantomar", "aprtojun", "maytojul", "augtooct")

full_condition_observed <- full_condition_observed %>%
  mutate(window = factor(window, levels = window_levels),
         season = factor(season, levels = c("Spring", "Fall")))

full_condition_predicted <- full_condition_predicted %>%
  mutate(window = factor(window, levels = window_levels),
         season = factor(season, levels = c("Spring", "Fall")))

# -----------------------------
# 3.27 Plot seasonal body condition relationships with environmental covariates
# -----------------------------

p_full_condition_relationships <- ggplot() +

  # horizonal line at zero
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    linewidth = 0.5
  ) +

  # 95% credible interval
  geom_ribbon(
    data = full_condition_predicted,
    aes(
      x = env_z,
      ymin = lower_95,
      ymax = upper_95
    ),
    alpha = 0.2
  ) +

  # Median fitted relationship
  geom_line(
    data = full_condition_predicted,
    aes(
      x = env_z,
      y = median
    ),
    linewidth = 1
  ) +

  # Observed annual values
  geom_point(
    data = full_condition_observed,
    aes(
      x = env_z,
      y = mean_log_Kn,
      colour = season
    ),
    alpha = 0.5,
    size = 1.8
  ) +

  facet_grid(
    rows = vars(index),
    cols = vars(window),
    labeller = labeller(
      window = c(
        jantomar = "Jan to Mar",
        aprtojun = "Apr to Jun",
        maytojul = "May to Jul",
        augtooct = "Aug to Oct"
      )
    ),
    scales = "free",
    # ncol = 3
  ) +

  labs(
    # title = "Seasonal body condition and environmental indices",
    # subtitle = 
    #   "Bayesian fitted relationships adjusted for nonlinear temporal trends",
    x = "Standardized environmental index",
    y = "Seasonal mean log condition"
  ) +

  theme_classic(base_size = 12)

p_full_condition_relationships

# save plot
ggsave(
  filename = "figures/seasonal_condition_environment_relationships.png",
  plot = p_full_condition_relationships,
  width = 12,
  height = 10,
  dpi = 300
)


# -----------------------------
# 3.28 Forest plot of posterior effects of environmental covariates on seasonal body condition
# -----------------------------

p_full_condition_forest <- full_condition_posterior_summary %>%
  mutate(
    window = factor(window, levels = c("augtooct", "maytojul", "aprtojun", "jantomar")),
    season = factor(season, levels = c("Spring", "Fall"))
  ) %>%
  ggplot(
    aes(
      x = median_pct,
      y = window,
      colour = season
    )
  ) +
  geom_vline(
    xintercept = 0,
    linetype = "dashed"
  ) +
  geom_errorbar(
    aes(
      xmin = lower_pct,
      xmax = upper_pct,
      colour = season
    ),
    height = 0.15
  ) +
  geom_point(
    size = 2.5
  ) +
  facet_wrap(
    ~ index,
    ncol = 2,
    scales = "free_x"
  ) +
  scale_y_discrete(
    labels = c(
      jantomar = "Jan to Mar",
      aprtojun = "Apr to Jun",
      maytojul = "May to Jul",
      augtooct = "Aug to Oct"
    )
  ) +
  labs(
    # title = "Environmental effects on annual body condition",
    # subtitle = "Posterior mean and 95% credible intervals",
    x = "Change in body condition (%) per 1-SD increase in index",
    y = NULL
  ) +
  theme_classic()

p_full_condition_forest

# save plot
ggsave(
  filename = "figures/seasonal_condition_environment_forest.png",
  plot = p_full_condition_forest,
  width = 10,
  height = 8,
  dpi = 300
)

