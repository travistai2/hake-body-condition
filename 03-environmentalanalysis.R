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
env_pdo <- format_covariate_month(pdo, "anomaly", "pdo_lagprev1month")
env_npgo <- format_covariate_month(npgo, "anomaly", "npgo_lagprev1month")
env_oni <- format_covariate_month(oni, "anomaly", "oni_lagprev1month")
env_npi <- format_covariate_month(npi_monthly , "value", "npi_lagprev1month")
env_mei <- format_covariate_month(mei , "anomaly", "mei_lagprev1month")

# combine variables
env_long <- bind_rows(env_pdo, env_npgo, env_oni, env_npi, env_mei)

##
# seasonal covariates of previous 1, 3, and 6 months - combined into one data set

# estimate prev 1 variables first
lag_months <- 1
env_long <- env_long %>%
  mutate(month = month + lag_months) %>%
  mutate(year = case_when(month > 12 ~ year + 1,
                          .default = year),
         month = case_when(month > 12 ~ (month - 12),
                           .default = month))

# make data wide
env_wide <- env_long %>%
  pivot_wider(names_from = covariate, values_from = value)

# use lag 1 variable names to then find the average across 3 and 6 months
lag1_vars <- c(
  "pdo_lagprev1month",
  "npgo_lagprev1month",
  "oni_lagprev1month",
  "npi_lagprev1month",
  "mei_lagprev1month"
)

# formatting data and making 3 and 6 month rolling averages
env_wide_rolling <- env_wide %>%
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

# join with body condition data
dat_env <- cond_monthly %>%
  left_join(env_wide_rolling, by = c("year", "month"))


# -----------------------------
# 3.2 Correlations with environmental covariates
# -----------------------------

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
# 1.10 Model-based predicted plots
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
    )
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

plot_env_posterior(
  model = env_models_bayes[["pdo_m_prev1"]],
  variable = "pdo_m_prev1"
)

# all plots in facet
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

env_posterior_draws <- imap_dfr(
  env_models_bayes,
  \(model, variable) {
    extract_env_posterior(
      model = model,
      variable = variable
    )
  }
)

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
      grepl("_prev1$", model) ~ "Previous month",
      grepl("_prev3$", model) ~ "Previous 3 months",
      grepl("_prev6$", model) ~ "Previous 6 months",
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


# percent change in body condition
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
    n_grid = 100
) {

  # Extract existing posterior slopes
  beta_z <- posterior_draws %>%
    filter(model == variable) %>%
    pull(beta)

  if (length(beta_z) == 0) {
    stop("No posterior draws found for ", variable)
  }

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
  sd_x <- sd(observed$x)

  # Convert standardized posterior slopes to original units
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
      title = toupper(variable),
      subtitle = "Posterior median slope and 95% credible interval",
      x = "Environmental index (original units)",
      y = "Mean log body condition"
    ) +
    theme_classic(base_size = 12)
}

env_linear_plots <- map(
  unique(env_posterior_draws$model),
  # unique(env_posterior_draws$model)[order(unique(env_posterior_draws$model))],
  plot_env_linear
)

env_linear_combined <- wrap_plots(
  env_linear_plots,
  ncol = 3
)

env_linear_combined




















##### -------------------------------------------------------------------- #####
##### -------------------------------------------------------------------- #####
##### -------------------------------------------------------------------- #####
##### -------------------------------------------------------------------- #####
##### -------------------------------------------------------------------- #####

# -----------------------------
# 6.1 Annual condition and environment
# -----------------------------


##### OPTION 1
# April-June condition, giving each month equal weight
annual_condition <- cond_monthly %>%
  filter(month %in% 4:6) %>%
  group_by(year) %>%
  summarise(
    mean_log_Kn = mean(mean_log_Kn, na.rm = TRUE),
    n_months = sum(!is.na(month)),
    n_fish = sum(n_fish),
    .groups = "drop"
  ) #%>%
  # filter(n_months == 3)

# April-June environmental conditions
# Use existing lagged indices
annual_environment <- env_wide_rolling %>%
  filter(month %in% (4:6 + 1)) %>% # add one as we have a lagfor each month
  group_by(year) %>%
  summarise(
    across(
      matches("^(pdo|npgo|oni|npi|mei)_m_prev[1]$"),
      ~ if (sum(!is.na(.x)) == 3) {
        mean(.x)
      } else {
        NA_real_
      }
    ),
    .groups = "drop"
  ) %>%
  rename_with(~ str_replace(.x, "_m_prev1$", "_annual_same"))


##### OPTION 2 -
# Aug-Oct condition, giving each month equal weight
annual_condition <- cond_monthly %>%
  filter(month %in% 8:10) %>%
  group_by(year) %>%
  summarise(
    mean_log_Kn = mean(mean_log_Kn, na.rm = TRUE),
    n_months = sum(!is.na(month)),
    n_fish = sum(n_fish),
    .groups = "drop"
  ) #%>%
# filter(n_months == 3)

# Aug-Oct environmental conditions
# Use existing lagged indices
annual_environment <- env_wide_rolling %>%
  filter(month %in% (8:10 + 1)) %>% # add one as we have a lagfor each month
  group_by(year) %>%
  summarise(
    across(
      matches("^(pdo|npgo|oni|npi|mei)_m_prev[1]$"),
      ~ if (sum(!is.na(.x)) == 3) {
        mean(.x)
      } else {
        NA_real_
      }
    ),
    .groups = "drop"
  ) %>%
  rename_with(~ str_replace(.x, "_m_prev1$", "_annual_same"))



# Join condition, environment, and age
annual_joint <- annual_condition %>%
  left_join(annual_environment, by = "year") %>%
  inner_join(
    annual_age_env %>%
      transmute(
        year = as.integer(Year),
        mean_age = mean_age
      ),
    by = "year"
  ) %>%
  arrange(year) %>%
  mutate(
    year_z = as.numeric(scale(year)),
    condition_z = as.numeric(scale(mean_log_Kn))
  )

glimpse(annual_joint)


# -----------------------------
# 6.2 plot age and condition together
# -----------------------------

annual_joint %>%
  select(year, mean_age, mean_log_Kn) %>%
  pivot_longer(
    cols = c(mean_age, mean_log_Kn),
    names_to = "response",
    values_to = "value"
  ) %>%
  group_by(response) %>%
  mutate(value_z = as.numeric(scale(value))) %>%
  ungroup() %>%
  ggplot(aes(year, value_z, colour = response)) +
  geom_line(alpha = 0.6) +
  geom_point(size = 1.5) +
  geom_smooth(
    method = "loess",
    se = FALSE,
    linewidth = 1
  ) +
  labs(
    x = "Year",
    y = "Standardized value",
    colour = NULL,
    title = "Temporal patterns in fish age and body condition"
  ) +
  theme_classic()




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



# -----------------------------
# 6.4 Environment effects on body condition
# -----------------------------

# Select environmental indices
annual_env_vars <- c(
  "pdo_annual_same",
  # "pdo_annual_prior",
  "npgo_annual_same",
  # "npgo_annual_prior",
  "oni_annual_same",
  # "oni_annual_prior",
  "npi_annual_same",
  # "npi_annual_prior",
  "mei_annual_same"
  # "mei_annual_prior"
)

# Fit one model per environmental variable
annual_condition_models <- set_names(
  map(annual_env_vars, function(v) {

    dat <- annual_joint %>%
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
  annual_env_vars
)

# Summarize environmental slopes
annual_condition_results <- imap_dfr(
  annual_condition_models,
  function(model, variable) {

    x <- fixef(model)["env_z", ]

    tibble(
      covariate = variable,
      estimate = x["Estimate"],
      lower = x["Q2.5"],
      upper = x["Q97.5"]
    )
  }
)

annual_condition_results


# -----------------------------
# 6.4.1 Extract posterior slopes
# -----------------------------

library(posterior)
library(ggdist)
library(patchwork)

annual_condition_draws <- purrr::imap_dfr(
  annual_condition_models,
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
    index = toupper(
      sub("_.*", "", covariate)
    )
  )

# Summarize posterior effects
annual_condition_posterior_summary <- annual_condition_draws %>%
  group_by(covariate, index) %>%
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

annual_condition_posterior_summary

# -----------------------------
# 6.4.2 Posterior effect size plot
# -----------------------------

annual_condition_draws <- annual_condition_draws %>%
  mutate(
    pct_change = 100 * (exp(beta) - 1),
    index = factor(
      index,
      levels = c("PDO", "NPGO", "ONI", "NPI", "MEI")
    )
  )

p_annual_condition_posterior <- ggplot(
  annual_condition_draws,
  aes(
    x = pct_change,
    y = index
  )
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
    title = "Environmental effects on annual body condition",
    subtitle = paste(
      "April-June environmental conditions;",
      "adjusted for temporal trends"
    ),
    x = "Change in body condition (%) per 1-SD increase in index",
    y = NULL
  ) +
  theme_classic(base_size = 12)

p_annual_condition_posterior



# -----------------------------
# 6.4.3 Bayesian environmental relationships
# -----------------------------

get_annual_condition_relationship <- function(
    model,
    variable,
    data = annual_joint,
    ndraws = 1000
) {

  # Original data used for this predictor
  dat <- data %>%
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
      mutate(covariate = variable),

    predicted = pred_summary %>%
      mutate(covariate = variable)
  )
}

# Run for all environmental variables
annual_condition_relationships <- imap(
  annual_condition_models,
  ~ get_annual_condition_relationship(
    model = .x,
    variable = .y
  )
)

# Combine observations
annual_condition_observed <- map_dfr(
  annual_condition_relationships,
  "observed"
)

# Combine predictions
annual_condition_predicted <- map_dfr(
  annual_condition_relationships,
  "predicted"
)


# -----------------------------
# 6.4.4 Plot all environmental relationships
# -----------------------------

p_annual_condition_relationships <- ggplot() +

  # 95% credible interval
  geom_ribbon(
    data = annual_condition_predicted,
    aes(
      x = env_z,
      ymin = lower_95,
      ymax = upper_95
    ),
    alpha = 0.2
  ) +

  # Median fitted relationship
  geom_line(
    data = annual_condition_predicted,
    aes(
      x = env_z,
      y = median
    ),
    linewidth = 1
  ) +

  # Observed annual values
  geom_point(
    data = annual_condition_observed,
    aes(
      x = env_z,
      y = mean_log_Kn
    ),
    alpha = 0.5,
    size = 1.8
  ) +

  facet_wrap(
    ~ covariate,
    scales = "free_x",
    ncol = 3,
    labeller = labeller(
      covariate = function(x) {
        toupper(sub("_annual_same$", "", x))
      }
    )
  ) +

  labs(
    title = "Annual body condition and environmental indices",
    subtitle = paste(
      "Bayesian fitted relationships adjusted",
      "for nonlinear temporal trends"
    ),
    x = "Standardized environmental index",
    y = "Annual mean log condition"
  ) +

  theme_classic(base_size = 12)

p_annual_condition_relationships


# -----------------------------
# 6.4.5 Combine plots
# -----------------------------

p_annual_condition_combined <-
  p_annual_condition_relationships /
  p_annual_condition_posterior +
  plot_layout(heights = c(2, 1))

p_annual_condition_combined

# ggsave(
#   filename = "annual_condition_environment_effects.png",
#   plot = p_annual_condition_combined,
#   width = 14,
#   height = 12,
#   dpi = 300
# )


# -----------------------------
# 6.4.6 Forest plots
# -----------------------------

p_annual_condition_forest <- annual_condition_results %>%
  mutate(
    index = toupper(
      sub("_.*", "", covariate)
    )
  ) %>%
  ggplot(
    aes(
      x = estimate,
      y = reorder(index, estimate)
    )
  ) +
  geom_vline(
    xintercept = 0,
    linetype = "dashed"
  ) +
  geom_errorbarh(
    aes(
      xmin = lower,
      xmax = upper
    ),
    height = 0.15
  ) +
  geom_point(
    size = 2.5
  ) +
  labs(
    title = "Environmental effects on annual body condition",
    subtitle = "Posterior mean and 95% credible intervals",
    x = "Change in mean log condition per 1-SD increase",
    y = NULL
  ) +
  theme_classic()

p_annual_condition_forest


p_annual_condition_combined <-
  p_annual_condition_relationships /
  p_annual_condition_forest +
  plot_layout(heights = c(2, 1))

p_annual_condition_combined












# -----------------------------
# 6.5 Joint model comparison - just PDO
# -----------------------------

dat_joint_pdo <- annual_joint %>%
  select(
    year,
    mean_age,
    condition_z,
    pdo_m_prev3
  ) %>%
  drop_na() %>%
  mutate(
    pdo_z = as.numeric(scale(pdo_m_prev3))
  )

# Null
m0 <- brm(
  mean_age ~ s(year, k = 8),
  data = dat_joint_pdo,
  family = student(),
  backend = "cmdstanr"
)

# Environment
m1 <- brm(
  mean_age ~ s(year, k = 8) + pdo_z,
  data = dat_joint_pdo,
  family = student(),
  backend = "cmdstanr"
)

# Condition
m2 <- brm(
  mean_age ~ s(year, k = 8) + condition_z,
  data = dat_joint_pdo,
  family = student(),
  backend = "cmdstanr"
)

# Combined
m3 <- brm(
  mean_age ~ s(year, k = 8) +
    pdo_z + condition_z,
  data = dat_joint_pdo,
  family = student(),
  backend = "cmdstanr"
)

# LOO comparison
loo_comparison <- loo_compare(
  loo(m0),
  loo(m1),
  loo(m2),
  loo(m3)
)

loo_comparison

# Examine posterior effects
fixef(m3)


### STEPPING THROUGH ALL VARIABLES
# -----------------------------
# 6.6.1 Environmental variables
# -----------------------------

library(dplyr)
library(purrr)
library(tidyr)
library(tibble)
library(brms)
library(posterior)
library(ggplot2)

annual_env_vars <- c(
  # "pdo_m_prev1",
  # "npgo_m_prev1",
  # "oni_m_prev1",
  # "npi_m_prev1",
  # "mei_m_prev1",
  #
  # "pdo_m_prev3",
  # "npgo_m_prev3",
  # "oni_m_prev3",
  # "npi_m_prev3",
  # "mei_m_prev3",

  "pdo_m_prev6",
  "npgo_m_prev6",
  "oni_m_prev6",
  "npi_m_prev6",
  "mei_m_prev6"
)

# Check availability
setdiff(annual_env_vars, names(annual_joint))




# -----------------------------
# 6.6.2 Model fitting function
# -----------------------------

fit_annual_joint_models <- function(
    variable,
    data = annual_joint
) {

  # Prepare data
  dat <- data %>%
    select(
      year,
      mean_age,
      condition_z,
      all_of(variable)
    ) %>%
    drop_na() %>%
    mutate(
      env_z = as.numeric(
        scale(.data[[variable]])
      )
    )

  if (nrow(dat) < 15) {
    stop("Too few observations for: ", variable)
  }

  # Define formulas
  formulas <- list(
    null = bf(
      mean_age ~ s(year, k = 8)
    ),

    environment = bf(
      mean_age ~ s(year, k = 8) + env_z
    ),

    condition = bf(
      mean_age ~ s(year, k = 8) + condition_z
    ),

    combined = bf(
      mean_age ~ s(year, k = 8) +
        env_z + condition_z
    )
  )

  # Fit models
  models <- map(
    formulas,
    function(form) {

      brm(
        formula = form,
        data = dat,
        family = student(),

        backend = "cmdstanr",

        chains = 4,
        cores = 4,
        iter = 3000,
        warmup = 1500,

        seed = 22,
        refresh = 0
      )
    }
  )

  list(
    variable = variable,
    data = dat,
    models = models
  )
}


# -----------------------------
# 6.6.3 Fit all environmental models
# -----------------------------

annual_joint_models <- annual_env_vars %>%
  set_names() %>%
  map(
    ~ fit_annual_joint_models(
      variable = .x,
      data = annual_joint
    )
  )

# Save models to avoid refitting
saveRDS(
  annual_joint_models,
  "annual_joint_models.rds"
)




# -----------------------------
# 6.6.4 Extract posterior effects
# -----------------------------

annual_joint_posterior <- imap_dfr(
  annual_joint_models,
  function(result, variable) {

    model <- result$models$combined

    draws <- as_draws_df(model)

    tibble(
      covariate = variable,

      env_beta = draws$b_env_z,

      condition_beta = draws$b_condition_z
    )
  }
)

# Summarize posterior distributions
annual_joint_summary <- annual_joint_posterior %>%
  group_by(covariate) %>%
  summarise(
    env_mean = mean(env_beta),
    env_lower = quantile(env_beta, 0.025),
    env_upper = quantile(env_beta, 0.975),
    prob_env_positive = mean(env_beta > 0),

    condition_mean = mean(condition_beta),
    condition_lower = quantile(condition_beta, 0.025),
    condition_upper = quantile(condition_beta, 0.975),
    prob_condition_positive = mean(condition_beta > 0),

    .groups = "drop"
  )

annual_joint_summary


# -----------------------------
# 6.6.5 Posterior coefficient plot
# -----------------------------

annual_joint_plot_data <- annual_joint_posterior %>%
  pivot_longer(
    cols = c(env_beta, condition_beta),
    names_to = "predictor",
    values_to = "beta"
  ) %>%
  mutate(
    predictor = recode(
      predictor,
      env_beta = "Environment",
      condition_beta = "Body condition"
    ),

    index = toupper(
      sub("_.*", "", covariate)
    ),

    exposure_window = case_when(
      grepl("prev1$", covariate) ~ "Previous month",
      grepl("prev3$", covariate) ~ "Previous 3 months",
      grepl("prev6$", covariate) ~ "Previous 6 months"
    )
  )

ggplot(
  annual_joint_plot_data,
  aes(
    x = beta,
    y = exposure_window
  )
) +
  ggdist::stat_halfeye(
    .width = c(0.50, 0.80, 0.95),
    point_interval = "median_qi"
  ) +
  geom_vline(
    xintercept = 0,
    linetype = "dashed"
  ) +
  facet_grid(
    predictor ~ index,
    scales = "free_x"
  ) +
  labs(
    title = "Environmental and condition associations with annual mean age",
    subtitle = "Posterior coefficients from combined models",
    x = "Change in annual mean age (years)",
    y = "Environmental exposure window"
  ) +
  theme_classic()



# -----------------------------
# 6.6.6 Calculate LOO
# -----------------------------

annual_joint_loo <- imap(
  annual_joint_models,
  function(result, variable) {

    loo_results <- map(
      result$models,
      ~ loo(
        .x,
        moment_match = FALSE
      )
    )

    comparison <- do.call(
      loo_compare,
      loo_results
    )

    list(
      variable = variable,
      loo = loo_results,
      comparison = comparison
    )
  }
)


condition_added_value <- imap_dfr(
  annual_joint_loo,
  function(result, variable) {

    comparison <- loo_compare(
      result$loo$combined,
      result$loo$environment
    )

    comparison <- as.data.frame(comparison)

    # Difference in ELPD: combined minus environment
    elpd_combined <- result$loo$combined$estimates[
      "elpd_loo", "Estimate"
    ]

    elpd_environment <- result$loo$environment$estimates[
      "elpd_loo", "Estimate"
    ]

    tibble(
      covariate = variable,
      elpd_gain_condition =
        elpd_combined - elpd_environment,

      se_diff = comparison[
        "environment",
        "se_diff"
      ]
    )
  }
)

condition_added_value

# -----------------------------
# 6.6.7 Summarize model comparison
# -----------------------------
annual_joint_loo_summary <- imap_dfr(
  annual_joint_loo,
  function(result, variable) {

    comparison <- as.data.frame(
      result$comparison
    )

    comparison %>%
      rownames_to_column("model") %>%
      mutate(
        covariate = variable
      ) %>%
      select(
        covariate,
        model,
        elpd_diff,
        se_diff
      )
  }
)

annual_joint_loo_summary


# plot LOO comparisons
ggplot(
  annual_joint_loo_summary,
  aes(
    x = elpd_diff,
    y = reorder(covariate, elpd_diff),
    colour = model
  )
) +
  geom_point(
    position = position_dodge(width = 0.6),
    size = 2
  ) +
  geom_errorbarh(
    aes(
      xmin = elpd_diff - se_diff,
      xmax = elpd_diff + se_diff
    ),
    position = position_dodge(width = 0.6),
    height = 0
  ) +
  geom_vline(
    xintercept = 0,
    linetype = "dashed"
  ) +
  labs(
    title = "Comparison of annual mean age models",
    subtitle = "ELPD differences relative to the best model for each index",
    x = "ELPD difference",
    y = "Environmental variable",
    colour = "Model"
  ) +
  theme_classic()


# -----------------------------
# 6.6.8 Condition add value to analysis
# -----------------------------

condition_added_value <- imap_dfr(
  annual_joint_loo,
  function(result, variable) {

    comparison <- loo_compare(
      result$loo$combined,
      result$loo$environment
    )

    comparison <- as.data.frame(comparison)

    # Difference in ELPD: combined minus environment
    elpd_combined <- result$loo$combined$estimates[
      "elpd_loo", "Estimate"
    ]

    elpd_environment <- result$loo$environment$estimates[
      "elpd_loo", "Estimate"
    ]

    tibble(
      covariate = variable,
      elpd_gain_condition =
        elpd_combined - elpd_environment,

      se_diff = comparison[
        "environment",
        "se_diff"
      ]
    )
  }
)

condition_added_value














