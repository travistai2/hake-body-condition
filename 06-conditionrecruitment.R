##### -------------------------------------------------------------------- #####
# Seasonal hake body condition and following-year recruitment
##### -------------------------------------------------------------------- #####

# Run section 2.5 of 02-lengthweightanalysis.R first to create cond_monthly.
# The condition index depends on the data source and length-weight model selected there.
# Set data_source to identify the input dataset; this does not filter the data.

library(dplyr)
library(tidyr)
library(ggplot2)
library(purrr)
library(broom)
library(brms)

if (!exists("data_source") || !is.character(data_source) ||
  length(data_source) != 1L || is.na(data_source) || !nzchar(data_source) ||
  data_source %in% c(".", "..") || grepl("[/\\\\]", data_source)) {
  stop("Set data_source to a single non-empty dataset name without path separators.")
}

condition_recruitment_output_dir <- file.path("figures", "hake_condition_recruitment", data_source)
if (!dir.exists(condition_recruitment_output_dir) &&
  !dir.create(condition_recruitment_output_dir, recursive = TRUE)) {
  stop("Could not create output directory: ", condition_recruitment_output_dir)
}

# -----------------------------
# 6.1 Seasonal body condition and recruitment
# -----------------------------

if (!exists("cond_monthly")) {
  stop("Run section 2.5 of 02-lengthweightanalysis.R to create cond_monthly.")
}

required_condition_columns <- c("year", "month", "mean_log_Kn", "n_fish")
missing_condition_columns <- setdiff(required_condition_columns, names(cond_monthly))
if (length(missing_condition_columns) > 0L) {
  stop("Missing columns in cond_monthly: ", paste(missing_condition_columns, collapse = ", "))
}

condition_monthly <- cond_monthly %>%
  ungroup() %>%
  select(all_of(required_condition_columns))

if (anyNA(condition_monthly$year) || anyNA(condition_monthly$month) ||
  any(!condition_monthly$month %in% 1:12)) {
  stop("Condition records must have a year and a month from 1 to 12.")
}
if (anyDuplicated(condition_monthly[c("year", "month")])) {
  stop("cond_monthly must contain at most one row per year and month.")
}
if (any(!is.finite(condition_monthly$n_fish) | condition_monthly$n_fish <= 0)) {
  stop("Each monthly condition estimate must have a positive fish count.")
}

# Match the seasonal condition summaries in section 3.21: equally weight
# available monthly mean log(Kn) values in April-June and August-October.
condition_seasonal <- condition_monthly %>%
  filter(month %in% c(4:6, 8:10)) %>%
  mutate(
    season = case_when(
      month %in% 4:6 ~ "Spring (Apr-Jun)",
      month %in% 8:10 ~ "Fall (Aug-Oct)"
    ),
    season = factor(season, levels = c("Spring (Apr-Jun)", "Fall (Aug-Oct)"))
  ) %>%
  group_by(year, season) %>%
  summarise(
    condition = if (any(is.finite(mean_log_Kn))) {
      mean(mean_log_Kn[is.finite(mean_log_Kn)])
    } else {
      NA_real_
    },
    n_months = sum(is.finite(mean_log_Kn)),
    n_fish = sum(n_fish[is.finite(mean_log_Kn)]),
    .groups = "drop"
  )

# As in script 05, omit recent recruitment estimates that may be poorly informed.
condition_recruitment <- as_tibble(pacea::hake_recruitment) %>%
  transmute(recruitment_year = as.integer(year), recruitment = median) %>%
  filter(!is.na(recruitment_year), recruitment_year < 2023, is.finite(recruitment))

if (nrow(condition_recruitment) == 0L || anyDuplicated(condition_recruitment$recruitment_year)) {
  stop("Recruitment must contain one estimate per year before 2023.")
}
if (any(condition_recruitment$recruitment <= 0)) {
  stop(
    "Recruitment must be positive for the natural log; check years: ",
    paste(condition_recruitment$recruitment_year[condition_recruitment$recruitment <= 0], collapse = ", ")
  )
}

condition_recruitment <- condition_recruitment %>%
  mutate(log_recruitment = log(recruitment))

# -----------------------------
# 6.2 Match condition to the following recruitment year
# -----------------------------

condition_recruitment_all <- crossing(
  condition_recruitment,
  season = factor(c("Spring (Apr-Jun)", "Fall (Aug-Oct)"),
    levels = c("Spring (Apr-Jun)", "Fall (Aug-Oct)")
  )
) %>%
  mutate(condition_year = recruitment_year - 1L) %>%
  left_join(condition_seasonal, by = c("condition_year" = "year", "season"))

message(
  "Seasonal condition unavailable for ",
  sum(!is.finite(condition_recruitment_all$condition)),
  " recruitment-year/season combinations."
)

condition_recruitment_data <- condition_recruitment_all %>%
  filter(is.finite(condition))

write.csv(
  condition_recruitment_data,
  file.path(condition_recruitment_output_dir, "analysis_data.csv"),
  row.names = FALSE
)

# -----------------------------
# 6.3 Linear models and observed relationships
# -----------------------------

condition_recruitment_group_checks <- condition_recruitment_data %>%
  group_by(season) %>%
  summarise(
    n = n(),
    condition_sd = sd(condition),
    log_recruitment_sd = sd(log_recruitment),
    .groups = "drop"
  ) %>%
  mutate(usable = n >= 10 & is.finite(condition_sd) & condition_sd > 0 &
    is.finite(log_recruitment_sd) & log_recruitment_sd > 0)

excluded_seasons <- condition_recruitment_group_checks %>% filter(!usable)
missing_seasons <- setdiff(levels(condition_recruitment_all$season), as.character(condition_recruitment_group_checks$season))
if (nrow(excluded_seasons) > 0L || length(missing_seasons) > 0L) {
  warning("Some seasons have insufficient, constant, or missing data for modeling.")
  if (nrow(excluded_seasons) > 0L) print(excluded_seasons)
  if (length(missing_seasons) > 0L) message("Seasons without matched data: ", paste(missing_seasons, collapse = ", "))
}

condition_recruitment_groups <- condition_recruitment_data %>%
  inner_join(
    condition_recruitment_group_checks %>% filter(usable) %>% select(season),
    by = "season"
  ) %>%
  group_by(season) %>%
  group_split()

if (length(condition_recruitment_groups) == 0L) {
  stop("No seasons have at least 10 usable years with variable condition and recruitment.")
}

condition_recruitment_linear_results <- map_dfr(condition_recruitment_groups, \(dat) {
  fit <- lm(log_recruitment ~ condition, data = dat)
  broom::tidy(fit, conf.int = TRUE) %>%
    mutate(season = dat$season[1], n = nrow(dat), r_squared = summary(fit)$r.squared)
}) %>%
  group_by(term) %>%
  mutate(p_adjust_bh = p.adjust(p.value, method = "BH")) %>%
  ungroup()

write.csv(
  condition_recruitment_linear_results,
  file.path(condition_recruitment_output_dir, "linear_model_results.csv"),
  row.names = FALSE
)
print(condition_recruitment_linear_results, n = Inf)

condition_recruitment_model_data <- bind_rows(condition_recruitment_groups)
condition_recruitment_linear_plot <- ggplot(
  condition_recruitment_model_data,
  aes(x = condition, y = log_recruitment)
) +
  geom_point(alpha = 0.7) +
  geom_smooth(method = "lm", formula = y ~ x, se = TRUE) +
  facet_wrap(~season) +
  labs(
    title = "Body condition and following-year hake recruitment",
    x = "Seasonal mean log(Kn)",
    y = "Natural log of age-0 recruitment (billions of fish)"
  ) +
  theme_bw()

ggsave(
  file.path(condition_recruitment_output_dir, "linear_condition_recruitment.png"),
  plot = condition_recruitment_linear_plot,
  width = 10,
  height = 5
)

# -----------------------------
# 6.4 Bayesian models and posterior summaries
# -----------------------------

# As in script 05, each regression standardizes both variables separately.
# The analysis uses point estimates of condition and recruitment and does not
# propagate uncertainty in the length-weight or stock assessment models.
condition_recruitment_bayesian_results <- map_dfr(condition_recruitment_groups, \(dat) {
  dat <- dat %>%
    mutate(x = as.numeric(scale(condition)), y = as.numeric(scale(log_recruitment)))
  label <- if (as.character(dat$season[1]) == "Spring (Apr-Jun)") "spring" else "fall"
  message("Fitting Bayesian model: ", label)

  fit <- brm(
    y ~ x,
    data = dat,
    family = gaussian(),
    prior = c(
      prior(normal(0, 1), class = "b"),
      prior(normal(0, 1), class = "Intercept"),
      prior(exponential(1), class = "sigma")
    ),
    chains = 4,
    cores = 4,
    iter = 2000,
    backend = "cmdstanr",
    seed = 1234,
    refresh = 0
  )

  saveRDS(fit, file.path(condition_recruitment_output_dir, paste0("bayesian_", label, ".rds")))
  as.data.frame(posterior_summary(fit, variable = c("b_Intercept", "b_x", "sigma"))) %>%
    tibble::rownames_to_column("parameter") %>%
    mutate(season = dat$season[1], n = nrow(dat))
})

write.csv(
  condition_recruitment_bayesian_results,
  file.path(condition_recruitment_output_dir, "bayesian_model_results.csv"),
  row.names = FALSE
)
print(condition_recruitment_bayesian_results, n = Inf)

# -----------------------------
# 6.5 Posterior distributions of condition coefficients
# -----------------------------

condition_recruitment_posterior_data <- map(condition_recruitment_groups, \(dat) {
  label <- if (as.character(dat$season[1]) == "Spring (Apr-Jun)") "spring" else "fall"
  model_file <- file.path(condition_recruitment_output_dir, paste0("bayesian_", label, ".rds"))
  if (!file.exists(model_file)) {
    stop("Missing Bayesian model: ", model_file, ". Run section 6.4 first.")
  }
  fit <- readRDS(model_file)
  if (!all(c("x", "y") %in% names(fit$data)) ||
    !isTRUE(all.equal(fit$data$x, as.numeric(scale(dat$condition)))) ||
    !isTRUE(all.equal(fit$data$y, as.numeric(scale(dat$log_recruitment))))) {
    stop("Saved model data do not match the current data for ", label, ". Rerun section 6.4.")
  }

  coefficient_draws <- as.numeric(posterior::as_draws_matrix(fit, variable = "b_x")[, "b_x"])
  fitted_values <- fitted(fit, re_formula = NA, probs = c(0.025, 0.975)) %>%
    as_tibble()
  log_mean <- mean(dat$log_recruitment)
  log_sd <- sd(dat$log_recruitment)

  list(
    model = fit,
    coefficients = tibble(season = dat$season[1], coefficient = coefficient_draws),
    fitted = bind_cols(dat, fitted_values) %>%
      mutate(
        Estimate = log_mean + log_sd * Estimate,
        Q2.5 = log_mean + log_sd * Q2.5,
        Q97.5 = log_mean + log_sd * Q97.5
      )
  )
})

condition_recruitment_coefficient_draws <- map_dfr(condition_recruitment_posterior_data, "coefficients")
condition_recruitment_coefficient_summary <- condition_recruitment_coefficient_draws %>%
  group_by(season) %>%
  summarise(
    median = median(coefficient),
    lower_95 = quantile(coefficient, 0.025),
    upper_95 = quantile(coefficient, 0.975),
    .groups = "drop"
  )

condition_recruitment_coefficient_plot <- ggplot(
  condition_recruitment_coefficient_draws,
  aes(x = coefficient)
) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_density(fill = "steelblue", alpha = 0.4) +
  facet_wrap(~season) +
  labs(
    title = "Posterior body-condition coefficients",
    subtitle = "Condition and log recruitment standardized separately within each model",
    x = "Standardized coefficient",
    y = "Posterior density"
  ) +
  theme_bw()

ggsave(
  file.path(condition_recruitment_output_dir, "posterior_condition_coefficients.png"),
  plot = condition_recruitment_coefficient_plot,
  width = 10,
  height = 5
)

# -----------------------------
# 6.6 Observed and posterior fitted recruitment
# -----------------------------

# Intervals describe expected values on the natural-log recruitment scale,
# not posterior predictive intervals for future observations.
condition_recruitment_posterior_fitted <- map_dfr(condition_recruitment_posterior_data, "fitted")

condition_recruitment_fitted_time_plot <- ggplot(
  condition_recruitment_posterior_fitted,
  aes(x = recruitment_year)
) +
  geom_ribbon(aes(ymin = Q2.5, ymax = Q97.5), alpha = 0.2) +
  geom_line(aes(y = Estimate)) +
  geom_point(aes(y = log_recruitment), alpha = 0.7) +
  facet_wrap(~season) +
  labs(
    title = "Observed and fitted recruitment following seasonal condition",
    subtitle = "Posterior expected value and 95% credible interval",
    x = "Recruitment year",
    y = "Natural log of age-0 recruitment (billions of fish)"
  ) +
  theme_bw()

condition_recruitment_observed_fitted_plot <- ggplot(
  condition_recruitment_posterior_fitted,
  aes(x = Estimate, y = log_recruitment)
) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  geom_segment(aes(x = Q2.5, xend = Q97.5, yend = log_recruitment), alpha = 0.35) +
  geom_point(alpha = 0.7) +
  facet_wrap(~season) +
  labs(
    title = "Observed versus fitted recruitment",
    subtitle = "Horizontal lines show 95% credible intervals for fitted values",
    x = "Fitted natural-log recruitment",
    y = "Observed natural-log recruitment"
  ) +
  theme_bw()

ggsave(
  file.path(condition_recruitment_output_dir, "posterior_fitted_time.png"),
  plot = condition_recruitment_fitted_time_plot,
  width = 10,
  height = 5
)
ggsave(
  file.path(condition_recruitment_output_dir, "posterior_observed_vs_fitted.png"),
  plot = condition_recruitment_observed_fitted_plot,
  width = 10,
  height = 5
)

# -----------------------------
# 6.7 Posterior recruitment trends and observed data
# -----------------------------

# Predictions describe the expected log recruitment, with a 95% credible
# interval for the trend rather than a predictive interval for observations.
condition_recruitment_posterior_trends <- map2_dfr(
  condition_recruitment_posterior_data,
  condition_recruitment_groups,
  \(posterior_data, dat) {
    condition_values <- seq(min(dat$condition), max(dat$condition), length.out = 100)
    prediction_data <- tibble(
      x = (condition_values - mean(dat$condition)) / sd(dat$condition)
    )
    trend <- fitted(
      posterior_data$model,
      newdata = prediction_data,
      re_formula = NA,
      probs = c(0.025, 0.975)
    ) %>%
      as_tibble()

    bind_cols(
      tibble(season = dat$season[1], condition = condition_values),
      trend
    ) %>%
      mutate(
        Estimate = mean(dat$log_recruitment) + sd(dat$log_recruitment) * Estimate,
        Q2.5 = mean(dat$log_recruitment) + sd(dat$log_recruitment) * Q2.5,
        Q97.5 = mean(dat$log_recruitment) + sd(dat$log_recruitment) * Q97.5
      )
  }
)

condition_recruitment_posterior_trend_plot <- ggplot(
  condition_recruitment_posterior_trends,
  aes(x = condition, y = Estimate)
) +
  geom_ribbon(aes(ymin = Q2.5, ymax = Q97.5), alpha = 0.2) +
  geom_line(linewidth = 0.8) +
  geom_point(
    data = condition_recruitment_model_data,
    aes(x = condition, y = log_recruitment),
    inherit.aes = FALSE,
    alpha = 0.7
  ) +
  facet_wrap(~season, scales = "free_x") +
  labs(
    title = "Body condition and following-year hake recruitment",
    subtitle = "Bayesian posterior mean trend and 95% credible interval",
    x = "Seasonal mean log(Kn)",
    y = "Natural log of age-0 recruitment (billions of fish)"
  ) +
  theme_bw()

ggsave(
  file.path(condition_recruitment_output_dir, "posterior_condition_recruitment_trend.png"),
  plot = condition_recruitment_posterior_trend_plot,
  width = 10,
  height = 5
)
