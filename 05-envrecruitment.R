##### -------------------------------------------------------------------- #####
##### -------------------------------------------------------------------- #####
##### -------------------------------------------------------------------- #####
##### -------------------------------------------------------------------- #####
##### -------------------------------------------------------------------- #####


##### -------------------------------------------------------------------- #####
# Hake recruitment and seasonal environmental indices
##### -------------------------------------------------------------------- #####

# Run section 3.1 of 03-environmentalanalysis.R first to create env_long.
# Set data_source to identify the input dataset; this does not filter the data.
# Recruitment is the pacea age-0 cohort estimate, in billions of fish.

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

recruitment_env_output_dir <- file.path("figures", "hake_recruitment", data_source)
if (!dir.exists(recruitment_env_output_dir) &&
  !dir.create(recruitment_env_output_dir, recursive = TRUE)) {
  stop("Could not create output directory: ", recruitment_env_output_dir)
}

# -----------------------------
# 5.1 Recruitment and seasonal environmental data
# -----------------------------

if (!exists("env_long")) {
  stop("Run section 3.1 of 03-environmentalanalysis.R to create env_long.")
}

indices <- c("pdo", "npgo", "oni", "npi", "mei")
season_levels <- c("Jan-Mar", "Apr-Jun", "Jul-Sep", "Oct-Dec")

required_env_columns <- c("year", "month", "covariate", "value")
missing_env_columns <- setdiff(required_env_columns, names(env_long))
if (length(missing_env_columns) > 0L) {
  stop("Missing columns in env_long: ", paste(missing_env_columns, collapse = ", "))
}

env_monthly <- env_long %>%
  ungroup() %>%
  filter(covariate %in% indices) %>%
  select(all_of(required_env_columns))

missing_indices <- setdiff(indices, unique(env_monthly$covariate))
if (length(missing_indices) > 0L) {
  stop("Missing environmental indices in env_long: ", paste(missing_indices, collapse = ", "))
}
if (anyNA(env_monthly$year) || anyNA(env_monthly$month) ||
  any(!env_monthly$month %in% 1:12)) {
  stop("Environmental records must have a year and a month from 1 to 12.")
}
if (anyDuplicated(env_monthly[c("year", "month", "covariate")])) {
  stop("env_long must contain at most one value per index, year, and month.")
}

# A seasonal mean is available only when all three monthly values are finite.
env_seasonal_recruitment <- env_monthly %>%
  mutate(
    season = factor(
      season_levels[(month - 1L) %/% 3L + 1L],
      levels = season_levels
    )
  ) %>%
  group_by(year, covariate, season) %>%
  summarise(
    index_value = if (n() == 3L && all(is.finite(value))) mean(value) else NA_real_,
    .groups = "drop"
  )

# Keep the draft's cutoff: recent estimates may be poorly informed by data.
recruitment <- as_tibble(pacea::hake_recruitment) %>%
  transmute(year = as.integer(year), recruitment = median) %>%
  filter(!is.na(year), year < 2023, is.finite(recruitment))

if (anyDuplicated(recruitment$year)) {
  stop("Recruitment must have one estimate per year.")
}
if (nrow(recruitment) == 0L) {
  stop("No recruitment estimates are available before 2023.")
}
if (any(recruitment$recruitment <= 0)) {
  stop(
    "Recruitment must be positive for the natural log; check years: ",
    paste(recruitment$year[recruitment$recruitment <= 0], collapse = ", ")
  )
}

recruitment <- recruitment %>%
  mutate(log_recruitment = log(recruitment))

# -----------------------------
# 5.2 Match previous environmental years to recruitment cohorts
# -----------------------------

# For cohort year Y, lag 1 uses seasons in Y-1, lag 2 in Y-2, and lag 3 in Y-3.
# No season from the recruitment year is included.
recruitment_env_all <- crossing(
  recruitment,
  covariate = indices,
  season = factor(season_levels, levels = season_levels),
  lag = 1:3
) %>%
  mutate(index_year = year - lag) %>%
  left_join(
    env_seasonal_recruitment,
    by = c("index_year" = "year", "covariate", "season")
  )

message(
  "Environmental values unavailable for ",
  sum(!is.finite(recruitment_env_all$index_value)),
  " recruitment/index/season/lag combinations."
)

recruitment_env <- recruitment_env_all %>%
  filter(is.finite(index_value))

write.csv(
  recruitment_env,
  file.path(recruitment_env_output_dir, "analysis_data.csv"),
  row.names = FALSE
)

# -----------------------------
# 5.3 Linear models and plots
# -----------------------------

# Use the same minimum sample size as the draft and report excluded groups.
recruitment_env_group_checks <- recruitment_env %>%
  group_by(covariate, season, lag) %>%
  summarise(
    n = n(),
    index_sd = sd(index_value),
    log_recruitment_sd = sd(log_recruitment),
    .groups = "drop"
  ) %>%
  mutate(usable = n >= 10 & is.finite(index_sd) & index_sd > 0 &
    is.finite(log_recruitment_sd) & log_recruitment_sd > 0)

excluded_groups <- recruitment_env_group_checks %>% filter(!usable)
if (nrow(excluded_groups) > 0L) {
  warning("Excluding ", nrow(excluded_groups), " groups with insufficient or constant data.")
  print(excluded_groups)
}

recruitment_env_groups <- recruitment_env %>%
  inner_join(
    recruitment_env_group_checks %>% filter(usable) %>% select(covariate, season, lag),
    by = c("covariate", "season", "lag")
  ) %>%
  group_by(covariate, season, lag) %>%
  group_split()

if (length(recruitment_env_groups) == 0L) {
  stop("No index-season-lag combinations have enough usable data.")
}

recruitment_env_linear_results <- map_dfr(recruitment_env_groups, \(dat) {
  fit <- lm(log_recruitment ~ index_value, data = dat)
  broom::tidy(fit, conf.int = TRUE) %>%
    mutate(
      covariate = dat$covariate[1],
      season = dat$season[1],
      lag = dat$lag[1],
      n = nrow(dat),
      r_squared = summary(fit)$r.squared
    )
}) %>%
  group_by(term) %>%
  mutate(p_adjust_bh = p.adjust(p.value, method = "BH")) %>%
  ungroup()

write.csv(
  recruitment_env_linear_results,
  file.path(recruitment_env_output_dir, "linear_model_results.csv"),
  row.names = FALSE
)
print(recruitment_env_linear_results, n = Inf)

recruitment_env_model_data <- bind_rows(recruitment_env_groups)
recruitment_env_linear_plots <- split(
  recruitment_env_model_data,
  recruitment_env_model_data$covariate
) %>%
  imap(\(dat, index) {
    ggplot(dat, aes(index_value, log_recruitment)) +
      geom_point(alpha = 0.7) +
      geom_smooth(method = "lm", formula = y ~ x, se = TRUE) +
      facet_grid(season ~ lag, labeller = label_both) +
      labs(
        title = paste(toupper(index), "and hake recruitment"),
        x = "Seasonal mean environmental index",
        y = "Natural log of age-0 recruitment (billions of fish)"
      ) +
      theme_bw()
  })

iwalk(recruitment_env_linear_plots, \(plot, index) {
  ggsave(
    file.path(recruitment_env_output_dir, paste0("linear_", index, ".png")),
    plot = plot,
    width = 10,
    height = 10
  )
})

# -----------------------------
# 5.4 Bayesian models and posterior summaries
# -----------------------------

# The separate regressions standardize the index and log recruitment; they use
# recruitment medians and do not propagate stock assessment uncertainty.
recruitment_env_bayesian_results <- map_dfr(recruitment_env_groups, \(dat) {
  dat <- dat %>%
    mutate(
      x = as.numeric(scale(index_value)),
      y = as.numeric(scale(log_recruitment))
    )

  label <- paste(dat$covariate[1], gsub("-", "_", dat$season[1]), dat$lag[1], sep = "_")
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

  saveRDS(fit, file.path(recruitment_env_output_dir, paste0("bayesian_", label, ".rds")))
  # png(file.path(recruitment_env_output_dir, paste0("bayesian_diagnostics_", label, ".png")),
  #   width = 1200, height = 800, res = 120
  # )
  # print(plot(fit, variable = "b_x"))
  # dev.off()

  as.data.frame(posterior_summary(fit, variable = c("b_Intercept", "b_x", "sigma"))) %>%
    tibble::rownames_to_column("parameter") %>%
    mutate(
      covariate = dat$covariate[1],
      season = dat$season[1],
      lag = dat$lag[1],
      n = nrow(dat)
    )
})

write.csv(
  recruitment_env_bayesian_results,
  file.path(recruitment_env_output_dir, "bayesian_model_results.csv"),
  row.names = FALSE
)
print(recruitment_env_bayesian_results, n = Inf)

# -----------------------------
# 5.5 Posterior distributions of environmental coefficients
# -----------------------------

# The fitted models use a separately standardized index and log recruitment
# for each index-season-lag combination.
recruitment_env_posterior_data <- map(recruitment_env_groups, \(dat) {
  label <- paste(dat$covariate[1], gsub("-", "_", dat$season[1]), dat$lag[1], sep = "_")
  model_file <- file.path(recruitment_env_output_dir, paste0("bayesian_", label, ".rds"))
  if (!file.exists(model_file)) {
    stop("Missing Bayesian model: ", model_file, ". Run section 5.4 first.")
  }
  fit <- readRDS(model_file)
  if (!all(c("x", "y") %in% names(fit$data)) ||
    !isTRUE(all.equal(fit$data$x, as.numeric(scale(dat$index_value)))) ||
    !isTRUE(all.equal(fit$data$y, as.numeric(scale(dat$log_recruitment))))) {
    stop("Saved model data do not match the current data for ", label, ". Rerun section 5.4.")
  }

  coefficient_draws <- as.numeric(posterior::as_draws_matrix(fit, variable = "b_x")[, "b_x"])
  fitted_values <- fitted(fit, re_formula = NA, probs = c(0.025, 0.975)) %>%
    as_tibble()
  log_mean <- mean(dat$log_recruitment)
  log_sd <- sd(dat$log_recruitment)

  list(
    coefficients = tibble(
      covariate = dat$covariate[1],
      season = dat$season[1],
      lag = dat$lag[1],
      coefficient = coefficient_draws
    ),
    fitted = bind_cols(
      dat %>% select(year, index_year, covariate, season, lag, index_value, log_recruitment),
      fitted_values
    ) %>%
      mutate(
        Estimate = log_mean + log_sd * Estimate,
        Q2.5 = log_mean + log_sd * Q2.5,
        Q97.5 = log_mean + log_sd * Q97.5
      )
  )
})

recruitment_env_coefficient_draws <- map_dfr(recruitment_env_posterior_data, "coefficients")
recruitment_env_coefficient_summary <- recruitment_env_coefficient_draws %>%
  group_by(covariate, season, lag) %>%
  summarise(
    median = median(coefficient),
    lower_95 = quantile(coefficient, 0.025),
    upper_95 = quantile(coefficient, 0.975),
    .groups = "drop"
  )

recruitment_env_coefficient_plots <- split(
  recruitment_env_coefficient_draws,
  recruitment_env_coefficient_draws$covariate
) %>%
  imap(\(dat, index) {
    ggplot(dat, aes(x = coefficient)) +
      geom_vline(xintercept = 0, linetype = "dashed") +
      geom_density(fill = "steelblue", alpha = 0.4) +
      facet_grid(season ~ lag, labeller = label_both) +
      labs(
        title = paste(toupper(index), "posterior environmental coefficients"),
        subtitle = "Index and log recruitment standardized separately within each model",
        x = "Standardized coefficient",
        y = "Posterior density"
      ) +
      theme_bw()
  })

iwalk(recruitment_env_coefficient_plots, \(plot, index) {
  ggsave(
    file.path(recruitment_env_output_dir, paste0("posterior_coefficients_", index, ".png")),
    plot = plot,
    width = 10,
    height = 10
  )
})

# -----------------------------
# 5.6 Observed and posterior fitted recruitment
# -----------------------------

# Fitted means and intervals are transformed back from standardized values to
# the natural-log recruitment scale. Intervals describe expected values, not
# the wider posterior predictive intervals for future observations.
recruitment_env_posterior_fitted <- map_dfr(recruitment_env_posterior_data, "fitted")

recruitment_env_fitted_time_plots <- split(
  recruitment_env_posterior_fitted,
  recruitment_env_posterior_fitted$covariate
) %>%
  imap(\(dat, index) {
    ggplot(dat, aes(x = year)) +
      geom_ribbon(aes(ymin = Q2.5, ymax = Q97.5), alpha = 0.2) +
      geom_line(aes(y = Estimate)) +
      geom_point(aes(y = log_recruitment), alpha = 0.7) +
      facet_grid(season ~ lag, labeller = label_both) +
      labs(
        title = paste(toupper(index), "observed and fitted recruitment"),
        subtitle = "Posterior expected value and 95% credible interval",
        x = "Recruitment year",
        y = "Natural log of age-0 recruitment (billions of fish)"
      ) +
      theme_bw()
  })

recruitment_env_observed_fitted_plots <- split(
  recruitment_env_posterior_fitted,
  recruitment_env_posterior_fitted$covariate
) %>%
  imap(\(dat, index) {
    ggplot(dat, aes(x = Estimate, y = log_recruitment)) +
      geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
      geom_segment(aes(x = Q2.5, xend = Q97.5, yend = log_recruitment), alpha = 0.35) +
      geom_point(alpha = 0.7) +
      facet_grid(season ~ lag, labeller = label_both) +
      labs(
        title = paste(toupper(index), "observed versus fitted recruitment"),
        subtitle = "Horizontal lines show 95% credible intervals for fitted values",
        x = "Fitted natural-log recruitment",
        y = "Observed natural-log recruitment"
      ) +
      theme_bw()
  })

iwalk(recruitment_env_fitted_time_plots, \(plot, index) {
  ggsave(
    file.path(recruitment_env_output_dir, paste0("posterior_fitted_time_", index, ".png")),
    plot = plot,
    width = 10,
    height = 10
  )
})

iwalk(recruitment_env_observed_fitted_plots, \(plot, index) {
  ggsave(
    file.path(recruitment_env_output_dir, paste0("posterior_observed_vs_fitted_", index, ".png")),
    plot = plot,
    width = 10,
    height = 10
  )
})
