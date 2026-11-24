##### -------------------------------------------------------------------- #####
##### -------------------------------------------------------------------- #####
##### -------------------------------------------------------------------- #####
##### -------------------------------------------------------------------- #####
##### -------------------------------------------------------------------- #####

# Annual mean age analysis with environmental variables

# -----------------------------
# 4.1 Summarize hake data into annual mean age and combine with environmental variables
# -----------------------------

# select data source for age analysis
data_source <- "LWAdat"
# data_source <- "usdat"
# data_source <- "LWAdat_can"

if (data_source == "LWAdat") {
  selected_data <- lwa
} else if (data_source == "usdat") {
  selected_data <- lwa
} else if (data_source == "LWAdat_can") {
  selected_data <- lwa %>% filter(Source %in% c("CAN_shoreside", "CAN_polish", "CAN_JV", "CAN_acoustic", "Acoustic Canada"))
}

# summarise annual mean age
annual_age <- selected_data %>%
  group_by(Year) %>%
  summarise(
    mean_age = mean(Age_yrs_numeric, na.rm = TRUE),
    median_age = median(Age_yrs_numeric, na.rm = TRUE),
    age_sd = sd(Age_yrs_numeric, na.rm = TRUE),
    age_q25 = quantile(Age_yrs_numeric, probs = 0.25, na.rm = TRUE),
    age_q75 = quantile(Age_yrs_numeric, probs = 0.75, na.rm = TRUE),
    age_iqr = IQR(Age_yrs_numeric, na.rm = TRUE),
    n_fish = sum(!is.na(Age_yrs_numeric)),
    n_sources = n_distinct(Source[!is.na(Source)]),
    .groups = "drop"
  ) %>%
  ungroup() %>%
  mutate(Year = as.numeric(as.character(Year)))

# summarise seasonal environmental variables
env_seasonal_summarise <- env_long %>%
  group_by(year, covariate) %>%
  summarise(
    jantomar = mean(value[month %in% 1:3], na.rm = TRUE),
    aprtojun = mean(value[month %in% 4:6], na.rm = TRUE),
    jultosep = mean(value[month %in% 7:9], na.rm = TRUE),
    octtodec = mean(value[month %in% 10:12], na.rm = TRUE),
    .groups = "drop"
  ) %>%
  pivot_longer(
    cols = c(jantomar, aprtojun, jultosep, octtodec),
    names_to = "season",
    values_to = "value"
  ) %>%
  mutate(covariate = paste(covariate, season, sep = "_")) %>%
  mutate(
    yearlag1 = year + 1,
    yearlag2 = year + 2,
    yearlag3 = year + 3
  ) %>%
  select(-season) %>%
  pivot_wider(names_from = covariate, values_from = value)

# join environmental data to current year and 1 year lag to annual mean age data
annual_age_env_current <- annual_age %>%
  left_join(env_seasonal_summarise %>% select(-yearlag1, -yearlag2, -yearlag3),
    by = c("Year" = "year")
  ) %>%
  mutate(scenario = "currentyear") %>%
  arrange(Year)

annual_age_env_lag1 <- annual_age %>%
  left_join(env_seasonal_summarise %>% select(-year, -yearlag2, -yearlag3),
    by = c("Year" = "yearlag1")
  ) %>%
  mutate(scenario = "lag1year") %>%
  arrange(Year)

annual_age_env_lag2 <- annual_age %>%
  left_join(env_seasonal_summarise %>% select(-year, -yearlag1, -yearlag3),
    by = c("Year" = "yearlag2")
  ) %>%
  mutate(scenario = "lag2year") %>%
  arrange(Year)

annual_age_env_lag3 <- annual_age %>%
  left_join(env_seasonal_summarise %>% select(-year, -yearlag1, -yearlag2),
    by = c("Year" = "yearlag3")
  ) %>%
  mutate(scenario = "lag3year") %>%
  arrange(Year)

# join data together
annual_age_env <- rbind(
  annual_age_env_current,
  annual_age_env_lag1,
  annual_age_env_lag2,
  annual_age_env_lag3
)

names(annual_age_env)

# -----------------------------
# 4.2 Plot mean sampled age through time and env variables over time
# -----------------------------

# plot env variables over time
p_env_variables_time <- annual_age_env %>%
  filter(scenario == "currentyear") %>%
  pivot_longer(
    cols = starts_with(c("pdo", "npgo", "oni", "bi", "npi", "mei")),
    names_to = "covariate",
    values_to = "value"
  ) %>%
  mutate(
    index = toupper(sub("_.*", "", covariate)),
    window = sub(".*_", "", covariate)
  ) %>%
  filter(window == "aprtojun") %>%
  ggplot(aes(x = Year, y = value)) +
  geom_point() +
  geom_smooth(method = "lm", se = TRUE, color = "blue") +
  facet_wrap(~index, scales = "free") +
  labs(
    title = "Environmental variables over time (Apr-Jun)",
    x = NULL,
    y = "Standardized value"
  ) +
  theme_classic()

p_env_variables_time

# save plot
ggsave(p_env_variables_time, filename = "figures/env_variables_over_time_spring_plot.png", width = 8, height = 5, units = "in", dpi = 300, overwrite = FALSE)

# mean age over time
p_annual_mean_age <- annual_age_env %>%
  filter(scenario == "currentyear") %>%
  ggplot(aes(x = Year, y = mean_age)) +
  geom_line(linewidth = 0.5) +
  geom_point(aes(size = n_fish), alpha = 0.5) +
  geom_smooth(method = "lm", se = TRUE, color = "blue") +
  scale_size_continuous(trans = "sqrt") +
  labs(
    title = "Annual mean age of sampled fish",
    x = NULL,
    y = "Mean age (years)",
    size = "Fish sampled"
  ) +
  theme_classic()

p_annual_mean_age

# save plot
ggsave(p_annual_mean_age, filename = paste0("figures/", data_source, "/annual_mean_age_plot.png"), width = 8, height = 5, units = "in", dpi = 300)

# -----------------------------
# 4.3 Plot mean age correlations with environmental variables
# -----------------------------

# function to plot mean age correlations with env variables
plot_mean_age_env_correlations <- function(data, scenario_label) {
  data %>%
    filter(scenario == scenario_label) %>%
    pivot_longer(
      cols = starts_with(c("pdo", "npgo", "oni", "bi", "npi", "mei")),
      names_to = "covariate",
      values_to = "value"
    ) %>%
    mutate(
      index = toupper(sub("_.*", "", covariate)),
      window = sub(".*_", "", covariate),
      window = recode(
        window,
        jantomar = "Jan to Mar",
        aprtojun = "Apr to Jun",
        jultosep = "Jul to Sep",
        octtodec = "Oct to Dec"
      ),
      window = factor(
        window,
        levels = c("Jan to Mar", "Apr to Jun", "Jul to Sep", "Oct to Dec")
      )
    ) %>%
    ggplot(aes(x = value, y = mean_age)) +
    geom_point(alpha = 0.5) +
    geom_smooth(method = "lm", se = TRUE, color = "blue") +
    facet_wrap(index ~ window, scales = "free", ncol = 4) +
    labs(
      title = paste("Annual mean age correlations with environmental variables (", scenario_label, ")", sep = ""),
      x = "Environmental variable (standardized)",
      y = "Mean age (years)"
    ) +
    theme_classic()
}

# plot current year correlations
plot_mean_age_env_correlations(annual_age_env, "currentyear")
ggsave(last_plot(), filename = paste0("figures/", data_source, "/annual_mean_age_env_cor_current.png"), width = 8, height = 8, units = "in", dpi = 300)

# plot lag 1 year correlations
plot_mean_age_env_correlations(annual_age_env, "lag1year")
ggsave(last_plot(), filename = paste0("figures/", data_source, "/annual_mean_age_env_cor_lag1.png"), width = 8, height = 8, units = "in", dpi = 300)

# plot lag 2 year correlations
plot_mean_age_env_correlations(annual_age_env, "lag2year")
ggsave(last_plot(), filename = paste0("figures/", data_source, "/annual_mean_age_env_cor_lag2.png"), width = 8, height = 8, units = "in", dpi = 300)

# plot lag 3 year correlations
plot_mean_age_env_correlations(annual_age_env, "lag3year")
ggsave(last_plot(), filename = paste0("figures/", data_source, "/annual_mean_age_env_cor_lag3.png"), width = 8, height = 8, units = "in", dpi = 300)

# correlations with one season at a time with current and next years mean age

cor_table <- cbind(
  # current year
  annual_age_env %>%
    filter(scenario == "currentyear") %>%
    select(mean_age, contains(c("pdo", "npgo", "oni", "bi", "npi", "mei"))) %>%
    cor(use = "complete.obs") %>%
    as.data.frame() %>%
    select(mean_age) %>%
    rename(mean_age_cor_currentyear = mean_age),
  # previous year
  annual_age_env %>%
    filter(scenario == "lag1year") %>%
    select(mean_age, contains(c("pdo", "npgo", "oni", "bi", "npi", "mei"))) %>%
    cor(use = "complete.obs") %>%
    as.data.frame() %>%
    select(mean_age) %>%
    rename(mean_age_cor_lag1year = mean_age),
  annual_age_env %>%
    filter(scenario == "lag2year") %>%
    select(mean_age, contains(c("pdo", "npgo", "oni", "bi", "npi", "mei"))) %>%
    cor(use = "complete.obs") %>%
    as.data.frame() %>%
    select(mean_age) %>%
    rename(mean_age_cor_lag2year = mean_age),
  annual_age_env %>%
    filter(scenario == "lag3year") %>%
    select(mean_age, contains(c("pdo", "npgo", "oni", "bi", "npi", "mei"))) %>%
    cor(use = "complete.obs") %>%
    as.data.frame() %>%
    select(mean_age) %>%
    rename(mean_age_cor_lag3year = mean_age)
) %>%
  rownames_to_column(var = "env_variable") %>%
  mutate(
    index = toupper(sub("_.*", "", env_variable)),
    window = sub(".*_", "", env_variable),
    window = factor(window, levels = c("jantomar", "aprtojun", "jultosep", "octtodec"))
  )

ggplot(cor_table, aes(x = window)) +
  geom_point(aes(y = mean_age_cor_currentyear, color = "Current Year")) +
  geom_point(aes(y = mean_age_cor_lag1year, color = "Lag 1 Year")) +
  geom_point(aes(y = mean_age_cor_lag2year, color = "Lag 2 Year")) +
  geom_point(aes(y = mean_age_cor_lag3year, color = "Lag 3 Year")) +
  facet_wrap(~index, scales = "free") +
  labs(
    x = "Environmental Variable",
    y = "Correlation with Mean Age",
    color = "Scenario"
  )

ggsave(last_plot(), filename = paste0("figures/", data_source, "/annual_mean_age_env_cor_comparison.png"), width = 8, height = 8, units = "in", dpi = 300)

# For all data:
# lag 2 year seems to generally have stronger correlations with mean age than the current or lag 1year
# lag 2 year effects on mean age could mean environmental effects on body condition affecting survival and reproduction
# lag 1 year effects on mean age could indicate delayed environmental impacts on survival and reproduction
# current year effects on mean age could reflect immediate environmental influences on survival (but not necessarily reproduction)

# -----------------------------
# 4.4 Function to run individual models for each environmental covariate
# -----------------------------

# function to fit annual age-environment models for a given environmental variable
fit_annual_age_env_model <- function(
  data,
  env_variable,
  year_lag,
  output_directory = "model-outputs/annual_mean_age_env",
  seed = 22
) {
  model_data <- data %>%
    filter(
      if_all(
        all_of(c("Year", "mean_age", env_variable)),
        ~ !is.na(.x)
      )
    )

  model_formula <- bf(
    as.formula(
      paste(
        "mean_age ~",
        env_variable,
        "+ s(Year, k = 8)"
      )
    )
  )

  age_priors <- c(
    # Prior for your population-level effects (env_variable)
    # Assumes an environmental shift won't drastically swing age by more than a few years
    prior(normal(0, 5), class = "b"),

    # Prior for the overall Intercept (your average population age)
    # Adjust the '10' to roughly match your species' average age
    prior(normal(8, 3), class = "Intercept"),

    # Prior for the standard deviation of the spline (Year smooth)
    # Keeps the year-to-year wiggliness from exploding into chaos
    prior(exponential(1), class = "sds"),

    # Prior for sigma (residual standard deviation)
    prior(exponential(1), class = "sigma"),

    # prior for nu (degrees of freedom for the Student's t distribution)
    prior(gamma(2, 0.1), class = "nu")
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
    prior = age_priors,
    chains = 4,
    iter = 3000,
    warmup = 1500,
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
        data_source,
        "_annual_mean_age_",
        env_variable,
        "_",
        year_lag
      )
    ),
    file_refit = "always",
    refresh = 500
  )
}

# -----------------------------
# 4.5 fit each env variable separately
# -----------------------------

# list of env variables
env_var_names <- c(
  "pdo_jantomar", "pdo_aprtojun", "pdo_jultosep", "pdo_octtodec",
  "npi_jantomar", "npi_aprtojun", "npi_jultosep", "npi_octtodec",
  "npgo_jantomar", "npgo_aprtojun", "npgo_jultosep", "npgo_octtodec",
  "oni_jantomar", "oni_aprtojun", "oni_jultosep", "oni_octtodec",
  "mei_jantomar", "mei_aprtojun", "mei_jultosep", "mei_octtodec"
)


# run model for current year
# scale env variables
annual_age_env_current_scaled <- annual_age_env_current %>%
  mutate(
    across(
      all_of(env_var_names),
      ~ as.vector(scale(.x)),
      .names = "{.col}"
    )
  )

# run model
annual_age_env_models_current <- env_var_names %>%
  set_names() %>%
  map(
    \(env_variable) {
      fit_annual_age_env_model(
        data = annual_age_env_current_scaled,
        env_variable = env_variable,
        year_lag = "current"
      )
    }
  )

# run model for lag 1 year
# scale env variables
annual_age_env_lag1_scaled <- annual_age_env_lag1 %>%
  mutate(
    across(
      all_of(env_var_names),
      ~ as.vector(scale(.x)),
      .names = "{.col}"
    )
  )

# run model
annual_age_env_models_lag1 <- env_var_names %>%
  set_names() %>%
  map(
    \(env_variable) {
      fit_annual_age_env_model(
        data = annual_age_env_lag1_scaled,
        env_variable = env_variable,
        year_lag = "lag1"
      )
    }
  )

# run model for lag 2 year
# scale env variables
annual_age_env_lag2_scaled <- annual_age_env_lag2 %>%
  mutate(
    across(
      all_of(env_var_names),
      ~ as.vector(scale(.x)),
      .names = "{.col}"
    )
  )

# run model
annual_age_env_models_lag2 <- env_var_names %>%
  set_names() %>%
  map(
    \(env_variable) {
      fit_annual_age_env_model(
        data = annual_age_env_lag2_scaled,
        env_variable = env_variable,
        year_lag = "lag2"
      )
    }
  )

# run model for lag 3 year
# scale env variables
annual_age_env_lag3_scaled <- annual_age_env_lag3 %>%
  mutate(
    across(
      all_of(env_var_names),
      ~ as.vector(scale(.x)),
      .names = "{.col}"
    )
  )

# run model
annual_age_env_models_lag3 <- env_var_names %>%
  set_names() %>%
  map(
    \(env_variable) {
      fit_annual_age_env_model(
        data = annual_age_env_lag3_scaled,
        env_variable = env_variable,
        year_lag = "lag3"
      )
    }
  )


# combine outputs and paste year_lag name for each variable name
names(annual_age_env_models_current) <- paste0(env_var_names, "_current")
names(annual_age_env_models_lag1) <- paste0(env_var_names, "_lag1")
names(annual_age_env_models_lag2) <- paste0(env_var_names, "_lag2")
names(annual_age_env_models_lag3) <- paste0(env_var_names, "_lag3")

annual_age_env_models_all <- c(
  annual_age_env_models_current,
  annual_age_env_models_lag1,
  annual_age_env_models_lag2,
  annual_age_env_models_lag3
)

# check model diagnostics
annual_age_env_diagnostics <- imap_dfr(
  annual_age_env_models_all,
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
# 4.6 extract possterior from each model
# -----------------------------

# function to extract across all models
extract_annual_env_draws <- function(
  model,
  env_variable
) {
  coefficient_name <- paste0(
    "b_",
    sub("^(.*)_(current|lag1|lag2|lag3)$", "\\1", env_variable)
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
  annual_age_env_models_all,
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
      grepl("^pdo", covariate) ~ "PDO",
      grepl("^npgo", covariate) ~ "NPGO",
      grepl("^oni", covariate) ~ "ONI",
      grepl("^bi", covariate) ~ "BI",
      grepl("^npi", covariate) ~ "NPI",
      grepl("^mei", covariate) ~ "MEI",
      TRUE ~ covariate
    ),
    exposure_window = case_when(
      grepl("_jantomar", covariate) ~
        "Jan to Mar",
      grepl("_aprtojun", covariate) ~
        "Apr to Jun",
      grepl("_jultosep", covariate) ~
        "Jul to Sep",
      grepl("_octtodec", covariate) ~
        "Oct to Dec",
      TRUE ~ "Other"
    ),
    exposure_window = factor(
      exposure_window,
      levels = c(
        "Oct to Dec",
        "Jul to Sep",
        "Apr to Jun",
        "Jan to Mar",
        "Other"
      )
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
    year_lag = case_when(
      grepl("_current$", covariate) ~ "current",
      grepl("_lag1$", covariate) ~ "lag1",
      grepl("_lag2$", covariate) ~ "lag2",
      grepl("_lag3$", covariate) ~ "lag3",
      TRUE ~ "Other"
    )
  )

# -----------------------------
# 4.7 summarise posterior coefficints
# -----------------------------

annual_age_env_summary <- annual_age_env_draws %>%
  group_by(
    covariate,
    index,
    exposure_window,
    year_lag
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

# output table of posterior coefficient summaries for each environmental variable
dir.create(paste0("model-outputs/", data_source), recursive = TRUE)
write.csv(
  annual_age_env_summary,
  file = paste0("model-outputs/", data_source, "/annual_age_env_summary.csv"),
  row.names = FALSE
)

# -----------------------------
# 4.8 plot all posterior coefficient distributions
# -----------------------------

# plot halfeye and save in object
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
    ~ year_lag + index,
    nrow = 4
    # scales = "free_x"
  ) +
  labs(
    title = "Previous year seasonal environment associations with annual mean age",
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

# save figure
ggsave(
  filename = paste0("figures/", data_source, "/annual_age_env_posterior_plot.png"),
  plot = annual_age_env_posterior_plot,
  width = 8,
  height = 6
)

# plot density distribution curve with year_lags on the same plot

annual_age_env_posterior_density_plot <- ggplot(
  annual_age_env_draws,
  aes(
    x = beta,
    y = exposure_window,
    colour = factor(year_lag)
  )
) +
  stat_slab(
    fill = NA
  ) +
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    linewidth = 0.6
  ) +
  facet_wrap(
    ~index,
    ncol = 3,
    scales = "free_x"
  ) +
  labs(
    title = "Posterior density of environmental effects on annual mean age",
    x = "Change in annual mean age per 1-SD increase",
    y = "Environmental exposure window"
  ) +
  theme_classic(base_size = 12)

annual_age_env_posterior_density_plot



# -----------------------------
# 4.9 forest plot of all posteriors
# -----------------------------

# function to plot forest plot, then done for each year_lag option
annual_age_env_forest_plot <- function(
  data = annual_age_env_summary,
  lag = "current"
) {
  data <- data %>%
    filter(
      year_lag == lag
    )

  ggplot(
    data,
    aes(
      x = mean_beta,
      y = exposure_window
    )
  ) +
    geom_vline(
      xintercept = 0,
      linetype = "dashed",
      linewidth = 0.6
    ) +
    geom_errorbar(
      aes(
        xmin = lower_beta,
        xmax = upper_beta
      ),
      width = 0.3
    ) +
    geom_point(
      size = 2.5
    ) +
    facet_wrap(
      vars(index),
      ncol = 2,
      scales = "free_x"
    ) +
    labs(
      title = paste0("Posterior environmental effects (", lag, " year) on annual mean age"),
      subtitle = "Points are means; lines are 95% credible intervals",
      x = "Change in annual mean age (years)",
      y = NULL
    ) +
    theme_classic(base_size = 12)
}

# plot current year_lag and save
annual_age_env_forest_plot(data = annual_age_env_summary, lag = "current")
ggsave(
  filename = paste0("figures/", data_source, "/annual_age_env_forest_plot_current.png"),
  plot = last_plot(),
  width = 8,
  height = 6
)

# plot lag1 year_lag
annual_age_env_forest_plot(data = annual_age_env_summary, lag = "lag1")
ggsave(
  filename = paste0("figures/", data_source, "/annual_age_env_forest_plot_lag1.png"),
  plot = last_plot(),
  width = 8,
  height = 6
)

# plot lag2 year_lag
annual_age_env_forest_plot(data = annual_age_env_summary, lag = "lag2")
ggsave(
  filename = paste0("figures/", data_source, "/annual_age_env_forest_plot_lag2.png"),
  plot = last_plot(),
  width = 8,
  height = 6
)

# plot lag3 year_lag
annual_age_env_forest_plot(data = annual_age_env_summary, lag = "lag3")
ggsave(
  filename = paste0("figures/", data_source, "/annual_age_env_forest_plot_lag3.png"),
  plot = last_plot(),
  width = 8,
  height = 6
)

# -----------------------------
# 4.10 Ribbon plot of posteriors
# -----------------------------

# order year lag and season windows chronologically in the summary data

annual_age_env_summary_ribbon <- annual_age_env_summary %>%
  mutate(
    year_lag = factor(
      year_lag,
      levels = c("lag3", "lag2", "lag1", "current")
    ),
    exposure_window = factor(
      exposure_window,
      levels = c("Jan to Mar", "Apr to Jun", "Jul to Sep", "Oct to Dec")
    ),
    time = as.numeric(as.factor(paste0(as.numeric(year_lag), as.numeric(exposure_window))))
  ) %>%
  arrange(year_lag, exposure_window)


str(annual_age_env_summary_ribbon)

ggplot(
  annual_age_env_summary_ribbon,
  aes(
    x = time,
    y = mean_beta
  )
) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    linewidth = 0.6
  ) +
  geom_ribbon(
    aes(
      ymin = lower_beta,
      ymax = upper_beta
    ),
    alpha = 0.3,
    colour = NA
  ) +
  geom_line(
    size = 1
  ) +
  geom_point(
    size = 2
  ) +
  facet_wrap(
    vars(index)
  ) +
  scale_x_continuous(
    breaks = unique(annual_age_env_summary_ribbon$time),
    labels = paste0(
      rep(levels(annual_age_env_summary_ribbon$year_lag),
        each = length(levels(annual_age_env_summary_ribbon$exposure_window))
      ),
      " - ",
      substr(levels(annual_age_env_summary_ribbon$exposure_window), 1, 3)
    )
  ) +
  labs(
    title = "Posterior environmental effects on annual mean age",
    subtitle = "Points are means; ribbons are 95% credible intervals",
    x = "Time",
    y = "Change in annual mean age (years) with +1 SD increase in environmental variable"
  ) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1)
  )

ggsave(
  filename = paste0("figures/", data_source, "/annual_age_env_ribbon_plot.png"),
  plot = last_plot(),
  width = 8,
  height = 6
)

# -----------------------------
# 4.11 testing leftover memory (residuals)
# -----------------------------

names(annual_age_env_models_all)

fit <- annual_age_env_models_all[[42]]

# 2. Extract residuals and make sure they are ordered by Year
res_df <- fit$data %>%
  mutate(
    residual = residuals(fit)[, "Estimate"]
  ) %>%
  arrange(Year)

# 3. Plot the Autocorrelation Function (ACF) of the residuals
acf(
  res_df$residual,
  main = "Residual Autocorrelation (Leftover Memory)"
)

# Check 1: Overall distribution (Does the model capture sharp cohort drops and old-age peaks?)
pp_check(fit, ndraws = 100) +
  ggtitle("Posterior Predictive Check: Density")

# Check 2: Time-series trajectory (Does the model track year-to-year shifts accurately?)
pp_check(fit, type = "ribbon", x = "Year") +
  ggtitle("Posterior Predictive Check: Mean Age Over Time")


##### -------------------------------------------------------------------- #####
##### -------------------------------------------------------------------- #####
##### -------------------------------------------------------------------- #####
##### -------------------------------------------------------------------- #####
##### -------------------------------------------------------------------- #####

# -----------------------------
# 4.21 Annual age environmental model with multiple lags and AR(3) autocorrelation in one model
# -----------------------------

# function to run annual age environmental model with multiple lags and AR(3) autocorrelation and a spline smoother of year reduced to k = 5
fit_annual_age_env_model_ar <- function(
  data,
  env_variable,
  data_source,
  output_directory = "model-outputs/annual_mean_age_env",
  seed = 22
) {
  # 1. Create lag 0, 1, and 2 columns chronologically before filtering NAs
  model_data <- data %>%
    arrange(Year) %>%
    mutate(
      "{env_variable}_lag0" := .data[[env_variable]],
      "{env_variable}_lag1" := lag(.data[[env_variable]], 1),
      "{env_variable}_lag2" := lag(.data[[env_variable]], 2),
      "{env_variable}_lag3" := lag(.data[[env_variable]], 3)
    ) %>%
    filter(
      if_all(
        all_of(c("Year", "mean_age", paste0(env_variable, "_lag0"), paste0(env_variable, "_lag1"), paste0(env_variable, "_lag2"))),
        ~ !is.na(.x)
      )
    )

  # 2. Combine all 4 lags, a broader spline (k = 5), and AR(3) autocorrelation
  model_formula <- bf(
    as.formula(
      paste0(
        "mean_age ~ ",
        env_variable, "_lag0 + ",
        env_variable, "_lag1 + ",
        env_variable, "_lag2 + ",
        env_variable, "_lag3 + ",
        "s(Year, k = 5) + ar(time = Year, p = 3)"
      )
    )
  )

  # 3. Updated priors (including tighter slope prior and AR prior)
  age_priors <- c(
    prior(normal(0, 1.5), class = "b"),
    prior(normal(0, 0.5), class = "ar"),
    prior(normal(8, 3), class = "Intercept"),
    prior(exponential(1), class = "sds"),
    prior(exponential(1), class = "sigma"),
    prior(gamma(2, 0.1), class = "nu")
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
    prior = age_priors,
    chains = 4,
    iter = 3000,
    warmup = 1500,
    backend = "cmdstanr",
    seed = seed,
    control = list(
      adapt_delta = 0.98,
      max_treedepth = 13
    ),
    save_pars = save_pars(
      all = TRUE
    ),
    file = file.path(
      output_directory,
      paste0(
        data_source,
        "_annual_mean_age_",
        env_variable,
        "_multilag_ar3"
      )
    ),
    file_refit = "always",
    refresh = 500
  )
}

# scaled data
head(annual_age_env_current_scaled)

# run model
annual_age_env_models_current_ar <- env_var_names %>%
  set_names() %>%
  map(
    \(env_variable) {
      fit_annual_age_env_model_ar(
        data = annual_age_env_current_scaled,
        env_variable = env_variable,
        data_source = data_source
      )
    }
  )


# -----------------------------
# 4.22 testing leftover memory (residuals) for new model
# -----------------------------

names(annual_age_env_models_current_ar)

fit <- annual_age_env_models_current_ar[[1]]

# 2. Extract residuals and make sure they are ordered by Year
res_df <- fit$data %>%
  mutate(
    residual = residuals(fit)[, "Estimate"]
  ) %>%
  arrange(Year)

# 3. Plot the Autocorrelation Function (ACF) of the residuals
acf(
  res_df$residual,
  main = "Residual Autocorrelation (Leftover Memory)"
)

# Check 1: Overall distribution (Does the model capture sharp cohort drops and old-age peaks?)
pp_check(fit, ndraws = 100) +
  ggtitle("Posterior Predictive Check: Density")

# Check 2: Time-series trajectory (Does the model track year-to-year shifts accurately?)
pp_check(fit, type = "ribbon", x = "Year") +
  ggtitle("Posterior Predictive Check: Mean Age Over Time")

# check model diagnostics
annual_age_env_diagnostics_ar <- imap_dfr(
  annual_age_env_models_current_ar,
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

annual_age_env_diagnostics_ar


# -----------------------------
# 4.23 extract possterior from each model
# -----------------------------

# function to extract across all models
extract_annual_env_draws <- function(
  model,
  env_variable
) {
  # Build the 3 expected brms parameter names for this environmental variable
  target_coefs <- paste0("b_", env_variable, "_lag", 0:3)

  model_draws <- as_draws_df(model)

  # Check that all 3 lag coefficients exist in the posterior draws
  missing_coefs <- setdiff(target_coefs, names(model_draws))
  if (length(missing_coefs) > 0) {
    stop(
      "Coefficient(s) not found: ",
      paste(missing_coefs, collapse = ", "),
      call. = FALSE
    )
  }

  # Extract and pivot into a tidy tibble matching your original output structure
  model_draws %>%
    as_tibble() %>%
    select(all_of(target_coefs)) %>%
    pivot_longer(
      cols = everything(),
      names_to = "covariate",
      names_prefix = "^b_",
      values_to = "beta"
    ) %>%
    mutate(
      env_variable = env_variable,
      lag = sub(".*_(lag[0-3])$", "\\1", covariate),
      .before = covariate
    )
}

# run function
annual_age_env_draws_ar <- imap_dfr(
  annual_age_env_models_current_ar,
  \(model, variable) {
    extract_annual_env_draws(
      model = model,
      env_variable = variable
    )
  }
)

# add labels for plotting
annual_age_env_draws_ar <- annual_age_env_draws_ar %>%
  mutate(
    index = case_when(
      grepl("^pdo", covariate) ~ "PDO",
      grepl("^npgo", covariate) ~ "NPGO",
      grepl("^oni", covariate) ~ "ONI",
      grepl("^bi", covariate) ~ "BI",
      grepl("^npi", covariate) ~ "NPI",
      grepl("^mei", covariate) ~ "MEI",
      TRUE ~ covariate
    ),
    exposure_window = case_when(
      grepl("_jantomar", covariate) ~
        "Jan to Mar",
      grepl("_aprtojun", covariate) ~
        "Apr to Jun",
      grepl("_jultosep", covariate) ~
        "Jul to Sep",
      grepl("_octtodec", covariate) ~
        "Oct to Dec",
      TRUE ~ "Other"
    ),
    exposure_window = factor(
      exposure_window,
      levels = c(
        "Oct to Dec",
        "Jul to Sep",
        "Apr to Jun",
        "Jan to Mar",
        "Other"
      )
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
    year_lag = case_when(
      grepl("_lag0$", covariate) ~ "lag0",
      grepl("_lag1$", covariate) ~ "lag1",
      grepl("_lag2$", covariate) ~ "lag2",
      grepl("_lag3$", covariate) ~ "lag3",
      TRUE ~ "Other"
    )
  )

# -----------------------------
# 4.24 summarise posterior coefficints
# -----------------------------

annual_age_env_summary_ar <- annual_age_env_draws_ar %>%
  group_by(
    covariate,
    index,
    exposure_window,
    year_lag
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
annual_age_env_summary_ar

# -----------------------------
# 4.25 plot all posterior coefficient distributions
# -----------------------------

# plot halfeye and save in object
annual_age_env_posterior_plot_ar <- ggplot(
  annual_age_env_draws_ar,
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
  # facet_wrap(
  #   ~ year_lag + index,
  #   nrow = 3
  # ) +
  facet_grid(
    year_lag ~ index,
    scales = "free_x"
  ) +
  labs(
    title = "Previous year seasonal environment associations with annual mean age",
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

annual_age_env_posterior_plot_ar

# save figure
ggsave(
  filename = paste0("figures/", data_source, "/annual_age_env_posterior_plot_ar.png"),
  plot = annual_age_env_posterior_plot_ar,
  width = 8,
  height = 6
)

# plot density distribution curve with year_lags on the same plot

annual_age_env_posterior_density_plot_ar <- ggplot(
  annual_age_env_draws_ar,
  aes(
    x = beta,
    y = exposure_window,
    colour = factor(year_lag)
  )
) +
  stat_slab(
    fill = NA
  ) +
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    linewidth = 0.6
  ) +
  facet_wrap(
    ~index,
    ncol = 3,
    scales = "free_x"
  ) +
  labs(
    title = "Posterior density of environmental effects on annual mean age",
    x = "Change in annual mean age per 1-SD increase",
    y = "Environmental exposure window"
  ) +
  theme_classic(base_size = 12)

annual_age_env_posterior_density_plot_ar

# save figure
ggsave(
  filename = paste0("figures/", data_source, "/annual_age_env_posterior_density_plot_ar.png"),
  plot = annual_age_env_posterior_density_plot_ar,
  width = 8,
  height = 6
)

# -----------------------------
# 4.26 forest plot of all posteriors
# -----------------------------

ggplot(
  annual_age_env_summary_ar,
  aes(
    x = mean_beta,
    y = exposure_window
  )
) +
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    linewidth = 0.6
  ) +
  geom_errorbar(
    aes(
      xmin = lower_beta,
      xmax = upper_beta
    ),
    width = 0.3
  ) +
  geom_point(
    size = 2.5
  ) +
  facet_grid(
    year_lag ~ index,
    scales = "free_x"
  ) +
  labs(
    title = paste0("Posterior environmental effects on annual mean age"),
    subtitle = "Points are means; lines are 95% credible intervals",
    x = "Change in annual mean age (years)",
    y = NULL
  ) +
  theme_classic(base_size = 12)

# save plot
ggsave(
  filename = paste0("figures/", data_source, "/annual_age_env_forest_plot_ar.png"),
  plot = last_plot(),
  width = 8,
  height = 6
)

# -----------------------------
# 4.27 Ribbon plot of posteriors
# -----------------------------

# order year lag and season windows chronologically in the summary data

annual_age_env_summary_ribbon_ar <- annual_age_env_summary_ar %>%
  mutate(
    year_lag = factor(
      year_lag,
      levels = c("lag2", "lag1", "lag0")
    ),
    exposure_window = factor(
      exposure_window,
      levels = c("Jan to Mar", "Apr to Jun", "Jul to Sep", "Oct to Dec")
    ),
    time = as.numeric(as.factor(paste0(as.numeric(year_lag), as.numeric(exposure_window))))
  ) %>%
  arrange(year_lag, exposure_window)


str(annual_age_env_summary_ribbon_ar)

ggplot(
  annual_age_env_summary_ribbon_ar,
  aes(
    x = time,
    y = mean_beta
  )
) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    linewidth = 0.6
  ) +
  geom_ribbon(
    aes(
      ymin = lower_beta,
      ymax = upper_beta
    ),
    alpha = 0.3,
    colour = NA
  ) +
  geom_line(
    size = 1
  ) +
  geom_point(
    size = 2
  ) +
  facet_wrap(
    vars(index)
  ) +
  scale_x_continuous(
    breaks = unique(annual_age_env_summary_ribbon$time),
    labels = paste0(
      rep(levels(annual_age_env_summary_ribbon$year_lag),
        each = length(levels(annual_age_env_summary_ribbon$exposure_window))
      ),
      " - ",
      substr(levels(annual_age_env_summary_ribbon$exposure_window), 1, 3)
    )
  ) +
  labs(
    title = "Posterior environmental effects on annual mean age",
    subtitle = "Points are means; ribbons are 95% credible intervals",
    x = "Time",
    y = "Change in annual mean age (years) with +1 SD increase in environmental variable"
  ) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1)
  )

# save figure
ggsave(
  filename = paste0("figures/", data_source, "/annual_age_env_ribbon_plot_ar.png"),
  plot = last_plot(),
  width = 8,
  height = 6
)

# -----------------------------
# 4.28 Cumulative lag effects (across lag0,lag1,lag2)for PDO, NPGO, MEI, ONI, and NPI
# -----------------------------

# PDO lag
# Check cumulative 4-year PDO effect (lag0 + lag1 + lag2 + lag3)
cumulative_pdo <- annual_age_env_models_current_ar[grep("^pdo_", names(annual_age_env_models_current_ar))] %>%
  imap(\(model, env_var) {
    draws <- as_draws_df(model)
    b0 <- draws[[paste0("b_", env_var, "_lag0")]]
    b1 <- draws[[paste0("b_", env_var, "_lag1")]]
    b2 <- draws[[paste0("b_", env_var, "_lag2")]]
    b3 <- draws[[paste0("b_", env_var, "_lag3")]]

    tibble(
      covariate = env_var,
      cumulative_beta = b0 + b1 + b2 + b3
    )
  }) %>%
  list_rbind() %>%
  group_by(covariate) %>%
  summarise(
    median_sum = median(cumulative_beta),
    lower_95 = quantile(cumulative_beta, 0.025),
    upper_95 = quantile(cumulative_beta, 0.975),
    prob_positive = mean(cumulative_beta > 0),
    prob_negative = mean(cumulative_beta < 0)
  )

# NPGO lag
# Check cumulative 4-year NPGO effect (lag0 + lag1 + lag2 + lag3)
cumulative_npgo <- annual_age_env_models_current_ar[grep("^npgo_", names(annual_age_env_models_current_ar))] %>%
  imap(\(model, env_var) {
    draws <- as_draws_df(model)
    b0 <- draws[[paste0("b_", env_var, "_lag0")]]
    b1 <- draws[[paste0("b_", env_var, "_lag1")]]
    b2 <- draws[[paste0("b_", env_var, "_lag2")]]
    b3 <- draws[[paste0("b_", env_var, "_lag3")]]

    tibble(
      covariate = env_var,
      cumulative_beta = b0 + b1 + b2 + b3
    )
  }) %>%
  list_rbind() %>%
  group_by(covariate) %>%
  summarise(
    median_sum = median(cumulative_beta),
    lower_95 = quantile(cumulative_beta, 0.025),
    upper_95 = quantile(cumulative_beta, 0.975),
    prob_positive = mean(cumulative_beta > 0),
    prob_negative = mean(cumulative_beta < 0)
  )

# MEI lag
# Check cumulative 4-year MEI effect (lag0 + lag1 + lag2 + lag3)
cumulative_mei <- annual_age_env_models_current_ar[grep("^mei_", names(annual_age_env_models_current_ar))] %>%
  imap(\(model, env_var) {
    draws <- as_draws_df(model)
    b0 <- draws[[paste0("b_", env_var, "_lag0")]]
    b1 <- draws[[paste0("b_", env_var, "_lag1")]]
    b2 <- draws[[paste0("b_", env_var, "_lag2")]]
    b3 <- draws[[paste0("b_", env_var, "_lag3")]]

    tibble(
      covariate = env_var,
      cumulative_beta = b0 + b1 + b2 + b3
    )
  }) %>%
  list_rbind() %>%
  group_by(covariate) %>%
  summarise(
    median_sum = median(cumulative_beta),
    lower_95 = quantile(cumulative_beta, 0.025),
    upper_95 = quantile(cumulative_beta, 0.975),
    prob_positive = mean(cumulative_beta > 0),
    prob_negative = mean(cumulative_beta < 0)
  )

# ONI lag
# Check cumulative 4-year ONI effect (lag0 + lag1 + lag2 + lag3)
cumulative_oni <- annual_age_env_models_current_ar[grep("^oni_", names(annual_age_env_models_current_ar))] %>%
  imap(\(model, env_var) {
    draws <- as_draws_df(model)
    b0 <- draws[[paste0("b_", env_var, "_lag0")]]
    b1 <- draws[[paste0("b_", env_var, "_lag1")]]
    b2 <- draws[[paste0("b_", env_var, "_lag2")]]
    b3 <- draws[[paste0("b_", env_var, "_lag3")]]

    tibble(
      covariate = env_var,
      cumulative_beta = b0 + b1 + b2 + b3
    )
  }) %>%
  list_rbind() %>%
  group_by(covariate) %>%
  summarise(
    median_sum = median(cumulative_beta),
    lower_95 = quantile(cumulative_beta, 0.025),
    upper_95 = quantile(cumulative_beta, 0.975),
    prob_positive = mean(cumulative_beta > 0),
    prob_negative = mean(cumulative_beta < 0)
  )

# NPI lag
# Check cumulative 4-year NPI effect (lag0 + lag1 + lag2 + lag3)
cumulative_npi <- annual_age_env_models_current_ar[grep("^npi_", names(annual_age_env_models_current_ar))] %>%
  imap(\(model, env_var) {
    draws <- as_draws_df(model)
    b0 <- draws[[paste0("b_", env_var, "_lag0")]]
    b1 <- draws[[paste0("b_", env_var, "_lag1")]]
    b2 <- draws[[paste0("b_", env_var, "_lag2")]]
    b3 <- draws[[paste0("b_", env_var, "_lag3")]]

    tibble(
      covariate = env_var,
      cumulative_beta = b0 + b1 + b2 + b3
    )
  }) %>%
  list_rbind() %>%
  group_by(covariate) %>%
  summarise(
    median_sum = median(cumulative_beta),
    lower_95 = quantile(cumulative_beta, 0.025),
    upper_95 = quantile(cumulative_beta, 0.975),
    prob_positive = mean(cumulative_beta > 0),
    prob_negative = mean(cumulative_beta < 0)
  )


# -----------------------------
# 4.29 Plot cumulative lag effects for PDO, NPGO, MEI, ONI, and NPI
# -----------------------------

# combine cumulative lag effects for all environmental variables
cumulative_all <- bind_rows(
  cumulative_pdo,
  cumulative_npgo,
  cumulative_mei,
  cumulative_oni,
  cumulative_npi
)

# create forest ggplot figure of cumulative effects, faceted by index
cumulative_all %>%
  mutate(
    index = case_when(
      grepl("^pdo", covariate) ~ "PDO",
      grepl("^npgo", covariate) ~ "NPGO",
      grepl("^oni", covariate) ~ "ONI",
      grepl("^bi", covariate) ~ "BI",
      grepl("^npi", covariate) ~ "NPI",
      grepl("^mei", covariate) ~ "MEI",
      TRUE ~ covariate
    ),
    exposure_window = case_when(
      grepl("_jantomar", covariate) ~
        "Jan to Mar",
      grepl("_aprtojun", covariate) ~
        "Apr to Jun",
      grepl("_jultosep", covariate) ~
        "Jul to Sep",
      grepl("_octtodec", covariate) ~
        "Oct to Dec",
      TRUE ~ "Other"
    ),
    exposure_window = factor(exposure_window, levels = c("Jan to Mar", "Apr to Jun", "Jul to Sep", "Oct to Dec", "Other"))
  ) %>%
  ggplot(aes(x = exposure_window, y = median_sum)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  geom_point() +
  geom_errorbar(aes(ymin = lower_95, ymax = upper_95), width = 0.2) +
  facet_wrap(~index, scales = "free_x") +
  theme_bw() +
  labs(
    x = "Exposure window",
    y = "Cumulative Effect"
  )


##### -------------------------------------------------------------------- #####
##### -------------------------------------------------------------------- #####
##### -------------------------------------------------------------------- #####
##### -------------------------------------------------------------------- #####
##### -------------------------------------------------------------------- #####

# Plot annual mean-age effects by index, exposure window, and lag, then evaluate fit and predictions.
# for the separate lag model

# -----------------------------
# 4.31 environmental effects plot - data organizing
# -----------------------------

# Model names include a time suffix; the predictor column names do not.
get_annual_model_context <- function(model_name) {
  if (!grepl("_(current|lag[1-3])$", model_name)) {
    stop(
      "Model name must end in _current or _lag1 through _lag3: ",
      model_name,
      call. = FALSE
    )
  }

  year_lag <- sub(
    "^.*_(current|lag[1-3])$",
    "\\1",
    model_name
  )
  env_variable <- sub(
    "_(current|lag[1-3])$",
    "",
    model_name
  )

  source_data <- switch(year_lag,
    current = annual_age_env_current_scaled,
    lag1 = annual_age_env_lag1_scaled,
    lag2 = annual_age_env_lag2_scaled,
    lag3 = annual_age_env_lag3_scaled,
    stop("Unsupported model time suffix: ", year_lag, call. = FALSE)
  )

  timing_label <- switch(year_lag,
    current = "Current year",
    lag1 = "Lag 1 year",
    lag2 = "Lag 2 years",
    lag3 = "Lag 3 years"
  )

  list(
    env_variable = env_variable,
    source_data = source_data,
    timing_label = timing_label
  )
}

get_model_data <- function(env_variable, source_data) {
  required_variables <- c(
    "Year",
    "mean_age",
    "n_fish",
    env_variable
  )
  missing_variables <- setdiff(required_variables, names(source_data))

  if (length(missing_variables) > 0L) {
    stop(
      "Missing variables in source_data: ",
      paste(missing_variables, collapse = ", "),
      call. = FALSE
    )
  }

  model_data <- source_data %>%
    filter(
      if_all(
        all_of(required_variables),
        ~ !is.na(.x)
      )
    ) %>%
    arrange(Year)

  if (nrow(model_data) == 0L) {
    stop(
      "No complete observations for environmental variable: ",
      env_variable,
      call. = FALSE
    )
  }

  model_data
}

get_annual_env_effect <- function(
  model,
  model_name,
  env_variable,
  source_data,
  prediction_length = 100
) {
  model_data <- get_model_data(
    env_variable = env_variable,
    source_data = source_data
  )
  env_values <- model_data[[env_variable]]

  prediction_data <- tibble(
    Year = rep(
      median(model_data$Year),
      prediction_length
    )
  )
  prediction_data[[env_variable]] <- seq(
    min(env_values),
    max(env_values),
    length.out = prediction_length
  )

  fitted_summary <- fitted(
    model,
    newdata = prediction_data,
    re_formula = NA,
    probs = c(0.025, 0.975)
  ) %>%
    as_tibble()

  bind_cols(
    tibble(env_value = prediction_data[[env_variable]]),
    fitted_summary
  ) %>%
    mutate(covariate = model_name)
}

# -----------------------------
# 4.32 function to plot one env relationship
# -----------------------------

plot_annual_age_env_effect <- function(
  model,
  model_name,
  env_variable,
  source_data
) {
  model_data <- get_model_data(
    env_variable = env_variable,
    source_data = source_data
  )

  effect_data <- get_annual_env_effect(
    model = model,
    model_name = model_name,
    env_variable = env_variable,
    source_data = source_data
  )

  model_context <- get_annual_model_context(model_name)
  index_label <- case_when(
    grepl("^pdo", env_variable) ~ "PDO",
    grepl("^npgo", env_variable) ~ "NPGO",
    grepl("^oni", env_variable) ~ "ONI",
    grepl("^bi", env_variable) ~ "BI",
    grepl("^npi", env_variable) ~ "NPI",
    grepl("^mei", env_variable) ~ "MEI",
    TRUE ~ env_variable
  )
  exposure_label <- case_when(
    grepl("_jantomar$", env_variable) ~ "Jan to Mar",
    grepl("_aprtojun$", env_variable) ~ "Apr to Jun",
    grepl("_jultosep$", env_variable) ~ "Jul to Sep",
    grepl("_octtodec$", env_variable) ~ "Oct to Dec",
    TRUE ~ env_variable
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
      title = paste(index_label, exposure_label),
      subtitle = model_context$timing_label,
      x = "Standardized environmental index",
      y = "Predicted annual mean age (years)",
      size = "Fish sampled"
    ) +
    theme_classic(base_size = 12)
}


# -----------------------------
# 4.33 effects plot for each env variable
# -----------------------------

annual_age_env_effect_plots <- imap(
  annual_age_env_models_all,
  \(model, model_name) {
    model_context <- get_annual_model_context(model_name)

    plot_annual_age_env_effect(
      model = model,
      model_name = model_name,
      env_variable = model_context$env_variable,
      source_data = model_context$source_data
    )
  }
)

combine_annual_age_env_plots_by_index <- function(plots_by_index) {
  imap(
    plots_by_index,
    \(index_plots, index) {
      index_plot <- wrap_plots(
        index_plots,
        ncol = 4,
        guides = "collect"
      ) +
        plot_annotation(
          title = toupper(index)
        )

      index_plot &
        theme(
          legend.position = "bottom"
        )
    }
  )
}

# Group by index. Within each index, columns are exposure windows and rows
# are current year through lag 3.
annual_age_env_effect_plots_by_index <- split(
  annual_age_env_effect_plots,
  factor(
    sub("_.*$", "", names(annual_age_env_effect_plots)),
    levels = c("pdo", "npgo", "oni", "npi", "mei")
  )
)

annual_age_env_effect_combined_by_index <-
  combine_annual_age_env_plots_by_index(
    annual_age_env_effect_plots_by_index
  )

walk(
  annual_age_env_effect_combined_by_index,
  print
)

# save each figure by index
walk2(
  annual_age_env_effect_combined_by_index,
  names(annual_age_env_effect_combined_by_index),
  ~ ggsave(
    filename = paste0("figures/", data_source, "/annual_age_env_effect_observed_with_posteriorpred_", .y, ".png"),
    plot = .x
  )
)


# -----------------------------
# 4.34 posterior predictions and model evaluation
# -----------------------------

get_annual_fitted_values <- function(
  model,
  model_name
) {
  model_context <- get_annual_model_context(model_name)
  model_data <- get_model_data(
    env_variable = model_context$env_variable,
    source_data = model_context$source_data
  )

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
    mutate(covariate = model_name)
}

# table of posterior values with error
annual_age_env_fitted <- imap_dfr(
  annual_age_env_models_all,
  \(model, model_name) {
    get_annual_fitted_values(
      model = model,
      model_name = model_name
    )
  }
)

# calculate residuals
annual_age_env_fitted <- annual_age_env_fitted %>%
  mutate(
    residual = mean_age - Estimate
  )

# observed and fitted values through time, grouped by index
annual_age_env_fitted_by_index <- split(
  annual_age_env_fitted,
  factor(
    sub("_.*$", "", annual_age_env_fitted$covariate),
    levels = c("pdo", "npgo", "oni", "npi", "mei")
  )
)

annual_age_env_time_series_plots_by_index <- imap(
  annual_age_env_fitted_by_index,
  \(index_data, index) {
    index_data$covariate <- factor(
      index_data$covariate,
      levels = names(annual_age_env_models_all)
    )

    ggplot(
      index_data,
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
        ~covariate,
        ncol = 4,
        labeller = labeller(covariate = \(x) gsub("_", " ", x))
      ) +
      labs(
        x = "Year",
        y = "Mean age",
        title = paste(
          "Observed and fitted annual mean age:",
          toupper(index)
        )
      ) +
      theme_classic()
  }
)

walk(
  annual_age_env_time_series_plots_by_index,
  print
)

# residuals through time, grouped by index
annual_age_env_residual_plots_by_index <- imap(
  annual_age_env_fitted_by_index,
  \(index_data, index) {
    index_data$covariate <- factor(
      index_data$covariate,
      levels = names(annual_age_env_models_all)
    )

    ggplot(
      index_data,
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
        ~covariate,
        ncol = 4,
        labeller = labeller(covariate = \(x) gsub("_", " ", x))
      ) +
      labs(
        x = "Year",
        y = "Observed - fitted mean age",
        title = paste("Residual patterns through time:", toupper(index))
      ) +
      theme_classic()
  }
)

walk(
  annual_age_env_residual_plots_by_index,
  print
)

# observed versus fitted, grouped by index
annual_age_env_observed_fitted_plots_by_index <- imap(
  annual_age_env_fitted_by_index,
  \(index_data, index) {
    index_data$covariate <- factor(
      index_data$covariate,
      levels = names(annual_age_env_models_all)
    )

    ggplot(
      index_data,
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
        ~covariate,
        ncol = 4,
        labeller = labeller(covariate = \(x) gsub("_", " ", x))
      ) +
      labs(
        x = "Fitted mean age",
        y = "Observed mean age",
        title = paste("Observed versus fitted annual mean age:", toupper(index))
      ) +
      theme_classic()
  }
)

walk(
  annual_age_env_observed_fitted_plots_by_index,
  print
)

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
# 4.35 function to plot observed vs fitted ages over time
# -----------------------------

plot_annual_fitted <- function(
  model,
  model_name
) {
  fitted_data <- get_annual_fitted_values(
    model = model,
    model_name = model_name
  )

  plot_title <- gsub("_", " ", model_name)

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
  annual_age_env_models_all,
  \(model, model_name) {
    plot_annual_fitted(
      model = model,
      model_name = model_name
    )
  }
)

# combine fitted-through-time plots by index
annual_age_env_fitted_model_plots_by_index <- split(
  annual_age_env_fitted_plots,
  factor(
    sub("_.*$", "", names(annual_age_env_fitted_plots)),
    levels = c("pdo", "npgo", "oni", "npi", "mei")
  )
)

annual_age_env_fitted_model_combined_by_index <-
  combine_annual_age_env_plots_by_index(
    annual_age_env_fitted_model_plots_by_index
  )

walk(
  annual_age_env_fitted_model_combined_by_index,
  print
)

# -----------------------------
# 4.36 posterior predictive intervals
# -----------------------------

# posterior predictive intervals
get_annual_predictive_summary <- function(
  model,
  model_name,
  ndraws = 2000
) {
  model_context <- get_annual_model_context(model_name)
  model_data <- get_model_data(
    env_variable = model_context$env_variable,
    source_data = model_context$source_data
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
    mutate(covariate = model_name)
}

set.seed(123)

annual_age_env_predictive <- imap_dfr(
  annual_age_env_models_all,
  \(model, model_name) {
    get_annual_predictive_summary(
      model = model,
      model_name = model_name,
      ndraws = 2000
    )
  }
)

annual_age_env_predictive

# posterior predictions with observations
plot_annual_predictive <- function(
  predictive_data,
  model_name
) {
  plot_title <- gsub("_", " ", model_name)

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
    \(predictive_data, model_name) {
      plot_annual_predictive(
        predictive_data = predictive_data,
        model_name = model_name
      )
    }
  )

annual_age_env_predictive_plots_by_index <- split(
  annual_age_env_predictive_plots,
  factor(
    sub("_.*$", "", names(annual_age_env_predictive_plots)),
    levels = c("pdo", "npgo", "oni", "npi", "mei")
  )
)

annual_age_env_predictive_combined_by_index <-
  combine_annual_age_env_plots_by_index(
    annual_age_env_predictive_plots_by_index
  )

walk(
  annual_age_env_predictive_combined_by_index,
  print
)

#####
#####
#####
#####
#####

# Evaluate the multi-lag autoregressive models and their posterior predictions.

# -----------------------------
# 4.41 autoregressive model data and conditional effects
# -----------------------------

get_annual_ar_model_data <- function(model, env_variable, age_data = annual_age) {
  model_data <- model$data
  lag_variables <- paste0(env_variable, "_lag", 0:3)
  required_variables <- c("Year", "mean_age", lag_variables)
  missing_variables <- setdiff(required_variables, names(model_data))
  missing_age_variables <- setdiff(c("Year", "mean_age", "n_fish"), names(age_data))

  if (length(missing_variables) > 0L) {
    stop("Missing variables in model data: ", paste(missing_variables, collapse = ", "), call. = FALSE)
  }
  if (length(missing_age_variables) > 0L) {
    stop("Missing variables in age data: ", paste(missing_age_variables, collapse = ", "), call. = FALSE)
  }
  if (anyDuplicated(age_data$Year)) {
    stop("Age data must have one row per year.", call. = FALSE)
  }

  age_rows <- match(model_data$Year, age_data$Year)
  if (anyNA(age_rows) || anyNA(age_data$n_fish[age_rows])) {
    stop("Missing fish counts for one or more model years.", call. = FALSE)
  }
  if (!isTRUE(all.equal(
    as.numeric(model_data$mean_age),
    as.numeric(age_data$mean_age[age_rows])
  ))) {
    stop("Age data do not match the observations used to fit the model.", call. = FALSE)
  }
  if (any(!complete.cases(model_data[, required_variables]))) {
    stop("Model data contain incomplete observations.", call. = FALSE)
  }

  model_data$n_fish <- age_data$n_fish[age_rows]
  model_data
}

get_annual_ar_effect <- function(model, env_variable, lag, model_data) {
  lag_variables <- paste0(env_variable, "_lag", 0:3)
  effect_variable <- lag_variables[lag + 1]
  reference_data <- as_tibble(model_data[1, c("Year", lag_variables)])
  reference_data$Year <- median(model_data$Year)
  for (variable in lag_variables) {
    reference_data[[variable]] <- median(model_data[[variable]])
  }

  effect <- brms::conditional_effects(
    model,
    effects = effect_variable,
    conditions = reference_data,
    re_formula = NA,
    robust = FALSE,
    prob = 0.95
  )[[1]]

  tibble(
    covariate = env_variable,
    lag = paste0("lag", lag),
    env_value = effect[[effect_variable]],
    Estimate = effect$estimate__,
    lower_95 = effect$lower__,
    upper_95 = effect$upper__
  )
}

annual_age_env_ar_model_data <- imap(
  annual_age_env_models_current_ar,
  \(model, variable) get_annual_ar_model_data(model, variable)
)

annual_age_env_ar_effects <- imap_dfr(
  annual_age_env_models_current_ar,
  \(model, variable) {
    map_dfr(
      0:3,
      \(lag) get_annual_ar_effect(model, variable, lag, annual_age_env_ar_model_data[[variable]])
    )
  }
)

annual_age_env_ar_effect_observations <- imap_dfr(
  annual_age_env_ar_model_data,
  \(model_data, variable) {
    map_dfr(
      0:3,
      \(lag) {
        env_value <- model_data[[paste0(variable, "_lag", lag)]]
        tibble(
          covariate = variable,
          lag = paste0("lag", lag),
          env_value = env_value,
          mean_age = model_data$mean_age,
          n_fish = model_data$n_fish
        )
      }
    )
  }
)

# -----------------------------
# 4.42 autoregressive conditional effect plots by index
# -----------------------------

annual_age_env_ar_effect_plots_by_index <- imap(
  split(
    annual_age_env_ar_effects,
    sub("_.*$", "", annual_age_env_ar_effects$covariate)
  ),
  \(index_effects, index) {
    index_effects$covariate <- factor(index_effects$covariate, levels = names(annual_age_env_models_current_ar))
    index_effects$lag <- factor(index_effects$lag, levels = paste0("lag", 0:3))
    index_observations <- annual_age_env_ar_effect_observations %>%
      filter(covariate %in% as.character(index_effects$covariate)) %>%
      mutate(
        covariate = factor(covariate, levels = names(annual_age_env_models_current_ar)),
        lag = factor(lag, levels = paste0("lag", 0:3))
      )

    ggplot(index_effects, aes(x = env_value, y = Estimate)) +
      geom_ribbon(aes(ymin = lower_95, ymax = upper_95), alpha = 0.2) +
      geom_line() +
      geom_point(
        data = index_observations,
        aes(x = env_value, y = mean_age, size = n_fish),
        inherit.aes = FALSE,
        alpha = 0.5
      ) +
      facet_grid(lag ~ covariate, scales = "free_x", labeller = label_both) +
      scale_size_continuous(trans = "sqrt") +
      labs(
        title = paste("Autoregressive model environmental effects:", toupper(index)),
        subtitle = "Other lags and year held at their medians; AR contribution excluded. Points are unadjusted observations.",
        x = "Standardized environmental index",
        y = "Predicted annual mean age (years)",
        size = "Fish sampled"
      ) +
      theme_classic()
  }
)

walk(annual_age_env_ar_effect_plots_by_index, print)

# -----------------------------
# 4.43 autoregressive fitted values and model evaluation
# -----------------------------

annual_age_env_ar_fitted <- imap_dfr(
  annual_age_env_models_current_ar,
  \(model, variable) {
    fitted_values <- fitted(model, re_formula = NA, probs = c(0.025, 0.975)) %>%
      as_tibble()

    bind_cols(
      annual_age_env_ar_model_data[[variable]] %>% select(Year, mean_age, n_fish),
      fitted_values
    ) %>%
      mutate(covariate = variable, residual = mean_age - Estimate)
  }
)

annual_age_env_ar_fit_summary <- annual_age_env_ar_fitted %>%
  group_by(covariate) %>%
  summarise(
    n = n(),
    RMSE = sqrt(mean(residual^2)),
    MAE = mean(abs(residual)),
    bias = mean(residual),
    .groups = "drop"
  ) %>%
  arrange(RMSE)

annual_age_env_ar_fit_summary

annual_age_env_ar_fitted_by_index <- split(
  annual_age_env_ar_fitted,
  sub("_.*$", "", annual_age_env_ar_fitted$covariate)
)

annual_age_env_ar_fitted_plots_by_index <- imap(
  annual_age_env_ar_fitted_by_index,
  \(index_data, index) {
    ggplot(index_data, aes(x = Year)) +
      geom_ribbon(aes(ymin = Q2.5, ymax = Q97.5), alpha = 0.2) +
      geom_line(aes(y = Estimate), linewidth = 0.8) +
      geom_point(aes(y = mean_age, size = n_fish), alpha = 0.65) +
      scale_size_continuous(trans = "sqrt") +
      facet_wrap(~covariate, ncol = 2) +
      labs(
        title = paste("Autoregressive model observed and fitted:", toupper(index)),
        subtitle = "Line is the posterior expected value; ribbon is the 95% credible interval",
        x = "Year",
        y = "Annual mean age (years)",
        size = "Fish sampled"
      ) +
      theme_classic()
  }
)

walk(annual_age_env_ar_fitted_plots_by_index, print)

annual_age_env_ar_residual_plots_by_index <- imap(
  annual_age_env_ar_fitted_by_index,
  \(index_data, index) {
    ggplot(index_data, aes(x = Year, y = residual)) +
      geom_hline(yintercept = 0, linetype = "dashed") +
      geom_point() +
      geom_smooth(method = "loess", se = FALSE) +
      facet_wrap(~covariate, ncol = 2) +
      labs(
        title = paste("Autoregressive model residuals:", toupper(index)),
        x = "Year",
        y = "Observed - fitted mean age"
      ) +
      theme_classic()
  }
)

walk(annual_age_env_ar_residual_plots_by_index, print)

annual_age_env_ar_observed_fitted_plots_by_index <- imap(
  annual_age_env_ar_fitted_by_index,
  \(index_data, index) {
    ggplot(index_data, aes(x = Estimate, y = mean_age)) +
      geom_point() +
      geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
      facet_wrap(~covariate, ncol = 2) +
      labs(
        title = paste("Autoregressive model observed versus fitted:", toupper(index)),
        x = "Fitted mean age",
        y = "Observed mean age"
      ) +
      theme_classic()
  }
)

walk(annual_age_env_ar_observed_fitted_plots_by_index, print)

# -----------------------------
# 4.44 autoregressive posterior predictive intervals
# -----------------------------

set.seed(123)

annual_age_env_ar_predictive <- imap_dfr(
  annual_age_env_models_current_ar,
  \(model, variable) {
    predictive_draws <- posterior_predict(model, re_formula = NA, ndraws = 2000)
    predictive_summary <- apply(
      predictive_draws,
      2,
      \(x) c(
        mean = mean(x),
        lower_50 = unname(quantile(x, 0.25)),
        upper_50 = unname(quantile(x, 0.75)),
        lower_95 = unname(quantile(x, 0.025)),
        upper_95 = unname(quantile(x, 0.975))
      )
    ) %>%
      t() %>%
      as_tibble()

    bind_cols(
      annual_age_env_ar_model_data[[variable]] %>% select(Year, mean_age, n_fish),
      predictive_summary
    ) %>%
      mutate(covariate = variable)
  }
)

annual_age_env_ar_predictive_plots_by_index <- imap(
  split(
    annual_age_env_ar_predictive,
    sub("_.*$", "", annual_age_env_ar_predictive$covariate)
  ),
  \(index_data, index) {
    ggplot(index_data, aes(x = Year)) +
      geom_ribbon(aes(ymin = lower_95, ymax = upper_95), alpha = 0.15) +
      geom_ribbon(aes(ymin = lower_50, ymax = upper_50), alpha = 0.30) +
      geom_line(aes(y = mean), linewidth = 0.9) +
      geom_point(aes(y = mean_age, size = n_fish), alpha = 0.70) +
      scale_size_continuous(trans = "sqrt") +
      facet_wrap(~covariate, ncol = 2) +
      labs(
        title = paste("Autoregressive model posterior predictions:", toupper(index)),
        subtitle = "Dark ribbon is the 50% predictive interval; light ribbon is the 95% predictive interval",
        x = "Year",
        y = "Annual mean age (years)",
        size = "Fish sampled"
      ) +
      theme_classic()
  }
)

walk(annual_age_env_ar_predictive_plots_by_index, print)
