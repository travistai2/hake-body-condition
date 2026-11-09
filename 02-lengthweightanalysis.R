#####
#
# Hake length weight analysis and body condition estimation
#
#


##### -------------------------------------------------------------------- #####
##### -------------------------------------------------------------------- #####
##### -------------------------------------------------------------------- #####
##### -------------------------------------------------------------------- #####
##### -------------------------------------------------------------------- #####

##### Calculate le cren's body condition

# Length-weight mixed effects analysis + Le Cren condition + pacea covariates

# -----------------------------
# 2.1 Clean data
# -----------------------------

# distribution of data by age
hist(lwa_dat$Age_yrs, breaks = 20, main = "Distribution of data by age", xlab = "Age (yrs)")
table(lwa_dat$Age_yrs)
# combine older years into one category of 19 or older 

# distribution of data by month
hist(lwa_dat$Month, breaks = 20, main = "Distribution of data by month", xlab = "Month")
table(lwa_dat$Month)

# currently using the length weight data that seems to have the most data (a combination of the other datasets
lwa <- lwa_dat %>%
  filter(
    !outlierL,
    !is.na(Weight_kg),
    !is.na(Length_cm),
    Weight_kg > 0,
    Length_cm > 0,
    Sex %in% c("F","M")
  ) %>%
  mutate(
    Sex = factor(Sex),
    Age_yrs = ifelse(Age_yrs >= 19, 19, Age_yrs),
    Age_yrs_numeric = Age_yrs, 
    Age_yrs = factor(Age_yrs),
    Source = factor(Source),
    Year = factor(Year),
    Month = factor(Month),
    year_month = interaction(Year, Month, drop = TRUE),
    log_weight = log(Weight_kg),
    log_length = log(Length_cm)
  )

# -----------------------------
# 2.2 Exploratory plots
# -----------------------------

# length weight relationship by sex
p_lw <- ggplot(lwa, aes(Length_cm, Weight_kg, colour = Sex)) +
  geom_point(alpha = 0.08, size = 0.5) +
  geom_smooth(method = "lm", formula = y ~ x, se = TRUE) +
  scale_x_log10() +
  scale_y_log10() +
  labs(
    x = "Length (cm)",
    y = "Weight (kg)",
    title = "Length-weight relationship by sex"
  ) +
  theme_bw()
p_lw

# length weight by age
p_age <- ggplot(lwa, aes(Length_cm, Weight_kg, colour = Age_yrs)) +
  geom_point(alpha = 0.08, size = 0.5) +
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE) +
#   scale_x_log10() +
#   scale_y_log10() +
  labs(
    x = "Length (cm)",
    y = "Weight (kg)",
    colour = "Age",
    title = "Length-weight relationship by age"
  ) +
  theme_bw()
p_age

# length weight by month
p_month <- ggplot(lwa, aes(Length_cm, Weight_kg, colour = Month)) +
  geom_point(alpha = 0.08, size = 0.5) +
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE) +
  scale_x_log10() +
  scale_y_log10() +
  labs(
    x = "Length (cm)",
    y = "Weight (kg)",
    colour = "Month",
    title = "Length-weight relationship by month"
  ) +
  theme_bw()
p_month

# length weight by month in separate panels
p_month_facet <- ggplot(lwa, aes(Length_cm, Weight_kg)) +
  geom_point(alpha = 0.08, size = 0.5) +
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE) +
  scale_x_log10() +
  scale_y_log10() +
  facet_wrap(~ Month) +
  labs(
    x = "Length (cm)",
    y = "Weight (kg)",
    title = "Length-weight relationship by month"
  ) +
  theme_bw()
p_month_facet

# length weight by month in separate panels with modelled line of overall length weight relationship as a dashed line
p_month_facet2 <- ggplot(lwa, aes(Length_cm, Weight_kg)) +
  geom_point(alpha = 0.08, size = 0.5) +
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE) +
  scale_x_log10() +
  scale_y_log10() +
  facet_wrap(~ Month) +
  geom_smooth(data = lwa %>% select(-Month), aes(Length_cm, Weight_kg), inherit.aes = FALSE, method = "lm", formula = y ~ x, se = FALSE, colour = "black", linetype = "dashed") +
  labs(
    x = "Length (cm)",
    y = "Weight (kg)",
    title = "Length-weight relationship by month"
  ) +
  theme_bw()
p_month_facet2

# rm(p_lw,p_age,p_month,p_month_facet,p_month_facet2)

# -----------------------------
# 2.3 Mixed effects length-weight model
# -----------------------------
# log(W) = log(a) + b log(L)

# modelling length weight with month and sex as fixed effects with interaction terms
m_lw1 <- glmmTMB(log_weight ~ log_length * Month * Sex, 
  data = lwa,
  family = gaussian(),
  control = glmmTMBControl(
    optCtrl = list(iter.max = 1e4, eval.max = 1e4)
  ))

# modelling length weight with just sex as a fixed effect
m_lw2 <- glmmTMB(log_weight ~ log_length * Sex, 
  data = lwa,
  family = gaussian(),
  control = glmmTMBControl(
    optCtrl = list(iter.max = 1e4, eval.max = 1e4)
  )) 

# modelling length weight with just sex as a fixed effect and month as a random effect
m_lw3 <- glmmTMB(log_weight ~ log_length * Sex + (1 | Month), 
  data = lwa,
  family = gaussian(),
  control = glmmTMBControl(
    optCtrl = list(iter.max = 1e4, eval.max = 1e4)
  )) 

# compare models with AIC and anova
AIC(m_lw1, m_lw2, m_lw3)
anova(m_lw1, m_lw2, m_lw3)

# -----------------------------
# 2.4 Predict expected weight and calculate Le Cren's condition
# -----------------------------

##### 2.4.1 Selected model includes month
# selected model - includes month
m_lw_selected <- m_lw1 

# Le Cren's relative condition factor:
# Kn = observed weight / expected weight

lwa_cond <- lwa %>%
  mutate(
    log_weight_hat = predict(m_lw_selected, type = "response"),
    weight_hat_kg = exp(log_weight_hat),
    Kn = Weight_kg / weight_hat_kg,
    log_Kn = log(Kn)
  )

summary(lwa_cond$Kn)

# histogram of body condition
p_kn_hist <- ggplot(lwa_cond, aes(Kn)) +
  geom_histogram(bins = 80) +
  geom_vline(xintercept = 1, linetype = "dashed") +
  coord_cartesian(xlim = quantile(lwa_cond$Kn, c(0.01, 0.99), na.rm = TRUE)) +
  labs(
    x = "Le Cren's relative condition factor, Kn",
    y = "Count",
    title = "Distribution of body condition"
  ) +
  theme_bw()
p_kn_hist

# mean body condition over time by month
p_kn_time <- lwa_cond %>%
  mutate(
    Year_num = as.numeric(as.character(Year)),
    Month_num = as.numeric(as.character(Month)),
    Month_name = lubridate::month(Month_num, label = TRUE)
  ) %>%
  group_by(Year_num, Month_name) %>%
  summarise(
    mean_Kn = mean(Kn, na.rm = TRUE),
    median_Kn = median(Kn, na.rm = TRUE),
    n = n(),
    .groups = "drop"
  ) %>%
  ggplot(aes(Year_num, mean_Kn)) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  geom_point(aes(size = n), alpha = 0.5) +
  geom_smooth(se = TRUE) +
  facet_wrap(~ Month_name, scales = "free_y") +
  labs(
    x = "Year",
    y = "Mean Kn",
    size = "n",
    title = "Mean body condition through time"
  ) +
  theme_bw()
p_kn_time

# mean body condition over time by month linear line fit
p_kn_time_linear <- lwa_cond %>%
  mutate(
    Year_num = as.numeric(as.character(Year)),
    Month_num = as.numeric(as.character(Month)),
    Month_name = lubridate::month(Month_num, label = TRUE)
  ) %>%
  group_by(Year_num, Month_name) %>%
  summarise(
    mean_Kn = mean(Kn, na.rm = TRUE),
    median_Kn = median(Kn, na.rm = TRUE),
    n = n(),
    .groups = "drop"
  ) %>%
  ggplot(aes(Year_num, mean_Kn)) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  geom_point(aes(size = n), alpha = 0.5) +
  geom_smooth(method = "lm", se = TRUE) +
  facet_wrap(~ Month_name, scales = "free_y") +
  labs(
    x = "Year",
    y = "Mean Kn",
    size = "n",
    title = "Mean body condition through time"
  ) +
  theme_bw()
p_kn_time_linear


# median body condition over time by month
p_kn_time_med <- lwa_cond %>%
  mutate(
    Year_num = as.numeric(as.character(Year)),
    Month_num = as.numeric(as.character(Month))
  ) %>%
  group_by(Year_num,Month_num) %>%
  summarise(
    mean_Kn = mean(Kn, na.rm = TRUE),
    median_Kn = median(Kn, na.rm = TRUE),
    n = n(),
    .groups = "drop"
  ) %>%
  ggplot(aes(Year_num, median_Kn)) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  geom_point(aes(size = n), alpha = 0.5) +
  geom_smooth(se = TRUE) +
  facet_wrap(~ Month_num, scales = "free_y") +
  labs(
    x = "Year",
    y = "Median Kn",
    size = "n",
    title = "Median body condition through time"
  ) +
  theme_bw()
p_kn_time_med

# distribution historgram of body condition
p_kn_hist

# mean body condition over time by month
p_kn_time
# p_kn_time_med # median condition over time by month

# save loess plot
ggsave(filename = "figures/mean_body_condition_through_time_modelmonth.png", plot = p_kn_time, width = 8, height = 6, dpi = 300)

# save linear plot
ggsave(filename = "figures/mean_body_condition_through_time_linear_modelmonth.png", plot = p_kn_time_linear, width = 8, height = 6, dpi = 300)

rm(p_kn_hist, p_kn_time, p_kn_time_med, p_kn_time_linear)
gc()

##### 2.4.2 Selected model excludes month
# selected model  - excludes month
m_lw_selected <- m_lw2

# Le Cren's relative condition factor:
# Kn = observed weight / expected weight

lwa_cond <- lwa %>%
  mutate(
    log_weight_hat = predict(m_lw_selected, type = "response"),
    weight_hat_kg = exp(log_weight_hat),
    Kn = Weight_kg / weight_hat_kg,
    log_Kn = log(Kn)
  )

summary(lwa_cond$Kn)

# histogram of body condition
p_kn_hist <- ggplot(lwa_cond, aes(Kn)) +
  geom_histogram(bins = 80) +
  geom_vline(xintercept = 1, linetype = "dashed") +
  coord_cartesian(xlim = quantile(lwa_cond$Kn, c(0.01, 0.99), na.rm = TRUE)) +
  labs(
    x = "Le Cren's relative condition factor, Kn",
    y = "Count",
    title = "Distribution of body condition"
  ) +
  theme_bw()
p_kn_hist

# mean body condition over time by month
p_kn_time <- lwa_cond %>%
  mutate(
    Year_num = as.numeric(as.character(Year)),
    Month_num = as.numeric(as.character(Month)),
    Month_name = lubridate::month(Month_num, label = TRUE)
  ) %>%
  group_by(Year_num, Month_name) %>%
  summarise(
    mean_Kn = mean(Kn, na.rm = TRUE),
    median_Kn = median(Kn, na.rm = TRUE),
    n = n(),
    .groups = "drop"
  ) %>%
  ggplot(aes(Year_num, mean_Kn)) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  geom_point(aes(size = n), alpha = 0.5) +
  geom_smooth(se = TRUE) +
  facet_wrap(~ Month_name, scales = "free_y") +
  labs(
    x = "Year",
    y = "Mean Kn",
    size = "n",
    title = "Mean body condition through time"
  ) +
  theme_bw()
p_kn_time

# mean body condition over time by month linear line fit
p_kn_time_linear <- lwa_cond %>%
  mutate(
    Year_num = as.numeric(as.character(Year)),
    Month_num = as.numeric(as.character(Month)),
    Month_name = lubridate::month(Month_num, label = TRUE)
  ) %>%
  group_by(Year_num, Month_name) %>%
  summarise(
    mean_Kn = mean(Kn, na.rm = TRUE),
    median_Kn = median(Kn, na.rm = TRUE),
    n = n(),
    .groups = "drop"
  ) %>%
  ggplot(aes(Year_num, mean_Kn)) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  geom_point(aes(size = n), alpha = 0.5) +
  geom_smooth(method = "lm", se = TRUE) +
  facet_wrap(~ Month_name, scales = "free_y") +
  labs(
    x = "Year",
    y = "Mean Kn",
    size = "n",
    title = "Mean body condition through time"
  ) +
  theme_bw()
p_kn_time_linear

# median body condition over time by month
p_kn_time_med <- lwa_cond %>%
  mutate(
    Year_num = as.numeric(as.character(Year)),
    Month_num = as.numeric(as.character(Month))
  ) %>%
  group_by(Year_num,Month_num) %>%
  summarise(
    mean_Kn = mean(Kn, na.rm = TRUE),
    median_Kn = median(Kn, na.rm = TRUE),
    n = n(),
    .groups = "drop"
  ) %>%
  ggplot(aes(Year_num, median_Kn)) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  geom_point(aes(size = n), alpha = 0.5) +
  geom_smooth(se = TRUE) +
  facet_wrap(~ Month_num, scales = "free_y") +
  labs(
    x = "Year",
    y = "Median Kn",
    size = "n",
    title = "Median body condition through time"
  ) +
  theme_bw()
p_kn_time_med

# distribution historgram of body condition
p_kn_hist

# mean body condition over time by month
p_kn_time
# p_kn_time_med # median condition over time by month

# save plot
ggsave(filename = "figures/mean_body_condition_through_time_modelnomonth.png", plot = p_kn_time, width = 8, height = 6, dpi = 300)

# save linear plot
ggsave(filename = "figures/mean_body_condition_through_time_linear_modelnomonth.png", plot = p_kn_time_linear, width = 8, height = 6, dpi = 300)

rm(p_kn_hist, p_kn_time, p_kn_time_med, p_kn_time_linear)
gc()

# -----------------------------
# 2.5 Prepare monthly condition index
# -----------------------------

# use selected model - m_lw1 = month included, m_lw2 = month excluded
m_lw_selected <- m_lw1
m_lw_selected <- m_lw2

# predict the expected weight and estimate condition
lwa_cond <- lwa %>%
  mutate(
    log_weight_hat = predict(m_lw_selected, type = "response"),
    weight_hat_kg = exp(log_weight_hat),
    Kn = Weight_kg / weight_hat_kg,
    log_Kn = log(Kn),
    year = as.integer(as.character(Year)),
    month = as.integer(as.character(Month))
  )

# summarize condition by month
cond_monthly <- lwa_cond %>%
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

# summarize condition by year
cond_annual <- lwa_cond %>%
  group_by(year) %>%
  summarise(
    mean_age = mean(as.numeric(as.character(Age_yrs)), na.rm = TRUE),
    mean_Kn = mean(Kn, na.rm = TRUE),
    median_Kn = median(Kn, na.rm = TRUE),
    mean_log_Kn = mean(log_Kn, na.rm = TRUE),
    sd_log_Kn = sd(log_Kn, na.rm = TRUE),
    se_log_Kn = sd_log_Kn / sqrt(n()),
    n_fish = n(),
    .groups = "drop"
  )

# mean monthly condition plot
cond_monthly_plot <- cond_monthly %>%
  mutate(date = make_date(year, month, 1))

p_cond_monthly <- ggplot(cond_monthly_plot,
             aes(x = date, y = mean_Kn)) +
  geom_hline(yintercept = 1, linetype = "dashed",
             colour = "grey50") +
  geom_point(alpha = 0.4, size = 1.5,
             colour = "steelblue") +
  geom_smooth(
    method = "lm",
    span = 0.3,
    se = TRUE,
    colour = "black",
    fill = "grey70",
    linewidth = 1
  ) +
  scale_x_date(date_breaks = "5 years",
               date_labels = "%Y") +
  labs(
    x = "Year",
    y = "Mean relative condition (Kn)",
    title = "Monthly mean body condition over time"
  ) +
  theme_classic(base_size = 13)

p_cond_monthly
ggsave(filename = "figures/mean_monthly_body_condition_over_time.png", plot = p_cond_monthly, width = 8, height = 6, dpi = 300)

# mean by year, across months
p_cond_monthly2 <- cond_monthly %>%
  mutate(
    month = factor(
      month,
      levels = 1:12,
      labels = month.abb
    )
  ) %>%
  ggplot(aes(x = month, y = mean_Kn)) +
  geom_hline(yintercept = 1, linetype = "dashed",
             colour = "grey50") +
  geom_boxplot(
    fill = "steelblue",
    alpha = 0.5,
    outlier.alpha = 0.3
  ) +
  geom_jitter(
    width = 0.15,
    alpha = 0.3,
    size = 1
  ) +
  labs(
    x = "Month",
    y = "Mean relative condition (Kn)",
    title = "Seasonal variation in body condition"
  ) +
  theme_classic(base_size = 13)

p_cond_monthly2
ggsave(filename = "figures/seasonal_variation_body_condition.png", plot = p_cond_monthly2, width = 8, height = 6, dpi = 300)  

# separating years and months
p_cond_monthly3 <- cond_monthly %>%
  mutate(month = factor(
    month,
    levels = 1:12,
    labels = month.abb
  )) %>%
  ggplot(aes(x = year, y = mean_Kn)) +
  geom_hline(yintercept = 1, linetype = "dashed",
             colour = "grey60") +
  geom_point(alpha = 0.5, size = 1) +
  geom_smooth(
    method = "loess",
    se = TRUE,
    colour = "steelblue"
  ) +
  facet_wrap(~ month, scales = "free_y", ncol = 3) +
  labs(
    x = "Year",
    y = "Mean relative condition (Kn)",
    title = "Temporal trends in body condition by month"
  ) +
  theme_classic(base_size = 12)

p_cond_monthly3
ggsave(filename = "figures/temporal_trends_body_condition_by_month.png", plot = p_cond_monthly3, width = 8, height = 6, dpi = 300)

# plotting annual mean condition over time
p_cond_annual <- ggplot(cond_annual, aes(x = year, y = mean_Kn)) +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "grey50") +
  geom_point(alpha = 0.5, size = 1.5, colour = "steelblue") +
  geom_smooth(method = "lm", span = 0.3, se = TRUE, colour = "black", fill = "grey70", linewidth = 1) +
  scale_x_continuous(breaks = seq(min(cond_annual$year), max(cond_annual$year), by = 5)) +
  labs(
    x = "Year",
    y = "Mean relative condition (Kn)",
    title = "Annual mean body condition over time"
  ) +
  theme_classic(base_size = 13)
p_cond_annual
ggsave(filename = "figures/annual_mean_body_condition_over_time.png", plot = p_cond


##### -------------------------------------------------------------------- #####
##### -------------------------------------------------------------------- #####
##### -------------------------------------------------------------------- #####
##### -------------------------------------------------------------------- #####
##### -------------------------------------------------------------------- #####

##### Bayesian analysis for length-weight relationship

# -----------------------------
# 2.11 Bayesian analysis packages
# -----------------------------

library(brms)
library(cmdstanr)
library(posterior)
library(tidybayes)

# parallel processing for cmdstanr
options(mc.cores = parallel::detectCores())

# Optional: check that cmdstan is available
cmdstanr::cmdstan_version()

# -----------------------------
# 2.12 Prepare data for Bayesian length-weight model
# -----------------------------

# redo variables for bayesian model
lwa_bayes <- lwa %>%
  mutate(
    year = as.integer(as.character(Year)),
    month = as.integer(as.character(Month)),

    # Explicit factor levels help when predicting
    Sex = droplevels(factor(Sex)),
    Source = droplevels(factor(Source)),
    Month_f = factor(
      month,
      levels = 1:12,
      labels = month.abb
    ),

    # Center log length to improve MCMC sampling.
    # Predictions are still returned on the original log-weight scale.
    log_length_mean = mean(log_length),
    log_length_c = log_length - log_length_mean,

    # Unique observation ID for joining predictions back to fish
    fish_id = row_number()
  )

log_length_center <- unique(lwa_bayes$log_length_mean)

stopifnot(length(log_length_center) == 1)

# -----------------------------
# 2.13 Priors
# -----------------------------

prior_intercept_mean <- median(lwa_bayes$log_weight)

priors_lw <- c(
  # Expected log weight at the centered length
  set_prior(
    paste0(
      "normal(",
      format(
        prior_intercept_mean,
        scientific = FALSE,
        digits = 15
      ),
      ", 1)"
    ),
    class = "Intercept"
  ),

  # Fixed effects, including length slope, sex effect,
  # month effects and length-by-sex interaction
  prior(
    normal(0, 1),
    class = "b"
  ),

  # Among-source variability - not needed for simpler LW model
  # prior(
  #   student_t(3, 0, 0.5),
  #   class = "sd"
  # ),

  # Individual residual variability in log weight
  prior(
    student_t(3, 0, 0.5),
    class = "sigma"
  )
)

priors_lw

get_prior(
  formula =
    log_weight ~ log_length_c * Sex,
  data = lwa_bayes,
  family = gaussian()
)

# -----------------------------
# 2.14 Bayesian length-weight model
# -----------------------------
#
# Deliberately excludes a Year random intercept so that annual
# weight-at-length variation remains in the condition index.
#
# Month controls for the average seasonal pattern.
# Month as a fixed or random effect takes a long time to run. Maybe consider alternative approach

form_lw_bayes <- bf(
  log_weight ~
    log_length_c * Sex
    # log_length_c * Sex * Month
    # log_length_c * Sex + (1 + log_length_c | Month)
)

m_lw_bayes <- brm(
  formula = form_lw_bayes,
  data = lwa_bayes,
  family = gaussian(),
  prior = priors_lw,

  backend = "cmdstanr",
  chains = 4,
  iter = 3000,
  warmup = 1500,
  seed = 2,

  control = list(
    adapt_delta = 0.95,
    max_treedepth = 12
  ),

  save_pars = save_pars(all = TRUE),
  file = "model-outputs/m_lw_bayes_sexmonth",
  file_refit = "on_change"
)

# -----------------------------
# 2.15 Model diagnostics
# -----------------------------

# show r-hats, ESS and coefficient estimates distributions
summary(m_lw_bayes)

plot(m_lw_bayes)

# These take a while to run so skipping these
# show density distribution of posterior predictions with observed data
# pp_check(
#   m_lw_bayes,
#   ndraws = 100,
#   type = "dens_overlay"
# )

# compare observed vs predicted as scatter plot
# pp_check(
#   m_lw_bayes,
#   ndraws = 100,
#   type = "scatter_avg"
# )

# variation explained by model
r2_lw <- bayes_R2(
  m_lw_bayes,
  ndraws = 500
)
r2_lw

# Fixed-effect posterior summaries
fixef(m_lw_bayes, probs = c(0.025, 0.5, 0.975))

# -----------------------------
# Sex-specific length slopes
# -----------------------------

lw_coef_draws <- as_draws_df(m_lw_bayes)

sex_levels <- levels(lwa_bayes$Sex)

# female are the reference
sex_levels

variables(lw_coef_draws)[
  stringr::str_detect(
    variables(lw_coef_draws),
    "log_length"
  )
]

# getting slope estimates (and distirbution) for each sex from bayesian draws
slope_draws <- lw_coef_draws %>%
  transmute(
    .draw,
    slope_F = b_log_length_c,
    slope_M = b_log_length_c + `b_log_length_c:SexM`
  )

# summarizing them
slope_summary <- slope_draws %>%
  pivot_longer(
    cols = starts_with("slope_"),
    names_to = "Sex",
    values_to = "b"
  ) %>%
  group_by(Sex) %>%
  summarise(
    median = median(b),
    lower = quantile(b, 0.025),
    upper = quantile(b, 0.975),
    .groups = "drop"
  )

# print slope summary
# female slope is greater than male slope
slope_summary

# -----------------------------
# 2.16 Construct posterior draws of monthly condition
# -----------------------------

# function to extract posterior distribution of monthly condition
make_monthly_condition_draws <- function(
    model,
    data,
    ndraws = 500,
    chunk_size = 10000,
    seed = 1837,
    include_source_effect = TRUE
) {

  set.seed(seed)

  # Select posterior draws once and use the same draw IDs
  # for every data chunk.
  total_draws <- posterior::ndraws(model)

  selected_draw_ids <- sort(
    sample(
      seq_len(total_draws),
      size = min(ndraws, total_draws),
      replace = FALSE
    )
  )

  # Integer index for year-month groups
  group_lookup <- data %>%
    distinct(year, month) %>%
    arrange(year, month) %>%
    mutate(group_id = row_number())

  data_indexed <- data %>%
    left_join(
      group_lookup,
      by = c("year", "month")
    )

  n_groups <- nrow(group_lookup)
  n_selected <- length(selected_draw_ids)

  # Accumulate:
  # 1. sum of observed log weights
  # 2. sum of posterior expected log weights
  # 3. number of fish
  observed_sum <- numeric(n_groups)
  expected_sum <- matrix(
    0,
    nrow = n_selected,
    ncol = n_groups
  )
  group_n <- integer(n_groups)

  chunk_starts <- seq(
    from = 1,
    to = nrow(data_indexed),
    by = chunk_size
  )

  for (chunk_number in seq_along(chunk_starts)) {

    row_start <- chunk_starts[chunk_number]
    row_end <- min(
      row_start + chunk_size - 1,
      nrow(data_indexed)
    )

    rows <- row_start:row_end
    chunk <- data_indexed[rows, ]

    message(
      "Processing chunk ",
      chunk_number,
      " of ",
      length(chunk_starts),
      ": rows ",
      row_start,
      "-",
      row_end
    )

    # Posterior draws of the expected log weight.
    #
    # re_formula = NULL:
    # include estimated Source random effects.
    #
    # re_formula = NA:
    # exclude all random effects and estimate condition relative
    # to the population-average source.
    mu_chunk <- posterior_linpred(
      model,
      newdata = chunk,
      transform = FALSE,
      re_formula = if (include_source_effect) NULL else NA,
      draw_ids = selected_draw_ids
    )

    group_id <- chunk$group_id

    # Add observed log weights and sample sizes by group
    observed_chunk <- rowsum(
      chunk$log_weight,
      group = group_id,
      reorder = FALSE
    )

    n_chunk <- rowsum(
      rep(1L, nrow(chunk)),
      group = group_id,
      reorder = FALSE
    )

    group_present <- as.integer(rownames(observed_chunk))

    observed_sum[group_present] <-
      observed_sum[group_present] +
      observed_chunk[, 1]

    group_n[group_present] <-
      group_n[group_present] +
      n_chunk[, 1]

    # Aggregate expected log weights within this chunk.
    #
    # rowsum() operates on rows, so transpose first:
    # fish x posterior draws.
    expected_chunk <- rowsum(
      t(mu_chunk),
      group = group_id,
      reorder = FALSE
    )

    expected_sum[, group_present] <-
      expected_sum[, group_present] +
      t(expected_chunk)
  }

  # Monthly mean log(Kn):
  #
  # mean(log(W observed) - log(W expected))
  #
  # = mean(log observed weight) -
  #   mean(expected log weight)
  mean_observed_log_weight <- observed_sum / group_n

  mean_expected_log_weight <- sweep(
    expected_sum,
    MARGIN = 2,
    STATS = group_n,
    FUN = "/"
  )

  mean_log_kn <- sweep(
    mean_expected_log_weight,
    MARGIN = 2,
    STATS = mean_observed_log_weight,
    FUN = function(expected, observed) observed - expected
  )

  # Convert posterior matrix to tidy format
  condition_draws <- as_tibble(mean_log_kn) %>%
    mutate(
      .draw_index = seq_len(n_selected),
      .draw = selected_draw_ids
    ) %>%
    pivot_longer(
      cols = starts_with("V"),
      names_to = "group_column",
      values_to = "mean_log_Kn"
    ) %>%
    mutate(
      group_id = as.integer(
        stringr::str_remove(group_column, "^V")
      )
    ) %>%
    select(
      .draw,
      .draw_index,
      group_id,
      mean_log_Kn
    ) %>%
    left_join(
      group_lookup,
      by = "group_id"
    ) %>%
    mutate(
      # Geometric-mean Kn for each posterior draw
      geometric_mean_Kn = exp(mean_log_Kn)
    ) %>%
    arrange(.draw_index, year, month)

  condition_draws
}

# calculate monthly condition draws from bayesian model
cond_monthly_draws <- make_monthly_condition_draws(
  model = m_lw_bayes,
  data = lwa_bayes,
  ndraws = 500,
  chunk_size = 10000,
  seed = 22,

  # TRUE adjusts expected weight for systematic source differences.
  include_source_effect = TRUE
)


# -----------------------------
# 2.17 Posterior summaries by year and month
# -----------------------------

cond_monthly_summary <- cond_monthly_draws %>%
  group_by(year, month) %>%
  summarise(
    median_log_Kn = median(mean_log_Kn),
    lower_log_Kn = quantile(mean_log_Kn, 0.025),
    upper_log_Kn = quantile(mean_log_Kn, 0.975),

    median_Kn = median(geometric_mean_Kn),
    lower_Kn = quantile(geometric_mean_Kn, 0.025),
    upper_Kn = quantile(geometric_mean_Kn, 0.975),

    probability_Kn_gt_1 = mean(geometric_mean_Kn > 1),

    .groups = "drop"
  ) %>%
  left_join(
    lwa_bayes %>%
      count(year, month, name = "n_fish"),
    by = c("year", "month")
  )

cond_monthly_summary


# plot monthly posterior condition index
ggplot(
  cond_monthly_summary,
  aes(
    x = as.Date(sprintf("%04d-%02d-15", year, month)),
    y = median_Kn
  )
) +
  geom_hline(
    yintercept = 1,
    linetype = "dashed"
  ) +
  geom_ribbon(
    aes(
      ymin = lower_Kn,
      ymax = upper_Kn
    ),
    alpha = 0.25
  ) +
  geom_line() +
  geom_point(
    aes(size = n_fish),
    alpha = 0.6
  ) +
  labs(
    x = NULL,
    y = "Posterior geometric mean Kn",
    size = "Number of fish",
    title = "Monthly posterior body-condition index",
    subtitle = "Ribbon shows the 95% credible interval"
  ) +
  theme_bw()

# seasonal patterns
ggplot(
  cond_monthly_summary,
  aes(
    x = year,
    y = median_Kn
  )
) +
  geom_hline(
    yintercept = 1,
    linetype = "dashed"
  ) +
  geom_ribbon(
    aes(
      ymin = lower_Kn,
      ymax = upper_Kn
    ),
    alpha = 0.25
  ) +
  geom_line() +
  facet_wrap(
    ~ month,
    ncol = 4
  ) +
  labs(
    x = "Year",
    y = "Posterior geometric mean Kn",
    title = "Body-condition anomalies by calendar month"
  ) +
  theme_bw()


# -----------------------------
# 2.18 Fish-level condition draws for validation
# -----------------------------

set.seed(22)

# sample from observed data used for bayes analysis
validation_fish <- lwa_bayes %>%
  slice_sample(n = min(2000, nrow(lwa_bayes)))

# posterior draws of expected log_weight
mu_validation <- posterior_linpred(
  m_lw_bayes,
  newdata = validation_fish,
  transform = FALSE,
  re_formula = NULL,
  ndraws = 500
)

# mutate the posterior draws into a dataframe with original observed data. calculated distribution of expected condition for each fish
fish_condition_draws <- as_tibble(mu_validation) %>%
  mutate(.draw = row_number()) %>%
  pivot_longer(
    cols = starts_with("V"),
    names_to = "fish_column",
    values_to = "expected_log_weight"
  ) %>%
  mutate(
    fish_row = as.integer(
      stringr::str_remove(fish_column, "^V")
    )
  ) %>%
  left_join(
    validation_fish %>%
      mutate(fish_row = row_number()) %>%
      select(
        fish_row,
        fish_id,
        year,
        month,
        Sex,
        Source,
        Length_cm,
        Weight_kg,
        log_weight
      ),
    by = "fish_row"
  ) %>%
  mutate(
    log_Kn = log_weight - expected_log_weight,
    Kn = exp(log_Kn)
  )

fish_condition_draws

# summarize each fish's posterior condition
fish_condition_summary <- fish_condition_draws %>%
  group_by(
    fish_row,
    fish_id,
    year,
    month,
    Sex,
    Source,
    Length_cm,
    Weight_kg
  ) %>%
  summarise(
    Kn = median(Kn),
    Kn_lower = quantile(Kn, 0.025),
    Kn_upper = quantile(Kn, 0.975),
    log_Kn = median(log_Kn),
    .groups = "drop"
  ) %>%
  arrange(year, month) %>%
  mutate(date = lubridate::make_date(year = year, month = month, day = 1))

# plot fish size and median body condition from posterior distribution
ggplot(
  fish_condition_summary,
  aes(x = Length_cm, y = Kn)
) +
  geom_point(alpha = 0.2) +
  geom_smooth() +
  geom_hline(yintercept = 1, linetype = 2)

# plot fish size and median body condition from posterior distribution - separated by sex
ggplot(
  fish_condition_summary,
  aes(x = Length_cm, y = Kn, colour = Sex)
) +
  geom_point(alpha = 0.15) +
  geom_smooth(
    method = "gam",
    formula = y ~ s(x),
    se = TRUE
  ) +
  geom_hline(
    yintercept = 1,
    linetype = "dashed"
  ) +
  labs(
    x = "Fork length (cm)",
    y = "Relative condition (Kn)",
    colour = "Sex"
  ) +
  theme_classic()

# log condition
ggplot(
  fish_condition_summary,
  aes(x = Length_cm, y = log_Kn, colour = Sex)
) +
  geom_point(alpha = 0.15) +
  geom_smooth(
    method = "gam",
    formula = y ~ s(x),
    se = TRUE
  ) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  theme_classic()

# check effects of sex on overall body condition
ggplot(
  fish_condition_summary,
  aes(x = Sex, y = Kn)
) +
  geom_boxplot() +
  geom_hline(yintercept = 1, linetype = 2)

# check effects of month
ggplot(
  fish_condition_summary,
  aes(x = factor(month), y = Kn)
) +
  geom_boxplot() +
  geom_hline(yintercept = 1, linetype = 2)

ggplot(
  fish_condition_summary,
  aes(x = Length_cm, y = Kn)
) +
  geom_point(alpha = 0.1) +
  geom_smooth(
    method = "gam",
    formula = y ~ s(x),
    se = FALSE
  ) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  facet_wrap(~ month) +
  labs(
    x = "Fork length (cm)",
    y = "Relative condition (Kn)"
  ) +
  theme_classic()

# check effects of sex and condition over time
ggplot(
  fish_condition_summary,
  aes(x = date, y = Kn, col = Sex)
) +
  geom_point(alpha = 0.15) +
  geom_smooth(
    method = "gam",
    formula = y ~ s(x),
    se = TRUE
  )

# check effects of source
ggplot(
  fish_condition_summary,
  aes(x = Source, y = Kn)
) +
  geom_boxplot() +
  geom_hline(yintercept = 1, linetype = 2) +
  coord_flip()

ggplot(
  fish_condition_summary,
  aes(x = Length_cm, y = Kn)
) +
  geom_point(alpha = 0.1) +
  geom_smooth(
    method = "gam",
    formula = y ~ s(x),
    se = FALSE
  ) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  facet_wrap(~ Source) +
  labs(
    x = "Fork length (cm)",
    y = "Relative condition (Kn)"
  ) +
  theme_classic()

