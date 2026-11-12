#####
#
# Hake length weight analysis and body condition
#
#


# -----------------------------
# 1.1 Read in data
# -----------------------------
#####


# US Hake weight-at-age data
# only has US at sea and US shore Source data
# has data from 1965 to 2025
us_dat <- read.csv("https://raw.githubusercontent.com/pacific-hake/hake-assessment/refs/heads/master/data-tables/us-weight-at-age.csv", header = TRUE)

write.csv(us_dat, "data/github_data/us-weight-at-age.csv", row.names = FALSE)

str(us_dat)
head(us_dat)
unique(us_dat$Source)
unique(us_dat$Year)


# length weight age data
# Has numerous sources
# has data from 1975 to 2021
lwa_dat <- read.csv("https://raw.githubusercontent.com/pacific-hake/hake-assessment/refs/heads/master/data-tables/length-weight-age/LWAdata_1975to2021.csv", header = TRUE)

write.csv(lwa_dat, "data/github_data/LWAdata_1975to2021.csv", row.names = FALSE)

str(lwa_dat)
head(lwa_dat)
unique(lwa_dat$Source)
unique(lwa_dat$Year)


# acoustic survey weight at age
# US and Can Acoustic sources
# 1995 to 2025, every 2-3 years
ac_dat <- read.csv("https://raw.githubusercontent.com/pacific-hake/hake-assessment/refs/heads/master/data-tables/survey-weight-at-age.csv", header = TRUE)

write.csv(ac_dat, "data/github_data/survey-weight-at-age.csv", row.names = FALSE)

str(ac_dat)
head(ac_dat)
unique(ac_dat$Source)
unique(ac_dat$Year)

# Biomass at age
#  Age is in columns, and the values represent the total biomass at each age
# These data are from Andy from Rebecca
# github directory

# # canada survey data
# can.dat <- read.csv(paste0("data/survey_biomass_age_proportions/", "Survey_Biomass_at_age_Canada.csv"), header = TRUE)

# # california survey data
# cal.dat <- read.csv(paste0("data/survey_biomass_age_proportions/", "Survey_Biomass_at_age_California.csv"), header = TRUE)

# # oregon survey data
# ore.dat <- read.csv(paste0("data/survey_biomass_age_proportions/", "Survey_Biomass_at_age_Oregon.csv"), header = TRUE)

# # total survey data
# tot.dat <- read.csv(paste0("data/survey_biomass_age_proportions/", "Survey_Biomass_at_Age_Total.csv"), header = TRUE)

# str(can.dat)
# str(cal.dat)
# str(ore.dat)
# str(tot.dat)

#####
# Checkout other data on github data/ folder

# # weight age sample sizes
# data_url <- paste0(
#     "https://github.com/pacific-hake/hake-assessment/blob/master/data/weight_age_sample_sizes_df.rda",
#     "?raw=true"
# )
# download.file(data_url, destfile = "data/github_data/weight_age_sample_sizes_df.rda")
# load("data/github_data/weight_age_sample_sizes_df.rda")
# weight_age_sample_sizes_df

# # weight at age df
# data_url <- paste0(
#     "https://github.com/pacific-hake/hake-assessment/blob/master/data/weight_at_age_df.rda",
#     "?raw=true"
# )
# download.file(data_url, destfile = "data/github_data/weight_at_age_df.rda")
# load("data/github_data/weight_at_age_df.rda")
# weight_at_age_df
# head(weight_at_age_df)

# # weight at age fn
# data_url <- paste0(
#     "https://github.com/pacific-hake/hake-assessment/blob/master/data/weight_at_age_fn.rda",
#     "?raw=true"
# )
# download.file(data_url, destfile = "data/github_data/weight_at_age_fn.rda")
# load("data/github_data/weight_at_age_fn.rda")
# head(weight_at_age_fn)

# # can ft bottom depth df
# data_url <- paste0(
#     "https://github.com/pacific-hake/hake-assessment/blob/master/data/can_ft_bottom_depth_df.rda",
#     "?raw=true"
# )
# download.file(data_url, destfile = "data/github_data/can_ft_bottom_depth_df.rda")
# load("data/github_data/can_ft_bottom_depth_df.rda")
# head(can_ft_bottom_depth_df)

# # can_ft_catch_by_month_df
# data_url <- paste0(
#     "https://github.com/pacific-hake/hake-assessment/blob/master/data/can_ft_catch_by_month_df.rda",
#     "?raw=true"
# )
# download.file(data_url, destfile = "data/github_data/can_ft_catch_by_month_df.rda")
# load("data/github_data/can_ft_catch_by_month_df.rda")
# head(can_ft_catch_by_month_df)


# # can_waa_fn
# data_url <- paste0(
#     "https://github.com/pacific-hake/hake-assessment/blob/master/data/can_waa_fn.rda",
#     "?raw=true"
# )
# download.file(data_url, destfile = "data/github_data/can_waa_fn.rda")
# load("data/github_data/can_waa_fn.rda")
# head(can_waa_fn)

# # recruit_dev_start_yr
# data_url <- paste0(
#     "https://github.com/pacific-hake/hake-assessment/blob/master/data/recruit_dev_start_yr.rda",
#     "?raw=true"
# )
# download.file(data_url, destfile = "data/github_data/recruit_dev_start_yr.rda")
# load("data/github_data/recruit_dev_start_yr.rda")
# head(recruit_dev_start_yr)

# # survey_by_country_df
# data_url <- paste0(
#     "https://github.com/pacific-hake/hake-assessment/blob/master/data/survey_by_country_df.rda",
#     "?raw=true"
# )
# download.file(data_url, destfile = "data/github_data/survey_by_country_df.rda")
# load("data/github_data/survey_by_country_df.rda")
# head(survey_by_country_df)

# # us_ms_age_df
# data_url <- paste0(
#     "https://github.com/pacific-hake/hake-assessment/blob/master/data/us_ms_age_df.rda",
#     "?raw=true"
# )
# download.file(data_url, destfile = "data/github_data/us_ms_age_df.rda")
# load("data/github_data/us_ms_age_df.rda")
# head(us_ms_age_df)

# # us_ss_age_df
# data_url <- paste0(
#     "https://github.com/pacific-hake/hake-assessment/blob/master/data/us_ss_age_df.rda",
#     "?raw=true"
# )
# download.file(data_url, destfile = "data/github_data/us_ss_age_df.rda")
# load("data/github_data/us_ss_age_df.rda")
# head(us_ss_age_df)

# # us_ti_ct_by_month_df
# data_url <- paste0(
#     "https://github.com/pacific-hake/hake-assessment/blob/master/data/us_ti_ct_by_month_df.rda",
#     "?raw=true"
# )
# download.file(data_url, destfile = "data/github_data/us_ti_ct_by_month_df.rda")
# load("data/github_data/us_ti_ct_by_month_df.rda")
# head(us_ti_ct_by_month_df)

# # none of these appear to be what I want - they don't have length data
# # remove objects from line 75
# rm(
#     weight_at_age_fn, can_ft_bottom_depth_df, can_ft_catch_by_month_df, can_waa_fn, recruit_dev_start_yr, survey_by_country_df, us_ms_age_df, us_ss_age_df, us_ti_ct_by_month_df, can_waa_fn,
#     weight_at_age_df, weight_age_sample_sizes_df
# )
