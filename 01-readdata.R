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

str(us_dat)
head(us_dat)
unique(us_dat$Source)
unique(us_dat$Year)

# length weight age data
# Has numerous sources
# has data from 1975 to 2021
lwa_dat <- read.csv("https://raw.githubusercontent.com/pacific-hake/hake-assessment/refs/heads/master/data-tables/length-weight-age/LWAdata_1975to2021.csv", header = TRUE)

str(lwa_dat)
head(lwa_dat)
unique(lwa_dat$Source)
unique(lwa_dat$Year)

# acoustic survey weight at age
# US and Can Acoustic sources
# 1995 to 2025, every 2-3 years
ac_dat <- read.csv("https://raw.githubusercontent.com/pacific-hake/hake-assessment/refs/heads/master/data-tables/survey-weight-at-age.csv", header = TRUE)

str(ac_dat)
head(ac_dat)
unique(ac_dat$Source)
unique(ac_dat$Year)

## NEED TO FIX BELOW AND GET DATA TO REPO
# Biomass at age
#  Age is in columns, and the values represent the total biomass at each age
# github directory
# git_dir <- "C:/github/hakeanalysis/data/"

# # canada survey data
# can.dat <- read.csv(paste0(git_dir, "Survey_Biomass_at_age_Canada.csv"), header = TRUE)

# # california survey data
# cal.dat <- read.csv(paste0(git_dir, "Survey_Biomass_at_age_California.csv"), header = TRUE)

# # oregon survey data
# ore.dat <- read.csv(paste0(git_dir, "Survey_Biomass_at_age_Oregon.csv"), header = TRUE)

# # total survey data
# tot.dat <- read.csv(paste0(git_dir, "Survey_Biomass_at_Age_Total.csv"), header = TRUE)

# str(can.dat)
# str(cal.dat)
# str(ore.dat)
# str(tot.dat)
