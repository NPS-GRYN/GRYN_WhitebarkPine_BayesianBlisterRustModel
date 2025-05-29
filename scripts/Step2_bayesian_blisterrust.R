### Blister Rust Hierarchical Occupancy Model ###
## Lydia M. Druin ##
# GRYN WBP Bayesian Blister Rust #

# Stored as "/bayesian_blisterrust.R" #
# Updated 03/21/25 #

### Purpose ###
# This code is for completing the analysis and running the model of HOccu2,
# the Bayesian hierarchical occupancy model described in Shanahan et al. 2021.
# Raw data is accessible in the "data/raw_data" folder and must be processed
# using the script "scripts/data_processing.R" before beginning!

### Note to User ###
# If you have already completed all model running and steps in this code and want
# to examine models, page down to the "Jump-to-Point". If you have already completed
# all model running and steps in this code and want to make figures/graphs,
# open the "scripts/model_plotting.R" script instead.

# For incorporating new data, be sure the input filenames and exported filenames
# are updated with your initials/date.

## Note ##
# This code and model were written by Wilson Wright (see Shanahan et al. 2021) and
# adapted by LMD in 2024.

# Set-up ####

### Clean Up Workspace ###
rm(list = ls()) # empty environment/loaded items in workspace--not required

gc() # release memory

### Load Libraries ###

install.packages(c("tidyverse", "rstan")) # only need to do once--but will not hurt to do every time

library(tidyverse)
library(rstan)

options(scipen = 999) # trick to force R to stay in numeric (not scientific) notation
options(dplyr.width = Inf, dplyr.print_min = 100) # used to set dimensions for viewing objects

options(mc.cores = parallel::detectCores()) # configure rstan to use available cores on your machine
rstan_options(auto_write = TRUE)

# Step 1. Load in Data ####

## Part A. Load in Source File ####

# Data developed with data_processing.R #
gye_obs <- readRDS("data/rds/GYE_WBP_Processed_AllObservers_BlisterRustData_LMDDec2024_late.rds")
tran_data <- readRDS("data/rds/GYE_WBP_Processed_TransectData_LMDDec2024.rds")

# Raw data, no need for updating UNLESS a new stand/transect is added #
stand_hiking <- read.csv("data/raw_data/stand_hiking.csv")

### NEEDS ATTENTION ####
# Load model--this can take some time
final_model <- stan_model("D:/R codes/final_model1LMD.stan")
# Are you incorporating data post-2023? Read the "README_BEFORE_RUNNING_MODEL.txt"
# file in the project folder. Make necessary changes before continuing!

## Part B. Arrange Data ####
gye_obs <- arrange(gye_obs, time, standID, transectID, treeID, obs_priority)
# This sorts the data by time period, then stand, then transect, then tree,
# and then finally observer priority--which makes repeated observations for a given
# tree within the same time period appear in consecutive rows in the dataframe,
# with the observer of the highest priority listed first.

## Part C. Format Tree-level Data ####
tree_data <- gye_obs %>%
  group_by(time, standID, transectID, siteID, treeID) %>%
  summarize(
    n_visits = length(inf_ind),
    naive_inf = first(inf_ind)
  )
# This creates a tree-level dataframe with the number of visits to each tree in a
# given time period (n_visits) and a naive infection indicator (naive_inf) based
# on whether blister rust was observed by the leading priority observer during
# each time period

tree_data <- arrange(tree_data, time, standID, transectID, treeID)

# Confirm the ordering is correct:
all.equal(
  rep(tree_data$time, times = tree_data$n_visits),
  gye_obs$time
)
# Equal TRUE? Yay!
all.equal(
  rep(tree_data$standID, times = tree_data$n_visits),
  gye_obs$standID
)
# Equal TRUE? Yay!
all.equal(
  rep(tree_data$transectID, times = tree_data$n_visits),
  gye_obs$transectID
)
# Equal TRUE? Yay!
all.equal(
  rep(tree_data$treeID, times = tree_data$n_visits),
  gye_obs$treeID
)
# Equal TRUE? Yay!

# Step 2. Assemble Data for Model ####

## Part A. Pull Together Covariates ####
# Recall (from the paper) that we are interested in the impact of covariates on
# both detection of blister rust and on infection prevalence.

# Observer Experience #
# Indicates if observer had at least one year of experience
obs_data <- gye_obs %>%
  group_by(obsID) %>%
  summarize(first_year = min(year)) %>%
  ungroup()

gye_obs <- left_join(gye_obs, obs_data, by = "obsID") # Join dataframes

gye_obs$experience_ind <- with(gye_obs, ifelse(first_year == year, 0, 1))
experience_ind <- gye_obs$experience_ind # Make a vector of this data

# Hiking Distance #
# Represents a measure in miles of the one-way hiking distance to each stand
stand_hike <- stand_hiking %>%
  select(
    "stand_number" = Stand,
    "total_distance" = OneWayTotalMiles
  )

gye_obs$stand_number <- floor(gye_obs$siteID)
gye_obs <- left_join(gye_obs, stand_hike, by = "stand_number") # Join dataframes
gye_obs <- select(
  gye_obs,
  -stand_number
) # Remove extra column

hike_distance <- gye_obs$total_distance # Make vector of this data

# Stand Elevation #
stand_elev_data <- tran_data %>%
  group_by(standID) %>%
  summarise(elev = mean(elev))
elev_scaled <- as.vector(scale(stand_elev_data$elev))
# Create "elev_scaled"--a measure of the scaled elevation (in meters) of each transect

# Transect Slope #

tran_slope <- tran_data$slope2 # Make vector of transect slope
# This is scaled and centered over a mean of zero
trans_id_new <- as.numeric(as.factor(gye_obs$siteID)) # Make vector of the site IDs

## Part B. Generate Tree Counts ####

n_trees <- nrow(tree_data) # Count of the number of individual trees surveyed, summed across all time periods

standID <- tran_data$standID # Make vector of all stands
n_stands <- max(standID) # How many stands in total are in the data?

n_trans <- nrow(tran_data) # How many transects in total are in the data?

tree_counts <- tree_data %>%
  group_by(time, standID, transectID, siteID) %>%
  summarize(n_trees = length(treeID))
# Create "n_trees"--a measure of the number of trees in each stand/transect by time period

tree_matrix <- model.matrix(~time, tree_counts) # Create matrix of trees through time
head(tree_matrix) # Note: time1 is equivalent to the intercept!

n_xcovs <- ncol(tree_matrix) # How many time periods are in the data?

## Part C. Generate Transect-level Counts ####

trans_id <- as.numeric(as.factor(tree_counts$siteID)) # Make vector of the site IDs

trans_ntrees <- tree_counts %>%
  spread(time, n_trees, fill = 0) %>%
  ungroup() %>%
  select(-standID, -transectID, -siteID) %>%
  as.matrix() # Creates a matrix of the counts of trees per each unique transect by time period
head(trans_ntrees) # time is the columns; siteID is the rows!

trans_year_id <- rep(1:nrow(tree_counts), times = tree_counts$n_trees)
# Create a list of 1 to y, where y is the newest time period (time5, in LMDs case).
# Each entry indicates the time period (from 1 to y) for every tree

n_trans_year <- length(trans_id) # How many transects have been visited since 2004?

## Part D. Generate Input Data from the Observations ####

n_tot <- nrow(gye_obs) # How many total observations are in the data?
n_obs <- max(gye_obs$obsID) # How many observers are in the data?

obsID <- gye_obs$obsID # Make vector of the observers

stand_all <- gye_obs$standID # Make vector of the stand IDs

gyeobs_matrix <- model.matrix(~time, data = gye_obs) # Create a matrix from comb_obs
head(gyeobs_matrix) # Note: time1 is equivalent to the intercept!

all_detections <- gye_obs$inf_ind # Make vector of the blister rust statuses

n_visits <- tree_data$n_visits # Make vector of all the number of visits per tree

tree_id <- rep(c(1:nrow(tree_data)), times = tree_data$n_visits)
# Make a list of tree IDs based upon the number of visits per tree

## Part E. Tree DBH ####

# Create dataframe of DBH data only
tree_dbh_data <- gye_obs %>%
  group_by(time, standID, transectID, siteID, treeID) %>%
  summarize(
    year = first(year),
    dbh1 = first(dbh1),
    dbh2 = first(dbh2),
    meas1 = first(meas_year1),
    meas2 = first(meas_year2)
  )

# Handle NAs in the DBH data
tree_dbh_data$dbh1[which(tree_dbh_data$dbh1 == -999)] <- NA # Assign NAs
tree_dbh_data$dbh1[which(is.na(tree_dbh_data$dbh1) == TRUE)] <- mean(tree_dbh_data$dbh1, na.rm = TRUE)
# If a tree has an earliest DBH value of NA, assign it the mean earliest DBH value instead
tree_dbh_data$dbh2[which(tree_dbh_data$dbh2 == -999)] <- NA # Assign NAs
tree_dbh_data$dbh2[which(is.na(tree_dbh_data$dbh2) == TRUE)] <- mean(tree_dbh_data$dbh2, na.rm = TRUE)
# If a tree has an latest DBH value of NA, assign it the mean latest DBH value instead

# Determine the difference between measurement years
tree_dbh_data <- tree_dbh_data %>%
  mutate(
    dif1 = abs(year - meas1),
    dif2 = abs(year - meas2)
  )
# Create values for the difference (dif1/2) between year surveyed and the earliest
# (meas1) or latest (meas2) year of DBH measurement

# Assign DBH value (either earliest or latestDBH)
tree_dbh_data$dbh3 <- ifelse(tree_dbh_data$dif1 < tree_dbh_data$dif2,
  tree_dbh_data$dbh1, tree_dbh_data$dbh2
)
# Voo Doo R-stats! Explanation below if interested:
# When the difference between the survey year and latest year of DBH measurement is
# larger than the difference between the survey year and earliest year of DBH measurement,
# take the value from dbh1 (earliest DBH); but, when that is smaller, take the
# value from dbh2 (latest DBH)--put the selected value into dbh3.
# The goal is to represent dbh3 as the most likely measurement of DBH for a given
# survey year. Recall that DBH is not measured every survey visit--trees are
# re-measured 10-15 years. This means if the earliest DBH was measured 16 years
# (say, 2004) from a given survey's year (say, 2020), but the latest DBH was measured
# 4 years (say, 2016) from a given survey's year, we will assign "dbh3" as equal
# to the DBH measurement from 2016, because that is more recent! Repeat for every
# year a tree was surveyed.

# Review the distributions of the data
summary(tree_dbh_data$dbh3)

# Scale and center the log DBH
tree_dbh <- as.vector(scale(log(tree_dbh_data$dbh3)))
# Create "tree_dbh"--a measure of the scaled DBH (in cm)

# Calculate the log mean DBH
mean_dbh <- as.vector(tapply(tree_dbh, tree_dbh_data$time, mean))

# Step 3. Run Model ####
## Part A. Define Data ####

model_data <- list(
  "n_trees" = n_trees,
  "n_xcovs" = n_xcovs,
  "xmat" = tree_matrix,
  "n_stands" = n_stands,
  "stand_elev" = elev_scaled,
  "n_trans" = n_trans,
  "stand_id" = standID,
  "trans_slope" = tran_slope,
  "trans_id" = trans_id,
  "trans_year_id" = trans_year_id,
  "n_trans_year" = n_trans_year,
  "n_tot" = n_tot,
  "n_obs" = n_obs,
  "obs_id" = obsID,
  "stand_id2" = stand_all,
  "vmat" = gyeobs_matrix,
  "n_visits" = n_visits,
  "dets" = all_detections,
  "mean_dbh" = mean_dbh,
  "tree_dbh" = tree_dbh,
  "tree_id" = tree_id,
  "trans_id2" = trans_id_new,
  "exp_ind" = experience_ind,
  "hike_dist" = hike_distance
)

# NOTE TO USER: You can save this data object! This is helpful for future use.
# Also, edit the initials/date to your own
saveRDS(model_data, "data/rds/GYE_WBP_ModelData_LMDDec2024_late.rds")

# JUMP-TO-POINT ####
# If you have completed all above steps and want to run the model with your data,
# You can load a previous save of the model data object now:
# Edit the initials/date to your own.
model_data <- readRDS("data/rds/GYE_WBP_ModelData_LMDDec2024_late.rds")

## Part B. Run Model ####

# NOTE TO USER: THIS TAKES A LONG TIME! AND A LOT OF COMPUTER MEMORY! #
# BE PREPARED TO LEAVE YOUR COMPUTER ALONE FOR A COUPLE HOURS. #
# Also: run the "saveRDS" line as well to automatically save the fitted model as
# soon as its done running! Edit initials/date to your own.

model_fit <- sampling(final_model, model_data, pars = c(
  "logit_psi", "log_psi",
  "log1m_psi"
), include = FALSE)
saveRDS(model_fit, "output/saved_fits/Model_Fit_LMDDec2024_late.rds")

# Step 4. Examine Fitted Model ####

# Note to user: if you have already fit the model and want to load a previous save
# of it, do so below:
# Edit initials/date to your own.
model_fit <- readRDS("output/saved_fits/Model_Fit_LMDDec2024_late.rds")

## Part A. Check Convergence through Rhat ####
# Assessing convergence is similar to checking "fit". Our model runs 4 chains for
# 2000 iterations each--a total of 8000 iterations. The goal is for these chains
# to converge & mix, meaning that over their 2000 iterations, each chain arrives
# at or follows very similar values. There are two ways to accomplish this
# below--both necessary to pass.

# In these model summaries, review the column "Rhat". This refers to the Gelman-
# Rubin (1992) statistic R-hat. Per Shanahan et al. (2021), we assess convergence
# as Rhat being below 1.05 for this model.

# First group of parameters:
summary(model_fit,
  pars = c(
    "theta0", "theta1", "theta_year", "alpha_dbh",
    "sigma_stand", "sigma_trans", "sigma_trans_year"
  ),
  use_cache = FALSE
)$summary

# Second group of parameters:
summary(model_fit,
  pars = c(
    "eta0", "eta1", "eta_year", "beta_exp", "beta_hike",
    "beta_dbh", "sigma_obs"
  ),
  use_cache = FALSE
)$summary

## Part B. Check Convergence through Traceplots ####
# Traceplots are the most visual way to check convergence
# Review the below plots--you don't want to see lines shooting off by themselves!

# First group of parameters:
traceplot(model_fit,
  pars = c(
    "theta0", "theta1", "theta_year", "alpha_dbh",
    "sigma_stand", "sigma_trans", "sigma_trans_year"
  ),
  ncol = 1
)

# Second group of parameters:
traceplot(model_fit, pars = c(
  "eta0", "eta1", "eta_year", "beta_exp", "beta_hike",
  "beta_dbh", "sigma_obs"
), ncol = 1)

## Part C. Review Prevalence Estimates ####
# This uses the mean--but Wilson Wright preferred/reported on the median--more on that later.

### NEEDS ATTENTION ####
# To accommodate post-2023 data, you will need to add "prev_t[i]" where
# [i] is equal to the new time period/s.
# IF YOU GET AN ERROR: "Unexpected string constant" it is because you need
# a comma between "prev_t5" and "prev_t[i]"

prev_estimates <- summary(model_fit, pars = c(
  "prev_t1", "prev_t2", "prev_t3",
  "prev_t4", "prev_t5"
), use_cache = FALSE)$summary %>%
  as.data.frame() %>%
  select(p_hat = 1, p_low = 4, p_upp = 8)
prev_estimates$time <- 1:model_data$n_xcovs
# These are our results (first glance)!
prev_estimates

# Step 5. Review Model Results of Infection Prevalence ####

## Part A. Recalculate the Median DBH by Time Period ####

# Calculate log mean and standard deviation of DBH
log_dbh_mean <- mean(log(tree_dbh_data$dbh3))
log_dbh_sd <- sd(log(tree_dbh_data$dbh3))

# Calculate the scaled log median tree DBH per time period
dbh_median <- log(as.vector(tapply(tree_dbh_data$dbh3, tree_dbh_data$time, median))) - log_dbh_mean / log_dbh_sd

## Part B. Extract Prevalence (psi) Parameters from the Fitted Model ####

post_theta0 <- rstan::extract(model_fit, pars = c("theta0"))$theta0
post_alpha1 <- rstan::extract(model_fit, pars = c("alpha_dbh"))$alpha_dbh
post_theta_year <- rstan::extract(model_fit, pars = c("theta_year"))$theta_year

## Part C. Calculate Overall Prevalence Estimates ####

### NEEDS ATTENTION ####
# To accommodate post-2023 data, you will need to add a new line of "prev_t[x]"
# where [x] is equal to the new time period/s. You will also need to update
# the following:
# "dbh_median[x], where [x] is equal to the new time period/s
# "post_theta_year[i, x-1], where [x-1] is equal to the new time period minus 1 (so, if time period 6, you would use 5!)


# Calculate prevalence by time period #
prev_t1 <- prev_t2 <- prev_t3 <- prev_t4 <- prev_t5 <- rep(NA, 4000)
for (i in 1:4000) {
  prev_t1[i] <- plogis(post_theta0[i] + post_alpha1[i] * dbh_median[1])
  prev_t2[i] <- plogis(post_theta0[i] + post_alpha1[i] * dbh_median[2] + post_theta_year[i, 1])
  prev_t3[i] <- plogis(post_theta0[i] + post_alpha1[i] * dbh_median[3] + post_theta_year[i, 2])
  prev_t4[i] <- plogis(post_theta0[i] + post_alpha1[i] * dbh_median[4] + post_theta_year[i, 3])
  prev_t5[i] <- plogis(post_theta0[i] + post_alpha1[i] * dbh_median[5] + post_theta_year[i, 4])
}

# Calculate Point-Estimates #

### NEEDS ATTENTION ####
# To accommodate post-2023 data, you will need to add "prev_t[i]" where
# [i] is equal to the new time period/s.
# IF YOU GET AN ERROR: "Unexpected string constant" it is because you need
# a comma between "(prev_t5)" and "(prev_t[i])"

# For the below code-chunk, copy+paste for each group: "mean(prev_t[i])",
# "quantile(prev_t[i], probs = 0.025)", "quantile(prev_t[i], probs = 0.975)",
# "quantile(prev_t[i], probs = 0.25)", "quantile(prev_t[i], probs = 0.75)",
# "sd(prev_t[i])/sqrt(length(prev_t[i]))", and "sd(prev_t[i])"

overall_prevalence <- data.frame(
  p_hat = c(
    mean(prev_t1), mean(prev_t2), mean(prev_t3), mean(prev_t4),
    mean(prev_t5)
  ),
  p_low = c(
    quantile(prev_t1, probs = 0.025), quantile(prev_t2, probs = 0.025),
    quantile(prev_t3, probs = 0.025), quantile(prev_t4, probs = 0.025),
    quantile(prev_t5, probs = 0.025)
  ),
  p_upp = c(
    quantile(prev_t1, probs = 0.975), quantile(prev_t2, probs = 0.975),
    quantile(prev_t3, probs = 0.975), quantile(prev_t4, probs = 0.975),
    quantile(prev_t5, probs = 0.975)
  ),
  p50_low = c(
    quantile(prev_t1, probs = 0.25), quantile(prev_t2, probs = 0.25),
    quantile(prev_t3, probs = 0.25), quantile(prev_t4, probs = 0.25),
    quantile(prev_t5, probs = 0.25)
  ),
  p50_upp = c(
    quantile(prev_t1, probs = 0.75), quantile(prev_t2, probs = 0.75),
    quantile(prev_t3, probs = 0.75), quantile(prev_t4, probs = 0.75),
    quantile(prev_t5, probs = 0.75)
  ),
  p_se = c(
    sd(prev_t1)/sqrt(length(prev_t1)), sd(prev_t2)/sqrt(length(prev_t2)),
    sd(prev_t3)/sqrt(length(prev_t3)), sd(prev_t4)/sqrt(length(prev_t4)),
    sd(prev_t5)/sqrt(length(prev_t5))
  ),
  p_sd = c(
    sd(prev_t1), sd(prev_t2), sd(prev_t3), sd(prev_t4), sd(prev_t5)
  ),
  time = 1:model_data$n_xcovs
)

overall_prevalence
# These estimates are calculated using the median DBH for each corresponding
# time period and are preferred over estimates calculated with the mean DBH.
# This is due to the spread of the data (sometimes, median is more representative
# of data than a mean) and is what Wilson reported in Shanahan et al. 2021.

# Step 6. Estimate Detection Probability ####

## Part A. Extract Detection (p) Parameters from the Fitted Model ####

post_eta0 <- rstan::extract(model_fit, pars = c("eta0"))$eta0
post_eta_year <- rstan::extract(model_fit, pars = c("eta_year"))$eta_year

## Part B. Calculate Detection Probability through Time ####

# Calculate probabilities

### NEEDS ATTENTION ####
# To accommodate post-2023 data, you will need to add "det_t[i] <-" (where
# [i] is equal to the new time period/s) between "det_t5 <-" and "rep(NA, 4000)"

# IF YOU GET AN ERROR: "Unexpected string constant" it is because you need
# a "<-" between everything
det_t1 <- det_t2 <- det_t3 <- det_t4 <- det_t5 <- rep(NA, 4000)

# For the below code-chunk, copy+paste a new line as the following:
# "det_t[x] <- plogis(post_eta0[i] + post_eta_year[i, x-1])" WHERE: x is equal
# to the new time period/s AND x-1 is equal to the new time period minus one!

for (i in 1:4000) {
  det_t1[i] <- plogis(post_eta0[i])
  det_t2[i] <- plogis(post_eta0[i] + post_eta_year[i, 1])
  det_t3[i] <- plogis(post_eta0[i] + post_eta_year[i, 2])
  det_t4[i] <- plogis(post_eta0[i] + post_eta_year[i, 3])
  det_t5[i] <- plogis(post_eta0[i] + post_eta_year[i, 4])
}
# Calculate point-estimates #

### NEEDS ATTENTION ####
# To accomodate post-2023 data you will need to:
# For the below code-chunk, copy+paste for each group: "mean(det_t[i])",
# "quantile(det_t[i], probs = 0.025)", "quantile(det_t[i], probs = 0.975)"
# IF YOU GET AN ERROR: "Unexpected string constant" it is because you need
# a comma between "(det_t5)" and "(det_t[i])"

detection_by_time <- data.frame(
  p_hat = c(
    mean(det_t1), mean(det_t2), mean(det_t3), mean(det_t4),
    mean(det_t5)
  ),
  p_low = c(
    quantile(det_t1, probs = 0.025), quantile(det_t2, probs = 0.025),
    quantile(det_t3, probs = 0.025), quantile(det_t4, probs = 0.025),
    quantile(det_t5, probs = 0.025)
  ),
  p_upp = c(
    quantile(det_t1, probs = 0.975), quantile(det_t2, probs = 0.975),
    quantile(det_t3, probs = 0.975), quantile(det_t4, probs = 0.975),
    quantile(det_t5, probs = 0.975)
  ),
  time = 1:n_xcovs
)
detection_by_time
# The first entry reports the overall detection probability and is not
# respective of observer priority. For example, the higher the priority,
# the higher likelihood the observer detects blister rust. For this reason,
# the overall estimate is biased lower (as it accommodates data regardless
# of priority). The second entry shows a big boost (adjustment) to the
# overall result and is not so informative. Anecdotally reference the later
# probabilities of detection to understand where detection probability sits
# relatively (and regarding higher priority observers)... Dec2024 is at ~83-85%

## Part C. Calculate Overall Detection Probability ####

# Calculate probabilities
det_overall <- rep(NA, 4000)

for (i in 1:4000) {
  det_overall[i] <- plogis(post_eta0[i])
}
# Calculate point-estimates #
overall_detection <- data.frame(
  p_hat = c(mean(det_overall)),
  p_low = c(quantile(det_overall, probs = 0.025)),
  p_upp = c(quantile(det_overall, probs = 0.975)),
  label = "Overall Detection Probability"
)
overall_detection

# Step 7. Export Results ####

# Don't forget to change initials/date to your own!

# Save Model Data, if not already:
saveRDS(model_data, "data/rds/GYE_WBP_ModelData_LMDDec2024_late.rds")

# Save Fitted Model, if not already:
saveRDS(model_fit, "output/saved_fits/Model_Fit_LMDDec2024_late.rds")

# Save Tree DBH data
write.csv(tree_dbh_data, "data/processed_data/TreeDBHdata_LMDDec2024_late.csv",
  row.names = FALSE
)

# Save Prevalence Estimates:
write.csv(overall_prevalence, "output/results/BlisterRust_PrevalenceEstimates_LMDDec2024.csv",
  row.names = FALSE
)

# Save Detection Estimates:
write.csv(overall_detection, "output/results/BlisterRust_DetectionEstimateOverall_LMDDec2024.csv",
  row.names = FALSE
)
write.csv(detection_by_time, "output/results/BlisterRust_DetectionEstimates_byTime_LMDDec2024.csv",
  row.names = FALSE
)

# END ----