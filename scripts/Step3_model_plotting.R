### Creating Plots/Figures from Model ###
## Lydia M. Druin ##
# GRYN WBP Bayesian Blister Rust #

# Stored as "/model_plotting.R" #
# Updated 03/21/25 #

### Purpose ###
# This code is for recreating figures as seen in Shanahan et al. 2021,
# the Bayesian hierarchical occupancy model of blister rust.
# Scripts "/scripts/data_processing.R" and "/scripts/bayesian_blisterrust.R"
# must be run first before being able to plot.

### Note to User ###
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
library(santoku)

options(scipen = 999) # trick to force R to stay in numeric (not scientific) notation
options(dplyr.width = Inf, dplyr.print_min = 100) # used to set dimensions for viewing objects

options(mc.cores = parallel::detectCores()) # configure rstan to use available cores on your machine
rstan_options(auto_write = TRUE)

# Step 1. Load in Data ####

## Part A. Load in Source File ####

# Note to user: be sure to adapt filenames to reflect how you named them!

# Data developed with data_processing.R #
gye_obs <- readRDS("data/rds/GYE_WBP_Processed_AllObservers_BlisterRustData_LMDDec2024_late.rds")
tran_data <- readRDS("data/rds/GYE_WBP_Processed_TransectData_LMDDec2024.rds")

# Objects developed with bayesian_blisterrust.R #
model_fit <- readRDS("output/saved_fits/Model_Fit_LMDDec2024_late.rds")
model_data <- readRDS("data/rds/GYE_WBP_ModelData_LMDDec2024_late.rds")

# Step 2. Replicate Figures from Shanahan et al. 2021 ####

# Note: Figure 1 is a map of stand locations by number of trees per stand and
# will not be replicated.

## Part A. Figure 2: Tree DBHs ####

# Load in tree DBH data from bayesian_blisterrust.R
tree_dbh_data <- read.csv("data/processed_data/TreeDBHdata_LMDDec2024_late.csv",
  header = TRUE
)

# Create bins for tree DBH
tree_dbh_data$dbh_cat <- as.factor(with(tree_dbh_data, ifelse(dbh3 < 0, "no data",
  ifelse(dbh3 <= 2.5, ">0 to 2.5",
    ifelse(dbh3 <= 10, ">2.5 to 10",
      ifelse(dbh3 <= 30, ">10 to 30", ">30")
    )
  )
)))

# Format data
tree_dbh_data$time <- factor(tree_dbh_data$time, levels = c(1:length(unique(tree_dbh_data$time))))
tree_dbh_data$dbh_cat <- factor(tree_dbh_data$dbh_cat, levels = c(">0 to 2.5", ">2.5 to 10", ">10 to 30", ">30"))

with(tree_dbh_data, table(time, dbh_cat))

# Make plot
figure2 <- ggplot(tree_dbh_data, aes(x = dbh_cat, group = time, fill = time)) +
  geom_bar(aes(fill = time), position = position_dodge(width = 0.75), width = 0.75) +
  theme_grey(base_size = 7) +
  theme(
    axis.title = element_text(size = 14), axis.text = element_text(size = 12),
    legend.position = "bottom", legend.title = element_blank(), legend.text = element_text(size = 12)
  ) +
  scale_fill_brewer(palette = "Set3") +
  xlab("DBH (cm) Size Class") +
  ylab("Number of Whitebark Pine Trees") +
  labs(title = "Number of Trees by Size Class during Time-Steps") +
  guides(fill = guide_legend(title = "Time-Step"))

figure2 # call plot; take a look!

# Save figure
ggsave("output/figures/BlisterRust_TreeDBHcategories_Figure2_LMDDec2024.png",
  plot = last_plot(), width = 8, height = 6.5, units = c("in")
)

## Part B. Figure 3: Blister rust estimates comparison  ####

# Save the results from the ratio estimator, located in
# "GRYN Whitebark Ratio Estimators/output/results/" to the file:
# "GRYN Whitebark Bayesian Blister Rust/output/results/" to complete this figure.

# Load in Ratio Estimator results and add 50% Confidence Intervals
ratio_est <- read.csv("output/results/RatioEst_BlisterRust_LMDDec2024.csv",
  header = TRUE
)
ratio_est$p50_low <- with(ratio_est, p_hat - (0.674 * p_se))
ratio_est$p50_upp <- with(ratio_est, p_hat + (0.674 * p_se))
ratio_est$p_sd = with(ratio_est, sqrt(p_var))

# Load in Bayesian estimates and add 50% Credible Intervals
bayes_est <- read.csv("output/results/BlisterRust_PrevalenceEstimates_LMDDec2024.csv",
  header = TRUE
)

# Combine the two methods
compare_est <- bind_rows(
  select(ratio_est, time = Time, p_hat, p_low, p50_low, p_upp, p50_upp, p_sd),
  select(bayes_est, time, p_hat, p_low, p50_low, p_upp, p50_upp, p_sd)
)
compare_est$method <- rep(c("Design-based Ratio Estimator", "Bayesian Hierarchical Occupancy Model"),
  each = length(unique(bayes_est$time))
)
compare_est$method = factor(compare_est$method, levels = c("Design-based Ratio Estimator", "Bayesian Hierarchical Occupancy Model"))

# Note on time: if you are including additional time periods beyond Time 5, you
# will need to manually add that in to the above line of code!
# In quotation marks, add "T[i] (years-years)" where i = time period #.
# OR--DONT RUN THE BELOW LINE OF CODE. The plot will still work! It just won't
# be formatted using the same x-axis labels as previous pubs/reports.
compare_est$time = rep(c("T1 (2004-2007)", "T2 (2008-2011)", "T3 (2012-2015)", "T4 (2016-2019)", "T5 (2020-2023)"), 2)

# Make figure
figure3 <- ggplot(
  data = compare_est,
  aes(x = factor(time), y = p_hat, group = method)) +
  geom_errorbar(aes(ymin = p_low, ymax = p_upp),
    width = 0, size = 0.71,
    position = position_dodge(width = 0.5)
  ) +
  geom_errorbar(aes(ymin = p50_low, ymax = p50_upp, color = method),
    width = 0, size = 2.25, position = position_dodge(width = 0.5)
  ) +
  geom_point(position = position_dodge(width = 0.5)) +
  theme_grey() +
  theme(
    axis.title = element_text(size = 14), axis.text = element_text(size = 12),
    legend.position = "bottom", legend.title = element_blank(), legend.text = element_text(size = 12)
  ) +
  scale_color_brewer(palette = "Set3") +
  xlab("Time-Step") +
  ylab("Prevalence point estimates") +
  labs(title = "Comparison of Blister Rust Infection Estimation Methods")
figure3 # call plot; take a look!

# Black bands indicate 95% intervals; colored bands indicate 50% intervals

# Save figure
ggsave("output/figures/BlisterRust_MethodsComparison_Figure3_LMDDec2024.png",
  plot = last_plot(), width = 8, height = 6.5, units = c("in")
)

## Part C. Figure 4: Posterior prevalence parameters ####

# Call posterior intercept parameters from blister rust prevalence ###

# NEEDS ATTENTION ####
# In the line "post_intercepts$param" you will make edits to incorporate
# new time periods! Specifically, after "Time 5 Adj.", add a comma! Then, add
# "Time [i] Adj.", enclosed by quotation marks, where [i] is the new time period/s

post_intercepts <- summary(model_fit,
  pars = c("theta0", "theta_year"),
  use_cache = FALSE
)$summary %>%
  as.data.frame() %>%
  select(1, 4, 8, 5, 7)
names(post_intercepts) <- c("p_mean", "p_low", "p_upp", "p50_low", "p50_upp")
post_intercepts$param <- c("Time 1", "Time 2 Adj.", "Time 3 Adj.", "Time 4 Adj.", "Time 5 Adj.")
post_intercepts$group <- c("Overall", rep("Adjustments", length(post_intercepts$param) - 1))
post_intercepts$group <- factor(post_intercepts$group, levels = c("Overall", "Adjustments"))

# Make figure
figure4 <- ggplot(post_intercepts, aes(x = factor(param), y = p_mean)) +
  geom_errorbar(aes(ymin = p_low, ymax = p_upp),
    width = 0, show.legend = FALSE,
    size = 0.71
  ) +
  geom_errorbar(aes(ymin = p50_low, ymax = p50_upp),
    width = 0, color = "#7570b3",
    size = 2.25
  ) +
  geom_point(show.legend = FALSE) +
  theme_grey() +
  theme(axis.title = element_text(size = 14), axis.text = element_text(size = 12)) +
  facet_grid(~group, scales = "free_x", space = "free_x") +
  geom_hline(yintercept = 0, lty = 2) +
  xlab("Intercept parameter") +
  ylab("Posterior distribution (logit scale)")
figure4 # call plot; take a look!

# Save figure
ggsave("output/figures/BlisterRust_PrevalenceParameters_Figure4_LMDDec2024.png",
  plot = last_plot(), width = 8, height = 6.5, units = c("in")
)

## Part D. Figure 5: Probability of infection by DBH ####

# Extract parameters from model
post_theta0 <- rstan::extract(model_fit, pars = c("theta0"))$theta0
post_alpha1 <- rstan::extract(model_fit, pars = c("alpha_dbh"))$alpha_dbh
post_theta_year <- rstan::extract(model_fit, pars = c("theta_year"))$theta_year

# Calculate log mean and standard deviation of DBH
log_dbh_mean <- mean(log(tree_dbh_data$dbh3))
log_dbh_sd <- sd(log(tree_dbh_data$dbh3))

dbh_possible <- c(1:length(unique(tree_dbh_data$dbh3))) # Number of possible DBH values
log_dbh_poss <- (log(dbh_possible) - log_dbh_mean) / log_dbh_sd # Create logged values of possible DBHs

# Calculate prevalence estimates by log DBH

### NEEDS ATTENTION ####
# To accommodate post-2023 data, you will need to add a new line of "dbh_prevX[,i]"
# where X in "dbh_prevX" is equal to the new time period/s (e.g., 6)
# ALSO RENAME, at the end: "post_theta_year[i, X-1]" where X-1 is equal to the
# new time period/s minus one (e.g., 6 minus 1 is 5)

dbh_prev1 <- dbh_prev2 <- dbh_prev3 <- dbh_prev4 <- dbh_prev5 <- matrix(NA, ncol = 4000, nrow = length(dbh_possible))
for (i in 1:4000) {
  dbh_prev1[, i] <- plogis(post_theta0[i] + post_alpha1[i] * log_dbh_poss)
  dbh_prev2[, i] <- plogis(post_theta0[i] + post_alpha1[i] * log_dbh_poss + post_theta_year[i, 1])
  dbh_prev3[, i] <- plogis(post_theta0[i] + post_alpha1[i] * log_dbh_poss + post_theta_year[i, 2])
  dbh_prev4[, i] <- plogis(post_theta0[i] + post_alpha1[i] * log_dbh_poss + post_theta_year[i, 3])
  dbh_prev5[, i] <- plogis(post_theta0[i] + post_alpha1[i] * log_dbh_poss + post_theta_year[i, 4])
}

### NEEDS ATTENTION ####
# To accommodate post-2023 data, you will need to:
# For the below code-chunk, copy+paste for each group: 
# "apply(dbh_prevX, 1, mean)", "apply(dbh_prevX, 1, quantile, probs = 0.025)",
# "apply(dbh_prevX, 1, quantile, probs = 0.975)"

# Format prevalence estimates into dataframe
dbh_summary <- data.frame(
  dbh = rep(dbh_possible, length(unique(tree_dbh_data$time))),
  p_mean = c(
    apply(dbh_prev1, 1, mean),
    apply(dbh_prev2, 1, mean),
    apply(dbh_prev3, 1, mean),
    apply(dbh_prev4, 1, mean),
    apply(dbh_prev5, 1, mean)
  ),
  p_low = c(
    apply(dbh_prev1, 1, quantile, probs = 0.025),
    apply(dbh_prev2, 1, quantile, probs = 0.025),
    apply(dbh_prev3, 1, quantile, probs = 0.025),
    apply(dbh_prev4, 1, quantile, probs = 0.025),
    apply(dbh_prev5, 1, quantile, probs = 0.025)
  ),
  p_upp = c(
    apply(dbh_prev1, 1, quantile, probs = 0.975),
    apply(dbh_prev2, 1, quantile, probs = 0.975),
    apply(dbh_prev3, 1, quantile, probs = 0.975),
    apply(dbh_prev4, 1, quantile, probs = 0.975),
    apply(dbh_prev5, 1, quantile, probs = 0.975)
  ),
  Time = rep(c(unique(tree_dbh_data$time)), each = length(dbh_possible))
)
dbh_summary$Time <- factor(dbh_summary$Time)

# Make plot
figure5 <- ggplot(dbh_summary, aes(x = dbh, y = p_mean, group = Time, color = Time)) +
  geom_ribbon(aes(ymin = p_low, ymax = p_upp, fill = Time, color = Time), alpha = 0.1, color = NA) +
  geom_line() +
  theme_grey() +
  theme(
    axis.title = element_text(size = 14), axis.text = element_text(size = 12),
    legend.position = "bottom", legend.title = element_blank(), legend.text = element_text(size = 12)
  ) +
  scale_color_brewer(palette = "Set3") +
  scale_fill_brewer(palette = "Set3") +
  xlab("Tree DBH (cm)") +
  ylab("Probability of Infection") +
  labs(title = "Likelihood of Blister Rust Infection by DBH during Time-Steps")
figure5 # call plot; take a look!

# Save figure
ggsave("output/figures/BlisterRust_Prevalence_by_DBH_Figure5_LMDDec2024.png",
  plot = last_plot(), width = 8, height = 6.5, units = c("in")
)

# Make plot comparing Time 1 with the newest Time Period
newest_time = length(unique(dbh_summary$Time))
subset_dbh = subset(dbh_summary, subset = Time %in% c(1, newest_time))
figure5_alt <- ggplot(subset_dbh, aes(x = dbh, y = p_mean, group = Time, color = Time)) +
  geom_ribbon(aes(ymin = p_low, ymax = p_upp, fill = Time, color = Time), alpha = 0.1, color = NA) +
  geom_line() +
  theme_grey() +
  theme(
    axis.title = element_text(size = 14), axis.text = element_text(size = 12),
    legend.position = "bottom", legend.title = element_blank(), legend.text = element_text(size = 12)
  ) +
  scale_color_brewer(palette = "Set3") +
  scale_fill_brewer(palette = "Set3") +
  xlab("Tree DBH (cm)") +
  ylab("Probability of Infection") +
  labs(title = "Likelihood of Blister Rust Infection by DBH during Time-Steps")
figure5_alt # call plot; take a look!

# Save figure
ggsave("output/figures/BlisterRust_Prevalence_by_DBH_Figure5_OldestNewestTimePeriods_LMDDec2024.png",
       plot = last_plot(), width = 8, height = 6.5, units = c("in")
)

# Step 3. Replicate Figures from Appendices ####

## Part A. Appendix 1 ####

# Figures from Appendix 1 involve the model-based ratio estimator and HLogis1.
# Because we do not utilize either of these two approaches beyond the paper
# itself, none of the figures are replicated in this code.

## Part B. Appendix 2 ####

# This appendix addresses model residuals for all the Bayesian approaches in
# the paper. Because we only utilize HOccu2, we will only replicate these
# residual plots.

### Figure S2.8: Binned occupancy residuals ####
# y = binned occupancy residuals, x = tree dbh
# y = binned occupancy residuals, x = log tree dbh
# should be 4000 points (1 for each posterior draw)
# bins have ~125 trees each
# Not able to be replicated at this time; LMD 03/24/25

### Figure S2.9: Stand coefficient residuals ####
# y = stand coefficient residuals, x = standardized stand elevation
# should be 4000 points (1 for each posterior draw)
# Not able to be replicated at this time; LMD 03/24/25

### Figure S2.10: Binned detection residuals ####
# y = binned detection residuals, x = tree dbh
# y = binned occupancy residuals, x = log tree dbh
# should be 4000 points (1 for each posterior draw)
# bins have ~50 observations each
# Not able to be replicated at this time; LMD 03/24/25

## Part C. Appendix 3 ####

### Figure S3: Posterior Detection Parameters ####

# Call posterior detection parameters ###

det_coefficients <- summary(model_fit,
  pars = c(
    "eta0", "eta_year",
    "beta_hike", "beta_exp", "eta1", "beta_dbh", "sigma_obs"
  ),
  use_cache = FALSE
)$summary %>%
  as.data.frame() %>%
  select(1, 4, 8, 5, 7)
names(det_coefficients) <- c("p_mean", "p_low", "p_upp", "p50_low", "p50_upp")
det_coefficients$param <- c(
  post_intercepts$param, "Hike Dist.", "Obs. Experience", "Slope",
  "Tree DBH", "Observer"
)
det_coefficients$param <- factor(det_coefficients$param, levels = c(
  post_intercepts$param, "Hike Dist.", "Obs. Experience", "Slope",
  "Tree DBH", "Observer"
))
det_coefficients$group <- c(rep("Intercept", length(unique(tree_dbh_data$time))), rep("Coefficient", 4), "Std. Dev.")
det_coefficients$group <- factor(det_coefficients$group, levels = c("Intercept", "Coefficient", "Std. Dev."))

# Make figure
figureS3 <- ggplot(det_coefficients, aes(x = param, y = p_mean)) +
  geom_errorbar(aes(ymin = p_low, ymax = p_upp),
    width = 0, show.legend = FALSE,
    size = 0.71
  ) +
  geom_errorbar(aes(ymin = p50_low, ymax = p50_upp),
    width = 0, color = "darkgreen",
    size = 2.25
  ) +
  geom_point(show.legend = FALSE) +
  theme_bw() +
  theme(axis.title = element_text(size = 14), axis.text = element_text(size = 12)) +
  facet_grid(~group, scales = "free_x", space = "free_x") +
  geom_hline(yintercept = 0, lty = 2) +
  xlab("Intercept parameter") +
  ylab("Posterior distribution (logit scale)")
figureS3 # call plot; take a look!

# Save figure
ggsave("output/figures/BlisterRust_DetectionParameters_Figure4_LMDDec2024.png",
  plot = last_plot(), width = 12, height = 6.5, units = c("in")
)

# Step 4. Miscellaneous Plots ####

## Part A. Infection by Time and DBH group ####

# Call infection parameters by DBH size and time

# NEEDS ATTENTION ####
# For incorporating post-2023 data, you will need to:
# After "small_t5" (the final entry under the "small" category), add in 
# "small_t[i]" were [i] is equal to the new time period/s, followed by a comma
# Repeat for "med", "large" and "xlarge" categories--EXCEPT do not add a comma
# after the "xlarge_t[i]" final addition.

size_summary <- summary(model_fit,
  pars = c(
    "small_t1", "small_t2", "small_t3", "small_t4", "small_t5",
    "med_t1", "med_t2", "med_t3", "med_t4", "med_t5",
    "large_t1", "large_t2", "large_t3", "large_t4", "large_t5",
    "xlarge_t1", "xlarge_t2", "xlarge_t3", "xlarge_t4", "xlarge_t5"
  ),
  use_cache = FALSE
) %>%
  as.data.frame() %>%
  select(1, 4, 8, 5, 7)
names(size_summary) <- c("p_mean", "p_low", "p_upp", "p50_low", "p50_upp")
size_summary$time <- rep(c(unique(tree_dbh_data$time)), 4)
size_summary$size <- c(rep("small", length(unique(tree_dbh_data$time))), 
                           rep("med", length(unique(tree_dbh_data$time))), 
                           rep("large", length(unique(tree_dbh_data$time))), 
                           rep("xlarge", length(unique(tree_dbh_data$time))))
size_summary$size <- factor(size_summary$size,
  levels = c("small", "med", "large", "xlarge")
)

# Make figure
figure_misc1_fixed <- ggplot(size_summary, aes(x = time, y = p_mean)) +
  geom_point() +
  geom_errorbar(aes(ymin = p_low, ymax = p_upp),
    width = 0, show.legend = FALSE,
    size = 0.71
  ) +
  geom_errorbar(aes(ymin = p50_low, ymax = p50_upp),
    width = 0, color = "darkgreen",
    size = 2.25
  ) +
  geom_point(show.legend = FALSE) +
  facet_grid(size ~ ., scale = "fixed") +
  theme_grey() +
  theme(axis.title = element_text(size = 14), axis.text = element_text(size = 12),
        legend.position = "bottom", legend.title = element_blank(), legend.text = element_text(size = 12))+
    xlab("Time-Step") +
  ylab("Blister Rust Infection Prevalence") +
  labs(title = "Greater Yellowstone Whitebark Pine Blister Rust Infection Prevalence")
figure_misc1_fixed # call plot; take a look!

# NOTE: This plot has a fixed (same) y-axis for all of the panels

# Save figure
ggsave("output/figures/BlisterRust_Infection_by_DBHTime_Figure_Misc1_SameYaxis_LMDDec2024.png",
  plot = last_plot(), width = 8, height = 6.5, units = c("in")
)

# NOTE: This plot has a y-axis that VARIES BY PANEL! It makes it easier to see
# the smaller values (like with small DBH category) but falsely distorts results
# AKA, the Time 5 result is not at 100% for every DBH category
figure_misc1 <- ggplot(size_summary, aes(x = time, y = p_mean)) +
  geom_point() +
  geom_errorbar(aes(ymin = p_low, ymax = p_upp),
                width = 0, show.legend = FALSE,
                size = 0.71
  ) +
  geom_errorbar(aes(ymin = p50_low, ymax = p50_upp),
                width = 0, color = "darkgreen",
                size = 2.25
  ) +
  geom_point(show.legend = FALSE) +
  facet_grid(size ~ ., scale = "free") +
  theme_grey() +
  theme(axis.title = element_text(size = 14), axis.text = element_text(size = 12),
        legend.position = "bottom", legend.title = element_blank(), legend.text = element_text(size = 12))+
  xlab("Time-Step") +
  ylab("Blister Rust Infection Prevalence") +
  labs(title = "Greater Yellowstone Whitebark Pine Blister Rust Infection Prevalence")
figure_misc1 # call plot; take a look!

# This plot has a custom y-axis based on each panel!

# Save figure
ggsave("output/figures/BlisterRust_Infection_by_DBHTime_Figure_Misc1_VariedYaxis_LMDDec2024.png",
       plot = last_plot(), width = 8, height = 6.5, units = c("in")
)

# END ----