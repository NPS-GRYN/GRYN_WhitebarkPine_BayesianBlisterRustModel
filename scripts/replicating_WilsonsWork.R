### WW Replication ###
## Lydia M. Druin ##
# GRYN WBP Bayesian Blister Rust #

# Stored as "/replicating_WilsonsWork.R" #
# Updated 03/21/25 #

### Purpose ###
# This code is for replicating WW's work from 2020 in the Bayesian hierarchical
# occupancy model described in Shanahan et al. 2021

### Note to User ###
# This code was originally written by Wilson Wright (11/19/19) and adapted by LMD.

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

# Note: the source file/s were processed using "/WilsonsWork-ErinCodefolder/DataProcessing.R"
# BELOW THIS POINT: ALL CODE IS WILSON'S with minor edits by LMD for function!

##Load data
comb_obs <- readRDS('WilsonsWork-ErinCodefolder/comb_obs.rds')
tran_data <- readRDS('WilsonsWork-ErinCodefolder/tran_data.rds')
stand_hiking <- read.csv('WilsonsWork-ErinCodefolder/stand_hiking.csv')

##Restructure data some
#Add time period factor
comb_obs$time <- factor(comb_obs$year)
levels(comb_obs$time) <- c("1", "1", "1", "1",
                           "2", "2", "2", "2",
                           "3", "3", "3", "3",
                           "4", "4", "4", "4")

#Arrange by time period, stand, then tree.
#This makes repeat observations for a tree within a time period appear
#in consecutive rows in the data frame
comb_obs <- arrange(comb_obs, time, stand_id, tran_id, tree_id)

#Create tree-level data frame and include number of visits per time period.
#note that this is really a tree by time period data frame.
#Also add a naive infection indicator based on the first observation,
#this should be modified to match which observers were prioritized on the
#multi-observer transects.
# LMD NOTE; THIS ABOVE NOTE FROM WILSON IS CRUCIAL--whether he intended to or not,
# the observations were NOT sorted by observer priority, and INSTEAD sorted by
# alphabetized last name of observer. THIS HAS BEEN CORRECTED IN 2024! And is
# a PRIMARY SOURCE OF DIFFERENCES IN RESULTS 2020 vs 2024 WORK!
tree_data <- comb_obs %>%
  group_by(time, stand_id, tran_id, site_id, tree_id) %>%
  summarise(n_visits = length(inf_ind),
            naive_inf = first(inf_ind))

#double check that ordering for everything is correct,
#all of these should be true if everything worked properly.
all.equal(rep(tree_data$time, times = tree_data$n_visits),
          comb_obs$time)
all.equal(rep(tree_data$stand_id, times = tree_data$n_visits),
          comb_obs$stand_id)
all.equal(rep(tree_data$tran_id, times = tree_data$n_visits),
          comb_obs$tran_id)
all.equal(rep(tree_data$tree_id, times = tree_data$n_visits),
          comb_obs$tree_id)

#Get some additional covariates for the (eventual) detection level 
#of the model. Many of these were included in the previous analysis.
obs_data <- comb_obs %>%
  group_by(obs_ID) %>%
  summarise(first_year = min(year)) %>%
  ungroup()

comb_obs <- left_join(comb_obs, obs_data, by = 'obs_ID')

comb_obs$exp_ind <- with(comb_obs, ifelse(first_year == year, 0, 1))

stand_hiking2 <- stand_hiking %>%
  select('stand_num' = Stand,
         'total' = OneWayTotalMiles)

comb_obs$stand_num <- floor(comb_obs$site_id)
comb_obs <- left_join(comb_obs, stand_hiking2, by = 'stand_num')
comb_obs <- select(comb_obs, -stand_num)
names(comb_obs)[20] <- 'total_hike'

####Setup data for Stan####
n_trees1 <- nrow(tree_data)

stand_id1 <- tran_data$stand_id
n_stands1 <- max(stand_id1)

n_trans1 <- nrow(tran_data)

stand_data1 <- tran_data %>%
  group_by(stand_id) %>%
  summarise(area_raw = first(area))
stand_data1$log_area_std <- as.vector(scale(log(stand_data1$area_raw)))

tree_counts1 <- tree_data %>%
  group_by(time, stand_id, tran_id, site_id) %>%
  summarise(n_trees = length(tree_id))

xmat1 <- model.matrix(~time, tree_counts1)
n_xcovs1 <- ncol(xmat1)
trans_id1 <- as.numeric(as.factor(tree_counts1$site_id))

trans_ntrees1 <- tree_counts1 %>%
  spread(time, n_trees, fill = 0) %>%
  ungroup() %>%
  select(4:7) %>%
  as.matrix()

trans_year_id <- rep(1:nrow(tree_counts1), times = tree_counts1$n_trees)
n_trans_year <- length(trans_id1)

n_tot <- nrow(comb_obs)
n_obs <- max(comb_obs$obs_ID)
obs_id <- comb_obs$obs_ID
stand_new <- comb_obs$stand_id
vmat <- model.matrix(~time, data = comb_obs)
dets_all <- comb_obs$inf_ind

n_visits <- tree_data$n_visits

#Final data pieces
stand_elev_data <- tran_data %>%
  group_by(stand_id) %>%
  summarise(elev = mean(elev))
stand_elev <- as.vector(scale(stand_elev_data$elev))

tree_dbh_data2 <- comb_obs %>%
  group_by(time, stand_id, tran_id, site_id, tree_id) %>%
  summarise(year = first(year),
            dbh1 = first(dbh1),
            dbh2 = first(dbh2),
            meas1 = first(meas_year1),
            meas2 = first(meas_year2))

tree_dbh_data2$dbh1[which(tree_dbh_data2$dbh1 == -999)] <- NA
tree_dbh_data2$dbh1[which(is.na(tree_dbh_data2$dbh1) == TRUE)] <- 
  mean(tree_dbh_data2$dbh1, na.rm = TRUE)

tree_dbh_data2$dbh2[which(tree_dbh_data2$dbh2 == -999)] <- NA
tree_dbh_data2$dbh2[which(is.na(tree_dbh_data2$dbh2) == TRUE)] <- 
  mean(tree_dbh_data2$dbh2, na.rm = TRUE)

tree_dbh_data2 <- tree_dbh_data2 %>%
  mutate(dif1 = abs(year - meas1),
         dif2 = abs(year - meas2))

tree_dbh_data2$dbh3 <- ifelse(tree_dbh_data2$dif1 < tree_dbh_data2$dif2,
                              tree_dbh_data2$dbh1, tree_dbh_data2$dbh2)

tree_dbh <- as.vector(scale(log(tree_dbh_data2$dbh3)))
mean_dbh <- as.vector(tapply(tree_dbh, tree_dbh_data2$time, mean))

tran_slope <- tran_data$slope2

tree_id <- rep(c(1:nrow(tree_data)), times = tree_data$n_visits)

trans_id_new <- as.numeric(as.factor(comb_obs$site_id))

exp_ind <- comb_obs$exp_ind
hike_dist <- comb_obs$total_hike

####Final model fit####
##In the final model we expand upon the basic model to include more
##of the covariates that are available for detection and prevalence.
##These include stand elevation, log tree dbh (can influence both),
##and transect slope for detection.

##Compile model
final_model1 <- stan_model('WilsonsWork-ErinCodefolder/final_model1.stan')

final_data1 <- list('n_trees' = n_trees1,
                    'n_xcovs' = n_xcovs1,
                    'xmat' = xmat1,
                    'n_stands' = n_stands1,
                    'stand_elev' = stand_elev,
                    'n_trans' = n_trans1,
                    'stand_id' = stand_id1,
                    'trans_slope' = tran_slope,
                    'trans_id' = trans_id1,
                    'trans_year_id' = trans_year_id,
                    'n_trans_year' = n_trans_year,
                    'n_tot' = n_tot,
                    'n_obs' = n_obs,
                    'obs_id' = obs_id,
                    'stand_id2' = stand_new,
                    'vmat' = vmat,
                    'n_visits' = n_visits,
                    'dets' = dets_all,
                    'mean_dbh' = mean_dbh,
                    'tree_dbh' = tree_dbh,
                    'tree_id' = tree_id,
                    'trans_id2' = trans_id_new,
                    'exp_ind' = exp_ind,
                    'hike_dist' = hike_dist)
# SAVE DATA, LMD:
saveRDS(final_data1, "WilsonsWork-ErinCodefolder/final_data1.rds")
saveRDS(tree_dbh_data2, "WilsonsWork-ErinCodefolder/tree_dbh_data.rds")

##Sample parameters
#final_fit1 <- sampling(final_model1, final_data1,
#                       pars = c('logit_psi', 'log_psi', 'log1m_psi'),
#                       include = FALSE)
final_fit1 = readRDS("WilsonsWork-ErinCodefolder/final_fit1.rds")

##checking convergence
summary(final_fit1, pars = c('theta0', 'theta1', 'theta_year',
                             'alpha_dbh',
                             'sigma_stand', 'sigma_trans',
                             'sigma_trans_year'),
        use_cache = FALSE)$summary
summary(final_fit1, pars = c('eta0', 'eta1', 'eta_year',
                             'beta_exp', 'beta_hike',
                             'beta_dbh',
                             'sigma_obs'),
        use_cache = FALSE)$summary

traceplot(final_fit1, pars = c('theta0', 'theta1', 'theta_year',
                               'alpha_dbh',
                               'sigma_stand', 'sigma_trans',
                               'sigma_trans_year'), ncol = 1)
traceplot(final_fit1, pars = c('eta0', 'eta1', 'eta_year',
                               'beta_exp', 'beta_hike',
                               'beta_dbh',
                               'sigma_obs'), ncol = 1)

##Look at the overall prevalence estimates (note that this is using the
#mean from each time period. I think using the median is better, see below).
final_est1 <- summary(final_fit1,
                      pars = c('prev_t1', 'prev_t2', 'prev_t3', 'prev_t4'),
                      use_cache = FALSE)$summary %>%
  as.data.frame() %>%
  select(p_hat = 1, p_low = 4, p_upp = 8)
final_est1$time = c(1, 2, 3, 4)

##Recalculate the "expected" tree size each time period using the median.
##This might make more sense because the median tree is invariant to the
##log transformation we are using for tree dbh.
logdbh_mean <- mean(log(tree_dbh_data2$dbh3))
logdbh_sd <- sd(log(tree_dbh_data2$dbh3))
dbh_median <- (log(as.vector(tapply(tree_dbh_data2$dbh3, tree_dbh_data2$time, median))) - logdbh_mean) / logdbh_sd

post_theta0 <- rstan::extract(final_fit1, pars = c('theta0'))$theta0
post_alpha1 <- rstan::extract(final_fit1, pars = c('alpha_dbh'))$alpha_dbh
post_theta_year <- rstan::extract(final_fit1, pars = c('theta_year'))$theta_year

prev2_t1 <- prev2_t2 <- prev2_t3 <- prev2_t4 <- rep(NA, 4000)
for(i in 1:4000){
  prev2_t1[i] <- plogis(post_theta0[i] + post_alpha1[i] * dbh_median[1])
  prev2_t2[i] <- plogis(post_theta0[i] + post_alpha1[i] * dbh_median[2] + post_theta_year[i, 1])
  prev2_t3[i] <- plogis(post_theta0[i] + post_alpha1[i] * dbh_median[3] + post_theta_year[i, 2])
  prev2_t4[i] <- plogis(post_theta0[i] + post_alpha1[i] * dbh_median[4] + post_theta_year[i, 3])
}

#Reexamine the overall prevalence estimates and make a new plot
##Look at the overall prevalence estimates
final_est2 <- data.frame(p_hat = c(mean(prev2_t1), mean(prev2_t2), mean(prev2_t3), mean(prev2_t4)),
                         p_low = c(quantile(prev2_t1, probs = 0.025), quantile(prev2_t2, probs = 0.025),
                                   quantile(prev2_t3, probs = 0.025), quantile(prev2_t4, probs = 0.025)),
                         p_upp = c(quantile(prev2_t1, probs = 0.975), quantile(prev2_t2, probs = 0.975),
                                   quantile(prev2_t3, probs = 0.975), quantile(prev2_t4, probs = 0.975)),
                         time = c(1, 2, 3, 4))
final_est2
write.csv(final_est2, "WilsonsWork-ErinCodefolder/overall_prevalence.csv")

# LMD: Detection
post_eta0 <- rstan::extract(model_fit, pars = c("eta0"))$eta0
post_eta_year <- rstan::extract(model_fit, pars = c("eta_year"))$eta_year
det_overall <- rep(NA, 4000)

for (i in 1:4000) {
  det_overall[i] <- plogis(post_eta0[i])
}
# Calculate point-estimates #
final_det <- data.frame(
  p_hat = c(mean(det_overall)),
  p_low = c(quantile(det_overall, probs = 0.025)),
  p_upp = c(quantile(det_overall, probs = 0.975)),
  label = "Overall Detection Probability"
)
final_det
write.csv(final_det, "WilsonsWork-ErinCodefolder/overall_detection.csv")

####Some plots####
size_summ <- summary(final_fit1,
                     pars = c('small_t1', 'small_t2', 'small_t3', 'small_t4',
                              'med_t1', 'med_t2', 'med_t3', 'med_t4',
                              'large_t1', 'large_t2', 'large_t3', 'large_t4',
                              'xlarge_t1', 'xlarge_t2', 'xlarge_t3', 'xlarge_t4'),
                     use_cache = FALSE) %>%
  as.data.frame() %>%
  select(1, 4, 8)
names(size_summ) <- c('p_mean', 'p_low', 'p_upp')
size_summ$time <- rep(c(1, 2, 3, 4), 4)
size_summ$size <- rep(c('small', 'med', 'large', 'xlarge'), each = 4)
size_summ$size <- factor(size_summ$size,
                         levels = c('small', 'med', 'large', 'xlarge'))

ggplot(size_summ, aes(x = time, y = p_mean)) + 
  geom_point() +
  geom_errorbar(aes(ymin = p_low, ymax = p_upp), width = 0) +
  theme_bw() +
  facet_grid(size ~ ., scale = "free") +
  xlab("Time Period") +
  ylab("Prevalence point estimates and 95% intervals") +
  ggtitle('Prevalence by time period and DBH size class')
ggsave("WilsonsWork-ErinCodefolder/Prevalence_byTime_andDBH.png",
       plot = last_plot(), width = 8, height = 6.5, units = c("in")
)

##Make a plot of estimated prevalence versus dbh for each time period
dbh_poss <- c(1:127)
ldbh_poss <- (log(dbh_poss) - logdbh_mean) / logdbh_sd

dbh_prev1 <- dbh_prev2 <- dbh_prev3 <- dbh_prev4 <- matrix(NA, ncol = 4000, nrow = 127)
for(i in 1:4000){
  dbh_prev1[, i] <- plogis(post_theta0[i] + post_alpha1[i] * ldbh_poss)
  dbh_prev2[, i] <- plogis(post_theta0[i] + post_alpha1[i] * ldbh_poss + post_theta_year[i, 1])
  dbh_prev3[, i] <- plogis(post_theta0[i] + post_alpha1[i] * ldbh_poss + post_theta_year[i, 2])
  dbh_prev4[, i] <- plogis(post_theta0[i] + post_alpha1[i] * ldbh_poss + post_theta_year[i, 3])
}

dbh_summ <- data.frame(dbh = rep(dbh_poss, 4),
                       p_mean = c(apply(dbh_prev1, 1, mean),
                                  apply(dbh_prev2, 1, mean),
                                  apply(dbh_prev3, 1, mean),
                                  apply(dbh_prev4, 1, mean)),
                       p_low  = c(apply(dbh_prev1, 1, quantile, probs = 0.025),
                                  apply(dbh_prev2, 1, quantile, probs = 0.025),
                                  apply(dbh_prev3, 1, quantile, probs = 0.025),
                                  apply(dbh_prev4, 1, quantile, probs = 0.025)),
                       p_upp  = c(apply(dbh_prev1, 1, quantile, probs = 0.975),
                                  apply(dbh_prev2, 1, quantile, probs = 0.975),
                                  apply(dbh_prev3, 1, quantile, probs = 0.975),
                                  apply(dbh_prev4, 1, quantile, probs = 0.975)),
                       Time   = rep(c(1, 2, 3, 4), each = 127))
dbh_summ$Time <- factor(dbh_summ$Time)

ggplot(dbh_summ, aes(x = dbh, y = p_mean, group = Time, color = Time)) +
  geom_ribbon(aes(ymin = p_low, ymax = p_upp, fill = Time), alpha = 0.1, color = NA) +
  geom_line() +
  theme_bw() +
  xlab('Tree DBH (cm)') +
  ylab('Probability of Infection') +
  ggtitle('Posterior mean and 95% interval for infection probability \nversus tree DBH in each time period')
ggsave("WilsonsWork-ErinCodefolder/DBH_PosteriorMeans.png",
       plot = last_plot(), width = 8, height = 6.5, units = c("in")
)

#Summary plot for the intercepts
plot(final_fit1, pars = c('theta0', 'theta_year')) +
  ggtitle('Posterior distribution of intercept parameters')
ggsave("WilsonsWork-ErinCodefolder/PosteriorMeans_Intercepts_rough.png",
       plot = last_plot(), width = 8, height = 6.5, units = c("in")
)

#better plot
post_intercepts <- summary(final_fit1, pars = c('theta0', 'theta_year'),
                           use_cache = FALSE)$summary %>%
  as.data.frame() %>%
  select(1, 4, 8)
names(post_intercepts) <- c('p_mean', 'p_low', 'p_upp')
post_intercepts$param <- c('Time1', 'Time2 Adj.', 'Time3 Adj.', 'Time4 Adj.')
post_intercepts$group <- c('Overall', 'Adjustments', 'Adjustments', 'Adjustments')
post_intercepts$group <- factor(post_intercepts$group, levels = c('Overall', 'Adjustments'))

ggplot(post_intercepts, aes(x = factor(param), y = p_mean, color = factor(param))) +
  geom_errorbar(aes(ymin = p_low, ymax = p_upp), width = 0, show.legend = FALSE) +
  geom_point(show.legend = FALSE) +
  theme_bw() +
  facet_grid(~group, scales = 'free_x', space = 'free_x') +
  geom_hline(yintercept = 0, lty = 2) +
  xlab('Parameter') +
  ylab('Posterior Summary (logit scale)') +
  ggtitle('Infection Probability Intercept Parameters')
ggsave("WilsonsWork-ErinCodefolder/PosteriorMeans_Intercepts_better.png",
       plot = last_plot(), width = 8, height = 6.5, units = c("in")
)
