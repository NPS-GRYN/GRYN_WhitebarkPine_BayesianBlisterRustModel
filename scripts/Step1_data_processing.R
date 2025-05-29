### Data Processing ###
## Lydia M. Druin ##
# GRYN WBP Bayesian Blister Rust #

# Stored as "/data_processing.R" #
# Updated 12/17/24 #

### Purpose ###
# This code is for processing raw data to be used in the Bayesian hierarchical
# occupancy model described in Shanahan et al. 2021. Raw data is accessible
# in the "data/raw_data" folder.

### Note to User ###
# If you have already completed all data processing steps within this code and want
# to run/examine models, open the "scripts/bayesian_blisterrust.R" file instead.

# For incorporating new data, be sure the input filenames and exported filenames
# are updated with your initials/date.

## Note ##
# This code was originally written by Wilson Wright (11/19/19) and adapted by LMD.

# Set-up ####

### Clean Up Workspace ###
rm(list = ls()) # empty environment/loaded items in workspace

gc() # release memory

### Load Libraries ###

install.packages("tidyverse") # only need to run this the first time
install.packages("santoku") # only need to run this the first time
# Installing packages every session won't hurt, however!

library(tidyverse)
library(santoku)

options(scipen = 999) # trick to force R to stay in numeric notation
options(dplyr.width = Inf, dplyr.print_min = 100)

# Step 1. Load in Data ####

## Part A. Load Source File ####

# Load in Survey Data #

sourcefile <- "GYE_WBP_AllObservers_CompleteCertifiedDataSet_ExtractFromMasterDatabase_July_23_2024_DO_NOT_DISTRIBUTE_LMDedits.csv"
# Note to user: If you are incorporating post-2023 data, this is a point where
# you need to change the above file name to *your* new file name! Make sure it is
# saved as a .csv file, and copy+paste your file name into the above line while
# retaining the quotation marks around your new file name!

# IMPORTANT INFO: Take a look at the "Raw Data Editing Log" in "data/raw_data" file.
# This file has descriptions of trees that we made changes to in 2024 in the raw data.
# You may need to make these changes yourself! Reach out to R Daley for further info.

path <- "data/raw_data"
sourcepath <- paste(path, sourcefile, sep = "/")
gye <- read.csv(sourcepath,
  header = TRUE, na.strings = c("-999", "-1998", "n/a", "na", "NA", "nodata")
)

glimpse(gye) # glimpse at the data--make sure it is what you expect!

# Load in Observer Priority Data #

sourcefile_Obs <- "GYE_WBP_ObserverPrecedence_byYear.csv"
# Note to user: This file is valid through 2023 and will need updated to incorporate
# observers in 2024-beyond.

sourcepath_Obs <- paste(path, sourcefile_Obs, sep = "/")
observer_priority <- read.csv(sourcepath_Obs, header = TRUE)

# Load in Transect and Stand Data #
# Note to user: these data DO NOT need to be updated UNLESS a new stand or transect is added
transect_data <- read.csv("data/raw_data/transect_data.csv")
stand_area <- read.csv("data/raw_data/StandArea.csv")

## Part B. Exclude Categories ####

# Remove Limber Pine #
with(gye, table(TreeSpecies))
gye <- subset(gye, TreeSpecies != "PIFL")
gye <- droplevels(gye)

# Remove BLM plots #
with(gye, table(AdminOrg))
gye <- subset(gye, AdminOrg != "BLM")
gye <- droplevels(gye)

## Part C. Create Unique Identifier for Transects and Trees ####

gye$Transect_unq <- as.numeric(with(gye, paste(Stand_ID, Transect_ID, sep = ".")))
gye$Tree_unq <- with(gye, paste(Stand_ID, Transect_ID, TreeID, sep = "."))

## Part D. Create Year Column for Year Tree Entered Monitoring ####

TreeAddDate2 <- strptime(as.character(gye$TreeAddDate), format = "%d-%b-%y")
gye$TreeAddYear <- format(TreeAddDate2, "%Y")

## Part E. Create Sampling Period Column ####

gye$SamplingPeriod <- chop_width(gye$SurveyYear,
  width = 4, labels =
    lbl_seq(start = "1")
)

with(gye, table(SurveyYear, SamplingPeriod))
# Look at the resulting table and confirm it is what you expect!

# Step 2. Tidy Data ####

## Part A. Remove Exception and Dead Trees ####

with(gye, table(TreeStatus))
idx <- subset(gye, TreeStatus != "L")
# Write these trees and their data to a csv file if you are interested:
# Note to user: update initials/date to your own
write.csv(idx, "data/review/Exception_and_Dead_Trees_LMD_Dec2024.csv", row.names = FALSE)

# Exception trees are not excluded from the overall data--this is because they
# are monitored before and/or after exception (meaning, unable to be found during
# a given monitoring visit, but are found and surveyed during other visits).
# For blister rust purposes, because a tree was unable to be surveyed during an
# "exception" visit, we will treat these surveys as "no data" and censor them.
# We also exclude all surveys of dead trees--only living tree surveys record
# data on blister rust infections.

gye <- subset(gye, TreeStatus == "L")
gye <- droplevels(gye)

## Part B. Exclude Beetle-Only Surveys ####

# Beetle-only surveys do not collect any data on blister rust.

with(gye, table(SurveyType))
gye <- subset(gye, SurveyType != "pine beetle only")
gye <- droplevels(gye)

## Part C. Create Binary Variables ####

### Beetle Infestation ###

with(gye, table(MtnPineBeetle, useNA = "always"))

gye$BeetlePresence <- ifelse(gye$MtnPineBeetle == "Y", 1, 0)

# Verify:
with(gye, table(BeetlePresence, MtnPineBeetle, useNA = "always"))
# All "BeetlePresence = 0" records are under "MtnPineBeetle = N"? ALL GOOD!

### Blister Rust Infection ###

with(gye, table(InfectionPresent, useNA = "always"))

gye$RustPresence <- ifelse(gye$InfectionPresent == "yes", 1, 0)

# Verify:
with(gye, table(RustPresence, InfectionPresent, useNA = "always"))
# All "RustPresence = 0" records are under "InfectionPresent = no"? ALL GOOD!

# Step 4. Refine Data ####

## Part A. Partition Out Multi-Observer Data ####
# This part is not necessary and only if you are interested in it, or want to
# export the data.

gye.multi <- subset(gye, MultiObserver == "yes")
gye.single <- subset(gye, MultiObserver == "no")

# Export these data if interested:
# Note to user: be sure to change file names to your initials/date!
write.csv(gye.single, "data/processed_data/SingleObserverData_LMDDec2024.csv", row.names = FALSE)
write.csv(gye.multi, "data/processed_data/MultiObserverData_LMDDec2024.csv", row.names = FALSE)

## Part B. Add Observer Priority to Data ####

# Tidy up the observer_priority dataframe
observer_priority$obs_year <- with(observer_priority, paste(LastName, Year, sep = "."))
priority <- observer_priority %>%
  select(
    obs_priority = Priority_Order,
    obs_year = obs_year
  )
glimpse(priority)

# Create the obs_year column in gye dataframe
gye$obs_year <- with(gye, paste(Observer, SurveyYear, sep = "."))

# Add the obs_priority column to gye so we know the order of priority!
gye <- left_join(gye, priority)
gye <- select(gye, -obs_year) # Remove extra column
glimpse(gye) # Take a look at the new column we just added!

## Part C. Select Columns to Keep and Rename ####

gye.tidy <- gye %>%
  select(
    siteID = Transect_unq,
    treeID = TreeID,
    year = SurveyYear,
    date = SurveyDate,
    time = SamplingPeriod,
    obs_name = Observer,
    obs_priority = obs_priority,
    dbh1 = EarliestDBH,
    height1 = EarliestHeightClass,
    dbh2 = LatestDBH,
    height2 = LatestHeightClass,
    meas_year1 = EarliestDBH_Height_Year,
    meas_year2 = LatestDBH_Height_Year,
    inf_ind = RustPresence,
    beetle_ind = BeetlePresence
  )

## Part D. Arrange Observations and Generate ID Numbers ####

gye.tidy <- arrange(gye.tidy, year, siteID, treeID, obs_priority)
gye.tidy[1:10, ] # look at first 10 records

# Create Stand ID #
# This creates a numeric number of each stand in order from 1 to 150
gye.tidy$standID <- gye.tidy$siteID %>%
  floor() %>%
  factor() %>%
  as.numeric()

# Create Transect [within stand] ID #
# This does the same as above, but of each transect within a stand, in order
# from 1 to 2 (generally each stand has no more than 2 transects within it)
gye.tidy$transectID_temp <- (gye.tidy$siteID - floor(gye.tidy$siteID)) * 10
gye.tidy <- gye.tidy %>%
  group_by(standID) %>%
  mutate(transectID = as.numeric(factor(transectID_temp))) %>%
  select(-transectID_temp) %>%
  ungroup()

# Create Observer ID #
# This does the same thing as above, but creates numeric indicators of each observer
gye.tidy$obsID <- gye.tidy$obs_name %>%
  factor() %>%
  as.numeric()

# Check Work #
gye.tidy %>%
  group_by(standID, transectID) %>%
  summarize(siteID = first(siteID)) %>%
  ungroup() %>%
  arrange(standID, transectID) %>%
  glimpse()
# Is standID and transectID in numerical order, from 1-n? Yes--ALL GOOD!

# Step 5. Process Transect and Stand Data ####

## Part A. Process Transect Data ####

# Create siteID variable to match gye data #
transect_data$siteID <- with(transect_data, paste(Stand_ID, Transect_ID, sep = ".")) %>%
  as.numeric()

# Create standID variable to match gye data #
transect_data$standID <- transect_data$siteID %>%
  floor() %>%
  factor() %>%
  as.numeric()

# Create transectID variable to match gye data #
transect_data$transectID_temp <- (transect_data$siteID -
  floor(transect_data$siteID)) * 10
transect_data <- transect_data %>%
  group_by(standID) %>%
  mutate(transectID = as.numeric(factor(transectID_temp))) %>%
  select(-transectID_temp) %>%
  ungroup()

# Check Work #
glimpse(arrange(transect_data, standID, transectID))
# Is standID and transectID in numerical order, from 1-n? Yes--ALL GOOD!

all.equal(
  arrange(transect_data, standID, transectID) %>%
    select(standID, transectID, siteID),
  gye.tidy %>%
    group_by(standID, transectID) %>%
    summarize(siteID = first(siteID)) %>%
    ungroup() %>%
    arrange(standID, transectID) %>%
    select(standID, transectID, siteID)
)
# Equals true? ALL GOOD!

## Part B. Process Stand Data ####

# Select Columns of Interest #
stand_area <- select(stand_area, Stand_ID, Transect_ID, area = StandAreaSQm)

# Create siteID variable to match gye and transect data #
stand_area$siteID <- with(stand_area, paste(Stand_ID, Transect_ID, sep = ".")) %>%
  as.numeric()

# Create standID variable to match gye and transect data #
stand_area$standID <- stand_area$siteID %>%
  floor() %>%
  factor() %>%
  as.numeric()

# Create transectID variable to match gye and transect data #
stand_area$transectID_temp <- (stand_area$siteID -
  floor(stand_area$siteID)) * 10
stand_area <- stand_area %>%
  group_by(standID) %>%
  mutate(transectID = as.numeric(factor(transectID_temp))) %>%
  select(-transectID_temp) %>%
  ungroup()

# Tidy Stand Data #
stand_area <- select(stand_area, -Stand_ID, -Transect_ID)

## Part C. Format Site Data ####

# Combine Site Data #
transect_data <- left_join(transect_data, stand_area,
  by = c("siteID", "standID", "transectID")
)

# Sort and Tidy Dataframe
# Note from WilsonWright: aspect was dropped because it wasn't found to be important before

transect_data <- transect_data %>%
  select(-Stand_ID, -Transect_ID) %>%
  arrange(standID, transectID) %>%
  select(standID, transectID, siteID, area,
    elev = elev_meters,
    slope = slope_degrees
  )

# Standardize (scale) Covariates #
transect_data$area2 <- as.vector(scale(transect_data$area))
transect_data$elev2 <- as.vector(scale(transect_data$elev))
transect_data$slope2 <- as.vector(scale(transect_data$slope))

# Step 6. Save Processed Data ####

# Note to user: edit initials/date to your own!
saveRDS(gye.tidy, "data/rds/GYE_WBP_Processed_AllObservers_BlisterRustData_LMDDec2024_late.rds")
saveRDS(transect_data, "data/rds/GYE_WBP_Processed_TransectData_LMDDec2024.rds")

# END ####

# Session Info ####
sessionInfo()
# R version 4.3.2 (2023-10-31 ucrt)
# Platform: x86_64-w64-mingw32/x64 (64-bit)
# Running under: Windows 10 x64 (build 19045)

# Matrix products: default


# locale:
# [1] LC_COLLATE=English_United States.utf8  LC_CTYPE=English_United States.utf8
# [3] LC_MONETARY=English_United States.utf8 LC_NUMERIC=C
# [5] LC_TIME=English_United States.utf8

# time zone: America/Chicago
# tzcode source: internal

# attached base packages:
# [1] stats     graphics  grDevices utils     datasets  methods   base

# other attached packages:
# [1] lubridate_1.9.2 forcats_1.0.0   stringr_1.5.1   dplyr_1.1.4     purrr_1.0.1
# [6] readr_2.1.4     tidyr_1.3.1     tibble_3.2.1    ggplot2_3.5.1   tidyverse_2.0.0

# loaded via a namespace (and not attached):
# [1] styler_1.10.3      utf8_1.2.4         generics_0.1.3     stringi_1.7.12
# [5] digest_0.6.33      hms_1.1.3          magrittr_2.0.3     grid_4.3.2
# [9] timechange_0.2.0   R.oo_1.26.0        R.cache_0.16.0     jsonlite_1.8.8
# [13] R.utils_2.12.3     pkgbuild_1.4.4     gridExtra_2.3      fansi_1.0.6
# [17] QuickJSR_1.1.3     scales_1.3.0       codetools_0.2-19   cli_3.6.1
# [21] rlang_1.1.1        R.methodsS3_1.8.2  munsell_0.5.1      withr_3.0.0
# [25] StanHeaders_2.32.6 tools_4.3.2        rstan_2.32.6       inline_0.3.19
# [29] parallel_4.3.2     tzdb_0.4.0         colorspace_2.1-0   vctrs_0.6.5
# [33] R6_2.5.1           matrixStats_1.0.0  stats4_4.3.2       lifecycle_1.0.4
# [37] pkgconfig_2.0.3    RcppParallel_5.1.7 pillar_1.9.0       gtable_0.3.5
# [41] loo_2.7.0          glue_1.6.2         Rcpp_1.0.11        tidyselect_1.2.1
# [45] rstudioapi_0.15.0  farver_2.1.1       labeling_0.4.3     compiler_4.3.2