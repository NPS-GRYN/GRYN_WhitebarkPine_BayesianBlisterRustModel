### Data Processing ###
## Lydia M. Druin ##
# GRYN WBP Summaries #

# Stored as "/dbh_investigation.R" #
# Updated 07/16/24 #

# TO DO ####

### Purpose ###
# This code is for calculating the change in DBH from earliest measurement to
# latest measurement. We do this to identify changes that don't make sense
# (e.g., huge increases in DBH or declines in DBH).

# Set-up ####

### Clean Up Workspace ###

rm(list = ls()) # empty environment/loaded items in workspace

gc() # release memory

### Load Libraries ###
install.packages("tidyverse") # only need to run this the first time
# Installing packages every session won't hurt, however!

library(tidyverse)

options(scipen = 999) # trick to force R to stay in numeric notation
options(dplyr.width = Inf, dplyr.print_min = 100)

# Step 1. Load in Data ####

## Part A. Load Source File ####

sourcefile <- "GYE_WBP_CompleteCertifiedDataSet_ExtractFromMasterDatabase_Feb_12_2023_DO_NOT_DISTRIBUTE_LMDedits.csv"
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

# Step 2. Identify Change in DBH ####

## Part A. Glance at Data ####

with(gye, table(EarliestDBHcat, LatestDBHcat))

## Part B. Make "DBHchange" column ####

gye$DBHchange = gye$LatestDBH - gye$EarliestDBH

with(gye, table(DBHchange)) # Review Data

## Part C. Identify Nonsensical Changes in DBH ####

DBHloss = subset(gye, DBHchange < 0)
DBH_hugegain = subset(gye, DBHchange > 6)

gye_DBH = rbind(DBHloss, DBH_hugegain)

# Step 3. Export Data for Review ####

write.csv(gye, "data/review/alldata_DBH_investigation.csv", row.names = FALSE)
write.csv(gye_DBH, "data/review/oddballdata_DBH_investigation.csv", row.names = FALSE)
