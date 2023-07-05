# Purpose: Analyze data related to the SF Clearwater potential velocity barrier
# 
# Authors: Mike Ackerman
# 
# Created: July 5, 2023
#   Modified:

# clear environment
rm(list = ls())

# load necessary packages
library(tidyverse)
library(here)

# load some data
load(here("data/derived_data/sf_cleartwater_passage_data.rda"))

# START HERE

sf_df = sf_comp_cth %>%
  mutate(site = str_sub(node,
                        start = 1,
                        end = 3))

# END SCRIPT
