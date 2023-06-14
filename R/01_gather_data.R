# Purpose: Gather data related to the SF Clearwater potential velocity barrier
# 
# Authors: Ryan N. Kinzer and Mike Ackerman
# 
# Created: May 31, 2023
# Modified: June 14, 2023

# clear environment
rm(list = ls())

# load required packages
list_of_packages = c("tidyverse", "PITcleanr", "here")
lapply(list_of_packages, library, character.only = T)

# get site configuration information from PTAGIS
# config = buildConfig()
# save(config, file = here("data/derived_data/config.rda"))
load(here("data/derived_data/config.rda"))

# get all observations from DART for adults at GRA and upstream
# dat = queryObsDART(species = "Steelhead",
#                    loc = "GRA",
#                    spawn_year = 2021)
# save(dat, file = here("data/derived_data/observations/sy2021_sthd_GRA_obs.rda"))
load(here("data/derived_data/observations/sy2021_sthd_GRA_obs.rda"))

# not certain why RK only had sy2021 above. sy2022 should also be complete. Any use for data prior to sy2021?

# isolate some tags and mark data
tags = dat %>%
  select(tag_id,
         contains("mark_"),
         contains("t_"),
         length,
         contains("rel_"),
         contains("trans_")) %>%
  distinct(tag_id, .keep_all = T)

# species code table to add species names to dat
spp_code = tibble(
  mark_species_name = c("Chinook", "Coho", "Steelhead", "Sockeye"),
  t_species = 1:4
)

# prep for compress
dart_obs = dat %>%
  # add species names
  left_join(spp_code, by = "t_species") %>%
  # recode obs_type
  mutate(obs_type = recode(
    obs_type,
    INT = "Observation",
    MRT = "Mark",
    REC = "Recapture"
  )) %>%
  # re-organize some columns
  select(
    tag_code = tag_id,
    mark_species_name,
    mark_rear_type_name = t_rear_type,
    event_type_name = obs_type,
    event_site_code_value = obs_site,
    event_date_time_value = obs_time,
    antenna_id = coil_id,
    antenna_group_configuration_value = config,
    everything()
  )

# and compress using PITcleanr
comp_obs = compress(dart_obs,
                    configuration = config,
                    ignore_event_vs_release = T) %>%
  left_join(tags, by = c("tag_code" = "tag_id"))

# isolate only SF Clearwater tag/mark data and observations
parent_child <- tribble(
  ~'parent', ~'child',
  #'GRA', 'SC1',
  'SC1', 'SC2',
  'SC2', 'SC3',
  'SC3', 'SC4',
  'SC4', 'CRA'
)

# add the parent_child nodes (i.e., arrays)
pc_nodes <- addParentChildNodes(parent_child, config)

# isolate using all observations, but only SF Clearwater nodes
sfclw_obs <- addDirection(comp_obs, pc_nodes)

sfclw_tags = tags %>%
  filter(tag_id %in%
           (sfclw_obs %>%
           distinct(tag_code) %>%
           pull(tag_code)))



# summarise
tmp = sfclw_obs %>%
  group_by(rel_site, node) %>%
  summarise(n_tags = n_distinct(tag_code))

newsome = sfclw_obs %>%
  filter(rel_site == "NEWSOC")

n_distinct(newsome$tag_code)
# only 1 fish impeded by barrier (~5000 detections), only 2 other fish reached SC3 and were only
# observed 1 time, 4 others had last obs at Sc1
# 15,900 tags, 2600 IDFG, 13300 from NPT - started in 2022 and 2023
8/15000

meadow = sfclw_obs =
  filter(rel_site == "MEAD2C")

sfclw_tags %>%
  mutate(rel_year = year(rel_date)) %>%
  group_by(rel_site, rel_year) %>%
  summarise(n_tags = n_distinct(tag_id))

# END SCRIPT

