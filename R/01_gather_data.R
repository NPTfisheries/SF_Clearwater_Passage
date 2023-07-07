# Purpose: Gather data related to the SF Clearwater potential velocity barrier
# 
# Authors: Mike Ackerman and Ryan N. Kinzer 
# 
# Created: May 31, 2023
# Modified: July 5, 2023

# clear environment
rm(list = ls())

# load necessary packages
library(tidyverse)
library(PITcleanr)
library(here)
library(magrittr)

#---------------------
# gather PIT-tag data

# get side configuration info from PTAGIS
# config = buildConfig()
# save(config, file = here("data/derived_data/config.rda"))
load(here("data/derived_data/config.rda"))

# GRA & south fork clearwater sites
sf_nodes = c("SC1",  # rkm 1; These rkms are from PTAGIS and I don't know their accuracy
             "SC2",  # rkm 2   
             "SC3",  # rkm 60
             "SC4",  # rkm 81
             "CRA")  # Crooked River IPTDS

sf_config = config %>%
  filter(node %in% c("GRA", sf_nodes))

# get all observations from DART for adults at GRA and upstream (includes newly and previously tagged fish)
sthd_2022_obs = queryObsDART(species = "Steelhead",
                             loc = "GRA",
                             spawn_year = 2022) %>%
  mutate(spawn_year = 2022) %>%
  group_by(tag_id) %>%
  filter(any(obs_site %in% sf_nodes))

sthd_2023_obs = queryObsDART(species = "Steelhead",
                             loc = "GRA",
                             spawn_year = 2023) %>%
  mutate(spawn_year = 2023) %>%
  group_by(tag_id) %>%
  filter(any(obs_site %in% sf_nodes))

# combine SY2022 and SY2023 and reduce to only those fish that were observed at a sf_node
sf_dart_obs = rbind(sthd_2022_obs, sthd_2023_obs)
rm(sthd_2022_obs, sthd_2023_obs)

# get lists of sy2022 and sy2023 steelhead that have been detected at sf_nodes
sf_2022_tags = sf_dart_obs %>%
  filter(spawn_year == 2022) %>%
  select(tag_id) %>%
  distinct()

sf_2023_tags = sf_dart_obs %>%
  filter(spawn_year == 2023) %>%
  select(tag_id) %>%
  distinct()

# write to .txt for upload to PTAGIS
write.table(sf_2022_tags,
            here("data/derived_data/sfclrwtr_sy2022_sthd_tag_list.txt"),
            quote = F,
            row.names = F,
            col.names = F,
            sep = "\t")

write.table(sf_2023_tags,
            here("data/derived_data/sfclrwtr_sy2023_sthd_tag_list.txt"),
            quote = F,
            row.names = F,
            col.names = F,
            sep = "\t")

# now query CTHs for sf_tags in PTAGIS

# read in CTH file(s)
sf_2022_cth = readCTH(here("data/derived_data/sfclrwtr_sy2022_sthd_cths.csv"))
sf_2023_cth = readCTH(here("data/derived_data/sfclrwtr_sy2023_sthd_cths.csv"))

# compress observations
sf_2022_comp_cth = compress(cth_file = sf_2022_cth,
                            file_type = "PTAGIS",   
                            max_minutes = NA,
                            configuration = config,
                            units = "days",
                            ignore_event_vs_release = TRUE) %>%
  mutate(spawn_year = 2022)

sf_2023_comp_cth = compress(cth_file = sf_2023_cth,
                            file_type = "PTAGIS",   
                            max_minutes = NA,
                            configuration = config,
                            units = "days",
                            ignore_event_vs_release = TRUE) %>%
  mutate(spawn_year = 2023)

# isolate only SF Clearwater cth observations
parent_child = tribble(~"parent", ~"child",
                       #"GRA", "SC1",
                       "SC1", "SC2",
                       "SC2", "SC3",
                       "SC3", "SC4",
                       "SC4", "CRA")

# add nodes i.e., arrays
pc_nodes = addParentChildNodes(parent_child, config)
pc_nodes_order = c("SC1", "SC2", "SC3", "SC4B0", "SC4A0", "CRAB0", "CRAA0")

# add directionality and indicate whether each detection should be kept
sf_2022_comp_cth = filterDetections(compress_obs = sf_2022_comp_cth,
                                    parent_child = pc_nodes,
                                    max_obs_date = NULL)
sf_2023_comp_cth = filterDetections(compress_obs = sf_2023_comp_cth,
                                    parent_child = pc_nodes,
                                    max_obs_date = NULL)

# merge compressed CTHs
sf_comp_cth = bind_rows(sf_2022_comp_cth, sf_2023_comp_cth) %>%
  filter(auto_keep_obs == T) %>%
  select(-direction,
         -auto_keep_obs,
         -user_keep_obs)

# clean environment
rm(sf_2022_comp_cth, sf_2023_comp_cth)

# convert compress ptagis cths into capture histories
sf_ch = sf_comp_cth %>%
  filter(node %in% pc_nodes_order) %>%
  mutate(node = factor(node,
                       levels = pc_nodes_order)) %>%
  select(tag_code, spawn_year, node) %>%
  distinct() %>%
  mutate(seen = 1) %>%
  # convert to wide format
  pivot_wider(names_from = node,
              values_from = seen,
              values_fill = 0,
              names_sort = T,
              names_expand = T) %>% # makes sure all nodes are included in ch's
  # unite(ch, 2:(length(pc_nodes_order) + 1), sep = "")
  left_join(sf_comp_cth %>%
              group_by(tag_code, spawn_year, node) %>%
              summarise(last_det = min(min_det)) %>%
              pivot_wider(names_from = node,
                          values_from = last_det,
                          names_glue = "{node}_first_det")) %>%
  left_join(sf_comp_cth %>%
              group_by(tag_code, spawn_year, node) %>%
              summarise(last_det = max(max_det)) %>%
              pivot_wider(names_from = node,
                          values_from = last_det,
                          names_glue = "{node}_last_det"))

# Get some biological data from LGR
trap_df = read_csv("C:/Git/SnakeRiverFishStatus/data/LGTrappingDB/LGTrappingDB_2023-06-27.csv")
sf_lgr_df = bind_rows(sf_2022_tags, sf_2023_tags) %>%
  distinct() %>%
  left_join(trap_df,
            by = c("tag_id" = "LGDNumPIT")) %>%
  select(tag_id,
         spawn_year = SpawnYear,
         lgr_collection_date = CollectionDate,
         srr = SRR,
         lgr_fl_mm = LGDFLmm,
         gen_sex = GenSex,
         pbt_by_hat = GenPBT_ByHat,
         pbt_rel_group = GenPBT_RGroup,
         bio_scale_final_age = BioScaleFinalAge,
         lgd_mark_ad = LGDMarkAD)

# a few tags w/ multiple records at LGR
sf_dups = sf_lgr_df %>%
  group_by(tag_id) %>%
  mutate(count = n()) %>%
  filter(count > 1) %>%
  slice(1, 5, 7) %>%
  select(-count)

# remove the few fish that have multiple records, then rbind those records we want to keep
sf_lgr_df %<>%
  filter(!tag_id %in% c("3DD.003D4941C2",
                        "3DD.003D82B197",
                        "3DD.003DEA1EBD")) %>%
  rbind(sf_dups)

# write out objects for analysis
save(sf_lgr_df,
     sf_dart_obs,
     sf_comp_cth,
     sf_ch,
     file = here("data/derived_data/sf_cleartwater_passage_data.rda"))

#------------------------
# gather stream gage data

# load necessary libraries
#library(remotes)
# install_github("DOI-USGS/dataRetrieval",
#                build_vignettes = TRUE, 
#                build_opts = c("--no-resave-data",
#                               "--no-manual"))
library(dataRetrieval)

# query stream gage data
sf_elk_gage_info = readNWISsite(13337500)                    # sf clearwater river nr Elk City, ID
sf_elk_daily_cfs = readNWISdv(siteNumbers = 13337500,        
                              parameterCd = "00060",    # mean daily cfs
                              startDate = "2021-07-01", 
                              endDate = "2023-06-30") %>%
  rename(daily_mean_cfs = X_00060_00003)
# Unfortunately, data is only available for the site through 10/17/2021

sf_stites_gage_info = readNWISsite(13338500)                 # sf clearwater river nr Stites, ID
sf_stites_daily_cfs = readNWISdv(siteNumbers = 13338500,        
                                parameterCd = "00060",    # mean daily cfs
                                startDate = "2021-07-01", 
                                endDate = "2023-06-30") %>%
  rename(daily_mean_cfs = X_00060_00003)
 
# write out stream gage data for analysis
save(sf_elk_gage_info,
     sf_elk_daily_cfs,
     sf_stites_gage_info,
     sf_stites_daily_cfs,
     file = here("data/derived_data/sf_clearwater_mean_daily_cfs.rda"))

### END SCRIPT
