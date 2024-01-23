# Purpose: Gather data related to the SF Clearwater potential velocity barrier
# 
# Authors: Mike Ackerman and Ryan N. Kinzer 
# 
# Created: May 31, 2023
#   Last Modified: January 23, 2024

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
# chinook, sy2022
chnk_2022_dart = queryObsDART(species = "Chinook",
                              loc = "GRA",
                              spawn_year = 2022) %>%
  mutate(species = "Chinook",
         spawn_year = 2022) %>%
  group_by(tag_id) %>%
  filter(any(obs_site %in% sf_nodes))

# chinook, sy2023
chnk_2023_dart = queryObsDART(species = "Chinook",
                              loc = "GRA",
                              spawn_year = 2023) %>%
  mutate(species = "Chinook",
         spawn_year = 2023) %>%
  group_by(tag_id) %>%
  filter(any(obs_site %in% sf_nodes))

# steelhead, sy2022
sthd_2022_dart = queryObsDART(species = "Steelhead",
                              loc = "GRA",
                              spawn_year = 2022) %>%
  mutate(species = "Steelhead",
         spawn_year = 2022) %>%
  group_by(tag_id) %>%
  filter(any(obs_site %in% sf_nodes))

# steelhead, sy2023
sthd_2023_dart = queryObsDART(species = "Steelhead",
                              loc = "GRA",
                              spawn_year = 2023) %>%
  mutate(species = "Steelhead",
         spawn_year = 2023) %>%
  group_by(tag_id) %>%
  filter(any(obs_site %in% sf_nodes))

# combine SY2022 and SY2023, both species, and only those fish observed at a sf_node
sf_dart_obs = rbind(sthd_2022_dart, 
                    sthd_2023_dart,
                    chnk_2022_dart,
                    chnk_2023_dart)
rm(sthd_2022_dart, sthd_2023_dart, chnk_2022_dart, chnk_2023_dart)

# create tag lists for fish observed at sf_nodes
sy2022_chnk_tags = sf_dart_obs %>%
  filter(species == "Chinook",
         spawn_year == 2022) %>%
  select(tag_id) %>%
  distinct() %>%
  pull() %>%
  write_lines(here("data/derived_data/tag_lists/sy2022_chnk_tag_list.txt"))

sy2023_chnk_tags = sf_dart_obs %>%
  filter(species == "Chinook",
         spawn_year == 2023) %>%
  select(tag_id) %>%
  distinct() %>%
  pull() %>%
  write_lines(here("data/derived_data/tag_lists/sy2023_chnk_tag_list.txt"))

sy2022_sthd_tags = sf_dart_obs %>%
  filter(species == "Steelhead",
         spawn_year == 2022) %>%
  select(tag_id) %>%
  distinct() %>%
  pull() %>%
  write_lines(here("data/derived_data/tag_lists/sy2022_sthd_tag_list.txt"))

sy2023_sthd_tags = sf_dart_obs %>%
  filter(species == "Steelhead",
         spawn_year == 2023) %>%
  select(tag_id) %>%
  distinct() %>%
  pull() %>%
  write_lines(here("data/derived_data/tag_lists/sy2023_sthd_tag_list.txt"))

# now query CTHs for sf_tags in PTAGIS

# read in CTH file(s)
sy2022_chnk_cth = readCTH(here("data/derived_data/cths/sy2022_chnk_cth.csv"))
sy2023_chnk_cth = readCTH(here("data/derived_data/cths/sy2023_chnk_cth.csv"))
sy2022_sthd_cth = readCTH(here("data/derived_data/cths/sy2022_sthd_cth.csv"))
sy2023_sthd_cth = readCTH(here("data/derived_data/cths/sy2023_sthd_cth.csv"))

# compress observations
sy2022_chnk_comp = compress(cth_file = sy2022_chnk_cth,
                            file_type = "PTAGIS",   
                            max_minutes = NA,
                            configuration = config,
                            units = "days",
                            ignore_event_vs_release = TRUE) %>%
  mutate(species = "Chinook",
         spawn_year = 2022)

sy2023_chnk_comp = compress(cth_file = sy2023_chnk_cth,
                            file_type = "PTAGIS",   
                            max_minutes = NA,
                            configuration = config,
                            units = "days",
                            ignore_event_vs_release = TRUE) %>%
  mutate(species = "Chinook",
         spawn_year = 2023)

sy2022_sthd_comp = compress(cth_file = sy2022_sthd_cth,
                            file_type = "PTAGIS",   
                            max_minutes = NA,
                            configuration = config,
                            units = "days",
                            ignore_event_vs_release = TRUE) %>%
  mutate(species = "Steelhead",
         spawn_year = 2022)

sy2023_sthd_comp = compress(cth_file = sy2023_sthd_cth,
                            file_type = "PTAGIS",   
                            max_minutes = NA,
                            configuration = config,
                            units = "days",
                            ignore_event_vs_release = TRUE) %>%
  mutate(species = "Steelhead",
         spawn_year = 2023)

# build SF Clearwater parent-child table
parent_child = tribble(~"parent", ~"child",
                       #"GRA", "SC1",
                       "SC1", "SC2",
                       "SC2", "SC3",
                       "SC3", "SC4",
                       "SC4", "CRA")

# add nodes i.e., arrays
pc_nodes = addParentChildNodes(parent_child, 
                               config)
pc_node_order = buildNodeOrder(parent_child = pc_nodes)

# add directionality and indicate whether each detection should be kept
sy2022_chnk_filter = filterDetections(compress_obs = sy2022_chnk_comp,
                                      parent_child = pc_nodes,
                                      max_obs_date = "20220915") %>%
  filter(auto_keep_obs == T) %>%
  select(species,
         spawn_year,
         everything(),
         -direction,
         -auto_keep_obs,
         -user_keep_obs)

sy2023_chnk_filter = filterDetections(compress_obs = sy2023_chnk_comp,
                                      parent_child = pc_nodes,
                                      max_obs_date = "20230915") %>%
  filter(auto_keep_obs == T) %>%
  select(species,
         spawn_year,
         everything(),
         -direction,
         -auto_keep_obs,
         -user_keep_obs)

sy2022_sthd_filter = filterDetections(compress_obs = sy2022_sthd_comp,
                                      parent_child = pc_nodes,
                                      max_obs_date = "20220531") %>%
  filter(auto_keep_obs == T) %>%
  select(species,
         spawn_year,
         everything(),
         -direction,
         -auto_keep_obs,
         -user_keep_obs)

sy2023_sthd_filter = filterDetections(compress_obs = sy2023_sthd_comp,
                                      parent_child = pc_nodes,
                                      max_obs_date = "20230531") %>%
  filter(auto_keep_obs == T) %>%
  select(species,
         spawn_year,
         everything(),
         -direction,
         -auto_keep_obs,
         -user_keep_obs)

# merge compressed CTHs
comp_filter = bind_rows(sy2022_chnk_filter,
                        sy2023_chnk_filter,
                        sy2022_sthd_filter,
                        sy2023_sthd_filter)

# write comp_filter
write_csv(comp_filter,
          file = "data/derived_data/sf_clearwater_filtered_detections.csv")

# estimate node efficiencies
sy2022_chnk_node_eff = estNodeEff(capHist_proc = sy2022_chnk_filter,
                                  node_order = pc_node_order) %>%
  mutate(species = "Chinook",
         spawn_year = 2022) %>%
  select(species,
         spawn_year,
         everything())

sy2023_chnk_node_eff = estNodeEff(capHist_proc = sy2023_chnk_filter,
                                  node_order = pc_node_order) %>%
  mutate(species = "Chinook",
         spawn_year = 2023) %>%
  select(species,
         spawn_year,
         everything())

sy2022_sthd_node_eff = estNodeEff(capHist_proc = sy2022_sthd_filter,
                                  node_order = pc_node_order) %>%
  mutate(species = "Steelhead",
         spawn_year = 2022) %>%
  select(species,
         spawn_year,
         everything())

sy2023_sthd_node_eff = estNodeEff(capHist_proc = sy2023_sthd_filter,
                                  node_order = pc_node_order) %>%
  mutate(species = "Steelhead",
         spawn_year = 2023) %>%
  select(species,
         spawn_year,
         everything())

# bind node efficiencies
node_eff = rbind(sy2022_sthd_node_eff,
                 sy2023_sthd_node_eff,
                 sy2022_chnk_node_eff,
                 sy2023_chnk_node_eff)

# write node efficiencies
write_csv(node_eff,
          file = paste0(here(), "/data/derived_data/sf_clearwater_node_efficiencies.csv"))

# convert filtered cths into capture histories
sy2022_chnk_ch = buildCapHist(filter_ch = sy2022_chnk_filter,
                              parent_child = parent_child,
                              configuration = config,
                              keep_cols = c("tag_code",
                                            "species",
                                            "spawn_year")) %>%
  mutate(species = "Chinook",
         spawn_year = 2022)

sy2023_chnk_ch = buildCapHist(filter_ch = sy2023_chnk_filter,
                              parent_child = parent_child,
                              configuration = config,
                              keep_cols = c("tag_code",
                                            "species",
                                            "spawn_year")) %>%
  mutate(species = "Chinook",
         spawn_year = 2023)

sy2022_sthd_ch = buildCapHist(filter_ch = sy2022_sthd_filter,
                              parent_child = parent_child,
                              configuration = config,
                              keep_cols = c("tag_code",
                                            "species",
                                            "spawn_year")) %>%
  mutate(species = "Steelhead",
         spawn_year = 2022)

sy2023_sthd_ch = buildCapHist(filter_ch = sy2023_sthd_filter,
                              parent_child = parent_child,
                              configuration = config,
                              keep_cols = c("tag_code",
                                            "species",
                                            "spawn_year")) %>%
  mutate(species = "Steelhead",
         spawn_year = 2023)

# define the capture history columns
ch_cols = defineCapHistCols(parent_child = parent_child,
                            configuration = config,
                            use_rkm = TRUE)

# get some biological data from LGR
trap_df = read_csv("C:/Git/SnakeRiverFishStatus/data/LGTrappingDB/LGTrappingDB_2023-11-20.csv")
sf_lgr_df = comp_filter %>%
  select(species,
         spawn_year,
         tag_code) %>%
  distinct() %>%
  # there is one repeat spawner; ignore for now
  #group_by(tag_code) %>%
  #mutate(count = n()) %>%
  #filter(count > 1)
  left_join(trap_df,
            by = c("tag_code" = "LGDNumPIT")) %>%
  select(tag_code,
         spawn_year = SpawnYear,
         lgr_collection_date = CollectionDate,
         srr = SRR,
         lgr_fl_mm = LGDFLmm,
         gen_sex = GenSex,
         pbt_by_hat = GenPBT_ByHat,
         pbt_rel_group = GenPBT_RGroup,
         bio_scale_final_age = BioScaleFinalAge,
         lgd_mark_ad = LGDMarkAD)

# write out objects for analysis
save(sf_lgr_df,
     comp_filter,
     node_eff,
     sy2022_chnk_ch,
     sy2023_chnk_ch,
     sy2022_sthd_ch,
     sy2023_sthd_ch,
     ch_cols,
     file = here("data/derived_data/sf_clearwater_passage_data.rda"))

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
                              endDate = "2023-12-31") %>%
  rename(daily_mean_cfs = X_00060_00003)
# Unfortunately, data is only available for the site through 10/17/2021

sf_stites_gage_info = readNWISsite(13338500)                 # sf clearwater river nr Stites, ID
sf_stites_daily_cfs = readNWISdv(siteNumbers = 13338500,        
                                parameterCd = "00060",    # mean daily cfs
                                startDate = "2021-07-01", 
                                endDate = "2023-12-31") %>%
  rename(daily_mean_cfs = X_00060_00003)
 
# write out stream gage data for analysis
save(sf_elk_gage_info,
     sf_elk_daily_cfs,
     sf_stites_gage_info,
     sf_stites_daily_cfs,
     file = here("data/derived_data/sf_clearwater_mean_daily_cfs.rda"))

### END SCRIPT
