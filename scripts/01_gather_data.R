# Purpose: Gather data for SF Clearwater related to:
#   1. potential velocity barrier evaluation,
#   2. run-timing.
# 
# Authors: Mike Ackerman and Ryan N. Kinzer 
# 
# Created: May 31, 2023
#   Last Modified: May 12, 2026

# clear environment
rm(list = ls())

# load necessary packages
library(PITcleanr)
library(tidyverse)
library(magrittr)
library(lubridate)
library(data.table)

#-------------------------------------------------
# retrieve site config info from PTAGIS, if needed

#config = buildConfig(node_assign = "site")
#save(config, file = "data/derived_data/config.rda")
load("data/derived_data/config.rda")

sf_clwr_sites = c(
  "SC1",  # rkm 1; These rkms are from PTAGIS and I don't know their accuracy
  "SC2",  # rkm 2   
  "SC3",  # rkm 60
  "SC4",  # rkm 81
  "CRA"   # Crooked River IPTDS
)

sf_clwr_config = config %>%
  filter(node %in% c("GRA", sf_clwr_sites))

sf_clwr_parent_child = tribble(
  ~parent, ~child,
  "SC1", "SC2",
  "SC2", "SC3",
  "SC3", "SC4",
  "SC4", "CRA"
)

#---------------------
# gather PIT-tag data

# output directories
dart_obs_dir = "data/derived_data/dart_obs/by_sy"
comp_obs_dir = "data/derived_data/comp_obs/by_sy"
dir.create(dart_obs_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(comp_obs_dir, recursive = TRUE, showWarnings = FALSE)

# set all species and spawn years available
species  = c("Chinook", "Steelhead")
spawn_yr = 2012:2026

sy_all = crossing(species, spawn_yr) %>%
  filter(!(species == "Chinook" & spawn_yr == 2026)) %>%
  mutate(
    spc_code = case_when(
      species == "Chinook"   ~ "chnk",
      species == "Steelhead" ~ "sthd"
    ),
    dart_file_name = paste0(spc_code, "_sy", spawn_yr, "_dart_obs.rds"),
    comp_file_name = paste0(spc_code, "_sy", spawn_yr, "_comp_obs.rds"),
    dart_file_path = file.path(dart_obs_dir, dart_file_name),
    comp_file_path = file.path(comp_obs_dir, comp_file_name)
  )

# set only species and spawn years to update/retrieve
#species_to_retrieve   = species
#spawn_yrs_to_retrieve = spawn_yr
species_to_retrieve   = "Steelhead"
spawn_yrs_to_retrieve = 2026

sy_retrieve = sy_all %>%
  filter(
    species  %in% species_to_retrieve,
    spawn_yr %in% spawn_yrs_to_retrieve
  )

# helper function: query DART, compress, filter, and save
compressDART_spc_yr = function(spc,
                               yr,
                               dart_file_path,
                               comp_file_path,
                               sites = sf_clwr_sites,
                               parent_child = sf_clwr_parent_child,
                               config_df = config) {
  
  max_obs_date = case_when(
    spc == "Chinook"   ~ paste0(yr, "0915"),
    spc == "Steelhead" ~ paste0(yr, "0531")
  )
  
  min_obs_datetime = case_when(
    spc == "Chinook"   ~ ymd_hms(paste0(yr,     "-03-01 01:00:00")),
    spc == "Steelhead" ~ ymd_hms(paste0(yr - 1, "-07-01 01:00:00"))
  )
  
  message("Retrieving and compressing observations for ", spc, ", spawn year ", yr)
  
  dart_out = compressDART(
    species                 = spc,
    loc                     = "GRA",
    spawn_year              = yr,
    configuration           = config_df,
    max_minutes             = NA,
    units                   = "days",
    ignore_event_vs_release = TRUE
  )
  
  dart_obs = dart_out$dart_obs %>%
    group_by(tag_code) %>%
    filter(any(event_site_code_value %in% sites)) %>%
    ungroup() %>%
    mutate(
      species    = spc,
      spawn_year = yr
    )
  
  saveRDS(dart_obs, dart_file_path)
  
  comp_obs = dart_out$compress_obs %>%
    group_by(tag_code) %>%
    filter(any(node %in% sites)) %>%
    ungroup() %>%
    mutate(
      species = spc,
      spawn_yr = yr
    ) %>%
    filter(min_det >= min_obs_datetime) %>%
    filterDetections(parent_child, max_obs_date) %>%
    filter(auto_keep_obs == TRUE) %>%
    select(
      species,
      spawn_yr,
      everything(),
      -direction,
      -auto_keep_obs,
      -user_keep_obs
    )
  
  saveRDS(comp_obs, comp_file_path)
  
  message("Saved raw DART observations: ",   basename(dart_file_path))
  message("Saved compressed observations: ", basename(comp_file_path))
  
  invisible(comp_obs)
}

# run compressDART_spc_yr over sy_retrieve
pwalk(
  sy_retrieve,
  \(species,
    spawn_yr,
    spc_code,
    dart_file_name,
    comp_file_name,
    dart_file_path,
    comp_file_path) {
    
    compressDART_spc_yr(
      spc = species,
      yr = spawn_yr,
      dart_file_path = dart_file_path,
      comp_file_path = comp_file_path
    )
  }
)

#------------------------
# LGTrapppingDB

# get some biological data from LGR
trap_df = read_csv("C:/Git/SnakeRiverFishStatus/data/LGTrappingDB/LGTrappingDB_2026-05-05.csv") %>%
  mutate(
    species = case_when(
      str_starts(SRR, "1") ~ "Chinook",
      str_starts(SRR, "3") ~ "Steelhead",
      TRUE ~ NA_character_
    ),
    spawn_yr = str_remove(SpawnYear, "^SY") %>% as.integer(),
    tag_code = LGDNumPIT
  ) %>%
  filter(species %in% c("Chinook", "Steelhead")) %>%
  # return just first record in case of multiple for given species, spawn_yr, and tag_code
  group_by(species, spawn_yr, tag_code) %>%
  slice(1) %>%
  ungroup()

sf_clwr_lgr_df = comp_obs_df %>%
  select(species,
         spawn_yr,
         tag_code) %>%
  distinct() %>%
  left_join(
    trap_df,
    by = c("species", "spawn_yr", "tag_code"),
    relationship = "many-to-one"
  ) %>%
  select(species,
         spawn_yr,
         tag_code,
         lgr_collection_date = CollectionDate,
         srr = SRR,
         lgr_fl_mm = LGDFLmm,
         gen_sex = GenSex,
         pbt_by_hat = GenPBT_ByHat,
         pbt_rel_group = GenPBT_RGroup,
         bio_scale_final_age = BioScaleFinalAge,
         lgd_mark_ad = LGDMarkAD)

save(sf_clwr_lgr_df, file = "data/derived_data/LGTrappingDB/sf_clearwater_lgtrappingdb.rda")

#------------------------
# IPTDS Environmental Probe Data

# load necessary libraries
library(fisheR)

# log into BioLogic database to retrieve API token
source("C:/Git/SnakeRiverIPTDS/keys/biologic_login.txt")

# set sf clearwater env probe sites and years
env_sites = c("SC1", "SC2", "SC4")
env_years = 2025:2026

# loop to request data from each site
for(s in env_sites) {
  for(y in env_years) {
    
    # log-in each time
    biologic_login(email, password)
    
    # set dates
    begin_dt <- paste0(y, "-01-01")
    end_dt   <- paste0(y, "-12-31")
    
    tryCatch({
      # pass API token to BioLogic; retrieve site environmental data
      env_df <- get_biologic_data(site = s,
                                  endpoint = "enviro",
                                  begin_dt = begin_dt,
                                  end_dt = end_dt)
      # save env_df, if it exists
      if (nrow(env_df) > 0) {
        env_df <- env_df %>%
          select(reader.site.slug,
                 parameter.slug,
                 parameter.units,
                 read_at,
                 value)
        
        save(env_df, file = paste0("data/derived_data/enviro", "/", s, "_", y, ".rda"))
        print(paste0("Environmental data saved for site ", s, ", year ", y, "."))
      } else {
        print(paste0("Environmental data does not exist for site ", s, ", year ", y, "."))
      }
    }, error = function(e) {
      # handle the error i.e., print an error message
      print(paste0("Error occurred for site ", s, ": ", conditionMessage(e)))
    })
  } # end env_years loop
} # end env_sites loop

#------------------------
# Stream Gage Data

# load necessary packages
# library(remotes)
# install_github("DOI-USGS/dataRetrieval",
#                build_vignettes = TRUE,
#                build_opts = c("--no-resave-data",
#                               "--no-manual"))
library(dataRetrieval)

# set start and end dates for data retrieval (elk city gage started in August 2002)
start_dt = "2003-01-01"
end_dt   = "2026-05-13"

# query stream gage data
sf_elk_gage_info = readNWISsite(13337500)                  # sf clearwater river nr Elk City, ID
sf_elk_daily_cfs = readNWISdv(siteNumbers = 13337500,        
                              parameterCd = "00060",       # mean daily cfs
                              startDate = start_dt, 
                              endDate = end_dt) %>%
  rename(daily_mean_cfs = X_00060_00003)
# unfortunately, data is only available for the site through 10/17/2021

sf_stites_gage_info = readNWISsite(13338500)                 # sf clearwater river nr Stites, ID
sf_stites_daily_cfs = readNWISdv(siteNumbers = 13338500,        
                                 parameterCd = "00060",    # mean daily cfs
                                 startDate = start_dt, 
                                 endDate = end_dt) %>%
  rename(daily_mean_cfs = X_00060_00003)

# merge items togeter
sf_gage_meta = rbind(sf_elk_gage_info, sf_stites_gage_info)
sf_gage_df = bind_rows(sf_elk_daily_cfs, sf_stites_daily_cfs)

# write out stream gage data for analysis
save(sf_gage_meta,
     sf_gage_df,
     file = "data/derived_data/enviro/sf_clearwater_mean_daily_cfs.rda")

### END SCRIPT