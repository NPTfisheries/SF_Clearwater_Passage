# Purpose: Gather data for SF Clearwater related to:
#   1. the potential velocity barrier evaluation,
#   2. run-timing.
# 
# Authors: Mike Ackerman and Ryan N. Kinzer 
# 
# Created: May 31, 2023
#   Last Modified: May 22, 2024

# clear environment
rm(list = ls())

# load necessary packages
library(tidyverse)
library(PITcleanr)
library(here)
library(magrittr)

#---------------------
# gather PIT-tag data

# get site configuration info from PTAGIS
# config = buildConfig(node_assign = "site")
# save(config, file = here("data/derived_data/config.rda"))
load(here("data/derived_data/config.rda"))

# GRA & south fork clearwater sites
sf_sites = c("SC1",  # rkm 1; These rkms are from PTAGIS and I don't know their accuracy
             "SC2",  # rkm 2   
             "SC3",  # rkm 60
             "SC4",  # rkm 81
             "CRA")  # Crooked River IPTDS

sf_config = config %>%
  filter(node %in% c("GRA", sf_sites))

# function to query DART observations for a given species and year for adults observed at LGR and within SF Clearwater
queryObsDART_spc_yr = function(spc, yr) {
  dart_df = queryObsDART(species = spc,
               loc = "GRA",
               spawn_year = yr) %>%
    # trim to adults observed in SF Clearwater
    group_by(tag_id) %>%
    filter(any(obs_site %in% sf_sites)) %>%
    mutate(species = spc,
           spawn_year = yr)
  
  print(paste0("Observations for ", spc, " adults in spawn year ", yr, " downloaded from DART."))
  return(dart_df)
}

# set species and years
species = c("Chinook", "Steelhead")
years = 2012:2024

# SKIP UNLESS DATA NEEDS TO BE UPDATED: get all Chinook & steelhead observations from DART for adults at 
# GRA and upstream (includes newly and previously tagged fish)

sy = crossing(species, years) %>%
  # incomplete data for SY2024 Chinook
  filter(!(species == "Chinook" & years == 2024))

dart_obs_list = map2(sy$species, sy$years, queryObsDART_spc_yr)
names(dart_obs_list) = paste0(sy$species, "_", sy$years)
save(dart_obs_list, file = here("data/derived_data/dart_observations/sy12-24_dart_obs.rda"))
# glimpse(dart_obs_list[["Steelhead_2024"]]) # view a single spc, yr result

# load dart_obs_list and convert to a data frame (rbindlist avoids issues with differing data types)
load(here("data/derived_data/dart_observations/sy12-24_dart_obs.rda"))
dart_obs_df = data.table::rbindlist(dart_obs_list)

# write out tag lists for each species and spawn year
for( spc in species ) {
  # set species code
  if(spc == "Chinook")   { spc_code = "chnk" }
  if(spc == "Steelhead") { spc_code = "sthd" }
  for ( yr in years ) {
    tag_list = dart_obs_df %>%
      filter(species == spc,
             spawn_year == yr) %>%
      select(tag_id) %>%
      distinct() %>%
      pull() %>%
      write_lines(paste0(here("output/tag_lists"), "/", spc_code, "_sy", yr, "_tag_list.txt"))
  } # end years loop
} # end species loop

# now query CTHs for tags in PTAGIS

# build parent-child table for SF Clearwater
parent_child = tribble(~"parent", ~"child",
                       #"GRA", "SC1",
                       "SC1", "SC2",
                       "SC2", "SC3",
                       "SC3", "SC4",
                       "SC4", "CRA")

# function to compress and filter detections for a given species and year
compressSpcYr = function(spc, yr) {
  # set species code
  if(spc == "Chinook")   { spc_code = "chnk" ; max_obs_date = paste0(yr, "0915") }
  if(spc == "Steelhead") { spc_code = "sthd" ; max_obs_date = paste0(yr, "0531") }
  
  file_path = here(paste0("data/derived_data/cths/", spc_code, "_sy", yr, ".csv"))
  cth = readCTH(file_path)
  
  comp_df = compress(cth_file = cth,
                     file_type = "PTAGIS",
                     max_minutes = NA,
                     configuration = config,
                     units = "days",
                     ignore_event_vs_release = TRUE) %>%
    mutate(species = spc,
           spawn_year = yr) %>%
    # filter to exclude detections prior to beginning of spawn year
    filter(case_when(
      spc == "Chinook"   ~ min_det >= ymd_hms(paste0(yr, "-03-01 01:00:00")),
      spc == "Steelhead" ~ min_det >= ymd_hms(paste0(yr-1, "-07-01 01:00:00"))
    )) %>%
    filterDetections(parent_child, max_obs_date) %>%
    filter(auto_keep_obs == T) %>%
    select(species,
           spawn_year,
           everything(),
           -direction,
           -auto_keep_obs,
           -user_keep_obs)
  print(paste0("Observations for ", spc, " and spawn year ", yr, " compressed and filtered!"))
  return(comp_df)
} # end compressSpcYr

# compress and filter all observations, combine into a list
comp_list = map2(sy$species, sy$years, compressSpcYr)
names(comp_list) = paste0(sy$species, "_", sy$years)
glimpse(comp_list[["Steelhead_2024"]])

# save compressed and filtered cths
save(comp_list, parent_child, file = here("data/derived_data/cths/sy12-24_compressed_filtered_obs.rda"))

#------------------------
# LGTrapppingDB

# get some biological data from LGR
trap_df = read_csv("C:/Git/SnakeRiverFishStatus/data/LGTrappingDB/LGTrappingDB_2024-05-21.csv")
sf_lgr_df = bind_rows(comp_list) %>%
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

save(sf_lgr_df, file = here("data/derived_data/LGTrappingDB/sf_clearwater_lgtrappingdb.rda"))

#------------------------
# IPTDS Environmental Probe Data

# load necessary libraries
library(fisheR)

# log into BioLogic database to retrieve API token
source("C:/Git/SnakeR_IPTDS/keys/biologic_login.txt")

# set sf clearwater env probe sites and years
env_sites = c("SC1", "SC2", "SC4")
env_years = 2021:2024

env_sites = "SC4"
env_year = 2017

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
        
        save(env_df, file = paste0(here("data/derived_data/enviro"), "/", s, "_", y, ".rda"))
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

# set start and end dates for data retrieval
start_dt = "2021-07-01"
end_dt   = "2024-12-31"

# query stream gage data
sf_elk_gage_info = readNWISsite(13337500)                  # sf clearwater river nr Elk City, ID
sf_elk_daily_cfs = readNWISdv(siteNumbers = 13337500,        
                              parameterCd = "00060",       # mean daily cfs
                              startDate = start_dt, 
                              endDate = end_dt) %>%
  rename(daily_mean_cfs = X_00060_00003)
# Unfortunately, data is only available for the site through 10/17/2021

sf_stites_gage_info = readNWISsite(13338500)                 # sf clearwater river nr Stites, ID
sf_stites_daily_cfs = readNWISdv(siteNumbers = 13338500,        
                                 parameterCd = "00060",    # mean daily cfs
                                 startDate = start_dt, 
                                 endDate = end_dt) %>%
  rename(daily_mean_cfs = X_00060_00003)

# write out stream gage data for analysis
save(sf_elk_gage_info,
     sf_elk_daily_cfs,
     sf_stites_gage_info,
     sf_stites_daily_cfs,
     file = here("data/derived_data/enviro/sf_clearwater_mean_daily_cfs.rda"))

### END SCRIPT