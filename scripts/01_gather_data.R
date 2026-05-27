# Purpose: Gather data for SF Clearwater related to:
#   1. potential velocity barrier evaluation,
#   2. run-timing.
# 
# Authors: Mike Ackerman and Ryan N. Kinzer 
# 
# Created: May 31, 2023
#   Last Modified: May 26, 2026

# clear environment
rm(list = ls())

# install GitHub packages, if needed
# remotes::install_github("DOI-USGS/dataRetrieval",
#                         build_vignettes = TRUE,
#                         build_opts = c("--no-resave-data",
#                                        "--no-manual"))
# remotes::install_github("ryankinzer/fisheR", ref = "master")

# load necessary packages
library(PITcleanr)
library(tidyverse)
library(magrittr)
library(lubridate)
library(data.table)
library(here)
library(dataRetrieval)
library(fisheR)

#---------
# metadata
update_dttm  = now(tzone = "America/Boise")
current_year = year(update_dttm)
save(update_dttm, file = here("data/derived_data/update_dttm.rda"))

#----------------
# set directories
dart_obs_dir = here("data/derived_data/dart_obs/by_sy")
comp_obs_dir = here("data/derived_data/comp_obs/by_sy")
enviro_dir   = here("data/derived_data/enviro")
lgr_out_dir  = here("data/derived_data/LGTrappingDB")

walk(
  c(dart_obs_dir, comp_obs_dir, enviro_dir, lgr_out_dir),
  dir.create,
  recursive    = TRUE,
  showWarnings = FALSE
)

#-------------------------------------------------
# site configuration
load("C:/Git/SnakeRiverFishStatus/data/configuration_files/site_config_LGR_20260116.rda") ; rm(crb_sites_sf, flowlines, sr_site_pops)

# add nodes to parent-child table
pc_nodes = parent_child %>%
  addParentChildNodes(.,  configuration = configuration)

#---------------------
# gather PIT-tag data
species  = c("Chinook", "Steelhead")

sy_all = crossing(
  species = species,
  spawn_yr = 2012:current_year
) %>%
  filter(!(species == "Chinook" & spawn_yr == current_year & as_date(update_dttm) < ymd(paste0(current_year, "-09-15")))) %>%
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

# to refresh everything:
#species_to_retrieve   = species
#spawn_yrs_to_retrieve = 2012:current_year

# or override, manually
species_to_retrieve   = "Steelhead"
spawn_yrs_to_retrieve = current_year

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
                               sites = c("SC1", "SC2", "SC3", "SC4", "CRA"),
                               nodes = c("SC1", "SC2", "SC3", "SC4_D", "SC4_U", "CRA"),
                               #parent_child = sf_clwr_parent_child,
                               #config_df = config,
                               parent_child = pc_nodes,
                               config_df = configuration
                               ) {
  
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
  
  kelt_sites = c("GRS", "GOA", "LMA", "IHR", "MCN", "JDA", "TDA", "BON")
  
  comp_obs = dart_out$compress_obs %>%
    group_by(tag_code) %>%
    filter(any(node %in% nodes)) %>%
    ungroup() %>%
    mutate(
      species = spc,
      spawn_yr = yr
    ) %>%
    filter(min_det >= min_obs_datetime) %>%
    filterDetections(pc_nodes, max_obs_date) %>%
    # deal with fish later observed as kelts so they don't get tossed out
    { 
      tmp = .
      
      if (spc == "Steelhead") {
        tmp = tmp %>%
          group_by(tag_code) %>%
          mutate(
            last_study_det = max(if_else(node %in% nodes, min_det, as_datetime(NA)), na.rm = TRUE),
            has_kelt_obs   = any(node %in% kelt_sites & min_det > last_study_det),
            auto_keep_obs  = case_when(
              node %in% nodes ~ TRUE,
              has_kelt_obs == TRUE & node %in% kelt_sites & min_det > last_study_det ~ FALSE,
              TRUE ~ auto_keep_obs
            )
          ) %>%
          ungroup() %>%
          select(-last_study_det, -has_kelt_obs)
      } 
      
      tmp
    } %>%
    filter(auto_keep_obs == TRUE) %>%
    # this filter removes any remaining fish determined to spawn outside of sf_clwr_sites
    group_by(tag_code) %>%
    filter(any(node %in% nodes)) %>%
    ungroup() %>%
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
      spc            = species,
      yr             = spawn_yr,
      dart_file_path = dart_file_path,
      comp_file_path = comp_file_path
    )
  }
)

# compile comp_obs objects for below
comp_obs_all = list.files(
  comp_obs_dir,
  pattern    = "_comp_obs\\.rds$",
  full.names = TRUE
) %>%
  map(readRDS) %>%
  bind_rows()

#-------------
# LGTrappingDB
lgr_file = list.files(
  path = "C:/Git/SnakeRiverFishStatus/data/LGTrappingDB/",
  pattern = "^LGTrappingDB_.*\\.csv$",
  full.names = TRUE
) %>%
  sort(decreasing = TRUE) %>%
  first()

message("Using LGTrappingDB file: ", basename(lgr_file))

trap_df = read_csv(lgr_file, show_col_types = F) %>%
  # trim to adults with a spawn year
  filter(
    LGDLifeStage == "RF",
    SpawnYear    != "None"
  ) %>%
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

sf_clwr_lgr_df = comp_obs_all %>%
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

save(
  sf_clwr_lgr_df, 
  file = file.path(lgr_out_dir, "sf_clearwater_lgtrappingdb.rda")
)

#-------------------------------
# IPTDS Environmental Probe Data

# log into BioLogic database to retrieve API token
source("C:/Git/SnakeRiverIPTDS/keys/biologic_login.txt")

# set sf clearwater env probe sites and years to retrieve
env_sites = c("SC1", "SC2", "SC4")
env_years = current_year

# loop to request data from each site
for(s in env_sites) {
  for(y in env_years) {
    
    # log-in each time
    biologic_login(email, password)
    
    # set dates
    begin_dt = paste0(y, "-01-01")
    end_dt   = paste0(y, "-12-31")
    
    tryCatch({
      # pass API token to BioLogic; retrieve site environmental data
      env_df <- get_biologic_data(site     = s,
                                  endpoint = "enviro",
                                  begin_dt = begin_dt,
                                  end_dt   = end_dt)
      # save env_df, if it exists
      if (nrow(env_df) > 0) {
        env_df = env_df %>%
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

# set start and end dates for data retrieval (elk city gage started in August 2002)
gage_start_dt = "2003-01-01"
gage_end_dt   = as.character(as_date(update_dttm) + days(1))

gage_sites = tribble(
  ~site_no,   ~monitoring_location_id,
  "13337500", "USGS-13337500",
  "13338500", "USGS-13338500"
)

sf_gage_meta = gage_sites %>%
  mutate(
    meta = map(site_no, \(x) {
      read_waterdata_monitoring_location(
        monitoring_location_number = x
      ) %>%
        sf::st_drop_geometry() %>%
        mutate(
          across(where(lubridate::is.POSIXt), as.character)
        )
    })
  ) %>%
  pull(meta) %>%
  bind_rows()

sf_gage_df = gage_sites %>%
  mutate(
    daily = map(monitoring_location_id, \(x) {
      read_waterdata_daily(
        monitoring_location_id = x,
        parameter_code = "00060",
        time = c(gage_start_dt, gage_end_dt)
      )
    })
  ) %>%
  pull(daily) %>%
  bind_rows() %>%
  select(
    monitoring_location_id,
    date = time,
    daily_mean_cfs = value,
    approval_status,
    qualifier
  ) %>%
  sf::st_drop_geometry() %>%
  mutate(
    site_no = str_remove(monitoring_location_id, "USGS-")
  )

# write out stream gage data for analysis
save(
  sf_gage_meta,
  sf_gage_df,
  file = file.path(enviro_dir, "sf_clearwater_mean_daily_cfs.rda")
)

message("Data update completed: ", update_dttm)

### END SCRIPT