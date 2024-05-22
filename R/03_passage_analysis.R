# Purpose: Analyze data related to the potential SF Clearwater velocity barrier
# 
# Authors: Mike Ackerman
# 
# Created: July 5, 2023
#   Modified: May 22, 2024

# clear environment
rm(list = ls())

# load necessary packages
library(tidyverse)
library(here)
library(PITcleanr)

#---------------------
# Compile Relevant Data

# load compressed and filtered PITcleanr observations
load(here("data/derived_data/cths/sy12-24_compressed_filtered_obs.rda"))

# load dart_obs_list and convert to a data frame (rbindlist avoids issues with differing data types)
load(here("data/derived_data/dart_observations/sy12-24_dart_obs.rda"))
dart_obs_df = data.table::rbindlist(dart_obs_list) ; rm(dart_obs_list)

# load environmental probe data
iptds_env_df = list.files(path = here("data/derived_data/enviro/"), pattern = "^SC.*\\.rda$", full.names = TRUE) %>%
  map_dfr(~ {
    get(load(.))
  })

# load stream gage data
load(here("data/derived_data/enviro/sf_clearwater_mean_daily_cfs.rda"))

#---------------------
# Site Detection Efficiencies

comp_df = bind_rows(comp_list) %>%
  filter(spawn_year %in% 2022:2024) ; rm(comp_list)

# build site order
sf_site_order = buildNodeOrder(parent_child)

# Function to estimate node efficiencies for a given species and year
estNodeEff_sy = function(spc, yr) {
  # filter data down to species and year
  obs_df = comp_df %>%
    filter(species == spc,
           spawn_year == yr)
  
  # estimate node efficiencies for given species and spawn year
  est_df = estNodeEff(capHist_proc = obs_df,
                      node_order = sf_site_order) %>%
    mutate(species = spc,
           spawn_year = yr) %>%
    select(species, 
           spawn_year,
           everything())
  
  print(paste0("Estimated node efficiencies for ", spc, " in spawn year ", yr))
  return(est_df)
}

# Estimate node efficiencies across species and spawn years
sy = crossing(species = comp_df$species, years = comp_df$spawn_year) %>%
  # incomplete data for SY2024 Chinook
  filter(!(species == "Chinook" & years == 2024))
node_est_list = map2(sy$species, sy$years, estNodeEff_sy)
names(node_est_list) = paste0(sy$species, "_", sy$years)
node_est_df = bind_rows(node_est_list) ; rm(node_est_list)

# write results
write_csv(node_est_df,
          file = paste0(here(), "/output/detection_probs/sf_clearwater_site_efficiencies.csv"))

#---------------------
# Site Detection Efficiencies

# load configuration file
load(here("data/derived_data/config.rda"))

# Function to convert cths into capture histories for a given species and year
buildCapHist_sy = function(spc, yr) {
  # filter data down to species and year
  obs_df = comp_df %>%
    filter(species == spc,
           spawn_year == yr)
  
  ch_df = buildCapHist(filter_ch = obs_df,
                       parent_child = parent_child,
                       configuration = config,
                       keep_cols = c("tag_code")) %>%
    mutate(species = spc,
           spawn_year = yr)
  
  print(paste0("Created capture histories for ", spc, " in spawn year ", yr))
  return(ch_df)
}

# Create capture histories for each species and spawn year
ch_list = map2(sy$species, sy$years, buildCapHist_sy)
names(ch_list) = paste0(sy$species, "_", sy$years)
ch_df = bind_rows(ch_list) ; rm(ch_list)

# define the capture history columns
ch_cols = defineCapHistCols(parent_child = parent_child,
                            configuration = config,
                            use_rkm = TRUE)

# save capture histories
save(ch_df,
     ch_cols,
     file = here("output/capture_histories/sy22-24_capture_histories.rda"))

# --------------------------
# Conversion Rates

# calculate conversion rates with standard error
conversion_df = node_est_df %>%
  group_by(species, spawn_year) %>%
  mutate(conv_rate = pmin(est_tags_at_node / lag(est_tags_at_node), 1),
         conv_rate_se = sqrt((1 / est_tags_at_node) + (1 / lag(est_tags_at_node))),
         conv_l95ci = pmax(conv_rate - 1.96 * conv_rate_se, 0),
         conv_u95ci = pmin(conv_rate + 1.96 * conv_rate_se, 1)) %>%
  ungroup() %>%
  select(species,
         spawn_year,
         node,
         obs_tags_at_node = tags_at_node,
         eff_est,
         eff_se,
         est_tags_at_node,
         conv_rate,
         conv_rate_se,
         conv_l95ci,
         conv_u95ci)
conversion_df

# summary of conversion rates SC3 -> SC4
conversion_df %>%
  filter(node == "SC4") %>%
  select(species, spawn_year, node, conv_rate, conv_l95ci, conv_u95ci) %>%
  mutate_at(vars(starts_with("conv")), ~ . * 100) %>%
  mutate(spawn_year = as.factor(spawn_year)) %>%
  ggplot(aes(x = spawn_year,
             y = conv_rate,
             fill = species)) +
  geom_col(position = "dodge") +
  geom_errorbar(aes(ymin = conv_l95ci, ymax = conv_u95ci),
                position = position_dodge(width = 0.9),
                width = 0.1) +
  theme_classic() +
  labs(x = "Spawn Year",
       y = "Conversion Rate (%)",
       fill = "Species",
       title = "Conversion Rate SC3 -> SC4")

# --------------------------
# Data Prep for Further Analysis

# detection by tag in wide format
comp_df_wide = comp_df %>%
  select(species,
         spawn_year,
         tag_code,
         node,
         min_det,
         travel_time) %>%
  # for instances where a tag has two compressed, get first detection
  group_by(species, spawn_year, tag_code, node) %>%
  slice(which.min(min_det)) %>%
  # need to convert datetime columns to character before pivot_wider
  mutate(min_det = as.character(min_det),
         travel_time = as.character(travel_time)) %>%
  # pivot wider
  group_by(species, spawn_year, tag_code) %>%
  pivot_wider(names_from = node,
              values_from = c(min_det, travel_time)) %>%
  ungroup()

# load relevant LGTrappingDB data
load(here("data/derived_data/LGTrappingDB/sf_clearwater_lgtrappingdb.rda"))

# environmental data from the SC4 probe
sc4_daily_depth = iptds_env_df %>%
  filter(reader.site.slug == "SC4",
         parameter.slug == "water_level") %>%
  mutate(value = as.numeric(value),
         date = date(read_at)) %>%
  group_by(reader.site.slug,
           parameter.slug,
           date) %>%
  summarise(mean = round(mean(value), 2)) %>%
  ungroup()

sf_df = ch_df %>%
  select(species, spawn_year, tag_code, cap_hist) %>%
  left_join(comp_df_wide) %>%
  # attach some useful data from LGTrappingDB
  left_join(sf_lgr_df %>%
              select(-spawn_year)) %>%
  # attach some useful mark and release information from DART
  left_join(dart_obs_df %>%
              select(tag_id,
                     mark_site,
                     mark_date,
                     rel_site,
                     rel_date) %>%
              distinct(),
            by = c("tag_code" = "tag_id")) %>%
  # join water level data from the SC4 probe, using date when fish was last detected at SC1, SC2, or SC3
  mutate(tmp = date(pmax(min_det_SC1, min_det_SC2, min_det_SC3, na.rm = T))) %>%
  left_join(sc4_daily_depth %>%
              filter(reader.site.slug == "SC4") %>%
              select(date,
                     sc4_avg_daily_depth_m = mean),
            by = c("tmp" = "date")) %>%
  # join water level data from stites usgs gage 13338500
  left_join(sf_stites_daily_cfs %>%
              select(Date,
                     daily_mean_cfs),
            by = c("tmp" = "Date")) %>%
  select(-tmp) %>%
  # did fish migrate past the potential velocity barrier?
  mutate(pass_sc3 = !is.na(min_det_SC3),
         success = !is.na(min_det_SC4) | !is.na(min_det_CRA)) %>%
  select(species,
         spawn_year,
         tag_code,
         cap_hist,
         pass_sc3,
         success,
         everything())

# --------------------------
# Release Groups

# query MRR sites
mrr_df = queryMRRMeta()

rg_df = sf_df %>%
  select(species,
         spawn_year,
         tag_code,
         rel_site) %>%
  left_join(mrr_df %>% 
              select(siteCode,
                     type,
                     rkm),
            by = c("rel_site" = "siteCode")) %>%
  mutate(release_location = case_when(
    rel_site == "NEWSOC" ~ "NEWSOC",
    type == "IntraDamReleaseSite" ~ "Unknown",
    TRUE ~ "Downstream")) %>%
  group_by(species, spawn_year, release_location) %>%
  summarise(n_tags = n(),
            .groups = "drop") %>%
  filter(release_location != "Unknown") %>%
  pivot_wider(names_from = release_location,
              values_from = n_tags,
              values_fill = 0)
rg_df


### END SCRIPT
