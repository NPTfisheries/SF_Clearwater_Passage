# Purpose: Analyze data related to the potential SF Clearwater velocity barrier
# 
# Authors: Mike Ackerman
# 
# Created: July 5, 2023
#   Modified: May 23, 2024

# clear environment
rm(list = ls())

# load necessary packages
library(tidyverse)
library(here)
library(PITcleanr)
library(janitor)

#---------------------
# Compile Relevant Data

# load compressed and filtered PITcleanr observations
load(here("data/derived_data/cths/sy12-24_compressed_filtered_obs.rda"))

# load dart_obs_list and convert to a data frame (rbindlist avoids issues with differing data types)
load(here("data/derived_data/dart_observations/sy12-24_dart_obs.rda"))
dart_obs_df = data.table::rbindlist(dart_obs_list) ; rm(dart_obs_list)

# load relevant LGTrappingDB data
load(here("data/derived_data/LGTrappingDB/sf_clearwater_lgtrappingdb.rda"))

# load environmental probe data
iptds_env_df = list.files(path = here("data/derived_data/enviro/"), pattern = "^SC.*\\.rda$", full.names = TRUE) %>%
  map_dfr(~ {
    get(load(.))
  })

# load stream gage data
load(here("data/derived_data/enviro/sf_clearwater_mean_daily_cfs.rda"))

# calculate ratios between stream gages to estimate recent elk city gage cfs
flow_bin = 500
flow_ratio = sf_gage_df %>%
  select(site_no,
         Date,
         daily_mean_cfs) %>%
  pivot_wider(names_from = site_no,
              values_from = daily_mean_cfs) %>%
  mutate(ratio = `13338500` / `13337500`,
         bin_13338500 = floor(`13338500` / flow_bin) * flow_bin) %>%
  group_by(bin_13338500) %>%
  summarise(avg_ratio = mean(ratio, na.rm = T)) %>%
  filter(!is.na(avg_ratio))

# create new gage data frame which estimates cfs at elk city gage for recent years
sf_gage_df2 = sf_gage_df %>%
  select(site_no,
         Date,
         daily_mean_cfs) %>%
  pivot_wider(names_from = site_no,
              values_from = daily_mean_cfs) %>%
  mutate(bin_13338500 = floor(`13338500` / flow_bin) * flow_bin) %>%
  left_join(flow_ratio) %>%
  mutate(`13337500` = if_else(is.na(`13337500`), `13338500` / avg_ratio, `13337500`)) %>%
  select(-bin_13338500, -avg_ratio) %>%
  pivot_longer(cols = c(`13337500`, `13338500`),
               names_to = "site_no",
               values_to = "daily_mean_cfs")

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
# Build Capture Histories

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
         # 90% CI; change to 1.96 for 95% CI
         conv_l90ci = pmax(conv_rate - 1.645 * conv_rate_se, 0),
         conv_u90ci = pmin(conv_rate + 1.645 * conv_rate_se, 1)) %>%
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
         conv_l90ci,
         conv_u90ci)
conversion_df

# summary of conversion rates, all sites
conversion_df %>%
  select(species, spawn_year, node, conv_rate, conv_l90ci, conv_u90ci) %>%
  mutate(spawn_year = as.factor(spawn_year)) %>%
  filter(node %in% c("SC2", "SC3", "SC4")) %>%
  ggplot(aes(x = node, y = conv_rate, fill = spawn_year)) +
  geom_col(position = "dodge") +
  geom_errorbar(aes(ymin = conv_l90ci, ymax = conv_u90ci),
                position = position_dodge(width = 0.9),
                width = 0.1) +
  theme_bw() +
  labs(x = "Site",
       y = "Conversion Rate (%)",
       fill = "Spawn Year") +
  facet_wrap(~species, nrow = 2)

# summary of conversion rates SC3 -> SC4
conversion_df %>%
  filter(node == "SC4") %>%
  select(species, spawn_year, node, conv_rate, conv_l90ci, conv_u90ci) %>%
  mutate_at(vars(starts_with("conv")), ~ . * 100) %>%
  mutate(spawn_year = as.factor(spawn_year)) %>%
  ggplot(aes(x = spawn_year,
             y = conv_rate,
             fill = species)) +
  geom_col(position = "dodge") +
  geom_errorbar(aes(ymin = conv_l90ci, ymax = conv_u90ci),
                position = position_dodge(width = 0.9),
                width = 0.1) +
  theme_bw() +
  labs(x = "Spawn Year",
       y = "Conversion Rate (%)",
       fill = "Species",
       title = "Conversion Rate SC3 -> SC4")

# --------------------------
# Data Prep for Further Analysis

# prep water level data from the sc4 probe
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

# compile observation, biological, environmental data, etc.
sf_df = comp_df %>%
  select(species,
         spawn_year,
         tag_code,
         node,
         min_det) %>%
  # for instances where a tag has two compressed observations at a single node, keep first detection
  group_by(species, spawn_year, tag_code, node) %>%
  slice(which.min(min_det)) %>%
  ungroup() %>%
  # pivot wider for the time-being to ease joining some data
  group_by(species, spawn_year, tag_code) %>%
  pivot_wider(names_from = node,
              values_from = c(min_det)) %>%
  ungroup() %>%
  # join capture histories
  left_join(ch_df) %>%
  # join some useful data from LGTrappingDB
  left_join(sf_lgr_df %>%
              select(-spawn_year)) %>%
  # join useful mark and release information from DART 
  left_join(dart_obs_df %>%
              select(tag_id,
                     mark_site,
                     rel_site,
                     t_rear_type) %>%
              distinct(),
            by = c("tag_code" = "tag_id")) %>%
  # join water level at SC4 using latest date that a fish was first observed at SC1, SC2, or SC3
  mutate(tmp_dt = date(pmax(SC1, SC2, SC3, na.rm = T))) %>%
  left_join(sc4_daily_depth %>%
              select(date,
                     sc4_avg_daily_depth_m = mean),
            by = c("tmp_dt" = "date")) %>%
  # join daily average cfs from stites usgs gage 13338500
  left_join(sf_gage_df2 %>%
              filter(site_no == 13338500) %>%
              select(Date,
                     daily_cfs_stites = daily_mean_cfs),
            by = c("tmp_dt" = "Date")) %>%
  # join daily average cfs from elk city usgs gage 13337500
  left_join(sf_gage_df2 %>%
              filter(site_no == 13337500) %>%
              select(Date,
                     daily_cfs_elk = daily_mean_cfs),
            by = c("tmp_dt" = "Date")) %>%
  select(-tmp_dt) %>%
  # did a fish pass SC3, SC4, and/or CRA?
  mutate(pass_sc3 = !is.na(SC3) | !is.na(SC4) | !is.na(CRA),
         success = !is.na(SC4) | !is.na(CRA)) %>%
  # some organizing
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

# summarize tags by species, spawn_year, and release location
rg_summ = sf_df %>%
  select(species,
         spawn_year,
         tag_code,
         rel_site) %>%
  left_join(mrr_df %>%
              select(siteCode,
                     type,
                     rkm),
            by = c("rel_site" = "siteCode")) %>%
  mutate(rel_loc = case_when(
    rel_site %in% c("REDP", "NEWSOC") ~ "Upstream",
    type == "IntraDamReleaseSite"     ~ "Unknown",
    TRUE                              ~ "Downstream")) %>%
  group_by(species, spawn_year, rel_loc) %>%
  summarise(n_tags = n(),
            .groups = "drop") %>%
  pivot_wider(names_from = rel_loc,
              values_from = n_tags,
              values_fill = 0)
rg_summ  

# steelhead passage, upstream vs. downstream releases
pass_by_rel_loc = sf_df %>%
  select(species,
         spawn_year,
         tag_code,
         pass_sc3,
         success,
         rel_site) %>%
  filter(species == "Steelhead") %>%
  left_join(mrr_df %>%
              select(siteCode,
                     type),
            by = c("rel_site" = "siteCode")) %>%
  mutate(rel_loc = case_when(
    rel_site %in% c("REDP", "NEWSOC") ~ "Upstream",
    type == "IntraDamReleaseSite"     ~ "Unknown",
    TRUE                              ~ "Downstream")) %>%
  # get only fish that at least arrived at SC3
  #filter(pass_sc3 == T | success == T) %>%
  left_join(conversion_df %>%
              select(species,
                     spawn_year,
                     node,
                     eff_est,
                     eff_se) %>%
              filter(species == "Steelhead",
                     node %in% c("SC3", "SC4")) %>%
              pivot_wider(names_from = node,
                          values_from = c(eff_est, eff_se))) %>%
  # calculate expansion rates
  mutate(sc3_exp = 1 / eff_est_SC3,
         sc4_exp = 1 / eff_est_SC4,
         sc3_exp_se = eff_se_SC3 / (eff_est_SC3^2),
         sc4_exp_se = eff_se_SC4 / (eff_est_SC4^2)) %>%
  select(-eff_est_SC3, -eff_est_SC4, -eff_se_SC3, -eff_se_SC4) %>%
  # summarise tags to sc3 and sc4
  group_by(rel_loc) %>%
  summarise(n_tags = n(),
            n_tags_sc3 = sum(pass_sc3, na.rm = T),
            n_tags_sc4 = sum(success, na.rm = T),
            est_tags_sc3 = sum(ifelse(pass_sc3, sc3_exp, 0), na.rm = T),
            est_tags_sc4 = sum(ifelse(success, sc4_exp, 0), na.rm = T),
            .groups = "drop") %>%
  # estimated conversion rates SC3 -> SC4
  mutate(conv_rate = est_tags_sc4 / est_tags_sc3,
         conv_rate_se = sqrt( (1 / est_tags_sc3) + (1 / est_tags_sc4) ),
         conv_l90ci = pmax(conv_rate - 1.645 * conv_rate_se, 0),
         conv_u90ci = pmin(conv_rate + 1.645 * conv_rate_se, 1))

rl_summ = pass_by_rel_loc %>%
  mutate_at(vars(starts_with("conv")), ~ . * 100) %>%
  mutate_at(vars(starts_with("conv"), starts_with("est")), ~ round(., 1)) %>%
  select(`Release Location` = rel_loc,
         `Obs Tags` = n_tags,
         `Obs Tags SC3` = n_tags_sc3,
         `Obs Tags SC4` = n_tags_sc4,
         `Est Tags SC3` = est_tags_sc3,
         `Est Tags SC4` = est_tags_sc4,
         `Conversion (%)` = conv_rate)
rl_summ

# --------------------------
# Flow

# passage by cfs at Elk City
pass_by_cfs = sf_df %>%
  select(species,
         spawn_year,
         tag_code,
         pass_sc3,
         success,
         rel_site,
         daily_cfs_elk) %>%
  mutate(cfs_last_det = if_else(daily_cfs_elk >= 600, "abv_600", "blw_600")) %>%
  left_join(mrr_df %>%
              select(siteCode,
                     type),
            by = c("rel_site" = "siteCode")) %>%
  mutate(rel_loc = case_when(
    rel_site %in% c("REDP", "NEWSOC") ~ "Upstream",
    type == "IntraDamReleaseSite"     ~ "Unknown",
    TRUE                              ~ "Downstream")) %>%
  select(-type) %>%
  # join conversion rates
  left_join(conversion_df %>%
              select(species,
                     spawn_year,
                     node,
                     eff_est,
                     eff_se) %>%
              filter(node %in% c("SC3", "SC4")) %>%
              pivot_wider(names_from = node,
                          values_from = c(eff_est, eff_se))) %>%
  # calculate expansion rates
  mutate(sc3_exp = 1 / eff_est_SC3,
         sc4_exp = 1 / eff_est_SC4,
         sc3_exp_se = eff_se_SC3 / (eff_est_SC3^2),
         sc4_exp_se = eff_se_SC4 / (eff_est_SC4^2)) %>%
  select(-eff_est_SC3, -eff_est_SC4, -eff_se_SC3, -eff_se_SC4) %>%
  group_by(species, rel_loc, cfs_last_det) %>%
  #group_by(species, cfs_last_det) %>%
  summarise(n_tags = n(),
            n_tags_sc3 = sum(pass_sc3, na.rm = T),
            n_tags_sc4 = sum(success, na.rm = T),
            est_tags_sc3 = sum(ifelse(pass_sc3, sc3_exp, 0), na.rm = T),
            est_tags_sc4 = sum(ifelse(success, sc4_exp, 0), na.rm = T),
            .groups = "drop") %>%
  # estimated conversion rates SC3 -> SC4
  mutate(conv_rate = pmin(est_tags_sc4 / est_tags_sc3, 1),
         conv_rate_se = sqrt( (1 / est_tags_sc3) + (1 / est_tags_sc4) ),
         conv_l90ci = pmax(conv_rate - 1.645 * conv_rate_se, 0),
         conv_u90ci = pmin(conv_rate + 1.645 * conv_rate_se, 1))

cfs_summ = pass_by_cfs %>%
  mutate_at(vars(starts_with("conv")), ~ . * 100) %>%
  mutate_at(vars(starts_with("conv"), starts_with("est")), ~ round(., 1)) %>%
  select(Species = species,
         `Release Location` = rel_loc,
         `Est CFS (Elk City)` = cfs_last_det,
         `Obs Tags` = n_tags,
         `Obs Tags SC3` = n_tags_sc3,
         `Obs Tags SC4` = n_tags_sc4,
         `Est Tags SC3` = est_tags_sc3,
         `Est Tags SC4` = est_tags_sc4,
         `Conversion (%)` = conv_rate)
cfs_summ

### END SCRIPT
