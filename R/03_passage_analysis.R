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

#---------------------
# Site Detection Efficiencies

comp_df = bind_rows(comp_list) %>%
  filter(spawn_year %in% 2022:2024)

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
ch_df = bind_rows(ch_list)

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
  mutate(conv_rate = est_tags_at_node / lag(est_tags_at_node),
         conv_rate_se = sqrt((1 / est_tags_at_node) + (1 / lag(est_tags_at_node)))) %>%
  ungroup() %>%
  select(species,
         spawn_year,
         node,
         obs_tags_at_node = tags_at_node,
         eff_est,
         eff_se,
         est_tags_at_node,
         conv_rate,
         conv_rate_se)
  
# river kilometers
site_rkms = tribble(~"site", ~"rkm",
                    "SC1", 1,
                    "SC2", 2,
                    "SC3", 60,
                    "SC4", 81,
                    "CRA", 94)

# --------------------------
# prepare data for analysis

# detection by tag in wide format
comp_filter_wide = comp_filter %>%
  # reset slots to start at 1
  # group_by(tag_code) %>%
  # mutate(slot = 1:n()) %>%
  # ungroup() %>%
  select(species,
         spawn_year,
         tag_code,
         node,
         min_det,
         travel_time) %>%
  # for instances where a tag has two detections at a single site, get first detection
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

# prepare data frame for analysis  
sf_df = rbind(sy2022_chnk_ch,
              sy2023_chnk_ch,
              sy2022_sthd_ch,
              sy2023_sthd_ch) %>%
  select(species, spawn_year, tag_code, cap_hist) %>%
  left_join(comp_filter_wide) %>%
  # attach some useful information from LGRTrappindDB
  left_join(sf_lgr_df %>%
              select(-spawn_year)) %>%
  # attach some useful mark and release information from DART
  left_join(sf_dart_obs %>%
              select(tag_id,
                     mark_site,
                     mark_date,
                     rel_site,
                     rel_date) %>%
              distinct(),
            by = c("tag_code" = "tag_id")) %>%
  # join water level data from the SC4 probe, using date when fish was last detected at SC1, SC2, or SC3
  mutate(tmp = date(pmax(min_det_SC1, min_det_SC2, min_det_SC3, na.rm = T))) %>%
  left_join(sc4_env %>%
              filter(metric == "water_level") %>%
              select(date,
                     sc4_water_level_m = mean),
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

#------------------
# passage by length
bins = seq(floor(min(sf_df$lgr_fl_mm, na.rm = T) / 50) * 50,
           ceiling(max(sf_df$lgr_fl_mm, na.rm = T) / 50) * 50,
           by = 50)
bin_labels = c("350-400",
               "400-450",
               "450-500",
               "500-550",
               "550-600",
               "600-650",
               "650-700",
               "700-750",
               "750-800",
               "800-850",
               "850-900",
               "900-950",
               "950-1000")

pass_by_fl = sf_df %>%
  select(species,
         spawn_year,
         tag_code,
         pass_sc3,
         success,
         lgr_fl_mm) %>%
  filter(!is.na(lgr_fl_mm)) %>%
  # get only fish that at least arrived at SC3
  filter(pass_sc3 == T | success == T) %>%
  mutate(fl_bin = cut(lgr_fl_mm,
                      breaks = bins,
                      labels = bin_labels)) %>%
  group_by(species, fl_bin) %>%
  summarize(n_tags = n(),
            success = sum(success)) %>%
  mutate(conv_rate = success / n_tags)
pass_by_fl

#------------------
# passage by saltwater age
pass_by_sw_age = sf_df %>%
  select(species,
         spawn_year,
         tag_code,
         pass_sc3,
         success,
         bio_scale_final_age) %>%
  # convert scale age to saltwater age
  mutate(sw_age = gsub(".*:","", bio_scale_final_age)) %>%
  mutate(sw_age = case_when(
    sw_age == "A" ~ NA,
    TRUE ~ sw_age)) %>%
  filter(!is.na(sw_age)) %>%
  # get only fish that at least arrived at SC3
  filter(pass_sc3 == T | success == T) %>%
  group_by(species, sw_age) %>%
  summarize(n_tags = n(),
            success = sum(success)) %>%
  mutate(conv_rate = success / n_tags)

#------------------
# fork length by age
fl_by_age = sf_df %>%
  select(species,
         lgr_fl_mm,
         bio_scale_final_age) %>%
  filter(!is.na(lgr_fl_mm), !is.na(bio_scale_final_age)) %>%
  mutate(sw_age = gsub(".*:","", bio_scale_final_age)) %>%
  mutate(sw_age = case_when(
    sw_age == "A" ~ NA,
    TRUE ~ sw_age)) %>%
  filter(!is.na(sw_age)) %>%
  ggplot(aes(x = sw_age,
             y = lgr_fl_mm)) +
  geom_boxplot() +
  theme_classic() +
  labs(x = "Saltwater Age",
       y = "FL (mm)") +
  facet_wrap(~species)

#------------------
# fork length by successful passage
fl_by_success = sf_df %>%
  select(species,
         tag_code,
         lgr_fl_mm,
         pass_sc3,
         success) %>%
  # get only fish that at least arrived at SC3
  filter(pass_sc3 == T | success == T) %>%
  filter(!is.na(lgr_fl_mm)) %>%
  ggplot(aes(x = success,
             y = lgr_fl_mm)) +
  geom_boxplot() +
  theme_classic() +
  labs(x = "Successful Passage to SC4",
       y = "FL (mm)") +
  facet_wrap(~species)

#------------------
# passage by SC4 water level
pass_by_wl = sf_df %>%
  select(species,
         spawn_year,
         tag_code,
         pass_sc3,
         success,
         sc4_water_level_m) %>%
  # get only fish that at least arrived at SC3
  filter(pass_sc3 == T | success == T) %>%
  filter(!is.na(sc4_water_level_m))

wl_bins = seq(floor(min(pass_by_wl$sc4_water_level_m, na.rm = T) / 0.2) * 0.2,
              ceiling(max(pass_by_wl$sc4_water_level_m, na.rm = T) / 0.2) * 0.2,
              by = 0.2)
wl_bin_labels = c("0.2-0.4",
                  "0.4-0.6",
                  "0.6-0.8",
                  "0.8-1.0",
                  "1.0-1.2",
                  "1.2-1.4",
                  "1.4-1.6") 

pass_by_wl2 = pass_by_wl %>%
  mutate(wl_bin = cut(sc4_water_level_m,
                      breaks = wl_bins,
                      labels = wl_bin_labels)) %>%
  group_by(species, wl_bin) %>%
  summarize(n_tags = n(),
            success = sum(success)) %>%
  mutate(conv_rate = success / n_tags)
  
#------------------
# water level by successful passage
wl_by_success = sf_df %>%
  select(species,
         tag_code,
         sc4_water_level_m,
         pass_sc3,
         success) %>%
  # get only fish that at least arrived at SC3
  filter(pass_sc3 == T | success == T) %>%
  filter(!is.na(sc4_water_level_m)) %>%
  ggplot(aes(x = success,
             y = sc4_water_level_m)) +
  geom_boxplot() +
  theme_classic() +
  labs(x = "Successful Passage to SC4",
       y = "Water Level (m) at SC4") +
  facet_wrap(~species)
wl_by_success

#------------------
# passage by Stites discharge
cfs_bins = seq(
  floor(min(sf_df$daily_mean_cfs, na.rm = T) / 1000) * 1000,
  ceiling(max(sf_df$daily_mean_cfs, na.rm = T) / 1000) * 1000,
  by = 1000)
tmp = paste0(cfs_bins, "-", lead(cfs_bins))
cfs_bin_labels = tmp[1:(length(tmp) - 1)]
  
pass_by_cfs = sf_df %>%
  select(species,
         tag_code,
         daily_mean_cfs,
         pass_sc3,
         success) %>%
  # get only fish that at least arrived at SC3
  filter(pass_sc3 == T | success == T) %>%
  filter(!is.na(daily_mean_cfs)) %>%
  mutate(cfs_bin = cut(daily_mean_cfs,
                       breaks = cfs_bins,
                       labels = cfs_bin_labels)) %>%
  group_by(species, cfs_bin) %>%
  summarize(n_tags = n(),
            success = sum(success)) %>%
  mutate(conv_rate = success / n_tags)

#------------------
# passage by release group
library(janitor)
tabyl(sf_df$rel_site)

pass_by_rg = sf_df %>%
  select(species,
         spawn_year,
         tag_code,
         pass_sc3,
         success,
         rel_site) %>%
  mutate(rel_group = case_when(
    rel_site == "CLWRSF" ~ "South Fork Clearwater River",
    rel_site == "MEAD2C" ~ "Meadow Creek",
    rel_site == "NEWSOC" ~ "Newsome Creek",
    TRUE ~ "Other")) %>%
  group_by(species, rel_group) %>%
  summarize(n_tags = n(),
            success = sum(success)) %>%
  mutate(conv_rate = success / n_tags) 
  
#------------------
# simple logistic regression
# fl_rg_data = sf_df %>%
#   filter(!rel_site %in% c("CLEARC", "PRDLD1"),
#          !is.na(lgr_fl_mm)) %>%
#   mutate(rel_group = case_when(
#     rel_site == "BONAFF" ~ "Unknown Release",
#     rel_site == "CLWRSF" ~ "South Fork Clearwater River",
#     rel_site == "COLR2"  ~ "Unknown Release",
#     rel_site == "LGRLDR" ~ "Unknown Release",
#     rel_site == "LGRRBR" ~ "Unknown Release",
#     rel_site == "LGRRRR" ~ "Unknown Release",
#     rel_site == "MEAD2C" ~ "Meadow Creek",
#     rel_site == "NEWSOC" ~ "Newsome Creek",
#     TRUE ~ rel_site)) %>%
#   mutate(pass_SC4 = if_else(pass_SC4 == T, 1, 0)) %>%
#   select(tag_code,
#          lgr_fl_mm,
#          pass_SC4,
#          rel_group)
# 
# # fit the logistic regression model
# model = glm(pass_SC4 ~ lgr_fl_mm + rel_group, data = fl_rg_data, family = binomial())
# 
# # create a sequence of x values for prediction
# x_pred = seq(min(fl_rg_data$lgr_fl_mm), max(fl_rg_data$lgr_fl_mm), length.out = 1000)
# 
# # create a data frame for prediction including groups
# pred_data = expand.grid(x = x_pred, group = unique(fl_rg_data$rel_group)) %>%
#   rename(lgr_fl_mm = x,
#          rel_group = group)
# 
# # predict the probabilities using the logistic regression model
# pred_probs = predict(model, newdata = pred_data, type = "response")
# 
# # combine the predicted data with the original data
# plot_data = cbind(pred_data, pred_prob = pred_probs)
# 
# # Plot the logistic regression curve and data points
# lr_p = ggplot() +
#   geom_point(aes(x = lgr_fl_mm,
#                  y = pass_SC4,
#                  color = rel_group),
#              data = fl_rg_data) +
#   geom_line(aes(x = lgr_fl_mm,
#                 y = pred_prob,
#                 color = rel_group),
#             data = plot_data,
#             linewidth = 0.5) +
#   labs(x = "FL (mm)",
#        y = "p(Successful Pass)",
#        color = "Release Group") +
#   theme_minimal()
# lr_p
# 
# # END SCRIPT
