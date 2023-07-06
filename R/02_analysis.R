# Purpose: Analyze data related to the potential SF Clearwater velocity barrier
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

# environmental data from the SC4 probe
sc4_raw = read_csv(here("data/raw_data/SC4_enviro_f62a18b7-c2a7-4df8-899a-f0f26ba9a7dc.csv"))
sc4_env = sc4_raw %>%
  rename(node = slug...1,
         metric = slug...2) %>%
  mutate(date = date(read_at)) %>%
  group_by(node, metric, date) %>%
  summarise(mean = round(mean(value), 2)) %>%
  ungroup()

# --------------------------
# prepare data for analysis
sf_df = sf_ch %>%
  rowwise() %>%
  mutate(SC1_SC2 = ifelse(any(c(SC1, SC2) == 1), 1, 0),
         SC3_SC4_CRA = ifelse(any(c(SC3, SC4B0, SC4A0, CRAB0, CRAA0) == 1), 1, 0),
         SC1_SC2_last_det = max(c_across(SC1_last_det:SC2_last_det), na.rm = T),
         SC3_SC4_CRA_first_det = min(c_across(SC3_first_det:SC4A0_first_det), na.rm = T),
         travel_time = difftime(SC3_SC4_CRA_first_det, SC1_SC2_last_det, units = "days")) %>%
  ungroup() %>%
  # create a column of capture histories
  unite(ch, SC1:CRAA0, sep = "") %>%
  select(tag_code,
         spawn_year,
         ch,
         SC1_SC2,
         SC3_SC4_CRA,
         SC1_SC2_last_det,
         SC3_SC4_CRA_first_det,
         travel_time) %>%
  # convert any cell containing "Inf" to NA
  mutate_all(~replace(., str_detect(., "Inf"), NA)) %>%
  # attach some useful information from lower granite trap
  left_join(sf_lgr_df %>% 
              select(-spawn_year),
            by = c("tag_code" = "tag_id")) %>%
  # attach some useful mark and release information from DART
  left_join(sf_dart_obs %>%
              select(tag_id,
                     mark_site,
                     mark_date,
                     rel_site,
                     rel_date) %>%
              distinct(),
            by = c("tag_code" = "tag_id")) %>%
  # join water level data from the SC4 probe
  mutate(tmp = date(SC1_SC2_last_det)) %>%
  left_join(sc4_env %>%
              filter(metric == "water_level") %>%
              select(date,
                     water_level_m = mean),
            by = c("tmp" = "date")) %>%
  select(-tmp) %>%
  # was passage successful?
  mutate(success = ifelse(SC3_SC4_CRA == 1, T, F))

# -----------------------
# passage by spawn year
pass_by_sy = sf_df %>%
  filter(SC1_SC2 == 1) %>%
  group_by(spawn_year) %>%
  summarize(n_tags = n(),
            n_tags_passed = sum(success),
            pct_passed = round((n_tags_passed / n_tags) * 100, 1),
            mean_travel_time_days = round(mean(travel_time, na.rm = T), 1),
            med_travel_time_days = round(median(travel_time, na.rm = T), 1)) %>%
  bind_rows(sf_df %>% 
              filter(SC1_SC2 == 1) %>%
              summarize(spawn_year = -9999,
                        n_tags = n(),
                        n_tags_passed = sum(success),
                        pct_passed = round((n_tags_passed / n_tags) * 100, 1),
                        mean_travel_time_days = round(mean(travel_time, na.rm = T), 1),
                        med_travel_time_days = round(median(travel_time, na.rm = T), 1))) %>%
  mutate(spawn_year = replace(spawn_year, spawn_year == -9999, "Totals"))
pass_by_sy

# -----------------------
# passage by release group, year
pass_by_rel_group = sf_df %>%
  filter(SC1_SC2 == 1) %>%
  mutate(rel_group = case_when(
    rel_site == "BONAFF" ~ "Dam Adults",
    rel_site == "CLEARC" ~ "Clear Creek",
    rel_site == "CLWRSF" ~ "South Fork Clearwater River",
    rel_site == "COLR2"  ~ "Columbia River",
    rel_site == "LGRLDR" ~ "Dam Adults",
    rel_site == "LGRRBR" ~ "LGR Juv Barged",
    rel_site == "LGRRRR" ~ "LGR Juv Return-to-River",
    rel_site == "MEAD2C" ~ "Meadow Creek",
    rel_site == "NEWSOC" ~ "Newsome Creek",
    rel_site == "PRDLD1" ~ "Dam Adults",
    TRUE ~ rel_site)) %>%
  group_by(rel_group, spawn_year) %>%
  summarize(n_tags = n(),
            n_tags_passed = sum(success),
            pct_passed = round((n_tags_passed / n_tags) * 100, 1),
            mean_travel_time_days = round(mean(travel_time, na.rm = T), 1),
            med_travel_time_days = round(median(travel_time, na.rm = T), 1)) %>%
  bind_rows(sf_df %>% 
              filter(SC1_SC2 == 1) %>%
              summarize(spawn_year = -9999,
                        n_tags = n(),
                        n_tags_passed = sum(success),
                        pct_passed = round((n_tags_passed / n_tags) * 100, 1),
                        mean_travel_time_days = round(mean(travel_time, na.rm = T), 1),
                        med_travel_time_days = round(median(travel_time, na.rm = T), 1))) %>%
  mutate(spawn_year = replace(spawn_year, spawn_year == -9999, "Totals"))
pass_by_rel_group  

# bar plot
pass_by_rel_group_p = pass_by_rel_group %>%
  select(rel_group, spawn_year, pct_passed) %>%
  filter(spawn_year != "Totals") %>%
  ggplot(aes(x = rel_group,
             y = pct_passed,
             fill = spawn_year)) +
  geom_bar(position = "dodge",
           stat = "identity") +
  theme_classic() +
  labs(x = "Release Group",
       y = "% Passed",
       fill = "Spawn Year") +
  ylim(0, 100) +
  theme(legend.position = "top",
        axis.text.x = element_text(angle = 45, vjust = 0.55))
pass_by_rel_group_p

#------------------
# passage by length

# create fl bins
bins = seq(500, 1000, by = 50)
bin_labels = c("500-550",
               "550-600",
               "600-650",
               "650-700",
               "700-750",
               "750-800",
               "800-850",
               "850-900",
               "900-950",
               "950-1000")

pass_by_length = sf_df %>%
  filter(SC1_SC2 == 1,
         !is.na(lgr_fl_mm)) %>%
  select(tag_code,
         success,
         travel_time,
         lgr_fl_mm) %>%
  mutate(fl_bin = cut(lgr_fl_mm,
                      breaks = bins,
                      labels = bin_labels)) %>%
  group_by(fl_bin) %>%
  summarize(n_tags = n(),
            n_tags_passed = sum(success),
            pct_passed = round((n_tags_passed / n_tags) * 100, 1),
            mean_travel_time_days = round(mean(travel_time, na.rm = T), 1),
            med_travel_time_days = round(median(travel_time, na.rm = T), 1))
pass_by_length

#------------------
# passage by ocean_age
pass_by_age = sf_df %>%
  filter(SC1_SC2 == 1,
         !is.na(bio_scale_final_age)) %>%
  select(tag_code,
         success,
         travel_time,
         bio_scale_final_age) %>%
  mutate(ocean_age = gsub(".*:","", bio_scale_final_age)) %>%
  mutate(ocean_age = case_when(
    ocean_age == "A" ~ NA,
    TRUE ~ ocean_age)) %>%
  filter(!is.na(ocean_age)) %>%
  group_by(ocean_age) %>%
  summarize(n_tags = n(),
            n_tags_passed = sum(success),
            pct_passed = round((n_tags_passed / n_tags) * 100, 1),
            mean_travel_time_days = round(mean(travel_time, na.rm = T), 1),
            med_travel_time_days = round(median(travel_time, na.rm = T), 1))
pass_by_age

fl_by_age = sf_df %>%
  select(lgr_fl_mm,
         bio_scale_final_age) %>%
  filter(!is.na(lgr_fl_mm), !is.na(bio_scale_final_age)) %>%
  mutate(ocean_age = gsub(".*:","", bio_scale_final_age)) %>%
  mutate(ocean_age = case_when(
    ocean_age == "A" ~ NA,
    TRUE ~ ocean_age)) %>%
  filter(!is.na(ocean_age)) %>%
  select(-bio_scale_final_age) %>%
  ggplot(aes(x = ocean_age,
             y = lgr_fl_mm)) +
  geom_boxplot() +
  theme_classic() +
  labs(x = "Ocean Age",
       y = "FL (mm)")
fl_by_age

fl_by_success = sf_df %>%
  select(lgr_fl_mm,
         success) %>%
  filter(!is.na(lgr_fl_mm)) %>%
  ggplot(aes(x = success,
             y = lgr_fl_mm)) +
  geom_boxplot() +
  theme_classic() +
  labs(x = "Successful Passage",
       y = "FL (mm)")
fl_by_success

#------------------
# passage by water level

# create fl bins
wl_bins = seq(0.3, 1.9, by = 0.2)
wl_bin_labels = c("0.3-0.5",
                  "0.5-0.7",
                  "0.7-0.9",
                  "0.9-1.1",
                  "1.1-1.3",
                  "1.3-1.5",
                  "1.5-1.7",
                  "1.7-1.9") 

pass_by_wl = sf_df %>%
  filter(SC1_SC2 == 1,
         !is.na(water_level_m)) %>%
  select(tag_code,
         success,
         travel_time,
         water_level_m) %>%
  mutate(wl_bin = cut(water_level_m,
                      breaks = wl_bins,
                      labels = wl_bin_labels)) %>%
  group_by(wl_bin) %>%
  summarize(n_tags = n(),
            n_tags_passed = sum(success),
            pct_passed = round((n_tags_passed / n_tags) * 100, 1),
            mean_travel_time_days = round(mean(travel_time, na.rm = T), 1),
            med_travel_time_days = round(median(travel_time, na.rm = T), 1))
pass_by_wl

# END SCRIPT
