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

# --------------------------
# load some data
load(here("data/derived_data/sf_cleartwater_passage_data.rda"))
load(here("data/derived_data/sf_clearwater_mean_daily_cfs.rda")) ; rm(sf_elk_daily_cfs, sf_elk_gage_info, sf_stites_gage_info)

# environmental data from the SC4 probe
sc4_env = read_csv(here("data/raw_data/SC4_enviro_f62a18b7-c2a7-4df8-899a-f0f26ba9a7dc.csv")) %>%
  rename(node = slug...1,
         metric = slug...2) %>%
  mutate(date = date(read_at)) %>%
  group_by(node, metric, date) %>%
  summarise(mean = round(mean(value), 2)) %>%
  ungroup()

# --------------------------
# prepare data for analysis
sf_df = sf_ch %>%
  # collapse some captures
  rowwise() %>%
  mutate(SC1_SC2 = ifelse(any(c(SC1, SC2) == 1), 1, 0),
         SC4_CRA = ifelse(any(c(SC4B0, SC4A0, CRAB0, CRAA0) == 1), 1, 0),
         SC1_SC2_last_det = max(c_across(SC1_last_det:SC2_last_det), na.rm = T),
         SC4_CRA_first_det = min(c_across(CRAA0_first_det:SC4A0_first_det), na.rm = T),
         SC1_SC2_to_SC3_tt = difftime(SC3_first_det, SC1_SC2_last_det, units = "days"),
         SC3_to_SC4_CRA_tt = difftime(SC4_CRA_first_det, SC3_last_det, units = "days")) %>%
  ungroup() %>%
  # convert strings containing "Inf" to NA
  mutate(across(SC1_SC2:SC3_to_SC4_CRA_tt, ~replace(., str_detect(., "Inf"), NA))) %>%
  #mutate(across(SC1_SC2:SC3_to_SC4_CRA_tt, ~if_else(str_detect(., "Inf"), NA, .))) %>%
  # create a column of capture histories
  mutate(SC3_tmp = SC3) %>%
  unite(ch, SC1:CRAA0, sep = "") %>%
  select(tag_code,
         spawn_year,
         ch,
         SC1_SC2,
         SC3 = SC3_tmp,
         SC4_CRA,
         SC1_SC2_last_det,
         SC3_first_det,
         SC3_last_det,
         SC4_CRA_first_det,
         SC1_SC2_to_SC3_tt,
         SC3_to_SC4_CRA_tt) %>%
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
  # join water level data from the SC4 probe, using date when fish was 
  # last detected at SC1 or SC2 (too few of detections at SC4)
  mutate(tmp = date(SC1_SC2_last_det)) %>%
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
  # did we detect fish passing SC3 or SC4?
  mutate(pass_SC3 = if_else(SC3 == 1 | SC4_CRA == 1, T, F),
         pass_SC4 = if_else(SC4_CRA == 1, T, F))

# -----------------------
# simple detection probs
sc1_sc2_p = sf_df %>%
  select(spawn_year,
         SC1_SC2,
         SC3,
         SC4_CRA) %>%
  mutate(pass_SC3 = if_else(SC3 == 1 | SC4_CRA == 1, T, F)) %>%
  filter(pass_SC3 == T) %>%
  group_by(spawn_year) %>%
  summarize(n_tags = n(),
            p = sum(SC1_SC2) / n())

sc3_p = sf_df %>%
  select(spawn_year, 
         SC3, 
         SC4_CRA) %>%
  filter(SC4_CRA == 1) %>%
  group_by(spawn_year) %>%
  summarize(n_tags = n(),
            p = sum(SC3) / n())

# -----------------------
# passage by spawn year
pass_by_sy = sf_df %>%
  group_by(spawn_year) %>%
  summarize(n_tags = n(),
            n_pass_SC3 = sum(pass_SC3),
            pct_pass_SC3 = round((n_pass_SC3 / n_tags)*100, 1),
            n_pass_SC4 = sum(pass_SC4),
            pct_SC3_pass_SC4 = round((n_pass_SC4 / n_pass_SC3)*100, 1),
            pct_pass_SC4 = round((n_pass_SC4 / n_tags)*100, 1)) %>%
  # add totals
  bind_rows(sf_df %>%
              summarize(spawn_year = -9999,
                        n_tags = n(),
                        n_pass_SC3 = sum(pass_SC3),
                        pct_pass_SC3 = round((n_pass_SC3 / n_tags)*100, 1),
                        n_pass_SC4 = sum(pass_SC4),
                        pct_SC3_pass_SC4 = round((n_pass_SC4 / n_pass_SC3)*100, 1),
                        pct_pass_SC4 = round((n_pass_SC4 / n_tags)*100, 1))) %>%
  mutate(spawn_year = replace(spawn_year, spawn_year == -9999, "Totals"))
pass_by_sy

# -----------------------
# passage by release group, year
pass_by_rg = sf_df %>%
  filter(!rel_site %in% c("CLEARC", "PRDLD1")) %>%
  mutate(rel_group = case_when(
    rel_site == "BONAFF" ~ "Unknown Release",
    rel_site == "CLWRSF" ~ "South Fork Clearwater River",
    rel_site == "COLR2"  ~ "Unknown Release",
    rel_site == "LGRLDR" ~ "Unknown Release",
    rel_site == "LGRRBR" ~ "Unknown Release",
    rel_site == "LGRRRR" ~ "Unknown Release",
    rel_site == "MEAD2C" ~ "Meadow Creek",
    rel_site == "NEWSOC" ~ "Newsome Creek",
    TRUE ~ rel_site)) %>%
  group_by(rel_group, spawn_year) %>%
  summarize(n_tags = n(),
            n_pass_SC3 = sum(pass_SC3),
            pct_pass_SC3 = round((n_pass_SC3 / n_tags)*100, 1),
            n_pass_SC4 = sum(pass_SC4),
            pct_SC3_pass_SC4 = round((n_pass_SC4 / n_pass_SC3)*100, 1),
            pct_pass_SC4 = round((n_pass_SC4 / n_tags)*100, 1)) %>%
  # add totals
  bind_rows(sf_df %>%
              filter(!rel_site %in% c("CLEARC", "PRDLD1")) %>%
              summarize(spawn_year = -9999,
                        n_tags = n(),
                        n_pass_SC3 = sum(pass_SC3),
                        pct_pass_SC3 = round((n_pass_SC3 / n_tags)*100, 1),
                        n_pass_SC4 = sum(pass_SC4),
                        pct_SC3_pass_SC4 = round((n_pass_SC4 / n_pass_SC3)*100, 1),
                        pct_pass_SC4 = round((n_pass_SC4 / n_tags)*100, 1))) %>%
  mutate(spawn_year = replace(spawn_year, spawn_year == -9999, "Totals"))

# bar plot
pass_by_rg_p = pass_by_rg %>%
  select(rel_group, spawn_year, pct_SC3_pass_SC4) %>%
  filter(spawn_year != "Totals") %>%
  ggplot(aes(x = rel_group,
             y = pct_SC3_pass_SC4,
             fill = spawn_year)) +
  geom_bar(position = "dodge",
           stat = "identity") +
  theme_classic() +
  labs(x = "Release Group",
       y = "% SC3 Passed SC4",
       fill = "Spawn Year") +
  ylim(0, 100) +
  theme(legend.position = "top",
        axis.text.x = element_text(angle = 45, vjust = 0.55))
pass_by_rg_p

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

pass_by_fl = sf_df %>%
  filter(!is.na(lgr_fl_mm)) %>%
  select(tag_code,
         pass_SC3,
         pass_SC4,
         lgr_fl_mm) %>%
  mutate(fl_bin = cut(lgr_fl_mm,
                      breaks = bins,
                      labels = bin_labels)) %>%
  group_by(fl_bin) %>%
  summarize(n_tags = n(),
            n_pass_SC3 = sum(pass_SC3),
            pct_pass_SC3 = round((n_pass_SC3 / n_tags)*100, 1),
            n_pass_SC4 = sum(pass_SC4),
            pct_SC3_pass_SC4 = round((n_pass_SC4 / n_pass_SC3)*100, 1),
            pct_pass_SC4 = round((n_pass_SC4 / n_tags)*100, 1))
pass_by_fl

#------------------
# passage by saltwater age
pass_by_sw_age = sf_df %>%
  filter(!is.na(bio_scale_final_age)) %>%
  select(tag_code,
         pass_SC3,
         pass_SC4,
         bio_scale_final_age) %>%
  mutate(sw_age = gsub(".*:","", bio_scale_final_age)) %>%
  mutate(sw_age = case_when(
    sw_age == "A" ~ NA,
    TRUE ~ sw_age)) %>%
  filter(!is.na(sw_age)) %>%
  group_by(sw_age) %>%
  summarize(n_tags = n(),
            n_pass_SC3 = sum(pass_SC3),
            pct_pass_SC3 = round((n_pass_SC3 / n_tags)*100, 1),
            n_pass_SC4 = sum(pass_SC4),
            pct_SC3_pass_SC4 = round((n_pass_SC4 / n_pass_SC3)*100, 1),
            pct_pass_SC4 = round((n_pass_SC4 / n_tags)*100, 1))
pass_by_sw_age

#------------------
# fork length by age
fl_by_age = sf_df %>%
  select(lgr_fl_mm,
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
       y = "FL (mm)")
fl_by_age

#------------------
# fork length by successful passage
fl_by_success = sf_df %>%
  select(tag_code,
         lgr_fl_mm,
         pass_SC3,
         pass_SC4) %>%
  filter(!is.na(lgr_fl_mm)) %>%
  ggplot(aes(x = pass_SC4,
             y = lgr_fl_mm)) +
  geom_boxplot() +
  theme_classic() +
  labs(x = "Successful Passage to SC4",
       y = "FL (mm)")
fl_by_success

#------------------
# passage by SC4 water level
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
  filter(!is.na(sc4_water_level_m)) %>%
  select(tag_code,
         pass_SC3,
         pass_SC4,
         sc4_water_level_m) %>%
  mutate(wl_bin = cut(sc4_water_level_m,
                      breaks = wl_bins,
                      labels = wl_bin_labels)) %>%
  group_by(wl_bin) %>%
  summarize(n_tags = n(),
            n_pass_SC3 = sum(pass_SC3),
            pct_pass_SC3 = round((n_pass_SC3 / n_tags)*100, 1),
            n_pass_SC4 = sum(pass_SC4),
            pct_SC3_pass_SC4 = round((n_pass_SC4 / n_pass_SC3)*100, 1),
            pct_pass_SC4 = round((n_pass_SC4 / n_tags)*100, 1))

pass_by_wl

#------------------
# passage by Stites discharge (cfs)
cfs_bins = seq(0, 10000, by = 250)
tmp = paste0(cfs_bins, "-", lead(cfs_bins))
cfs_bin_labels = tmp[1:(length(tmp) - 1)]

pass_by_cfs = sf_df %>%
  filter(!is.na(daily_mean_cfs)) %>%
  select(tag_code,
         pass_SC3,
         pass_SC4,
         daily_mean_cfs) %>%
  mutate(cfs_bin = cut(daily_mean_cfs,
                       breaks = cfs_bins,
                       labels = cfs_bin_labels)) %>%
  group_by(cfs_bin) %>%
  summarize(n_tags = n(),
            n_pass_SC3 = sum(pass_SC3),
            pct_pass_SC3 = round((n_pass_SC3 / n_tags)*100, 1),
            n_pass_SC4 = sum(pass_SC4),
            pct_SC3_pass_SC4 = round((n_pass_SC4 / n_pass_SC3)*100, 1),
            pct_pass_SC4 = round((n_pass_SC4 / n_tags)*100, 1))
pass_by_cfs

#------------------
# simple logistic regression
fl_rg_data = sf_df %>%
  filter(!rel_site %in% c("CLEARC", "PRDLD1"),
         !is.na(lgr_fl_mm)) %>%
  mutate(rel_group = case_when(
    rel_site == "BONAFF" ~ "Unknown Release",
    rel_site == "CLWRSF" ~ "South Fork Clearwater River",
    rel_site == "COLR2"  ~ "Unknown Release",
    rel_site == "LGRLDR" ~ "Unknown Release",
    rel_site == "LGRRBR" ~ "Unknown Release",
    rel_site == "LGRRRR" ~ "Unknown Release",
    rel_site == "MEAD2C" ~ "Meadow Creek",
    rel_site == "NEWSOC" ~ "Newsome Creek",
    TRUE ~ rel_site)) %>%
  mutate(pass_SC4 = if_else(pass_SC4 == T, 1, 0)) %>%
  select(tag_code,
         lgr_fl_mm,
         pass_SC4,
         rel_group)

# fit the logistic regression model
model = glm(pass_SC4 ~ lgr_fl_mm + rel_group, data = fl_rg_data, family = binomial())

# create a sequence of x values for prediction
x_pred = seq(min(fl_rg_data$lgr_fl_mm), max(fl_rg_data$lgr_fl_mm), length.out = 1000)

# create a data frame for prediction including groups
pred_data = expand.grid(x = x_pred, group = unique(fl_rg_data$rel_group)) %>%
  rename(lgr_fl_mm = x,
         rel_group = group)

# predict the probabilities using the logistic regression model
pred_probs = predict(model, newdata = pred_data, type = "response")

# combine the predicted data with the original data
plot_data = cbind(pred_data, pred_prob = pred_probs)

# Plot the logistic regression curve and data points
lr_p = ggplot() +
  geom_point(aes(x = lgr_fl_mm,
                 y = pass_SC4,
                 color = rel_group),
             data = fl_rg_data) +
  geom_line(aes(x = lgr_fl_mm,
                y = pred_prob,
                color = rel_group),
            data = plot_data,
            linewidth = 0.5) +
  labs(x = "FL (mm)",
       y = "p(Successful Pass)",
       color = "Release Group") +
  theme_minimal()
lr_p

# END SCRIPT
