# Purpose: Summarize steelhead run timing in South Fork Clearwater River
# 
# Authors: Mike Ackerman
# 
# Created: May 20, 2024
#   Last Modified: May 21, 2025

# clear environment
rm(list = ls())

# load necessary packages
library(tidyverse)
library(here)
library(PITcleanr)
library(janitor)
library(data.table)

#---------------------
# Compile Relevant Data

# load dart_obs_list and convert to a data frame (rbindlist avoids issues with differing data types)
load(here("data/derived_data/dart_observations/sy12-25_dart_obs.rda"))
dart_obs_df = data.table::rbindlist(dart_obs_list) %>%
  # only need steelhead for this analysis
  filter(species == "Steelhead")

# load compressed and filtered PITcleanr observations
load(here("data/derived_data/cths/sy12-25_compressed_filtered_obs.rda"))
comp_df = bind_rows(comp_list) %>%
  filter(species == "Steelhead")

#---------------------
# Begin Run-Timing Summaries

# get rears from dart_obs_df
rears = dart_obs_df %>%
  select(tag_id,
         spawn_year,
         # what is the source of t_rear_type?
         t_rear_type) %>%
  distinct()

# compile into single data frame and add rear
sf_obs = comp_df %>%
  # add rear from dart_obs_df
  left_join(dart_obs_df %>%
              select(tag_id,
                     spawn_year,
                     t_rear_type) %>%
              distinct(),
            by = c("tag_code" = "tag_id",
                   "spawn_year")) %>%
  select(spawn_year,
         tag_code,
         node,
         min_det,
         n_dets,
         rear = t_rear_type) %>%
  mutate(date = format(min_det, "%m-%d")) %>%
  # only 3 records in spawn_year 2012
  filter(spawn_year != 2012) %>%
  # only 11 records for CRA across 3 year
  filter(node != "CRA")
  
# number of rear by year
sf_obs %>% 
  tabyl(spawn_year, rear)

# number per node by year
sf_obs %>%
  tabyl(spawn_year, node)

# summarize tags per day by spawn_year, node, and rear
tags_per_day = sf_obs %>%
  # remove records with unknown rear
  filter(rear != "U") %>%
  group_by(spawn_year,
           node,
           rear,
           date) %>%
  summarize(n_tags = n_distinct(tag_code),
            .groups = "drop")

# create sequence of dates of interest
run_dates = seq.Date(from = as.Date("2025-01-01"),
                     to = as.Date("2025-05-31"),
                     by = "day") %>%
  # convert to %m-%d format
  format("%m-%d")

# create blank data frame with all unique combinations 
run_df = expand.grid(
  spawn_year = unique(tags_per_day$spawn_year),
  node = unique(tags_per_day$node),
  rear = unique(tags_per_day$rear),
  date = run_dates) %>%
  left_join(tags_per_day,
            by = c("spawn_year", "node", "rear", "date")) %>%
  # replace NAs with 0
  mutate(n_tags = replace_na(n_tags, 0))

# summarize totals for each date
totals = run_df %>%
  group_by(node, rear, date) %>%
  summarise(n_tags = sum(n_tags),
            .groups = "drop") %>%
  mutate(spawn_year = "Total")

# combine run_df and totals
run_df2 = rbind(run_df, totals) %>%
  arrange(spawn_year,
          node,
          rear,
          date) %>%
  # remove any groups that have no data
  group_by(spawn_year,
           node,
           rear) %>%
  filter(any(n_tags != 0)) %>%
  # add cumulative tags by spawn_year, rear, and node
  mutate(c_tags = cumsum(n_tags),
         max_tags = max(c_tags),
         p_tags = c_tags / max_tags) %>%
  ungroup()

# run timing date quantiles
sf_run_quants = run_df2 %>%
  group_by(spawn_year, 
           node, 
           rear) %>%
  mutate(q10 = date[which.min(abs(p_tags - 0.10))],
         q25 = date[which.min(abs(p_tags - 0.25))],
         q50 = date[which.min(abs(p_tags - 0.50))],
         q75 = date[which.min(abs(p_tags - 0.75))],
         q90 = date[which.min(abs(p_tags - 0.90))]) %>%
  distinct(spawn_year, node, rear, .keep_all = TRUE) %>%
  select(-date,
         -n_tags,
         -c_tags,
         -p_tags,
         n_tags = max_tags) %>%
  arrange(spawn_year, node, rear)
  
# write run timing quantiles table
write_csv(sf_run_quants,
          here("output/run_timing/sf_clearwater_sthd_run_timing_quantiles.csv"))

sites = unique(run_df$node)
plot_list = list()
for(s in sites) {
  
  site_df = run_df2 %>%
    filter(node == s) %>%
    mutate(date = as.Date(date, format = "%m-%d")) %>%
    # remove combos with small sample sizes
    group_by(spawn_year, node, rear) %>%
    filter(any(max_tags >= 20)) %>%
    ungroup() %>%
    # recode rears
    mutate(rear = recode(rear,
                         "H" = "Hatchery",
                         "W" = "Natural"))
   
    n_yrs = length(unique(site_df$spawn_year))
  
  site_p = site_df %>%
    ggplot(aes(x = date,
               y = p_tags,
               group = spawn_year)) +
    geom_line(aes(color = spawn_year,
                  linetype = spawn_year,
                  alpha = spawn_year == "Total")) +
    scale_color_manual(values = c(rep("gray20", n_yrs - 1), "black")) +
    scale_linetype_manual(values = c(rep("dotted", n_yrs - 1), "solid")) +
    scale_alpha_manual(values = c(rep(0.5, n_yrs - 1), 1)) +
    scale_x_date(date_breaks = "7 days", date_labels = "%b %d") +
    facet_wrap(~ rear, ncol = 1) +
    labs(x = NULL, y = "Proportion of PIT Tags", title = paste0("Cumulative Proportion of PIT Tags Passing ", s)) +
    theme_bw() +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = -45, hjust = 0.5, vjust = 0.5))
  
  plot_list[[s]] = site_p
 
}
all_sites_p = gridExtra::marrangeGrob(plot_list, nrow = 1, ncol = 1)
ggsave(paste0(here("output/run_timing/sf_clearwater_sthd_run_time.pdf")),
       all_sites_p,
       width = 8.5,
       height = 14,
       units = "in")

### END SCRIPT