rm(list=ls())

## data collation
packages <- c("data.table", "tidyverse", "sf", "sp", "spdep",'ggplot2','extRemes',
              "geofacet", "ggpubr", "ggthemes", 'readxl', 'lubridate','dplyr','readr',"extRemes","ggcorrplot")
lapply(packages, library, character.only = TRUE)

#######################Mitigation potential in the future##################################


models_to_keep <- c('ACCESS-ESM1-5', 'BCC-CSM2-MR', 'CMCC-CM2-SR5', 'CMCC-ESM2', 
                    'CNRM-CM6-1', 'EC-Earth3', 'FGOALS-g3', 'GFDL-CM4', 
                    'MIROC6', 'MPI-ESM1-2-HR')


# 1. Filter observation and historical simulation data (1980–2020)
precbase0 <- precbase %>%
  filter(year(date) >= 1981 & year(date) <= 2020)

his0 <- his %>%
  filter(year(date) >= 1981 & year(date) <= 2020 & model %in% models_to_keep)

# 2. Apply quantile mapping correction to future precipitation projections
files <- list.files(path = "/01 prec/", 
                    pattern = "CityDailyPrecipitation_AllModels_AllScenarios\\d+.csv", 
                    full.names = TRUE)

for (file in files) {
  precx <- read_csv(file) %>%
    filter(model %in% models_to_keep)
  
  prec_future_corrected <- precx %>%
    group_by(citycode, model, scenario) %>%
    do({
      obs_data <- precbase0 %>% filter(citycode == .$citycode[1]) %>% pull(prec)
      sim_data <- his0 %>% filter(citycode == .$citycode[1] & model == .$model[1]) %>% pull(prec)
      
      fit_qm <- fitQmapQUANT(obs = obs_data, mod = sim_data, qstep = 0.01)
      future_sim <- .$prec
      future_calibrated <- doQmapQUANT(future_sim, fit_qm)
      
      data.frame(citycode = .$citycode[1], date = .$date, prec = future_calibrated)
    }) %>%
    ungroup()
  
  prec_future_corrected <- merge(prec_future_corrected, return_levels3, by = "citycode", all.x = TRUE)
  
  test <- prec_future_corrected %>%
    filter(prec >= 50 | prec > return_level) %>%
    distinct(citycode, model, scenario, date, .keep_all = TRUE)
  
  file_name <- tools::file_path_sans_ext(basename(file))
  output_path <- paste0("/01 prec/correct/", file_name, "_corrected.csv")
  write_csv(test, output_path)
}

# 3. Merge all corrected projection files
directory <- "/01 prec/correct/"
corrected_files <- list.files(directory, full.names = TRUE)
corrected_data <- do.call(rbind, lapply(corrected_files, read.csv))

corrected_data <- corrected_data %>%
  rename(prec = prec_future_calibrated) %>%
  select(-return_level, -return_period)

# 4. Assign return periods based on return level thresholds
combined <- corrected_data %>%
  left_join(return_levels, by = "citycode") %>%
  arrange(citycode, return_period) %>%
  group_by(citycode, date, scenario, model) %>%
  summarize(
    prec = first(prec),
    return_period = case_when(
      prec < min(return_level, na.rm = TRUE) ~ NA_real_,
      prec >= max(return_level, na.rm = TRUE) ~ max(return_period, na.rm = TRUE),
      TRUE ~ {
        idx <- which(return_level[-length(return_level)] <= prec & prec < return_level[-1])
        if (length(idx) > 0) return_period[idx] else NA_real_
      }
    ),
    .groups = "drop"
  )

# 5. Detect consecutive precipitation events by city
combined$date <- as.Date(combined$date)
combined <- combined %>%
  arrange(scenario, model, citycode, date) %>%
  group_by(citycode) %>%
  mutate(event1 = cumsum(c(1, diff(date)) > 1)) %>%
  ungroup()

# 6. Cluster events across adjacent and nearby cities
scenarios_models <- unique(combined %>% select(scenario, model))
out_dir <- "/01 prec/size/"

for (row in seq_len(nrow(scenarios_models))) {
  current_scenario <- scenarios_models$scenario[row]
  current_model <- scenarios_models$model[row]
  
  subset_data <- combined %>% filter(scenario == current_scenario, model == current_model)
  years <- unique(year(subset_data$date))
  
  for (year_i in years) {
    year_data <- subset_data %>% filter(year(date) == year_i)
    year_data$event2 <- NA_integer_
    event_counter <- 1
    
    for (i in seq_len(nrow(year_data))) {
      if (is.na(year_data$event2[i])) {
        current_event1 <- year_data$event1[i]
        current_city <- year_data$citycode[i]
        
        queue <- year_data %>% filter(citycode == current_city, event1 == current_event1) %>%
          select(date, citycode) %>% distinct()
        checked <- data.frame(date = as.Date(character()), citycode = numeric())
        
        while (nrow(queue) > 0) {
          current <- queue[1, ]; queue <- queue[-1, ]
          if (nrow(checked %>% filter(date == current$date & citycode == current$citycode)) > 0) next
          
          checked <- rbind(checked, current)
          current_event_rows <- year_data %>% filter(citycode == current$citycode, event1 == current_event1)
          year_data$event2[which(year_data$citycode %in% current_event_rows$citycode & year_data$event1 %in% current_event_rows$event1)] <- event_counter
          
          adjacent <- adjacency_map %>% filter(citycode == current$citycode) %>% pull(adjacent_citycode)
          nearby <- dis %>% filter(ocitycode == current$citycode, distance <= 200) %>% pull(dcitycode)
          related_cities <- unique(c(adjacent, nearby))
          
          for (related_city in related_cities) {
            related_event1_rows <- year_data %>%
              filter(citycode == related_city, date == current$date, is.na(event2)) %>%
              select(date, citycode) %>% distinct()
            
            if (nrow(related_event1_rows) > 0) {
              queue <- unique(rbind(queue, related_event1_rows))
            }
          }
        }
        
        event_counter <- event_counter + 1
      }
    }
    
    write.csv(year_data, paste0(out_dir, current_scenario, "_", current_model, "_", year_i, ".csv"), row.names = FALSE)
  }
}

# 7. Combine event files and compute affected city count
files_final <- list.files("/01 prec/size/", full.names = TRUE)
future <- do.call(rbind, lapply(files_final, read.csv))

future2 <- future %>%
  group_by(scenario, model, year, event2) %>%
  summarize(affectedsize = n_distinct(citycode), .groups = "drop")

future <- merge(future, future2, by = c("scenario", "model", "year", "event2"))

    
##8. calculate affected cities in each event 
future2<-future %>%
  group_by(scenario,model,year,event2) %>%
  summarize(affectedsize = n_distinct(citycode), .groups = "drop")

future<-merge(future,future2,by=c('scenario',"model","year","event2"))

##9. percentage of city-days in multi-city events. fig s4
future_filtered <- future %>%
  filter(affectedsize > 2)

proportion_df <- future %>%
  group_by(scenario, year, model) %>%
  summarize(
    prop = mean(affectedsize > 2), 
    .groups = "drop"
  )

ratio_summary_result <- proportion_df %>%
  group_by(scenario, year) %>%
  summarize(
    row_mean = mean(prop, na.rm = TRUE),
    row_25th = quantile(prop, 0.25, na.rm = TRUE),
    row_75th = quantile(prop, 0.75, na.rm = TRUE),
    .groups = "drop"
  )

p <- ggplot(ratio_summary_result, aes(x = year, color = scenario, fill = scenario)) +
  geom_ribbon(aes(ymin = row_25th, ymax = row_75th), alpha = 0.2, color = NA) +
  geom_line(aes(y = row_mean), size = 0.5) +
  theme_bw() +
  theme(
    text = element_text(family = "Arial", size = 7),
    legend.position = "none",
    axis.title = element_text(size = 7),
    axis.text = element_text(size = 6),
    plot.title = element_text(hjust = 0.5),
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", size = 0.5)
  ) +
  labs(
    x = "Year",
    y = "Proportion of Models with Affected Size > 2"
  ) +
  scale_x_continuous(limits = c(2025, 2100), breaks = seq(2025, 2100, by = 25)) +
  coord_cartesian(expand = FALSE)



# Future Rainfall Event Analysis under Climate Scenarios--fig 4

## 1. Event count by city across models and years
future_event2_city <- future %>%
  group_by(scenario, model, year, citycode) %>%
  summarize(
    event2_count = n_distinct(event2),
    row_count = n(),
    .groups = "drop"
  )

## 2. Take model average
future_event2_city1 <- future_event2_city %>%
  group_by(scenario, year, citycode) %>%
  summarize(
    event2_count = mean(event2_count),
    row_count = mean(row_count),
    .groups = "drop"
  )

## 3. Sum across years
future_event2_city2 <- future_event2_city1 %>%
  group_by(scenario, citycode) %>%
  summarize(
    event2_count = sum(event2_count),
    row_count = sum(row_count),
    .groups = "drop"
  )

## 4. Summary statistics by scenario
summary(future_event2_city2$event2_count[future_event2_city2$scenario=='ssp245'])
summary(future_event2_city2$row_count[future_event2_city2$scenario=='ssp245'])

summary(future_event2_city2$event2_count[future_event2_city2$scenario=='ssp585'])
summary(future_event2_city2$row_count[future_event2_city2$scenario=='ssp585'])

## 5. Calculate yearly statistics for event2 count and day count
future_event2_counts <- future %>%
  group_by(scenario, model, year) %>%
  summarize(
    event2_count = n_distinct(event2),
    row_count = n(),
    .groups = "drop"
  )

## 6. Compute future ratios compared to baseline
future_summary_result_with_ratios <- future_event2_counts %>%
  mutate(
    row_ratio = (row_count - 450) / 450,
    event_ratio = (event2_count - 102) / 102
  )

future_summary_result <- future_summary_result_with_ratios %>%
  group_by(scenario, year) %>%
  summarise(
    row_mean = mean(row_ratio, na.rm = TRUE),
    row_25th = quantile(row_ratio, 0.25, na.rm = TRUE),
    row_75th = quantile(row_ratio, 0.75, na.rm = TRUE),
    event_mean = mean(event_ratio, na.rm = TRUE),
    event_25th = quantile(event_ratio, 0.25, na.rm = TRUE),
    event_75th = quantile(event_ratio, 0.75, na.rm = TRUE),
    .groups = "drop"
  )

## 7. Visualization - Future Daily City Count
p <- ggplot(future_summary_result, aes(x = year, color = scenario, fill = scenario)) +
  geom_ribbon(aes(ymin = row_25th, ymax = row_75th), alpha = 0.2, color = NA) +
  geom_line(aes(y = row_mean), size = 0.5) +
  theme_bw() +
  theme(
    text = element_text(family = "Arial", size = 7),
    legend.position = "none",
    axis.title = element_text(size = 7),
    axis.text = element_text(size = 6),
    plot.title = element_text(hjust = 0.5),
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", size = 0.5)
  ) +
  labs(x = "Year", y = "Change Ratio of Extreme Rainfall Day (%)") +
  scale_x_continuous(limits = c(2025, 2100), breaks = seq(2025, 2100, by = 25)) +
  coord_cartesian(expand = FALSE)

## 8. Visualization - Future Event Count
p <- ggplot(future_summary_result, aes(x = year, color = scenario, fill = scenario)) +
  geom_ribbon(aes(ymin = event_25th, ymax = event_75th), alpha = 0.2, color = NA) +
  geom_line(aes(y = event_mean), size = 0.5) +
  theme_bw() +
  theme(
    text = element_text(family = "Arial", size = 7),
    legend.position = "none",
    axis.title = element_text(size = 7),
    axis.text = element_text(size = 6),
    plot.title = element_text(hjust = 0.5),
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", size = 0.5)
  ) +
  labs(x = "Year", y = "Change Ratio of Extreme Rainfall Event (%)") +
  scale_x_continuous(limits = c(2025, 2100), breaks = seq(2025, 2100, by = 25)) +
  coord_cartesian(expand = FALSE)

## 9. Export summary result
write_csv(future_summary_result, "~/Desktop/04 heat&warning/08 result/fig/fig4 spfuturenumberchange.csv")

# (Note: Multi-city event analysis and plots will be added separately)





##Analysis of Future benefit from intercity connectivity -- fig 5  

# 1. Read input data
future <- read_csv("future.csv")                # future event records
events0 <- read_csv("events0.csv")              # baseline event records in 2023
windspeed <- read_csv("windspeed.csv")           # average windspeed at different precipitation level (return_period) from 19821-2020
max_cpolicy0 <- read_csv("max_cpolicy0.csv")    # maximum city-level typhoon warning from 2021-2023
socioeco <- read_csv("socioeco.csv")            # socioeconomic data
network_iso <- read_csv("network_iso.csv")      # isolated network features
network_weak <- read_csv("network_weak.csv")    # weak-control network features

# 2. Select peak precipitation per event
futurepro <- future %>%
  left_join(windspeed, by = c("citycode", "return_period" = "returnperiod")) %>%
  left_join(max_cpolicy0, by = "citycode") %>%
  group_by(citycode, scenario, model, year, event2) %>%
  slice_max(prec, n = 1) %>%
  ungroup()

futurepro$total_prec<-12.2512+1.1032*futurepro$prec

# 3. Compute total_prec and policy tier
futurepro <- futurepro %>%
  mutate(
    total_prec = 12.2512 + 1.1032 * prec,
    max_Rpolicy = case_when(
      total_prec <  50 ~ 0,
      prec <         50 ~ 1,
      total_prec < 100 ~ 2,
      TRUE             ~ 3
    ),
    return_period = coalesce(return_period, 1)
  )

# 4. Merge socioeconomic and network data
socioeco <- socioeco %>%
  select(city_code, mean_elevation, proportion_10, proportion_60,
         damdes, gdppercapita, pop)
futurepro <- futurepro %>%
  inner_join(socioeco, by = c("citycode" = "city_code")) %>%
  inner_join(network_iso, by = "citycode") %>%
  rename_with(~ sub("\.x$", "", .x, perl = TRUE), ends_with(".x"))

futurepro$tergdp<-futurepro$tergdp/100
futurepro$total_prec<-futurepro$total_prec/100
futurepro$mean_elevation<-futurepro$mean_elevation/1000

# 5. Fit baseline regression model on historical events
events1 <- events0 %>%
  group_by(citycode, event2) %>%
  slice_max(prec, n = 1) %>%
  ungroup() %>%
  mutate(total_prec = 12.2512 + 1.1032 * prec) %>%
  inner_join(socioeco, by = c("citycode" = "city_code")) %>%
  inner_join(network_iso, by = "citycode") %>%
  mutate(return_period = coalesce(return_period, 1))

events1$tergdp.x<-events1$tergdp.x/100
events1$total_prec<-events1$total_prec/100
events1$mean_elevation.x<-events1$mean_elevation.x/1000

events1$max_Cpolicy<-events1$Cpolicy
events1$max_Rpolicy<-events1$Rpolicy

events1iso <- events1 %>%
  left_join(
    iso %>% select(citycode, interaction_unevenness(mad), economic-tier_mobility_concentration, rich_neighbor_preference), 
    by = "citycode"
  ) %>%
  mutate(
    interaction_unevenness(mad) = coalesce(interaction_unevenness(mad).y, interaction_unevenness(mad).x),
    economic-tier_mobility_concentration = coalesce(economic-tier_mobility_concentration.y, economic-tier_mobility_concentration.x),
    rich_neighbor_preference = coalesce(rich_neighbor_preference.y, rich_neighbor_preference.x)
  ) %>%
  select(-interaction_unevenness(mad).x, -interaction_unevenness(mad).y, 
         -economic-tier_mobility_concentration.x, -economic-tier_mobility_concentration.y, 
         -rich_neighbor_preference.x, -rich_neighbor_preference.y)


overall$tergdp<-overall$tergdp/100
overall$total_prec<-overall$total_prec/100
overall$mean_elevation<-overall$mean_elevation/1000
overall$unique_city_count
hlm_model_test1 <- lm(
  total_intrac ~ total_prec + max_wind_diff + mean_elevation +proportion_10+proportion_60+
    proportion_90 + damdes + max_Rpolicy + max_Cpolicy  +mobility_per_capita+interaction_unevenness(mad)+economic-tier_mobility_concentration+rich_neighbor_preference,
  data = overall)

##  2023 observed network
events1$predicted_total_intrac2 <- predict(hlm_model_test1, newdata = events1)

##  isolated network
events1iso$predicted_total_intrac0 <- predict(hlm_model_test1, newdata = events1iso)


events1iso$netpredit<-events1iso$predicted_total_intrac2-events1iso$predicted_total_intrac0
events1iso$finalgdp<-events1iso$predicted_total_intrac2*events1iso$pop*events1iso$gdppercapita*(100-events1iso$primgdp)/(10000000000*365)
events1iso$savegdp<-events1iso$netpredit*events1iso$pop*events1iso$gdppercapita*(100-events1iso$primgdp)/(10000000000*365)

events1iso$overgdp<-events1iso$pop*events1iso$gdppercapita/(100000000)

events1isosta<-events1iso%>%
  group_by(citycode)%>%
  summarise( total_savegdp = sum(savegdp, na.rm = TRUE),
             overgdp = first(overgdp),
             savegdpper = 100 * total_savegdp / overgdp)

events1isosta$total_savegdpus<-events1isosta$total_savegdp/70.467 ##change from 100 million RMB to  billion us


# 6. future predicted benefits

##  2023 observed network
futurepro1$predicted_total_intrac2 <- predict(hlm_model_test1, newdata = futurepro1)

##  isolated network
futureproiso <- futurepro1 %>%
  left_join(
    iso %>% select(citycode, interaction_unevenness(mad), economic-tier_mobility_concentration, rich_neighbor_preference), 
    by = "citycode"
  ) %>%
  mutate(
    interaction_unevenness(mad) = coalesce(interaction_unevenness(mad).y, interaction_unevenness(mad).x),
    economic-tier_mobility_concentration = coalesce(economic-tier_mobility_concentration.y, economic-tier_mobility_concentration.x),
    rich_neighbor_preference = coalesce(rich_neighbor_preference.y, rich_neighbor_preference.x)
  ) %>%
  select(-interaction_unevenness(mad).x, -interaction_unevenness(mad).y, 
         -economic-tier_mobility_concentration.x, -economic-tier_mobility_concentration.y, 
         -rich_neighbor_preference.x, -rich_neighbor_preference.y)

futureproiso$predicted_total_intrac0 <- predict(hlm_model_test1, newdata = futureproiso)

futureproiso$netpredit<-futureproiso$predicted_total_intrac2-futureproiso$predicted_total_intrac0

futureproiso$finalgdp<-futureproiso$predicted_total_intrac2*futureproiso$pop*futureproiso$gdppercapita*(100-futureproiso$primgdp)/(10000000000*365)

futureproiso$isogdp<-futureproiso$predicted_total_intrac0*futureproiso$pop*futureproiso$gdppercapita*(100-futureproiso$primgdp)/(10000000000*365)

futureproiso$savegdp<-futureproiso$netpredit*futureproiso$pop*futureproiso$gdppercapita*(100-futureproiso$primgdp)/(10000000000*365)

futureproiso_sta1<-futureproiso%>%
  group_by(scenario,model,year,citycode)%>%
  summarise(netpredit=sum(netpredit),
            predicted_total_intrac0=sum(predicted_total_intrac0),
            predicted_total_intrac2=sum(predicted_total_intrac2),
            affectgdp=sum(finalgdp),
            isogdp=sum(isogdp),
            savegdp=sum(savegdp))


futureproiso_sta2 <- futureproiso_sta1 %>%
  group_by(scenario, year, citycode) %>%
  summarise(
    netpredit_mean = mean(netpredit, na.rm = TRUE),
    netpredit_25th = quantile(netpredit, 0.25, na.rm = TRUE),
    netpredit_75th = quantile(netpredit, 0.75, na.rm = TRUE),
    
    affectedgdp_mean = mean(affectgdp, na.rm = TRUE),
    affectedgdp_25th = quantile(affectgdp, 0.25, na.rm = TRUE),
    affectedgdp_75th = quantile(affectgdp, 0.75, na.rm = TRUE),
    
    isogdp_mean = mean(isogdp, na.rm = TRUE),
    isogdp_25th = quantile(isogdp, 0.25, na.rm = TRUE),
    isogdp_75th = quantile(isogdp, 0.75, na.rm = TRUE),
    
    savegdp_mean = mean(savegdp, na.rm = TRUE),
    savegdp_25th = quantile(savegdp, 0.25, na.rm = TRUE),
    savegdp_75th = quantile(savegdp, 0.75, na.rm = TRUE),
    .groups = "drop"
  )


# aggregate across cities
futureproiso_result <- futureproiso_sta2 %>%
  group_by(scenario, year) %>%
  summarise(
    netpredit = sum(netpredit_mean, na.rm = TRUE),
    netpredit_25th = sum(netpredit_25th, na.rm = TRUE),
    netpredit_75th = sum(netpredit_75th, na.rm = TRUE),
    
    affectedgdp = sum(affectedgdp_mean, na.rm = TRUE),
    affectedgdp_25th = sum(affectedgdp_25th, na.rm = TRUE),
    affectedgdp_75th = sum(affectedgdp_75th, na.rm = TRUE),
    
    isogdp = sum(isogdp_mean, na.rm = TRUE),
    isogdp_25th = sum(isogdp_25th, na.rm = TRUE),
    isogdp_75th = sum(isogdp_75th, na.rm = TRUE),
    
    savegdp = sum(savegdp_mean, na.rm = TRUE),
    savegdp_25th = sum(savegdp_25th, na.rm = TRUE),
    savegdp_75th = sum(savegdp_75th, na.rm = TRUE),
    .groups = "drop"
  )

# 7. fig 5

futureproiso_result1 <- futureproiso_result %>%
  group_by(scenario, year) %>%
  mutate(
    netpredit_ratio =100* (netpredit - sum(events1iso$netpredit))/sum(events1iso$netpredit),
    netpredit_25th_ratio = 100*(netpredit_25th -sum(events1iso$netpredit))/sum(events1iso$netpredit),
    netpredit_75th_ratio = 100*(netpredit_75th -sum(events1iso$netpredit))/sum(events1iso$netpredit),
    
    savegdp_ratio1 = 100*(savegdp -sum(events1iso$savegdp))/sum(events1iso$savegdp),
    savegdp_25th_ratio1 =100* (savegdp_25th-sum(events1iso$savegdp)) / sum(events1iso$savegdp),
    savegdp_75th_ratio1 = 100*(savegdp_75th-sum(events1iso$savegdp)) / sum(events1iso$savegdp),
  )


p <- ggplot(futureproiso_result1, aes(x = year, color = scenario, fill = scenario)) +
  geom_ribbon(aes(ymin = netpredit_25th_ratio, ymax = netpredit_75th_ratio), 
              alpha = 0.2, color = NA) +

  geom_line(aes(y = netpredit_ratio), size = 0.4) +
  theme_bw(base_family = "sans") +
  theme(
    text = element_text(family = "sans", size = 7),
    legend.position = "none",
    axis.title = element_text(size = 7),
    axis.text = element_text(size = 6),
    plot.title = element_text(hjust = 0.5),
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", size = 0.5)
  ) +
  labs(
    x = "Year",
    y = "Resilience improvement(%)"
  ) +
  scale_x_continuous(limits = c(2025, 2100), breaks = seq(2025, 2100, by = 25)) +
  scale_y_continuous(limits = c(0, 350), breaks = seq(0, 340, by = 50)) +
  coord_cartesian(expand = FALSE)


p1 <- ggplot(futureproiso_result1, aes(x = year, color = scenario, fill = scenario)) +
  geom_ribbon(aes(ymin = savegdp_25th_ratio1, ymax = savegdp_75th_ratio1), 
              alpha = 0.2, color = NA) +
  geom_line(aes(y = savegdp_ratio1), size = 0.4) +
  theme_bw(base_family = "Arial") +
  theme(
    text = element_text(family = "Arial", size = 7),
    legend.position = "none",
    axis.title = element_text(size = 7),
    axis.text = element_text(size = 6),
    plot.title = element_text(hjust = 0.5),
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", size = 0.5)
  ) +
  labs(
    x = "Year",
    y = "Avoided economic loss (%)"
  ) +
  scale_x_continuous(limits = c(2025, 2100), breaks = seq(2025, 2100, by = 25)) +
  scale_y_continuous(limits = c(-25, 140), breaks = seq(0, 125, by = 25)) +
  coord_cartesian(expand = FALSE)


p2 <- ggplot(futureproiso_result1, aes(x = year, color = scenario, fill = scenario)) +
  geom_ribbon(aes(ymin = savegdp_25th/70.467, ymax = savegdp_75th/70.467), 
              alpha = 0.2, color = NA) +
  geom_line(aes(y = savegdp/70.467), size = 0.4) +
  theme_bw(base_family = "Arial") +
  theme(
    text = element_text(family = "Arial", size = 7),
    legend.position = "none",
    axis.title = element_text(size = 7),
    axis.text = element_text(size = 6),
    plot.title = element_text(hjust = 0.5),
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", size = 0.5)
  ) +
  labs(
    x = "Year",
    y = "Avoided economic loss (billion USD)"
  ) +
  scale_x_continuous(limits = c(2025, 2100), breaks = seq(2025, 2100, by = 25)) +
  scale_y_continuous(limits = c(5, 20), breaks = seq(5, 20, by = 5)) +
  coord_cartesian(expand = FALSE)

# summary(futureproiso_result1$savegdp/70.467)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 8.678   9.823  10.594  10.773  11.201  15.014 
# 
# summary(futureproiso_result1$savegdp/70.467[futureproiso_result1$scenario=='ssp245'])
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
#   8.678   8.847   9.015   9.015   9.184   9.352     150 
  
data2<-futureproiso_result1%>%
  filter(scenario=='ssp585')
summary(data2$savegdp/70.467)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 8.965  10.157  11.109  11.277  12.531  15.014 


### compared to connectivity in weak-control network are exactly the same with above, with intercity connectivity features changed to those in weak-control network


##########

##Future benefit from intercity connectivity of different city groups-- extended data fig 3 


### Stratify the results by city population size

futureproiso_resultcg <- futureproiso_resultcg %>%
  mutate(
    pop_category = case_when(
      pop >= 10000000 ~ "mega",
      pop > 5000000&pop <=10000000 ~ "large",
      pop > 1000000&pop <=5000000 ~ "median",
      pop <= 1000000 ~ "small",
    )
  )

futureproweak_resultcg <- futureproweak_resultcg %>%
  mutate(
    pop_category = case_when(
      pop >= 10000000 ~ "mega",
      pop > 5000000&pop <=10000000 ~ "large",
      pop > 1000000&pop <=5000000 ~ "median",
      pop <= 1000000 ~ "small",
    )
  )

futureproiso_resultcg0 <- futureproiso_resultcg %>%
  group_by(scenario, year,pop_category) %>%
  summarise(
    affectedgdp = sum(affectedgdp_mean, na.rm = TRUE),
    affectedgdp_25th = sum(affectedgdp_25th, na.rm = TRUE),
    affectedgdp_75th = sum(affectedgdp_75th, na.rm = TRUE),
    
    isogdp = sum(isogdp_mean, na.rm = TRUE),
    isogdp_25th = sum(isogdp_25th, na.rm = TRUE),
    isogdp_75th = sum(isogdp_75th, na.rm = TRUE),
    
    savegdp = sum(savegdp_mean, na.rm = TRUE),
    savegdp_25th = sum(savegdp_25th, na.rm = TRUE),
    savegdp_75th = sum(savegdp_75th, na.rm = TRUE),
    .groups = "drop"
  )

futureproiso_resultcg1 <- futureproiso_resultcg0 %>%
  group_by(scenario, year,pop_category) %>%
  mutate(
    savegdp_ratio1 = 100*(savegdp)/gdp,
    affectedgdp_ratio1 = 100*(affectedgdp )/gdp,
  )

futureproiso_resultcg10 <- futureproiso_resultcg1 %>%
  group_by(scenario, pop_category) %>%
  summarise(
    savegdp_ratio1     = mean(savegdp_ratio1,     na.rm = TRUE),
    affectedgdp_ratio1 = mean(affectedgdp_ratio1, na.rm = TRUE)
  ) %>%
  ungroup()



