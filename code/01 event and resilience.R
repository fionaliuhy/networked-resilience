rm(list=ls())

# -------------------------------------------------------------------------
# 0. Load Packages and Settings
# -------------------------------------------------------------------------
packages <- c("data.table", "tidyverse", "sf", "sp", "spdep",
              "dlnm", "tsModel", "hydroGOF","RColorBrewer", 
              "geofacet", "ggpubr", "ggthemes","extRemes", 
              "gghalves", "extrafont", "relaimpo", "effectsize")

# Load packages
lapply(packages, library, character.only = TRUE)

# Set options
options(digits = 8)  
options(scipen = 999)  

# Font settings for plots
# font_import(pattern = "Arial") # Run once if needed
loadfonts(device = "pdf")

# -------------------------------------------------------------------------
# Step 1: Identify Extreme Rainfall Events
# -------------------------------------------------------------------------

# 1.1 Load Data and Initial Filter
data <- read_csv("/Volumes/LaCie/00 data/data2223.csv")

# Filter for extreme events (Precipitation >= 50mm or Return Period >= 3 years)
events0 <- data %>%
  filter((prec >= 50 | returnperiod >= 3)) 

# 1.2 Temporal Grouping
# Temporally contiguous rainfall in the same city is classified as a single event
events1 <- events0 %>%
  arrange(citycode, date) %>%  
  group_by(citycode) %>%
  mutate(
    date_diff = date - lag(date, default = first(date) - 1), 
    group_id = cumsum(date_diff > 1) 
  ) %>%
  group_by(citycode, group_id) %>%
  slice_max(prec, n = 1, with_ties = FALSE) %>%  
  ungroup() %>%
  dplyr::select(-date_diff, -group_id) 

# 1.3 Spatial Grouping
# Spatially contiguous events within the same time window are grouped
adjacency_distance <- read_csv("/Volumes/LaCie/04 mobility/adjacency_distance.csv")
dis <- read_csv("/Volumes/LaCie/04 mobility/distance.csv")

# Prepare adjacency map
adjacency_map <- adjacency_distance %>%
  distinct(ocitycode, adjacent_codes) %>%
  separate_rows(adjacent_codes, sep = ",") %>%
  mutate(adjacent_citycode = as.numeric(adjacent_codes)) %>%
  rename(citycode = ocitycode)

# Mark temporal events
events0 <- events0 %>%
  arrange(citycode, date) %>%
  group_by(citycode) %>%
  mutate(
    event1 = cumsum(c(1, diff(date)) > 1)  
  ) %>%
  ungroup()

# Initialize spatial clustering
events0$event2 <- NA_integer_
event_counter <- 1

# 

# Flood-fill algorithm to identify spatio-temporal events
for (i in seq_len(nrow(events0))) {
  if (is.na(events0$event2[i])) {
    current_event1 <- events0$event1[i]
    current_city <- events0$citycode[i]
    
    same_event1_records <- events0 %>%
      filter(citycode == current_city, event1 == current_event1)
    
    queue <- same_event1_records %>%
      dplyr::select(date, citycode) %>%
      distinct()
    
    checked <- data.frame(date = as.Date(character()), citycode = numeric())
    
    while (nrow(queue) > 0) {
      current <- queue[1, ]
      queue <- queue[-1, ]
      
      # Skip if already checked
      if (nrow(checked %>% filter(date == current$date & citycode == current$citycode)) > 0) {
        next
      }
      checked <- rbind(checked, current)
      
      # Assign Event ID
      current_city_event1 <- events0 %>%
        filter(citycode == current$citycode) %>%
        filter(event1 == events0$event1[events0$citycode == current$citycode & events0$date == current$date])
      
      events0$event2[events0$citycode %in% current_city_event1$citycode & 
                       events0$event1 %in% current_city_event1$event1] <- event_counter
      
      # Find neighbors (Adjacent or within 200km)
      adjacent_cities <- adjacency_map %>%
        filter(citycode == current$citycode) %>%
        pull(adjacent_citycode)
      
      close_cities <- dis %>%
        filter(ocitycode == current$citycode & distance <= 200) %>%
        pull(dcitycode)
      
      related_cities <- unique(c(adjacent_cities, close_cities))
      
      # Add neighbors to queue if they have an event on the same day
      for (related_city in related_cities) {
        related_records <- events0 %>%
          filter(citycode == related_city, date == current$date, is.na(event2))
        
        if (nrow(related_records) > 0) {
          related_event1_records <- events0 %>%
            filter(citycode == related_city, event1 == related_records$event1[1], is.na(event2)) %>%
            dplyr::select(date, citycode) %>%
            distinct()
          
          queue <- rbind(queue, related_event1_records)
          queue <- unique(queue)
        }
      }
    }
    event_counter <- event_counter + 1
  }
}

# 1.4 Descriptive Statistics and Plots (Fig 1)
times <- events0 %>%
  filter(year == 2023) %>%
  group_by(citycode) %>%
  summarize(event_count = n(), .groups = "drop")
fwrite(times, '/Volumes/LaCie/new fig/fig 1a.csv', row.names = F)

# Count unique cities per event
event_counts <- events0 %>%
  group_by(event2) %>%
  summarise(unique_city_count = n_distinct(citycode), .groups = "drop")

com <- merge(events0, event_counts, by = 'event2')

# Count distribution of event sizes
event2_counts <- com %>%
  filter(year == 2023) %>%
  group_by(unique_city_count) %>%
  summarize(
    event2_count = n_distinct(event2),  
    row_count = n(),      
    .groups = "drop"
  )
fwrite(event2_counts, '/Volumes/LaCie/new fig/fig 1b.csv', row.names = F)

# Plot: Unique City Count vs Event Count (Fig 1c part 1)
p1 <- ggplot(event2_counts, aes(x = unique_city_count, y = event2_count)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  labs(x = "Unique City Count", y = "Event2 Count") +
  theme_classic() +
  theme(
    text = element_text(family = "Arial", size = 7),
    axis.title = element_text(size = 6.5),
    axis.text = element_text(size = 6.5),  
    plot.title = element_text(hjust = 0.5)
  ) + coord_flip()

ggsave(filename = "/Volumes/LaCie/new fig/fig1c1.pdf", plot = p1, width = 3.2, height = 5, units = "cm", device = cairo_pdf)

# Plot: Extreme rainfall duration distribution (Fig 1d)
event_days <- events0 %>%
  filter(year == 2023) %>%
  group_by(citycode, event2) %>%
  summarize(extreme_day_count = n(), .groups = "drop")

p3 <- ggplot(event_days, aes(x = extreme_day_count)) +
  geom_histogram(aes(y = after_stat(density)), binwidth = 1, color = "black", fill = "lightblue") +
  labs(x = "Extreme Rainfall Days", y = "Frequency") +
  theme_classic() +
  theme(
    text = element_text(family = "Arial", size = 7),
    axis.title = element_text(size = 6.5),
    axis.text = element_text(size = 6.5),  
    plot.title = element_text(hjust = 0.5)
  ) + coord_flip()

ggsave(filename = "/Volumes/LaCie/new fig/fig1d.pdf", plot = p3, width = 6, height = 3, units = "cm", device = cairo_pdf)


# -------------------------------------------------------------------------
# Step 2: Calculate Mobility Changes Relative to Baseline
# -------------------------------------------------------------------------

# 2.1 Filter Events for Analysis
# Exclude holidays and non-severe events
events <- events0 %>%
  filter((prec >= 50 | returnperiod >= 3) & holiday == 0) %>% 
  arrange(citycode, date) %>%  
  group_by(citycode) %>%
  mutate(event_group = cumsum(c(1, diff(as.Date(date)) > 1))) %>%
  group_by(citycode, event_group) %>%
  summarise(
    date = min(date),  # Keep start date
    returnperiod = first(returnperiod), 
    prec = max(prec),
    event2 = event2,
    .groups = "drop"
  ) %>%
  mutate(eventid = row_number()) 

data1 <- merge(data, events, by = c('citycode', "date", "returnperiod"), all.x = T)

# 2.2 Exclude Compound Hazards (14-day gap rule)
events <- data1 %>%
  mutate(date = as.Date(date)) %>%
  arrange(citycode, date) %>%
  filter(!is.na(eventid)) %>%
  group_by(citycode, eventid, event2) %>%
  summarise(start_date = min(date), .groups = "drop") %>%
  arrange(citycode, start_date) %>%
  group_by(citycode) %>%
  mutate(
    prev_start_date = lag(start_date),
    gap_prev = as.integer(start_date - prev_start_date)  
  ) %>%
  ungroup() %>%
  filter(is.na(prev_start_date) | gap_prev > 14)

# 2.3 Expand Time Window (-14 to +28 days)
data0 <- events %>%
  rowwise() %>%
  mutate(date_seq = list(seq(start_date - 14, start_date + 28, by = "day"))) %>%
  tidyr::unnest(cols = c(date_seq)) %>% 
  rename(date = date_seq) %>%
  ungroup() %>%
  mutate(date = as.Date(date)) %>%
  left_join(data %>% mutate(date = as.Date(date)), by = c("citycode", "date")) %>%
  group_by(citycode, eventid) %>%
  mutate(exposure = as.integer(date - start_date)) %>% # Relative day index
  ungroup()

# 2.4 Calculate Baseline
# Load NPI (Stringency Index) Data
npi <- read_csv("/Volumes/LaCie/10 NPI/ox/china npi.csv")
npi$Date <- as.Date(as.character(npi$Date), format = "%Y%m%d")
data0 <- merge(data0, npi, by.x = c('prov_code', 'date'), by.y = c('prov_code', 'Date'), all.x = T)
data0 <- data0 %>% dplyr::select(-...1, -RegionName)
data0$StringencyIndex_Average <- ifelse(is.na(data0$StringencyIndex_Average), 0, data0$StringencyIndex_Average)

# Define Baseline (Pre-event, no rain, no holiday)
base <- data0 %>%
  filter(prec <= 2 & exposure < 0 & holiday == 0) %>%
  group_by(citycode, eventid, event2) %>%
  summarise(
    windbase = mean(wind_speed),
    tempbase = mean(max_temp),
    npibase = mean(StringencyIndex_Average),
    # Differentiate baseline between weekdays (1-5) and weekends (0,6)
    base_day_6_and_0 = ifelse(
      !is.na(mean(intra[day %in% c(6, 0)], na.rm = TRUE)), 
      mean(intra[day %in% c(6, 0)], na.rm = TRUE),
      mean(intra[day %in% 1:5], na.rm = TRUE)
    ),
    base_day_1_to_5 = ifelse(
      !is.na(mean(intra[day %in% 1:5], na.rm = TRUE)),
      mean(intra[day %in% 1:5], na.rm = TRUE),
      mean(intra[day %in% c(6, 0)], na.rm = TRUE)
    ),
    .groups = "drop"
  )

data0 <- data0 %>%
  mutate(selected_base = ifelse(day %in% 1:5, "base_day_1_to_5",
                                ifelse(day %in% c(6, 0), "base_day_6_and_0", NA)))

data0 <- data0 %>%
  left_join(base %>% pivot_longer(cols = starts_with("base_day"), 
                                  names_to = "selected_base", 
                                  values_to = "base_value"), 
            by = c("citycode", "eventid", "selected_base"))

# Calculate percentage change in mobility (intrac)
data0$intrac <- data0$intra / data0$base_value - 1

# Merge event characteristics back
data0 <- merge(data0, event_day_count, by.x = c('citycode', 'event2.x'), by.y = c('citycode', 'event2'))
data0 <- merge(data0, event_counts, by.x = c('event2.x'), by.y = c('event2'), all.x = T)

# 2.5 Aggregate and Plot Trajectories
avg_all <- data0 %>%
  filter(year == 2023) %>%
  group_by(exposure) %>%
  summarise(
    mean_intrac = mean(intrac, na.rm = TRUE),
    q25 = quantile(intrac, 0.25, na.rm = TRUE),
    q75 = quantile(intrac, 0.75, na.rm = TRUE)
  )

# Calculate subgroups for plotting (1-day vs 2-day duration)
avg_d1 <- data0 %>% filter(year == 2023, event_day_count == 1) %>%
  group_by(exposure) %>%
  summarise(mean_intrac = mean(intrac, na.rm = T), q25 = quantile(intrac, 0.25, na.rm = T), q75 = quantile(intrac, 0.75, na.rm = T))

avg_d2 <- data0 %>% filter(year == 2023, event_day_count == 2) %>%
  group_by(exposure) %>%
  summarise(mean_intrac = mean(intrac, na.rm = T), q25 = quantile(intrac, 0.25, na.rm = T), q75 = quantile(intrac, 0.75, na.rm = T))

# Plot Mobility Trajectory (Fig 1d_2)
p_traj <- ggplot(data0, aes(x = exposure, y = intrac, group = interaction(citycode, eventid))) +
  geom_line(color = "blue", alpha = 0.05) + 
  
  # All Data
  geom_ribbon(data = avg_all, aes(x = exposure, ymin = q25, ymax = q75), fill = "gray70", alpha = 0.3, inherit.aes = FALSE) +
  geom_line(data = avg_all, aes(x = exposure, y = mean_intrac), color = "black", size = 0.6, inherit.aes = FALSE) +
  
  # 1-Day Events
  geom_ribbon(data = avg_d1, aes(x = exposure, ymin = q25, ymax = q75), fill = "green", alpha = 0.3, inherit.aes = FALSE) +
  geom_line(data = avg_d1, aes(x = exposure, y = mean_intrac), color = "green4", size = 0.6, inherit.aes = FALSE) +
  
  # 2-Day Events
  geom_ribbon(data = avg_d2, aes(x = exposure, ymin = q25, ymax = q75), fill = "orange", alpha = 0.3, inherit.aes = FALSE) +
  geom_line(data = avg_d2, aes(x = exposure, y = mean_intrac), color = "orange3", size = 0.6, inherit.aes = FALSE) +
  
  labs(x = "Exposure (day)", y = "Baseline (%)") +
  xlim(-7, 7) + ylim(-0.52, 0.27) +
  theme_bw() +
  theme(text = element_text(family = "Arial", size = 7), plot.title = element_text(hjust = 0.5))

ggsave(filename = "/Volumes/LaCie/new fig/fig1d_2.pdf", plot = p_traj, width = 7.5, height = 5, units = "cm", device = cairo_pdf)


# 2.6 Inflow Analysis (Supplementary)
inflowbase <- data0 %>%
  filter(prec <= 2 & exposure < 0 & holiday == 0) %>%
  group_by(citycode, eventid) %>%
  summarise(
    base_day_6_and_0 = ifelse(!is.na(mean(inflow[day %in% c(6, 0)], na.rm = T)), mean(inflow[day %in% c(6, 0)], na.rm = T), mean(inflow[day %in% 1:5], na.rm = T)),
    base_day_1_to_5 = ifelse(!is.na(mean(inflow[day %in% 1:5], na.rm = T)), mean(inflow[day %in% 1:5], na.rm = T), mean(inflow[day %in% c(6, 0)], na.rm = T))
  )

datainflow <- data0 %>%
  left_join(inflowbase %>% pivot_longer(cols = starts_with("base_day"), names_to = "selected_base", values_to = "base_value"),
            by = c("citycode", "eventid", "selected_base"))

datainflow$inflowc <- datainflow$inflow / datainflow$base_value.y - 1

# Plot Inflow Change
plot_datainflow <- datainflow %>%
  filter(year == 2023) %>%
  group_by(exposure) %>%
  summarise(
    mean_inflowc = 100 * median(inflowc, na.rm = TRUE),
    lower_bound = 100 * quantile(inflowc, 0.25, na.rm = TRUE),
    upper_bound = 100 * quantile(inflowc, 0.75, na.rm = TRUE)
  )

p_inflow <- ggplot(plot_datainflow, aes(x = exposure, y = mean_inflowc)) +
  geom_ribbon(aes(ymin = lower_bound, ymax = upper_bound), fill = "#2C3E50", alpha = 0.25) +
  geom_line(color = "#2C3E50", size = 0.8) +
  labs(x = "Exposure day", y = "Inflow change relative to baseline (%)") +
  xlim(-7, 8) + theme_bw() +
  theme(text = element_text(family = "Arial", size = 7))

ggsave(filename = "/Volumes/LaCie/new fig/inflowchange.pdf", plot = p_inflow, width = 5, height = 5, units = "cm", device = cairo_pdf)

# -------------------------------------------------------------------------
# Step 3: Calculate Mobility Resilience Metrics
# -------------------------------------------------------------------------
# 

data2 <- data0[data0$exposure >= -3, ]

result0 <- data2 %>%
  group_by(citycode, event2.x) %>%
  summarise(
    base = mean(base_value),
    date_exposure_zero = min(date[exposure == 0], na.rm = TRUE),
    
    # Identify the day when mobility first falls below the baseline
    date_before = {
      d1 <- date[exposure <= 0 & exposure >= -3 & intrac <= 0]
      d2 <- date[exposure <= 0 & exposure >= -3 & intrac > 0]
      cand1 <- if (length(d1) > 0) min(d1) else as.Date(NA)
      cand2 <- if (length(d2) > 0) max(d2) else as.Date(NA)
      max(cand1, cand2, na.rm = TRUE)
    },
    
    # Identify the day of recovery
    date_after = min(date[exposure > 0 & exposure <= 25 & intrac >= 0], na.rm = TRUE),
    
    # Impacted period duration
    period = ifelse(!is.na(date_before) & !is.na(date_after), as.numeric(date_after - date_before), NA),
    
    # Recovery time
    recovery = ifelse(!is.na(date_before) & !is.na(date_after), as.numeric(date_after - date_exposure_zero), NA),
    
    # Maximum impact (Max drop)
    max_intrac = ifelse(!is.na(date_before) & !is.na(date_after), min(intrac[date >= date_before & date < date_after], na.rm = TRUE), NA),
    
    # Total performance loss (Cumulative drop)
    total_intrac = ifelse(!is.na(date_before) & !is.na(date_after), sum(intrac[date >= date_before & date < date_after], na.rm = TRUE), NA),
    
    # Sensitivity analysis version of Total Loss
    total_intrac2 = ifelse(!is.na(date_before) & !is.na(date_after), sum(intrac[date >= date_exposure_zero & date < date_after], na.rm = TRUE), NA),
    
    # Event characteristics during impact
    max_prec = ifelse(!is.na(date_before) & !is.na(date_after), max(prec[date >= date_before & date < date_after], na.rm = TRUE), NA),
    
    # Environmental changes
    max_wind_diff = ifelse(!is.na(date_before) & !is.na(date_after), max(wind_speed[date >= date_before & date < date_after] - windbase[date >= date_before & date < date_after], na.rm = TRUE), NA)
  ) %>%
  ungroup()

# Filter valid results
result <- result0[!is.na(result0$base), ]
result <- result %>% filter(is.finite(period))
overall <- result %>% filter(total_intrac < 0)

# Merge metadata
overall <- merge(overall, event_day_count, by = c('citycode', "event2.x"), by.y = c('citycode', 'event2'), all.x = TRUE)
overall <- merge(overall, event_counts, by = c('event2.x'), by.y = c('event2'), all.x = TRUE)
overall$year <- year(overall$date_before)

# 3.1 Group Comparisons and Visualization
overall_long1 <- bind_rows(
  overall %>% mutate(group = "All"),
  overall %>% filter(event_day_count == 1) %>% mutate(group = "D1"),
  overall %>% filter(event_day_count > 1) %>% mutate(group = "D>1")
) %>% filter(year == 2023)

overall_long1$total_intrac <- (-1) * overall_long1$total_intrac # Convert loss to positive for plotting

# Plot: Distribution of Total Mobility Loss (Fig 1e)
p_box <- ggplot(overall_long1, aes(x = group, y = total_intrac, fill = group)) +
  geom_half_violin(side = "r", alpha = 0.5, color = NA, width = 1.0) +
  geom_boxplot(width = 0.15, outlier.shape = NA, alpha = 0.8) +
  stat_summary(fun = mean, geom = "point", shape = 21, size = 1.5, fill = "red", color = "black") +
  labs(x = "", y = "Total Intracity Mobility Loss") +
  coord_cartesian(ylim = c(0, 1.2)) +
  theme_bw() +
  theme(text = element_text(family = "Arial", size = 6.5), legend.position = "none")

ggsave(filename = "/Volumes/LaCie/new fig/fig1e.pdf", plot = p_box, width = 6.5, height = 5, units = "cm", device = cairo_pdf)

# 3.2 Regression and Relative Importance Analysis
# Calculate derived metrics
overall <- overall %>% filter(year == 2023)
overall$shape <- (-1) * overall$total_intrac2 / ((overall$max_intrac * overall$recovery) / 2)
overall$total_intrac <- (-1) * overall$total_intrac # Make positive for log
overall$max_intrac <- abs(overall$max_intrac)

# Log transformation for regression
overall$logtpl <- log(overall$total_intrac)
overall$logmi <- log(overall$max_intrac)
overall$logshape <- log(overall$shape)
overall$logrecovery <- log(overall$recovery)

# Linear Models
model <- lm(logshape ~ logmi + logrecovery, data = overall)
summary(model)

model2 <- lm(logtpl ~ logmi + logrecovery + logshape, data = overall)
summary(model2)

# Relative Importance (LMG)
metrics <- calc.relimp(model2, type = "lmg", rela = TRUE)

lmg_results <- data.frame(
  Factor = c("Maximum Impact (Resistance)", "Recovery Time (Rapidity)", "Trajectory Shape (Efficiency)"),
  Contribution = c(metrics$lmg["logmi"] * 100, 
                   metrics$lmg["logrecovery"] * 100, 
                   metrics$lmg["logshape"] * 100)
)
print(lmg_results)

# 3.3 Statistical Tests (Wilcoxon)
group1 <- overall$total_intrac[overall$event_day_count == 1]
group2 <- overall$total_intrac[overall$event_day_count > 1]

wilcox_test_result <- wilcox.test(group1, group2)
print(wilcox_test_result)
print(cohens_d(group1, group2))