rm(list=ls())

# -------------------------------------------------------------------------
# 1. Setup and Library Loading
# -------------------------------------------------------------------------
packages <- c("data.table", "tidyverse", "sf", "sp", "spdep",
              "dlnm", "tsModel", "hydroGOF","RColorBrewer", 
              "geofacet", "ggpubr", "ggthemes","extRemes")

# Load packages
lapply(packages, library, character.only = TRUE)

options(digits = 8)  
options(scipen = 999) 

# -------------------------------------------------------------------------
# 2. Data Loading
# -------------------------------------------------------------------------
# Load main OD (Origin-Destination) data for 2023
data_raw <- read_csv("/Volumes/LaCie/04 mobility/2023/2023.csv") 

# Load aggregated time-series data
data_ts <- read_csv("/Volumes/LaCie/00 data/data2223.csv")

# -------------------------------------------------------------------------
# 3. Case Study 1: Zhengzhou (City Code: 410100)
# -------------------------------------------------------------------------

# --- 3.1 Time Series Visualization (Single-City Event) ---
data_sub <- data_ts %>%
  filter(citycode == 410100) %>%
  filter(date >= as.Date("2023-07-29"), date <= as.Date("2023-08-11")) %>%
  select(date, inflow, outflow, netflow) %>%
  pivot_longer(cols = c(inflow, outflow, netflow), names_to = "flow_type", values_to = "mobility")

p1 <- ggplot(data_sub, aes(x = date, y = mobility, color = flow_type)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 1.5) +
  geom_vline(xintercept = as.Date("2023-08-05"), linetype = "dashed", color = "black", linewidth = 0.6) +
  scale_color_manual(
    values = c(inflow = "#D55E00", outflow = "#009E73", netflow = "#0072B2"),
    labels = c(inflow = "Inflow", outflow = "Outflow", netflow = "Net flow"),
    name = ""
  ) +
  theme_bw(base_size = 12) +
  labs(x = "Date", y = "Mobility index", title = "Zhengzhou (Single-City Event)") +
  theme(text = element_text(family = "Arial", size = 7), plot.title = element_text(hjust = 0.5))

print(p1)

# --- 3.2 Statistical Test (Wilcoxon) for Event 1 ---
event_date <- as.Date("2023-08-05")
data_test <- data_ts %>%
  filter(citycode == 410100) %>%
  filter(date > event_date - 7, date <= event_date + 6, date != event_date) %>%
  mutate(period = ifelse(date < event_date, "Before", "After"))

results_zz1 <- data_test %>%
  select(inflow, outflow, netflow, period) %>%
  pivot_longer(cols = c(inflow, outflow, netflow), names_to = "flow_type", values_to = "mobility") %>%
  group_by(flow_type) %>%
  summarise(
    p_value = wilcox.test(mobility ~ period, exact = FALSE)$p.value,
    median_before = mean(mobility[period == "Before"], na.rm = TRUE),
    median_after  = mean(mobility[period == "After"],  na.rm = TRUE)
  )
print(results_zz1)

# --- 3.3 Time Series Visualization (Multi-City Event) ---
data_sub_multi <- data_ts %>%
  filter(citycode == 410100) %>%
  filter(date >= as.Date("2023-08-19"), date <= as.Date("2023-09-01")) %>%
  select(date, inflow, outflow, netflow) %>%
  pivot_longer(cols = c(inflow, outflow, netflow), names_to = "flow_type", values_to = "mobility")

p2 <- ggplot(data_sub_multi, aes(x = date, y = mobility, color = flow_type)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 1.5) +
  geom_vline(xintercept = as.Date("2023-08-26"), linetype = "dashed", color = "black", linewidth = 0.6) +
  scale_color_manual(
    values = c(inflow = "#D55E00", outflow = "#009E73", netflow = "#0072B2"),
    labels = c(inflow = "Inflow", outflow = "Outflow", netflow = "Net flow")
  ) +
  theme_bw(base_size = 12) +
  labs(x = "Date", y = "Intercity mobility", title = "Zhengzhou (Multi-City Event)")

print(p2)

# --- 3.4 Statistical Test for Event 2 ---
event_date <- as.Date("2023-08-26")
data_test <- data_ts %>%
  filter(citycode == 410100) %>%
  filter(date > event_date - 7, date <= event_date + 6, date != event_date) %>%
  mutate(period = ifelse(date < event_date, "Before", "After"))

results_zz2 <- data_test %>%
  select(inflow, outflow, netflow, period) %>%
  pivot_longer(cols = c(inflow, outflow, netflow), names_to = "flow_type", values_to = "mobility") %>%
  group_by(flow_type) %>%
  summarise(
    p_value = wilcox.test(mobility ~ period, exact = FALSE)$p.value,
    median_before = median(mobility[period == "Before"], na.rm = TRUE),
    median_after  = median(mobility[period == "After"],  na.rm = TRUE)
  )
print(results_zz2)



# -------------------------------------------------------------------------
# 4. Comparative Analysis: Single vs. Multi-City Events
# -------------------------------------------------------------------------
# Filter main OD data for Zhengzhou
od_zz <- data_raw %>%
  filter(ocitycode %in% c(410100) | dcitycode %in% c(410100)) %>%
  mutate(mobility = replace_na(mobility, 0))

# --- Zhengzhou Comparisons ---

# Baseline and Change Calculation (Single Event)
base_410100 <- od_zz %>%
  filter(dcitycode == 410100, time >= "20230722" & time <= "20230804") %>%
  group_by(ocitycode) %>%
  summarize(base = mean(mobility), normalratio = base / 704.53805)

in_410100s <- od_zz %>%
  filter(dcitycode == 410100, time >= "20230805" & time <= "20230807") %>%
  group_by(ocitycode) %>%
  summarize(avg_mobility = mean(mobility), .groups = "drop") %>%
  left_join(base_410100, by = "ocitycode") %>%
  mutate(
    mobility_change = avg_mobility - base,
    changeratio = mobility_change / 857.31936
  )

# Normalize to percentage
in_410100s$normalratio <- in_410100s$normalratio * 100
in_410100s$changeratio <- in_410100s$changeratio * 100

# Baseline and Change Calculation (Multi Event)
base_410100_m <- od_zz %>%
  filter(dcitycode == 410100, time >= "20230812" & time <= "20230825") %>%
  group_by(ocitycode) %>%
  summarize(base = mean(mobility), normalratio = base / 704.53805)

in_410100m <- od_zz %>%
  filter(dcitycode == 410100, time >= "20230826" & time <= "20230828") %>%
  group_by(ocitycode) %>%
  summarize(avg_mobility = mean(mobility), .groups = "drop") %>%
  left_join(base_410100_m, by = "ocitycode") %>%
  mutate(
    mobility_change = avg_mobility - base,
    changeratio = mobility_change / 704.53805
  )

in_410100m$normalratio <- in_410100m$normalratio * 100
in_410100m$changeratio <- in_410100m$changeratio * 100

# Visualization: Comparison Scatter Plot (Zhengzhou)
combined_data_zz <- bind_rows(
  in_410100m %>% mutate(source = "in_410100m"),
  in_410100s %>% mutate(source = "in_410100s")
)

# Plotting Comparison
ggplot(combined_data_zz, aes(x = normalratio, y = changeratio)) +
  geom_point(data = subset(combined_data_zz, source == "in_410100s"), color = "#1b7837", size = 2, alpha = 0.6) +
  geom_point(data = subset(combined_data_zz, source == "in_410100m"), color = "#0571b0", size = 2, alpha = 0.6) +
  geom_smooth(data = subset(combined_data_zz, source == "in_410100s"), method = "lm", se = TRUE, color = "#1b7837") +
  geom_smooth(data = subset(combined_data_zz, source == "in_410100m"), method = "lm", se = TRUE, color = "#0571b0") +
  labs(x = "Normal Ratio (%)", y = "Change Ratio (%)", title = "Zhengzhou Comparison") +
  theme_bw()

# Correlation between Single and Multi Event Baselines (Zhengzhou)
df_zz <- data.frame(
  x = in_410100s$normalratio,
  y = in_410100m$normalratio
)

figs21 <- ggplot(df_zz, aes(x = x, y = y)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, color = "blue") +
  labs(x = "Baseline inflow in single-city event", y = "Baseline inflow in multi-city event", title = "Zhengzhou Baseline Correlation") +
  theme_bw() + theme(text = element_text(size = 7))

cor.test(in_410100s$normalratio, in_410100m$normalratio)

