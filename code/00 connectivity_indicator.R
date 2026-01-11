rm(list=ls())

# -------------------------------------------------------------------------
# 1. Load Packages
# -------------------------------------------------------------------------
library(tidyverse)
library(igraph)
library(readr)
library(sf)
library(spdep)

# -------------------------------------------------------------------------
# 2. Define Helper Functions
# -------------------------------------------------------------------------

# Function to calculate Gini coefficient manually
gini_manual <- function(x) {
  x <- sort(na.omit(x))
  n <- length(x)
  G <- sum((2 * (1:n) - n - 1) * x)
  # Handle case where sum(x) is 0 to avoid division by zero
  if (sum(x) == 0) return(0)
  return(G / (n * sum(x)))
}

# -------------------------------------------------------------------------
# 3. Load and Preprocess Data
# -------------------------------------------------------------------------

# --- 3.1 Load Main Mobility Data ---
# Replace with your specific year/file path
data_raw <- read_csv("/Volumes/LaCie/04 mobility/2023/2023.csv") 

# Filter dates (exclude spring festival)
data_filtered <- data_raw %>%
  filter(time > 20230215 | time < 20230107)

# Basic cleaning: Remove NAs and self-loops
data_clean <- data_filtered %>%
  filter(!is.na(ocitycode) & !is.na(dcitycode) & !is.na(mobility)) %>%
  filter(ocitycode != dcitycode) %>%
  mutate(
    ocitycode = as.numeric(ocitycode),
    dcitycode = as.numeric(dcitycode)
  )

# --- 3.2 Load Auxiliary Data ---
# Distance Data
dis <- read_csv("/Volumes/LaCie/04 mobility/distance.csv")

# Socio-economic Data (GDP and Population)
socieco <- read_csv('/Volumes/LaCie/07 socioeco/soci.csv')
socieco <- socieco[, c(1, 9, 12, 4)]
colnames(socieco) <- c('citycode', 'gdppercapita', "level", "pop")
socieco$citycode <- as.numeric(socieco$citycode)

# Adjacency Data
adjacency_distance <- read_csv("/Volumes/LaCie/04 mobility/adjacency_distance.csv")

# --- 3.3 Merge and Aggregate ---
# Merge distance first
data_merged <- data_clean %>%
  left_join(dis, by = c("ocitycode", "dcitycode"))

# Aggregate mobility between city pairs (summing over time)
# Note: Adjust scaling factor (e.g., *324/143) here if comparing different time window lengths
data_agg <- data_merged %>%
  group_by(ocitycode, dcitycode) %>%
  summarise(
    mobility = sum(mobility, na.rm = TRUE),
    distance = first(distance), # Assuming distance doesn't change
    drive_distance = first(drive_distance),
    .groups = "drop"
  )

# Merge Socio-economic data for Origin (x) and Destination (y)
data_final <- data_agg %>%
  left_join(socieco, by = c("ocitycode" = "citycode")) %>%
  left_join(socieco, by = c("dcitycode" = "citycode"), suffix = c(".x", ".y")) %>%
  mutate(
    gdppercapita.x = as.numeric(gdppercapita.x),
    gdppercapita.y = as.numeric(gdppercapita.y),
    pop.x = as.numeric(pop.x)
  )

# -------------------------------------------------------------------------
# 4. Calculation of 15 Metrics
# -------------------------------------------------------------------------

# --- Part A: Network Centrality (Graph Theory) ---
# Create Graph Object
g <- graph_from_data_frame(
  data_final[, c("ocitycode", "dcitycode", "mobility", "drive_distance")], 
  directed = TRUE
)

# 1. Pagerank
pagerank_scores <- page_rank(g, weights = E(g)$mobility)$vector

# 2. Eigenvector Centrality
eigenvector_scores <- eigen_centrality(g, weights = log1p(E(g)$mobility))$vector

# 3. Closeness Centrality (Distance based)
# Using inverse of distance or direct distance depending on interpretation. 
# Standard interpretation for spatial interaction:
closeness_scores <- closeness(g, weights = E(g)$drive_distance, normalized = TRUE)

df_centrality <- data.frame(
  city = as.numeric(names(pagerank_scores)),
  `Pagerank` = pagerank_scores,
  `Eigenvector centrality` = eigenvector_scores,
  `Closeness centrality` = closeness_scores,
  check.names = FALSE
)

# --- Part B: Basic Mobility & Weighted Averages ---

# 4. Mobility per capita
# 10. Mobility-weighted distance
# 5. Mobility-weighted GDP per capita
df_basic <- data_final %>%
  group_by(ocitycode) %>%
  summarise(
    # Mobility per capita
    total_mobility = sum(mobility, na.rm = TRUE),
    population = first(pop.x),
    `Mobility per capita` = total_mobility / population, 
    
    # Mobility-weighted distance
    # Note: Averaging the weighted distance of outgoing and incoming (via join later) or just outgoing.
    # Based on your code "distance_preferenceo", this is usually outgoing average:
    `Mobility-weighted distance` = sum(mobility * distance, na.rm = TRUE) / sum(mobility, na.rm = TRUE),
    
    # Mobility-weighted GDP per capita (Origin perspective)
    w_gdp_out = sum(mobility * gdppercapita.y, na.rm = TRUE) / sum(mobility, na.rm = TRUE),
    .groups = "drop"
  )

# Calculate Destination perspective for Weighted GDP to average them (as per your logic)
df_basic_dest <- data_final %>%
  group_by(dcitycode) %>%
  summarise(
    w_gdp_in = sum(mobility * gdppercapita.x, na.rm = TRUE) / sum(mobility, na.rm = TRUE),
    .groups = "drop"
  )

# Merge and calculate final Weighted GDP
df_basic <- df_basic %>%
  left_join(df_basic_dest, by = c("ocitycode" = "dcitycode")) %>%
  mutate(
    `Mobility-weighted GDP per capita` = (w_gdp_out + w_gdp_in) / 2
  ) %>%
  select(ocitycode, `Mobility per capita`, `Mobility-weighted distance`, `Mobility-weighted GDP per capita`)


# --- Part C: Inequality & Concentration (Within Group) ---

# 12. Intercity interaction unevenness (MAD)
# 13. Maximum connected ratio
# 14. Intercity interaction inequality (Gini)
df_inequality <- data_final %>%
  group_by(ocitycode) %>%
  mutate(
    total_mob = sum(mobility),
    ratio = mobility / total_mob
  ) %>%
  summarise(
    median_ratio = median(ratio, na.rm = TRUE),
    max_ratio = max(ratio, na.rm = TRUE),
    
    `Intercity interaction unevenness` = mean(abs(ratio - median_ratio), na.rm = TRUE),
    `Maximum connected ratio` = abs(max_ratio - median_ratio),
    `Intercity interaction inequality` = gini_manual(ratio),
    .groups = "drop"
  )

# --- Part D: Thresholds & Administrative Linkages ---

# 15. Within province linkage
df_province <- data_final %>%
  mutate(
    oprov = substr(ocitycode, 1, 2),
    dprov = substr(dcitycode, 1, 2)
  ) %>%
  group_by(ocitycode) %>%
  summarise(
    within_mob = sum(mobility[oprov == dprov], na.rm = TRUE),
    total_mob = sum(mobility, na.rm = TRUE),
    `Within province linkage` = ifelse(total_mob > 0, within_mob / total_mob, NA),
    .groups = "drop"
  )

# 11. Elite city linkage (Ratio of flow with high GDP cities)
gdp_threshold <- quantile(as.numeric(socieco$gdppercapita), 0.75, na.rm = TRUE)

df_elite <- data_final %>%
  # We need to look at both In and Out flow for "Linkage" usually, 
  # but here we aggregate by Origin (O) and Destination (D) separately then combine
  pivot_longer(cols = c(ocitycode, dcitycode), names_to = "type", values_to = "citycode") %>%
  group_by(citycode) %>%
  summarise(
    # If type is ocitycode, we look at gdppercapita.y (destination gdp)
    # If type is dcitycode, we look at gdppercapita.x (origin gdp)
    target_gdp = ifelse(type == "ocitycode", gdppercapita.y, gdppercapita.x),
    flow = mobility
  ) %>%
  # Filter only valid rows created by pivoting
  filter(!is.na(citycode)) %>%
  group_by(citycode) %>%
  summarise(
    elite_flow = sum(flow[target_gdp > gdp_threshold], na.rm = TRUE),
    total_flow = sum(flow, na.rm = TRUE),
    `Elite city linkage` = ifelse(total_flow > 0, elite_flow / total_flow, 0),
    .groups = "drop"
  )

# --- Part E: Neighborhood / Adjacency Metrics ---

# Prepare Adjacency Data
adj_clean <- adjacency_distance %>%
  distinct(ocitycode, adjacent_codes) %>%
  mutate(adjacent_codes = strsplit(as.character(adjacent_codes), ",")) %>%
  unnest_longer(adjacent_codes) %>%
  mutate(adjacent_codes = as.numeric(adjacent_codes))

# Join Main Data with Adjacency to get flows specifically to neighbors
# We join data_final where dcitycode is the neighbor
data_neighbors <- adj_clean %>%
  left_join(data_final, by = c("ocitycode" = "ocitycode", "adjacent_codes" = "dcitycode"))

# 6. Near-neighbor preference
# 7. Economic standing
# 8. Economic bias
# 9. Neighbor linkage

df_neighbors <- data_neighbors %>%
  group_by(ocitycode) %>%
  summarise(
    # -- Data prep --
    total_neighbor_mobility = sum(mobility, na.rm = TRUE),
    avg_neighbor_gdp = mean(gdppercapita.y, na.rm = TRUE), # Average GDP of neighbors
    avg_neighbor_dist = mean(drive_distance, na.rm = TRUE), # Average distance to neighbors
    self_gdp = mean(gdppercapita.x, na.rm = TRUE),
    
    # Weighted preferences
    pref_gdp = sum(mobility * gdppercapita.y, na.rm = TRUE) / sum(mobility, na.rm = TRUE),
    pref_dist = sum(mobility * drive_distance, na.rm = TRUE) / sum(mobility, na.rm = TRUE),
    
    .groups = "drop"
  ) %>%
  mutate(
    `Near-neighbor preference` = pref_dist / avg_neighbor_dist,
    `Economic standing` = pref_gdp / avg_neighbor_gdp,
    `Economic bias` = self_gdp / pref_gdp
  )

# For Neighbor Linkage (Total flow to neighbors / Total flow everywhere)
# We need Total Flow from df_basic and Neighbor Flow from above
df_neighbor_linkage <- df_neighbors %>%
  select(ocitycode, total_neighbor_mobility) %>%
  left_join(df_basic %>% select(ocitycode, total_mobility = `Mobility per capita`), by = "ocitycode") %>%
  # Note: `Mobility per capita` col in df_basic currently holds ratio, let's re-fetch total
  left_join(data_final %>% group_by(ocitycode) %>% summarise(tot = sum(mobility)), by="ocitycode") %>%
  mutate(
    `Neighbor linkage` = total_neighbor_mobility / tot
  ) %>%
  select(ocitycode, `Neighbor linkage`)

# -------------------------------------------------------------------------
# 5. Merge All Metrics
# -------------------------------------------------------------------------

final_metrics <- df_centrality %>%
  left_join(df_basic, by = c("city" = "ocitycode")) %>%
  left_join(df_inequality, by = c("city" = "ocitycode")) %>%
  left_join(df_province, by = c("city" = "ocitycode")) %>%
  left_join(df_elite, by = c("city" = "citycode")) %>%
  left_join(df_neighbors %>% select(ocitycode, `Near-neighbor preference`, `Economic standing`, `Economic bias`), 
            by = c("city" = "ocitycode")) %>%
  left_join(df_neighbor_linkage, by = c("city" = "ocitycode"))

# -------------------------------------------------------------------------
# 6. Final Cleaning and Export
# -------------------------------------------------------------------------

# Select exactly the 15 metrics + City ID
final_output <- final_metrics %>%
  select(
    city,
    `Closeness centrality`,
    `Mobility per capita`,
    `Pagerank`,
    `Eigenvector centrality`,
    `Mobility-weighted GDP per capita`,
    `Near-neighbor preference`,
    `Economic standing`,
    `Economic bias`,
    `Neighbor linkage`,
    `Mobility-weighted distance`,
    `Elite city linkage`,
    `Intercity interaction unevenness`,
    `Maximum connected ratio`,
    `Intercity interaction inequality`,
    `Within province linkage`
  )

# Handle potential Infinite values or NaNs from divisions
final_output <- final_output %>%
  mutate(across(where(is.numeric), ~ ifelse(is.infinite(.), NA, .))) %>%
  mutate(across(where(is.numeric), ~ ifelse(is.nan(.), NA, .)))

# Export
write_csv(final_output, "/Volumes/LaCie/04 mobility/network_metrics_final.csv")

print("Calculation complete. File saved.")
head(final_output)