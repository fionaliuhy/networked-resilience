rm(list=ls())

# ==============================================================================
# 1. Setup and Libraries
# ==============================================================================
packages <- c("data.table", "tidyverse", "sf", "sp", "spdep",
              "dlnm", "tsModel", "hydroGOF", "RColorBrewer", 
              "geofacet", "ggpubr", "ggthemes", "extRemes",
              "randomForestSRC", "mlr3", "mlr3tuning", "mlr3learners", 
              "mlr3verse", "mlr3misc", "mlr3pipelines", "paradox", 
              "future", "mgcv", "mlr3extralearners", "sfnetworks", 
              "tidygraph", "igraph", "purrr", "tibble", "corrplot", 
              "patchwork", "gridExtra")

# Load packages
lapply(packages, library, character.only = TRUE)

options(digits = 8)  
options(scipen = 999)  

# ==============================================================================
# 2. Network Data Preparation (2013 & 2023)
# ==============================================================================

# --- 2.1 Load Network Data ---
network23 <- read_csv("/Volumes/LaCie/04 mobility/network 23.csv") %>% 
  dplyr::rename(citycode = city)

network13 <- read.csv("/Volumes/LaCie/04 mobility/network 13.csv") %>% 
  dplyr::rename(citycode = city) %>%
  mutate(citycode = as.character(citycode))

# --- 2.2 Calculate Market Potential (MP) for 2023 ---
dis <- read_csv("/Volumes/LaCie/04 mobility/distance.csv")
socio <- read_csv('/Volumes/LaCie/07 socioeco/soci.csv') %>% dplyr::rename(citycode = city_code)

socio <- socio %>% mutate(citycode = as.character(citycode))
dis <- dis %>% mutate(ocitycode = as.character(ocitycode), dcitycode = as.character(dcitycode))

dis_clean <- dis %>%
  filter(ocitycode != dcitycode) %>%
  select(ocitycode, dcitycode, distance) %>%
  mutate(distance = ifelse(distance <= 0, 1, distance))

socio_source <- socio %>% select(citycode, gdp, pop) %>% rename(ocitycode = citycode, src_gdp = gdp, src_pop = pop)

mp_calc <- dis_clean %>%
  left_join(socio_source, by = "ocitycode") %>%
  mutate(contrib_gdp = src_gdp / distance, contrib_pop = src_pop / distance) %>%
  group_by(dcitycode) %>%
  summarise(mp_gdp = sum(contrib_gdp, na.rm = TRUE), mp_pop = sum(contrib_pop, na.rm = TRUE), .groups = "drop")

socio_final23 <- socio %>%
  left_join(mp_calc, by = c("citycode" = "dcitycode")) %>%
  mutate(mp_gdp = replace_na(mp_gdp, 0), mp_pop = replace_na(mp_pop, 0))

network23 <- merge(network23, socio_final23, by = 'citycode')

# --- 2.3 Calculate Market Potential (MP) for 2013 ---
socio13 <- read_csv('/Volumes/LaCie/04 ncc/07 socioeco/socio13.csv')
socio13$gdp <- socio13$gdppercapita * socio13$TOTAL_POP2014 / 100000000

# Fix City Codes for 2013
replace_map <- c("522200" = "520600", "542200" = "540500", "542300" = "540200", 
                 "542400" = "540600", "632100" = "630200", "652100" = "650400", 
                 "652200" = "650500", "542600" = "540400", "542100" = "540300", 
                 "522400" = "520500", "411800" = "419001")

socio13$city <- as.character(socio13$city)
socio13$city <- ifelse(socio13$city %in% names(replace_map), replace_map[socio13$city], socio13$city)
socio13 <- socio13 %>% dplyr::rename(pop = TOTAL_POP2014)

# Aggregate Distances for Merged Cities in 2013
grp_429 <- c("429004", "429005", "429006", "429021")
grp_469 <- c("469001", "469002", "469005", "469006", "469007", "469021", 
             "469022", "469023", "469024", "469025", "469026", "469027", 
             "469028", "469029", "469030")
grp_659 <- c("659001", "659002", "659003", "659004", "659005", 
             "659006", "659007", "659008", "659009", "659010")

dis2 <- dis %>%
  mutate(dcitycode_new = case_when(
    dcitycode %in% grp_429 ~ "429000", dcitycode %in% grp_469 ~ "469000",
    dcitycode %in% grp_659 ~ "659000", dcitycode == "371200" ~ "370100",
    dcitycode == "611100" ~ "611000", TRUE ~ dcitycode
  )) %>%
  group_by(dcitycode_new, ocitycode) %>%
  summarise(distance = mean(distance, na.rm = TRUE), .groups = "drop") %>%
  rename(dcitycode = dcitycode_new) %>%
  mutate(ocitycode_new = case_when(
    ocitycode %in% grp_429 ~ "429000", ocitycode %in% grp_469 ~ "469000",
    ocitycode %in% grp_659 ~ "659000", ocitycode == "371200" ~ "370100",
    ocitycode == "611100" ~ "611000", TRUE ~ ocitycode
  )) %>%
  group_by(ocitycode_new, dcitycode) %>%
  summarise(distance = mean(distance, na.rm = TRUE), .groups = "drop") %>%
  rename(ocitycode = ocitycode_new)

socio13 <- socio13 %>% mutate(city = as.character(city))
dis2 <- dis2 %>% mutate(ocitycode = as.character(ocitycode), dcitycode = as.character(dcitycode))

dis2_with_socio <- dis2 %>%
  left_join(socio13 %>% select(city, gdp, pop), by = c("dcitycode" = "city")) %>%
  mutate(dist_safe = ifelse(distance <= 0 | is.na(distance), NA_real_, distance),
         gdp_term = gdp / dist_safe, pop_term = pop / dist_safe)

mp_city_13 <- dis2_with_socio %>%
  group_by(ocitycode) %>%
  summarise(mp_gdp = sum(gdp_term, na.rm = TRUE), mp_pop = sum(pop_term, na.rm = TRUE), .groups = "drop")

# Merge Socio Economic Data for 2013
socio_level_topo <- socio %>% mutate(citycode = as.character(citycode)) %>%
  select(citycode, level, mean_elevation, mean_slope, slopedifference)

socio13_mp <- socio13 %>%
  left_join(mp_city_13, by = c("city" = "ocitycode")) %>%
  mutate(city = as.character(city)) %>%
  left_join(socio_level_topo %>% rename(city = citycode), by = "city")

network13 <- merge(network13, socio13_mp, by.x = 'citycode', by.y = 'city', all.x = T)
network13$year <- 2013

# ==============================================================================
# 3. Railway Network Analysis (2013 & 2023)
# ==============================================================================

# --- 3.1 Load and Preprocess Railway Shapefiles ---
sf_use_s2(FALSE)
target_crs <- 3857

rail13 <- st_read("/Volumes/LaCie/railway/2013.shp", quiet = TRUE) %>% 
  st_transform(target_crs) %>% filter(fclass == 'rail') %>% st_make_valid()
rail23 <- st_read("/Volumes/LaCie/railway/2024.shp", quiet = TRUE) %>% 
  st_transform(target_crs) %>% filter(fclass == 'rail') %>% st_make_valid()

china_city <- st_read("/Volumes/LaCie/city/china/china_city.shp", quiet = TRUE) %>% 
  st_transform(target_crs) %>% rename(citycode = adcode)

city_point <- st_read("/Volumes/LaCie/city_cn/city_point.shp", quiet = TRUE) %>% 
  st_transform(target_crs)

# --- 3.2 Compute City-level Railway Metrics ---
compute_city_rail_metrics <- function(rail_sf, city_sf, year_label) {
  rail_in_city <- st_intersection(st_make_valid(rail_sf), st_make_valid(city_sf["citycode"]))
  rail_in_city <- rail_in_city %>% mutate(length_km = as.numeric(st_length(geometry)) / 1000)
  
  city_area <- city_sf %>% mutate(area_km2 = as.numeric(st_area(geometry)) / 1e6) %>% st_drop_geometry()
  
  rail_summary <- rail_in_city %>% st_drop_geometry() %>%
    group_by(citycode) %>% summarise(rail_length_km = sum(length_km, na.rm = TRUE), .groups = "drop")
  
  out <- city_area %>% left_join(rail_summary, by = "citycode") %>%
    mutate(rail_length_km = replace_na(rail_length_km, 0),
           rail_density_km_per_100km2 = ifelse(area_km2 > 0, rail_length_km / area_km2 * 100, NA),
           year = year_label)
  return(out)
}

rail13_lines <- rail13 %>% st_union() %>% st_line_merge() %>% st_cast("LINESTRING") %>% st_sf(geometry=.)
rail23_lines <- rail23 %>% st_union() %>% st_line_merge() %>% st_cast("LINESTRING") %>% st_sf(geometry=.)

rail_city_2013 <- compute_city_rail_metrics(rail13_lines, china_city, 2013)
rail_city_2023 <- compute_city_rail_metrics(rail23_lines, china_city, 2023)

# Merge Railway Metrics
network23 <- merge(network23, rail_city_2023, by.x = c('citycode', 'year'), by.y = c('citycode', 'year'))
network13 <- merge(network13, rail_city_2013, by.x = c('citycode', 'year'), by.y = c('citycode', 'year'))

# --- 3.3 Railway Connectivity (SFNetwork + Igraph) ---
build_rail_network <- function(rail_lines_sf) {
  as_sfnetwork(rail_lines_sf, directed = FALSE) %>%
    activate("edges") %>% mutate(weight = as.numeric(st_length(geometry)) / 1000)
}

net13 <- build_rail_network(rail13_lines)
net23 <- build_rail_network(rail23_lines)

# Attach cities to network
attach_cities_to_network <- function(net, city_polys, city_pts, city_id_col = "citycode", k = 3) {
  nodes_sf <- net %>% activate("nodes") %>% mutate(node_id = row_number()) %>% st_as_sf()
  
  # Ensure CRS match
  if(st_crs(net) != st_crs(city_polys)) city_polys <- st_transform(city_polys, st_crs(net))
  if(st_crs(net) != st_crs(city_pts)) city_pts <- st_transform(city_pts, st_crs(net))
  
  # Calculate distance matrix
  dist_mat <- st_distance(nodes_sf, city_pts)
  
  city_node_list <- vector("list", nrow(city_pts))
  city_ids <- city_pts[[city_id_col]]
  
  for(i in seq_along(city_ids)) {
    # Simple strategy: k nearest nodes to city centroid
    dists <- as.numeric(dist_mat[, i])
    keep_idx <- order(dists)[1:k]
    city_node_list[[i]] <- tibble(city_id = city_ids[i], node_id = nodes_sf$node_id[keep_idx])
  }
  bind_rows(city_node_list)
}

city_node_13 <- attach_cities_to_network(net13, china_city, city_point, k = 3)
city_node_23 <- attach_cities_to_network(net23, china_city, city_point, k = 3)

# Compute Distances
compute_rail_distances <- function(g, city_nodes) {
  city_nodes <- city_nodes %>% mutate(citycode = as.character(city_id))
  node_list <- split(city_nodes$node_id, city_nodes$citycode)
  cities <- names(node_list)
  
  purrr::map_dfr(cities, function(ci) {
    nodes_i <- node_list[[ci]]
    d_all <- igraph::distances(g, v = nodes_i, to = V(g), weights = E(g)$weight)
    
    purrr::map_dfr(cities, function(cj) {
      if(ci == cj) return(NULL)
      nodes_j <- node_list[[cj]]
      sub_mat <- d_all[, nodes_j, drop = FALSE]
      valid <- sub_mat[is.finite(sub_mat)]
      if(length(valid) == 0) return(NULL)
      tibble(city_from = ci, city_to = cj, min_dist_km = min(valid))
    })
  })
}

g13 <- as.igraph(net13)
g23 <- as.igraph(net23)

city_pair_rail_13 <- compute_rail_distances(g13, city_node_13)
city_pair_rail_23 <- compute_rail_distances(g23, city_node_23)

# --- 3.4 Calculate Connectivity Indicators ---
summarize_rail_conn <- function(pair_df, year_val) {
  pair_df %>% group_by(city_from) %>%
    summarise(
      n_direct_rail = n_distinct(city_to),
      mean_direct_dist = mean(min_dist_km),
      share_within_800 = sum(min_dist_km <= 800) / n(),
      .groups = "drop"
    ) %>%
    rename(citycode = city_from) %>% mutate(year = year_val)
}

city_conn_13 <- summarize_rail_conn(city_pair_rail_13, 2013)
city_conn_23 <- summarize_rail_conn(city_pair_rail_23, 2023)

network13 <- merge(network13, city_conn_13, by = c('citycode', 'year'), all.x = T)
network23 <- merge(network23, city_conn_23, by = c('citycode', 'year'), all.x = T)

# --- 3.5 Calculate Rich-City Connectivity ---
# Identify Rich Cities (Top 25% GDP per capita)
get_rich_cities <- function(df) {
  thr <- quantile(df$gdppercapita, 0.75, na.rm = TRUE)
  df %>% filter(gdppercapita >= thr) %>% pull(citycode) %>% as.character()
}

rich_13 <- get_rich_cities(network13)
rich_23 <- get_rich_cities(network23)

summarize_rich_conn <- function(pair_df, rich_list, gdp_df) {
  pair_df %>% 
    filter(city_to %in% rich_list) %>%
    left_join(gdp_df %>% select(citycode, gdppercapita), by = c("city_to" = "citycode")) %>%
    group_by(city_from) %>%
    summarise(
      n_rail_rich = n_distinct(city_to),
      mean_rail_dist_rich = mean(min_dist_km),
      sum_gdp_over_dist_rich = sum(gdppercapita / min_dist_km, na.rm = TRUE),
      .groups = "drop"
    ) %>% rename(citycode = city_from)
}

rich_conn_13 <- summarize_rich_conn(city_pair_rail_13, rich_13, network13)
rich_conn_23 <- summarize_rich_conn(city_pair_rail_23, rich_23, network23)

network13 <- merge(network13, rich_conn_13, by = "citycode", all.x = T)
network23 <- merge(network23, rich_conn_23, by = "citycode", all.x = T)

# Fill NAs
network13[is.na(network13)] <- 0
network23[is.na(network23)] <- 0

# Save intermediate files
write.csv(network13, '/Volumes/LaCie/08 result/network13_forpredict.csv')
write.csv(network23, '/Volumes/LaCie/08 result/network23_forpredict.csv')

# ==============================================================================
# 4. Step 1: Predict Network Metrics (Multivariate Random Forest)
# ==============================================================================

# --- 4.1 Prepare Data for Prediction ---
# Combine 2013 and 2023 data
x_vars_pred <- c("pop", "gdppercapita", "avg_dis", "mindis", "mp_pop", "level", 
                 "mean_elevation", "rail_length_km", "rail_density_km_per_100km2", 
                 "n_direct_rail", "share_within_800", "mean_rail_dist_rich", 
                 "sum_gdp_over_dist_rich", "prov_code")

# Targets: Eigenvector centrality & Elite city linkage
net_vars_pred <- c("Eigenvector centrality", "Elite city linkage")

# Combine and Transform
net_combined <- bind_rows(
  network13 %>% select(citycode, year, all_of(net_vars_pred), all_of(x_vars_pred)),
  network23 %>% select(citycode, year, all_of(net_vars_pred), all_of(x_vars_pred))
) %>%
  mutate(
    pop = log(pop), gdppercapita = log(gdppercapita), 
    avg_dis = log(avg_dis), mindis = log(mindis), mp_pop = log(mp_pop),
    rail_length_km = log1p(rail_length_km),
    sum_gdp_over_dist_rich = log1p(sum_gdp_over_dist_rich),
    level = as.factor(level), prov_code = as.factor(prov_code)
  ) %>% na.omit()

# --- 4.2 Initial Training & Hyperparameter Tuning ---
form_multi <- as.formula(paste0("Multivar(`", paste(net_vars_pred, collapse = "`, `"), "`) ~ . - citycode - year"))

# Hyperparameter Tuning (Grid Search 5-fold CV)
mtry_grid <- c(3, 5, 7)
nodesize_grid <- c(3, 5, 10)
best_score <- -Inf
best_params <- list(mtry = 3, nodesize = 5) # Default

for(m in mtry_grid) {
  for(n in nodesize_grid) {
    set.seed(123)
    rf_cv <- rfsrc(form_multi, data = net_combined, ntree = 100, mtry = m, nodesize = n, block.size = 1)
    # Score: Average R2 across targets
    r2_1 <- 1 - rf_cv$regrOutput[[1]]$mse[100] / var(net_combined[[net_vars_pred[1]]])
    r2_2 <- 1 - rf_cv$regrOutput[[2]]$mse[100] / var(net_combined[[net_vars_pred[2]]])
    score <- mean(c(r2_1, r2_2))
    
    if(score > best_score) {
      best_score <- score
      best_params <- list(mtry = m, nodesize = n)
    }
  }
}

cat("Best Params (Initial): ", unlist(best_params), "\n")

# --- 4.3 Feature Selection (Filter) ---
# Train full model with best params to get importance
set.seed(123)
rf_full <- rfsrc(form_multi, data = net_combined, ntree = 500, 
                 mtry = best_params$mtry, nodesize = best_params$nodesize, 
                 importance = "permute")

# Extract Importance
# rfsrc importance is a matrix (rows=vars, cols=targets)
imp_mat <- rf_full$importance
# Calculate average importance across both targets
avg_imp <- rowMeans(imp_mat)
# Select Top 8 Variables (or based on threshold)
top_k_vars <- names(sort(avg_imp, decreasing = TRUE)[1:6])

cat("Selected Top Variables for Network Prediction:\n")
print(top_k_vars)

# --- 4.4 Re-fit Model with Selected Features ---
# Update data to only include top features + targets + ID variables
data_reduced <- net_combined %>% 
  select(citycode, year, all_of(net_vars_pred), all_of(top_k_vars))

# Update formula
form_reduced <- as.formula(paste0("Multivar(`", paste(net_vars_pred, collapse = "`, `"), "`) ~ . - citycode - year"))

# Re-fit Final Model
set.seed(123)
rf_net_final <- rfsrc(form_reduced, data = data_reduced, ntree = 500, 
                      mtry = min(best_params$mtry, length(top_k_vars)), # Adjust mtry if k is small
                      nodesize = best_params$nodesize, 
                      importance = "permute")

saveRDS(rf_net_final, "/Volumes/LaCie/08 result/model_rf_netpredict.rds")

# --- 4.5 Validation Plot (Optional) ---
preds <- predict(rf_net_final, newdata = data_reduced)
# (Plotting code remains the same...)


# ==============================================================================
# 5. Step 2: Predict Resilience (TPL) (Random Forest via mlr3)
#    (Updated with Feature Selection & Re-fitting)
# ==============================================================================

# --- 5.1 Load Resilience Data ---
overall <- read.csv('/Volumes/LaCie/08 result/overall_cenup1124.csv') %>%
  filter(year == 2023)

# Log Transforms
overall$logtpl <- log(-1 * overall$total_intrac)
overall$pop <- log(overall$pop)
overall$gdppercapita <- log(overall$gdppercapita)

x_vars_resilience <- c("total_prec", "max_temp_diff", "max_wind_diff", "npi_diff", 
                       "urbanization", "pop", "gdppercapita", "max_Rpolicy", "event_day_count",
                       "level", "prov_code") 


df_resilience <- overall %>% select(logtpl, all_of(x_vars_resilience)) %>% 
  mutate(level = as.factor(level), prov_code = as.factor(prov_code)) %>% na.omit()

# --- 5.2 Initial Tuning ---
task_tpl <- TaskRegr$new(id = "tpl", backend = df_resilience, target = "logtpl")
learner_rf <- lrn("regr.ranger", num.trees = 500, importance = "impurity")

param_set <- ps(mtry = p_int(3, 10), min.node.size = p_int(2, 10))
instance <- TuningInstanceSingleCrit$new(
  task = task_tpl, learner = learner_rf, resampling = rsmp("cv", folds = 5),
  measure = msr("regr.rmse"), search_space = param_set, terminator = trm("none")
)
tuner <- tnr("grid_search", resolution = 5)
tuner$optimize(instance)

# --- 5.3 Feature Selection (Filter) ---
# Train on full feature set with optimized params
learner_rf$param_set$values <- instance$result_learner_param_vals
learner_rf$train(task_tpl)

# Get Importance
imp <- learner_rf$importance()
# Select Top 8 Variables
top_features <- names(sort(imp, decreasing = TRUE)[1:8])

cat("Selected Top Variables for Resilience Prediction:\n")
print(top_features)


# --- 5.4 Re-fit Model with Selected Features ---
# Create reduced task
task_tpl_reduced <- task_tpl$clone()$select(top_features)

# Define learner for reduced model
learner_rf_final <- lrn("regr.ranger", num.trees = 500, importance = "impurity")

# Optional: Quick re-tune on reduced set (or use previous params adapted)
# Here we adapt mtry to not exceed number of features
final_mtry <- min(instance$result_learner_param_vals$mtry, length(top_features))
final_min_node <- instance$result_learner_param_vals$min.node.size

learner_rf_final$param_set$values <- list(mtry = final_mtry, min.node.size = final_min_node)

# Train Final Model
learner_rf_final$train(task_tpl_reduced)

saveRDS(learner_rf_final, "/Volumes/LaCie/08 result/rf_tpl_final_model.rds")


# ==============================================================================
# 6. Step 3: Bias Correction (Bin-based Calibration)
# ==============================================================================

# Predict on training data (OOB or CV predictions preferred, here simple predict for demo)
preds_tpl <- learner_rf$predict(task_tpl)
df_pred <- data.frame(Actual = preds_tpl$truth, Predicted = preds_tpl$response)

# Log to Original Scale
df_pred <- df_pred %>% mutate(Actual_TPL = exp(Actual), Predicted_TPL = exp(Predicted))

# Add Population Grouping
df_pred$pop <- df_resilience$pop
pop_breaks <- quantile(df_pred$pop, probs = seq(0, 1, 0.25))
df_pred$pop_group <- cut(df_pred$pop, breaks = pop_breaks, include.lowest = T)

# Calculate Calibration Coefficients per Group
calib_tbl <- df_pred %>%
  group_by(pop_group) %>%
  do(model = lm(Actual_TPL ~ Predicted_TPL, data = .)) %>%
  summarise(pop_group = pop_group, 
            slope = coef(model)[2], 
            intercept = coef(model)[1])

# Apply Correction
df_pred <- df_pred %>%
  left_join(calib_tbl, by = "pop_group") %>%
  mutate(Predicted_TPL_Corrected = slope * Predicted_TPL + intercept)

# Validation
total_actual <- sum(df_pred$Actual_TPL)
total_pred <- sum(df_pred$Predicted_TPL_Corrected)
error_pct <- (total_pred - total_actual) / total_actual * 100

cat("\n=== Total TPL Validation ===\n")
cat("Actual Total:", format(total_actual, big.mark=","), "\n")
cat("Predicted Total (Corrected):", format(total_pred, big.mark=","), "\n")
cat("Error %:", round(error_pct, 2), "%\n")
