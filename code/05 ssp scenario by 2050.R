rm(list=ls())

# ==============================================================================
# 1. Setup and Libraries
# ==============================================================================
packages <- c("dplyr", "randomForest", "data.table", "tidyr", "stats", "readr", 
              "mlr3", "mlr3learners", "xgboost", "ggplot2")

# Load packages
lapply(packages, library, character.only = TRUE)

options(digits = 8)
options(scipen = 999)

# Common Paths
base_path <- "/Volumes/LaCie/"

# ==============================================================================
# 2. Load Future Climate and Socioeconomic Data
# ==============================================================================

# --- 2.1 Future Climate Data ---
files_final <- list.files(paste0(base_path, "01 prec/size"), full.names = TRUE)
future <- do.call(rbind, lapply(files_final, read.csv))

# Calculate unique city count per event
future_event_count <- future %>%
  group_by(scenario, model, year, event2) %>%
  summarise(unique_city_count = n_distinct(citycode), .groups = "drop")

future <- merge(future, future_event_count, by = c('scenario', 'model', 'year', 'event2'), all.x = T)

# Summarize precipitation metrics (used later in MC Stage 2)
future1 <- future %>%
  group_by(citycode, scenario, year, model, event2, unique_city_count) %>%
  summarise(
    total_prec = sum(prec, na.rm = TRUE),
    max_prec = max(prec, na.rm = TRUE),
    .groups = "drop"
  )

# --- 2.2 Future Population and GDP ---
pop <- read.csv(paste0(base_path, "future ssp/ssp_pop_old.csv")) %>%
  rename(citycode = city_code) %>%
  mutate(citycode = as.character(ifelse(citycode == "5e+05", "500000", citycode)))

gdp <- read.csv(paste0(base_path, "future ssp/GDP_city_filled.csv")) %>%
  rename(citycode = city_code) %>%
  mutate(citycode = as.character(ifelse(citycode == "5e+05", "500000", citycode))) %>%
  select(citycode, year, GDP, scenario)

# Merge Pop and GDP
pop <- merge(pop, gdp, by = c('citycode', 'year', 'scenario'))
pop$gdppercapita <- 100000000 * pop$GDP / pop$pop
pop$GDP <- pop$GDP * 100

# --- 2.3 Baseline Network Data (2023) ---
stat23 <- read_csv(paste0(base_path, "08 result/network23_forpredict.csv")) %>%
  mutate(citycode = as.character(ifelse(citycode == "5e+05", "500000", citycode))) %>%
  select(citycode, pop, urbanization, old65, gdppercapita, level, prov_code,
         mp_pop, richshare, n_direct_rail, sum_gdp_over_dist_rich, 
         rail_length_km, rail_density_km_per_100km2)

# ==============================================================================
# 3. Bias Correction for Socioeconomic Projections
# ==============================================================================
# Calculate bias based on 2023 actuals vs projections and apply to future years

# 1. Calculate 2023 Bias
pop_2023_baseline <- pop %>%
  filter(year == 2023) %>%
  select(citycode, scenario, pop_2023_proj = pop, gdp_2023_proj = gdppercapita)

bias_table <- pop_2023_baseline %>%
  left_join(stat23 %>% select(citycode, pop_real = pop, gdp_real = gdppercapita), by = "citycode") %>%
  mutate(
    diff_pop = pop_real - pop_2023_proj,
    diff_gdp = gdp_real - gdp_2023_proj
  ) %>%
  select(citycode, scenario, diff_pop, diff_gdp)

# 2. Apply Bias Correction
pop_corrected <- pop %>%
  left_join(bias_table, by = c("citycode", "scenario")) %>%
  mutate(
    popc = pmax(0, pop + diff_pop),
    gdppercapitac = pmax(0, gdppercapita + diff_gdp)
  ) %>%
  select(-diff_pop, -diff_gdp)

# ==============================================================================
# 4. Calculate Future Market Potential (MP)
# ==============================================================================

dis <- read_csv(paste0(base_path, "04 mobility/distance.csv")) %>%
  mutate(ocitycode = as.character(ifelse(ocitycode == "5e+05", "500000", ocitycode)),
         dcitycode = as.character(ifelse(dcitycode == "5e+05", "500000", dcitycode)))

dis_clean <- dis %>%
  filter(ocitycode != dcitycode) %>%
  select(ocitycode, dcitycode, distance) %>%
  mutate(distance = ifelse(distance <= 0, 1, distance))

# Prepare Source Data
socio_source <- pop_corrected %>%
  select(citycode, scenario, year, popc) %>%
  rename(ocitycode = citycode, src_pop = popc)

# Calculate MP
mp_calc <- dis_clean %>%
  left_join(socio_source, by = "ocitycode") %>%
  mutate(contrib_pop = src_pop / distance) %>%
  group_by(dcitycode, scenario, year) %>%
  summarise(mp_pop = sum(contrib_pop, na.rm = TRUE), .groups = "drop")

# Merge MP back to Population Data
popnew <- merge(pop_corrected, mp_calc, by.x = c('citycode', 'scenario', 'year'), by.y = c('dcitycode', 'scenario', 'year'), all.x = T)

# Merge static variables from 2023
stat1 <- stat23 %>%
  select(citycode, level, prov_code, richshare, n_direct_rail, sum_gdp_over_dist_rich, rail_length_km)

popnew$gdp_day <- (popnew$gdppercapitac * popnew$popc / 1000000) / 365
popnew <- merge(popnew, stat1, by = 'citycode', all.x = T)

# ==============================================================================
# 5. Urbanization Projection
# ==============================================================================

# Load Provincial Urbanization Targets
urban1 <- read.csv(paste0(base_path, "future ssp/DATA_Provincial_Urbanization_Rate 2/urbanShareSSP2.csv")) %>%
  rename(year = X) %>%
  pivot_longer(cols = -year, names_to = "Province_English", values_to = "urbanization") %>%
  mutate(scenario = 'ssp245') %>% filter(year > 2019)

urban2 <- read.csv(paste0(base_path, "future ssp/DATA_Provincial_Urbanization_Rate 2/urbanShareSSP5.csv")) %>%
  rename(year = X) %>%
  pivot_longer(cols = -year, names_to = "Province_English", values_to = "urbanization") %>%
  mutate(scenario = 'ssp585') %>% filter(year > 2019)

urban <- rbind(urban1, urban2) %>%
  mutate(year = as.numeric(year), urbanization = as.numeric(urbanization))

# 2020 Baseline Urbanization
urban20 <- stat23 %>% select(citycode, urbanization) %>% rename(urbanization_2020 = urbanization)

# Function to calculate city-level urbanization iteratively
calculate_city_urbanization <- function(popnew, urban, urban20) {
  # Pre-processing
  df_main <- popnew %>%
    mutate(year = as.numeric(year), popc = as.numeric(popc)) %>%
    arrange(scenario, citycode, year) %>%
    left_join(urban %>% rename(urban_rate_prov = urbanization), by = c("Province_English", "scenario", "year")) %>%
    left_join(urban20, by = "citycode")
  
  df_main$urban_city <- NA_real_
  df_main <- df_main %>%
    mutate(urban_city = ifelse(year == 2020, urbanization_2020, urban_city)) %>%
    select(-urbanization_2020)
  
  years <- sort(unique(df_main$year))
  future_years <- years[years > 2020]
  
  for (yr in future_years) {
    # T-1 Data
    df_prev <- df_main %>% filter(year == (yr - 1)) %>%
      select(citycode, scenario, Province_English, pop_prev = popc, urban_rate_prev = urban_city) %>%
      mutate(city_urban_pop_prev = pop_prev * urban_rate_prev)
    
    # T-1 Provincial Sum
    prov_sum_prev <- df_prev %>%
      group_by(Province_English, scenario) %>%
      summarise(prov_urban_pop_prev = sum(city_urban_pop_prev, na.rm = TRUE), .groups = "drop")
    
    # T Data
    df_curr <- df_main %>% filter(year == yr) %>%
      select(citycode, scenario, Province_English, pop_curr = popc, urban_rate_prov_curr = urban_rate_prov)
    
    # T Provincial Target
    prov_sum_curr <- df_curr %>%
      group_by(Province_English, scenario) %>%
      summarise(
        prov_total_pop_curr = sum(pop_curr, na.rm = TRUE),
        urban_rate_prov_curr = mean(urban_rate_prov_curr, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      mutate(prov_urban_pop_target = prov_total_pop_curr * urban_rate_prov_curr)
    
    # Calculation
    calc_df <- df_curr %>%
      left_join(df_prev %>% select(citycode, scenario, city_urban_pop_prev, urban_rate_prev), by = c("citycode", "scenario")) %>%
      left_join(prov_sum_prev, by = c("Province_English", "scenario")) %>%
      left_join(prov_sum_curr %>% select(Province_English, scenario, prov_urban_pop_target), by = c("Province_English", "scenario")) %>%
      mutate(
        delta_prov_urban = prov_urban_pop_target - prov_urban_pop_prev,
        share = if_else(prov_urban_pop_prev > 0, city_urban_pop_prev / prov_urban_pop_prev, 0),
        city_urban_pop_curr = city_urban_pop_prev + (delta_prov_urban * share),
        raw_urban_rate = if_else(pop_curr > 0, city_urban_pop_curr / pop_curr, 0),
        final_urban_rate = pmin(pmax(pmax(raw_urban_rate, urban_rate_prev, na.rm = TRUE), 0), 1)
      )
    
    update_dat <- calc_df %>% mutate(year = yr) %>% select(citycode, scenario, year, final_urban_rate)
    
    # Update main dataframe
    df_main <- df_main %>%
      left_join(update_dat, by = c("citycode", "scenario", "year")) %>%
      mutate(urban_city = coalesce(final_urban_rate, urban_city)) %>%
      select(-final_urban_rate)
  }
  return(df_main)
}

result <- calculate_city_urbanization(popnew, urban, urban20)
result$urbanization <- result$urban_city
result$n_direct_rail <- replace_na(result$n_direct_rail, 0)
result$sum_gdp_over_dist_rich <- replace_na(result$sum_gdp_over_dist_rich, 0)

# Merge Industry Structure Data (Lambda)
industry_city_ssp <- read.csv(paste0(base_path, "future ssp/industry_city_ssp.csv")) %>%
  select(citycode, year, scenario, lambda1, lambda2, lambda3)

result0 <- merge(result, industry_city_ssp, by = c('citycode', 'year', 'scenario'), all.x = T)

# Log Transformations for Prediction Model
result0$pop <- log(result0$popc)
result0$gdppercapita <- log(result0$gdppercapitac)
result0$mp_pop <- log(result0$mp_pop)
result0$rail_length_km <- log1p(result0$rail_length_km)
result0$sum_gdp_over_dist_rich <- log1p(result0$sum_gdp_over_dist_rich)

# Save intermediate result
write.csv(result0, paste0(base_path, 'future ssp/ssp_overall0.csv'), row.names = FALSE)

# ==============================================================================
# 6. Network Growth Simulation (Scenario Simulation)
# ==============================================================================

# Parameters
target_n <- 365       # Target direct rail connections
target_rail1 <- 0.02  # Standard annual growth rate
target_rail2 <- 0.05  # Catch-up annual growth rate
thr_len <- 785.45     # Threshold for rail length
thr_gdp <- 7945.7     # Threshold for GDP potential
min_rail_seed <- 1
min_gdp_seed <- 1671

# Simulation Logic (Catch-up Mode)
result2 <- result0 %>%
  group_by(citycode, scenario) %>%
  group_modify(~{
    df <- arrange(.x, year)
    yrs <- df$year
    
    if (any(yrs == 2023)) {
      idx_2023 <- which(yrs == 2023)[1]
      
      # 1) Direct Rail Connections (Linear Growth to 2035)
      start_n <- replace_na(df$n_direct_rail[idx_2023], 0)
      idx_23_35 <- which(yrs >= 2023 & yrs <= 2035)
      
      if (length(idx_23_35) > 0 && start_n < target_n) {
        n_seq <- seq(from = start_n, to = target_n, length.out = length(idx_23_35))
        df$n_direct_rail[idx_23_35] <- as.integer(round(n_seq))
      }
      
      idx_after_35 <- which(yrs > 2035)
      if (length(idx_after_35) > 0 && start_n < target_n) {
        df$n_direct_rail[idx_after_35] <- target_n
      }
      
      # 2) Rich Share (Linear Growth to 1 by 2050)
      start_rich <- replace_na(df$richshare[idx_2023], 0)
      idx_23_50 <- which(yrs >= 2023 & yrs <= 2050)
      
      if (length(idx_23_50) > 0 && start_rich < 1) {
        rich_seq <- seq(from = start_rich, to = 1, length.out = length(idx_23_50))
        df$richshare[idx_23_50] <- pmin(pmax(rich_seq, 0), 1)
      }
      idx_after_50 <- which(yrs > 2050)
      if (length(idx_after_50) > 0) df$richshare[idx_after_50] <- 1
      
      # 3) Rail Length (Segmented Exponential Growth)
      base_rail <- replace_na(df$rail_length_km[idx_2023], 0)
      start_val_rail <- ifelse(base_rail == 0, min_rail_seed, base_rail)
      idx_future <- which(yrs > 2023)
      
      if (length(idx_future) > 0) {
        if (start_val_rail < thr_len) {
          # Catch-up mode
          years_from_23 <- yrs[idx_future] - 2023
          fast_curve <- start_val_rail * (1 + target_rail2)^years_from_23
          catch_idx <- which(fast_curve >= thr_len)[1]
          
          if (is.na(catch_idx)) {
            df$rail_length_km[idx_future] <- fast_curve
          } else {
            df$rail_length_km[idx_future[1:catch_idx]] <- fast_curve[1:catch_idx]
            if (catch_idx < length(idx_future)) {
              idx_rem <- idx_future[(catch_idx + 1):length(idx_future)]
              years_rem <- yrs[idx_rem] - yrs[idx_future[catch_idx]]
              df$rail_length_km[idx_rem] <- fast_curve[catch_idx] * (1 + target_rail1)^years_rem
            }
          }
        } else {
          # Standard mode
          df$rail_length_km[idx_future] <- start_val_rail * (1 + target_rail1)^(yrs[idx_future] - 2023)
        }
      }
      
      # 4) GDP Potential (Segmented Exponential Growth)
      base_gdp <- replace_na(df$sum_gdp_over_dist_rich[idx_2023], 0)
      start_val_gdp <- ifelse(base_gdp == 0, min_gdp_seed, base_gdp)
      
      if (length(idx_future) > 0) {
        if (start_val_gdp < thr_gdp) {
          years_from_23 <- yrs[idx_future] - 2023
          fast_curve <- start_val_gdp * (1 + target_rail2)^years_from_23
          catch_idx <- which(fast_curve >= thr_gdp)[1]
          
          if (is.na(catch_idx)) {
            df$sum_gdp_over_dist_rich[idx_future] <- fast_curve
          } else {
            df$sum_gdp_over_dist_rich[idx_future[1:catch_idx]] <- fast_curve[1:catch_idx]
            if (catch_idx < length(idx_future)) {
              idx_rem <- idx_future[(catch_idx + 1):length(idx_future)]
              years_rem <- yrs[idx_rem] - yrs[idx_future[catch_idx]]
              df$sum_gdp_over_dist_rich[idx_rem] <- fast_curve[catch_idx] * (1 + target_rail1)^years_rem
            }
          }
        } else {
          df$sum_gdp_over_dist_rich[idx_future] <- start_val_gdp * (1 + target_rail1)^(yrs[idx_future] - 2023)
        }
      }
    }
    df
  }) %>% ungroup()

# Apply log1p to simulated growth variables
result2$rail_length_km <- log1p(result2$rail_length_km)
result2$sum_gdp_over_dist_rich <- log1p(result2$sum_gdp_over_dist_rich)

write.csv(result2, paste0(base_path, 'future ssp/ssp_overall2.csv'), row.names = FALSE)

# ==============================================================================
# 7. Predict Network Metrics & Calibrate (Random Forest)
# ==============================================================================

# Load Models and Params
model2 <- readRDS(paste0(base_path, "08 result/model_rf_netpredict.rds"))
dml_stab <- read.csv(paste0(base_path, "08 result/updatecen/log/multi_dml_en_tpl_rep_5factor.csv"))
network23 <- read_csv(paste0(base_path, "04 mobility/network 23 new2.csv")) %>% dplyr::rename(citycode = city)

# Parameters for Stratified Calibration
pop_breaks <- c(-Inf, 13.815511, 15.424948, 16.118096, Inf)
group_params <- tibble(
  pop_group = 1:4,
  eig_slope = c(0.656, 0.944, 0.882, 1.219), eig_int = c(0.021, 0.007, 0.042, -0.146),
  wr_slope = c(1.362, 1.265, 1.275, 1.307), wr_int = c(-0.098, -0.119, -0.133, -0.16)
)

# Extract DML Coefficients
beta_mat <- dml_stab %>% filter(D %in% c("eigenvector", "withrich")) %>%
  select(rep, D, coef) %>% pivot_wider(names_from = D, values_from = coef) %>%
  arrange(rep) %>% select(eigenvector, withrich) %>% as.matrix()
n_beta <- nrow(beta_mat)

# --- 7.1 Predict Raw Network Metrics ---
pred_obj <- predict(model2, newdata = result2)

raw_pred_net <- result2 %>%
  mutate(
    citycode = as.character(citycode), scenario = as.character(scenario),
    eigen_pred_raw = pred_obj$regrOutput$eigenvector$predicted,
    withrich_pred_raw = pred_obj$regrOutput$withrich$predicted,
    pop_group = cut(pop, breaks = pop_breaks, labels = FALSE, include.lowest = TRUE)
  ) %>%
  select(year, scenario, citycode, pop_group, eigen_pred_raw, withrich_pred_raw)

# --- 7.2 Calculate 2023 Residual Bias ---
network23_use <- network23 %>% mutate(citycode = as.character(citycode)) %>%
  select(citycode, eigen_23 = eigenvector, withrich_23 = withrich)

bias_calc <- raw_pred_net %>%
  filter(year == 2023) %>% distinct(citycode, scenario, .keep_all = TRUE) %>%
  left_join(network23_use, by = "citycode") %>%
  left_join(group_params, by = "pop_group") %>%
  mutate(
    bias_eig = eigen_23 - (eig_int + eig_slope * eigen_pred_raw),
    bias_wr = withrich_23 - (wr_int + wr_slope * withrich_pred_raw)
  ) %>%
  select(citycode, scenario, bias_eig, bias_wr)

# --- 7.3 Apply Calibration to All Years ---
future_net_base <- raw_pred_net %>%
  left_join(group_params, by = "pop_group") %>%
  left_join(bias_calc, by = c("citycode", "scenario")) %>%
  left_join(network23_use, by = "citycode") %>%
  mutate(
    bias_eig = replace_na(bias_eig, 0), bias_wr = replace_na(bias_wr, 0),
    eigen_mean_calib = pmin(1, pmax(0, eig_int + eig_slope * eigen_pred_raw + bias_eig)),
    withrich_mean_calib = pmin(1, pmax(0, wr_int + wr_slope * withrich_pred_raw + bias_wr))
  ) %>%
  filter(!is.na(eigen_23), !is.na(withrich_23))

# ==============================================================================
# 8. Monte Carlo Stage 1: Relative TPL Change
# ==============================================================================

set.seed(2025)
B <- 1000
eig_rmse_mean <- 0.056; eig_rmse_sd <- 0.007
wr_rmse_mean <- 0.181; wr_rmse_sd <- 0.004

city_rel_results <- vector("list", B)

cat("Starting Stage 1 Monte Carlo...\n")
for (b in 1:B) {
  # Sample parameters
  idx_beta <- sample(1:n_beta, 1)
  beta_b <- beta_mat[idx_beta, ]
  rmse_eig_b <- max(rnorm(1, eig_rmse_mean, eig_rmse_sd), 1e-6)
  rmse_wr_b <- max(rnorm(1, wr_rmse_mean, wr_rmse_sd), 1e-6)
  
  df_b <- future_net_base
  n_row <- nrow(df_b)
  
  # Add noise
  eps_eig <- rnorm(n_row, 0, rmse_eig_b)
  eps_wr <- rnorm(n_row, 0, rmse_wr_b)
  
  df_b <- df_b %>%
    mutate(
      eigen_sim_b = pmin(1, pmax(0, eigen_mean_calib + eps_eig)),
      withrich_sim_b = pmin(1, pmax(0, withrich_mean_calib + eps_wr)),
      raw_d_eigen = eigen_sim_b - eigen_23,
      raw_d_wr = withrich_sim_b - withrich_23,
      raw_dlog_tpl = beta_b["eigenvector"] * raw_d_eigen + beta_b["withrich"] * raw_d_wr
    )
  
  # Zero-alignment (Relative to 2023 simulation)
  base_noise_23 <- df_b %>% filter(year == 2023) %>%
    select(citycode, scenario, base_bias = raw_dlog_tpl) %>% distinct()
  
  df_b <- df_b %>%
    left_join(base_noise_23, by = c("citycode", "scenario")) %>%
    mutate(
      final_dlog_tpl = raw_dlog_tpl - replace_na(base_bias, 0),
      rel_tpl = exp(final_dlog_tpl) - 1,
      rel_pct = rel_tpl * 100
    )
  
  # Store City-level results
  city_rel_results[[b]] <- df_b %>% select(citycode, year, scenario, rel_pct) %>% mutate(iter = b)
}

city_rel_mc <- bind_rows(city_rel_results)

# Calculate Median Trajectory for Stage 2
city_median_trajectory <- city_rel_mc %>%
  group_by(citycode, year, scenario) %>%
  filter(year >= 2023 & year <= 2060) %>%
  summarise(
    median_rel_pct = median(rel_pct, na.rm = TRUE),
    .groups = "drop"
  )
write.csv(city_median_trajectory, paste0(base_path, '08 result/final_city_tpl_relativechange_slowpri.csv'), row.names = FALSE)

# ==============================================================================
# 9. Monte Carlo Stage 2: Economic Outcome
# ==============================================================================

# --- 9.1 Data Setup ---
tplmodel <- readRDS(paste0(base_path, "08 result/rf_top8_final_model.rds"))
city_rel_med <- city_median_trajectory

# TPL Calibration Params
tpl_group_params <- tibble(
  pop_group = 1:4, tpl_slope = c(1.172, 1.284, 1.311, 1.679), tpl_int = c(0.029, 0.002, -0.004, 0.006)
)
tpl_rmse_mean <- 0.831; tpl_rmse_sd <- 0.216

# Prepare Socio Data
socio_data <- result2 %>%
  mutate(
    citycode = as.character(citycode), scenario = as.character(scenario),
    gdp_day = GDP / 365,
    pop_group = cut(pop, breaks = pop_breaks, labels = FALSE, include.lowest = TRUE)
  )
pop_group_23 <- socio_data %>% filter(year == 2023) %>% select(citycode, pop_group_23 = pop_group) %>% distinct()

# Prepare Climate Data
future_climate <- future1 %>%
  mutate(citycode = as.character(citycode), scenario = as.character(scenario)) %>%
  select(year, scenario, citycode, model, total_prec, max_prec, unique_city_count, max_wind_diff, max_temp_diff)

# Build Master Data
future_master <- future_climate %>%
  inner_join(socio_data, by = c("year", "scenario", "citycode")) %>%
  left_join(pop_group_23, by = "citycode")

# Fill missing columns for TPL model
missing_vars <- setdiff(c("old65", "urbanization"), names(future_master))
if(length(missing_vars) > 0) future_master[missing_vars] <- 0

# Predict Log TPL (Baseline)
future_master$pred_log_tpl_raw <- predict(tplmodel, newdata = future_master)

# Merge Params & Network Effect
future_master <- future_master %>%
  left_join(tpl_group_params, by = "pop_group") %>%
  left_join(city_rel_med, by = c("citycode", "year", "scenario")) %>%
  mutate(rel_med = replace_na(median_rel_pct / 100, 0)) %>%
  filter(!is.na(pred_log_tpl_raw), !is.na(lambda1), !is.na(gdp_day))

# Pre-calculate GDP Sums for denominator
gdp_base_overall <- socio_data %>% left_join(pop_group_23, by="citycode") %>%
  group_by(year, scenario) %>% summarise(gdp_total = sum(GDP, na.rm=T), .groups="drop")

gdp_base_pop <- socio_data %>% left_join(pop_group_23, by="citycode") %>%
  group_by(year, scenario, pop_group_23) %>% summarise(gdp_total_group = sum(GDP, na.rm=T), .groups="drop")

# --- 9.2 Simulation ---
dt_master <- as.data.table(future_master)
climate_models <- unique(dt_master$model)
set.seed(2025)
B <- 1000
results_list <- vector("list", B)
results_list_pop <- vector("list", B)

cat("Starting Stage 2 Monte Carlo...\n")
pb <- txtProgressBar(min = 0, max = B, style = 3)

for (b in 1:B) {
  # Sample uncertainty
  selected_model <- sample(climate_models, 1)
  rmse_tpl_b <- max(rnorm(1, tpl_rmse_mean, tpl_rmse_sd), 1e-6)
  
  df_b <- copy(dt_master[model == selected_model])
  if (nrow(df_b) == 0) next
  
  n_rows <- nrow(df_b)
  eps_tpl <- rnorm(n_rows, 0, rmse_tpl_b)
  
  # Calculate Baseline TPL (Calibrated)
  df_b[, tpl_log_raw := pred_log_tpl_raw + eps_tpl]
  df_b[, tpl_pred_lev := exp(tpl_log_raw)]
  df_b[, tpl_base_val := pmax(0, pmin(365, tpl_int + tpl_slope * tpl_pred_lev))]
  
  # Sample Lambda
  l_min <- pmin(df_b$lambda1, df_b$lambda2, na.rm=TRUE); l_min[is.na(l_min)] <- 0
  l_max <- pmax(df_b$lambda1, df_b$lambda2, na.rm=TRUE); l_max[is.na(l_max)] <- 0
  lambda_vec <- runif(n_rows, l_min, l_max)
  
  # Calculate Outcome
  df_b[, tpl_delta := tpl_base_val * rel_med] # Change due to network
  df_b[, eco_outcome := tpl_delta * gdp_day * lambda_vec]
  
  # Aggregate Overall
  agg_res <- df_b[, .(total_outcome = sum(eco_outcome, na.rm=T)), by = .(year, scenario)]
  agg_res[, iter := b]
  results_list[[b]] <- agg_res
  
  # Aggregate by Pop Group
  agg_res_pop <- df_b[, .(total_outcome = sum(eco_outcome, na.rm=T)), by = .(year, scenario, pop_group_23)]
  agg_res_pop[, iter := b]
  results_list_pop[[b]] <- agg_res_pop
  
  setTxtProgressBar(pb, b)
}
close(pb)

# ==============================================================================
# 10. Final Summaries
# ==============================================================================

# Overall Summary
final_mc_df <- rbindlist(results_list, use.names = TRUE, fill = TRUE)
summary_outcome <- final_mc_df %>%
  filter(year >= 2023 & year <= 2050) %>%
  group_by(year, scenario) %>%
  summarise(
    outcome_mean = mean(total_outcome, na.rm=T),
    outcome_lower = quantile(total_outcome, 0.025, na.rm=T),
    outcome_upper = quantile(total_outcome, 0.975, na.rm=T),
    .groups = "drop"
  ) %>%
  left_join(gdp_base_overall, by = c("year", "scenario")) %>%
  mutate(outcome_pct = 100 * outcome_mean / gdp_total)

write.csv(summary_outcome, paste0(base_path, "08 result/final_all_ecoloss_raildeve_percent.csv"), row.names = FALSE)

# Pop Group Summary
final_mc_df_pop <- rbindlist(results_list_pop, use.names = TRUE, fill = TRUE)
summary_outcome_pop <- final_mc_df_pop %>%
  filter(year >= 2023 & year <= 2050) %>%
  group_by(year, scenario, pop_group_23) %>%
  summarise(
    outcome_mean = mean(total_outcome, na.rm=T),
    outcome_lower = quantile(total_outcome, 0.025, na.rm=T),
    outcome_upper = quantile(total_outcome, 0.975, na.rm=T),
    .groups = "drop"
  ) %>%
  left_join(gdp_base_pop, by = c("year", "scenario", "pop_group_23")) %>%
  mutate(outcome_pct_group = 100 * outcome_mean / gdp_total_group)

write.csv(summary_outcome_pop, paste0(base_path, "08 result/finalpopgroup_ecoloss_raildeve_percent.csv"), row.names = FALSE)

# Plotting
ggplot(summary_outcome_pop, aes(x = year, y = outcome_pct_group, color = factor(pop_group_23), fill = factor(pop_group_23))) +
  geom_ribbon(aes(ymin = 100 * outcome_lower / gdp_total_group, ymax = 100 * outcome_upper / gdp_total_group), alpha = 0.2, color = NA) +
  geom_line(size = 1) +
  facet_wrap(~ scenario) +
  theme_minimal() +
  labs(title = "Outcome as Share of Group GDP", y = "Outcome (%)", color = "Pop Group", fill = "Pop Group")