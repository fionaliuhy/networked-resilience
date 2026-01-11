rm(list=ls())

# ==============================================================================
# 1. Setup and Libraries
# ==============================================================================
packages <- c("data.table", "tidyverse", "fixest", "glmnet", 
              "randomForest", "DoubleML", "mlr3", "mlr3learners")

# Load packages
lapply(packages, library, character.only = TRUE)

options(digits = 8)  
options(scipen = 999)  

# Define output path
base_path <- "/Volumes/LaCie/08 result/updatecen/log/"

# ==============================================================================
# 2. Data Loading and Preprocessing
# ==============================================================================

# Load main dataset (Assuming column names are already the 15 English indicators)
overall <- read.csv('/Volumes/LaCie/08 result/overall_cenup1124.csv')

# Define the 15 Network Variables (Treatments)
net_vars <- c(
  "Closeness centrality", "Mobility per capita", "Pagerank", "Eigenvector centrality",
  "Mobility-weighted GDP per capita", "Near-neighbor preference", "Economic standing",
  "Economic bias", "Neighbor linkage", "Mobility-weighted distance", "Elite city linkage",
  "Intercity interaction unevenness", "Maximum connected ratio", 
  "Intercity interaction inequality", "Within province linkage"
)

# --- Transformations ---
# Log transformations for specific skewed variables
overall$`Mobility-weighted distance` <- log(overall$`Mobility-weighted distance`)
overall$`Mobility-weighted GDP per capita` <- log(overall$`Mobility-weighted GDP per capita`)

# Outcome Variable Construction (Total Performance Loss - TPL)
overall$shape <- (-1) * overall$total_intrac2 / ((overall$max_intrac * overall$recovery2) / 2)

# Log transformation of outcomes and controls
overall$logtpl <- log(overall$total_intrac)
overall$logmi <- log(overall$max_intrac)
overall$logrecovery <- log(overall$recovery)
overall$logshape <- log(overall$shape)

overall$pop <- log(overall$pop)
overall$gdppercapita <- log(overall$gdppercapita)
overall$schoolnum <- log1p(overall$schoolnum)
overall$mean_nightlight_y <- log(overall$mean_nightlight_y)
overall$roaddes <- log(overall$roaddes)
overall$roadper <- log(overall$roadper)

# Define Variables
y_var <- "logtpl"

control_vars <- c(
  "total_prec", "max_temp_diff", "max_wind_diff", "npi_diff", "max_prec",
  "mean_elevation", "urbanization", "old65", "tergdp", "pop",
  "gdppercapita", "mean_nightlight_y",
  "proportion_10", "proportion_30", "proportion_50", "proportion_60",
  "proportion_80", "proportion_90",
  "roaddes", "roadper", "schoolnum", "damdes", "ndvi", "unique_city_count", "event_day_count",
  "max_Rpolicy", "max_Cpolicy"
)

# ==============================================================================
# 3. Residualization (Removing Fixed Effects)
# ==============================================================================
# Regress variables on City and Year Fixed Effects and store residuals

df_fe <- overall %>% select(citycode, year, all_of(c(y_var, net_vars, control_vars)))
df_res <- df_fe %>% select(citycode, year)

vars_to_demean <- c(y_var, net_vars, control_vars)

for (v in vars_to_demean) {
  # Use backticks to handle spaces in variable names
  fml <- as.formula(paste0("`", v, "` ~ 1 | citycode + year"))
  m <- feols(fml, data = df_fe)
  df_res[[paste0(v, "_res")]] <- resid(m)
}

# Define residual variable names
control_res_vars <- paste0(control_vars, "_res")
y_res_var        <- paste0(y_var, "_res")
net_res_vars     <- paste0(net_vars, "_res")

# ==============================================================================
# 4. Diagnostic Analysis (Post-Residualization)
# ==============================================================================
# Check if controls (X) can predict residuals of Y or D using Elastic Net and RF.

cat("\n=== Diagnostic Analysis ===\n")

# Prepare data matrix
diag_df <- df_res %>% select(all_of(c(y_res_var, net_res_vars, control_res_vars))) %>% na.omit()
X_mat <- as.matrix(diag_df[, control_res_vars])

# Function: CV R2 for Elastic Net
get_en_r2 <- function(y, X) {
  tryCatch({
    cv_fit <- cv.glmnet(X, y, alpha = 0.5)
    preds <- predict(cv_fit, newx = X, s = "lambda.min")
    return(1 - sum((y - preds)^2) / sum((y - mean(y))^2))
  }, error = function(e) NA)
}

# Function: OOB R2 for Random Forest
get_rf_r2 <- function(y, X) {
  tryCatch({
    rf_fit <- randomForest(x = X, y = y, ntree = 100)
    return(mean(rf_fit$rsq)) # Returns vector of rsq by tree, take mean or last
  }, error = function(e) NA)
}

# Run Diagnostics
diag_results <- list()
targets <- c(y_res_var, net_res_vars)

for (target in targets) {
  y_vec <- diag_df[[target]]
  r2_en <- get_en_r2(y_vec, X_mat)
  r2_rf <- get_rf_r2(y_vec, X_mat)
  
  diag_results[[length(diag_results) + 1]] <- data.frame(
    Variable = target,
    R2_ElasticNet = r2_en,
    R2_RandomForest = r2_rf
  )
}

diag_table <- bind_rows(diag_results)
print(diag_table)
write.csv(diag_table, paste0(base_path, "diagnostic_residuals_r2.csv"), row.names = FALSE)

# ==============================================================================
# 5. Step 1: Single-Treatment DML (Elastic Net)
###for four resilience metrics， each metrics took 30min
# ==============================================================================
cat("\n=== Step 1: Single-Treatment DML (Elastic Net) ===\n")


set.seed(123)
rep_times <- 100
stab_list_single_en <- list()

for (d_name in net_vars) {
  d_res_var <- paste0(d_name, "_res")
  
  # Prepare DoubleML data object
  df_dml <- df_res %>% select(all_of(c(y_res_var, d_res_var, control_res_vars))) %>% as.data.table()
  dml_data <- DoubleMLData$new(data = df_dml, y_col = y_res_var, d_cols = d_res_var, x_cols = control_res_vars)
  
  for (r in 1:rep_times) {
    set.seed(1000 * which(net_vars == d_name) + r)
    ml_g <- lrn("regr.glmnet", alpha = 0.5, standardize = TRUE)
    ml_m <- lrn("regr.glmnet", alpha = 0.5, standardize = TRUE)
    
    dml_plr <- DoubleMLPLR$new(data = dml_data, ml_g = ml_g, ml_m = ml_m, n_folds = 10, score = "partialling out")
    dml_plr$fit()
    
    stab_list_single_en[[length(stab_list_single_en) + 1]] <- data.frame(
      D = d_name, rep = r, coef = as.numeric(dml_plr$coef), se = as.numeric(dml_plr$se), method = "Single_EN"
    )
  }
  cat("Finished Single EN:", d_name, "\n")
}

res_single_en <- bind_rows(stab_list_single_en)
write.csv(res_single_en, paste0(base_path, "step1_single_dml_en.csv"), row.names = FALSE)

# ==============================================================================
# 6. Step 2: Placebo Test (Single EN + Shuffle)
# ==============================================================================
cat("\n=== Step 2: Placebo Test (Shuffled Treatment) ===\n")

set.seed(999)
stab_list_placebo <- list()

for (d_name in net_vars) {
  d_res_var <- paste0(d_name, "_res")
  temp_data <- df_res %>% select(all_of(c(y_res_var, d_res_var, control_res_vars)))
  
  for (r in 1:rep_times) {
    # Shuffle Treatment
    current_permuted_D <- sample(temp_data[[d_res_var]])
    df_dml_placebo <- temp_data
    df_dml_placebo[[d_res_var]] <- current_permuted_D
    df_dml_placebo <- as.data.table(df_dml_placebo)
    
    dml_data_r <- DoubleMLData$new(data = df_dml_placebo, y_col = y_res_var, d_cols = d_res_var, x_cols = control_res_vars)
    set.seed(2000 * which(net_vars == d_name) + r)
    
    ml_g <- lrn("regr.glmnet", alpha = 0.5, standardize = TRUE)
    ml_m <- lrn("regr.glmnet", alpha = 0.5, standardize = TRUE)
    
    dml_plr <- DoubleMLPLR$new(data = dml_data_r, ml_g = ml_g, ml_m = ml_m, n_folds = 10, score = "partialling out")
    dml_plr$fit()
    
    stab_list_placebo[[length(stab_list_placebo) + 1]] <- data.frame(
      D = d_name, rep = r, coef = as.numeric(dml_plr$coef), se = as.numeric(dml_plr$se), method = "Placebo"
    )
  }
  cat("Finished Placebo:", d_name, "\n")
}

res_placebo <- bind_rows(stab_list_placebo)
write.csv(res_placebo, paste0(base_path, "step2_placebo_dml_en.csv"), row.names = FALSE)

# ==============================================================================
# 7. Step 3: Single-Treatment DML (Random Forest)
### each metrics took 2-3h
# ==============================================================================
cat("\n=== Step 3: Single-Treatment DML (Random Forest) ===\n")

set.seed(456)
stab_list_single_rf <- list()

for (d_name in net_vars) {
  d_res_var <- paste0(d_name, "_res")
  df_dml <- df_res %>% select(all_of(c(y_res_var, d_res_var, control_res_vars))) %>% as.data.table()
  dml_data <- DoubleMLData$new(data = df_dml, y_col = y_res_var, d_cols = d_res_var, x_cols = control_res_vars)
  
  for (r in 1:rep_times) {
    set.seed(3000 * which(net_vars == d_name) + r)
    
    # Random Forest Learners
    ml_g <- lrn("regr.ranger", num.trees = 200, min.node.size = 30, max.depth = 5)
    ml_m <- lrn("regr.ranger", num.trees = 200, min.node.size = 30, max.depth = 5)
    
    dml_plr <- DoubleMLPLR$new(data = dml_data, ml_g = ml_g, ml_m = ml_m, n_folds = 10, score = "partialling out")
    
    tryCatch({
      dml_plr$fit()
      stab_list_single_rf[[length(stab_list_single_rf) + 1]] <- data.frame(
        D = d_name, rep = r, coef = as.numeric(dml_plr$coef), se = as.numeric(dml_plr$se), method = "Single_RF"
      )
    }, error = function(e) cat("Error in RF rep", r, "\n"))
  }
  cat("Finished Single RF:", d_name, "\n")
}

res_single_rf <- bind_rows(stab_list_single_rf)
write.csv(res_single_rf, paste0(base_path, "step3_single_dml_rf.csv"), row.names = FALSE)

# ==============================================================================
# 8. Step 4: Multi-Treatment DML (Elastic Net)
###for four resilience metrics， each metrics took 30min
# ==============================================================================
cat("\n=== Step 4: Multi-Treatment DML (Elastic Net) ===\n")

df_dml_multi <- df_res %>% select(all_of(c(y_res_var, net_res_vars, control_res_vars))) %>% as.data.table()
dml_data_multi <- DoubleMLData$new(data = df_dml_multi, y_col = y_res_var, d_cols = net_res_vars, x_cols = control_res_vars)

set.seed(789)
stab_list_multi <- vector("list", rep_times)

for (r in seq_len(rep_times)) {
  set.seed(4000 + r)
  ml_l <- lrn("regr.cv_glmnet", alpha = 0.5, s = "lambda.min", standardize = TRUE)
  ml_m <- lrn("regr.cv_glmnet", alpha = 0.5, s = "lambda.min", standardize = TRUE)
  
  dml_plr <- DoubleMLPLR$new(data = dml_data_multi, ml_l = ml_l, ml_m = ml_m, n_folds = 10, score = "partialling out")
  dml_plr$fit()
  
  d_names <- gsub("_res$", "", names(dml_plr$coef))
  stab_list_multi[[r]] <- data.frame(
    D = d_names, rep = r, coef = as.numeric(dml_plr$coef), se = as.numeric(dml_plr$se), method = "Multi_EN"
  )
  if(r %% 10 == 0) cat("Multi-DML Rep:", r, "\n")
}

res_multi_en <- bind_rows(stab_list_multi)
write.csv(res_multi_en, paste0(base_path, "step4_multi_dml_en_coeffs.csv"), row.names = FALSE)

# ==============================================================================
# 9. Step 5: Scenario / Counterfactual Analysis (Raw + Winsorized)
# ==============================================================================
cat("\n=== Step 5: Scenario Analysis ===\n")

# --- 9.1 Load Only Specific Scenarios ---
scen_path <- "/Volumes/LaCie04 mobility/"

# Assuming these files already contain the 15 English Variable Names
scenario2 <- read.csv(paste0(scen_path, "network 20 new2.csv")) %>% mutate(scenario = "iso")
scenario3 <- read.csv(paste0(scen_path, "network 21 new2.csv")) %>% mutate(scenario = "weak")
scenario5 <- read.csv(paste0(scen_path, "network 23 new2.csv")) %>% mutate(scenario = "s23") # Baseline
scenario6 <- read.csv(paste0(scen_path, "network 24 new2.csv")) %>% mutate(scenario = "s24")
scenario9 <- read.csv(paste0(scen_path, "network 13 new4 filled.csv")) %>% mutate(scenario = "eco2")

# Combine list
scen_list <- list(scenario2, scenario3, scenario5, scenario6, scenario9)
names(scen_list) <- c("iso", "weak", "s23", "s24", "eco2")
baseline_net <- scen_list[["s23"]]

# --- 9.2 Prepare Data ---
# Baseline TPL (2023)
baseline_tpl_df <- overall %>%
  filter(year == 2023) %>%
  mutate(citycode = as.character(citycode), tpl_base_abs = total_intrac) %>%
  select(citycode, tpl_base_abs)

# Coefficients Matrix
beta_by_rep <- res_multi_en %>% select(D, rep, coef) %>% pivot_wider(names_from = D, values_from = coef)

# --- 9.3 Calculation Function (Raw + Winsorized) ---
run_scenario_calc <- function(scen_name, scen_df, baseline_net, baseline_tpl, beta_by_rep, net_vars, winsorize = FALSE) {
  
  df_proc <- scen_df %>%
    mutate(citycode = as.character(citycode)) %>%
    select(citycode, all_of(net_vars)) %>%
    inner_join(baseline_net %>% mutate(citycode = as.character(citycode)) %>% select(citycode, all_of(net_vars)), 
               by = "citycode", suffix = c("_scen", "_base")) %>%
    inner_join(baseline_tpl, by = "citycode")
  
  if (nrow(df_proc) == 0) return(NULL)
  
  # Winsorization (Optional)
  if (winsorize) {
    limits <- baseline_net %>%
      summarise(across(all_of(net_vars), list(min = ~min(., na.rm=T), max = ~max(., na.rm=T))))
    
    for (v in net_vars) {
      v_min <- limits[[paste0(v, "_min")]]
      v_max <- limits[[paste0(v, "_max")]]
      df_proc[[paste0(v, "_scen")]] <- pmax(v_min, pmin(df_proc[[paste0(v, "_scen")]], v_max))
    }
  }
  
  # Calculate Delta D
  for (v in net_vars) {
    df_proc[[paste0("d_", v)]] <- df_proc[[paste0(v, "_scen")]] - df_proc[[paste0(v, "_base")]]
  }
  
  # Monte Carlo Simulation
  beta_mat <- as.matrix(beta_by_rep[, net_vars])
  n_reps <- nrow(beta_mat)
  results <- list()
  total_tpl_baseline <- sum(df_proc$tpl_base_abs, na.rm = TRUE)
  
  for (i in 1:n_reps) {
    beta_vec <- beta_mat[i, ]
    delta_log_tpl <- numeric(nrow(df_proc))
    
    for(v in net_vars) {
      delta_log_tpl <- delta_log_tpl + (beta_vec[v] * df_proc[[paste0("d_", v)]])
    }
    
    city_rel_change <- exp(delta_log_tpl) - 1
    city_abs_delta <- df_proc$tpl_base_abs * city_rel_change
    
    results[[i]] <- data.frame(
      scenario = scen_name,
      type = ifelse(winsorize, "Winsorized", "Raw"),
      rep = i,
      mean_rel_change = mean(city_rel_change, na.rm = TRUE),
      total_delta_tpl = sum(city_abs_delta, na.rm = TRUE),
      total_ratio_change = sum(city_abs_delta, na.rm = TRUE) / total_tpl_baseline
    )
  }
  bind_rows(results)
}

# --- 9.4 Execute Both Analyses ---
res_raw <- bind_rows(lapply(names(scen_list), function(nm) {
  run_scenario_calc(nm, scen_list[[nm]], baseline_net, baseline_tpl_df, beta_by_rep, net_vars, winsorize = FALSE)
}))

res_win <- bind_rows(lapply(names(scen_list), function(nm) {
  run_scenario_calc(nm, scen_list[[nm]], baseline_net, baseline_tpl_df, beta_by_rep, net_vars, winsorize = TRUE)
}))

all_scen_results <- bind_rows(res_raw, res_win)

# --- 9.5 Summarize ---
final_summary <- all_scen_results %>%
  pivot_longer(cols = c("mean_rel_change", "total_delta_tpl", "total_ratio_change"), names_to = "metric") %>%
  group_by(scenario, type, metric) %>%
  summarise(
    est = mean(value),
    ci_low = quantile(value, 0.025),
    ci_high = quantile(value, 0.975),
    .groups = "drop"
  )

write.csv(final_summary, paste0(base_path, "step5_scenario_summary_final.csv"), row.names = FALSE)
print(final_summary)
