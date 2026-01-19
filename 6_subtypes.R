# ============================================================================
# TYPHOON IMPACT ANALYSIS - FIXED PARAMETERS VERSION
# COMPLETE ANALYSIS: FIGURES 1 TO 16 (NO OMISSIONS)
# ============================================================================
# Sys.setenv('R_MAX_VSIZE' = 64 * 1024^3)

library(rstan)
library(ggplot2)
library(dplyr)
library(lubridate)
library(bayesplot)
library(tidyr)
library(reshape2)
library(patchwork)
library(grid)
library(gridExtra)
library(gtable)
library(knitr)
library(scales)

options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

# ============================================================================
# 1. DATA LOADING AND PREPROCESSING
# ============================================================================
cat("=== 1. Loading and preprocessing data ===\n")

# 请确保路径正确
data_path <- "/Users/chenjiaqi/Desktop/COVID-19_HK/typhoon/HK_ILI_COVID_Sep.csv"
raw_data <- read.csv(data_path)
raw_data$date <- as.Date(raw_data$date, format = "%Y/%m/%d")

#start_date <- as.Date("2025-09-13")
start_date <- as.Date("2025-06-01")
#start_date <- as.Date("2025-08-02")
#start_date <- as.Date("2022-12-31")
typhoon_start_date <- as.Date("2025-09-27")
#typhoon_start_date <- as.Date("2025-10-11")
#typhoon_start_date <- as.Date("2025-12-06")
end_date <- as.Date("2025-12-13")


fitting_data <- raw_data[raw_data$date >= start_date & raw_data$date < typhoon_start_date, ]
validation_data <- raw_data[raw_data$date >= typhoon_start_date & raw_data$date <= end_date, ]

fitting_data <- fitting_data[complete.cases(fitting_data), ]
validation_data <- validation_data[complete.cases(validation_data), ]

strain_names <- c("B", "H3", "H1", "COVID", "RSV", "HFMD")
fitting_data <- fitting_data[, c("date", strain_names)]
validation_data <- validation_data[, c("date", strain_names)]

cat("Fitting period:", as.character(range(fitting_data$date)), "\n")
cat("Validation period:", as.character(range(validation_data$date)), "\n")
cat("Total fitting weeks:", nrow(fitting_data), "\n")
cat("Strains:", paste(strain_names, collapse=", "), "\n\n")

# ============================================================================
# 2. BASIC PARAMETERS & ATTACK RATE PERIODS
# ============================================================================
T_weeks <- nrow(fitting_data)
T_weeks_validation <- nrow(validation_data)
T_weeks_forecast <- 25
typhoon_weeks <- 1
N_strains <- 6
T_weeks_total <- T_weeks + T_weeks_forecast

cat("=== Basic Parameters ===\n")
cat("T_weeks (historical):", T_weeks, "\n")
cat("T_weeks_forecast:", T_weeks_forecast, "\n\n")

# Define Attack Rate Periods
cat("=== Defining Attack Rate Periods ===\n")
period1_start_date <- as.Date("2025-11-29")
#period1_start_date <- as.Date("2025-09-20")
#period1_end_date <- as.Date("2025-09-27")
period1_end_date <- as.Date("2025-12-06") 
#period2_additional_weeks <- 3 
period2_additional_weeks <- 0 
period1_start_week <- which(fitting_data$date == period1_start_date)

if (period1_end_date <= max(fitting_data$date)) {
  period1_end_week <- which(fitting_data$date == period1_end_date)
} else {
  period1_end_week <- as.numeric(difftime(period1_end_date, min(fitting_data$date), units = "weeks")) + 1
}

period2_end_week <- period1_end_week + period2_additional_weeks

cat("Period 1: Week", period1_start_week, "to", period1_end_week, "\n")
cat("Period 2: Week", period1_start_week, "to", period2_end_week, "\n\n")

# ============================================================================
# 3. STAN DATA PREPARATION (WITH FIXED PARAMETERS)
# ============================================================================
cat("=== Setting up Fixed Parameters based on Literature ===\n")

# 1. Define Fixed Parameters (Rates per day) [Source: Literature]
# Order: B, H3, H1, COVID, RSV, HFMD

# Sigma = 1 / Incubation Period
# Flu (2d), COVID (3.4d), RSV (4.4d), HFMD (4.5d)

# Prepare Matrices
cases_matrix <- as.matrix(fitting_data[, strain_names])
cases_matrix[cases_matrix < 0] <- 0
cases_matrix <- round(cases_matrix)

validation_matrix <- as.matrix(validation_data[, strain_names])
validation_matrix[validation_matrix < 0] <- NA
validation_matrix <- round(validation_matrix)

# 2. Update Stan Data List
stan_data <- list(
  T_weeks = T_weeks,
  T_weeks_forecast = T_weeks_forecast,
  N_strains = N_strains,
  cases = cases_matrix,
  num_knots = 8,
  spline_degree = 3,
  population = 7524000,
  typhoon_weeks = typhoon_weeks,
  child_ratio = 0.029,
  #child_ratio = 29,
  
  # Attack rate periods
  period1_start_week = period1_start_week,
  period1_end_week = period1_end_week,
  period2_additional_weeks = period2_additional_weeks
  
  # NEW: Pass Fixed Parameters to Stan
  #fixed_sigma = fixed_sigma,
  #fixed_gamma = fixed_gamma,
  #fixed_mu = fixed_mu
)

cat("Stan data prepared with fixed parameter vectors.\n\n")

# ============================================================================
# 4. MODEL COMPILATION AND FITTING
# ============================================================================
cat("=== Compiling Stan model ===\n")
# 指向修改后的 .stan 文件 (fixed version)
model_file <- "/Users/chenjiaqi/SPH Dropbox/Jackie Chen/Shared RESV_HK_Typhoon/Fig_rds/6_subtypes_last_epi.stan"
model <- stan_model(model_file)

cat("\n=== Starting MCMC sampling ===\n")
# Uncomment to run new fit:


fit <- sampling(
  model,
  data = stan_data,
  iter = 2000,
  warmup = 1000,
  chains = 4,
  thin = 1,
  cores = 4,
  control = list(adapt_delta = 0.95, max_treedepth = 15)
)


saveRDS(fit, file = "/Users/chenjiaqi/SPH Dropbox/Jackie Chen/rds/fitting_result/initial_value_09_13.rds")





#fit <- readRDS("/Users/chenjiaqi/SPH Dropbox/Jackie Chen/rds/fitting_result/initial_value_08_02.rds")



# ============================================================================
# EXTRACT RESULTS FROM SEIR MODEL (NO CROSS-IMMUNITY)
# ============================================================================
cat("\n=== Extracting results ===\n")

pred_cases <- rstan::extract(fit, pars = "pred_cases")$pred_cases
forecast_cases <- rstan::extract(fit, pars = "forecast_cases")$forecast_cases
forecast_cases_typhoon <- rstan::extract(fit, pars = "forecast_cases_typhoon")$forecast_cases_typhoon
R_eff <- rstan::extract(fit, pars = "R_eff")$R_eff
R_eff_scenarios <- rstan::extract(fit, pars = "R_eff_scenarios")$R_eff_scenarios
R0_t <- rstan::extract(fit, pars = "R0_t")$R0_t
reduction_typhoon <- rstan::extract(fit, pars = "reduction_typhoon_period")$reduction_typhoon_period
avg_weekly_typhoon <- rstan::extract(fit, pars = "avg_weekly_typhoon")$avg_weekly_typhoon
avg_weekly_recovery <- rstan::extract(fit, pars = "avg_weekly_recovery")$avg_weekly_recovery

# Extended scenarios
cases_extended <- rstan::extract(fit, pars = "cases_extended_typhoon")$cases_extended_typhoon
incidence_per10k_extended <- rstan::extract(fit, pars = "incidence_per10k_extended")$incidence_per10k_extended
reduction_extended <- rstan::extract(fit, pars = "reduction_extended_typhoon")$reduction_extended_typhoon

# Attack rates
attack_rate_baseline <- rstan::extract(fit, pars = "attack_rate_baseline")$attack_rate_baseline
attack_rate_scenarios <- rstan::extract(fit, pars = "attack_rate_scenarios")$attack_rate_scenarios
attack_rate_diff <- rstan::extract(fit, pars = "attack_rate_diff")$attack_rate_diff

# Period-specific attack rates
attack_rate_period1_baseline <- rstan::extract(fit, pars = "attack_rate_period1_baseline")$attack_rate_period1_baseline
attack_rate_period1_scenarios <- rstan::extract(fit, pars = "attack_rate_period1_scenarios")$attack_rate_period1_scenarios
attack_rate_period1_extended <- rstan::extract(fit, pars = "attack_rate_period1_extended")$attack_rate_period1_extended

attack_rate_period2_baseline <- rstan::extract(fit, pars = "attack_rate_period2_baseline")$attack_rate_period2_baseline
attack_rate_period2_scenarios <- rstan::extract(fit, pars = "attack_rate_period2_scenarios")$attack_rate_period2_scenarios
attack_rate_period2_extended <- rstan::extract(fit, pars = "attack_rate_period2_extended")$attack_rate_period2_extended

cat("\n=== Checking extracted dimensions ===\n")
cat("pred_cases dimensions:", dim(pred_cases), "\n")
cat("forecast_cases dimensions:", dim(forecast_cases), "\n")
cat("forecast_cases_typhoon dimensions:", dim(forecast_cases_typhoon), "\n")
cat("R_eff dimensions:", dim(R_eff), "\n")
cat("R_eff_scenarios dimensions:", dim(R_eff_scenarios), "\n")
cat("R0_t dimensions:", dim(R0_t), "\n")
cat("reduction_typhoon dimensions:", dim(reduction_typhoon), "\n")
cat("cases_extended dimensions:", dim(cases_extended), "\n")
cat("reduction_extended dimensions:", dim(reduction_extended), "\n")
cat("attack_rate_baseline dimensions:", dim(attack_rate_baseline), "\n")
cat("attack_rate_scenarios dimensions:", dim(attack_rate_scenarios), "\n")
cat("attack_rate_diff dimensions:", dim(attack_rate_diff), "\n")

# Verify attack rate period dimensions
cat("\n=== Attack Rate Period Dimensions ===\n")
cat("attack_rate_period1_baseline dimensions:", dim(attack_rate_period1_baseline), "\n")
cat("attack_rate_period1_scenarios dimensions:", dim(attack_rate_period1_scenarios), "\n")
cat("attack_rate_period1_extended dimensions:", dim(attack_rate_period1_extended), "\n")
cat("attack_rate_period2_baseline dimensions:", dim(attack_rate_period2_baseline), "\n")
cat("attack_rate_period2_scenarios dimensions:", dim(attack_rate_period2_scenarios), "\n")
cat("attack_rate_period2_extended dimensions:", dim(attack_rate_period2_extended), "\n")

# Check if Period 1 extends into forecast
if (period1_end_week > T_weeks) {
  cat("\n✓ Period 1 extends into forecast period\n")
  cat("  Historical component: weeks", period1_start_week, "to", T_weeks, "\n")
  cat("  Forecast component: weeks", T_weeks + 1, "to", period1_end_week, "\n")
  cat("  → Attack rates will VARY by scenario in forecast portion\n")
} else {
  cat("\n✓ Period 1 entirely within historical data\n")
  cat("  → Attack rates will be IDENTICAL across scenarios\n")
}

if (period2_end_week > T_weeks) {
  cat("\n✓ Period 2 extends into forecast period\n")
  cat("  → Attack rates will VARY by scenario\n")
} else {
  cat("\n✓ Period 2 entirely within historical data\n")
}

# ============================================================================
# DEFINE SCENARIO LABELS
# ============================================================================
typhoon_scenarios <- c(
  "0% (Baseline)",
  # 3 days
  "3d -50%", "3d -30%", "3d -20%", "3d -10%", "3d -5%",
  "3d +5%", "3d +10%", "3d +20%", "3d +30%", "3d +50%",
  # 7 days
  "7d -50%", "7d -30%", "7d -20%", "7d -10%", "7d -5%",
  "7d +5%", "7d +10%", "7d +20%", "7d +30%", "7d +50%",
  # 10 days
  "10d -50%", "10d -30%", "10d -20%", "10d -10%", "10d -5%",
  "10d +5%", "10d +10%", "10d +20%", "10d +30%", "10d +50%",
  # 14 days
  "14d -50%", "14d -30%", "14d -20%", "14d -10%", "14d -5%",
  "14d +5%", "14d +10%", "14d +20%", "14d +30%", "14d +50%",
  # Shifted timing (7 days, -50%)
  "7d -50% E7", "7d -50% E5", "7d -50% E3",
  "7d -50% D3", "7d -50% D5", "7d -50% D7"
)

# ============================================================================
# CALCULATE STATISTICS
# ============================================================================
cat("\n=== Calculating statistics ===\n")

pred_mean <- apply(pred_cases, c(2,3), mean)
pred_median <- apply(pred_cases, c(2,3), median)
pred_lower <- apply(pred_cases, c(2,3), quantile, probs = 0.025)
pred_upper <- apply(pred_cases, c(2,3), quantile, probs = 0.975)

forecast_mean <- apply(forecast_cases, c(2,3), mean)
forecast_median <- apply(forecast_cases, c(2,3), median)
forecast_lower <- apply(forecast_cases, c(2,3), quantile, probs = 0.025)
forecast_upper <- apply(forecast_cases, c(2,3), quantile, probs = 0.975)

# Typhoon forecast summaries (4D: scenarios × weeks × strains)
forecast_typhoon_mean <- apply(forecast_cases_typhoon, c(2,3,4), mean)
forecast_typhoon_median <- apply(forecast_cases_typhoon, c(2,3,4), median)
forecast_typhoon_lower <- apply(forecast_cases_typhoon, c(2,3,4), quantile, probs = 0.025)
forecast_typhoon_upper <- apply(forecast_cases_typhoon, c(2,3,4), quantile, probs = 0.975)

# R_eff scenarios summaries (4D: scenarios × weeks × strains)
R_eff_scenarios_mean <- apply(R_eff_scenarios, c(2,3,4), mean)
R_eff_scenarios_lower <- apply(R_eff_scenarios, c(2,3,4), quantile, probs = 0.025)
R_eff_scenarios_upper <- apply(R_eff_scenarios, c(2,3,4), quantile, probs = 0.975)

# Original attack rates
attack_rate_baseline_mean <- apply(attack_rate_baseline, 2, mean)
attack_rate_baseline_lower <- apply(attack_rate_baseline, 2, quantile, probs = 0.025)
attack_rate_baseline_upper <- apply(attack_rate_baseline, 2, quantile, probs = 0.975)

attack_rate_scenarios_mean <- apply(attack_rate_scenarios, c(2,3), mean)
attack_rate_scenarios_lower <- apply(attack_rate_scenarios, c(2,3), quantile, probs = 0.025)
attack_rate_scenarios_upper <- apply(attack_rate_scenarios, c(2,3), quantile, probs = 0.975)

attack_rate_diff_mean <- apply(attack_rate_diff, c(2,3), mean)
attack_rate_diff_lower <- apply(attack_rate_diff, c(2,3), quantile, probs = 0.025)
attack_rate_diff_upper <- apply(attack_rate_diff, c(2,3), quantile, probs = 0.975)

# Period-specific attack rate statistics
# Period 1
attack_rate_period1_baseline_mean <- apply(attack_rate_period1_baseline, 2, mean)
attack_rate_period1_baseline_lower <- apply(attack_rate_period1_baseline, 2, quantile, probs = 0.025)
attack_rate_period1_baseline_upper <- apply(attack_rate_period1_baseline, 2, quantile, probs = 0.975)

attack_rate_period1_scenarios_mean <- apply(attack_rate_period1_scenarios, c(2,3), mean)
attack_rate_period1_scenarios_lower <- apply(attack_rate_period1_scenarios, c(2,3), quantile, probs = 0.025)
attack_rate_period1_scenarios_upper <- apply(attack_rate_period1_scenarios, c(2,3), quantile, probs = 0.975)

attack_rate_period1_extended_mean <- apply(attack_rate_period1_extended, c(2,3), mean)
attack_rate_period1_extended_lower <- apply(attack_rate_period1_extended, c(2,3), quantile, probs = 0.025)
attack_rate_period1_extended_upper <- apply(attack_rate_period1_extended, c(2,3), quantile, probs = 0.975)

# Period 2
attack_rate_period2_baseline_mean <- apply(attack_rate_period2_baseline, 2, mean)
attack_rate_period2_baseline_lower <- apply(attack_rate_period2_baseline, 2, quantile, probs = 0.025)
attack_rate_period2_baseline_upper <- apply(attack_rate_period2_baseline, 2, quantile, probs = 0.975)

attack_rate_period2_scenarios_mean <- apply(attack_rate_period2_scenarios, c(2,3), mean)
attack_rate_period2_scenarios_lower <- apply(attack_rate_period2_scenarios, c(2,3), quantile, probs = 0.025)
attack_rate_period2_scenarios_upper <- apply(attack_rate_period2_scenarios, c(2,3), quantile, probs = 0.975)

attack_rate_period2_extended_mean <- apply(attack_rate_period2_extended, c(2,3), mean)
attack_rate_period2_extended_lower <- apply(attack_rate_period2_extended, c(2,3), quantile, probs = 0.025)
attack_rate_period2_extended_upper <- apply(attack_rate_period2_extended, c(2,3), quantile, probs = 0.975)

cat("Attack rate statistics calculated successfully\n\n")

# Extended scenarios mean
incidence_per10k_mean <- apply(incidence_per10k_extended, c(2,3), mean)
reduction_extended_mean <- apply(reduction_extended, c(2,3), mean)

cat("Statistics calculated successfully\n\n")

# ============================================================================
# PREPARE VISUALIZATION DATA
# ============================================================================
library(ggplot2)
library(dplyr)
library(tidyr)

theme_set(theme_minimal(base_size = 11))

forecast_dates <- seq(typhoon_start_date, by = "week", length.out = T_weeks_forecast)
all_dates <- c(fitting_data$date, forecast_dates)

# Date markers for visualization
typhoon_start_date_line <- max(fitting_data$date)
typhoon_end_date_line <- max(fitting_data$date) + 7
typhoon_period_end <- max(fitting_data$date) + 7
recovery_start_date <- max(fitting_data$date) + 7

cat("\n=== Date markers ===\n")
cat("Last fitting date:", as.character(max(fitting_data$date)), "\n")
cat("First forecast date:", as.character(min(forecast_dates)), "\n")
cat("Typhoon start date (data split):", as.character(typhoon_start_date), "\n")
cat("Typhoon start line (visualization):", as.character(typhoon_start_date_line), "\n")
cat("Typhoon end (red line):", as.character(typhoon_end_date_line), "\n")
cat("Recovery start:", as.character(recovery_start_date), "\n\n")

# Prepare data frames
results_df <- data.frame(
  week = rep(1:T_weeks, N_strains),
  date = rep(fitting_data$date, N_strains),
  strain = factor(rep(strain_names, each = T_weeks), levels = strain_names),
  observed = as.vector(cases_matrix),
  predicted = as.vector(pred_median),
  pred_mean = as.vector(pred_mean),
  pred_lower = as.vector(pred_lower),
  pred_upper = as.vector(pred_upper)
)

last_fit_point <- results_df %>%
  group_by(strain) %>%
  slice_tail(n = 1) %>%
  ungroup()

forecast_baseline_raw <- data.frame(
  week = rep((T_weeks + 1):(T_weeks + T_weeks_forecast), N_strains),
  date = rep(forecast_dates, N_strains),
  strain = factor(rep(strain_names, each = T_weeks_forecast), levels = strain_names),
  observed = NA,
  predicted = as.vector(forecast_median),
  pred_mean = as.vector(forecast_mean),
  pred_lower = as.vector(forecast_lower),
  pred_upper = as.vector(forecast_upper)
)

forecast_baseline <- rbind(
  last_fit_point %>% select(week, date, strain, observed, predicted, pred_mean, pred_lower, pred_upper),
  forecast_baseline_raw
)

# Prepare typhoon scenario forecasts
typhoon_forecast_list <- list()

for (typhoon_idx in 1:47) {
  for (strain_idx in 1:N_strains) {
    last_point <- last_fit_point %>% filter(strain == strain_names[strain_idx])
    
    scenario_forecast <- data.frame(
      week = (T_weeks + 1):(T_weeks + T_weeks_forecast),
      date = forecast_dates,
      strain = strain_names[strain_idx],
      observed = NA,
      predicted = forecast_typhoon_median[typhoon_idx, , strain_idx],
      pred_mean = forecast_typhoon_mean[typhoon_idx, , strain_idx],
      pred_lower = forecast_typhoon_lower[typhoon_idx, , strain_idx],
      pred_upper = forecast_typhoon_upper[typhoon_idx, , strain_idx],
      scenario = typhoon_scenarios[typhoon_idx]
    )
    
    last_point$scenario <- typhoon_scenarios[typhoon_idx]
    scenario_with_overlap <- rbind(
      last_point %>% select(week, date, strain, observed, predicted, pred_mean, pred_lower, pred_upper, scenario),
      scenario_forecast
    )
    
    typhoon_forecast_list[[length(typhoon_forecast_list) + 1]] <- scenario_with_overlap
  }
}

typhoon_forecast_df <- do.call(rbind, typhoon_forecast_list)
typhoon_forecast_df$strain <- factor(typhoon_forecast_df$strain, levels = strain_names)

# Validation data (if available)
validation_plot_df <- data.frame(
  date = rep(validation_data$date, N_strains),
  strain = factor(rep(strain_names, each = T_weeks_validation), levels = strain_names),
  observed = as.vector(validation_matrix)
)
validation_plot_df <- validation_plot_df[!is.na(validation_plot_df$observed), ]

cat("Visualization data prepared\n\n")

# ============================================================================
# FIGURE 1: Historical Fit + Baseline Forecast + Validation Data
# ============================================================================
cat("\n=== Creating Figure 1: Fit, Baseline Forecast, and Validation ===\n")

fit_and_forecast_df <- rbind(
  cbind(results_df, period = "Fitted"),
  cbind(forecast_baseline, period = "Forecast")
)

p1_fit_forecast <- ggplot(fit_and_forecast_df, aes(x = date)) +
  annotate("rect",
           xmin = typhoon_start_date_line,
           xmax = typhoon_end_date_line,
           ymin = -Inf, ymax = Inf,
           fill = "pink", alpha = 0.15) +
  annotate("rect",
           xmin = recovery_start_date,
           xmax = max(forecast_dates),
           ymin = -Inf, ymax = Inf,
           fill = "lightgreen", alpha = 0.1) +
  geom_ribbon(data = filter(fit_and_forecast_df, period == "Fitted"),
              aes(ymin = pred_lower, ymax = pred_upper),
              alpha = 0.25, fill = "#377EB8") +
  geom_line(data = filter(fit_and_forecast_df, period == "Fitted"),
            aes(y = predicted), color = "#377EB8", size = 1.3) +
  geom_ribbon(data = filter(fit_and_forecast_df, period == "Forecast"),
              aes(ymin = pred_lower, ymax = pred_upper),
              alpha = 0.25, fill = "#E41A1C") +
  geom_line(data = filter(fit_and_forecast_df, period == "Forecast"),
            aes(y = predicted), color = "#E41A1C", size = 1.3) +
  geom_point(data = filter(fit_and_forecast_df, !is.na(observed), period == "Fitted"),
             aes(y = observed), color = "black", size = 1.2, alpha = 0.6) +
  geom_point(data = validation_plot_df,
             aes(x = date, y = observed),
             color = "gold", fill = "yellow", shape = 23, size = 1.2, stroke = 0.8) +
  geom_vline(xintercept = typhoon_start_date_line,
             linetype = "dashed", color = "gray40", size = 0.8) +
  geom_vline(xintercept = typhoon_end_date_line,
             linetype = "dashed", color = "red", size = 0.8) +
  facet_wrap(~strain, scales = "free_y", ncol = 2, nrow = 3) +
  scale_y_continuous(trans = "sqrt") +
  scale_x_date(date_labels = "%Y-%m", date_breaks = "4 months") +
  labs(
    title = "SEIR Model: Historical Fit and Baseline Forecast (No Cross-Immunity)",
    subtitle = "Black dots: Observed (fitted) | Yellow diamonds: Validation data | Blue: Model fit | Red: Baseline forecast\nGray dashed: Typhoon start | Red dashed: Typhoon end | Each strain evolves independently",
    x = NULL,
    y = "Weekly Cases/Hospitalizations (√ scale)"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "gray30", fill = NA, size = 0.5),
    panel.background = element_rect(fill = "white", color = NA),
    strip.text = element_text(size = 11, face = "bold"),
    strip.background = element_rect(fill = "gray95", color = "gray30", size = 0.5),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
    plot.title = element_text(size = 13, face = "bold"),
    plot.subtitle = element_text(size = 8.5, color = "gray30")
  )

print(p1_fit_forecast)

# ============================================================================
# ADDITIONAL NOTE ON MODEL CHANGES
# ============================================================================
cat("\n=== MODEL MODIFICATIONS SUMMARY ===\n")
cat("✓ SEIRS → SEIR: Removed immunity waning (mu parameters)\n")
cat("✓ Cross-immunity: Removed all cross-protection between strains\n")
cat("✓ Independence: Each strain now evolves independently\n")
cat("✓ State transitions: S → E → I → R (no R → S)\n")
cat("✓ R_eff calculation: Simplified to R0 × S for each strain\n\n")


















# ============================================================================
# 提取9月20日的β值和SEIR状态值（含95% CI）
# ============================================================================

# 确定目标日期对应的天数索引
target_date <- as.Date("2025-09-27")
day_index <- as.numeric(difftime(target_date, start_date, units = "days")) + 1

cat("目标日期:", as.character(target_date), "\n")
cat("对应天数索引:", day_index, "\n\n")

# 提取transmission_rate (β = transmission_rate)
transmission_rate <- rstan::extract(fit, pars = "transmission_rate")$transmission_rate

# 提取states
states <- rstan::extract(fit, pars = "states")$states

# 计算β的统计量
beta_results <- data.frame(Strain = strain_names)
beta_results$Beta <- sapply(1:N_strains, function(i) {
  sprintf("%.4f (%.4f – %.4f)",
          mean(transmission_rate[, i, day_index]),
          quantile(transmission_rate[, i, day_index], probs = 0.025),
          quantile(transmission_rate[, i, day_index], probs = 0.975))
})

# 提取SEIR状态
seir_results <- data.frame(Strain = strain_names)

for (i in 1:N_strains) {
  s_idx <- 1
  e_idx <- 2 + 3*(i-1)
  i_idx <- 3 + 3*(i-1)
  r_idx <- 4 + 3*(i-1)
  
  seir_results$S[i] <- sprintf("%.6f (%.6f – %.6f)",
                               mean(states[, day_index, s_idx]),
                               quantile(states[, day_index, s_idx], probs = 0.025),
                               quantile(states[, day_index, s_idx], probs = 0.975))
  
  seir_results$E[i] <- sprintf("%.6f (%.6f – %.6f)",
                               mean(states[, day_index, e_idx]),
                               quantile(states[, day_index, e_idx], probs = 0.025),
                               quantile(states[, day_index, e_idx], probs = 0.975))
  
  seir_results$I[i] <- sprintf("%.6f (%.6f – %.6f)",
                               mean(states[, day_index, i_idx]),
                               quantile(states[, day_index, i_idx], probs = 0.025),
                               quantile(states[, day_index, i_idx], probs = 0.975))
  
  seir_results$R[i] <- sprintf("%.6f (%.6f – %.6f)",
                               mean(states[, day_index, r_idx]),
                               quantile(states[, day_index, r_idx], probs = 0.025),
                               quantile(states[, day_index, r_idx], probs = 0.975))
}

# 输出结果
cat("=== β值 (2025-09-27) ===\n")
print(beta_results, row.names = FALSE)

cat("\n=== SEIR状态值 (2025-09-27) ===\n")
print(seir_results, row.names = FALSE)