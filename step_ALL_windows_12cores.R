################################################################################
# Step Count Feature Engineering - ALL Windows (12 Cores - FULLY PARALLEL)
# Filename: step_ALL_windows_12cores_CORRECTED.R
# Purpose: Calculate dynamic features from step count data across multiple time windows
#          Time windows based on EMA data (ema_day) for each mlife_id
#          Windows END 4 hours before current time
# Time Windows: 
#   - 1 hour: (5 hours before to 4 hours before)
#   - 4 hours: (8 hours before to 4 hours before)
#   - 1 day: (1 day + 4 hours before to 4 hours before)
#   - 7 days: (7 days + 4 hours before to 4 hours before)
# Metrics: Mean, MAD, P90-P10, Entropy
# Parallel Processing: ALL sections use 12 cores
# Author: Research Team
# Date: 2026
################################################################################

# =============================================================================
# SECTION 1: SETUP AND DEPENDENCIES
# =============================================================================

library(tidyverse)
library(data.table)
library(slider)
library(lubridate)
library(parallel)
library(future.apply)
library(progressr)

# Set number of cores for parallel processing
n_cores <- max(1, detectCores() - 1)
cat("Using", n_cores, "cores for parallel processing\n")

# =============================================================================
# SECTION 2: HELPER FUNCTIONS
# =============================================================================

## 2.1 Set Entropy (Shannon Entropy) ----
set_entropy <- function(x, min_n = 2) {
  # Remove NA values
  x <- x[!is.na(x)]
  
  # Check if there are enough non-NA values
  if (length(x) < min_n) {
    return(NA_real_)
  }
  
  # Calculate probability distribution
  p <- table(x) / length(x)
  
  # Calculate Shannon entropy in bits (log2)
  -sum(p * log2(p))
}

## 2.2 P90-P10 Range ----
p90_p10 <- function(x) {
  # Remove NA values
  x <- x[!is.na(x)]
  
  # Check if there are enough data points
  if (length(x) < 3) {
    return(NA_real_)
  }
  
  # Calculate quantiles
  p10 <- quantile(x, 0.1, names = FALSE, type = 7, na.rm = TRUE)
  p90 <- quantile(x, 0.9, names = FALSE, type = 7, na.rm = TRUE)
  
  # Check for NA values in calculated quantiles
  if (is.na(p10) || is.na(p90)) {
    return(NA_real_)
  }
  
  # Check if p90 equals p10 (no variation)
  if (p90 == p10) {
    return(0)
  }
  
  # Return the range
  as.numeric(p90 - p10)
}

## 2.3 Median Absolute Deviation (MAD) ----
calculate_mad <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) < 2) return(NA_real_)
  mad(x, na.rm = TRUE)
}

## 2.4 Mean ----
calculate_mean <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) < 1) return(NA_real_)
  mean(x, na.rm = TRUE)
}

# =============================================================================
# SECTION 3: LOAD DATA
# =============================================================================

cat("\n=== Loading step count data ===\n")

# Load the cleaned step count data
steps <- fread("steps_clean.csv")

cat("Loaded step data:\n")
cat("  Rows:", nrow(steps), "\n")
cat("  Unique UIDs:", length(unique(steps$uid)), "\n")
cat("  Date range:", as.character(min(steps$time_bin)), "to", 
    as.character(max(steps$time_bin)), "\n")

# Ensure time_bin is POSIXct
steps <- steps %>%
  mutate(time_bin = as.POSIXct(time_bin, tz = "UTC")) %>%
  arrange(uid, time_bin)

cat("\n=== Loading EMA data ===\n")

# Load EMA features to get time windows per participant
ema <- fread("ema_features.csv")

cat("Loaded EMA data:\n")
cat("  Rows:", nrow(ema), "\n")
cat("  Unique mlife_ids:", length(unique(ema$mlife_id)), "\n")

# Check if ema_day column exists
if (!"ema_day" %in% names(ema)) {
  stop("ERROR: ema_day column not found in ema_features.csv")
}

# Extract unique mlife_id and their ema_day ranges
ema_windows <- ema %>%
  select(mlife_id, ema_day) %>%
  group_by(mlife_id) %>%
  summarise(
    min_ema_day = min(ema_day, na.rm = TRUE),
    max_ema_day = max(ema_day, na.rm = TRUE),
    n_ema_days = n_distinct(ema_day),
    .groups = 'drop'
  )

cat("\nEMA window summary:\n")
print(head(ema_windows))
cat("  Participants with EMA data:", nrow(ema_windows), "\n")

# Add ema_day to step data
# First, we need to create ema_day in steps based on calendar days
steps <- steps %>%
  mutate(calendar_date = as.Date(time_bin, tz = "UTC"))

# Create a mapping from step uid to ema mlife_id
# Assuming uid in steps corresponds to mlife_id in ema
# If there's a different mapping, adjust accordingly
steps <- steps %>%
  rename(mlife_id = uid)

# Add ema_day to steps by matching calendar date ranges per participant
# For each participant, assign ema_day based on chronological order of calendar dates
steps <- steps %>%
  group_by(mlife_id) %>%
  mutate(
    # Create a day number based on unique dates
    ema_day = as.integer(factor(calendar_date, levels = unique(sort(calendar_date))))
  ) %>%
  ungroup()

cat("\nStep data with ema_day:\n")
cat("  Rows:", nrow(steps), "\n")
cat("  Unique mlife_ids:", length(unique(steps$mlife_id)), "\n")
cat("  EMA day range:", min(steps$ema_day, na.rm = TRUE), "to", 
    max(steps$ema_day, na.rm = TRUE), "\n")

# =============================================================================
# SECTION 4: CALCULATE ROLLING WINDOWS BASED ON EMA_DAY
#            ALL WINDOWS END 4 HOURS BEFORE CURRENT TIME
# =============================================================================

cat("\n=== Calculating features with ema_day-based windows ===\n")
cat("Window definitions (all end 4 hours before current time):\n")
cat("  1-hour:  5 hours before → 4 hours before\n")
cat("  4-hour:  8 hours before → 4 hours before\n")
cat("  1-day:   28 hours before → 4 hours before\n")
cat("  7-day:   172 hours before → 4 hours before\n\n")

# Split data by participant for parallel processing
uid_list <- unique(steps$mlife_id)
data_split <- split(steps, steps$mlife_id)

# Setup parallel backend
plan(multisession, workers = n_cores)

# Progress bar
handlers(global = TRUE)

# Calculate all features in parallel for each participant
cat("\nProcessing all time windows in parallel...\n")

steps_features_list <- with_progress({
  p <- progressor(along = uid_list)
  
  future_lapply(data_split, function(participant_data) {
    p(sprintf("mlife_id: %s", unique(participant_data$mlife_id)))
    
    participant_data <- participant_data %>%
      arrange(time_bin)
    
    # Initialize feature columns
    participant_data <- participant_data %>%
      mutate(
        # 1-hour windows (5h before to 4h before)
        step_mean_1h = NA_real_,
        step_mad_1h = NA_real_,
        step_p90p10_1h = NA_real_,
        
        # 4-hour windows (8h before to 4h before)
        step_mean_4h = NA_real_,
        step_mad_4h = NA_real_,
        step_p90p10_4h = NA_real_,
        step_entropy_4h = NA_real_,
        
        # 1-day windows (28h before to 4h before)
        step_mean_1d = NA_real_,
        step_mad_1d = NA_real_,
        step_p90p10_1d = NA_real_,
        step_entropy_1d = NA_real_,
        
        # 7-day windows (172h before to 4h before)
        step_mean_7d = NA_real_,
        step_mad_7d = NA_real_,
        step_p90p10_7d = NA_real_,
        step_entropy_7d = NA_real_
      )
    
    # Calculate features for each row
    for (i in seq_len(nrow(participant_data))) {
      current_time <- participant_data$time_bin[i]
      current_ema_day <- participant_data$ema_day[i]
      
      # Define the end of all windows (4 hours before current time)
      window_end <- current_time - hours(4)
      
      # 1-hour window: 5 hours before to 4 hours before
      window_start_1h <- current_time - hours(5)
      window_1h <- participant_data %>%
        filter(time_bin > window_start_1h & time_bin <= window_end)
      
      if (nrow(window_1h) > 0) {
        participant_data$step_mean_1h[i] <- calculate_mean(window_1h$total_step_count)
        participant_data$step_mad_1h[i] <- calculate_mad(window_1h$total_step_count)
        participant_data$step_p90p10_1h[i] <- p90_p10(window_1h$total_step_count)
      }
      
      # 4-hour window: 8 hours before to 4 hours before
      window_start_4h <- current_time - hours(8)
      window_4h <- participant_data %>%
        filter(time_bin > window_start_4h & time_bin <= window_end)
      
      if (nrow(window_4h) > 0) {
        participant_data$step_mean_4h[i] <- calculate_mean(window_4h$total_step_count)
        participant_data$step_mad_4h[i] <- calculate_mad(window_4h$total_step_count)
        participant_data$step_p90p10_4h[i] <- p90_p10(window_4h$total_step_count)
        participant_data$step_entropy_4h[i] <- set_entropy(window_4h$total_step_count, min_n = 3)
      }
      
      # 1-day window: 28 hours before to 4 hours before (24h + 4h)
      window_start_1d <- current_time - hours(28)
      window_1d <- participant_data %>%
        filter(time_bin > window_start_1d & time_bin <= window_end)
      
      if (nrow(window_1d) > 0) {
        participant_data$step_mean_1d[i] <- calculate_mean(window_1d$total_step_count)
        participant_data$step_mad_1d[i] <- calculate_mad(window_1d$total_step_count)
        participant_data$step_p90p10_1d[i] <- p90_p10(window_1d$total_step_count)
        participant_data$step_entropy_1d[i] <- set_entropy(window_1d$total_step_count, min_n = 5)
      }
      
      # 7-day window: 172 hours before to 4 hours before (168h + 4h = 7*24 + 4)
      window_start_7d <- current_time - hours(172)
      window_7d <- participant_data %>%
        filter(time_bin > window_start_7d & time_bin <= window_end)
      
      if (nrow(window_7d) > 0) {
        participant_data$step_mean_7d[i] <- calculate_mean(window_7d$total_step_count)
        participant_data$step_mad_7d[i] <- calculate_mad(window_7d$total_step_count)
        participant_data$step_p90p10_7d[i] <- p90_p10(window_7d$total_step_count)
        participant_data$step_entropy_7d[i] <- set_entropy(window_7d$total_step_count, min_n = 10)
      }
    }
    
    return(participant_data)
  }, future.seed = TRUE)
})

# Combine results
steps_features <- bind_rows(steps_features_list)

cat("\nAll features calculated based on time windows ending 4 hours before current time\n")

# =============================================================================
# SECTION 5: FEATURE SUMMARY AND QUALITY CHECKS
# =============================================================================

cat("\n=== Feature Summary ===\n")

# List all features created
feature_cols <- c(
  # 1-hour features
  "step_mean_1h", "step_mad_1h", "step_p90p10_1h",
  
  # 4-hour features
  "step_mean_4h", "step_mad_4h", "step_p90p10_4h", "step_entropy_4h",
  
  # 1-day features
  "step_mean_1d", "step_mad_1d", "step_p90p10_1d", "step_entropy_1d",
  
  # 7-day features
  "step_mean_7d", "step_mad_7d", "step_p90p10_7d", "step_entropy_7d"
)

# Calculate completeness for each feature
completeness <- steps_features %>%
  summarise(across(all_of(feature_cols), 
                   ~sprintf("%.1f%%", 100 * mean(!is.na(.x)))))

cat("\nFeature Completeness:\n")
print(t(completeness))

# Summary statistics
cat("\nSummary Statistics:\n")
summary_stats <- steps_features %>%
  select(all_of(feature_cols)) %>%
  summary()

print(summary_stats)

# Check for issues
cat("\nData Quality Checks:\n")

# Check for infinite values
inf_check <- steps_features %>%
  summarise(across(all_of(feature_cols), 
                   ~sum(is.infinite(.x))))

if (any(inf_check > 0)) {
  cat("WARNING: Infinite values detected:\n")
  print(inf_check[, inf_check > 0])
} else {
  cat("✓ No infinite values detected\n")
}

# Check for NaN values
nan_check <- steps_features %>%
  summarise(across(all_of(feature_cols), 
                   ~sum(is.nan(.x))))

if (any(nan_check > 0)) {
  cat("WARNING: NaN values detected:\n")
  print(nan_check[, nan_check > 0])
} else {
  cat("✓ No NaN values detected\n")
}

# =============================================================================
# SECTION 6: SAVE RESULTS
# =============================================================================

cat("\n=== Saving results ===\n")

# Save full dataset with all features
output_file <- "step_features_all_windows_lag4h.csv"
fwrite(steps_features, output_file)
cat("✓ Saved complete dataset to:", output_file, "\n")
cat("  Rows:", nrow(steps_features), "\n")
cat("  Columns:", ncol(steps_features), "\n")

# Save feature-only dataset
features_only <- steps_features %>%
  select(mlife_id, time_bin, ema_day, calendar_date, total_step_count, all_of(feature_cols))

features_file <- "step_features_only_lag4h.csv"
fwrite(features_only, features_file)
cat("✓ Saved features-only dataset to:", features_file, "\n")

# Save summary statistics
summary_df <- data.frame(
  feature = feature_cols,
  mean = sapply(steps_features[, feature_cols, with = FALSE], mean, na.rm = TRUE),
  sd = sapply(steps_features[, feature_cols, with = FALSE], sd, na.rm = TRUE),
  min = sapply(steps_features[, feature_cols, with = FALSE], min, na.rm = TRUE),
  max = sapply(steps_features[, feature_cols, with = FALSE], max, na.rm = TRUE),
  missing_pct = sapply(steps_features[, feature_cols, with = FALSE], function(x) 100 * mean(is.na(x)))
)

summary_file <- "step_features_summary_lag4h.csv"
fwrite(summary_df, summary_file)
cat("✓ Saved feature summary to:", summary_file, "\n")

# =============================================================================
# SECTION 7: VISUALIZATION (OPTIONAL)
# =============================================================================

cat("\n=== Creating visualizations ===\n")

# Example: Distribution of mean step counts across time windows
plot_data <- steps_features %>%
  select(mlife_id, time_bin, step_mean_1h, step_mean_4h, step_mean_1d, step_mean_7d) %>%
  pivot_longer(
    cols = starts_with("step_mean"),
    names_to = "window",
    values_to = "mean_steps"
  ) %>%
  mutate(
    window = factor(
      window,
      levels = c("step_mean_1h", "step_mean_4h", "step_mean_1d", "step_mean_7d"),
      labels = c("1 hour", "4 hours", "1 day", "7 days")
    )
  )

p <- ggplot(plot_data, aes(x = window, y = mean_steps, fill = window)) +
  geom_violin(alpha = 0.7) +
  geom_boxplot(width = 0.2, outlier.alpha = 0.3) +
  scale_y_log10() +
  labs(
    title = "Distribution of Mean Step Counts Across Time Windows",
    subtitle = "All windows end 4 hours before current time",
    x = "Time Window",
    y = "Mean Step Count (log scale)",
    fill = "Window"
  ) +
  theme_minimal() +
  theme(legend.position = "none")

ggsave("step_features_distribution_lag4h.png", p, width = 10, height = 6, dpi = 300)
cat("✓ Saved distribution plot to: step_features_distribution_lag4h.png\n")

# =============================================================================
# SECTION 8: COMPLETION MESSAGE
# =============================================================================

cat("\n" , rep("=", 80), "\n", sep = "")
cat("STEP COUNT FEATURE ENGINEERING COMPLETE (FULLY PARALLEL)\n")
cat(rep("=", 80), "\n", sep = "")

cat("\nProcessing Summary:\n")
cat("  Total participants:", length(unique(steps_features$mlife_id)), "\n")
cat("  Total observations:", nrow(steps_features), "\n")
cat("  Features created:", length(feature_cols), "\n")
cat("  Time windows: 1h, 4h, 1d, 7d (ALL ending 4 hours before current time)\n")
cat("  Metrics: Mean, MAD, P90-P10, Entropy\n")
cat("  Cores used:", n_cores, "\n")

cat("\nWindow Definitions:\n")
cat("  1-hour:  (current_time - 5h) to (current_time - 4h)\n")
cat("  4-hour:  (current_time - 8h) to (current_time - 4h)\n")
cat("  1-day:   (current_time - 28h) to (current_time - 4h)\n")
cat("  7-day:   (current_time - 172h) to (current_time - 4h)\n")

cat("\nOutput Files:\n")
cat("  1.", output_file, "- Complete dataset with all features\n")
cat("  2.", features_file, "- Features only\n")
cat("  3.", summary_file, "- Feature summary statistics\n")
cat("  4. step_features_distribution_lag4h.png - Visualization\n")

cat("\n✓ All processing complete!\n\n")

# Clean up parallel backend
plan(sequential)
