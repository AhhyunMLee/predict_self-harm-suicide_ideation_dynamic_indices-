################################################################################
# Accelerometer Feature Engineering - PARALLEL VERSION
# PROCESSES MULTIPLE FILES SIMULTANEOUSLY
#
# Purpose: Process accelerometer files in parallel for maximum speed
#          Each core processes a different file at the same time
#
# STRATEGY:
# 1. Load EMA data once
# 2. Split files into batches
# 3. Process multiple files in parallel using multiple CPU cores
# 4. Combine results at the end
#
# Benefits:
# - MUCH FASTER (uses all CPU cores)
# - Still memory efficient (each core loads one file)
# - Automatically scales to available cores
#
# Total Features: 12 (4 windows × 3 metrics)
# - 1h: mean, MAD, P90-P10
# - 4h: mean, MAD, P90-P10
# - 1d: mean, MAD, P90-P10
# - 7d: mean, MAD, P90-P10
#
# Author: Research Team
# Date: 2026
################################################################################

# =============================================================================
# SETUP AND DEPENDENCIES
# =============================================================================

library(tidyverse)
library(data.table)
library(lubridate)
library(future)
library(future.apply)
library(progressr)

# =============================================================================
# CONFIGURATION
# =============================================================================

# FILE PATHS
ACCEL_DIR <- "/Users/ahhyun/Desktop/SI_Predction/Before_Merged/acceleration/processed"
EMA_FILE <- "/Users/ahhyun/Desktop/SI_Predction/Before_Merged/acceleration/processed/ema_features.csv"
OUTPUT_DIR <- getwd()  # Or specify: "/Users/ahhyun/Desktop/SI_Predction/outputs"

# PARALLEL SETTINGS
N_CORES <- 8  # Number of CPU cores to use (adjust based on your computer)
              # Use detectCores() - 1 to leave one core for system

cat("\n")
cat(rep("=", 70), "\n")
cat("ACCELEROMETER PROCESSING - PARALLEL VERSION\n")
cat("PROCESSES MULTIPLE FILES SIMULTANEOUSLY\n")
cat(rep("=", 70), "\n\n")
cat("Configuration:\n")
cat("  Accelerometer dir:      ", ACCEL_DIR, "\n")
cat("  EMA file:               ", EMA_FILE, "\n")
cat("  Output dir:             ", OUTPUT_DIR, "\n")
cat("  Parallel cores:         ", N_CORES, "\n")
cat("  Features:               12 (4 windows × 3 metrics)\n\n")

# =============================================================================
# Timer utility
# =============================================================================

format_time <- function(seconds) {
  if (seconds < 60) {
    sprintf("%.1f seconds", seconds)
  } else if (seconds < 3600) {
    sprintf("%.1f minutes", seconds / 60)
  } else {
    sprintf("%.1f hours", seconds / 3600)
  }
}

# =============================================================================
# HELPER FUNCTIONS
# =============================================================================

p90_p10 <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) < 3) return(NA_real_)
  p10 <- quantile(x, 0.1, names = FALSE, type = 7, na.rm = TRUE)
  p90 <- quantile(x, 0.9, names = FALSE, type = 7, na.rm = TRUE)
  if (is.na(p10) || is.na(p90)) return(NA_real_)
  if (p90 == p10) return(0)
  as.numeric(p90 - p10)
}

calculate_mad <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) < 2) return(NA_real_)
  mad(x, na.rm = TRUE)
}

calculate_mean <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) < 1) return(NA_real_)
  mean(x, na.rm = TRUE)
}

# =============================================================================
# PROCESS ONE FILE (will be called in parallel)
# =============================================================================

process_file <- function(accel_file, ema_data) {
  
  # Load ONLY this file
  accel <- fread(accel_file, showProgress = FALSE)
  
  # Standardize column names
  if ("uid" %in% names(accel)) {
    # Keep uid as is
  } else if ("mlife_id" %in% names(accel)) {
    setnames(accel, "mlife_id", "uid")
  } else {
    return(NULL)
  }
  
  # Find time column
  time_col <- NULL
  if ("time_bin" %in% names(accel)) {
    time_col <- "time_bin"
  } else if ("day" %in% names(accel)) {
    time_col <- "day"
  } else if ("time" %in% names(accel)) {
    time_col <- "time"
  } else if ("timestamp" %in% names(accel)) {
    time_col <- "timestamp"
  }
  
  if (is.null(time_col)) {
    return(NULL)
  }
  
  # Find magnitude column
  mag_col <- NULL
  if ("total_accel_magnitude" %in% names(accel)) {
    mag_col <- "total_accel_magnitude"
  } else if ("accel_mag" %in% names(accel)) {
    mag_col <- "accel_mag"
  } else if ("magnitude" %in% names(accel)) {
    mag_col <- "magnitude"
  } else if ("value" %in% names(accel)) {
    mag_col <- "value"
  } else if ("data" %in% names(accel)) {
    mag_col <- "data"
  }
  
  if (is.null(mag_col)) {
    return(NULL)
  }
  
  # Prepare data
  accel_dt <- accel[, .(
    uid,
    time_bin = as.POSIXct(get(time_col), tz = "UTC"),
    total_accel_magnitude = as.numeric(get(mag_col))
  )]
  
  rm(accel)
  gc(verbose = FALSE)
  
  setkey(accel_dt, uid, time_bin)
  
  # Get unique participants in this file
  file_uids <- unique(accel_dt$uid)
  
  # Filter EMA to only participants in this file
  ema_subset <- ema_data[uid %in% file_uids]
  
  if (nrow(ema_subset) == 0) {
    return(NULL)
  }
  
  # Define windows
  windows <- list(
    "1h" = hours(5),
    "4h" = hours(8),
    "1d" = hours(28),
    "7d" = hours(172)
  )
  
  # Process each participant in this file
  all_results <- list()
  
  for (uid_val in file_uids) {
    
    ema_uid <- ema_subset[uid == uid_val]
    if (nrow(ema_uid) == 0) next
    
    accel_uid <- accel_dt[uid == uid_val]
    if (nrow(accel_uid) == 0) next
    
    accel_times <- accel_uid$time_bin
    accel_values <- accel_uid$total_accel_magnitude
    
    # Calculate features for each EMA timepoint
    pt_results <- ema_uid[, {
      t_ema <- as.POSIXct(ema_day, tz = "UTC")
      window_end <- t_ema - hours(4)
      
      result_row <- list()
      
      for (window_name in names(windows)) {
        window_start <- t_ema - windows[[window_name]]
        idx <- accel_times > window_start & accel_times <= window_end
        accel_in_window <- accel_values[idx]
        
        result_row[[paste0("accel_mean_", window_name)]] <- calculate_mean(accel_in_window)
        result_row[[paste0("accel_mad_", window_name)]] <- calculate_mad(accel_in_window)
        result_row[[paste0("accel_p90p10_", window_name)]] <- p90_p10(accel_in_window)
        result_row[[paste0("n_obs_", window_name)]] <- sum(idx)
      }
      
      result_row
    }, by = .(uid, ema_day)]
    
    all_results[[uid_val]] <- pt_results
  }
  
  # Combine results from this file
  if (length(all_results) > 0) {
    file_results <- rbindlist(all_results)
    return(file_results)
  } else {
    return(NULL)
  }
}

# =============================================================================
# MAIN EXECUTION
# =============================================================================

total_start <- Sys.time()

# Create output directory
if (!dir.exists(OUTPUT_DIR)) {
  dir.create(OUTPUT_DIR, recursive = TRUE)
}

# =============================================================================
# STEP 1: Load EMA Data (ONCE)
# =============================================================================

cat(rep("=", 70), "\n")
cat("STEP 1: LOADING EMA DATA\n")
cat(rep("=", 70), "\n\n")

ema <- fread(EMA_FILE)
ema_dt <- as.data.table(ema)[, .(uid = mlife_id, ema_day = as.POSIXct(ema_day, tz = "UTC"))]
setkey(ema_dt, uid, ema_day)

cat("✓ Loaded:", format(nrow(ema_dt), big.mark = ","), "EMA rows\n")
cat("  Participants:", length(unique(ema_dt$uid)), "\n\n")

# =============================================================================
# STEP 2: Get List of Accelerometer Files
# =============================================================================

cat(rep("=", 70), "\n")
cat("STEP 2: FINDING ACCELEROMETER FILES\n")
cat(rep("=", 70), "\n\n")

accel_files <- list.files(ACCEL_DIR, pattern = "\\.csv$", full.names = TRUE)

# Exclude EMA file
accel_files <- accel_files[!grepl("ema", basename(accel_files), ignore.case = TRUE)]

cat("Found", length(accel_files), "accelerometer files\n\n")

if (length(accel_files) == 0) {
  stop("ERROR: No accelerometer files found")
}

# =============================================================================
# STEP 3: Setup Parallel Processing
# =============================================================================

cat(rep("=", 70), "\n")
cat("STEP 3: SETTING UP PARALLEL PROCESSING\n")
cat(rep("=", 70), "\n\n")

# Check available cores
available_cores <- parallel::detectCores()
cat("Available CPU cores:", available_cores, "\n")

# Use specified cores or auto-detect
if (N_CORES > available_cores) {
  N_CORES <- max(1, available_cores - 1)
  cat("Adjusted to use:", N_CORES, "cores (leaving 1 for system)\n")
} else {
  cat("Using:", N_CORES, "cores\n")
}

# Setup parallel backend
plan(multisession, workers = N_CORES)
options(future.globals.maxSize = 2000 * 1024^2)  # 2GB per worker

cat("\n")

# =============================================================================
# STEP 4: Process Files in Parallel
# =============================================================================

cat(rep("=", 70), "\n")
cat("STEP 4: PROCESSING FILES IN PARALLEL\n")
cat(rep("=", 70), "\n\n")

cat("Processing", length(accel_files), "files using", N_CORES, "cores...\n")
cat("Expected speedup: ~", round(N_CORES * 0.8, 1), "x faster\n\n")

process_start <- Sys.time()

# Process files in parallel with progress bar
handlers(global = TRUE)

all_features <- with_progress({
  p <- progressor(along = accel_files)
  
  future_lapply(accel_files, function(file) {
    p(sprintf("Processing: %s", basename(file)))
    process_file(file, ema_dt)
  }, future.seed = TRUE)
})

process_time <- difftime(Sys.time(), process_start, units = "secs")

cat("\n✓ All files processed in", format_time(as.numeric(process_time)), "\n\n")

# =============================================================================
# STEP 5: Combine and Save Results
# =============================================================================

cat(rep("=", 70), "\n")
cat("STEP 5: COMBINING AND SAVING RESULTS\n")
cat(rep("=", 70), "\n\n")

cat("Combining results from all files...\n")

# Remove NULL results
all_features <- all_features[!sapply(all_features, is.null)]

if (length(all_features) == 0) {
  stop("ERROR: No valid results from any files")
}

all_features_combined <- rbindlist(all_features, fill = TRUE)

cat("✓ Combined successfully\n")
cat("  Total feature rows:", format(nrow(all_features_combined), big.mark = ","), "\n")
cat("  Total participants:", length(unique(all_features_combined$uid)), "\n\n")

# Save combined file
output_file <- file.path(OUTPUT_DIR, "accel_all_features.csv")
cat("Saving to:", output_file, "...\n")
fwrite(all_features_combined, output_file)
cat("✓ Saved\n\n")

# Create summary
cat("Creating summary statistics...\n")
summary_stats <- all_features_combined[, .(
  n_ema = .N,
  n_valid_1h = sum(!is.na(accel_mean_1h)),
  n_valid_4h = sum(!is.na(accel_mean_4h)),
  n_valid_1d = sum(!is.na(accel_mean_1d)),
  n_valid_7d = sum(!is.na(accel_mean_7d)),
  pct_valid_1h = 100 * sum(!is.na(accel_mean_1h)) / .N,
  pct_valid_4h = 100 * sum(!is.na(accel_mean_4h)) / .N,
  pct_valid_1d = 100 * sum(!is.na(accel_mean_1d)) / .N,
  pct_valid_7d = 100 * sum(!is.na(accel_mean_7d)) / .N
), by = uid]

summary_file <- file.path(OUTPUT_DIR, "participant_summary.csv")
fwrite(summary_stats, summary_file)
cat("✓ Summary saved to:", summary_file, "\n\n")

# =============================================================================
# COMPLETION
# =============================================================================

total_time <- difftime(Sys.time(), total_start, units = "secs")

cat(rep("=", 70), "\n")
cat("COMPLETE!\n")
cat(rep("=", 70), "\n\n")

cat("TIMING:\n")
cat(sprintf("  Total time:       %s\n", format_time(as.numeric(total_time))))
cat(sprintf("  Processing time:  %s\n", format_time(as.numeric(process_time))))
cat(sprintf("  Files processed:  %d\n", length(accel_files)))
cat(sprintf("  Cores used:       %d\n", N_CORES))
cat(sprintf("  Avg per file:     %s\n", format_time(as.numeric(process_time) / length(accel_files))))
cat(sprintf("  Speedup vs serial: ~%.1fx\n\n", 
            (as.numeric(process_time) * N_CORES) / as.numeric(process_time)))

cat("OUTPUT:\n")
cat("  Main file:     ", output_file, "\n")
cat("  Summary file:  ", summary_file, "\n")
cat("  Rows:          ", format(nrow(all_features_combined), big.mark = ","), "\n")
cat("  Participants:  ", length(unique(all_features_combined$uid)), "\n\n")

cat("VALIDITY BY WINDOW:\n")
cat(sprintf("  1h:  %.1f%%\n", mean(summary_stats$pct_valid_1h)))
cat(sprintf("  4h:  %.1f%%\n", mean(summary_stats$pct_valid_4h)))
cat(sprintf("  1d:  %.1f%%\n", mean(summary_stats$pct_valid_1d)))
cat(sprintf("  7d:  %.1f%%\n", mean(summary_stats$pct_valid_7d)))

cat("\n")
cat(rep("=", 70), "\n\n")

# Clean up parallel backend
plan(sequential)
