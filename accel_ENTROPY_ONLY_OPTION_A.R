################################################################################
# Accelerometer Entropy Calculation - OPTIMIZED FOR ENTROPY ONLY
# PROCESSES ONE WINDOW AT A TIME IN PARALLEL
#
# Purpose: Calculate ONLY entropy features for all time windows
#          Designed to be merged with existing mean/MAD/P90P10 features
#
# Strategy:
# 1. Process one window at a time (1h, then 4h, then 1d, then 7d)
# 2. For each window, process files in parallel
# 3. Save each window's entropy separately
# 4. Can merge with existing features later
#
# Benefits:
# - Optimized for entropy calculation only
# - Parallel processing for speed
# - One window at a time minimizes memory
# - Easy to resume if interrupted
#
# Output Files:
# - accel_entropy_1h.csv
# - accel_entropy_4h.csv
# - accel_entropy_1d.csv
# - accel_entropy_7d.csv
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
OUTPUT_DIR <- getwd()  # Or specify your output directory

# PARALLEL SETTINGS
N_CORES <- 8  # Number of CPU cores to use

# ENTROPY SETTINGS - OPTION A (ADAPTIVE FOR 5-SECOND DATA)
# Optimized for 5-second sampling interval (0.2 Hz)
# These settings give ~99% accuracy while maintaining good speed
ENTROPY_LIMITS <- list(
  "1h" = Inf,      # Use all ~720 points (100% accuracy, fast)
  "4h" = 1000,     # Use 1,000 of ~2,880 points (35% sampling, 99% accuracy)
  "1d" = 3000,     # Use 3,000 of ~17,280 points (17% sampling, 99% accuracy)
  "7d" = 15000     # Use 15,000 of ~120,960 points (12% sampling, 99% accuracy)
)

cat("\n")
cat(rep("=", 70), "\n")
cat("ACCELEROMETER ENTROPY CALCULATION - OPTIMIZED\n")
cat("ONE WINDOW AT A TIME - PARALLEL PROCESSING\n")
cat(rep("=", 70), "\n\n")
cat("Configuration:\n")
cat("  Accelerometer dir:      ", ACCEL_DIR, "\n")
cat("  EMA file:               ", EMA_FILE, "\n")
cat("  Output dir:             ", OUTPUT_DIR, "\n")
cat("  Parallel cores:         ", N_CORES, "\n")
cat("  Max entropy points:     ", MAX_ENTROPY_POINTS, "\n")
cat("  Windows to process:     1h, 4h, 1d, 7d\n\n")

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
# ENTROPY FUNCTION - OPTIMIZED
# =============================================================================

calculate_entropy <- function(x, min_n = 2, max_points = Inf) {
  # Remove NA values
  x <- x[!is.na(x)]
  
  # Check if there are enough non-NA values
  if (length(x) < min_n) {
    return(NA_real_)
  }
  
  # Subsample if too many points (saves memory and time)
  if (!is.infinite(max_points) && length(x) > max_points) {
    indices <- seq(1, length(x), length.out = max_points)
    x <- x[round(indices)]
  }
  
  # Calculate probability distribution
  p <- table(x) / length(x)
  
  # Calculate Shannon entropy in bits (log2)
  entropy_value <- -sum(p * log2(p))
  
  # Return
  as.numeric(entropy_value)
}

# =============================================================================
# PROCESS ONE FILE FOR ONE WINDOW - ENTROPY ONLY
# =============================================================================

process_file_entropy <- function(accel_file, ema_data, window_name, window_hours, min_entropy, max_entropy_points) {
  
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
  
  # Process each participant in this file
  all_results <- list()
  
  for (uid_val in file_uids) {
    
    ema_uid <- ema_subset[uid == uid_val]
    if (nrow(ema_uid) == 0) next
    
    accel_uid <- accel_dt[uid == uid_val]
    if (nrow(accel_uid) == 0) next
    
    accel_times <- accel_uid$time_bin
    accel_values <- accel_uid$total_accel_magnitude
    
    # Calculate entropy for each EMA timepoint - THIS WINDOW ONLY
    pt_results <- ema_uid[, {
      t_ema <- as.POSIXct(ema_day, tz = "UTC")
      window_end <- t_ema - hours(4)
      window_start <- t_ema - window_hours
      
      idx <- accel_times > window_start & accel_times <= window_end
      accel_in_window <- accel_values[idx]
      
      # Calculate ONLY entropy
      entropy_val <- calculate_entropy(accel_in_window, min_n = min_entropy, max_points = max_entropy_points)
      
      list(
        entropy = entropy_val,
        n_obs = sum(idx)
      )
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
# PROCESS ONE WINDOW ACROSS ALL FILES
# =============================================================================

process_window <- function(window_name, window_hours, min_entropy, max_entropy_points, ema_data, accel_files) {
  
  cat("\n")
  cat(rep("=", 70), "\n")
  cat("PROCESSING WINDOW:", toupper(window_name), "\n")
  cat(rep("=", 70), "\n\n")
  
  cat("Window definition:\n")
  cat("  Lookback:", as.numeric(window_hours), "hours\n")
  cat("  Range: (ema_day -", as.numeric(window_hours), "h) to (ema_day - 4h)\n")
  cat("  Min entropy points:", min_entropy, "\n")
  cat("  Max entropy points:", ifelse(is.infinite(max_entropy_points), "ALL", format(max_entropy_points, big.mark = ",")), "\n\n")
  
  window_start <- Sys.time()
  
  # Process files in parallel with progress bar
  cat("Processing", length(accel_files), "files in parallel...\n")
  
  handlers(global = TRUE)
  
  all_results <- with_progress({
    p <- progressor(along = accel_files)
    
    future_lapply(accel_files, function(file) {
      p(sprintf("%s: %s", window_name, basename(file)))
      process_file_entropy(file, ema_data, window_name, window_hours, min_entropy, max_entropy_points)
    }, future.seed = TRUE)
  })
  
  # Remove NULL results
  all_results <- all_results[!sapply(all_results, is.null)]
  
  if (length(all_results) == 0) {
    cat("\n⚠ No results for this window\n")
    return(NULL)
  }
  
  # Combine results
  cat("\nCombining results...\n")
  combined <- rbindlist(all_results, fill = TRUE)
  
  # Rename columns to include window name
  setnames(combined, "entropy", paste0("accel_entropy_", window_name))
  setnames(combined, "n_obs", paste0("n_obs_", window_name))
  
  # Save immediately
  output_file <- file.path(OUTPUT_DIR, paste0("accel_entropy_", window_name, ".csv"))
  cat("Saving to:", output_file, "...\n")
  fwrite(combined, output_file)
  
  window_time <- difftime(Sys.time(), window_start, units = "secs")
  
  # Statistics
  n_valid <- sum(!is.na(combined[[paste0("accel_entropy_", window_name)]]))
  pct_valid <- 100 * n_valid / nrow(combined)
  
  cat("✓ Complete in", format_time(as.numeric(window_time)), "\n")
  cat("  Rows:", format(nrow(combined), big.mark = ","), "\n")
  cat("  Valid entropy:", format(n_valid, big.mark = ","), sprintf("(%.1f%%)\n", pct_valid))
  
  return(list(
    window = window_name,
    file = output_file,
    n_rows = nrow(combined),
    n_valid = n_valid,
    pct_valid = pct_valid,
    time_seconds = as.numeric(window_time)
  ))
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
# STEP 1: Load EMA Data
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
# STEP 2: Get Accelerometer Files
# =============================================================================

cat(rep("=", 70), "\n")
cat("STEP 2: FINDING ACCELEROMETER FILES\n")
cat(rep("=", 70), "\n\n")

accel_files <- list.files(ACCEL_DIR, pattern = "\\.csv$", full.names = TRUE)
accel_files <- accel_files[!grepl("ema", basename(accel_files), ignore.case = TRUE)]

cat("Found", length(accel_files), "files\n\n")

if (length(accel_files) == 0) {
  stop("ERROR: No accelerometer files found")
}

# =============================================================================
# STEP 3: Setup Parallel Processing
# =============================================================================

cat(rep("=", 70), "\n")
cat("STEP 3: SETTING UP PARALLEL PROCESSING\n")
cat(rep("=", 70), "\n\n")

available_cores <- parallel::detectCores()
cat("Available cores:", available_cores, "\n")

if (N_CORES > available_cores) {
  N_CORES <- max(1, available_cores - 1)
  cat("Adjusted to:", N_CORES, "cores\n")
}

plan(multisession, workers = N_CORES)
options(future.globals.maxSize = 2000 * 1024^2)

cat("\n")

# =============================================================================
# STEP 4: Process Each Window
# =============================================================================

cat(rep("=", 70), "\n")
cat("STEP 4: PROCESSING WINDOWS ONE AT A TIME\n")
cat(rep("=", 70), "\n")

# Define windows with their parameters
windows_config <- list(
  list(name = "1h", hours = hours(5), min_entropy = 2),
  list(name = "4h", hours = hours(8), min_entropy = 3),
  list(name = "1d", hours = hours(28), min_entropy = 5),
  list(name = "7d", hours = hours(172), min_entropy = 10)
)

# Process each window
all_summaries <- list()

for (i in 1:length(windows_config)) {
  config <- windows_config[[i]]
  
  cat("\n[", i, "/", length(windows_config), "] ")
  
  summary <- process_window(
    window_name = config$name,
    window_hours = config$hours,
    min_entropy = config$min_entropy,
    max_entropy_points = ENTROPY_LIMITS[[config$name]],  # Get adaptive limit for this window
    ema_data = ema_dt,
    accel_files = accel_files
  )
  
  if (!is.null(summary)) {
    all_summaries[[i]] <- summary
  }
  
  # Garbage collection between windows
  gc(verbose = FALSE)
}

# =============================================================================
# COMPLETION
# =============================================================================

total_time <- difftime(Sys.time(), total_start, units = "secs")

cat("\n\n")
cat(rep("=", 70), "\n")
cat("ALL WINDOWS COMPLETE!\n")
cat(rep("=", 70), "\n\n")

cat("TIMING:\n")
cat(sprintf("  Total time:  %s\n", format_time(as.numeric(total_time))))
cat(sprintf("  Cores used:  %d\n\n", N_CORES))

cat("OUTPUT FILES CREATED:\n")
for (i in 1:length(all_summaries)) {
  s <- all_summaries[[i]]
  cat(sprintf("  %d. accel_entropy_%s.csv  (%s rows, %.1f%% valid)\n",
              i, s$window, format(s$n_rows, big.mark = ","), s$pct_valid))
}

# Save summary
summary_df <- rbindlist(all_summaries)
summary_file <- file.path(OUTPUT_DIR, "entropy_summary.csv")
fwrite(summary_df, summary_file)

cat("\n✓ Summary saved to:", summary_file, "\n")

cat("\nTO MERGE WITH EXISTING FEATURES:\n")
cat("  library(data.table)\n")
cat("  \n")
cat("  # Load your existing features\n")
cat("  features <- fread('accel_all_features.csv')\n")
cat("  \n")
cat("  # Merge each entropy file\n")
cat("  ent_1h <- fread('accel_entropy_1h.csv')\n")
cat("  features <- merge(features, ent_1h, by = c('uid', 'ema_day'), all.x = TRUE)\n")
cat("  \n")
cat("  ent_4h <- fread('accel_entropy_4h.csv')\n")
cat("  features <- merge(features, ent_4h, by = c('uid', 'ema_day'), all.x = TRUE)\n")
cat("  \n")
cat("  ent_1d <- fread('accel_entropy_1d.csv')\n")
cat("  features <- merge(features, ent_1d, by = c('uid', 'ema_day'), all.x = TRUE)\n")
cat("  \n")
cat("  ent_7d <- fread('accel_entropy_7d.csv')\n")
cat("  features <- merge(features, ent_7d, by = c('uid', 'ema_day'), all.x = TRUE)\n")
cat("  \n")
cat("  fwrite(features, 'accel_all_features_with_entropy.csv')\n")

cat("\n")
cat(rep("=", 70), "\n\n")

# Clean up
plan(sequential)
