################################################################################
# Step Count Feature Engineering - TIMED VERSION (EMA-based)
# Purpose: Calculate dynamic features from step count data using EMA timepoints
#          Windows END 4 hours before ema_day time (similar to HR processing)
# Time Windows: 
#   - 1 hour: (5 hours before to 4 hours before)
#   - 4 hours: (8 hours before to 4 hours before)
#   - 1 day: (28 hours before to 4 hours before)
#   - 7 days: (172 hours before to 4 hours before)
# Metrics: Mean, MAD, P90-P10, Entropy
# Author: Research Team
# Date: 2026
################################################################################

# =============================================================================
# SECTION 1: SETUP AND DEPENDENCIES
# =============================================================================

library(tidyverse)
library(data.table)
library(lubridate)
library(future)
library(furrr)

# =============================================================================
# CONFIGURATION
# =============================================================================

N_CORES <- 12  # Number of cores for parallel processing

cat("\n")
cat(rep("=", 70), "\n")
cat("STEP COUNT PROCESSING - TIMED VERSION (EMA-based)\n")
cat(rep("=", 70), "\n\n")
cat("Configuration:\n")
cat("  Cores:                  ", N_CORES, "\n\n")

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
# SECTION 2: HELPER FUNCTIONS
# =============================================================================

## Set Entropy (Shannon Entropy) ----
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

## P90-P10 Range ----
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

## Median Absolute Deviation (MAD) ----
calculate_mad <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) < 2) return(NA_real_)
  mad(x, na.rm = TRUE)
}

## Mean ----
calculate_mean <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) < 1) return(NA_real_)
  mean(x, na.rm = TRUE)
}

## Calculate all metrics at once ----
calculate_all_metrics <- function(x, min_entropy = 2) {
  x_clean <- x[!is.na(x)]
  n <- length(x_clean)
  
  if (n == 0) {
    return(list(mean = NA_real_, MAD = NA_real_, P90P10 = NA_real_, entropy = NA_real_))
  }
  
  m <- mean(x_clean)
  mad_val <- if (n >= 2) mad(x_clean, na.rm = TRUE) else NA_real_
  p90p10_val <- p90_p10(x_clean)
  entropy <- set_entropy(x_clean, min_n = min_entropy)
  
  list(mean = m, MAD = mad_val, P90P10 = p90p10_val, entropy = entropy)
}

# =============================================================================
# Process all windows for a single UID
# =============================================================================

process_all_windows_for_uid <- function(uid_val, ema_uid, step_batch) {
  
  step_uid <- step_batch[uid == uid_val]
  
  if (nrow(step_uid) == 0) {
    empty_result <- data.table(uid = uid_val, ema_day = ema_uid$ema_day)
    for (suffix in c("1h", "4h", "1d", "7d")) {
      empty_result[, paste0("step_mean_", suffix) := NA_real_]
      empty_result[, paste0("step_mad_", suffix) := NA_real_]
      empty_result[, paste0("step_p90p10_", suffix) := NA_real_]
      if (suffix %in% c("4h", "1d", "7d")) {
        empty_result[, paste0("step_entropy_", suffix) := NA_real_]
      }
      empty_result[, paste0("n_obs_", suffix) := 0L]
    }
    return(empty_result)
  }
  
  # Get step times and values
  step_times <- as.POSIXct(step_uid$time_bin, tz = "UTC")
  step_values <- step_uid$total_step_count
  
  # Define windows (lookback periods from current time)
  windows <- list(
    "1h" = hours(5),   # 5 hours before to 4 hours before
    "4h" = hours(8),   # 8 hours before to 4 hours before
    "1d" = hours(28),  # 28 hours before to 4 hours before (24h + 4h)
    "7d" = hours(172)  # 172 hours before to 4 hours before (168h + 4h)
  )
  
  # Minimum entropy thresholds per window
  min_entropy <- list("1h" = 2, "4h" = 3, "1d" = 5, "7d" = 10)
  
  # Process each EMA timepoint
  results <- ema_uid[, {
    t_ema <- as.POSIXct(ema_day, tz = "UTC")
    window_end <- t_ema - hours(4)
    
    result_row <- list()
    
    for (suffix in names(windows)) {
      window_start <- t_ema - windows[[suffix]]
      idx <- step_times > window_start & step_times <= window_end
      steps_in_window <- step_values[idx]
      
      metrics <- calculate_all_metrics(steps_in_window, 
                                       min_entropy = min_entropy[[suffix]])
      
      result_row[[paste0("step_mean_", suffix)]] <- metrics$mean
      result_row[[paste0("step_mad_", suffix)]] <- metrics$MAD
      result_row[[paste0("step_p90p10_", suffix)]] <- metrics$P90P10
      
      # Only add entropy for 4h, 1d, 7d windows (not 1h)
      if (suffix %in% c("4h", "1d", "7d")) {
        result_row[[paste0("step_entropy_", suffix)]] <- metrics$entropy
      }
      
      result_row[[paste0("n_obs_", suffix)]] <- sum(idx)
    }
    
    result_row
  }, by = .(uid, ema_day)]
  
  return(results)
}

# =============================================================================
# MAIN PROCESSING
# =============================================================================

total_start <- Sys.time()

cat("Setting up parallel processing...\n")
plan(multisession, workers = N_CORES)
options(future.globals.maxSize = 4000 * 1024^2)
cat("✓ Using", N_CORES, "cores\n\n")

# =============================================================================
# Load Data
# =============================================================================

cat("STEP 1: Loading data...\n")
load_start <- Sys.time()

# Load EMA data with ema_day column
ema <- fread("ema_features.csv")

# Load step count data
steps <- fread("steps_clean.csv")

# Prepare EMA data - select uid (mlife_id) and ema_day
ema_dt <- as.data.table(ema)[, .(uid = mlife_id, ema_day = as.POSIXct(ema_day, tz = "UTC"))]
setkey(ema_dt, uid, ema_day)

# Prepare step data
# Rename uid column if needed to match ema (assuming uid in steps = mlife_id in ema)
if ("uid" %in% names(steps)) {
  step_dt <- as.data.table(steps)[, .(uid, time_bin = as.POSIXct(time_bin, tz = "UTC"), 
                                       total_step_count)]
} else if ("mlife_id" %in% names(steps)) {
  step_dt <- as.data.table(steps)[, .(uid = mlife_id, time_bin = as.POSIXct(time_bin, tz = "UTC"), 
                                       total_step_count)]
} else {
  stop("ERROR: Cannot find uid or mlife_id column in steps data")
}

setkey(step_dt, uid, time_bin)

load_time <- difftime(Sys.time(), load_start, units = "secs")
cat("✓ Data loaded in", format_time(as.numeric(load_time)), "\n")
cat("  EMA rows:", format(nrow(ema_dt), big.mark = ","), "\n")
cat("  Step rows:", format(nrow(step_dt), big.mark = ","), "\n")
cat("  Participants:", length(unique(ema_dt$uid)), "\n\n")

# =============================================================================
# Process in Batches
# =============================================================================

cat("STEP 2: Processing all participants...\n\n")
process_start <- Sys.time()

all_uids <- unique(ema_dt$uid)
batch_size <- 50
n_batches <- ceiling(length(all_uids) / batch_size)

all_results <- list()
batch_times <- numeric(n_batches)

for (batch in 1:n_batches) {
  batch_start <- Sys.time()
  
  start_idx <- (batch - 1) * batch_size + 1
  end_idx <- min(batch * batch_size, length(all_uids))
  batch_uids <- all_uids[start_idx:end_idx]
  
  cat(sprintf("Batch %d/%d (UIDs %d-%d)... ", batch, n_batches, start_idx, end_idx))
  
  step_batch <- step_dt[uid %in% batch_uids]
  ema_batch <- ema_dt[uid %in% batch_uids]
  
  batch_results <- future_map_dfr(
    batch_uids, 
    function(uid_val) {
      ema_uid <- ema_batch[uid == uid_val]
      process_all_windows_for_uid(uid_val, ema_uid, step_batch)
    },
    .progress = FALSE,
    .options = furrr_options(seed = TRUE)
  )
  
  all_results[[batch]] <- batch_results
  rm(step_batch, ema_batch, batch_results)
  gc(verbose = FALSE)
  
  batch_time <- difftime(Sys.time(), batch_start, units = "secs")
  batch_times[batch] <- as.numeric(batch_time)
  
  avg_time_per_batch <- mean(batch_times[1:batch])
  remaining_batches <- n_batches - batch
  est_remaining <- avg_time_per_batch * remaining_batches
  
  cat(sprintf("%s (ETA: %s)\n", 
              format_time(as.numeric(batch_time)),
              format_time(est_remaining)))
}

process_time <- difftime(Sys.time(), process_start, units = "secs")
cat("\n✓ Processing complete in", format_time(as.numeric(process_time)), "\n\n")

# =============================================================================
# Combine and Save
# =============================================================================

cat("STEP 3: Combining results...\n")
combine_start <- Sys.time()

step_all_windows <- rbindlist(all_results)
rm(all_results)
gc(verbose = FALSE)

combine_time <- difftime(Sys.time(), combine_start, units = "secs")
cat("✓ Combined in", format_time(as.numeric(combine_time)), "\n\n")

cat("STEP 4: Saving results...\n")
save_start <- Sys.time()

output_file <- "step_all_windows_ema_based.csv"
fwrite(step_all_windows, output_file)

save_time <- difftime(Sys.time(), save_start, units = "secs")
cat("✓ Saved in", format_time(as.numeric(save_time)), "\n\n")

# =============================================================================
# Summary
# =============================================================================

total_time <- difftime(Sys.time(), total_start, units = "secs")

cat(rep("=", 70), "\n")
cat("PROCESSING COMPLETE!\n")
cat(rep("=", 70), "\n\n")

cat("TIMING BREAKDOWN:\n")
cat(sprintf("  Data loading:    %s\n", format_time(as.numeric(load_time))))
cat(sprintf("  Processing:      %s\n", format_time(as.numeric(process_time))))
cat(sprintf("  Combining:       %s\n", format_time(as.numeric(combine_time))))
cat(sprintf("  Saving:          %s\n", format_time(as.numeric(save_time))))
cat(sprintf("  --------------------------------\n"))
cat(sprintf("  TOTAL TIME:      %s\n\n", format_time(as.numeric(total_time))))

cat("OUTPUT:\n")
cat("  File:", output_file, "\n")
cat("  Rows:", format(nrow(step_all_windows), big.mark = ","), "\n")
cat("  Columns:", ncol(step_all_windows), "\n\n")

cat("FEATURE VALIDITY:\n")
feature_cols <- grep("step_(mean|mad|p90p10|entropy)", names(step_all_windows), value = TRUE)
for (col in feature_cols) {
  n_valid <- sum(!is.na(step_all_windows[[col]]))
  pct <- 100 * n_valid / nrow(step_all_windows)
  cat(sprintf("  %-25s: %6d / %6d (%.1f%%)\n", 
              col, n_valid, nrow(step_all_windows), pct))
}

cat("\nOBSERVATIONS PER WINDOW:\n")
for (suffix in c("1h", "4h", "1d", "7d")) {
  n_obs_col <- paste0("n_obs_", suffix)
  total_obs <- sum(step_all_windows[[n_obs_col]], na.rm = TRUE)
  avg_obs <- mean(step_all_windows[[n_obs_col]], na.rm = TRUE)
  cat(sprintf("  %s: total=%s, avg=%.1f\n", 
              suffix, 
              format(total_obs, big.mark = ","), 
              avg_obs))
}

cat("\nWINDOW DEFINITIONS:\n")
cat("  All windows END 4 hours before ema_day time\n")
cat("  1-hour:  (ema_day - 5h) to (ema_day - 4h)\n")
cat("  4-hour:  (ema_day - 8h) to (ema_day - 4h)\n")
cat("  1-day:   (ema_day - 28h) to (ema_day - 4h)\n")
cat("  7-day:   (ema_day - 172h) to (ema_day - 4h)\n\n")

cat(rep("=", 70), "\n")
cat("✓ ALL DONE!\n")
cat(rep("=", 70), "\n\n")

plan(sequential)
