# ============================================================================
# HR Data: TIMED VERSION - 14 cores, MAX=500 entropy points
# ============================================================================

library(tidyverse)
library(lubridate)
library(data.table)
library(pracma)
library(future)
library(furrr)

# ============================================================================
# CONFIGURATION
# ============================================================================

CALCULATE_ENTROPY <- TRUE
MAX_ENTROPY_POINTS <- 500  # Your requested setting
N_CORES <- 14              # Your requested cores

cat("\n")
cat(rep("=", 70), "\n")
cat("HR PROCESSING - TIMED VERSION\n")
cat(rep("=", 70), "\n\n")
cat("Configuration:\n")
cat("  Cores:                  ", N_CORES, "\n")
cat("  Max entropy points:     ", MAX_ENTROPY_POINTS, "\n")
cat("  Calculate entropy:      ", CALCULATE_ENTROPY, "\n\n")

# ============================================================================
# Timer utility
# ============================================================================

format_time <- function(seconds) {
  if (seconds < 60) {
    sprintf("%.1f seconds", seconds)
  } else if (seconds < 3600) {
    sprintf("%.1f minutes", seconds / 60)
  } else {
    sprintf("%.1f hours", seconds / 3600)
  }
}

# ============================================================================
# Setup
# ============================================================================

total_start <- Sys.time()

cat("Setting up parallel processing...\n")
plan(multisession, workers = N_CORES)
options(future.globals.maxSize = 4000 * 1024^2)
cat("✓ Using", N_CORES, "cores\n\n")

# ============================================================================
# Helper Functions
# ============================================================================

calculate_sample_entropy <- function(x, m = 2, r = 0.2) {
  if (!CALCULATE_ENTROPY) return(NA_real_)
  if (length(x) < 10 || all(is.na(x))) return(NA_real_)
  x_clean <- x[!is.na(x)]
  if (length(x_clean) < 10) return(NA_real_)
  
  # Subsample for speed
  if (length(x_clean) > MAX_ENTROPY_POINTS) {
    indices <- seq(1, length(x_clean), length.out = MAX_ENTROPY_POINTS)
    x_clean <- x_clean[round(indices)]
  }
  
  sd_val <- sd(x_clean, na.rm = TRUE)
  if (is.na(sd_val) || sd_val == 0) return(NA_real_)
  
  tolerance <- r * sd_val
  if (tolerance == 0 || is.na(tolerance)) return(NA_real_)
  
  tryCatch({
    ent <- pracma::sample_entropy(x_clean, edim = m, r = tolerance)
    if (is.null(ent) || is.na(ent) || is.nan(ent) || is.infinite(ent)) {
      NA_real_
    } else {
      as.numeric(ent)
    }
  }, error = function(e) NA_real_)
}

calculate_all_metrics <- function(x) {
  x_clean <- x[!is.na(x)]
  n <- length(x_clean)
  
  if (n == 0) {
    return(list(mean = NA_real_, MAD = NA_real_, P90P10 = NA_real_, entropy = NA_real_))
  }
  
  m <- mean(x_clean)
  mad_val <- mad(x_clean, center = median(x_clean), constant = 1)
  
  p90p10 <- if (n >= 3) {
    p10 <- quantile(x_clean, 0.10, names = FALSE)
    p90 <- quantile(x_clean, 0.90, names = FALSE)
    cur <- x_clean[n]
    if (is.na(p90) || is.na(p10) || is.na(cur) || p90 == p10) {
      NA_real_
    } else {
      (cur - p10) / (p90 - p10)
    }
  } else {
    NA_real_
  }
  
  entropy <- calculate_sample_entropy(x_clean)
  
  list(mean = m, MAD = mad_val, P90P10 = p90p10, entropy = entropy)
}

process_all_windows_for_uid <- function(uid_val, ema_uid, hr_batch) {
  
  hr_uid <- hr_batch[uid == uid_val]
  
  if (nrow(hr_uid) == 0) {
    empty_result <- data.table(uid = uid_val, ema_time = ema_uid$ema_time)
    for (suffix in c("1h", "4h", "1d", "7d")) {
      empty_result[, paste0("mean_", suffix, "_hr") := NA_real_]
      empty_result[, paste0("MAD_", suffix, "_hr") := NA_real_]
      empty_result[, paste0("P90P10_pos_", suffix, "_hr") := NA_real_]
      empty_result[, paste0("entropy_", suffix, "_hr") := NA_real_]
      empty_result[, paste0("n_obs_", suffix) := 0L]
    }
    return(empty_result)
  }
  
  hr_times <- as.POSIXct(hr_uid$time_sensor, tz = "UTC")
  year(hr_times) <- 2021
  hr_values <- hr_uid$value
  
  windows <- list(
    "1h" = hours(5),
    "4h" = hours(8),
    "1d" = hours(28),
    "7d" = hours(172)
  )
  
  results <- ema_uid[, {
    t_ema <- as.POSIXct(ema_time, tz = "UTC")
    year(t_ema) <- 2021
    window_end <- t_ema - hours(4)
    
    result_row <- list()
    
    for (suffix in names(windows)) {
      window_start <- t_ema - windows[[suffix]]
      idx <- hr_times > window_start & hr_times <= window_end
      hr_in_window <- hr_values[idx]
      
      metrics <- calculate_all_metrics(hr_in_window)
      
      result_row[[paste0("mean_", suffix, "_hr")]] <- metrics$mean
      result_row[[paste0("MAD_", suffix, "_hr")]] <- metrics$MAD
      result_row[[paste0("P90P10_pos_", suffix, "_hr")]] <- metrics$P90P10
      result_row[[paste0("entropy_", suffix, "_hr")]] <- metrics$entropy
      result_row[[paste0("n_obs_", suffix)]] <- sum(idx)
    }
    
    result_row
  }, by = .(uid, ema_time)]
  
  return(results)
}

# ============================================================================
# Load Data
# ============================================================================

cat("STEP 1: Loading data...\n")
load_start <- Sys.time()

ema <- fread("ema_features.csv")
hr <- fread("hr_clean.csv")

ema_dt <- as.data.table(ema)[, .(uid = mlife_id, ema_time = as.POSIXct(ema_day, tz = "UTC"))]
setkey(ema_dt, uid, ema_time)

if ("data" %in% names(hr) && "day" %in% names(hr)) {
  hr_dt <- as.data.table(hr)[, .(uid, time_sensor = as.POSIXct(day, tz = "UTC"), value = as.numeric(data))]
} else {
  stop("ERROR: Need 'data' and 'day' columns")
}
setkey(hr_dt, uid, time_sensor)

load_time <- difftime(Sys.time(), load_start, units = "secs")
cat("✓ Data loaded in", format_time(as.numeric(load_time)), "\n")
cat("  EMA rows:", format(nrow(ema_dt), big.mark = ","), "\n")
cat("  HR rows:", format(nrow(hr_dt), big.mark = ","), "\n")
cat("  Participants:", length(unique(ema_dt$uid)), "\n\n")

# ============================================================================
# Process in Batches
# ============================================================================

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
  
  hr_batch <- hr_dt[uid %in% batch_uids]
  ema_batch <- ema_dt[uid %in% batch_uids]
  
  batch_results <- future_map_dfr(
    batch_uids, 
    function(uid_val) {
      ema_uid <- ema_batch[uid == uid_val]
      process_all_windows_for_uid(uid_val, ema_uid, hr_batch)
    },
    .progress = FALSE,
    .options = furrr_options(seed = TRUE)
  )
  
  all_results[[batch]] <- batch_results
  rm(hr_batch, ema_batch, batch_results)
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

# ============================================================================
# Combine and Save
# ============================================================================

cat("STEP 3: Combining results...\n")
combine_start <- Sys.time()

hr_all_windows <- rbindlist(all_results)
rm(all_results)
gc(verbose = FALSE)

combine_time <- difftime(Sys.time(), combine_start, units = "secs")
cat("✓ Combined in", format_time(as.numeric(combine_time)), "\n\n")

cat("STEP 4: Saving results...\n")
save_start <- Sys.time()

output_file <- "hr_all_windows_all_metrics.csv"
fwrite(hr_all_windows, output_file)

save_time <- difftime(Sys.time(), save_start, units = "secs")
cat("✓ Saved in", format_time(as.numeric(save_time)), "\n\n")

# ============================================================================
# Summary
# ============================================================================

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
cat("  Rows:", format(nrow(hr_all_windows), big.mark = ","), "\n")
cat("  Columns:", ncol(hr_all_windows), "\n\n")

cat("ENTROPY VALIDITY:\n")
entropy_cols <- grep("entropy", names(hr_all_windows), value = TRUE)
for (col in entropy_cols) {
  n_valid <- sum(!is.na(hr_all_windows[[col]]))
  pct <- 100 * n_valid / nrow(hr_all_windows)
  cat(sprintf("  %-20s: %6d / %6d (%.1f%%)\n", 
              col, n_valid, nrow(hr_all_windows), pct))
}

cat("\nOBSERVATIONS PER WINDOW:\n")
for (suffix in c("1h", "4h", "1d", "7d")) {
  n_obs_col <- paste0("n_obs_", suffix)
  total_obs <- sum(hr_all_windows[[n_obs_col]], na.rm = TRUE)
  avg_obs <- mean(hr_all_windows[[n_obs_col]], na.rm = TRUE)
  cat(sprintf("  %s: total=%s, avg=%.1f\n", 
              suffix, 
              format(total_obs, big.mark = ","), 
              avg_obs))
}

cat("\n")
cat(rep("=", 70), "\n")
cat("✓ ALL DONE!\n")
cat(rep("=", 70), "\n\n")

plan(sequential)
