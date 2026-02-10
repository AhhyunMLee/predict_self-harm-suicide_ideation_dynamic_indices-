# ============================================================================
# ACCELEROMETER Data: All 4 Metrics for 1h, 4h, 1d, 7d Windows
# ULTRA-OPTIMIZED for 20GB Data - Single-pass + 4-hour lag + 12 Cores
# ============================================================================

library(tidyverse)
library(lubridate)
library(data.table)
library(pracma)
library(future)
library(furrr)

# ============================================================================
# CRITICAL SETTINGS FOR LARGE DATA
# ============================================================================

cat("Configuring for accelerometer data...\n")

# Set up parallel processing
n_cores <- 14  # Using 14 cores as requested
plan(multisession, workers = n_cores)
options(future.globals.maxSize = 4000 * 1024^2)  # 4GB max per worker

# data.table settings for speed
setDTthreads(n_cores)

cat("✓ Optimized parallel processing configured (", n_cores, " cores)\n")
cat("  Memory limit per worker: 4GB\n")

# ============================================================================
# Helper Functions
# ============================================================================

# SPEED TOGGLE: Set to FALSE to skip entropy (10x faster but no entropy features)
CALCULATE_ENTROPY <- TRUE  # Set to TRUE for full entropy calculation
MAX_ENTROPY_POINTS <- 500  # Using 500 points as requested

calculate_sample_entropy <- function(x, m = 2, r = 0.2) {
  if (!CALCULATE_ENTROPY) return(NA_real_)  # Skip if disabled
  if (length(x) < 10 || all(is.na(x))) return(NA_real_)
  x_clean <- x[!is.na(x)]
  if (length(x_clean) < 10) return(NA_real_)
  
  # SPEED OPTIMIZATION: Subsample large datasets to MAX_ENTROPY_POINTS
  if (length(x_clean) > MAX_ENTROPY_POINTS) {
    indices <- seq(1, length(x_clean), length.out = MAX_ENTROPY_POINTS)
    x_clean <- x_clean[round(indices)]
  }
  
  tolerance <- r * sd(x_clean, na.rm = TRUE)
  if (tolerance == 0) return(NA_real_)
  tryCatch({
    pracma::sample_entropy(x_clean, edim = m, r = tolerance)
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

# ============================================================================
# CHUNKED DATA LOADING - KEY FOR 20GB DATA
# ============================================================================

load_accel_by_chunks <- function(accel_file, chunk_size = 5000000) {
  
  cat("\nLoading accelerometer data in chunks...\n")
  
  # Check if file exists
  if (!file.exists(accel_file)) {
    stop(sprintf("ERROR: File '%s' not found!\nPlease check the filename and path.", accel_file))
  }
  
  # First, read header to get actual column names
  cat("Reading file header to detect column names...\n")
  header_sample <- fread(accel_file, nrows = 1)
  actual_cols <- names(header_sample)
  cat("Available columns:", paste(actual_cols, collapse = ", "), "\n")
  
  # Map to expected column names
  uid_col <- NULL
  time_col <- NULL
  value_col <- NULL
  
  # Try to find UID column
  uid_patterns <- c("uid", "mlife_id", "participant_id", "id", "user_id")
  for (pattern in uid_patterns) {
    matches <- grep(pattern, actual_cols, ignore.case = TRUE, value = TRUE)
    if (length(matches) > 0) {
      uid_col <- matches[1]
      break
    }
  }
  
  # Try to find time column
  time_patterns <- c("day", "date", "time", "timestamp", "datetime", "dat")
  for (pattern in time_patterns) {
    matches <- grep(pattern, actual_cols, ignore.case = TRUE, value = TRUE)
    if (length(matches) > 0) {
      time_col <- matches[1]
      break
    }
  }
  
  # Try to find value column (magnitude)
  value_patterns <- c("magnitude", "mag", "value", "accel", "acceleration")
  for (pattern in value_patterns) {
    matches <- grep(pattern, actual_cols, ignore.case = TRUE, value = TRUE)
    if (length(matches) > 0) {
      value_col <- matches[1]
      break
    }
  }
  
  # Check if all columns found
  if (is.null(uid_col)) {
    stop(sprintf("ERROR: Could not find UID column. Available columns: %s\nPlease specify uid_col manually.", 
                 paste(actual_cols, collapse = ", ")))
  }
  if (is.null(time_col)) {
    stop(sprintf("ERROR: Could not find time/date column. Available columns: %s\nPlease specify time_col manually.", 
                 paste(actual_cols, collapse = ", ")))
  }
  if (is.null(value_col)) {
    stop(sprintf("ERROR: Could not find magnitude/value column. Available columns: %s\nPlease specify value_col manually.", 
                 paste(actual_cols, collapse = ", ")))
  }
  
  cat("✓ Using columns:\n")
  cat("  UID:", uid_col, "\n")
  cat("  Time:", time_col, "\n")
  cat("  Value:", value_col, "\n\n")
  
  # Get total rows (with better error handling)
  wc_cmd <- sprintf("wc -l < '%s'", accel_file)
  total_rows_result <- tryCatch({
    as.numeric(system(wc_cmd, intern = TRUE)) - 1
  }, error = function(e) {
    # Fallback: try to count rows using R (slower but works)
    cat("  Using R to count rows (this may take a moment)...\n")
    nrow(fread(accel_file, select = 1, showProgress = FALSE))
  })
  
  if (is.na(total_rows_result) || total_rows_result <= 0) {
    stop("ERROR: Could not determine file size. File may be empty or corrupted.")
  }
  
  total_rows <- total_rows_result
  cat("Total rows:", format(total_rows, big.mark = ","), "\n")
  
  # Calculate chunks
  n_chunks <- ceiling(total_rows / chunk_size)
  cat("Processing in", n_chunks, "chunks of", format(chunk_size, big.mark = ","), "rows\n")
  
  accel_list <- list()
  
  for (i in 1:n_chunks) {
    skip_rows <- (i - 1) * chunk_size
    
    cat(sprintf("  Chunk %d/%d...", i, n_chunks))
    
    chunk <- fread(
      accel_file,
      skip = skip_rows + 1,  # +1 to skip header on first read
      nrows = chunk_size,
      select = c(uid_col, time_col, value_col),  # Use actual column names
      showProgress = FALSE
    )
    
    # Rename to standard names
    setnames(chunk, 
             old = c(uid_col, time_col, value_col),
             new = c("uid", "day", "magnitude"))
    
    # Convert day to POSIXct
    chunk[, day := as.POSIXct(day, tz = "UTC")]
    
    # Set keys for fast lookup
    setkey(chunk, uid, day)
    
    accel_list[[i]] <- chunk
    
    cat(sprintf(" %s rows loaded\n", format(nrow(chunk), big.mark = ",")))
    
    # Clear memory every 5 chunks
    if (i %% 5 == 0) gc()
  }
  
  cat("\nCombining chunks...\n")
  accel_dt <- rbindlist(accel_list)
  rm(accel_list)
  gc()
  
  cat("✓ Total accelerometer data loaded:", format(nrow(accel_dt), big.mark = ","), "rows\n")
  
  return(accel_dt)
}

# ============================================================================
# OPTIMIZED: Process ALL windows at once per UID (SINGLE-PASS)
# ============================================================================

process_all_windows_for_uid <- function(uid_val, ema_uid, accel_batch) {
  
  accel_uid <- accel_batch[uid == uid_val]
  
  if (nrow(accel_uid) == 0) {
    # Return empty results for all windows
    empty_result <- data.table(
      uid = uid_val,
      ema_time = ema_uid$ema_time
    )
    
    # Add columns for all windows
    for (suffix in c("1h", "4h", "1d", "7d")) {
      empty_result[, paste0("mean_", suffix, "_accel") := NA_real_]
      empty_result[, paste0("MAD_", suffix, "_accel") := NA_real_]
      empty_result[, paste0("P90P10_pos_", suffix, "_accel") := NA_real_]
      empty_result[, paste0("entropy_", suffix, "_accel") := NA_real_]
      empty_result[, paste0("n_obs_", suffix) := 0L]
    }
    
    return(empty_result)
  }
  
  # Pre-compute year-agnostic times ONCE
  accel_times <- as.POSIXct(accel_uid$day, tz = "UTC")
  year(accel_times) <- 2021
  accel_values <- accel_uid$magnitude
  
  # Define windows (all end 4 hours before current time)
  # 1h: 5 hours before to 4 hours before
  # 4h: 8 hours before to 4 hours before
  # 1d: 28 hours before to 4 hours before
  # 7d: 172 hours before to 4 hours before
  
  windows <- list(
    "1h" = hours(5),
    "4h" = hours(8),
    "1d" = hours(28),
    "7d" = hours(172)
  )
  
  # Process all EMA times for this participant - ALL WINDOWS AT ONCE
  results <- ema_uid[, {
    t_ema <- as.POSIXct(ema_time, tz = "UTC")
    year(t_ema) <- 2021
    window_end <- t_ema - hours(4)  # All windows end 4 hours before
    
    # Calculate all windows at once
    result_row <- list()
    
    for (suffix in names(windows)) {
      window_start <- t_ema - windows[[suffix]]
      
      # Extract accel values in this window
      idx <- accel_times > window_start & accel_times <= window_end
      accel_in_window <- accel_values[idx]
      
      # Calculate metrics
      metrics <- calculate_all_metrics(accel_in_window)
      
      # Store results
      result_row[[paste0("mean_", suffix, "_accel")]] <- metrics$mean
      result_row[[paste0("MAD_", suffix, "_accel")]] <- metrics$MAD
      result_row[[paste0("P90P10_pos_", suffix, "_accel")]] <- metrics$P90P10
      result_row[[paste0("entropy_", suffix, "_accel")]] <- metrics$entropy
      result_row[[paste0("n_obs_", suffix)]] <- sum(idx)
    }
    
    result_row
  }, by = .(uid, ema_time)]
  
  return(results)
}

# ============================================================================
# Memory-Efficient Batch Processing with Single-Pass Windows
# ============================================================================

process_all_windows_batched <- function(ema_dt, accel_dt, batch_size = 50) {
  
  cat("\n=== Processing ALL windows in SINGLE PASS (memory-efficient batching) ===\n")
  
  all_uids <- unique(ema_dt$uid)
  n_uids <- length(all_uids)
  
  cat("Total participants:", n_uids, "\n")
  cat("Batch size:", batch_size, "UIDs per batch\n")
  cat("Using", n_cores, "cores + chunked processing\n")
  cat("Max entropy points:", MAX_ENTROPY_POINTS, "\n")
  cat("\nWindow definitions (all end 4 hours before current time):\n")
  cat("  1-hour:  5 hours before → 4 hours before\n")
  cat("  4-hour:  8 hours before → 4 hours before\n")
  cat("  1-day:   28 hours before → 4 hours before\n")
  cat("  7-day:   172 hours before → 4 hours before\n\n")
  
  # Process UIDs in batches to avoid memory overflow
  n_batches <- ceiling(n_uids / batch_size)
  
  all_results <- list()
  batch_times <- numeric(n_batches)
  
  format_time <- function(seconds) {
    if (seconds < 60) sprintf("%.1f sec", seconds)
    else if (seconds < 3600) sprintf("%.1f min", seconds / 60)
    else sprintf("%.1f hrs", seconds / 3600)
  }
  
  for (batch in 1:n_batches) {
    batch_start <- Sys.time()
    
    start_idx <- (batch - 1) * batch_size + 1
    end_idx <- min(batch * batch_size, n_uids)
    batch_uids <- all_uids[start_idx:end_idx]
    
    cat(sprintf("Batch %d/%d (UIDs %d-%d)... ", batch, n_batches, start_idx, end_idx))
    
    # Subset accel data for this batch only
    accel_batch <- accel_dt[uid %in% batch_uids]
    
    # Process this batch in parallel - ALL WINDOWS AT ONCE PER UID
    batch_results <- future_map_dfr(batch_uids, function(uid_val) {
      ema_uid <- ema_dt[uid == uid_val]
      process_all_windows_for_uid(uid_val, ema_uid, accel_batch)
    }, .progress = FALSE, .options = furrr_options(seed = TRUE))
    
    all_results[[batch]] <- batch_results
    
    # Clear batch from memory
    rm(accel_batch, batch_results)
    gc()
    
    batch_time <- difftime(Sys.time(), batch_start, units = "secs")
    batch_times[batch] <- as.numeric(batch_time)
    
    # Calculate ETA
    avg_time_per_batch <- mean(batch_times[1:batch])
    remaining_batches <- n_batches - batch
    est_remaining <- avg_time_per_batch * remaining_batches
    
    cat(sprintf("%s (ETA: %s)\n", 
                format_time(as.numeric(batch_time)),
                format_time(est_remaining)))
  }
  
  cat("\nCombining all batches...\n")
  final_results <- rbindlist(all_results)
  rm(all_results)
  gc()
  
  cat("✓ All windows processed in single pass!\n")
  return(final_results)
}

# ============================================================================
# Step 1: Load EMA data (small, loads fully)
# ============================================================================

cat("\nStep 1: Loading EMA data...\n")

# Load EMA data
ema <- fread("ema_features.csv")

# Prepare EMA data
ema_dt <- as.data.table(ema)[, .(uid = mlife_id, ema_time = as.POSIXct(ema_day, tz = "UTC"))]
setkey(ema_dt, uid, ema_time)

cat("✓ EMA data prepared\n")
cat("  EMA rows:", format(nrow(ema_dt), big.mark = ","), "\n")
cat("  Unique participants:", uniqueN(ema_dt$uid), "\n")

# ============================================================================
# Step 2: Load Accelerometer data (20GB - using chunking)
# ============================================================================

cat("\n", rep("=", 60), "\n")
cat("LOADING ACCELEROMETER DATA\n")
cat(rep("=", 60), "\n")

# Check what accelerometer files are available
cat("\nLooking for accelerometer data files...\n")
possible_files <- c(
  "merged_accelerometer.csv",
  "accel_clean.csv",
  "accelerometer.csv",
  "accel.csv"
)

accel_file <- NULL
for (f in possible_files) {
  if (file.exists(f)) {
    accel_file <- f
    cat("✓ Found:", f, "\n")
    break
  }
}

if (is.null(accel_file)) {
  cat("\nERROR: No accelerometer file found!\n")
  cat("Please specify your accelerometer filename.\n")
  cat("Example: accel_file <- 'your_accel_file.csv'\n\n")
  stop("Accelerometer file not found. Please update the filename.")
}

cat("\nUsing accelerometer file:", accel_file, "\n")

# OPTION 1: Load in chunks (recommended for large data)
accel_dt <- load_accel_by_chunks(accel_file, chunk_size = 5000000)

# OPTION 2: If you have enough RAM (32GB+), load directly (uncomment to use)
# cat("\nLoading accelerometer data directly (requires sufficient RAM)...\n")
# accel_dt <- fread(accel_file, 
#                   select = c("uid", "day", "magnitude"),
#                   showProgress = TRUE)
# accel_dt[, day := as.POSIXct(day, tz = "UTC")]
# setkey(accel_dt, uid, day)
# cat("✓ Loaded:", format(nrow(accel_dt), big.mark = ","), "rows\n")

cat("\n✓ Accelerometer data ready\n")
cat("  Accel rows:", format(nrow(accel_dt), big.mark = ","), "\n")
cat("  Unique participants:", uniqueN(accel_dt$uid), "\n")

# ============================================================================
# Step 3: Process ALL windows in SINGLE PASS (OPTIMIZED)
# ============================================================================

cat("\n", rep("=", 60), "\n")
cat("Processing ALL 4 Windows for ACCELEROMETER Data\n")
cat("SINGLE-PASS OPTIMIZATION: All windows calculated together\n")
cat(rep("=", 60), "\n")

# Start total timer
total_start <- Sys.time()

# Process all windows at once (3-4x faster than separate passes)
accel_all_windows <- process_all_windows_batched(ema_dt, accel_dt, batch_size = 50)

# Calculate total processing time
total_time <- difftime(Sys.time(), total_start, units = "secs")

# ============================================================================
# Step 4: Save results
# ============================================================================

cat("\n=== Saving results ===\n")

# Save combined file with all windows
output_file <- "accel_all_windows_all_metrics.csv"
fwrite(accel_all_windows, output_file)
cat("✓ Saved complete dataset to:", output_file, "\n")
cat("  Rows:", format(nrow(accel_all_windows), big.mark = ","), "\n")
cat("  Columns:", ncol(accel_all_windows), "\n")

# Also save individual window files for compatibility
for (suffix in c("1h", "4h", "1d", "7d")) {
  window_cols <- c("uid", "ema_time", 
                   paste0("mean_", suffix, "_accel"),
                   paste0("MAD_", suffix, "_accel"),
                   paste0("P90P10_pos_", suffix, "_accel"),
                   paste0("entropy_", suffix, "_accel"),
                   paste0("n_obs_", suffix))
  
  window_file <- paste0("accel_all_4metrics_", suffix, ".csv")
  fwrite(accel_all_windows[, ..window_cols], window_file)
  cat("✓ Saved:", window_file, "\n")
}

# ============================================================================
# Step 5: Summary statistics
# ============================================================================

cat("\n=== Summary Statistics ===\n")

# Calculate completeness for each feature
feature_cols <- grep("^(mean|MAD|P90P10|entropy)_", names(accel_all_windows), value = TRUE)

completeness <- accel_all_windows[, lapply(.SD, function(x) {
  sprintf("%.1f%%", 100 * mean(!is.na(x)))
}), .SDcols = feature_cols]

cat("\nFeature Completeness:\n")
print(t(completeness))

# Show sample data
cat("\nSample from combined dataset (first 10 rows):\n")
print(head(as.data.frame(accel_all_windows), 10))

# Count observations per window
cat("\nObservations per window:\n")
for (suffix in c("1h", "4h", "1d", "7d")) {
  n_obs_col <- paste0("n_obs_", suffix)
  total_obs <- sum(accel_all_windows[[n_obs_col]], na.rm = TRUE)
  avg_obs <- mean(accel_all_windows[[n_obs_col]], na.rm = TRUE)
  cat(sprintf("  %s: total=%s, avg=%.1f\n", 
              suffix, 
              format(total_obs, big.mark = ","), 
              avg_obs))
}

# ============================================================================
# Step 6: Completion message
# ============================================================================

cat("\n", rep("=", 60), "\n")
cat("ALL ACCELEROMETER PROCESSING COMPLETE!\n")
cat(rep("=", 60), "\n")

# Format time function
format_time <- function(seconds) {
  if (seconds < 60) sprintf("%.1f seconds", seconds)
  else if (seconds < 3600) sprintf("%.1f minutes", seconds / 60)
  else sprintf("%.1f hours", seconds / 3600)
}

cat("\n*** TOTAL PROCESSING TIME:", format_time(as.numeric(total_time)), "***\n\n")

cat("Files created:\n")
cat("  1. accel_all_windows_all_metrics.csv (COMBINED - all windows)\n")
cat("  2. accel_all_4metrics_1h.csv\n")
cat("  3. accel_all_4metrics_4h.csv\n")
cat("  4. accel_all_4metrics_1d.csv\n")
cat("  5. accel_all_4metrics_7d.csv\n")

cat("\nEach file contains:\n")
cat("  - uid, ema_time\n")
cat("  - mean_{window}_accel\n")
cat("  - MAD_{window}_accel\n")
cat("  - P90P10_pos_{window}_accel\n")
cat("  - entropy_{window}_accel\n")
cat("  - n_obs_{window}\n")

cat("\nWindow Definitions:\n")
cat("  1-hour:  (ema_time - 5h) to (ema_time - 4h)\n")
cat("  4-hour:  (ema_time - 8h) to (ema_time - 4h)\n")
cat("  1-day:   (ema_time - 28h) to (ema_time - 4h)\n")
cat("  7-day:   (ema_time - 172h) to (ema_time - 4h)\n")

cat("\nConfiguration used:\n")
cat("  ✓ Cores:", n_cores, "\n")
cat("  ✓ Max entropy points:", MAX_ENTROPY_POINTS, "\n")
cat("  ✓ SINGLE-PASS processing (all windows calculated together)\n")
cat("  ✓ 4-hour lag applied to all windows\n")

# Close parallel
plan(sequential)

cat("\n✓ ALL DONE!\n")
