################################################################################
# Digital Phenotype Data Processing Pipeline
# Purpose: Process EMA and digital phenotyping data for mental health analysis
# Author: Research Team
# Date: 2026
################################################################################

# =============================================================================
# SECTION 1: PACKAGE INSTALLATION AND LOADING
# =============================================================================

# List of required packages
required_packages <- c(
  "tidyverse", "jsonlite", "data.table", "lubridate", "geosphere",
  "DT", "pracma", "ggplot2", "future.apply", "furrr", "slider",
  "TSEntropies", "bit64", "vroom", "parallel", "dbscan", "hms"
)

# Install missing packages
install_if_missing <- function(packages) {
  new_packages <- packages[!(packages %in% installed.packages()[,"Package"])]
  if(length(new_packages)) install.packages(new_packages)
}

install_if_missing(required_packages)

# Load libraries
library(tidyverse)
library(jsonlite)
library(data.table)
library(lubridate)
library(geosphere)
library(DT)
library(pracma)
library(ggplot2)
library(future.apply)
library(furrr)
library(slider)
library(TSEntropies)
library(bit64)
library(vroom)
library(parallel)
library(dbscan)

# =============================================================================
# SECTION 2: EMA DATA PROCESSING
# =============================================================================

## 2.1 Load and Select Columns ----
ema <- read.csv("3._TD_cohort_EMA_FINAL.csv")

ema_selected <- ema %>%
  select(
    mlife_id, compliance, date_rep, ema_time, ema_timestamp,RespTime_12, 
    starts_with("PHQ."), starts_with("GAD.")
  ) %>%
  # Rename PHQ.# to q# and GAD.# to g#
  rename_with(~ str_replace(.x, "^PHQ\\.", "q"), starts_with("PHQ.")) %>%
  rename_with(~ str_replace(.x, "^GAD\\.", "g"), starts_with("GAD.")) %>%
  mutate(
    # Extract first number from RespTime_0
    first_timestamp = str_extract(as.character(RespTime_12), "(?<=;)[0-9.]+"),
    
    # Convert Unix timestamp to datetime
    ema_day = as.POSIXct(as.numeric(first_timestamp), 
                         origin = "1970-01-01", tz = "UTC")
  )

## 2.2 Fill Missing EMA Timepoints ----
# Parse dates
ema_selected$date_rep <- as.Date(ema_selected$date_rep, format = "%m/%d/%Y")

# Define EMA schedule (morning, afternoon, evening)
ema_levels <- c("morning", "afternoon", "evening")
slot_map <- setNames(1:3, ema_levels)

# Create complete grid for each participant
ema_selected$ema_time <- as.character(ema_selected$ema_time)
ids <- unique(ema_selected$mlife_id)

grid_list <- lapply(ids, function(id) {
  dsub <- ema_selected[ema_selected$mlife_id == id, ]
  dsub <- dsub[!is.na(dsub$date_rep) & !is.na(dsub$ema_time), ]
  
  if(nrow(dsub) == 0) return(NULL)
  
  day0 <- min(dsub$date_rep)
  day_index <- as.integer(dsub$date_rep - day0)
  slot_num <- unname(slot_map[dsub$ema_time])
  obs_index <- day_index * 3 + slot_num
  full_index <- min(obs_index):max(obs_index)
  full_day_index <- (full_index - 1) %/% 3
  full_slot_num <- (full_index - 1) %% 3 + 1
  
  data.frame(
    mlife_id = id,
    date_rep = day0 + full_day_index,
    ema_time = ema_levels[full_slot_num],
    stringsAsFactors = FALSE
  )
})

grid <- do.call(rbind, grid_list)

# Merge with observed data
ema_complete <- merge(
  grid, ema_selected,
  by = c("mlife_id", "date_rep", "ema_time"),
  all.x = TRUE
)

# Restore compliance information
comp_map <- ema_selected %>%
  filter(!is.na(compliance)) %>%
  distinct(mlife_id, compliance)

ema_complete$compliance <- ifelse(
  is.na(ema_complete$compliance),
  comp_map$compliance[match(ema_complete$mlife_id, comp_map$mlife_id)],
  ema_complete$compliance
)

## 2.3 Create DateTime Column ----
ema_complete <- ema_complete %>%
  mutate(
    ema_time = factor(ema_time, levels = ema_levels, ordered = TRUE),
    date_rep = as.Date(date_rep),
    dt = as_datetime(date_rep, tz = "UTC") + hm(as.character(ema_timestamp_hhmm))
  ) %>%
  arrange(mlife_id, date_rep, ema_time) %>%
  select(-date_rep, -ema_timestamp_hhmm) %>%
  rename(uid = mlife_id, ema_day = dt)

## 2.4 Calculate Dynamic Features ----

# Helper functions
set_entropy <- function(x, min_n = 2) {
  x <- x[!is.na(x)]
  if (length(x) < min_n) return(NA_real_)
  p <- table(x) / length(x)
  -sum(p * log2(p))
}

p90_p10 <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) < 2) return(NA_real_)
  quantile(x, 0.9, type = 7) - quantile(x, 0.1, type = 7)
}

# Calculate rolling features for all PHQ and GAD items
vars <- c(paste0("q", 1:9), "g1", "g2")

ema_features <- ema_complete %>%
  group_by(mlife_id) %>%
  mutate(
    across(
      all_of(vars),
      list(
        # 1-day window (3 timepoints)
        mean_d1 = ~slide_dbl(.x, ~mean(.x, na.rm = TRUE), 
                             .before = 2, .complete = TRUE),
        mad_d1 = ~slide_dbl(.x, ~mad(.x, na.rm = TRUE), 
                            .before = 2, .complete = TRUE),
        sentropy_d1 = ~slide_dbl(.x, ~set_entropy(.x, min_n = 2), 
                                 .before = 2, .complete = TRUE),
        p90p10_d1 = ~slide_dbl(.x, p90_p10, 
                               .before = 2, .complete = TRUE),
        
        # 7-day window (21 timepoints)
        mean_d7 = ~slide_dbl(.x, ~mean(.x, na.rm = TRUE), 
                             .before = 20, .complete = TRUE),
        mad_d7 = ~slide_dbl(.x, ~mad(.x, na.rm = TRUE), 
                            .before = 20, .complete = TRUE),
        sentropy_d7 = ~slide_dbl(.x, ~set_entropy(.x, min_n = 2), 
                                 .before = 20, .complete = TRUE),
        p90p10_d7 = ~slide_dbl(.x, p90_p10, 
                               .before = 20, .complete = TRUE)
      ),
      .names = "{.col}_{.fn}"
    )
  ) %>%
  ungroup()

# Save EMA features
write.csv(ema_features, "ema_features.csv", row.names = FALSE)

# =============================================================================
# SECTION 3: DIGITAL PHENOTYPE DATA - MERGE INDIVIDUAL FILES
# =============================================================================

## 3.1 Merge Heart Rate Files ----
main_folder <- "/Users/ahhyun/Desktop/SI_Predction"
folder_name <- "hr"
folder_path <- file.path(main_folder, folder_name)

start_time <- Sys.time()
csv_files <- list.files(folder_path, pattern = "\\.csv$", full.names = TRUE)

cat("Found", length(csv_files), "CSV files\n")

# Filter valid files
valid_files <- csv_files[sapply(csv_files, function(f) {
  test <- tryCatch(fread(f, nrows = 1, showProgress = FALSE), error = function(e) NULL)
  !is.null(test) && nrow(test) > 0
})]

cat("Processing", length(valid_files), "valid files\n")

# Read and combine files
if (length(valid_files) == 1) {
  hr <- fread(valid_files[1], showProgress = FALSE)
} else {
  num_cores <- min(14, length(valid_files))
  cl <- makeCluster(num_cores)
  clusterEvalQ(cl, library(data.table))
  
  data_list <- parLapply(cl, valid_files, function(f) fread(f, showProgress = FALSE))
  stopCluster(cl)
  
  hr <- rbindlist(data_list, fill = TRUE)
}

# Clean up
if ("V1" %in% names(hr)) hr[, V1 := NULL]
hr <- as.data.frame(hr)

cat("Combined HR data:", nrow(hr), "rows\n")
write.csv(hr, file.path(folder_path, "hr.csv"), row.names = FALSE)

# =============================================================================
# SECTION 4: LOAD AND CLEAN DIGITAL PHENOTYPE DATA
# =============================================================================

setwd("/Users/ahhyun/Desktop/SI_Predction/Final")

## 4.1 Load All Data Files ----
accelerometer <- vroom("acceleration.csv")
activity <- fread("activity.csv")
battery <- fread("battery.csv")
ema_features <- fread("ema_features.csv")
gps <- fread("gps.csv")
hr <- fread("hr.csv")
rr <- fread("rr.csv")
steps <- fread("steps.csv")
stress <- fread("stress.csv")
convo <- fread("convo.csv")
sms <- fread("smslog.csv")
call <- fread("calllog.csv")

## 4.2 Remove V1 Column from All DataFrames ----
data_frames <- c("activity", "battery", "convo", "gps", "hr", "rr", "steps", "stress")

for (df_name in data_frames) {
  df <- get(df_name)
  if ("V1" %in% names(df)) {
    df <- df %>% select(-V1)
    assign(df_name, df)
  }
}

## 4.3 Extract GPS Coordinates ----
clean <- gsub('""', '"', gps$data)
parsed <- lapply(clean, fromJSON)

get_val <- function(x, name) if (!is.null(x[[name]])) x[[name]] else NA

gps <- gps %>%
  mutate(
    alt = sapply(parsed, get_val, "ALTITUDE"),
    lon = sapply(parsed, get_val, "LONGITUDE"),
    acc = sapply(parsed, get_val, "ACCURACY"),
    lat = sapply(parsed, get_val, "LATITUDE"),
    speed = sapply(parsed, get_val, "SPEED"),
    bearing = sapply(parsed, get_val, "BEARING")
  )

## 4.4 Extract Conversation Start/End Times ----
clean <- gsub('""', '"', convo$event_data, fixed = TRUE)
parsed <- lapply(clean, fromJSON, simplifyVector = FALSE)

get_val <- function(x, name) {
  v <- x[[name]]
  if (is.null(v)) return(NA_character_)
  as.character(v)
}

convo <- convo %>%
  mutate(
    start = as.integer64(sapply(parsed, get_val, "CONVERSATION_START")),
    end = as.integer64(sapply(parsed, get_val, "CONVERSATION_END"))
  ) %>%
  rename(day = event_time) %>%
  select(-`_id`)

## 4.5 Extract Step Counts ----
steps_clean <- steps %>%
  mutate(step_count = as.numeric(str_extract(data, "(?<=\\[)[^,]+"))) %>%
  select(uid, day, step_count)

## 4.6 Remove Invalid Dates (1970 timestamps) ----
steps <- steps %>% filter(year(day) >= 2000)
stress <- stress %>% filter(year(day) >= 2000)

write.csv(stress, "stress_clean.csv", row.names = FALSE)

# =============================================================================
# SECTION 5: AGGREGATE EVENT-BASED DATA TO 10-MINUTE INTERVALS
# =============================================================================

## 5.1 Conversation Time (10-minute bins) ----
convo <- convo %>%
  mutate(
    end_time = as.POSIXct(as.numeric(end) / 1000, origin = "1970-01-01", tz = "UTC"),
    start_time = as.POSIXct(as.numeric(start) / 1000, origin = "1970-01-01", tz = "UTC")
  )

convo_per10min <- convo %>%
  group_by(uid) %>%
  reframe(Date = seq.POSIXt(min(start_time), max(end_time), by = "1 min")) %>%
  left_join(
    convo %>%
      mutate(Date = map2(start_time, end_time, ~seq.POSIXt(.x, .y, by = "1 min"))) %>%
      unnest(Date) %>%
      mutate(on_phone = 1) %>%
      select(uid, Date, on_phone),
    by = c("uid", "Date")
  ) %>%
  mutate(
    on_phone = ifelse(is.na(on_phone), 0, 1),
    time_bin = floor_date(Date, "10 minutes")
  ) %>%
  group_by(uid, time_bin) %>%
  summarise(data = sum(on_phone), .groups = "drop")

write.csv(convo_per10min, "conv_clean.csv", row.names = FALSE)

## 5.2 SMS Activity (10-minute bins) ----
sms_per10min <- sms %>%
  group_by(uid) %>%
  reframe(
    time_bin = seq.POSIXt(
      floor_date(min(time), "10 minutes"),
      floor_date(max(time), "10 minutes"),
      by = "10 min"
    )
  ) %>%
  left_join(
    sms %>%
      mutate(time_bin = floor_date(time, "10 minutes")) %>%
      group_by(uid, time_bin) %>%
      summarise(total_body_length = sum(body_len, na.rm = TRUE), .groups = "drop"),
    by = c("uid", "time_bin")
  ) %>%
  mutate(total_body_length = ifelse(is.na(total_body_length), 0, total_body_length))

write.csv(sms_per10min, "sms_clean.csv", row.names = FALSE)

## 5.3 Call Duration (10-minute bins) ----
call_per10min <- call %>%
  group_by(uid) %>%
  reframe(
    time_bin = seq.POSIXt(
      floor_date(min(time), "10 minutes"),
      floor_date(max(time), "10 minutes"),
      by = "10 min"
    )
  ) %>%
  left_join(
    call %>%
      mutate(time_bin = floor_date(time, "10 minutes")) %>%
      group_by(uid, time_bin) %>%
      summarise(total_duration = sum(duration, na.rm = TRUE), .groups = "drop"),
    by = c("uid", "time_bin")
  ) %>%
  mutate(total_duration = ifelse(is.na(total_duration), 0, total_duration))

write.csv(call_per10min, "call_clean.csv", row.names = FALSE)

## 5.4 Step Counts (10-minute bins) ----
steps_per10min <- steps_clean %>%
  group_by(uid) %>%
  reframe(
    time_bin = seq.POSIXt(
      floor_date(min(day), "10 minutes"),
      floor_date(max(day), "10 minutes"),
      by = "10 min"
    )
  ) %>%
  left_join(
    steps_clean %>%
      mutate(time_bin = floor_date(day, "10 minutes")) %>%
      group_by(uid, time_bin) %>%
      summarise(total_step_count = sum(step_count, na.rm = TRUE), .groups = "drop"),
    by = c("uid", "time_bin")
  ) %>%
  mutate(total_step_count = ifelse(is.na(total_step_count), 0, total_step_count))

write.csv(steps_per10min, "steps_clean.csv", row.names = FALSE)

# =============================================================================
# SECTION 6: GPS-BASED FEATURES
# =============================================================================

## 6.1 Calculate Movement Distance ----
gps_movement <- gps %>%
  arrange(uid, day) %>%
  group_by(uid) %>%
  mutate(
    distance_m = distHaversine(
      cbind(lag(lon), lag(lat)),
      cbind(lon, lat)
    ),
    time_diff_min = as.numeric(difftime(day, lag(day), units = "mins")),
    distance_m = ifelse(is.na(distance_m), 0, distance_m)
  ) %>%
  ungroup()

write.csv(gps_movement, "gps_clean_1.csv", row.names = FALSE)

## 6.2 Calculate Time at Home ----

# Subset to nighttime for home detection
loc_df <- gps %>%
  mutate(
    time = as.POSIXct(day),
    time_only = hms::as_hms(time)
  )

night_start <- hms::as_hms("22:00:00")
night_end <- hms::as_hms("07:00:00")

loc_night <- loc_df %>%
  filter(time_only >= night_start | time_only <= night_end)

# Function to find home locations using DBSCAN
find_home_locations <- function(df, uid_val, eps_meters = 50, min_samples = 10) {
  participant_data <- df %>% filter(uid == uid_val)
  if (nrow(participant_data) == 0) return(NULL)
  
  coords <- as.matrix(participant_data[, c("lat", "lon")])
  eps_degrees <- eps_meters / 111000
  db_result <- dbscan(coords, eps = eps_degrees, minPts = min_samples)
  
  participant_data %>%
    mutate(cluster = db_result$cluster) %>%
    filter(cluster > 0) %>%
    group_by(cluster) %>%
    summarise(
      lat = mean(lat),
      lon = mean(lon),
      n_points = n(),
      .groups = 'drop'
    ) %>%
    arrange(desc(n_points))
}

# Helper functions
calc_distance <- function(lat1, lon1, lat2, lon2) {
  distHaversine(c(lon1, lat1), c(lon2, lat2))
}

is_near_home <- function(lat, lon, home_coords, threshold_meters = 100) {
  if (is.null(home_coords) || nrow(home_coords) == 0) return(FALSE)
  distances <- sapply(1:nrow(home_coords), function(i) {
    calc_distance(lat, lon, home_coords$lat[i], home_coords$lon[i])
  })
  any(distances <= threshold_meters)
}

# Find home locations for all participants
unique_uids <- unique(loc_night$uid)
home_coords_list <- list()

cat("Finding home locations for", length(unique_uids), "participants\n")

for (i in seq_along(unique_uids)) {
  uid_val <- unique_uids[i]
  if (i %% 10 == 0) cat("Progress:", i, "/", length(unique_uids), "\n")
  
  home_coords <- find_home_locations(loc_night, uid_val)
  if (!is.null(home_coords) && nrow(home_coords) > 0) {
    home_coords_list[[uid_val]] <- home_coords
  }
}

cat("Detected homes for", length(home_coords_list), "participants\n")

# Calculate time at home per 10-minute interval
calculate_time_at_home <- function(df, uid_val, home_coords, threshold_meters = 100) {
  participant_data <- df %>%
    filter(uid == uid_val) %>%
    rowwise() %>%
    mutate(is_home = is_near_home(lat, lon, home_coords, threshold_meters)) %>%
    ungroup() %>%
    mutate(time_interval = floor_date(time, "10 minutes"))
  
  participant_data %>%
    group_by(uid, time_interval) %>%
    summarise(
      total_points = n(),
      points_at_home = sum(is_home),
      prop_at_home = mean(is_home),
      minutes_at_home = prop_at_home * 10,
      minutes_away = (1 - prop_at_home) * 10,
      .groups = 'drop'
    )
}

# Process all participants
all_time_summaries <- list()
uids_to_process <- names(home_coords_list)

for (i in seq_along(uids_to_process)) {
  uid_val <- uids_to_process[i]
  if (i %% 10 == 0) cat("Processing:", i, "/", length(uids_to_process), "\n")
  
  time_summary <- calculate_time_at_home(loc_df, uid_val, home_coords_list[[uid_val]])
  if (!is.null(time_summary)) all_time_summaries[[uid_val]] <- time_summary
}

time_at_home_df <- bind_rows(all_time_summaries)

# Save results
write.csv(time_at_home_df, "gps_clean_2.csv", row.names = FALSE)

cat("\nTime at home calculation complete!\n")
cat("Total intervals:", nrow(time_at_home_df), "\n")

# =============================================================================
# SECTION 7: ACCELEROMETER DATA PROCESSING
# =============================================================================

## 7.1 Process Individual Accelerometer Files ----
setwd("/Users/ahhyun/Desktop/SI_Predction/Before_Merged/acceleration")

process_accelerometer_file <- function(file_path) {
  if (file.size(file_path) == 0) return(NULL)
  
  dt <- tryCatch(
    fread(file_path, showProgress = FALSE, fill = TRUE),
    error = function(e) NULL
  )
  
  if (is.null(dt) || nrow(dt) == 0 || !("data" %in% names(dt))) return(NULL)
  
  # Extract x, y, z from data column
  s <- gsub("\\[|\\]", "", as.character(dt$data))
  xyz <- tstrsplit(s, ",", fixed = TRUE, fill = NA_character_)
  
  if (length(xyz) < 3) return(NULL)
  
  dt[, `:=`(
    x = as.numeric(xyz[[1]]),
    y = as.numeric(xyz[[2]]),
    z = as.numeric(xyz[[3]])
  )]
  
  dt[, magnitude := sqrt(x^2 + y^2 + z^2)]
  dt
}

# Setup parallel processing
file_list <- list.files(pattern = "\\.csv$", full.names = TRUE)
out_dir <- file.path(getwd(), "processed")
dir.create(out_dir, showWarnings = FALSE)

plan(multisession, workers = max(1, availableCores() - 1))

# Process files
handlers(global = TRUE)
out_files <- with_progress({
  p <- progressor(along = file_list)
  
  future_lapply(file_list, function(fp) {
    p(basename(fp))
    out_file <- file.path(out_dir, paste0(tools::file_path_sans_ext(basename(fp)), 
                                          "_processed.csv"))
    
    if (file.exists(out_file) && file.info(out_file)$size > 0) return(out_file)
    
    dt <- process_accelerometer_file(fp)
    if (is.null(dt)) return(NA_character_)
    
    fwrite(dt, out_file)
    out_file
  }, future.seed = TRUE)
})

cat("Processed:", sum(!is.na(out_files)), "of", length(out_files), "files\n")

## 7.2 Merge Processed Accelerometer Files ----
folder_path <- file.path(getwd(), "processed")
file_list <- list.files(folder_path, pattern = "\\.(csv|txt|tsv)$", full.names = TRUE)

cat("Merging", length(file_list), "accelerometer files\n")

# Read only needed columns
read_selected_columns <- function(file_path) {
  tryCatch(
    fread(file_path, select = c("day", "uid", "magnitude"), showProgress = FALSE),
    error = function(e) {
      warning("Error reading:", file_path)
      NULL
    }
  )
}

accel <- rbindlist(lapply(file_list, read_selected_columns), fill = TRUE)

cat("Combined accelerometer data:", nrow(accel), "rows\n")
cat("Unique UIDs:", length(unique(accel$uid)), "\n")

write.csv(accel, "accelcsv", row.names = FALSE)

# =============================================================================
# SECTION 8: DATA QUALITY CHECKS
# =============================================================================

## 8.1 Check for Invalid Dates ----
check_1970_dates <- function(df, date_col = "day") {
  df %>%
    mutate(year = year(!!sym(date_col))) %>%
    filter(year == 1970) %>%
    summarise(unique_uids = n_distinct(uid)) %>%
    pull(unique_uids)
}

cat("\nData Quality Checks:\n")
cat("UIDs with 1970 dates in stress:", check_1970_dates(stress), "\n")

## 8.2 Check for Outlier Values ----
check_outliers <- function(df, col_name, threshold) {
  df %>%
    summarise(
      total_rows = n(),
      rows_above_threshold = sum(!!sym(col_name) >= threshold, na.rm = TRUE),
      percentage = (sum(!!sym(col_name) >= threshold, na.rm = TRUE) / n()) * 100
    )
}

stress_outliers <- check_outliers(stress, "data", 4294967294)
cat("Stress outlier check:\n")
print(stress_outliers)

cat("\n=== Data Processing Complete ===\n")

# =============================================================================
# SECTION 9: Combine with ema and digital phenotype data
#.            & calculate the dynamic indices
# =============================================================================

# print(head(ema),5) #mlife_id, dt
# print(head(step),5) #v1, total_step_count, time_bin, uid
# print(head(accel),5) #data, day magnnitude
# print(head(hr),5) #data, day, uid

# step (didn't calculate entropy for 1 hr window due to the lack of data)
source("step_TIMED_ema_based.R")

source("simple_1h_test.R")
# Summary Statistics:
#   step_mean_1h       step_mad_1h       step_p90p10_1h    
# Min.   :    0.00   Min.   :    0.00   Min.   :    0.00  
# 1st Qu.:    0.00   1st Qu.:    0.00   1st Qu.:    0.00  
# Median :    3.50   Median :    0.00   Median :   10.50  
# Mean   :   26.22   Mean   :   14.99   Mean   :   53.16  
# 3rd Qu.:   28.67   3rd Qu.:   12.60   3rd Qu.:   68.00  
# Max.   :15265.00   Max.   :11315.94   Max.   :12212.00  
# NA's   :7248       NA's   :7550       NA's   :7852      
#   step_mean_4h        step_mad_4h       step_p90p10_4h    
#  Min.   :0.000e+00   Min.   :    0.00   Min.   :    0.00  
#  1st Qu.:9.167e-01   1st Qu.:    0.00   1st Qu.:    0.00  
#  Median :1.038e+01   Median :    0.00   Median :   30.70  
#  Mean   :2.631e+01   Mean   :   10.92   Mean   :   72.07  
#  3rd Qu.:3.438e+01   3rd Qu.:    0.00   3rd Qu.:   98.30  
#  Max.   :1.526e+04   Max.   :11315.94   Max.   :12212.00  
#  NA's   :7248        NA's   :7550       NA's   :7852      
# step_entropy_4h   step_mean_1d       step_mad_1d       
# Min.   :0.0000   Min.   :    0.00   Min.   :    0.000  
# 1st Qu.:0.2499   1st Qu.:   10.69   1st Qu.:    0.000  
# Median :1.2220   Median :   21.13   Median :    0.000  
# Mean   :1.5341   Mean   :   26.44   Mean   :    1.136  
# 3rd Qu.:2.5806   3rd Qu.:   36.04   3rd Qu.:    0.000  
# Max.   :4.5850   Max.   :15265.00   Max.   :11315.944  
# NA's   :7852     NA's   :7248       NA's   :7550       
#  step_p90p10_1d     step_entropy_1d  step_mean_7d     
#  Min.   :    0.00   Min.   :0.000   Min.   :    0.00  
#  1st Qu.:   31.90   1st Qu.:1.611   1st Qu.:   14.03  
#  Median :   66.70   Median :2.327   Median :   23.24  
#  Mean   :   80.68   Mean   :2.248   Mean   :   26.56  
#  3rd Qu.:  110.70   3rd Qu.:2.990   3rd Qu.:   35.47  
#  Max.   :12212.00   Max.   :5.987   Max.   :15265.00  
#  NA's   :7852       NA's   :8456    NA's   :7248      
# step_mad_7d        step_p90p10_7d     step_entropy_7d
# Min.   :0.000e+00   Min.   :    0.00   Min.   :0.000  
# 1st Qu.:0.000e+00   1st Qu.:   41.30   1st Qu.:2.121  
# Median :0.000e+00   Median :   72.30   Median :2.806  
# Mean   :2.238e-01   Mean   :   79.98   Mean   :2.720  
# 3rd Qu.:0.000e+00   3rd Qu.:  108.30   3rd Qu.:3.440  
# Max.   :1.132e+04   Max.   :12212.00   Max.   :5.987  
# NA's   :7550        NA's   :7852       NA's   :9966 

# hr
# 1h window:  75,200 × 0.05 sec  = 3,760 sec    = 1.0 hour
# 4h window:  75,200 × 0.8 sec   = 60,160 sec   = 16.7 hours
# 1d window:  75,200 × 30 sec    = 2,256,000 sec = 626.7 hours
# 7d window:  75,200 × 15 min    = 67,680,000 sec = 18,800 hours (!)

# Observations per window:
# 1h: total=6281681, avg=83.5
# 4h: total=25086445, avg=333.4
# 1d: total=149993118, avg=1993.5
# 7d: total=1005739441, avg=13367.1

# TOTAL WITHOUT PARALLELIZATION: ~19,443 hours = 810 days (!!)

# used MAX entropy point = 500
source("hr_TIMED_14cores_500pts.R")

# accel

source("accel_PARALLEL.R")

#entropy 
# sampling for efficiency
# Window   Points Available   Need for 99%        Recommended Setting
# 1h       720                All 720             Inf (use all)
# 4h       2,880              ~500–1,000          1,000
# 1d       17,280              ~2,000–3,000        3,000
# 7d       120,960             ~10,000–15,000      15,000


source("accel_ENTROPY_ONLY_OPTION_A.R")

# TO MERGE WITH EXISTING FEATURES:
library(data.table)

cat("\nMERGING ALL FEATURES...\n\n")

# Load base
features <- fread('accel_all_features.csv')
features <- unique(features, by = c('uid', 'ema_day'))  # Remove duplicates
cat("1. Loaded accel features:", nrow(features), "rows\n")

# Merge accel entropy
ent_1h <- fread('accel_entropy_1h.csv')
ent_1h <- unique(ent_1h, by = c('uid', 'ema_day'))
features <- merge(features, ent_1h, by = c('uid', 'ema_day'), all.x = TRUE)

ent_4h <- fread('accel_entropy_4h.csv')
ent_4h <- unique(ent_4h, by = c('uid', 'ema_day'))
features <- merge(features, ent_4h, by = c('uid', 'ema_day'), all.x = TRUE)

ent_1d <- fread('accel_entropy_1d.csv')
ent_1d <- unique(ent_1d, by = c('uid', 'ema_day'))
features <- merge(features, ent_1d, by = c('uid', 'ema_day'), all.x = TRUE)

ent_7d <- fread('accel_entropy_7d.csv')
ent_7d <- unique(ent_7d, by = c('uid', 'ema_day'))
features <- merge(features, ent_7d, by = c('uid', 'ema_day'), all.x = TRUE)
cat("2. Merged accel entropy:", ncol(features), "columns\n")

# Merge EMA
ema <- fread('ema_features.csv')
setnames(ema, "mlife_id", "uid")
ema[, ema_day := as.POSIXct(ema_day, tz = "UTC")]
ema <- unique(ema, by = c('uid', 'ema_day'))
features <- merge(features, ema, by = c('uid', 'ema_day'), all.x = TRUE)
cat("3. Merged EMA:", ncol(features), "columns\n")

# Merge HR
hr <- fread('hr_all_windows_all_metrics.csv')
setnames(hr, "ema_time", "ema_day")
hr[, ema_day := as.POSIXct(ema_day, tz = "UTC")]
hr <- unique(hr, by = c('uid', 'ema_day'))
features <- merge(features, hr, by = c('uid', 'ema_day'), all.x = TRUE)
cat("4. Merged HR:", ncol(features), "columns\n")

# Merge steps
steps <- fread('step_all_windows_ema_based.csv')
if ("mlife_id" %in% names(steps)) setnames(steps, "mlife_id", "uid")
steps[, ema_day := as.POSIXct(ema_day, tz = "UTC")]
steps <- unique(steps, by = c('uid', 'ema_day'))
features <- merge(features, steps, by = c('uid', 'ema_day'), all.x = TRUE)
cat("5. Merged steps:", ncol(features), "columns\n")

# Save
fwrite(features, 'all_features_combined.csv')
cat("\n✓ DONE! Saved to: all_features_combined.csv\n")
cat("  Rows:", nrow(features), "\n")
cat("  Columns:", ncol(features), "\n\n")
