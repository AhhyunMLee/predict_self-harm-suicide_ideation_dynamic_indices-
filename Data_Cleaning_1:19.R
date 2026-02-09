#download packages
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
library(pracma) 
library(dplyr) 
library(tidyr)
library(lubridate)
library(dplyr)
library(slider)
library(lubridate)
library(TSEntropies)
library(data.table)
library(jsonlite)  
library(bit64)
library(vroom)


###### Data Cleaning ######

###### 1. EMA #### ######

# load the file
ema <- read.csv("3._TD_cohort_EMA_FINAL.csv")
View(ema)

# select the columns that will be used for the study
# change the column name (PHQ.# -> q#, gad.# -> g#)
# make a new column (ema_timestamp_hhmm: right version of ema_timestamp)
# some of the data is in a wrong format (ema_timestamp: some of them has date+time of a day)
ema_selected <- ema %>%
  select(
    mlife_id, compliance, date_rep, ema_time, ema_timestamp,
    `PHQ.1`,`PHQ.2`,`PHQ.3`,`PHQ.4`,`PHQ.5`,`PHQ.6`,`PHQ.7`,`PHQ.8`,`PHQ.9`,
    `GAD.1`, `GAD.2`
  ) %>%
  rename_with(~ str_replace(.x, "^PHQ\\.", "q"), starts_with("PHQ.")) %>%
  rename_with(~ str_replace(.x, "^GAD\\.", "g"), starts_with("GAD.")) %>%
  mutate(
    ema_ts_chr = as.character(unlist(ema_timestamp)),
    ema_hm = str_extract(ema_ts_chr, "\\b\\d{1,2}:\\d{1,2}\\b"),
    ema_timestamp_hhmm = if_else(
      !is.na(ema_hm),
      sprintf(
        "%02d:%02d",
        as.integer(sub(":.*", "", ema_hm)),
        as.integer(sub(".*:", "", ema_hm))
      ),
      NA_character_
    )
  ) %>%
  select(-ema_ts_chr, -ema_hm) 
View(ema_selected) # need to check out the new column (ema_timestamp_hhmm) and compare with the wrong column (ema_timestamp)
summary(ema_selected)

# delete the wrong column (ema_timestamp)
ema_selected <- ema_selected %>% select(-ema_timestamp)
View(ema_selected)

# add NA rows for missing data
# parse date (m/d/Y)
ema_selected$date_rep <- as.Date(ema_selected$date_rep, format = "%m/%d/%Y")

# define EMA order (M → A → E)
ema_levels <- c("morning", "afternoon", "evening")
slot_map   <- setNames(1:3, ema_levels)

ema_selected$ema_time <- as.character(ema_selected$ema_time)

ids <- unique(ema_selected$mlife_id)

grid_list <- lapply(ids, function(id) {
  
  dsub <- ema_selected[ema_selected$mlife_id == id, ]
  dsub <- dsub[!is.na(dsub$date_rep) & !is.na(dsub$ema_time), ]
  
  day0 <- min(dsub$date_rep)
  
  day_index <- as.integer(dsub$date_rep - day0)
  slot_num  <- unname(slot_map[dsub$ema_time])
  
  obs_index <- day_index * 3 + slot_num
  
  full_index <- min(obs_index):max(obs_index)
  
  full_day_index <- (full_index - 1) %/% 3
  full_slot_num  <- (full_index - 1) %% 3 + 1
  
  data.frame(
    mlife_id = id,
    date_rep = day0 + full_day_index,
    ema_time = ema_levels[full_slot_num],
    stringsAsFactors = FALSE
  )
})

grid <- do.call(rbind, grid_list)

# merge observed data onto completed grid
out <- merge(
  grid,
  ema_selected,
  by = c("mlife_id", "date_rep", "ema_time"),
  all.x = TRUE
)

# keep compliance per participant
comp_map <- ema_selected[!is.na(ema_selected$compliance),
                         c("mlife_id", "compliance")]
comp_map <- comp_map[!duplicated(comp_map$mlife_id), ]

out$compliance <- ifelse(
  is.na(out$compliance),
  comp_map$compliance[match(out$mlife_id, comp_map$mlife_id)],
  out$compliance
)

# enforce correct ordering
out$ema_time <- factor(out$ema_time, levels = ema_levels, ordered = TRUE)
out <- out[order(out$mlife_id, out$date_rep, out$ema_time), ]

View(out)
ema_final <- out #change the data frame name

View(ema_final)

# combine date_rep + ema_timestamp_hhmm into a POSIXct datetime
ema_final2 <- ema_final %>%
  mutate(
    date_rep = as.Date(date_rep),
    dt = as_datetime(date_rep, tz = "UTC") + hm(as.character(ema_timestamp_hhmm))
  ) 
View(ema_final2)

# drop date_rep, ema_timestamp_hhmm
ema_final2 <- ema_final2 %>%
  select(-date_rep, -ema_timestamp_hhmm)
View(ema_final2)

# get dynamic indices

# helpers 
# Set entropy (Shannon entropy of the value distribution), in bits (log2)
set_entropy <- function(x, min_n = 2) {
  x <- x[!is.na(x)]
  if (length(x) < min_n) return(NA_real_)
  
  p <- table(x) / length(x)
  -sum(p * log2(p))
}

p90_p10 <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) < 2) return(NA_real_)
  as.numeric(
    stats::quantile(x, 0.9, names = FALSE, type = 7) -
      stats::quantile(x, 0.1, names = FALSE, type = 7)
  )
}

# quick check on your example vector:
set_entropy(c(1,2,NA,50,0,22,30,NA,0,24,2,NA,0,5,NA,0,2,18,1,0,NA))
# expected: ~2.852217


# choose the EMA variable you want features for (example: q1) 
# mad: (R) median absolute deviation (scaled by default constant)

vars <- c(paste0("q", 1:9), "g1", "g2")

ema_features <- ema_final2 %>%
  group_by(mlife_id) %>%
  mutate(
    across(
      all_of(vars),
      list(
        mean_d1     = ~slide_dbl(.x, ~mean(.x, na.rm = TRUE),
                                 .before = 2, .complete = TRUE),
        mad_d1      = ~slide_dbl(.x, ~mad(.x, na.rm = TRUE),
                                 .before = 2, .complete = TRUE),
        sentropy_d1 = ~slide_dbl(.x, ~set_entropy(.x, min_n = 2),
                                 .before = 2, .complete = TRUE),
        p90p10_d1   = ~slide_dbl(.x, p90_p10,
                                 .before = 2, .complete = TRUE),
        
        mean_d7      = ~slide_dbl(.x, ~mean(.x, na.rm = TRUE),
                                  .before = 20, .complete = TRUE),
        mad_d7       = ~slide_dbl(.x, ~mad(.x, na.rm = TRUE),
                                  .before = 20, .complete = TRUE),
        sentropy_d7  = ~slide_dbl(.x, ~set_entropy(.x, min_n = 2),
                                  .before = 20, .complete = TRUE),
        p90p10_d7    = ~slide_dbl(.x, p90_p10,
                                  .before = 20, .complete = TRUE)
      ),
      .names = "{.col}_{.fn}"
    )
  ) %>%
  ungroup()

View(ema_features)


# save the data frame
write.csv(ema_features, "ema_features.csv", row.names = FALSE)


###### 2. Digital Phenotype ####

# --- Merge data ---########
library(data.table)
library(parallel)

# 1) Define paths (HR)
main_folder <- "/Users/ahhyun/Desktop/SI_Predction"
folder_name <- "hr"
folder_path <- file.path(main_folder, folder_name)

# 2) Start timing
start_time <- Sys.time()

# 3) List CSV files
csv_files <- list.files(
  path = folder_path,
  pattern = "\\.csv$",
  full.names = TRUE,
  ignore.case = TRUE
)

if (length(csv_files) == 0) {
  stop("No CSV files found in: ", folder_path)
}

cat("\nFound", length(csv_files), "CSV files in", folder_name, "folder\n")

# 4) Filter out empty/unreadable files
cat("Checking for empty/unreadable files...\n")
valid_files <- character(0)

for (i in seq_along(csv_files)) {
  f <- csv_files[i]
  
  test_data <- tryCatch(
    fread(f, nrows = 1, showProgress = FALSE),
    error = function(e) NULL
  )
  
  if (is.null(test_data) || nrow(test_data) == 0) {
    cat("SKIPPED:", basename(f), "- empty or unreadable\n")
  } else {
    valid_files <- c(valid_files, f)
  }
  
  if (i %% 10 == 0 || i == length(csv_files)) {
    cat(sprintf("  Checked: %d/%d files (%.1f%%)\n",
                i, length(csv_files), (i / length(csv_files)) * 100))
  }
}

if (length(valid_files) == 0) {
  stop("All files are empty or unreadable!")
}

cat("\nProcessing", length(valid_files), "valid files...\n")

# 5) Read + combine
if (length(valid_files) == 1) {
  cat("Reading 1 file...\n")
  hr <- fread(valid_files[1], showProgress = FALSE)
  cat("Progress: 100%\n")
} else {
  num_cores <- min(14, length(valid_files))
  cat("Using", num_cores, "cores for parallel processing\n\n")
  
  cl <- makeCluster(num_cores)
  clusterEvalQ(cl, library(data.table))
  
  total_files <- length(valid_files)
  chunk_size <- ceiling(total_files / 10)
  data_list <- list()
  
  cat("Reading files in parallel...\n")
  
  for (chunk_start in seq(1, total_files, by = chunk_size)) {
    chunk_end <- min(chunk_start + chunk_size - 1, total_files)
    chunk_files <- valid_files[chunk_start:chunk_end]
    
    chunk_data <- parLapply(cl, chunk_files, function(file) {
      fread(file, showProgress = FALSE)
    })
    
    data_list <- c(data_list, chunk_data)
    
    progress_pct <- (chunk_end / total_files) * 100
    cat(sprintf("  Progress: %d/%d files (%.1f%%)\n",
                chunk_end, total_files, progress_pct))
  }
  
  stopCluster(cl)
  
  cat("\nCombining data...\n")
  hr <- rbindlist(data_list, fill = TRUE)
}

# 6) Remove V1 column if it exists
if ("V1" %in% names(hr)) {
  hr[, V1 := NULL]
}

# 7) Convert to data.frame if you really want (optional)
hr <- as.data.frame(hr)

# 8) End timing + summary
end_time <- Sys.time()
elapsed <- end_time - start_time

cat("\n========================================\n")
cat("SUCCESS!\n")
cat("Combined", nrow(hr), "rows,", ncol(hr), "columns\n")
cat("Time:", round(elapsed, 2), attr(elapsed, "units"), "\n")
cat("========================================\n")

# 9) Check unique UIDs (only if uid exists)
if ("uid" %in% names(hr)) {
  cat("\nUnique UIDs:", length(unique(hr$uid)), "\n")
} else {
  cat("\nNote: No 'uid' column found.\n")
}

# 10) Write output
write.csv(hr, file.path(folder_path, "hr.csv"), row.names = FALSE)

    
    
    
    
    
# --- load data  -- ####
    
setwd("/Users/ahhyun/Desktop/SI_Predction/Final")

accelerometer <- vroom("acceleration.csv") 
activity     <- fread("activity.csv")
battery      <- fread("battery.csv")
ema_features <- fread("ema_features.csv")
gps          <- fread("gps.csv")
hr           <- fread("hr.csv") 
rr           <- fread("rr.csv") #?
steps        <- fread("steps.csv")
stress       <- fread("stress.csv")
convo        <- fread("convo.csv")
sms          <- fread("smslog.csv")
call         <- fread("calllog.csv")

# --- data cleaning --- ####

colnames(ema_features) 
#  "mlife_id", "ema_time", "dt" (only indicated the one that I will use)
#  change mlife_id to uid - done
#  dt to ema_day - done

# data inspection
View(activity) #v1, data[#,#], day (event base), uid
# delete v1 - done
########## check out how to understand data
#{'still': 0, 'unknown': 1, 'tilting': 2, 'on foot': 3, 'on bike': 4, 'in vehicle': 5, 'walking': 6, 'running': 7}
View(convo) #v1, _id, event_data, event_data [start: #, end: #], event_time (event base), uid
# delete v1 - done
# delete _id - done
# change event_time to day - done
# extract start/end time of conversation - done
View(gps) #v1, data, day (every 15 mins), uid
# delete v1 - done
# data: extract altitude, longitude, accuracy, latitude, speed, bearing (what is it?) - done
View(hr) # data, day (every 10 seconds), uid
View(battery) #v1, data[#], day (hourly), uid
# delete v1 - done
View(rr) # v1, data, day (every minutes), uid
# delete v1 - done

View(steps) #v1, data [#,#], day (every 30 seconds), uid
# delete v1 - done
########### data [#,#]: make sure how to interpret
#[step count, duration, total steps, start time]

# day (why some of the data is in 1970?) - drop
View(stress) #v1, data, day (mot periodic), uid
# delete v1 - done
########### data: outlier for x - over 100?
# day (why some of the data is in 1970?) drop

# 1. delete V1 data
dfs <- c("activity", "battery", "convo",
         "gps", "hr", "rr", "steps", "stress")

for (d in dfs) {
  assign(
    d,
    get(d) %>% select(-any_of("V1"))
  )
}

# 2. extract number from data  (gps, conv, step counts)

# gps
# Fix doubled quotes
clean <- gsub('""', '"', gps$data)

# Parse each row
parsed <- lapply(clean, fromJSON)

# Helper to safely extract (returns NA if missing)
get_val <- function(x, name) if (!is.null(x[[name]])) x[[name]] else NA

# Create new columns
gps$alt  <- sapply(parsed, get_val, "ALTITUDE")
gps$lon <- sapply(parsed, get_val, "LONGITUDE")
gps$acc  <- sapply(parsed, get_val, "ACCURACY")
gps$lat  <- sapply(parsed, get_val, "LATITUDE")
gps$speed     <- sapply(parsed, get_val, "SPEED")
gps$bearing   <- sapply(parsed, get_val, "BEARING")

View(gps)

# conv
clean  <- gsub('""', '"', convo$event_data, fixed = TRUE)
parsed <- lapply(clean, fromJSON, simplifyVector = FALSE)

get_val <- function(x, name) {
  v <- x[[name]]
  if (is.null(v)) return(NA_character_)
  as.character(v)   # <- key: keep as string
}

convo$start <- as.integer64(sapply(parsed, get_val, "CONVERSATION_START"))
convo$end   <- as.integer64(sapply(parsed, get_val, "CONVERSATION_END"))
View(convo)

# change the column name
# drop rows with odd day values (start with 1970)
ema_features <- ema_features %>%
  rename(uid = mlife_id) %>%
  rename(ema_day = dt)
colnames(ema_features) 

convo <- convo %>%
  rename(day = event_time) %>%
  select(-"_id")

# step counts 
View(steps) #[step count, duration, total steps, start time] →only use the step count 
steps_clean <- steps %>%
  mutate(
    step_count = as.numeric(str_extract(data, "(?<=\\[)[^,]+"))
  ) %>%
  select(uid, day, step_count)
View(steps_clean)

# drop rows with outlier (steps, stress - drop rows with 1970)
steps <- steps %>%
  dplyr::filter(lubridate::year(day) >= 2000)

stress <- stress %>%
  dplyr::filter(lubridate::year(day) >= 2000)
fwrite(stress, "stress_clean.csv") # save data


# 3. Data calculation (1): event based timing -> every 10 mins (Convo, SMS, Calls, steps) #######
# change the event based digital phenotype data to regular distribution

# conversation
View(convo) #start(UTC) #end(UTC) #day
# change to UTC
convo$end_time <- as.POSIXct(as.numeric(convo$end) / 1000, origin = "1970-01-01", tz = "UTC")
convo$start_time <- as.POSIXct(as.numeric(convo$start) / 1000, origin = "1970-01-01", tz = "UTC")
summary(convo)

library(data.table)
library(parallel)

# The following code is adapted from:
# https://stackoverflow.com/questions/57245771/how-do-i-fill-time-sequence-based-on-start-and-end-time-in-r

# seq.POSIXt() creates a sequence of times (e.g., 2021-10-25 17:17:01 CEST, 2021-10-25 17:17:02 CEST, 2021-10-25 17:17:03 CEST). This function is vectorized, which means that it is directly performed on the whole vector (= column).
vec_seq <- Vectorize(seq.POSIXt, vectorize.args = c("from", "to"))

library(dplyr)
library(lubridate)

convo_per10min <- convo %>%
  # First create per-minute data
  group_by(uid) %>%
  reframe(
    Date = seq.POSIXt(min(start_time), max(end_time), by = "1 min")
  ) %>%
  left_join(
    convo %>% 
      mutate(Date = map2(start_time, end_time, ~seq.POSIXt(.x, .y, by = "1 min"))) %>%
      unnest(Date) %>%
      mutate(on_phone = 1) %>%  # Mark as on phone
      select(uid, Date, on_phone),
    by = c("uid", "Date")
  ) %>%
  # Replace NA with 0
  mutate(on_phone = ifelse(is.na(on_phone), 0, 1)) %>%
  # Create 10-minute bins and count
  mutate(time_bin = floor_date(Date, "10 minutes")) %>%
  group_by(uid, time_bin) %>%
  summarise(data = sum(on_phone), .groups = "drop")

View(convo_per10min)
summary(convo_per10min)
hist(convo_per10min$data)
fwrite(convo_per10min, "conv_clean.csv") # save data

# SMS
View(sms)

sms_per10min <- sms %>%
  # First create complete 10-minute sequence for each uid
  group_by(uid) %>%
  reframe(
    time_bin = seq.POSIXt(
      floor_date(min(time), "10 minutes"),  # Start at first 10-min mark
      floor_date(max(time), "10 minutes"),  # End at last 10-min mark
      by = "10 min"
    )
  ) %>%
  # Join with actual SMS data aggregated by 10-minute bins
  left_join(
    sms %>%
      mutate(time_bin = floor_date(time, "10 minutes")) %>%
      group_by(uid, time_bin) %>%
      summarise(total_body_length = sum(body_len, na.rm = TRUE), .groups = "drop"),
    by = c("uid", "time_bin")
  ) %>%
  # Replace NA with 0 (no messages in that window)
  mutate(total_body_length = ifelse(is.na(total_body_length), 0, total_body_length))

View(sms_per10min)
fwrite(sms_per10min, "sms_clean.csv") # save data

#call
View(call)

call_per10min <- call %>%
  # First create complete 10-minute sequence for each uid
  group_by(uid) %>%
  reframe(
    time_bin = seq.POSIXt(
      floor_date(min(time), "10 minutes"),  # Start at first 10-min mark
      floor_date(max(time), "10 minutes"),  # End at last 10-min mark
      by = "10 min"
    )
  ) %>%
  # Join with actual call data aggregated by 10-minute bins
  left_join(
    call %>%
      mutate(time_bin = floor_date(time, "10 minutes")) %>%
      group_by(uid, time_bin) %>%
      summarise(total_duration = sum(duration, na.rm = TRUE), .groups = "drop"),
    by = c("uid", "time_bin")
  ) %>%
  # Replace NA with 0 (no calls in that window)
  mutate(total_duration = ifelse(is.na(total_duration), 0, total_duration))

View(call_per10min)
fwrite(call_per10min, "call_clean.csv") # save data

# steps
steps_per10min <- steps_clean %>%
  # Create complete 10-minute sequence for each uid
  group_by(uid) %>%
  reframe(
    time_bin = seq.POSIXt(
      floor_date(min(day), "10 minutes"),
      floor_date(max(day), "10 minutes"),
      by = "10 min"
    )
  ) %>%
  # Join with actual step data aggregated by 10-minute bins
  left_join(
    steps_clean %>%
      mutate(time_bin = floor_date(day, "10 minutes")) %>%
      group_by(uid, time_bin) %>%
      summarise(total_step_count = sum(step_count, na.rm = TRUE), .groups = "drop"),
    by = c("uid", "time_bin")
  ) %>%
  # Replace NA with 0 (no steps in that window)
  mutate(total_step_count = ifelse(is.na(total_step_count), 0, total_step_count))

View(steps_per10min)
fwrite(steps_per10min, "steps_clean.csv") # save data

# # 3. Data calculation: time-based data (accelerometer, GPS)
#
setwd("/Users/ahhyun/Desktop/SI_Predction/Final")
library(vroom)
accelerometer <- vroom("acceleration.csv") 
summa(accelerometer)
library(data.table)

# Ensure it's a data.table
setDT(accelerometer)

# Check current size
cat("Dataset size:", format(object.size(accelerometer), units = "GB"), "\n")
cat("Dimensions:", nrow(accelerometer), "rows x", ncol(accelerometer), "cols\n")

# OPTION 1: Direct aggregation (try this first - often fastest)
# This does everything in one pass without copies
cat("Starting aggregation...\n")
system.time({
  results <- accelerometer[, .(
    x = mean(x),
    y = mean(y),
    z = mean(z),
    day = min(day)
  ), by = .(uid, time_5sec = floor(as.numeric(day) / 5) * 5)]
})

cat("Results size:", format(object.size(results), units = "MB"), "\n")
cat("Results rows:", nrow(results), "\n")

#gps 
# distance
# interval: between 10 mins to 15 mins, distance as meter
# Calculate distance between consecutive GPS points
gps_movement <- gps %>%
  arrange(uid, day) %>%
  group_by(uid) %>%
  mutate(
    # Calculate distance from previous point (in meters)
    distance_m = distHaversine(
      cbind(lag(longitude), lag(latitude)),
      cbind(longitude, latitude)
    ),
    # Time difference from previous point (in minutes)
    time_diff_min = as.numeric(difftime(day, lag(day), units = "mins")),
    # Set first point distance to 0
    distance_m = ifelse(is.na(distance_m), 0, distance_m)
  ) %>%
  ungroup()
View(gps_movement)
fwrite(gps_movement, "gps_clean_1.csv") #save_data

# time at home ######
# get a home location 
# calculate time to spend time at home
# calculate by every 10 mins

View(gps)

# Dawson's code 

# ===============================================================================
# Import Libraries & Read Data
# ===============================================================================

# Load required libraries
library(tidyverse)    # For data manipulation
library(lubridate)    # For date/time handling
library(dbscan)       # For DBSCAN clustering
library(geosphere)    # For geographic distance calculations

# Read data #gps

# Display data overview
loc_df <- gps
loc_df$time <- gps$day
print(head(loc_df))
cat("Data shape:", nrow(loc_df), "rows x", ncol(loc_df), "columns\n")

# REMOVED THE SUBSET - Process ALL participants now
cat("Number of unique UIDs:", length(unique(loc_df$uid)), "\n")

# ===============================================================================
# Subset to Nighttime for Home Detection
# ===============================================================================

cat("\nFiltering to nighttime hours (10pm-7am)...\n")

# Convert 'time' to datetime format
loc_df <- loc_df %>%
  mutate(time = as.POSIXct(time))

# Define nighttime range (10pm–7am)
night_start <- hms::as_hms("22:00:00")
night_end <- hms::as_hms("07:00:00")

# Extract time component
loc_df <- loc_df %>%
  mutate(time_only = hms::as_hms(time))

# Filter data to nighttime hours (10pm–7am) for home detection
loc_night <- loc_df %>%
  filter(time_only >= night_start | time_only <= night_end)

cat("Nighttime data shape:", nrow(loc_night), "rows\n")
cat("Unique UIDs in nighttime data:", length(unique(loc_night$uid)), "\n")

# ===============================================================================
# DBSCAN Clustering for Home Location Identification
# ===============================================================================

# Function to apply DBSCAN to find home locations
find_home_locations <- function(df, uid_val, eps_meters = 50, min_samples = 10) {
  # Filter data for specific participant
  participant_data <- df %>%
    filter(uid == uid_val)
  
  # Return NULL if no data
  if (nrow(participant_data) == 0) {
    return(NULL)
  }
  
  # Extract coordinates
  coords <- participant_data %>%
    select(lat, lon) %>%
    as.matrix()
  
  # Convert eps from meters to degrees (approximate)
  eps_degrees <- eps_meters / 111000
  
  # Apply DBSCAN
  db_result <- dbscan(coords, eps = eps_degrees, minPts = min_samples)
  
  # Add cluster labels to participant data
  participant_data <- participant_data %>%
    mutate(cluster = db_result$cluster)
  
  # Find centroids of each cluster (excluding noise points labeled as 0)
  home_coords <- participant_data %>%
    filter(cluster > 0) %>%
    group_by(cluster) %>%
    summarise(
      lat = mean(lat),
      lon = mean(lon),
      n_points = n(),
      .groups = 'drop'
    ) %>%
    arrange(desc(n_points))
  
  return(home_coords)
}

# Function to calculate distance between two points
calc_distance <- function(lat1, lon1, lat2, lon2) {
  distHaversine(c(lon1, lat1), c(lon2, lat2))
}

# Function to check if a point is near any home location
is_near_home <- function(lat, lon, home_coords, threshold_meters = 100) {
  if (is.null(home_coords) || nrow(home_coords) == 0) {
    return(FALSE)
  }
  
  distances <- sapply(1:nrow(home_coords), function(i) {
    calc_distance(lat, lon, home_coords$lat[i], home_coords$lon[i])
  })
  
  return(any(distances <= threshold_meters))
}

# ===============================================================================
# Process All Participants - Find Home Locations
# ===============================================================================

cat("\n=== Finding home locations for ALL participants ===\n")

# Get unique participant IDs
unique_uids <- unique(loc_night$uid)
cat("Processing", length(unique_uids), "participants...\n")

# Initialize home coordinates list
home_coords_list <- list()

# Process each participant with progress tracking
for (i in seq_along(unique_uids)) {
  uid_val <- unique_uids[i]
  
  # Print progress every 10 participants
  if (i %% 10 == 0) {
    cat("Progress:", i, "/", length(unique_uids), "participants processed\n")
  }
  
  home_coords <- find_home_locations(loc_night, uid_val, eps_meters = 50, min_samples = 10)
  
  if (!is.null(home_coords) && nrow(home_coords) > 0) {
    home_coords_list[[uid_val]] <- home_coords
  }
}

cat("\nCompleted home location detection!\n")
cat("Participants with detected homes:", length(home_coords_list), "/", length(unique_uids), "\n")
cat("Participants without detected homes:", length(unique_uids) - length(home_coords_list), "\n")

# ===============================================================================
# Calculate Time at Home for Each 10-Minute Interval
# ===============================================================================

cat("\n=== Calculating time at home in 10-minute intervals ===\n")

# Function to process a single participant
calculate_time_at_home <- function(df, uid_val, home_coords, threshold_meters = 100) {
  
  # Filter data for this participant
  participant_data <- df %>%
    filter(uid == uid_val)
  
  if (nrow(participant_data) == 0) {
    return(NULL)
  }
  
  # Check if each point is at home
  participant_data <- participant_data %>%
    rowwise() %>%
    mutate(is_home = is_near_home(lat, lon, home_coords, threshold_meters)) %>%
    ungroup()
  
  # Create 10-minute intervals
  # Floor time to nearest 10-minute interval
  participant_data <- participant_data %>%
    mutate(
      time_interval = floor_date(time, "10 minutes")
    )
  
  # Calculate time spent at home per 10-minute interval
  time_summary <- participant_data %>%
    group_by(uid, time_interval) %>%
    summarise(
      total_points = n(),
      points_at_home = sum(is_home),
      points_away = sum(!is_home),
      prop_at_home = mean(is_home),
      .groups = 'drop'
    ) %>%
    mutate(
      # Assuming GPS samples are relatively evenly distributed within the interval
      # Calculate approximate minutes at home (this is a rough estimate)
      # For more accurate results, you'd need to know the exact sampling frequency
      minutes_at_home = prop_at_home * 10,
      minutes_away = (1 - prop_at_home) * 10
    )
  
  return(time_summary)
}

# Process all participants with home locations
all_time_summaries <- list()

uids_to_process <- names(home_coords_list)
cat("Processing time at home for", length(uids_to_process), "participants with detected homes...\n")

for (i in seq_along(uids_to_process)) {
  uid_val <- uids_to_process[i]
  
  # Print progress every 10 participants
  if (i %% 10 == 0) {
    cat("Progress:", i, "/", length(uids_to_process), "participants processed\n")
  }
  
  home_coords <- home_coords_list[[uid_val]]
  
  time_summary <- calculate_time_at_home(
    loc_df, 
    uid_val, 
    home_coords, 
    threshold_meters = 100
  )
  
  if (!is.null(time_summary)) {
    all_time_summaries[[uid_val]] <- time_summary
  }
}

cat("\nCompleted time at home calculations!\n")

# Combine all summaries into one dataframe
time_at_home_df <- bind_rows(all_time_summaries)

print(head(time_at_home_df))

write.csv(time_at_home_df, "gps_clean_2.csv")
# ===============================================================================
# Display and Save Results
# ===============================================================================

# Display results
cat("\n=== Time at Home Summary (10-minute intervals) ===\n")
print(head(time_at_home_df, 20))
cat("\nTotal intervals:", nrow(time_at_home_df), "\n")
cat("Unique UIDs in results:", length(unique(time_at_home_df$uid)), "\n")

# Summary statistics by participant
participant_summary <- time_at_home_df %>%
  group_by(uid) %>%
  summarise(
    total_intervals = n(),
    total_minutes_at_home = sum(minutes_at_home),
    total_minutes_away = sum(minutes_away),
    avg_minutes_at_home_per_interval = mean(minutes_at_home),
    prop_time_at_home = sum(minutes_at_home) / (sum(minutes_at_home) + sum(minutes_away)),
    .groups = 'drop'
  ) %>%
  mutate(
    total_hours_at_home = total_minutes_at_home / 60,
    total_hours_away = total_minutes_away / 60
  )

cat("\n=== Summary by Participant ===\n")
print(head(participant_summary, 20))
cat("Total participants in summary:", nrow(participant_summary), "\n")

# Daily summary
daily_summary <- time_at_home_df %>%
  mutate(date = as.Date(time_interval)) %>%
  group_by(uid, date) %>%
  summarise(
    total_minutes_at_home = sum(minutes_at_home),
    total_minutes_away = sum(minutes_away),
    num_intervals = n(),
    .groups = 'drop'
  ) %>%
  mutate(
    hours_at_home = total_minutes_at_home / 60,
    hours_away = total_minutes_away / 60,
    prop_time_at_home = total_minutes_at_home / (total_minutes_at_home + total_minutes_away)
  )

cat("\n=== Daily Summary (first 20 rows) ===\n")
print(head(daily_summary, 20))

# ===============================================================================
# Create summary of participants WITHOUT detected homes
# ===============================================================================

all_uids <- unique(loc_df$uid)
uids_with_homes <- names(home_coords_list)
uids_without_homes <- setdiff(all_uids, uids_with_homes)

if (length(uids_without_homes) > 0) {
  cat("\n=== Participants WITHOUT detected home locations ===\n")
  cat("Count:", length(uids_without_homes), "\n")
  
  # Create dataframe of participants without homes
  participants_no_home <- data.frame(
    uid = uids_without_homes,
    has_home_location = FALSE,
    reason = "Insufficient nighttime GPS clusters"
  )
  
  # Get nighttime data counts for these participants
  nighttime_counts <- loc_night %>%
    filter(uid %in% uids_without_homes) %>%
    group_by(uid) %>%
    summarise(nighttime_points = n(), .groups = 'drop')
  
  participants_no_home <- participants_no_home %>%
    left_join(nighttime_counts, by = "uid")
  
  print(head(participants_no_home, 20))
  
  # Save this information
  write_csv(participants_no_home, "participants_without_home_locations.csv")
  cat("\nSaved to: participants_without_home_locations.csv\n")
}

# ===============================================================================
# Save Results to CSV
# ===============================================================================

cat("\n=== Saving results to CSV files ===\n")

# Save 10-minute interval data
write_csv(time_at_home_df, "time_at_home_10min_intervals_ALL.csv")
cat("✓ 10-minute interval data saved to: time_at_home_10min_intervals_ALL.csv\n")
cat("  Rows:", nrow(time_at_home_df), "\n")
cat("  Unique UIDs:", length(unique(time_at_home_df$uid)), "\n")

# Save participant summary
write_csv(participant_summary, "time_at_home_participant_summary_ALL.csv")
cat("✓ Participant summary saved to: time_at_home_participant_summary_ALL.csv\n")
cat("  Participants:", nrow(participant_summary), "\n")

# Save daily summary
write_csv(daily_summary, "time_at_home_daily_summary_ALL.csv")
cat("✓ Daily summary saved to: time_at_home_daily_summary_ALL.csv\n")
cat("  Rows:", nrow(daily_summary), "\n")

cat("\n=== Processing Complete! ===\n")
cat("Original data: ", length(unique(gps$uid)), "unique UIDs\n")
cat("With detected homes:", length(unique(time_at_home_df$uid)), "UIDs\n")
cat("Without detected homes:", length(uids_without_homes), "UIDs\n")
cat("\nFiles created:\n")
cat("1. time_at_home_10min_intervals_ALL.csv - Detailed 10-minute interval data\n")
cat("2. time_at_home_participant_summary_ALL.csv - Overall summary per participant\n")
cat("3. time_at_home_daily_summary_ALL.csv - Daily summary per participant\n")
cat("4. participants_without_home_locations.csv - List of participants without detected homes\n")

# View the results
View(time_at_home_df)
cat("\nFinal check - Unique UIDs in results:", length(unique(time_at_home_df$uid)), "\n")




####################################
# Data inspection #

# 1. Count UIDs with 1970 year in date data
# Assuming date column is named 'date_column' and uid column is 'uid'

uids_with_1970 <- stress %>%
  filter(format(as.Date(day), "%Y") == "1970") %>%
  summarise(unique_uids = n_distinct(uid))

print(paste("Number of unique UIDs with 1970 dates:", uids_with_1970$unique_uids))
#302
# For timestamp column - checking 1970 dates
uids_with_1970 <- stress %>%
  mutate(year = as.integer(format(as.POSIXct(day, origin="1970-01-01"), "%Y"))) %>%
  filter(year == 1970) %>%
  summarise(unique_uids = n_distinct(uid))


# 2. Count and percentage of rows with value 4294967294
# Assuming the column is named 'your_column'

rows_with_value <- stress %>%
  summarise(
    total_rows = n(),
    rows_with_4294967294 = sum(data == 4294967294, na.rm = TRUE),
    percentage = (sum(data == 4294967294, na.rm = TRUE) / n()) * 100
  )

print(paste("Total rows:", rows_with_value$total_rows))
print(paste("Rows with 4294967294:", rows_with_value$rows_with_4294967294))
print(paste("Percentage:", round(rows_with_value$percentage, 2), "%")) #7.15 %

# change this to NA?

# 3. Unlock frequncy?

unlock <- fread("unlock.csv")
length(unique(unlock$uid))

library(dplyr)
library(stringr)

# Faster method: Count commas + 1 (assuming non-empty arrays)
unlock_with_counts <- unlock %>%
  mutate(
    num_days_with_data = str_count(data, ",") + 1
  )

# Distribution
distribution <- unlock_with_counts %>%
  group_by(num_days_with_data) %>%
  summarise(
    count_uids = n(),
    percentage = (n() / nrow(unlock_with_counts)) * 100
  ) %>%
  arrange(num_days_with_data)

print(distribution)

# Summary stats
summary_stats <- unlock_with_counts %>%
  summarise(
    total_uids = n(),
    min_count = min(num_days_with_data),
    max_count = max(num_days_with_data),
    mean_count = mean(num_days_with_data),
    median_count = median(num_days_with_data)
  )

print(summary_stats)

# Visualization
library(ggplot2)
ggplot(distribution, aes(x = num_days_with_data, y = count_uids)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  labs(
    title = "Distribution of Array Length per UID",
    x = "Number of Values",
    y = "Number of UIDs"
  ) +
  theme_minimal()

print(head(unlock))

####### clean accelerate one by one #######
install.packages(c("data.table", "future.apply", "progressr"))
library(data.table)
library(future)
library(future.apply)
library(progressr)

setwd("/Users/ahhyun/Desktop/SI_Predction/Before_Merged/acceleration")
file_list <- list.files(pattern="\\.csv$", full.names = TRUE)

out_dir <- file.path(getwd(), "processed")
dir.create(out_dir, showWarnings = FALSE)

# macOS-friendly parallel
plan(multisession, workers = max(1, availableCores() - 1))

# atomic write (prevents half-written files)
atomic_write_csv <- function(dt, out_file) {
  tmp <- paste0(out_file, ".tmp_", Sys.getpid(), "_", as.integer(runif(1, 1, 1e9)))
  fwrite(dt, tmp)
  file.rename(tmp, out_file)
}

process_accelerometer_file_fast <- function(file_path) {
  if (file.size(file_path) == 0) return(NULL)
  
  dt <- tryCatch(
    fread(file_path, showProgress = FALSE, fill = TRUE, warn = FALSE),
    error = function(e) NULL
  )
  if (is.null(dt) || nrow(dt) == 0) return(NULL)
  if (!("data" %in% names(dt))) return(NULL)
  
  # data like "[#,#,#,#]" -> remove brackets
  s <- gsub("\\[|\\]", "", as.character(dt[["data"]]))
  
  # split once, take first 3 values only (ignore 4th+)
  xyz <- tstrsplit(s, ",", fixed = TRUE, fill = NA_character_)
  if (length(xyz) < 3) return(NULL)
  
  x <- as.numeric(xyz[[1]])
  y <- as.numeric(xyz[[2]])
  z <- as.numeric(xyz[[3]])
  
  dt[, `:=`(x = x, y = y, z = z)]
  dt[, magnitude := sqrt(x*x + y*y + z*z)]
  
  dt
}

handlers(global = TRUE)
out_files <- with_progress({
  p <- progressor(along = file_list)
  
  future_lapply(file_list, function(fp) {
    p(basename(fp))
    
    out_file <- file.path(
      out_dir,
      paste0(tools::file_path_sans_ext(basename(fp)), "_processed.csv")
    )
    
    # resume: skip already processed outputs
    if (file.exists(out_file) && file.info(out_file)$size > 0) return(out_file)
    
    dt <- process_accelerometer_file_fast(fp)
    if (is.null(dt)) return(NA_character_)
    
    atomic_write_csv(dt, out_file)
    out_file
  }, future.seed = TRUE)
})

# optional: quick summary
cat("Done.\nProcessed:", sum(!is.na(out_files)), "of", length(out_files), "\n")
cat("Outputs in:", out_dir, "\n")

failed_files <- file_list[is.na(out_files)]


length(failed_files)
failed_files # empty files

###### Merge the files ###########



# Efficient Data Merging Script for Large Files
# Extracts only day, uid, and magnitude columns and combines them

# Install required package if not already installed
if (!require("data.table")) {
  install.packages("data.table")
}

library(data.table)

# Set your folder path here
folder_path <- "/Users/ahhyun/Desktop/SI_Predction/Before_Merged/acceleration/processed"

# Get list of all data files in the folder
# Adjust the pattern based on your file type (e.g., "*.csv", "*.txt", "*.tsv")
file_list <- list.files(path = folder_path, 
                        pattern = "\\.(csv|txt|tsv)$", 
                        full.names = TRUE)

print(paste("Found", length(file_list), "files to process"))

# Function to read only specific columns efficiently
read_selected_columns <- function(file_path) {
  # Using fread for fast reading of large files
  # select parameter reads only specified columns
  tryCatch({
    dt <- fread(file_path, 
                select = c("day", "uid", "magnitude"),
                showProgress = FALSE)
    return(dt)
  }, error = function(e) {
    warning(paste("Error reading file:", file_path, "-", e$message))
    return(NULL)
  })
}

# Read and combine all files efficiently
print("Reading and merging files...")

# Use rbindlist which is very efficient for combining data.tables
combined_data <- rbindlist(lapply(file_list, read_selected_columns), 
                           fill = TRUE)

print(paste("Total rows in combined data:", nrow(combined_data)))
print(paste("Columns:", paste(names(combined_data), collapse = ", ")))

# Display summary
print("Data summary:")
print(summary(combined_data))

# Save the merged data
output_file <- "merged_data.csv"
fwrite(combined_data, output_file, showProgress = TRUE)

print(paste("Merged data saved to:", output_file))

# Optional: Clear memory
# rm(combined_data)
# gc()
      
length(unique(combined_data$uid))
head(combined_data)
