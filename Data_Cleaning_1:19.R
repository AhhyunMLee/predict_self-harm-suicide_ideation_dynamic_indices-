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
cat("\nWrote:", f
    
    
    
    
    
# --- load data  -- ####
    
setwd("/Users/ahhyun/Desktop/SI_Predction/Final")
    
# acceleration <- fread("acceleration.csv") #taking too long to download # will work on it later
activity     <- fread("activity.csv")
battery      <- fread("battery.csv")
ema_features <- fread("ema_features.csv")
gps          <- fread("gps.csv")
hr           <- fread("hr.csv") 
rr           <- fread("rr.csv") 
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

# 2. extract number from data 

# Fix doubled quotes
clean <- gsub('""', '"', gps$data)

# Parse each row
parsed <- lapply(clean, fromJSON)

# Helper to safely extract (returns NA if missing)
get_val <- function(x, name) if (!is.null(x[[name]])) x[[name]] else NA

# Create new columns
gps$altitude  <- sapply(parsed, get_val, "ALTITUDE")
gps$longitude <- sapply(parsed, get_val, "LONGITUDE")
gps$accuracy  <- sapply(parsed, get_val, "ACCURACY")
gps$latitude  <- sapply(parsed, get_val, "LATITUDE")
gps$speed     <- sapply(parsed, get_val, "SPEED")
gps$bearing   <- sapply(parsed, get_val, "BEARING")

View(gps)

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

steps <- steps %>%
  dplyr::filter(lubridate::year(day) >= 2000)

stress <- stress %>%
  dplyr::filter(lubridate::year(day) >= 2000)

# 3. Data calculation (Convo, SMS, Calls)
# change the event based digital phenotype data to regular distribution

# conversation
View(Convo) #start(UTC) #end(UTC) #day
convo$end_time <- as.POSIXct(as.numeric(convo$end) / 1000, origin = "1970-01-01", tz = "UTC")
convo$start_time <- as.POSIXct(as.numeric(convo$start) / 1000, origin = "1970-01-01", tz = "UTC")
View(convo)
summary(convo) #found 1 NA values


library(data.table)
library(pbapply)
library(parallel)

library(data.table)
library(pbapply)
library(parallel)

dt <- as.data.table(convo)
dt[, `:=`(start_time = as.POSIXct(start_time),
          end_time = as.POSIXct(end_time))]

uid_list <- split(dt, by = "uid")

# Set up cluster
cl <- makeCluster(12)
clusterEvalQ(cl, library(data.table))

# Expand events with progress bar
results <- pblapply(uid_list, function(x) {
  tryCatch({
    # Expand events
    events <- rbindlist(lapply(1:nrow(x), function(i) {
      if (x$start_time[i] >= x$end_time[i]) {
        return(data.table(
          uid = x$uid[i],
          time = x$start_time[i],
          event_source = x$event_source[i]
        ))
      }
      data.table(
        uid = x$uid[i],
        time = seq(x$start_time[i], x$end_time[i], by = "1 sec"),
        event_source = x$event_source[i]
      )
    }))
    
    # Create complete timeline with NAs for gaps
    complete_timeline <- data.table(
      uid = x$uid[1],
      time = seq(min(events$time), max(events$time), by = "1 sec")
    )
    
    # Merge to fill in event_source (NAs where no event)
    result <- merge(complete_timeline, events, by = c("uid", "time"), all.x = TRUE)
    
    return(result)
  }, error = function(e) {
    data.table(uid = character(), time = as.POSIXct(character()), event_source = character())
  })
}, cl = cl)

stopCluster(cl)

# Combine results
expanded_data <- rbindlist(results)
expanded_data <- rbindlist(results)