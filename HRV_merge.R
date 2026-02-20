#!/usr/bin/env Rscript
# ============================================================================
# HRV File Merger for R (FIXED VERSION)
# ============================================================================
# Purpose: Merges Daily HRV Summary and HRV Details files from multiple participants
# Input:   Folder with uid#### subfolders, each containing CSV files
# Output:  2 master CSV files combining all participants
#
# FIXES IN THIS VERSION:
#   1) Classifies files by FILENAME (not by first line inside the file)
#   2) Reads CSVs WITHOUT skipping the first row (since your files start with headers)
#   3) Optionally searches recursively for CSVs (in case there are nested folders)
#
# USAGE:
#   1. Update base_folder and output_folder paths below
#   2. Run in R: source("merge_hrv_fixed.R")
#   3. Or run from terminal: Rscript merge_hrv_fixed.R
#
# REQUIRED PACKAGES:
#   install.packages(c("dplyr", "readr", "purrr"))
# ============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(purrr)
})

# ============================================================================
# CONFIGURATION - UPDATE THESE PATHS FOR YOUR DATA
# ============================================================================
base_folder   <- "/Users/ahhyun/Desktop/Elizabeth_HRV_all"  # Folder containing uid#### subfolders
output_folder <- "/Users/ahhyun/Desktop/Elizabeth_HRV_all"                    # Where to save merged files

# ============================================================================
# FUNCTIONS
# ============================================================================

#' Extract UID from folder name
#' @param folder_name Character string like "uid1022"
#' @return Character string like "1022"
extract_uid <- function(folder_name) {
  gsub("^uid", "", folder_name, ignore.case = TRUE)
}

#' Read HRV CSV file and add uid column
#' NOTE: Your files already start with column headers, so we DO NOT skip rows.
#' @param file_path Path to CSV file
#' @param uid Participant UID
#' @return Data frame with uid column added
read_hrv_file <- function(file_path, uid) {
  df <- read_csv(file_path, show_col_types = FALSE)
  df$uid <- uid
  df
}

#' Merge all HRV files from all participants
#' @param base_folder Path to folder containing uid#### subfolders
#' @param output_folder Path to save output files
merge_all_hrv_files <- function(base_folder, output_folder) {
  
  # Create output folder if it doesn't exist
  if (!dir.exists(output_folder)) {
    dir.create(output_folder, recursive = TRUE)
  }
  
  # Find all uid folders
  all_folders <- list.dirs(base_folder, full.names = TRUE, recursive = FALSE)
  uid_folders <- all_folders[grepl("^uid[0-9]{4}$", basename(all_folders), ignore.case = TRUE)]
  uid_folders <- sort(uid_folders)
  
  cat("============================================================\n")
  cat("HRV File Merger - Master Files (R) [FIXED]\n")
  cat("============================================================\n\n")
  cat(sprintf("Base folder:   %s\n", base_folder))
  cat(sprintf("Output folder: %s\n\n", output_folder))
  cat(sprintf("Found %d participant folder(s)\n\n", length(uid_folders)))
  
  # Initialize lists to store all data
  all_daily_summary <- list()
  all_hrv_details   <- list()
  
  # Process each participant folder
  for (uid_folder in uid_folders) {
    folder_name <- basename(uid_folder)
    uid <- extract_uid(folder_name)
    
    cat(sprintf("Processing %s (UID: %s)\n", folder_name, uid))
    
    # Get all CSV files in this folder (recursive in case files are nested)
    all_files <- list.files(
      uid_folder,
      pattern = "\\.csv$",
      full.names = TRUE,
      recursive = TRUE,
      ignore.case = TRUE
    )
    
    if (length(all_files) == 0) {
      cat("  No CSV files found\n\n")
      next
    }
    
    # Categorize files by FILENAME (more reliable than reading file first line)
    daily_summary_files <- all_files[
      grepl("Daily Heart Rate Variability Summary", basename(all_files), ignore.case = TRUE)
    ]
    
    details_files <- all_files[
      grepl("Heart Rate Variability Details", basename(all_files), ignore.case = TRUE)
    ]
    
    # Read Daily Summary files
    if (length(daily_summary_files) > 0) {
      cat(sprintf("  Found %d Daily Summary file(s)\n", length(daily_summary_files)))
      for (file in daily_summary_files) {
        df <- read_hrv_file(file, uid)
        all_daily_summary[[length(all_daily_summary) + 1]] <- df
      }
    } else {
      cat("  Found 0 Daily Summary file(s)\n")
    }
    
    # Read Details files
    if (length(details_files) > 0) {
      cat(sprintf("  Found %d Details file(s)\n", length(details_files)))
      for (file in details_files) {
        df <- read_hrv_file(file, uid)
        all_hrv_details[[length(all_hrv_details) + 1]] <- df
      }
    } else {
      cat("  Found 0 Details file(s)\n")
    }
    
    cat("\n")
  }
  
  wrote_anything <- FALSE
  
  # Merge all Daily Summary data
  if (length(all_daily_summary) > 0) {
    cat("============================================================\n")
    cat("Merging ALL Daily HRV Summary data...\n")
    
    merged_daily <- bind_rows(all_daily_summary)
    
    # If timestamp exists, reorder and sort
    if ("timestamp" %in% names(merged_daily)) {
      other_cols <- setdiff(names(merged_daily), c("uid", "timestamp"))
      merged_daily <- merged_daily %>%
        select(uid, timestamp, all_of(other_cols)) %>%
        arrange(uid, timestamp)
    } else {
      merged_daily <- merged_daily %>% select(uid, everything())
    }
    
    output_file <- file.path(output_folder, "all_daily_hrv_summary.csv")
    write_csv(merged_daily, output_file)
    
    cat(sprintf("✓ Saved: %s\n", output_file))
    cat(sprintf("  Total rows: %d\n", nrow(merged_daily)))
    cat(sprintf("  Unique participants: %d\n", n_distinct(merged_daily$uid)))
    wrote_anything <- TRUE
  } else {
    cat("No Daily Summary files were merged.\n")
  }
  
  # Merge all HRV Details data
  if (length(all_hrv_details) > 0) {
    cat("\nMerging ALL HRV Details data...\n")
    
    merged_details <- bind_rows(all_hrv_details)
    
    if ("timestamp" %in% names(merged_details)) {
      other_cols <- setdiff(names(merged_details), c("uid", "timestamp"))
      merged_details <- merged_details %>%
        select(uid, timestamp, all_of(other_cols)) %>%
        arrange(uid, timestamp)
    } else {
      merged_details <- merged_details %>% select(uid, everything())
    }
    
    output_file <- file.path(output_folder, "all_hrv_details.csv")
    write_csv(merged_details, output_file)
    
    cat(sprintf("✓ Saved: %s\n", output_file))
    cat(sprintf("  Total rows: %d\n", nrow(merged_details)))
    cat(sprintf("  Unique participants: %d\n", n_distinct(merged_details$uid)))
    wrote_anything <- TRUE
  } else {
    cat("No Details files were merged.\n")
  }
  
  cat("\n============================================================\n")
  if (wrote_anything) {
    cat("Merge completed! Files saved in your output folder.\n")
  } else {
    cat("Merge completed, but no matching files were found to merge.\n")
    cat("Check that filenames include:\n")
    cat("  - 'Daily Heart Rate Variability Summary'\n")
    cat("  - 'Heart Rate Variability Details'\n")
  }
  cat("============================================================\n")
}

# ============================================================================
# MAIN EXECUTION
# ============================================================================
merge_all_hrv_files(base_folder, output_folder)
