#!/usr/bin/env Rscript
# 00_data_retrieval.R (v2)
# Phase 1: Download all multi-omics muscle data from NASA OSDR S3 + GEO + USDA

suppressPackageStartupMessages({
  library(httr)
  library(jsonlite)
  library(readxl)
  library(xml2)
})

data_dir <- "/workspace/osdr_pipeline/data"
dir.create(file.path(data_dir, "osdr_raw"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(data_dir, "osdr_metadata"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(data_dir, "geo_raw"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(data_dir, "usda"), showWarnings = FALSE, recursive = TRUE)

get_s3_listing <- function(osd_id, max_keys = 1000) {
  all_keys <- c()
  marker <- ""
  repeat {
    url <- sprintf("https://nasa-osdr.s3.amazonaws.com/?prefix=%s/&max-keys=%d&marker=%s",
                   osd_id, max_keys, URLencode(marker, reserved = TRUE))
    tryCatch({
      resp <- GET(url, timeout(30))
      content_xml <- content(resp, "text")
      doc <- read_xml(content_xml)
      xml_ns_strip(doc)
      keys <- xml_text(xml_find_all(doc, "//Key"))
      all_keys <- c(all_keys, keys)
      truncated <- xml_text(xml_find_all(doc, "//IsTruncated"))
      if (length(truncated) > 0 && truncated[1] == "true" && length(keys) > 0) {
        marker <- keys[length(keys)]
      } else {
        break
      }
    }, error = function(e) {
      cat(sprintf("  S3 listing failed for %s: %s\n", osd_id, conditionMessage(e)))
      break
    })
  }
  return(all_keys)
}

cat("=== Downloading OSDR datasets ===\n")

datasets_to_download <- c(
  "OSD-401", "OSD-21", "OSD-111", "OSD-135",
  "OSD-684", "OSD-530", "OSD-569", "OSD-571",
  "OSD-370", "OSD-195",
  "OSD-53",
  "OSD-488", "OSD-662", "OSD-660", "OSD-664",
  "OSD-652", "OSD-654", "OSD-661", "OSD-663",
  "OSD-691",
  "OSD-139", "OSD-353"
)

for (osd_id in datasets_to_download) {
  cat(sprintf("\n--- %s ---\n", osd_id))
  dest_dir <- file.path(data_dir, "osdr_raw", osd_id)
  dir.create(dest_dir, showWarnings = FALSE, recursive = TRUE)

  keys <- get_s3_listing(osd_id)
  if (length(keys) == 0) {
    cat(sprintf("  No files found for %s\n", osd_id))
    next
  }

  # Filter for processed files (xlsx, csv, txt, zip) - skip raw fastq, bam, tar
  processed_keys <- keys[grepl("\\.(xlsx|csv|txt|zip)$", keys, ignore.case = TRUE)]
  # Also include .results files (RSEM) for RNA-seq
  processed_keys <- c(processed_keys, keys[grepl("\\.results$", keys, ignore.case = TRUE)])
  # Exclude large bam files
  processed_keys <- processed_keys[!grepl("\\.bam", processed_keys, ignore.case = TRUE)]
  # Exclude fastq
  processed_keys <- processed_keys[!grepl("\\.fastq", processed_keys, ignore.case = TRUE)]

  cat(sprintf("  Total keys: %d, processed: %d\n", length(keys), length(processed_keys)))

  for (key in processed_keys) {
    url <- paste0("https://nasa-osdr.s3.amazonaws.com/", key)
    filename <- basename(key)
    if (filename == "") next
    dest <- file.path(dest_dir, filename)
    if (file.exists(dest) && file.info(dest)$size > 100) {
      cat(sprintf("  SKIP: %s\n", filename))
      next
    }
    tryCatch({
      head_resp <- HEAD(url, timeout(15))
      size <- as.numeric(headers(head_resp)$`content-length` %||% 0)
      if (!is.na(size) && size > 500e6) {
        cat(sprintf("  SKIP (too large %.0f MB): %s\n", size / 1e6, filename))
        next
      }
      download.file(url, dest, quiet = TRUE, mode = "wb")
      cat(sprintf("  OK: %s (%.1f MB)\n", filename, file.info(dest)$size / 1e6))
    }, error = function(e) {
      cat(sprintf("  FAIL: %s - %s\n", filename, conditionMessage(e)))
    })
  }
}

# ============================================================
# OSDR metadata via REST API
# ============================================================
cat("\n=== Downloading OSDR metadata ===\n")
for (osd_id in datasets_to_download) {
  url <- sprintf("https://visualization.osdr.nasa.gov/biodata/api/v2/dataset/%s/", osd_id)
  dest <- file.path(data_dir, "osdr_metadata", paste0(osd_id, "_metadata.json"))
  if (file.exists(dest) && file.info(dest)$size > 100) next
  tryCatch({
    resp <- GET(url, timeout(30))
    writeLines(content(resp, "text"), dest)
    cat(sprintf("  OK: %s\n", osd_id))
  }, error = function(e) cat(sprintf("  FAIL: %s\n", osd_id)))
}

# ============================================================
# GEO: GSE213808
# ============================================================
cat("\n=== Downloading GEO GSE213808 ===\n")
geo_dir <- file.path(data_dir, "geo_raw")
geo_files <- list(
  "GSE213808_TGB_022_3group_N6_totalcount_quantilenormalized.xlsx" = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE213nnn/GSE213808/suppl/GSE213808_TGB_022_3group_N6_totalcount_quantilenormalized.xlsx",
  "GSE213808_Mus_musculus_mm10_GeneAnnotation.xlsx" = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE213nnn/GSE213808/suppl/GSE213808_Mus_musculus_mm10_GeneAnnotation.xlsx"
)
for (fname in names(geo_files)) {
  dest <- file.path(geo_dir, fname)
  if (file.exists(dest) && file.info(dest)$size > 100) next
  tryCatch({
    download.file(geo_files[[fname]], dest, quiet = TRUE, mode = "wb")
    cat(sprintf("  OK: %s (%.1f MB)\n", fname, file.info(dest)$size / 1e6))
  }, error = function(e) cat(sprintf("  FAIL: %s\n", fname)))
}

# ============================================================
# USDA FoodData Central
# ============================================================
cat("\n=== Downloading USDA FoodData Central ===\n")
usda_dir <- file.path(data_dir, "usda")
# Try multiple URL patterns
usda_urls <- list(
  "sr_legacy.zip" = c(
    "https://fdc.nal.usda.gov/fdc-datasets/SR_Legacy_csv_2021-10-28.zip",
    "https://fdc.nal.usda.gov/fdc-datasets/SR_Legacy_csv_2024-04-19.zip"
  ),
  "foundation.zip" = c(
    "https://fdc.nal.usda.gov/fdc-datasets/Foundation_Food_csv_2021-10-28.zip",
    "https://fdc.nal.usda.gov/fdc-datasets/Foundation_Food_csv_2024-04-19.zip"
  )
)
for (fname in names(usda_urls)) {
  dest <- file.path(usda_dir, fname)
  if (file.exists(dest) && file.info(dest)$size > 100) next
  for (url in usda_urls[[fname]]) {
    tryCatch({
      download.file(url, dest, quiet = TRUE, mode = "wb")
      if (file.exists(dest) && file.info(dest)$size > 100) {
        cat(sprintf("  OK: %s (%.1f MB)\n", fname, file.info(dest)$size / 1e6))
        unzip(dest, exdir = file.path(usda_dir, sub("\\.zip$", "", fname)))
        break
      }
    }, error = function(e) {})
  }
  if (!file.exists(dest) || file.info(dest)$size < 100) {
    cat(sprintf("  FAIL: %s (all URLs failed)\n", fname))
  }
}

# ============================================================
# Summary
# ============================================================
cat("\n=== Download Summary ===\n")
for (osd_id in datasets_to_download) {
  dest_dir <- file.path(data_dir, "osdr_raw", osd_id)
  if (dir.exists(dest_dir)) {
    files <- list.files(dest_dir)
    total_size <- sum(file.size(file.path(dest_dir, files)), na.rm = TRUE)
    cat(sprintf("  %s: %d files, %.1f MB\n", osd_id, length(files), total_size / 1e6))
  }
}
cat("\nDone.\n")
