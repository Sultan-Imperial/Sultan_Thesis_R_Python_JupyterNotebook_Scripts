#!/usr/bin/env Rscript

## =============================================================
## Arm-level CNA aggregation from segment-level CNA data
##
## - Input 1: segment-level CNA file from Data Portal (cBioPortal)
## - Input 2: chromosome arm boundary file (p/q lengths per chr)
## - Output : arm-level CNA summary (one row per sample × chr × arm)
##
## Author: Sultan N. Alharbi
## =============================================================

suppressPackageStartupMessages({
  library(optparse)
  library(GenomicRanges)
  library(dplyr)
})

## -------------------------------------------------------------
## Function: buildPQArms
## -------------------------------------------------------------

#' Build chromosome p and q arm ranges from arm-length table
#'
#' @param arm_df A data.frame containing at least:
#'   - chr: chromosome name (e.g. "chr1")
#'   - p_length: length of the p arm in base pairs
#'   - q_length: length of the q arm in base pairs
#'
#' @return A GRanges object with one range per chromosome arm
#'   (p and q), with metadata column `arm` ("p" or "q").
#'
buildPQArms <- function(arm_df) {
  # p arms: from 1 to p_length
  p_arms <- data.frame(
    chr       = arm_df$chr,
    arm       = "p",
    arm_start = 1L,
    arm_end   = arm_df$p_length
  )
  
  # q arms: from p_length + 1 to p_length + q_length
  q_arms <- data.frame(
    chr       = arm_df$chr,
    arm       = "q",
    arm_start = arm_df$p_length + 1L,
    arm_end   = arm_df$p_length + arm_df$q_length
  )
  
  arm_df_out <- rbind(p_arms, q_arms)
  
  GRanges(
    seqnames = arm_df_out$chr,
    ranges   = IRanges(start = arm_df_out$arm_start,
                       end   = arm_df_out$arm_end),
    arm      = arm_df_out$arm
  )
}

## -------------------------------------------------------------
## Function: compute Arm-Level CNA
## -------------------------------------------------------------

#' Compute arm-level CNA from segment-level CNA data
#'
#' @param seg_df A data.frame of segment-level CNA data containing:
#'   - sample: sample identifier
#'   - chr: chromosome (e.g. "chr1")
#'   - start, end: genomic segment coordinates (1-based)
#'   - seg.mean: log2 copy-number ratio for the segment
#'   - num.mark: number of markers/probes in the segment
#'
#' @param arms_df A data.frame with chromosome arm lengths, as
#'   required by \code{buildPQArms()}:
#'   - chr, p_length, q_length
#'
#' @return A data.frame with one row per sample × chr × arm,
#'   containing:
#'   - sample, chr, arm
#'   - mean_bpWeighted     (length-weighted CNA)
#'   - total_markers       (weighted sum of num.mark)
#'   - total_bpOverlap     (total bp overlapping the arm)
#'   - top_event           ("GAIN", "LOSS", or "NEUT")
#'   - top_event_bp        (bp supporting the dominant event)
#'   - top_event_fraction  (top_event_bp / total_bpOverlap)
#'
computeArmLevelCNA <- function(seg_df, arms_df) {
  
  # -----------------------------------------------------------
  # 1. Classify segments as GAIN / LOSS / NEUTRAL
  # -----------------------------------------------------------
  seg_df$event <- ifelse(
    seg_df$seg.mean >=  0.3, "GAIN",
    ifelse(seg_df$seg.mean <= -0.3, "LOSS", "NEUT")
  )
  
  # -----------------------------------------------------------
  # 2. Convert segments and arms to GRanges
  # -----------------------------------------------------------
  seg_gr <- GRanges(
    seqnames = seg_df$chr,
    ranges   = IRanges(start = seg_df$start,
                       end   = seg_df$end),
    sample   = seg_df$sample,
    event    = seg_df$event,
    seg.mean = seg_df$seg.mean,
    num.mark = seg_df$num.mark
  )
  
  arms_gr <- buildPQArms(arms_df)
  
  # -----------------------------------------------------------
  # 3. Find and quantify overlaps between segments and arms
  # -----------------------------------------------------------
  hits <- findOverlaps(seg_gr, arms_gr)
  
  overlap_ranges <- pintersect(
    seg_gr[queryHits(hits)],
    arms_gr[subjectHits(hits)]
  )
  
  overlap_df <- data.frame(
    sample      = seg_gr$sample[queryHits(hits)],
    chr         = as.character(seqnames(seg_gr)[queryHits(hits)]),
    arm         = arms_gr$arm[subjectHits(hits)],
    seg_mean    = seg_gr$seg.mean[queryHits(hits)],
    num_mark    = seg_gr$num.mark[queryHits(hits)],
    event       = seg_gr$event[queryHits(hits)],
    overlap_len = width(overlap_ranges),
    seg_len     = width(seg_gr[queryHits(hits)])
  )
  
  # -----------------------------------------------------------
  # 4. Compute length-weighted segment contributions per arm
  # -----------------------------------------------------------
  overlap_df$weight        <- overlap_df$overlap_len / overlap_df$seg_len
  overlap_df$weighted_mean <- overlap_df$seg_mean * overlap_df$weight
  
  # Arm-level summary of CNA
  arm_summary <- overlap_df %>%
    dplyr::group_by(sample, chr, arm) %>%
    dplyr::summarise(
      mean_bpWeighted = sum(weighted_mean) / sum(weight),
      total_markers   = sum(num_mark * weight),
      total_bpOverlap = sum(overlap_len),
      .groups         = "drop"
    )
  
  # -----------------------------------------------------------
  # 5. Determine dominant event per arm (GAIN/LOSS/NEUT)
  # -----------------------------------------------------------
  event_summary <- overlap_df %>%
    dplyr::group_by(sample, chr, arm, event) %>%
    dplyr::summarise(
      event_bp = sum(overlap_len),
      .groups  = "drop"
    ) %>%
    dplyr::group_by(sample, chr, arm) %>%
    dplyr::slice_max(event_bp, n = 1, with_ties = FALSE) %>%
    dplyr::ungroup()
  
  # -----------------------------------------------------------
  # 6. Merge summaries and compute dominant event fraction
  # -----------------------------------------------------------
  final_summary <- arm_summary %>%
    dplyr::left_join(
      event_summary %>% dplyr::select(sample, chr, arm, event, event_bp),
      by = c("sample", "chr", "arm")
    ) %>%
    dplyr::mutate(
      top_event          = event,
      top_event_bp       = event_bp,
      top_event_fraction = event_bp / total_bpOverlap
    ) %>%
    dplyr::select(
      sample, chr, arm,
      mean_bpWeighted, total_markers, total_bpOverlap,
      top_event, top_event_bp, top_event_fraction
    )
  
  as.data.frame(final_summary)
}

## -------------------------------------------------------------
## Command-line interface
## -------------------------------------------------------------

option_list <- list(
  make_option(
    c("-s", "--segments"),
    type    = "character",
    help    = "Path to segment-level CNA file (e.g. cBioPortal .seg or .txt)",
    metavar = "FILE"
  ),
  make_option(
    c("-a", "--arms"),
    type    = "character",
    help    = "Path to chromosome arm boundary file (e.g. pq_arm_boundaries.csv)",
    metavar = "FILE"
  ),
  make_option(
    c("-o", "--output"),
    type    = "character",
    default = "arm_level_results.csv",
    help    = "Output CSV file for arm-level CNA summary [default: %default]",
    metavar = "FILE"
  )
)

# Only parse options / run main if script is called non-interactively
if (!interactive()) {
  opt <- parse_args(OptionParser(option_list = option_list))
  
  if (is.null(opt$segments) || is.null(opt$arms)) {
    stop("Please provide both --segments and --arms file paths.", call. = FALSE)
  }
  
  message("Reading segment file: ", opt$segments)
  seg_df <- read.delim(opt$segments, header = TRUE, stringsAsFactors = FALSE)
  
  message("Reading arm boundary file: ", opt$arms)
  arms_df <- read.csv(opt$arms, header = TRUE, stringsAsFactors = FALSE)
  
  message("Computing arm-level CNA summary...")
  arm_res <- computeArmLevelCNA(seg_df, arms_df)
  
  message("Writing results to: ", opt$output)
  write.csv(arm_res, opt$output, row.names = FALSE)
  
  message("Done.")
}

## -------------------------------------------------------------
## Example interactive usage (for R console / RStudio)
## -------------------------------------------------------------
# In an interactive R session:
#
#   source("segmentation_to_armLevel.R")
#
#   seg_df  <- read.delim("TCGA_combined_study_segments_233_samples.seg",
#                         header = TRUE, stringsAsFactors = FALSE)
#   arms_df <- read.csv("pq_arm_boundaries.csv",
#                       header = TRUE, stringsAsFactors = FALSE)
#
#   arm_res <- computeArmLevelCNA(seg_df, arms_df)
#   write.csv(arm_res, "arm_level_results.csv", row.names = FALSE)
#
# From the command line:
#
#   Rscript compute_arm_level_cna.R \
#       --segments TCGA_combined_study_segments_233_samples.seg \
#       --arms pq_arm_boundaries.csv \
#       --output arm_level_results.csv
##
## =============================================================
