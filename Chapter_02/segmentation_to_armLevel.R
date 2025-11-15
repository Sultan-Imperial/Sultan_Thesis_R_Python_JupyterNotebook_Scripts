#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(GenomicRanges)
  library(dplyr)
})

buildPQArms <- function(armDF) {
  # Simplified version without rowwise
  p_arms <- data.frame(
    chr = armDF$chr,
    arm = "p",
    arm_start = 1,
    arm_end = armDF$p_length
  )
  
  q_arms <- data.frame(
    chr = armDF$chr,
    arm = "q",
    arm_start = armDF$p_length + 1,
    arm_end = armDF$p_length + armDF$q_length
  )
  
  armDF_out <- rbind(p_arms, q_arms)
  
  arms_gr <- GRanges(
    seqnames = armDF_out$chr,
    ranges = IRanges(start=armDF_out$arm_start, end=armDF_out$arm_end),
    arm = armDF_out$arm
  )
  
  return(arms_gr)
}

computeArmLevelCNA <- function(segDF, armsDF) {
  # Add event classification
  segDF$event <- ifelse(segDF$seg.mean >= 0.3, "GAIN",
                        ifelse(segDF$seg.mean <= -0.3, "LOSS", "NEUT"))
  
  # Create GRanges objects
  seg_gr <- GRanges(
    seqnames = segDF$chr,
    ranges = IRanges(start=segDF$start, end=segDF$end),
    sample = segDF$sample,
    event = segDF$event,
    seg.mean = segDF$seg.mean,
    num.mark = segDF$num.mark
  )
  
  arms_gr <- buildPQArms(armsDF)
  
  # Find overlaps
  hits <- findOverlaps(seg_gr, arms_gr)
  
  # Calculate overlaps
  overlap_ranges <- pintersect(seg_gr[queryHits(hits)], arms_gr[subjectHits(hits)])
  
  # Create overlap dataframe
  overlap_df <- data.frame(
    sample = seg_gr$sample[queryHits(hits)],
    chr = as.character(seqnames(seg_gr)[queryHits(hits)]),
    arm = arms_gr$arm[subjectHits(hits)],
    seg_mean = seg_gr$seg.mean[queryHits(hits)],
    num_mark = seg_gr$num.mark[queryHits(hits)],
    event = seg_gr$event[queryHits(hits)],
    overlap_len = width(overlap_ranges),
    seg_len = width(seg_gr[queryHits(hits)])
  )
  
  # Calculate weighted summaries
  overlap_df$weight <- overlap_df$overlap_len / overlap_df$seg_len
  overlap_df$weighted_mean <- overlap_df$seg_mean * overlap_df$weight
  
  # Summarize by arm
  arm_summary <- overlap_df %>%
    group_by(sample, chr, arm) %>%
    summarize(
      mean_bpWeighted = sum(weighted_mean) / sum(weight),
      total_markers = sum(num_mark * weight),
      total_bpOverlap = sum(overlap_len),
      .groups = "drop"
    )
  
  # Calculate dominant event
  event_summary <- overlap_df %>%
    group_by(sample, chr, arm, event) %>%
    summarize(
      event_bp = sum(overlap_len),
      .groups = "drop"
    ) %>%
    group_by(sample, chr, arm) %>%
    slice_max(event_bp, n = 1) %>%
    ungroup()
  
  # Merge summaries
  final_summary <- arm_summary %>%
    left_join(
      event_summary %>% select(sample, chr, arm, event, event_bp),
      by = c("sample", "chr", "arm")
    ) %>%
    mutate(
      top_event = event,
      top_event_bp = event_bp,
      top_event_fraction = event_bp / total_bpOverlap
    ) %>%
    select(
      sample, chr, arm, 
      mean_bpWeighted, total_markers, total_bpOverlap,
      top_event, top_event_bp, top_event_fraction
    )
  
  return(as.data.frame(final_summary))
}

# Interactive usage example
if (interactive()) {
  # Example usage:
  # segDF <- read.delim("your_cbioportal_file.seg", header=TRUE)
  # armDF <- read.csv("pq_arm_boundaries.csv", header=TRUE)
  # results <- computeArmLevelCNA(segDF, armDF)
  # write.csv(results, "arm_level_results.csv", row.names=FALSE)
}
########Or if you prefer to run it interactively in R:#########

# Source the script
## source("cbioportal_to_armLevel.R")

# Read your input files
## segDF <- read.delim("your_cbioportal_file.seg", header=TRUE, stringsAsFactors=FALSE)
## armDF <- read.csv("pq_arm_boundaries.csv", header=TRUE, stringsAsFactors=FALSE)

# Run the analysis
## armRes <- computeArmLevelCNA(segDF, armDF)

# Write the results
## write.csv(armRes, "output_arm_level.csv", row.names=FALSE)
