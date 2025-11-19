#!/usr/bin/env python3
"""
Author: Sultan N. Alharbi

Title: Cell-Free DNA Fragmentomics and Tumour Fraction as Diagnostic Biomarkers for Hepatocellular Carcinoma

PhD Thesis — Imperial College London, 2026

Description:
Fragment-size feature extraction pipeline for cfDNA sWGS data.

This script:

1. Runs Picard's CollectInsertSizeMetrics on per-sample BAM files
   (or consumes pre-computed Picard metrics files).
2. Parses insert size histograms.
3. Bins fragment lengths into user-defined ranges, e.g.
   P(1_30), P(31_60), ...
4. Computes per-sample proportions for each bin:
   proportion(bin) = count_in_bin / total_fragment_count.
5. Outputs a sample × feature matrix as CSV/TSV.

Author: Sultan N. Alharbi
Python: 3.10+
"""

from __future__ import annotations

import argparse
import csv
import logging
import os
import subprocess
from collections import Counter
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Tuple, Optional

import numpy as np
import pandas as pd

# ---------------------------------------------------------------------
# Logging configuration
# ---------------------------------------------------------------------

LOG = logging.getLogger(__name__)


def configure_logging(verbosity: int = 1) -> None:
    """Configure global logging level and format."""
    level = logging.INFO if verbosity == 1 else logging.DEBUG if verbosity > 1 else logging.WARNING
    logging.basicConfig(
        level=level,
        format="%(asctime)s [%(levelname)s] %(name)s - %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )


# ---------------------------------------------------------------------
# Data structures
# ---------------------------------------------------------------------

@dataclass
class InsertSizeBin:
    """A single fragment size bin definition."""
    start: int  # inclusive
    end: int    # inclusive

    @property
    def label(self) -> str:
        """Return a human-readable feature label, e.g. P(1_30)."""
        return f"P({self.start}_{self.end})"

    def contains(self, size: int) -> bool:
        """Return True if an insert size falls within this bin."""
        return self.start <= size <= self.end


# ---------------------------------------------------------------------
# Picard integration
# ---------------------------------------------------------------------

def run_picard_collect_insert_size_metrics(
    bam_path: Path,
    picard_jar: Path,
    metrics_path: Path,
    histogram_path: Path,
    java_opts: str = "-Xmx4g",
    picard_opts: Optional[List[str]] = None,
) -> None:
    """
    Run Picard CollectInsertSizeMetrics for a single BAM file.

    Parameters
    ----------
    bam_path : Path
        Input BAM file path.
    picard_jar : Path
        Path to Picard jar file.
    metrics_path : Path
        Path to insert size metrics output file (text).
    histogram_path : Path
        Path to histogram PDF output.
    java_opts : str
        JVM options, e.g. "-Xmx4g".
    picard_opts : list of str, optional
        Additional Picard options, e.g. ["M=0.5"].
    """
    picard_opts = picard_opts or []
    cmd = [
        "java",
        *java_opts.split(),
        "-jar",
        str(picard_jar),
        "CollectInsertSizeMetrics",
        f"I={bam_path}",
        f"O={metrics_path}",
        f"H={histogram_path}",
        "M=0.5",  # default minimum percent to include in histogram
        *picard_opts,
    ]
    LOG.info("Running Picard CollectInsertSizeMetrics: %s", " ".join(map(str, cmd)))
    subprocess.run(cmd, check=True)


def parse_picard_insert_metrics(metrics_path: Path) -> Dict[int, int]:
    """
    Parse Picard CollectInsertSizeMetrics metrics file and return a dict of insert_size -> count.

    The Picard metrics file contains a header and a HISTOGRAM section with tab-delimited columns:
    insert_size <TAB> count <TAB> ...

    Parameters
    ----------
    metrics_path : Path
        Path to Picard metrics text file.

    Returns
    -------
    size_counts : dict
        Mapping from insert_size (int) to count (int).
    """
    size_counts: Dict[int, int] = {}
    in_histogram = False

    with open(metrics_path, "r") as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith("#"):
                # Picard comments; skip.
                continue
            if line.startswith("## HISTOGRAM"):
                in_histogram = True
                # Next line will be header.
                continue
            if not in_histogram:
                # Still in metrics header; skip.
                continue

            # Histogram section: header row, then data rows.
            parts = line.split("\t")
            if parts[0] == "insert_size":
                # header row within histogram; skip
                continue
            try:
                insert_size = int(parts[0])
                count = int(parts[1])
            except (ValueError, IndexError):
                continue

            # Ignore non-positive sizes
            if insert_size <= 0:
                continue

            size_counts[insert_size] = size_counts.get(insert_size, 0) + count

    return size_counts


# ---------------------------------------------------------------------
# Binning logic
# ---------------------------------------------------------------------

def generate_bins(bin_size: int, max_size: int) -> List[InsertSizeBin]:
    """
    Create contiguous size bins from 1 to max_size inclusive, each of width bin_size.

    Example (bin_size=30, max_size=300) -> bins:
      [1–30], [31–60], [61–90], ..., [271–300]

    Parameters
    ----------
    bin_size : int
        Width of each bin (in base pairs).
    max_size : int
        Maximum insert size to consider. Sizes above this are ignored or can be
        handled separately.

    Returns
    -------
    bins : list of InsertSizeBin
    """
    if bin_size <= 0:
        raise ValueError("bin_size must be positive")
    if max_size <= 0:
        raise ValueError("max_size must be positive")

    bins: List[InsertSizeBin] = []
    start = 1
    while start <= max_size:
        end = min(start + bin_size - 1, max_size)
        bins.append(InsertSizeBin(start=start, end=end))
        start = end + 1
    return bins


def bin_insert_sizes(
    size_counts: Dict[int, int],
    bins: List[InsertSizeBin],
) -> Tuple[Dict[str, int], int]:
    """
    Bin insert sizes into predefined ranges.

    Parameters
    ----------
    size_counts : dict
        Mapping from insert_size (int) -> count (int).
    bins : list of InsertSizeBin
        Bin definitions.

    Returns
    -------
    bin_counts : dict
        Mapping from bin_label -> count.
    total_count : int
        Total number of fragments counted across all bins.
    """
    bin_counts = Counter()
    total_count = 0

    for size, count in size_counts.items():
        # Find first bin that contains this size.
        for b in bins:
            if b.contains(size):
                bin_counts[b.label] += count
                total_count += count
                break
        # Sizes > max_size are currently ignored; could be added as overflow bin.

    return dict(bin_counts), total_count


def compute_bin_proportions(
    bin_counts: Dict[str, int],
    total_count: int,
) -> Dict[str, float]:
    """
    Convert bin counts to proportions.

    proportion(bin) = count(bin) / total_count.

    Parameters
    ----------
    bin_counts : dict
        Bin label -> count.
    total_count : int
        Total count across bins.

    Returns
    -------
    proportions : dict
        Bin label -> proportion.
    """
    if total_count <= 0:
        # Avoid division by zero; return zeros for all bins.
        return {label: 0.0 for label in bin_counts.keys()}

    return {label: count / total_count for label, count in bin_counts.items()}


# ---------------------------------------------------------------------
# High-level per-sample processing
# ---------------------------------------------------------------------

def sample_name_from_bam(bam_path: Path) -> str:
    """Derive a sample ID from a BAM filename (stem without extension)."""
    return bam_path.stem


def process_sample_from_bam(
    bam_path: Path,
    picard_jar: Path,
    bin_size: int,
    max_size: int,
    temp_dir: Optional[Path] = None,
    keep_metrics: bool = False,
) -> Tuple[str, Dict[str, float]]:
    """
    Run Picard on a BAM and compute fragment-size bin proportions.

    Parameters
    ----------
    bam_path : Path
        Input BAM.
    picard_jar : Path
        Picard jar path.
    bin_size : int
        Bin size (bp).
    max_size : int
        Max size considered; sizes > max_size are ignored.
    temp_dir : Path, optional
        Directory to write Picard metrics + histogram files.
    keep_metrics : bool
        If True, keep Picard metrics and histogram files; otherwise remove.

    Returns
    -------
    sample_id : str
        Sample identifier.
    proportions : dict
        Feature_name -> proportion.
    """
    sample_id = sample_name_from_bam(bam_path)
    temp_dir = temp_dir or Path("insert_size_metrics")
    temp_dir.mkdir(parents=True, exist_ok=True)

    metrics_path = temp_dir / f"{sample_id}_insert_metrics.txt"
    hist_path = temp_dir / f"{sample_id}_insert_histogram.pdf"

    run_picard_collect_insert_size_metrics(
        bam_path=bam_path,
        picard_jar=picard_jar,
        metrics_path=metrics_path,
        histogram_path=hist_path,
    )

    size_counts = parse_picard_insert_metrics(metrics_path)
    bins = generate_bins(bin_size=bin_size, max_size=max_size)
    bin_counts, total = bin_insert_sizes(size_counts, bins)
    proportions = compute_bin_proportions(bin_counts, total)

    if not keep_metrics:
        try:
            metrics_path.unlink()
        except FileNotFoundError:
            pass
        try:
            hist_path.unlink()
        except FileNotFoundError:
            pass

    return sample_id, proportions


def process_sample_from_metrics(
    metrics_path: Path,
    bin_size: int,
    max_size: int,
) -> Tuple[str, Dict[str, float]]:
    """
    Compute bin proportions given an existing Picard metrics file.

    This is useful for testing or when Picard has already been run.

    Parameters
    ----------
    metrics_path : Path
        Picard metrics file path.
    bin_size : int
        Bin size.
    max_size : int
        Max size considered.

    Returns
    -------
    sample_id : str
        Sample identifier (derived from filename).
    proportions : dict
        Feature_name -> proportion.
    """
    sample_id = metrics_path.stem.replace("_insert_metrics", "")
    size_counts = parse_picard_insert_metrics(metrics_path)
    bins = generate_bins(bin_size=bin_size, max_size=max_size)
    bin_counts, total = bin_insert_sizes(size_counts, bins)
    proportions = compute_bin_proportions(bin_counts, total)
    return sample_id, proportions


# ---------------------------------------------------------------------
# Feature matrix construction
# ---------------------------------------------------------------------

def build_feature_matrix(
    sample_to_features: Dict[str, Dict[str, float]]
) -> pd.DataFrame:
    """
    Convert sample -> feature mapping into a DataFrame.

    Rows: samples.
    Columns: feature names (bin labels).
    """
    df = pd.DataFrame.from_dict(sample_to_features, orient="index")
    df.index.name = "sample_id"
    df = df.sort_index(axis=0).sort_index(axis=1)
    return df


# ---------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------

def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description=(
            "Extract cfDNA fragment-size features from Picard CollectInsertSizeMetrics.\n"
            "Can either run Picard on BAMs or consume existing metrics files."
        )
    )
    subparsers = parser.add_subparsers(dest="mode", required=True)

    # Mode 1: BAM input
    p_bam = subparsers.add_parser(
        "from-bam",
        help="Run Picard on BAM files and compute fragment-size features.",
    )
    p_bam.add_argument(
        "--bam-list",
        required=True,
        help="Text file with one BAM path per line.",
    )
    p_bam.add_argument(
        "--picard-jar",
        required=True,
        help="Path to Picard jar file.",
    )
    p_bam.add_argument(
        "--bin-size",
        type=int,
        default=30,
        help="Size of fragment-length bins in base pairs (default: 30).",
    )
    p_bam.add_argument(
        "--max-size",
        type=int,
        default=600,
        help="Maximum insert size to consider (default: 600 bp).",
    )
    p_bam.add_argument(
        "--temp-dir",
        type=str,
        default="insert_size_metrics",
        help="Directory to store Picard metrics and histograms.",
    )
    p_bam.add_argument(
        "--keep-metrics",
        action="store_true",
        help="Keep Picard metrics/histogram files instead of deleting them.",
    )
    p_bam.add_argument(
        "--output",
        required=True,
        help="Output CSV file for fragment-size feature matrix.",
    )
    p_bam.add_argument(
        "-v",
        "--verbose",
        action="count",
        default=1,
        help="Increase log verbosity level.",
    )

    # Mode 2: pre-computed metrics
    p_met = subparsers.add_parser(
        "from-metrics",
        help="Compute fragment-size features from existing Picard metrics files.",
    )
    p_met.add_argument(
        "--metrics-list",
        required=True,
        help="Text file with one Picard metrics path per line.",
    )
    p_met.add_argument(
        "--bin-size",
        type=int,
        default=30,
        help="Size of fragment-length bins in base pairs (default: 30).",
    )
    p_met.add_argument(
        "--max-size",
        type=int,
        default=600,
        help="Maximum insert size to consider (default: 600 bp).",
    )
    p_met.add_argument(
        "--output",
        required=True,
        help="Output CSV file for fragment-size feature matrix.",
    )
    p_met.add_argument(
        "-v",
        "--verbose",
        action="count",
        default=1,
        help="Increase log verbosity level.",
    )

    return parser.parse_args()


def main() -> None:
    args = parse_args()
    configure_logging(args.verbose)

    if args.mode == "from-bam":
        bam_list_path = Path(args.bam_list)
        picard_jar = Path(args.picard_jar)
        temp_dir = Path(args.temp_dir)

        with open(bam_list_path) as fh:
            bam_files = [Path(line.strip()) for line in fh if line.strip()]

        LOG.info("Found %d BAM files to process.", len(bam_files))

        sample_to_features: Dict[str, Dict[str, float]] = {}
        for bam_path in bam_files:
            if not bam_path.exists():
                LOG.error("BAM file does not exist: %s", bam_path)
                continue
            sample_id, feats = process_sample_from_bam(
                bam_path=bam_path,
                picard_jar=picard_jar,
                bin_size=args.bin_size,
                max_size=args.max_size,
                temp_dir=temp_dir,
                keep_metrics=args.keep_metrics,
            )
            sample_to_features[sample_id] = feats

    elif args.mode == "from-metrics":
        metrics_list_path = Path(args.metrics_list)
        with open(metrics_list_path) as fh:
            metrics_files = [Path(line.strip()) for line in fh if line.strip()]

        LOG.info("Found %d metrics files to process.", len(metrics_files))

        sample_to_features = {}
        for metrics_path in metrics_files:
            if not metrics_path.exists():
                LOG.error("Metrics file does not exist: %s", metrics_path)
                continue
            sample_id, feats = process_sample_from_metrics(
                metrics_path=metrics_path,
                bin_size=args.bin_size,
                max_size=args.max_size,
            )
            sample_to_features[sample_id] = feats
    else:
        raise RuntimeError(f"Unknown mode: {args.mode}")

    df = build_feature_matrix(sample_to_features)
    output_path = Path(args.output)
    df.to_csv(output_path, index=True)
    LOG.info("Saved fragment-size feature matrix to %s", output_path)


if __name__ == "__main__":
    main()
