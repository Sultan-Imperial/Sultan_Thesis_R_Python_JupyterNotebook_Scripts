#!/usr/bin/env python

"""
Author: Sultan N. Alharbi

Title: Cell-Free DNA Fragmentomics and Tumour Fraction as Diagnostic Biomarkers for Hepatocellular Carcinoma

PhD Thesis — Imperial College London, 2026

Description:
cfDNA 5' end motif frequency pipeline.

This script reads one or more aligned BAM files and computes the frequency
of all possible k-mer motifs (default: 4-mer) at the 5' ends of cfDNA fragment
ends.

Key features
------------
- Handles BAM input.
- Filters reads by mapping quality and SAM flags.
- Correctly identifies the 5' end of each aligned fragment end:
    * Forward-strand reads: first aligned base.
    * Reverse-strand reads: last aligned base; the motif is reported as the
      reverse complement of the aligned query sequence so that motifs are always
      oriented 5'→3' with respect to the physical fragment.
- Ignores soft-clipped bases at read ends (uses first/last aligned base).
- Produces counts and normalised frequencies for all 4^k motifs.
- Optional multiprocessing across samples.
"""

import argparse
import itertools
import multiprocessing as mp
from collections import Counter
from pathlib import Path
from typing import Dict, Iterable, List, Tuple

import numpy as np
import pandas as pd
import pysam

VALID_BASES = ("A", "C", "G", "T")


# ----------------------------------------------------------------------
# Utility functions
# ----------------------------------------------------------------------
def generate_all_motifs(k: int = 4) -> List[str]:
    """Generate all possible k-mer motifs over A,C,G,T in lexicographic order."""
    return ["".join(p) for p in itertools.product(VALID_BASES, repeat=k)]


def reverse_complement(seq: str) -> str:
    """Return the reverse complement of a DNA sequence."""
    comp = str.maketrans("ACGTacgt", "TGCAtgca")
    return seq.translate(comp)[::-1]


def is_valid_dna(sequence: str) -> bool:
    """Check that a motif is composed only of canonical A/C/G/T bases."""
    return len(sequence) > 0 and all(base in VALID_BASES for base in sequence)


# ----------------------------------------------------------------------
# Core motif extraction logic
# ----------------------------------------------------------------------
def iter_5prime_motifs_from_bam(
    bam_path: Path,
    motif_length: int = 4,
    min_mapq: int = 30,
    include_duplicates: bool = False,
    include_secondary: bool = False,
    include_supplementary: bool = False,
) -> Iterable[str]:
    """
    Yield 5' end motifs (as 5'→3' fragment-oriented sequences) from an aligned BAM file.

    Parameters
    ----------
    bam_path : Path
        Input BAM or  file.
    motif_length : int
        Length of motif to extract (typically 4).
    min_mapq : int
        Minimum mapping quality (MAPQ) to keep a read.
    include_duplicates : bool
        If False, reads flagged as PCR/optical duplicates are excluded.
    include_secondary : bool
        If False, secondary alignments are excluded.
    include_supplementary : bool
        If False, supplementary alignments are excluded.

    Yields
    ------
    motif : str
        The k-mer motif at the 5' fragment end, oriented 5'→3'.
    """
    bam_path = Path(bam_path)
    mode = "rb"  # autodetect BAM
    with pysam.AlignmentFile(bam_path, mode) as bam:
        for read in bam.fetch(until_eof=True):
            # 1) Basic filtering
            if read.is_unmapped:
                continue
            if read.mapping_quality < min_mapq:
                continue
            if read.is_qcfail:
                continue
            if not include_duplicates and read.is_duplicate:
                continue
            if not include_secondary and read.is_secondary:
                continue
            if not include_supplementary and read.is_supplementary:
                continue
            if read.query_sequence is None:
                continue

            # Require enough aligned bases to cover the motif
            if read.query_alignment_length is None or read.query_alignment_length < motif_length:
                continue

            # 2) Extract motif in query space, accounting for strand and soft clipping
            if not read.is_reverse:
                # Forward-strand: 5' end = first aligned base
                start_idx = read.query_alignment_start
                end_idx = start_idx + motif_length
                motif = read.query_sequence[start_idx:end_idx]
            else:
                # Reverse-strand: 5' end = last aligned base
                end_idx = read.query_alignment_end
                start_idx = end_idx - motif_length
                if start_idx < 0:
                    continue
                subseq = read.query_sequence[start_idx:end_idx]
                motif = reverse_complement(subseq)

            if len(motif) != motif_length:
                # Can occur near contig boundaries or very short alignments
                continue
            if not is_valid_dna(motif):
                continue

            yield motif


def count_end_motifs(
    bam_path: Path,
    motif_length: int = 4,
    min_mapq: int = 30,
    include_duplicates: bool = False,
    include_secondary: bool = False,
    include_supplementary: bool = False,
) -> Tuple[str, Counter]:
    """
    Count 5' end motifs for a single BAM file.

    Returns
    -------
    sample_name : str
        Sample identifier inferred from file name stem.
    motif_counts : collections.Counter
        Counter mapping motif -> count for this sample.
    """
    bam_path = Path(bam_path)
    sample_name = bam_path.stem

    motifs = iter_5prime_motifs_from_bam(
        bam_path=bam_path,
        motif_length=motif_length,
        min_mapq=min_mapq,
        include_duplicates=include_duplicates,
        include_secondary=include_secondary,
        include_supplementary=include_supplementary,
    )
    motif_counts = Counter(motifs)
    return sample_name, motif_counts


def counts_to_dataframe(
    sample_counts: Dict[str, Counter],
    motif_length: int = 4,
) -> pd.DataFrame:
    """
    Convert per-sample motif Counters to a tidy DataFrame with counts and frequencies.

    Parameters
    ----------
    sample_counts : dict
        Mapping sample_name -> Counter(motif -> count).
    motif_length : int
        Length of motifs (e.g. 4).

    Returns
    -------
    df : pd.DataFrame
        Columns:
        - sample
        - motif
        - count
        - frequency (proportion of ends with that motif)
    """
    all_motifs = generate_all_motifs(motif_length)
    records = []

    for sample, counter in sample_counts.items():
        total = sum(counter.values())
        if total == 0:
            for motif in all_motifs:
                records.append(
                    {"sample": sample, "motif": motif, "count": 0, "frequency": 0.0}
                )
            continue

        for motif in all_motifs:
            c = counter.get(motif, 0)
            freq = c / total
            records.append(
                {"sample": sample, "motif": motif, "count": c, "frequency": freq}
            )

    df = pd.DataFrame.from_records(records)
    return df


def process_single_bam(args) -> Tuple[str, Counter]:
    """Helper for multiprocessing: unpack tuple of arguments."""
    return count_end_motifs(*args)


def run_pipeline(
    bam_files: List[str],
    motif_length: int = 4,
    min_mapq: int = 30,
    include_duplicates: bool = False,
    include_secondary: bool = False,
    include_supplementary: bool = False,
    n_workers: int = 1,
    output_path: str = "end_motif_frequencies.tsv",
) -> pd.DataFrame:
    """
    Run the end motif counting pipeline on a list of BAM files.

    Parameters
    ----------
    bam_files : list of str
        Paths to input BAM files.
    motif_length : int
        Length of motif (typically 4).
    min_mapq : int
        Minimum mapping quality.
    include_duplicates : bool
        Include duplicate-marked reads if True.
    include_secondary : bool
        Include secondary alignments if True.
    include_supplementary : bool
        Include supplementary alignments if True.
    n_workers : int
        Number of parallel worker processes. Use 1 to disable multiprocessing.
    output_path : str
        Output TSV file where combined counts and frequencies are written.

    Returns
    -------
    df : pd.DataFrame
        Combined motif counts and frequencies for all samples.
    """
    bam_files = [str(Path(b)) for b in bam_files]

    if n_workers is None or n_workers < 1:
        n_workers = 1

    sample_counts: Dict[str, Counter] = {}

    worker_args = [
        (
            Path(b),
            motif_length,
            min_mapq,
            include_duplicates,
            include_secondary,
            include_supplementary,
        )
        for b in bam_files
    ]

    if n_workers == 1:
        for args in worker_args:
            sample, counts = process_single_bam(args)
            sample_counts[sample] = counts
    else:
        with mp.Pool(processes=n_workers) as pool:
            for sample, counts in pool.imap_unordered(process_single_bam, worker_args):
                sample_counts[sample] = counts

    df = counts_to_dataframe(sample_counts, motif_length=motif_length)

    # Sanity check: counts sum to total ends per sample
    check_df = df.groupby("sample")["count"].sum().reset_index(name="total_count")
    for _, row in check_df.iterrows():
        sample = row["sample"]
        expected = sum(sample_counts[sample].values())
        assert (
            row["total_count"] == expected
        ), f"Count mismatch for sample {sample}: df={row['total_count']}, counter={expected}"

    # Save to TSV
    df.to_csv(output_path, sep="\t", index=False)
    return df


# ----------------------------------------------------------------------
# CLI and tests
# ----------------------------------------------------------------------
def parse_args() -> argparse.Namespace:
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description=(
            "Compute 5' end motif counts and frequencies from BAM files.\n"
            "Motifs are reported in the 5'→3' orientation of the fragment."
        )
    )
    parser.add_argument(
        "-i",
        "--input",
        nargs="+",
        required=False,
        help="One or more BAM files.",
    )
    parser.add_argument(
        "-l",
        "--bam-list",
        type=str,
        required=False,
        help="Optional text file with one BAM path per line.",
    )
    parser.add_argument(
        "-k",
        "--motif-length",
        type=int,
        default=4,
        help="Motif length (default: 4).",
    )
    parser.add_argument(
        "--min-mapq",
        type=int,
        default=30,
        help="Minimum mapping quality to include a read (default: 30).",
    )
    parser.add_argument(
        "--include-duplicates",
        action="store_true",
        help="Include reads marked as PCR/optical duplicates.",
    )
    parser.add_argument(
        "--include-secondary",
        action="store_true",
        help="Include secondary alignments.",
    )
    parser.add_argument(
        "--include-supplementary",
        action="store_true",
        help="Include supplementary alignments.",
    )
    parser.add_argument(
        "-o",
        "--output",
        type=str,
        default="end_motif_frequencies.tsv",
        help="Output TSV file (default: end_motif_frequencies.tsv).",
    )
    parser.add_argument(
        "-j",
        "--jobs",
        type=int,
        default=1,
        help="Number of parallel worker processes (default: 1).",
    )
    parser.add_argument(
        "--run-tests",
        action="store_true",
        help="Run built-in unit tests on a small synthetic BAM file and exit.",
    )

    return parser.parse_args()


# -------------------------------
# Testing utilities
# -------------------------------
def create_synthetic_bam(bam_path: Path, motif_length: int = 4) -> None:
    """
    Create a small synthetic BAM file for testing 5' end motif extraction.

    Contains:
    - One forward read with known 5' motif.
    - One reverse read where the motif is obtained by reverse complement.

    All reads are mapped with high MAPQ and no soft clipping.
    """
    header = {
        "HD": {"VN": "1.0"},
        "SQ": [{"LN": 1000, "SN": "chr1"}],
    }

    with pysam.AlignmentFile(bam_path, "wb", header=header) as bam:
        # Forward read: sequence "AAAACCCC" -> 5' motif should be "AAAA"
        a = pysam.AlignedSegment()
        a.query_name = "read_fwd"
        a.query_sequence = "AAAACCCC"
        a.flag = 0  # mapped, forward
        a.reference_id = 0
        a.reference_start = 100
        a.mapping_quality = 60
        a.cigar = ((0, len(a.query_sequence)),)  # 0 = M
        a.next_reference_id = -1
        a.next_reference_start = -1
        a.template_length = 0
        a.query_qualities = pysam.qualitystring_to_array("FFFFFFFF")
        bam.write(a)

        # Reverse read: sequence "GGGGTTTT" aligned to reverse strand.
        # Its 5' end motif should be reverse complement of last k bases.
        b = pysam.AlignedSegment()
        b.query_name = "read_rev"
        b.query_sequence = "GGGGTTTT"
        b.flag = 16  # mapped, reverse strand
        b.reference_id = 0
        b.reference_start = 200
        b.mapping_quality = 60
        b.cigar = ((0, len(b.query_sequence)),)
        b.next_reference_id = -1
        b.next_reference_start = -1
        b.template_length = 0
        b.query_qualities = pysam.qualitystring_to_array("FFFFFFFF")
        bam.write(b)

    pysam.index(str(bam_path))


def run_unit_tests(tmp_dir: Path) -> None:
    """Run a minimal set of unit tests to validate motif extraction."""
    bam_path = tmp_dir / "synthetic.bam"
    create_synthetic_bam(bam_path)

    sample, counts = count_end_motifs(
        bam_path=bam_path,
        motif_length=4,
        min_mapq=30,
        include_duplicates=False,
        include_secondary=False,
        include_supplementary=False,
    )

    assert sample == "synthetic", f"Unexpected sample name: {sample}"

    # Expected motifs:
    # read_fwd: "AAAACCCC" (forward) -> motif "AAAA"
    # read_rev: "GGGGTTTT" (reverse) -> last 4 bases "TTTT", reverse complement "AAAA"
    assert counts["AAAA"] == 2, f"Expected 2 occurrences of 'AAAA', found {counts['AAAA']}"
    total = sum(counts.values())
    assert total == 2, f"Expected 2 total motifs, found {total}"

    print("All unit tests passed.")


def main() -> None:
    args = parse_args()

    if args.run_tests:
        tmp_dir = Path("test_output")
        tmp_dir.mkdir(exist_ok=True)
        run_unit_tests(tmp_dir)
        return

    bam_files: List[str] = []
    if args.input:
        bam_files.extend(args.input)
    if args.bam_list:
        with open(args.bam_list) as f:
            bam_files.extend(line.strip() for line in f if line.strip())

    if not bam_files:
        raise SystemExit("No BAM files provided. Use --input or --bam-list.")

    df = run_pipeline(
        bam_files=bam_files,
        motif_length=args.motif_length,
        min_mapq=args.min_mapq,
        include_duplicates=args.include_duplicates,
        include_secondary=args.include_secondary,
        include_supplementary=args.include_supplementary,
        n_workers=args.jobs,
        output_path=args.output,
    )

    print(
        f"Written motif counts and frequencies for {df['sample'].nunique()} samples "
        f"to {args.output}"
    )


if __name__ == "__main__":
    main()
