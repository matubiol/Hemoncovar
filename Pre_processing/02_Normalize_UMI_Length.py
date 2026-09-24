#!/usr/bin/env python3
"""
Normalizes dual UMI segment lengths PER SAMPLE, only when necessary.

Some sequencing runs have inconsistent UMI segment lengths (e.g. mixing 6bp
and 7bp) within the same FASTQ, which breaks fgbio's GroupReadsByUmi/
CorrectUmis (both require all UMIs to be the same length). Other runs are
perfectly consistent (e.g. always 7bp+7bp), and truncating those would
throw away real UMI information.

This script:
  1. Scans R1 once to detect the distinct UMI segment lengths present.
  2. If lengths are already consistent, copies the FASTQs through unchanged
     (fast binary copy, no per-read rewriting, no information lost).
  3. If lengths are inconsistent, truncates each UMI segment independently
     to the minimum length observed FOR THAT SAMPLE (not a hardcoded value)
     -- so different batches/runs are each truncated only as much as their
     own data actually requires.
"""
import gzip
import shutil
from collections import Counter
from typing import Optional


def detect_umi_segment_lengths(fastq_path: str) -> Optional[tuple]:
    """Scan all read headers in a FASTQ and return
    (is_consistent, min_len1, min_len2, lengths_seen)."""
    lengths_seen = Counter()

    with gzip.open(fastq_path, "rt") as fin:
        for i, line in enumerate(fin):
            if i % 4 != 0:
                continue
            name = line.split(" ", 1)[0]
            name_fields = name.split(":")
            if len(name_fields) < 8:
                continue
            umi_segments = name_fields[7].split("+")
            if len(umi_segments) != 2:
                continue
            lengths_seen[(len(umi_segments[0]), len(umi_segments[1]))] += 1

    if not lengths_seen:
        return None  # No parseable UMI found -- nothing to normalize

    is_consistent = len(lengths_seen) == 1
    min_len1 = min(l1 for l1, _ in lengths_seen)
    min_len2 = min(l2 for _, l2 in lengths_seen)
    return is_consistent, min_len1, min_len2, dict(lengths_seen)


def truncate_fastq(in_path: str, out_path: str, len1: int, len2: int) -> None:
    with gzip.open(in_path, "rt") as fin, gzip.open(out_path, "wt") as fout:
        for i, line in enumerate(fin):
            if i % 4 == 0:
                parts = line.rstrip("\n").split(" ", 1)
                name = parts[0]
                rest = " " + parts[1] if len(parts) > 1 else ""
                name_fields = name.split(":")
                if len(name_fields) >= 8:
                    umi_segments = name_fields[7].split("+")
                    if len(umi_segments) == 2:
                        umi_segments[0] = umi_segments[0][:len1]
                        umi_segments[1] = umi_segments[1][:len2]
                        name_fields[7] = "+".join(umi_segments)
                        name = ":".join(name_fields)
                line = name + rest + "\n"
            fout.write(line)


def copy_fastq(in_path: str, out_path: str) -> None:
    """Fast binary copy -- used when no truncation is needed."""
    with gzip.open(in_path, "rb") as fin, gzip.open(out_path, "wb") as fout:
        shutil.copyfileobj(fin, fout)


if __name__ == "__main__":
    r1_in, r2_in = snakemake.input.r1, snakemake.input.r2
    r1_out, r2_out = snakemake.output.r1, snakemake.output.r2

    with open(snakemake.log[0], "w") as log:
        result = detect_umi_segment_lengths(r1_in)

        if result is None:
            log.write("No dual UMI found in read names -- copying through unchanged.\n")
            copy_fastq(r1_in, r1_out)
            copy_fastq(r2_in, r2_out)
        else:
            is_consistent, min_len1, min_len2, lengths_seen = result
            log.write(f"UMI segment lengths observed: {lengths_seen}\n")

            if is_consistent:
                log.write(
                    f"Lengths already consistent ({min_len1}+{min_len2}) -- "
                    "copying through unchanged, no truncation applied.\n"
                )
                copy_fastq(r1_in, r1_out)
                copy_fastq(r2_in, r2_out)
            else:
                log.write(
                    f"Inconsistent lengths detected -- truncating both "
                    f"segments to ({min_len1}+{min_len2}) for this sample.\n"
                )
                truncate_fastq(r1_in, r1_out, min_len1, min_len2)
                truncate_fastq(r2_in, r2_out, min_len1, min_len2)