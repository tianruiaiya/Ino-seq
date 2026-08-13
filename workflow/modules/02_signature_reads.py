#!/usr/bin/env python3
"""Ino-seq Step 01 mutation/read statistics.

Select ABE-associated reads from a UMI-deduplicated BAM, generate read-level
mutation statistics, and write the selected reads to an indexed BAM file.
"""

import argparse
import sys

import pysam


def calculate_max_consecutive_bases(sequence):
    """Return the maximum run length of an identical base in *sequence*."""
    if not sequence:
        return 0

    max_count = 1
    current_count = 1
    for i in range(1, len(sequence)):
        if sequence[i] == sequence[i - 1]:
            current_count += 1
        else:
            max_count = max(max_count, current_count)
            current_count = 1
    return max(max_count, current_count)


def collect_mutations(read, fasta, reference_base):
    """Collect substitutions of the specified reference base in CIGAR M blocks."""
    ref_start = read.reference_start
    ref_end = read.reference_end
    query_seq = read.query_sequence
    ref_seq = fasta.fetch(reference=read.reference_name, start=ref_start, end=ref_end)

    mutations = []
    ref_pos = ref_start
    query_pos = 0
    for op, length in read.cigartuples:
        if op == 0:  # M: match or mismatch
            for i in range(length):
                rbase_i = ref_seq[ref_pos - ref_start + i]
                qbase_i = query_seq[query_pos + i]
                if rbase_i.upper() == reference_base and qbase_i.upper() != reference_base:
                    mutations.append(f"{ref_pos + i + 1}:{rbase_i}>{qbase_i}")
            ref_pos += length
            query_pos += length
    return mutations


def write_selected_read(stats, outbam, read, fasta, location, direction, rbase, qbase, reference_base):
    mutations = collect_mutations(read, fasta, reference_base)
    mutation_str = ";".join(mutations) if mutations else "None"
    max_consecutive = calculate_max_consecutive_bases(read.query_sequence)

    flank_start = max(0, location - 21)
    flank_end = location + 20
    flank_seq = fasta.fetch(reference=read.reference_name, start=flank_start, end=flank_end)

    # Convert the candidate site to the workflow reporting coordinate.
    reported_location = location + 1 if direction == "-" else location - 1
    stats.write(
        f"{read.query_name}\t{read.reference_name}\t{read.reference_start + 1}\t{read.reference_end}\t"
        f"{reported_location}\t{direction}\t{rbase}\t{qbase}\t{mutation_str}\t{flank_seq}\t{max_consecutive}\n"
    )
    outbam.write(read)


def run(output_prefix, bam_path, reference_fasta, location_qual_threshold=30):
    outfile = output_prefix + "_end.bam"
    stats_file = output_prefix + ".end"

    with (
        pysam.AlignmentFile(bam_path) as bam,
        pysam.FastaFile(reference_fasta) as fasta,
        pysam.AlignmentFile(outfile, "wb", template=bam) as outbam,
        open(stats_file, "w") as stats,
    ):
        stats.write(
            "query_name\tref_name\tref_start\tref_end\tlocation\tcleavage_direction\t"
            "rbase\tqbase\tmutations\tflank_seq\tPoly\n"
        )

        for read in bam:
            # Apply read-level selection filters.
            if (
                read.is_unmapped
                or read.is_duplicate
                or "S" in read.cigarstring
                or "H" in read.cigarstring
                or read.mapping_quality < 30
            ):
                continue
            if read.is_read1:
                continue

            if read.is_forward:
                location = read.reference_start + 1
                try:
                    rbase = fasta.fetch(reference=read.reference_name, start=location, end=location + 1)
                    qbase = read.query[1:2]

                    if read.query_qualities is not None and len(read.query_qualities) > 1:
                        start_qual = read.query_qualities[0]
                        location_qual = read.query_qualities[1]
                        qual_check = start_qual > 20 and location_qual > location_qual_threshold
                    else:
                        qual_check = False

                    if rbase.upper() == "T" and qbase.upper() == "C" and qual_check:
                        write_selected_read(
                            stats, outbam, read, fasta, location, "-", rbase, qbase, "T"
                        )
                # Preserve the historical behavior: report a malformed read and
                # continue processing the remainder of the BAM.
                except Exception as exc:  # noqa: BLE001
                    print(read, file=sys.stderr)
                    print(f"error occurred at:{read.reference_name}:{location + 1}", file=sys.stderr)
                    print(exc, file=sys.stderr)
            else:
                location = read.reference_end
                try:
                    rbase = fasta.fetch(reference=read.reference_name, start=location - 2, end=location - 1)
                    qbase = read.query[-2:-1]

                    if read.query_qualities is not None and len(read.query_qualities) >= 2:
                        end_qual = read.query_qualities[-1]
                        location_qual = read.query_qualities[-2]
                        qual_check = end_qual > 20 and location_qual > location_qual_threshold
                    else:
                        qual_check = False

                    if rbase.upper() == "A" and qbase.upper() == "G" and qual_check:
                        write_selected_read(
                            stats, outbam, read, fasta, location, "+", rbase, qbase, "A"
                        )
                # Preserve the historical per-read error handling here as well.
                except Exception as exc:  # noqa: BLE001
                    print(read, file=sys.stderr)
                    print(f"error occurred at:{read.reference_name}:{location - 1}", file=sys.stderr)
                    print(exc, file=sys.stderr)

    pysam.index(outfile)


def parse_args():
    parser = argparse.ArgumentParser(
        description="Identify Ino-seq ABE signature reads and generate read-level statistics."
    )
    parser.add_argument("output_prefix")
    parser.add_argument("bam")
    parser.add_argument("reference_fasta")
    parser.add_argument("location_qual_threshold", nargs="?", type=int, default=30)
    return parser.parse_args()


def main():
    args = parse_args()
    run(args.output_prefix, args.bam, args.reference_fasta, args.location_qual_threshold)


if __name__ == "__main__":
    main()
