#!/usr/bin/env python3
"""Align an sgRNA sequence to candidate-site reference windows.

The fuzzy-matching and bulge-scoring behavior is retained from the historical
``identify.py`` implementation, including searches on both strands and the
weighted score ``substitutions + 3 * (insertions + deletions)``.
"""

from __future__ import annotations

import argparse
import logging
from pathlib import Path

import pyfaidx
import regex
from Bio.Data import IUPACData

from workflow.lib.io_utils import ensure_parent

LOGGER = logging.getLogger("inoseq.align_sgrna")

ALIGNMENT_COLUMNS = [
    "WindowSequence",
    "Site_SubstitutionsOnly.Sequence",
    "Site_SubstitutionsOnly.NumSubstitutions",
    "Site_SubstitutionsOnly.Strand",
    "Site_SubstitutionsOnly.Start",
    "Site_SubstitutionsOnly.End",
    "Site_GapsAllowed.Sequence",
    "Site_GapsAllowed.Length",
    "Site_GapsAllowed.Score",
    "Site_GapsAllowed.Substitutions",
    "Site_GapsAllowed.Insertions",
    "Site_GapsAllowed.Deletions",
    "Site_GapsAllowed.Strand",
    "Site_GapsAllowed.Start",
    "Site_GapsAllowed.End",
    "RealignedTargetSequence",
]


def reverse_complement(sequence: str) -> str:
    table = str.maketrans("ACGTacgt", "TGCATGCA")
    return sequence.translate(table)[::-1]


def regex_from_sequence(
    sequence: str, lookahead: bool = True, indels: int = 1, errors: int = 7
) -> tuple[str, str]:
    values = IUPACData.ambiguous_dna_values
    pattern = "".join(f"[{values[base]}]" for base in sequence.upper())
    if lookahead:
        pattern = f"(?b:{pattern})"
    standard = pattern + f"{{s<={errors}}}"
    gapped = pattern + (
        f"{{i<={indels},d<={indels},s<={errors},3i+3d+1s<={errors}}}"
    )
    return standard, gapped


def extended_pattern(sequence: str, indels: int = 1, errors: int = 7) -> str:
    extended = {
        "N": "[ATCGN]",
        "-": "[ATCGN]",
        "Y": "[CTY]",
        "R": "[AGR]",
        "W": "[ATW]",
        "S": "[CGS]",
        "A": "A",
        "T": "T",
        "C": "C",
        "G": "G",
        "V": "[ACG]",
        "M": "[AC]",
        "K": "[GT]",
        "H": "[ATC]",
        "B": "[GTC]",
        "D": "[GAT]",
    }
    pattern = "".join(extended[base] for base in sequence.upper())
    return f"(?b:{pattern})" + (
        f"{{i<={indels},d<={indels},s<={errors},3i+3d+1s<={errors}}}"
    )


def realigned_sequences(target: str, alignment: regex.Match, errors: int = 7) -> tuple[object, str]:
    match_sequence = alignment.group()
    substitutions, insertions, deletions = alignment.fuzzy_counts
    realigned_fuzzy = (substitutions, max(0, insertions - 1), max(0, deletions - 1))

    if insertions:
        if target.find("N") > len(target) / 2:
            target_realignments = [
                target[: index + 1] + "-" + target[index + 1 :]
                for index in range(target.find("N") + 1)
            ]
        else:
            target_realignments = [
                target[:index] + "-" + target[index:]
                for index in range(target.find("N"), len(target))
            ]
    else:
        target_realignments = [target]

    selected_target: object = None
    selected_offtarget = ""
    for candidate in target_realignments:
        if deletions:
            match_realignments = [
                match_sequence[: index + 1] + "-" + match_sequence[index + 1 :]
                for index in range(len(match_sequence) - 1)
            ]
            match_patterns = [
                match_sequence[: index + 1]
                + candidate[index + 1]
                + match_sequence[index + 1 :]
                for index in range(len(match_sequence) - 1)
            ]
        else:
            match_realignments = match_patterns = [match_sequence]

        # Compatibility freeze: the historical positional call passed
        # ``errors`` into the indel-limit argument while leaving errors=7.
        pattern = extended_pattern(candidate, indels=errors, errors=7)
        for match_pattern, match_alignment in zip(match_patterns, match_realignments):
            realignment = regex.search(pattern, match_pattern, regex.BESTMATCH)
            if realignment and realignment.fuzzy_counts == realigned_fuzzy:
                selected_target = candidate
                selected_offtarget = match_alignment
    return selected_target, selected_offtarget


def align_sequences(target: str, window: str, max_score: int = 8) -> list[object]:
    window = window.upper()
    standard_pattern, gapped_pattern = regex_from_sequence(target, errors=max_score)
    standard_alignments = [
        ("+", regex.search(standard_pattern, window, regex.BESTMATCH)),
        ("-", regex.search(standard_pattern, reverse_complement(window), regex.BESTMATCH)),
    ]
    gapped_alignments = [
        ("+", regex.search(gapped_pattern, window, regex.BESTMATCH)),
        ("-", regex.search(gapped_pattern, reverse_complement(window), regex.BESTMATCH)),
    ]

    chosen_standard = None
    chosen_standard_strand = ""
    lowest_mismatch = max_score + 1
    for strand, match in standard_alignments:
        if match is not None and match.fuzzy_counts[0] < lowest_mismatch:
            chosen_standard = match
            chosen_standard_strand = strand
            lowest_mismatch = match.fuzzy_counts[0]

    chosen_gapped = None
    chosen_gapped_strand = ""
    lowest_distance_score = 100
    for strand, match in gapped_alignments:
        if match is None:
            continue
        substitutions, insertions, deletions = match.fuzzy_counts
        if not (insertions or deletions):
            continue
        distance_score = substitutions + 3 * (insertions + deletions)
        edit_distance = substitutions + insertions + deletions
        if distance_score < lowest_distance_score and edit_distance < lowest_mismatch:
            chosen_gapped = match
            chosen_gapped_strand = strand
            lowest_distance_score = distance_score

    if chosen_standard:
        standard_values: list[object] = [
            chosen_standard.group(),
            chosen_standard.fuzzy_counts[0],
            chosen_standard_strand,
            chosen_standard.start(),
            chosen_standard.end(),
        ]
    else:
        standard_values = ["", "", "", "", ""]

    gapped_values: list[object] = ["", "", "", "", "", "", "", "", "", "none"]
    if chosen_gapped:
        realigned_target, realigned_offtarget = realigned_sequences(
            target, chosen_gapped, max_score
        )
        if realigned_offtarget:
            substitutions, insertions, deletions = chosen_gapped.fuzzy_counts
            gapped_values = [
                realigned_offtarget,
                len(chosen_gapped.group()),
                substitutions + 3 * (insertions + deletions),
                substitutions,
                insertions,
                deletions,
                chosen_gapped_strand,
                chosen_gapped.start(),
                chosen_gapped.end(),
                realigned_target,
            ]
        else:
            chosen_gapped_strand = ""
    return [*standard_values, *gapped_values]


def get_sequence(genome: pyfaidx.Fasta, chromosome: str, start: int, end: int) -> str:
    return str(genome[chromosome][start:end])


def _absolute_bounds(
    strand: object, start: object, end: object, window_start: int, window_end: int
) -> tuple[object, object]:
    if strand == "+":
        return int(start) + window_start, int(end) + window_start
    if strand == "-":
        return window_end - int(end), window_end - int(start)
    return "", ""


def process_candidates(
    input_file: str | Path,
    reference_fasta: str | Path,
    sgrna: str,
    output_file: str | Path,
    window_size: int = 25,
    max_score: int = 8,
) -> int:
    output_path = ensure_parent(output_file)
    genome = pyfaidx.Fasta(str(reference_fasta))
    processed = 0
    try:
        with open(input_file) as source, open(output_path, "w") as destination:
            header = source.readline().rstrip("\n")
            if not header:
                raise ValueError(f"empty candidate file: {input_file}")
            destination.write(header + "\t" + "\t".join(ALIGNMENT_COLUMNS) + "\n")

            for line_number, raw_line in enumerate(source, start=2):
                line = raw_line.rstrip("\n")
                if not line.strip():
                    continue
                fields = line.split("\t")
                if len(fields) < 8:
                    raise ValueError(f"{input_file}:{line_number}: expected at least 8 columns")
                chromosome = fields[0]
                start = int(fields[1])
                end = int(fields[2])
                max_position = int(fields[5]) if fields[5] else 0
                center = (start + end) // 2 if fields[6].lower() == "no" else max_position
                window_start = center - window_size
                window_end = center + window_size

                try:
                    window_sequence = get_sequence(
                        genome, chromosome, window_start, window_end
                    )
                    values = align_sequences(sgrna, window_sequence, max_score)
                    standard_start, standard_end = _absolute_bounds(
                        values[2], values[3], values[4], window_start, window_end
                    )
                    gapped_start, gapped_end = _absolute_bounds(
                        values[11], values[12], values[13], window_start, window_end
                    )
                    output_values = [
                        window_sequence,
                        values[0],
                        values[1],
                        values[2],
                        standard_start,
                        standard_end,
                        values[5],
                        values[6],
                        values[7],
                        values[8],
                        values[9],
                        values[10],
                        values[11],
                        gapped_start,
                        gapped_end,
                        values[14],
                    ]
                except Exception as exc:  # noqa: BLE001 - preserve historical blank-row behavior.
                    LOGGER.warning(
                        "Could not analyze %s:%s-%s: %s",
                        chromosome,
                        window_start,
                        window_end,
                        exc,
                    )
                    output_values = [""] * len(ALIGNMENT_COLUMNS)
                destination.write(line + "\t" + "\t".join(map(str, output_values)) + "\n")
                processed += 1
    finally:
        genome.close()
    return processed


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True)
    parser.add_argument("--ref", required=True, help="Reference-genome FASTA")
    parser.add_argument("--sgrna", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--window", type=int, default=25)
    parser.add_argument("--max-score", type=int, default=8)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    logging.basicConfig(level=logging.INFO, format="[%(levelname)s] %(message)s")
    rows = process_candidates(
        args.input, args.ref, args.sgrna, args.output, args.window, args.max_score
    )
    print(f"[05 align] {rows} candidate(s) -> {args.output}")


if __name__ == "__main__":
    main()
