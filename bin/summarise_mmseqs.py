#!/usr/bin/env python
# /// script
# requires-python = ">=3.13,<3.14"
# dependencies = [
#     "polars==1.35.1",
# ]
# ///
"""Provide a command line tool to generate stacked plot from taxonomies."""

import argparse
from pathlib import Path

import polars as pl
import polars.selectors as cs

TAXONOMIC_LEVELS = ["phylum", "class", "order", "family", "genus", "species"]
PIDENT_CUTOFF = 0.6
MINIMUM_TAXA_PERCENTAGE = 0.01


def main(argv=None):
    args = parse_args(argv)
    last_tax_level = "species"  # TODO: make this configurable
    tax_levels = TAXONOMIC_LEVELS[: TAXONOMIC_LEVELS.index(last_tax_level) + 1]

    ldf = load_mmseqs_file(args.mmseqs_tsv, args.min_identity, args.min_alignment_length, tax_levels)

    #
    # Case 1: Unique assignments
    #
    ldf_unique_matches = ldf.unique(subset="query", keep="none").select(
        ["query", "taxname", *TAXONOMIC_LEVELS[:-1]]
    )
    ldf_taxonomy = ldf.select(["taxname", *tax_levels]).unique(keep="first")

    # From reads with multiple alignments, keep only the first alignment
    uniq_ids = set(ldf_unique_matches.select("query").collect().to_series())
    ldf_first_match = (
        ldf.filter(~pl.col("query").is_in(uniq_ids))
        .sort(["pident", "tcov"], descending=True)
        .unique(subset=["query"], keep="first")
        .rename(lambda col: f"{col}_first" if col in ["target", "pident", "taxname", *tax_levels] else col)
    )

    # All reads having multiple alignments, compared against the first alignment
    ldf_multi_match = (
        ldf_first_match.join(ldf, on="query", how="left")
        .with_columns(pident_diff=(pl.col("pident_first") - pl.col("pident")).round(2))
        .filter(pl.col("pident_diff") < PIDENT_CUTOFF)
        .with_columns(
            is_unique=pl.len().over("query") == 1,
        )
    )

    #
    # Case 2: Secondary alignments with a difference with first alignment higher than PIDENT_CUTOFF
    #         In these cases keep the first alignment only
    #
    ldf_invalid_sec_matches = (
        ldf_multi_match.filter(pl.col("is_unique"))
        .select("query", cs.matches(r"*_first$"))
        .rename(lambda col: col.replace("_first", ""))
        .select("query", "taxname", *TAXONOMIC_LEVELS[:-1])
    )

    # Reads with valid secondary alignments
    ldf_valid_sec_matches = ldf_multi_match.filter(~pl.col("is_unique"))

    #
    # Case 3: Reads with secondary alignments with different genus than the first alignment
    #
    ldf_valid_sec_alns_diff_genus = ldf_valid_sec_matches.filter(
        (pl.col("genus") != pl.col("genus_first")).any().over("query")
    )

    ldf_mixed_genus = (
        ldf_valid_sec_alns_diff_genus.group_by("query")
        .agg(cs.by_name(*TAXONOMIC_LEVELS[:-1]).unique().sort().str.join(" / "))
        .with_columns(
            taxname=pl.col("genus") + pl.lit(" (mixed species)"),
        )
        # remove those reads that have alignments with more than two genus
        .filter(pl.col("genus").str.count_matches("/").lt(2))
        .select("query", "taxname", *TAXONOMIC_LEVELS[:-1])
    )

    # Update taxonomy with mixed genus entries
    ldf_taxonomy = pl.concat(
        [
            ldf_taxonomy.filter(~pl.col("taxname").str.contains("mixed")).select("taxname", *TAXONOMIC_LEVELS[:-1]),
            ldf_mixed_genus.select("taxname", *TAXONOMIC_LEVELS[:-1]),
        ]
    )

    #
    # Case 4: Reads with secondary alignments with same genus and same identity
    #         In these cases keep the first alignment only
    #
    ldf_valid_sec_matches_same_genus_zero = (
        ldf_valid_sec_matches.filter(
            (pl.col("genus") == pl.col("genus_first")).all().over("query"),
            (pl.col("species") != pl.col("species_first")) if last_tax_level == "species" else pl.lit(True),
        )
        .filter(pl.col("pident_diff") == 0)
        .unique(subset="query", keep="first")
        .select("query", "taxname", *TAXONOMIC_LEVELS[:-1])
        .with_columns(taxname=pl.col("genus") + pl.lit(" (mixed species)"))
    )

    #
    # Case 5: Reads with secondary alignments with same genus but different identity
    #
    same_genus_zero_ids = set(ldf_valid_sec_matches_same_genus_zero.select("query").collect().to_series())
    ldf_valid_sec_matches_same_genus_non_zero = (
        ldf_valid_sec_matches.filter(~pl.col("query").is_in(same_genus_zero_ids))
        .select("query", "taxname", *[f"{col}_first" for col in TAXONOMIC_LEVELS[:-1]])
        .rename(lambda col: col.replace("_first", ""))
    )

    ldf_final = pl.concat(
        [
            ldf_unique_matches,
            ldf_invalid_sec_matches,
            ldf_valid_sec_matches_same_genus_zero,
            ldf_valid_sec_matches_same_genus_non_zero,
            ldf_mixed_genus,
        ]
    ).rename(mapping={"taxname": "species"})

    #
    # Summarize abundances and write output files
    #

    # Collect once to avoid executing the full pipeline for each output
    df_final = ldf_final.collect()
    df_taxonomy = ldf_taxonomy.unique(keep="first").collect()

    for tax_level in tax_levels:
        df_abundance = (
            df_final.lazy()
            .group_by(tax_level)
            .len("count")
            .sort("count", descending=True)
            .with_columns(
                perc=(pl.col("count") / pl.col("count").sum() * 100).round(2),
            )
            # only keep taxa with at least 0.01% abundance and recalculate percentages
            .filter(pl.col("perc") >= MINIMUM_TAXA_PERCENTAGE)
            .with_columns(perc=(pl.col("count") / pl.col("count").sum() * 100).round(2))
        )

        # Add metadata columns dynamically
        if args.metadata:
            df_abundance = df_abundance.with_columns(
                **{key: pl.lit(value) for key, value in args.metadata.items()}
            )

        # Use sample name from metadata for output filename, or fallback to "output"
        sample_name = args.metadata.get("sample", "output") if args.metadata else "output"
        df_abundance.sink_csv(f"{sample_name}_{tax_level}.csv")

    # Use sample name from metadata for output filename, or fallback to "output"
    sample_name = args.metadata.get("sample", "output") if args.metadata else "output"

    Path("taxlineage").mkdir(exist_ok=True)
    df_taxonomy.write_csv(f"taxlineage/{sample_name}_taxlineage.csv")

    # Picrust input
    (
        df_final.select("query")
        .with_columns(pl.lit(1).alias(sample_name))
        .write_csv(f"reads_{sample_name}.tsv", separator="\t", include_header=True)
    )


def load_mmseqs_file(
    file_path: str, min_identity: float, min_alignment_length: int, tax_levels: list[str] = TAXONOMIC_LEVELS
) -> pl.LazyFrame:
    """Load MMseqs file and filter by identity and alignment length."""
    colnames = ["query", "target", "pident", "tcov", "alnlen", "taxname", "taxlineage"]

    ldf = (
        pl.scan_csv(file_path, separator="\t", has_header=False, new_columns=colnames)
        .filter((pl.col("pident") >= (min_identity) * 100) & (pl.col("alnlen") >= min_alignment_length))
        .with_columns(
            pident=pl.col("pident").round(1),
            tcov=pl.col("tcov").round(3),
        )
        # extract taxonomic levels from taxlineage
        .with_columns(
            *[
                pl.col("taxlineage").str.extract(rf".*;{level[0]}_([^;]*)(?:;.*)?", 1).alias(level)
                for level in tax_levels
            ],
        )
        # fill null values with parent taxonomic levels
        .with_columns(
            *[
                pl.col(level).fill_null(pl.lit(f"{plevel}:") + pl.col(plevel))
                for level, plevel in zip(tax_levels[1:], tax_levels)
            ]
        )
        .sort(["pident", "tcov"], descending=True)
        # if a read has more than one alignment to the same taxon, keep the first only
        .unique(subset=["query", "taxname"], keep="first")
    )
    return ldf


def parse_metadata(metadata_str: str) -> dict[str, str]:
    """Parse metadata string in format key=value,key2=value2 into a dictionary."""
    if not metadata_str:
        return {}

    metadata = {}
    for item in metadata_str.split(","):
        if "=" not in item:
            raise argparse.ArgumentTypeError(f"Invalid metadata format: {item}. Expected key=value")
        key, value = item.split("=", 1)
        metadata[key.strip()] = value.strip()
    return metadata


def parse_args(argv=None):
    parser = argparse.ArgumentParser(description="Merge blast out file with lineages and summary")

    parser.add_argument(
        "-t",
        "--mmseqs-tsv",
        required=True,
        metavar="FILE",
        help="MMseqs2 output file in TSV format containing taxonomic assignments",
    )
    parser.add_argument(
        "-i",
        "--min-identity",
        required=False,
        default=0.95,
        metavar="FLOAT",
        type=float,
        help="Minimum sequence identity threshold (0-1) for filtering alignments",
    )
    parser.add_argument(
        "-a",
        "--min-alignment-length",
        required=False,
        default=1000,
        metavar="INT",
        type=int,
        help="Minimum alignment length in base pairs for filtering alignments",
    )
    parser.add_argument(
        "-m",
        "--metadata",
        required=False,
        default=None,
        metavar="KEY=VALUE",
        type=parse_metadata,
        help=(
            "Metadata columns to add to output. Format: key=value,key2=value2 (e.g., sample=S1,group=control)"
        ),
    )

    return parser.parse_args(argv)


if __name__ == "__main__":
    main()
