#!/usr/bin/env python3
"""
Grab phenotype discovery GWAS sumstats for the MVP traits listed in the cleaned
phenotype file. For each MVP trait, the script finds matching entries in the
table_of_contents files, extracts the corresponding files from the .tar
archives, and writes them to ancestry/ICD10_category specific directories. A
simple summary (per ancestry) and a log of missing phecodes are also produced.
"""

from __future__ import annotations

import glob
import logging
import os
import tarfile
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple

import pandas as pd


REPO_ROOT = Path(__file__).resolve().parents[1]
PHENO_PATH = REPO_ROOT / "data/4-integrate_MVP_UKB/2025_casecount/MVP_UKB_final_cleaned_pheno_df.csv"
MVP_GWAS_DIR = REPO_ROOT / "data/5-MVP_download"
OUTPUT_DIR = REPO_ROOT / "data/6-grab_phenotype_discovery_GWAS"
ANCESTRIES = {"AFR", "AMR", "EAS", "EUR", "META"}
LOGGER_NAME = "grab_phenotype_discovery_gwas"


def setup_logging(output_dir: Path) -> logging.Logger:
    """Initialise file + console logging."""
    output_dir.mkdir(parents=True, exist_ok=True)
    log_path = output_dir / "grab_phenotype_discovery_GWAS.log"
    logger = logging.getLogger(LOGGER_NAME)
    logger.setLevel(logging.INFO)

    if logger.handlers:
        return logger

    formatter = logging.Formatter("%(asctime)s [%(levelname)s] %(message)s", "%Y-%m-%d %H:%M:%S")

    fh = logging.FileHandler(log_path)
    fh.setLevel(logging.INFO)
    fh.setFormatter(formatter)
    logger.addHandler(fh)

    ch = logging.StreamHandler()
    ch.setLevel(logging.INFO)
    ch.setFormatter(formatter)
    logger.addHandler(ch)

    logger.info("Logging initialised. Log file: %s", log_path)
    return logger


def parse_table_of_contents(toc_path: Path, logger: logging.Logger) -> pd.DataFrame:
    """Parse a single table_of_contents.txt into a dataframe."""
    rows: List[Dict[str, object]] = []
    tar_path = Path(str(toc_path).replace(".table_of_contents.txt", ""))

    with toc_path.open("r") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line or line.startswith("d"):
                continue
            parts = line.split(maxsplit=5)
            if len(parts) < 6:
                logger.warning("Skipping unparsable line in %s: %r", toc_path, line)
                continue

            member_path = parts[5]
            filename = os.path.basename(member_path)
            tokens = filename.split(".")
            phecode = next((tok for tok in tokens if tok.startswith("Phe_")), None)
            ancestry = next((tok for tok in tokens if tok in ANCESTRIES), None)

            rows.append(
                {
                    "member_path": member_path,
                    "filename": filename,
                    "phecode": phecode,
                    "ancestry": ancestry,
                    "toc_path": str(toc_path),
                    "tar_path": str(tar_path),
                    "is_metadata": filename.endswith("metadata.txt"),
                }
            )

    if not rows:
        logger.warning("No entries parsed from %s", toc_path)
        return pd.DataFrame(
            columns=[
                "member_path",
                "filename",
                "phecode",
                "ancestry",
                "toc_path",
                "tar_path",
                "is_metadata",
            ]
        )

    return pd.DataFrame(rows)


def load_all_table_of_contents(mvp_dir: Path, logger: logging.Logger) -> pd.DataFrame:
    """Load and concatenate all table_of_contents files in the MVP directory."""
    toc_files = sorted(glob.glob(str(mvp_dir / "*.table_of_contents.txt")))
    logger.info("Found %d table_of_contents files in %s", len(toc_files), mvp_dir)

    frames: List[pd.DataFrame] = []
    for toc_file in toc_files:
        df = parse_table_of_contents(Path(toc_file), logger)
        if not df.empty:
            frames.append(df)

    if not frames:
        logger.warning("No table_of_contents data available; returning empty dataframe.")
        return pd.DataFrame()

    combined = pd.concat(frames, ignore_index=True)
    combined = combined.dropna(subset=["phecode"])
    logger.info("Aggregated %d entries across all table_of_contents files.", len(combined))
    return combined


def extract_member(
    tar_path: Path,
    member_path: str,
    dest_dir: Path,
    logger: logging.Logger,
) -> Optional[Path]:
    """Extract a single member from the tar into dest_dir; return destination path or None."""
    if not tar_path.exists():
        logger.warning("Tar file not found: %s", tar_path)
        return None

    dest_dir.mkdir(parents=True, exist_ok=True)
    dest_path = dest_dir / os.path.basename(member_path)

    if dest_path.exists():
        logger.info("Skipping existing file: %s", dest_path)
        return dest_path

    try:
        with tarfile.open(tar_path, "r:*") as tar:
            try:
                member = tar.getmember(member_path)
            except KeyError:
                logger.warning("Member %s not found in %s", member_path, tar_path)
                return None

            if not member.isfile():
                logger.warning("Member %s in %s is not a regular file; skipping.", member_path, tar_path)
                return None

            with tar.extractfile(member) as src, dest_path.open("wb") as dst:
                if src is None:
                    logger.warning("Could not extract %s from %s", member_path, tar_path)
                    return None
                for chunk in iter(lambda: src.read(1024 * 1024), b""):
                    dst.write(chunk)

        logger.info("Extracted %s to %s", member_path, dest_path)
        return dest_path
    except tarfile.ReadError as exc:
        logger.error("Failed to open tar %s: %s", tar_path, exc)
        return None


def build_trait_groups(pheno: pd.DataFrame) -> Dict[str, List[str]]:
    """Group MVP_Trait values by ICD10_category."""
    trait_df = pheno[["ICD10_category", "MVP_Trait"]].dropna()
    trait_df = trait_df[trait_df["MVP_Trait"].astype(str).str.startswith("Phe_")]
    trait_df = trait_df.drop_duplicates()

    grouped: Dict[str, List[str]] = {}
    for category, group in trait_df.groupby("ICD10_category"):
        grouped[category] = sorted(group["MVP_Trait"].unique())
    return grouped


def write_summary_matrices(
    ancestry_hits: Dict[str, Set[Tuple[str, str]]],
    categories: List[str],
    traits: List[str],
    output_dir: Path,
    logger: logging.Logger,
) -> None:
    """Create ancestry-specific CSVs indicating which phecodes were found per ICD10 category."""
    for ancestry, hits in ancestry_hits.items():
        matrix = pd.DataFrame(False, index=categories, columns=traits)
        for category, trait in hits:
            if category in matrix.index and trait in matrix.columns:
                matrix.loc[category, trait] = True

        summary_path = output_dir / f"{ancestry}_icd10_phecode_matrix.csv"
        matrix.to_csv(summary_path)
        logger.info("Saved ancestry summary: %s", summary_path)


def main() -> None:
    logger = setup_logging(OUTPUT_DIR)

    if not PHENO_PATH.exists():
        logger.error("Phenotype file not found: %s", PHENO_PATH)
        return

    pheno = pd.read_csv(PHENO_PATH)
    logger.info("Loaded phenotype table with %d rows from %s", len(pheno), PHENO_PATH)

    toc_entries = load_all_table_of_contents(MVP_GWAS_DIR, logger)
    if toc_entries.empty:
        logger.error("No table_of_contents entries to search; exiting.")
        return

    trait_groups = build_trait_groups(pheno)
    categories = sorted(trait_groups.keys())
    all_traits = sorted(pheno["MVP_Trait"].dropna().astype(str).unique())

    missing: List[Tuple[str, str]] = []
    ancestry_hits: Dict[str, Set[Tuple[str, str]]] = defaultdict(set)
    extracted_records: List[Dict[str, object]] = []

    for category, traits in trait_groups.items():
        logger.info("Processing ICD10 category %s with %d MVP traits.", category, len(traits))
        for trait in traits:
            trait_matches = toc_entries[toc_entries["phecode"] == trait]
            if trait_matches.empty:
                missing.append((category, trait))
                logger.warning("No files found for %s in category %s", trait, category)
                continue

            for _, row in trait_matches.iterrows():
                ancestry = row["ancestry"] or "UNKNOWN"
                dest_dir = OUTPUT_DIR / ancestry / category
                dest_path = extract_member(Path(row["tar_path"]), row["member_path"], dest_dir, logger)
                if dest_path:
                    ancestry_hits[ancestry].add((category, trait))
                    extracted_records.append(
                        {
                            "ICD10_category": category,
                            "MVP_Trait": trait,
                            "ancestry": ancestry,
                            "tar_path": row["tar_path"],
                            "member_path": row["member_path"],
                            "dest_path": str(dest_path),
                            "is_metadata": row["is_metadata"],
                        }
                    )

    if missing:
        missing_df = pd.DataFrame(missing, columns=["ICD10_category", "MVP_Trait"])
        missing_path = OUTPUT_DIR / "missing_phecodes.csv"
        missing_df.to_csv(missing_path, index=False)
        logger.info("Recorded %d missing phecodes to %s", len(missing), missing_path)
    else:
        logger.info("All phecodes were found in the table_of_contents files.")

    if extracted_records:
        extracted_df = pd.DataFrame(extracted_records)
        extracted_log_path = OUTPUT_DIR / "extracted_files.csv"
        extracted_df.to_csv(extracted_log_path, index=False)
        logger.info("Wrote extracted file log to %s", extracted_log_path)
    else:
        logger.warning("No files were extracted.")

    if ancestry_hits:
        write_summary_matrices(ancestry_hits, categories, all_traits, OUTPUT_DIR, logger)
    else:
        logger.warning("No ancestry hits to summarise.")

    logger.info(
        "Finished. Extracted %d files across %d categories. Missing %d phecodes.",
        len(extracted_records),
        len(categories),
        len(missing),
    )


if __name__ == "__main__":
    main()
