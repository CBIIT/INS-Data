"""
gather_dbgap_subset_clones.py
2026-03-13 ZD

This script creates cloned dataset entries for dbGaP studies that are
also available through other NCI data repositories (e.g. CRDC-GC, GDC, CTDC).

For each dbGaP study tagged with a storage distribution, this script:
    1. Clones the row from the final curated clean dbGaP output
    2. Replaces `dataset_source_repo` with the subset repository key
    3. Generates a new deterministic UUID5 based on the new repo + phs accession
    4. Blanks out `dataset_storage_distribution` (avoids redundant self-reference)

All other field values (including dataset_source_id and dataset_source_url,
which still point to dbGaP) are preserved unchanged.

This script is intended to be run AFTER package_output_data.py has produced
the final curated clean dbGaP TSV. It does not modify the source file.

Input:
    - config.DBGAP_OUTPUT_CURATED_CLEANED  (final packaged curated dbGaP TSV)
    - Subset CSVs in config.DBGAP_SUBSET_INPUT_DIR

Output:
    - config.DBGAP_SUBSET_CLONES_OUTPUT_PATH  (TSV of cloned subset datasets)
"""

import os
import sys
import uuid

import pandas as pd

# Append the project's root directory to the Python path
# This allows for importing config when running as part of main.py or alone
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import config


# Deterministic UUID5 namespace (same as used throughout the dbGaP pipeline)
UUID_NAMESPACE = uuid.UUID('12345678-1234-5678-1234-567812345678')


def load_subset_phs(subset_key: str) -> set:
    """Load a subset CSV and return the set of short phs accessions.

    Args:
        subset_key (str): Key from DBGAP_STORAGE_DISTRIBUTION_MAP
            (e.g. 'CRDC-GC', 'GDC', 'CTDC'). Used to build the filename.

    Returns:
        set: Short phs accession strings (e.g. {'phs002790', 'phs001287'})
    """

    filepath = os.path.join(
        config.DBGAP_SUBSET_INPUT_DIR,
        f"{subset_key}_subset_{config.DBGAP_CSV_VERSION}.csv"
    )

    if not os.path.exists(filepath):
        print(f"  WARNING: Subset file not found, skipping: {filepath}")
        return set()

    subset_df = pd.read_csv(filepath, keep_default_na=False)

    if 'accession' not in subset_df.columns:
        print(f"  WARNING: No 'accession' column in {filepath}, skipping.")
        return set()

    # Strip versioning to get short phs (e.g. phs002790.v9.p3 -> phs002790)
    short_phs = set(
        subset_df['accession'].astype(str).str.split('.').str[0]
    )

    return short_phs


def generate_clone_uuid(repo_key: str, source_id: str) -> str:
    """Generate a deterministic UUID5 for a cloned dataset entry.

    Uses the same namespace as the main dbGaP pipeline but with the
    subset repository key instead of 'dbGaP', producing a unique
    but reproducible UUID for each repo+phs combination.

    Args:
        repo_key (str): subset repository key (e.g. 'CRDC-GC')
        source_id (str): The dataset_source_id (phs accession)

    Returns:
        str: UUID5 string
    """

    return str(uuid.uuid5(
        UUID_NAMESPACE,
        '||'.join([repo_key, source_id])
    ))


def build_subset_clones(dbgap_df: pd.DataFrame) -> pd.DataFrame:
    """Build cloned dataset rows for all configured subset repositories.

    For each repository key in DBGAP_STORAGE_DISTRIBUTION_MAP:
        1. Load the corresponding subset CSV to get relevant phs accessions
        2. Filter the dbGaP dataframe to matching rows
        3. Clone those rows with updated repo, UUID, and blanked distribution

    Args:
        dbgap_df (pd.DataFrame): The final curated clean dbGaP dataset TSV.
            Must contain 'dataset_source_id', 'dataset_source_repo',
            'dataset_uuid', and 'dataset_storage_distribution' columns.

    Returns:
        pd.DataFrame: Concatenated cloned rows for all subset repositories.
    """

    all_clones = []

    for repo_key, label in config.DBGAP_STORAGE_DISTRIBUTION_MAP.items():

        print(f"\n  Processing: {repo_key} ({label})")

        # Get phs accessions for this subset
        subset_phs = load_subset_phs(repo_key)
        if not subset_phs:
            continue

        # Filter dbGaP rows that match this subset
        mask = dbgap_df['dataset_source_id'].isin(subset_phs)
        matched_df = dbgap_df[mask].copy()

        if len(matched_df) == 0:
            print(f"    No matching rows found in curated clean TSV.")
            continue

        # Clone and relabel
        matched_df['dataset_source_repo'] = repo_key
        matched_df['dataset_uuid'] = matched_df['dataset_source_id'].apply(
            lambda phs: generate_clone_uuid(repo_key, phs)
        )
        matched_df['dataset_storage_distribution'] = ''

        all_clones.append(matched_df)

        print(f"    Cloned {len(matched_df)} datasets as '{repo_key}' entries.")

    # Combine all clones into a single dataframe
    if all_clones:
        clones_df = pd.concat(all_clones, ignore_index=True)
    else:
        # Return empty df with same columns if no clones produced
        clones_df = pd.DataFrame(columns=dbgap_df.columns)

    return clones_df


def gather_dbgap_subset_clones():
    """Main function: read curated clean dbGaP TSV, clone tagged rows
    per subset repository, and export as a single TSV.

    Reads:
        config.DBGAP_OUTPUT_CURATED_CLEANED

    Writes:
        config.DBGAP_SUBSET_CLONES_OUTPUT_PATH
    """

    print(f"\n---\nDATASETS - dbGaP Subset Clones:\n"
          f"Cloning tagged dbGaP datasets as subset repository entries...\n---")

    # Read the final curated clean dbGaP output
    input_path = config.DBGAP_OUTPUT_CURATED_CLEANED
    if not os.path.exists(input_path):
        raise FileNotFoundError(
            f"Curated clean dbGaP file not found: {input_path}\n"
            f"Run package_output_data.py first to generate this file."
        )

    dbgap_df = pd.read_csv(input_path, sep='\t', dtype=str,
                           keep_default_na=False)
    print(f"\nLoaded {len(dbgap_df)} dbGaP datasets from {input_path}")

    # Validate required columns
    required_cols = {'dataset_source_id', 'dataset_source_repo',
                     'dataset_uuid', 'dataset_storage_distribution'}
    missing = required_cols - set(dbgap_df.columns)
    if missing:
        raise KeyError(f"Missing required columns: {missing}")

    # Build clones for all subset repositories
    clones_df = build_subset_clones(dbgap_df)

    # Summary
    print(f"\n{'='*60}")
    print(f"  dbGaP Subset Clone Summary")
    print(f"{'='*60}")
    print(f"  Source dbGaP rows:     {len(dbgap_df):>6}")
    print(f"  Total clones created:  {len(clones_df):>6}")
    print(f"{'─'*60}")
    if len(clones_df) > 0:
        for repo, count in clones_df['dataset_source_repo'].value_counts().items():
            print(f"    {repo:<12}  {count:>6} clones")
    print(f"{'='*60}")

    # Export
    output_path = config.DBGAP_SUBSET_CLONES_OUTPUT_PATH
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    clones_df.to_csv(output_path, sep='\t', index=False, encoding='utf-8')

    print(f"\nSuccess! Subset clones saved to {output_path}\n")

    return clones_df


# Run module as a standalone script when called directly
if __name__ == "__main__":

    print(f"Running {os.path.basename(__file__)} as standalone module...")

    gather_dbgap_subset_clones()
