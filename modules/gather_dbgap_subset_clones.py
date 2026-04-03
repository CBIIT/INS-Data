"""
gather_dbgap_subset_clones.py
2026-03-13 ZD

This script creates cloned dataset entries for dbGaP studies that are
also available through other NCI data repositories (e.g. CRDC, GDC).

For each dbGaP study tagged with a storage distribution, this script:
    1. Clones the row from the final curated clean dbGaP output
    2. Replaces `dataset_source_repo` with the target repository name
       (looked up via config.DBGAP_SUBSET_REPO_MAP)
    3. Generates a new deterministic UUID5 based on the target repo + phs
    4. Preserves `dataset_storage_distribution` on the clone

Multiple subset keys can map to the same target repository. For example,
both the CRDC and CTDC subset CSVs produce clones labelled "CRDC".
When the same phs accession appears in more than one subset that maps to
the same repo, only one clone is produced (deduplicated).

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
            (e.g. 'CRDC', 'GDC', 'CTDC'). Used to build the filename.

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
    target repository name instead of 'dbGaP', producing a unique
    but reproducible UUID for each repo+phs combination.

    Args:
        repo_key (str): Target repository name (e.g. 'CRDC', 'GDC')
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
        3. Clone those rows with updated repo, UUID, and preserved distribution

    Multiple subset keys may map to the same target repo (via
    DBGAP_SUBSET_REPO_MAP). When this happens the phs sets are merged
    and only one clone per phs is produced for that repo.

    Args:
        dbgap_df (pd.DataFrame): The final curated clean dbGaP dataset TSV.
            Must contain 'dataset_source_id', 'dataset_source_repo',
            'dataset_uuid', and 'dataset_storage_distribution' columns.

    Returns:
        pd.DataFrame: Concatenated cloned rows for all subset repositories.
    """

    repo_map = config.DBGAP_SUBSET_REPO_MAP

    # Gather phs sets per target repo (merge keys that share the same repo)
    repo_phs: dict[str, set] = {}
    for subset_key in config.DBGAP_STORAGE_DISTRIBUTION_MAP:
        target_repo = repo_map[subset_key]

        print(f"\n  Loading subset: {subset_key}  →  target repo: {target_repo}")

        phs = load_subset_phs(subset_key)
        if not phs:
            continue

        print(f"    {len(phs)} accessions loaded from {subset_key} CSV")

        repo_phs.setdefault(target_repo, set()).update(phs)

    # Build clones per target repo (deduplicated)
    all_clones = []

    for target_repo, phs_set in repo_phs.items():

        print(f"\n  Cloning for target repo '{target_repo}' "
              f"({len(phs_set)} unique accessions)")

        mask = dbgap_df['dataset_source_id'].isin(phs_set)
        matched_df = dbgap_df[mask].copy()

        if len(matched_df) == 0:
            print(f"    No matching rows found in curated clean TSV.")
            continue

        # Clone and relabel — keep dataset_storage_distribution as-is
        matched_df['dataset_source_repo'] = target_repo
        matched_df['dataset_uuid'] = matched_df['dataset_source_id'].apply(
            lambda phs: generate_clone_uuid(target_repo, phs)
        )

        all_clones.append(matched_df)

        print(f"    Cloned {len(matched_df)} datasets as '{target_repo}' entries.")

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
