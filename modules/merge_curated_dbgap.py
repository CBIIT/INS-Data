"""
merge_curated_dbgap.py
2026-03-11 ZD

This script merges a previously curated dbGaP datasets TSV with a newly
gathered dbGaP datasets TSV. The short `dataset_source_id` identifies a study,
while `dataset_source_accession` includes its dbGaP version and participant
set (for example, phs002790.v7.p1).

        - Unchanged accessions keep the previously curated row.
        - Updated accessions use the newly gathered row for curator review.
        - New studies use the newly gathered row for curator review.
        - Old-only studies keep the previously curated row for curator review.

The merged output is saved as `dbgap_datasets_merged.tsv` in the current
dbGaP output directory defined in config.py.
"""

import os
import uuid

import pandas as pd

import config


# Key column used to identify unique studies across old and new datasets
MERGE_KEY = 'dataset_source_id'

# Full dbGaP accession used to detect source release changes
ACCESSION_COL = 'dataset_source_accession'

STATUS_ORDER = ['updated', 'new', 'old_only', 'unchanged']


def load_dbgap_tsv(filepath: str) -> pd.DataFrame:
    """Load a dbGaP datasets TSV and perform basic validation.

    Args:
        filepath (str): Path to a tab-separated dbGaP datasets file.

    Returns:
        pd.DataFrame: Loaded dataframe.

    Raises:
        FileNotFoundError: If the filepath does not exist.
        KeyError: If the expected merge key column is missing.
    """

    if not os.path.exists(filepath):
        raise FileNotFoundError(f"File not found: {filepath}")

    df = pd.read_csv(filepath, sep='\t', dtype=str, keep_default_na=False)

    required_cols = {MERGE_KEY, ACCESSION_COL}
    missing_cols = required_cols - set(df.columns)
    if missing_cols:
        raise KeyError(f"Expected columns missing from {filepath}: "
                       f"{sorted(missing_cols)}. "
                       f"Columns found: {list(df.columns)}")

    duplicate_ids = df.loc[df[MERGE_KEY].duplicated(keep=False), MERGE_KEY]
    if not duplicate_ids.empty:
        raise ValueError(f"Duplicate '{MERGE_KEY}' values found in "
                         f"{filepath}: {sorted(duplicate_ids.unique())}")

    return df


def merge_dbgap_datasets(old_curated_path: str,
                         new_gathered_path: str,
                         output_path: str) -> pd.DataFrame:
    """Merge an old curated dbGaP TSV with a new gathered dbGaP TSV.

    Row selection rules:
        1. Shared studies with unchanged full accessions keep old curation.
        2. Shared studies with changed full accessions use new gathered data.
        3. New studies use new gathered data.
        4. Old-only studies retain old curated data.

    Args:
        old_curated_path (str): Filepath of the previously curated TSV.
        new_gathered_path (str): Filepath of the newly gathered TSV.
        output_path (str): Filepath for the merged output TSV.

    Returns:
        pd.DataFrame: The merged dataframe.
    """

    # Load both datasets
    print(f"Loading old curated dbGaP file: {old_curated_path}")
    old_df = load_dbgap_tsv(old_curated_path)

    print(f"Loading new gathered dbGaP file: {new_gathered_path}")
    new_df = load_dbgap_tsv(new_gathered_path)

    # Regenerate deterministic UUID5 values on the old curated data
    # This ensures UUIDs are consistent with the current uuid5 scheme
    # without modifying the original old curated TSV on disk
    NAMESPACE = uuid.UUID('12345678-1234-5678-1234-567812345678')
    old_df['dataset_uuid'] = old_df.apply(
        lambda row: str(uuid.uuid5(
            NAMESPACE,
            '||'.join([str(row['dataset_source_repo']),
                       str(row[MERGE_KEY])]))),
        axis=1)
    print(f"Regenerated deterministic UUID5 values for "
          f"{len(old_df)} old curated rows.")

    # Validate that both files share the same column schema
    if list(old_df.columns) != list(new_df.columns):
        old_set = set(old_df.columns)
        new_set = set(new_df.columns)
        only_old_cols = old_set - new_set
        only_new_cols = new_set - old_set
        print(f"\nWARNING: Column mismatch detected between old and new files.")
        if only_old_cols:
            print(f"  Columns only in old: {only_old_cols}")
        if only_new_cols:
            print(f"  Columns only in new: {only_new_cols}")
        print(f"  Proceeding with union of all columns. "
              f"Missing values will be blank.\n")

    # Identify the three groups using set operations on the merge key
    old_ids = set(old_df[MERGE_KEY])
    new_ids = set(new_df[MERGE_KEY])

    shared_ids = old_ids & new_ids
    old_only_ids = old_ids - new_ids
    new_only_ids = new_ids - old_ids

    old_accessions = old_df.set_index(MERGE_KEY)[ACCESSION_COL].to_dict()
    new_accessions = new_df.set_index(MERGE_KEY)[ACCESSION_COL].to_dict()
    updated_ids = {
        phs for phs in shared_ids
        if old_accessions[phs].strip() != new_accessions[phs].strip()
    }
    unchanged_ids = shared_ids - updated_ids

    # --- Build the merged dataframe ---

    groups = []
    for status, source_df, ids in [
            ('updated', new_df, updated_ids),
            ('new', new_df, new_only_ids),
            ('old_only', old_df, old_only_ids),
            ('unchanged', old_df, unchanged_ids)]:
        group_df = source_df[source_df[MERGE_KEY].isin(ids)].copy()
        group_df['curation_status'] = status
        group_df['previous_source_accession'] = group_df[MERGE_KEY].map(
            old_accessions).fillna('')
        group_df['current_source_accession'] = group_df[MERGE_KEY].map(
            new_accessions).fillna('')
        groups.append(group_df)

    merged_df = pd.concat(groups, ignore_index=True)
    merged_df['curation_status'] = pd.Categorical(
        merged_df['curation_status'], categories=STATUS_ORDER, ordered=True)
    merged_df.sort_values(
        ['curation_status', MERGE_KEY], inplace=True, ignore_index=True)
    merged_df['curation_status'] = merged_df['curation_status'].astype(str)

    # --- Console summary ---
    print(f"\n{'='*60}")
    print(f"  dbGaP Merge Summary")
    print(f"{'='*60}")
    print(f"  Old curated rows:       {len(old_df):>6}")
    print(f"  New gathered rows:      {len(new_df):>6}")
    print(f"{'─'*60}")
    print(f"  Updated (use new):      {len(updated_ids):>6}")
    print(f"  New:                    {len(new_only_ids):>6}")
    print(f"  Old-only:               {len(old_only_ids):>6}")
    print(f"  Unchanged (keep old):   {len(unchanged_ids):>6}")
    print(f"{'─'*60}")
    print(f"  MERGED TOTAL:           {len(merged_df):>6}")
    print(f"{'='*60}\n")

    # Warn about old-only rows
    if old_only_ids:
        print(f"NOTE: {len(old_only_ids)} study(ies) found only in the old "
              f"curated file (not in new search results). These have been "
              f"retained in the merged output:")
        for phs in sorted(old_only_ids):
            print(f"  - {phs}")
        print()

    # --- Export ---
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    merged_df.to_csv(output_path, sep='\t', index=False, encoding='utf-8')
    print(f"Merged dbGaP file saved to {output_path}\n")

    return merged_df


# Run module as a standalone script when called directly
if __name__ == "__main__":

    print(f"\nRunning {os.path.basename(__file__)} as standalone module...\n")

    # Default paths from config
    old_curated_path = config.DBGAP_MERGED_OLD_CURATED_PATH
    new_gathered_path = config.DBGAP_OUTPUT_PATH
    output_path = config.DBGAP_MERGED_OUTPUT_PATH

    merge_dbgap_datasets(old_curated_path, new_gathered_path, output_path)
