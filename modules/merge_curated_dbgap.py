"""
merge_curated_dbgap.py
2026-03-11 ZD

This script merges a previously curated dbGaP datasets TSV with a newly
gathered dbGaP datasets TSV. The merge uses `dataset_source_id` (phs
accession) as the key:

    - Rows present in the OLD curated file are kept as-is (preserving
      hand-edited values and existing UUIDs).
    - Rows present only in the NEW file are appended (these are newly
      added studies that have not yet been curated).
    - Rows present only in the OLD file are retained with a warning
      (they may have been removed from the latest dbGaP search results
      but should not be silently dropped from the curated set).

After merging, the script detects title changes between old and new for
shared studies and exports a review CSV so the user can decide whether
to accept the new title or keep the old one.  If a reviewed CSV is found
from a prior run, the user's decisions are applied automatically.

The merged output is saved as `dbgap_datasets_merged.tsv` in the current
dbGaP output directory defined in config.py.
"""

import os
import sys
import uuid

import pandas as pd

# Append the project's root directory to the Python path
# This allows for importing config when running as part of main.py or alone
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import config


# Key column used to identify unique studies across old and new datasets
MERGE_KEY = 'dataset_source_id'

# Column to review for changes between old and new
TITLE_COL = 'dataset_title'


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

    if MERGE_KEY not in df.columns:
        raise KeyError(f"Expected column '{MERGE_KEY}' not found in "
                       f"{filepath}. Columns found: {list(df.columns)}")

    return df


def detect_title_changes(old_df: pd.DataFrame,
                         new_df: pd.DataFrame) -> pd.DataFrame:
    """Identify studies where the title differs between old and new datasets.

    Args:
        old_df (pd.DataFrame): Old curated dataset.
        new_df (pd.DataFrame): New gathered dataset.

    Returns:
        pd.DataFrame: Review table with columns:
            dataset_source_id, old_title, new_title, use_new_title
    """

    # Build lookup dicts keyed on the merge key
    old_titles = old_df.set_index(MERGE_KEY)[TITLE_COL].to_dict()
    new_titles = new_df.set_index(MERGE_KEY)[TITLE_COL].to_dict()

    # Find shared IDs with differing titles
    changes = []
    for phs in sorted(set(old_titles) & set(new_titles)):
        old_val = str(old_titles[phs]).strip()
        new_val = str(new_titles[phs]).strip()
        if old_val != new_val:
            changes.append({
                MERGE_KEY: phs,
                'old_title': old_val,
                'new_title': new_val,
                'use_new_title': ''   # blank = user has not reviewed yet
            })

    review_df = pd.DataFrame(changes)
    return review_df


def export_review_csv(review_df: pd.DataFrame, review_path: str) -> None:
    """Export the title review table as a CSV for manual editing.

    The user should fill in the `use_new_title` column:
        - 'yes' or 'y' → accept the new title
        - 'no', 'n', or blank → keep the old (curated) title

    Args:
        review_df (pd.DataFrame): Output of detect_title_changes.
        review_path (str): Filepath for the review CSV.
    """

    os.makedirs(os.path.dirname(review_path), exist_ok=True)
    review_df.to_csv(review_path, index=False, encoding='utf-8')
    print(f"Title review CSV exported to {review_path}")
    print(f"  → {len(review_df)} title change(s) detected.")
    print(f"  → Edit the 'use_new_title' column (yes/no) and re-run to apply.\n")


def apply_title_decisions(merged_df: pd.DataFrame,
                          review_path: str) -> pd.DataFrame:
    """Apply user decisions from a reviewed title-change CSV.

    For rows where use_new_title is 'yes' or 'y', the title in the merged 
    dataframe is replaced with the new title from the review CSV.  All other
    rows keep their current (old curated) title.

    Args:
        merged_df (pd.DataFrame): The merged dataframe (old titles by default).
        review_path (str): Path to the reviewed CSV with user decisions.

    Returns:
        pd.DataFrame: Updated merged dataframe with accepted title changes.
    """

    review_df = pd.read_csv(review_path, dtype=str, keep_default_na=False)

    # Validate expected columns
    required_cols = {MERGE_KEY, 'old_title', 'new_title', 'use_new_title'}
    if not required_cols.issubset(set(review_df.columns)):
        missing = required_cols - set(review_df.columns)
        raise KeyError(f"Review CSV is missing required columns: {missing}")

    # Filter to only rows where the user accepted the new title
    accepted = review_df[
        review_df['use_new_title'].str.strip().str.lower().isin(['yes', 'y'])
    ]

    if len(accepted) == 0:
        print(f"No title changes accepted in review CSV.\n")
        return merged_df

    # Apply accepted new titles
    accept_map = dict(zip(accepted[MERGE_KEY], accepted['new_title']))
    update_count = 0
    for idx, row in merged_df.iterrows():
        phs = row[MERGE_KEY]
        if phs in accept_map:
            merged_df.at[idx, TITLE_COL] = accept_map[phs]
            update_count += 1

    print(f"Applied {update_count} accepted title change(s) from review CSV.\n")

    return merged_df


def merge_dbgap_datasets(old_curated_path: str,
                         new_gathered_path: str,
                         output_path: str,
                         review_path: str) -> pd.DataFrame:
    """Merge an old curated dbGaP TSV with a new gathered dbGaP TSV.

    Precedence rules:
        1. Rows in the old curated file are always kept as-is.
        2. New rows (in new but not in old) are appended at the bottom.
        3. Old-only rows (in old but not in new) are retained with a 
           console warning.

    After merging, title changes are detected and either:
        - exported as a review CSV (first run), or
        - applied from the reviewed CSV (subsequent run).

    Args:
        old_curated_path (str): Filepath of the previously curated TSV.
        new_gathered_path (str): Filepath of the newly gathered TSV.
        output_path (str): Filepath for the merged output TSV.
        review_path (str): Filepath for the title review CSV.

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

    # --- Build the merged dataframe ---

    # 1. Start with ALL rows from the old curated file (shared + old-only)
    merged_df = old_df.copy()
    merged_df['curation_status'] = 'previously_curated'

    # 2. Append rows that exist only in the new file (at the end)
    new_only_df = new_df[new_df[MERGE_KEY].isin(new_only_ids)].copy()
    new_only_df['curation_status'] = 'new_study'
    merged_df = pd.concat([merged_df, new_only_df], ignore_index=True)

    # --- Console summary ---
    print(f"\n{'='*60}")
    print(f"  dbGaP Merge Summary")
    print(f"{'='*60}")
    print(f"  Old curated rows:       {len(old_df):>6}")
    print(f"  New gathered rows:      {len(new_df):>6}")
    print(f"{'─'*60}")
    print(f"  Shared (kept from old): {len(shared_ids):>6}")
    print(f"  New-only (added):       {len(new_only_ids):>6}")
    print(f"  Old-only (retained):    {len(old_only_ids):>6}")
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

    # --- Title change review ---
    print(f"Checking for title changes between old and new datasets...")
    review_df = detect_title_changes(old_df, new_df)

    if len(review_df) == 0:
        print(f"No title changes detected.\n")
    elif os.path.exists(review_path):
        # A reviewed CSV already exists — apply the user's decisions
        print(f"Found existing review CSV: {review_path}")
        merged_df = apply_title_decisions(merged_df, review_path)
    else:
        # First run — export the review CSV for user editing
        export_review_csv(review_df, review_path)

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
    review_path = config.DBGAP_MERGE_REVIEW_PATH

    merge_dbgap_datasets(old_curated_path, new_gathered_path,
                         output_path, review_path)
