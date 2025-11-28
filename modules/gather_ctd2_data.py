"""
gather_ctd2_data.py
2025-11-26 ZD

This module defines primary function gather_ctd2_data that will process a CSV
of curated datasets from the Cancer Target Discovery and Development (CTD^2) 
Network. 
The intermediate ctd2_datasets.csv will be finalized with 
`package_output_data.py` into a TSV ready for loading within INS.

Note: The "gather" filename is used for consistency with other modules that handle
input data even though this module isn't pulling from automated external sources.

"""

import os
import sys
import uuid

import pandas as pd


# Append the project's root directory to the Python path
# This allows for importing config when running as part of main.py or alone
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import config

def get_composite_uuid5(df: pd.DataFrame, 
                        fields: list, 
                        uuid_col: str = 'dataset_uuid') -> pd.DataFrame:
    """
    Generate a UUID5 for each row in the DataFrame based on a composite of 
    specified fields. The UUID5 is deterministic and reproducible for the 
    same field values. Raises an error if duplicate UUIDs are generated.

    Args:
        df (pd.DataFrame): Input DataFrame.
        fields (list): List of column names to combine for the UUID5 'name'.
        uuid_col (str): Name of the output column for the UUIDs.

    Returns:
        pd.DataFrame: DataFrame with a new column containing the UUID5 values.
    """
    
    # Use a fixed namespace UUID for reproducibility (can be any valid UUID)
    NAMESPACE = uuid.UUID('12345678-1234-5678-1234-567812345678')

    def row_to_uuid(row):
        name = '||'.join(str(row[field]) for field in fields)
        return uuid.uuid5(NAMESPACE, name)

    df[uuid_col] = df.apply(row_to_uuid, axis=1)

    # Check for duplicate UUIDs and raise error if found
    duplicate_uuids = df[uuid_col][df[uuid_col].duplicated()]
    if not duplicate_uuids.empty:
        raise ValueError(
            f"Duplicate UUID5 values found in {uuid_col} column! "
            f"Duplicates:\n{duplicate_uuids.to_list()}\n"
            f"Check your input data for non-unique combinations of {fields}."
        )

    return df



def gather_ctd2_data():
    """Process a CSV of CTD^2 dataset metadata into an intermediate CSV for INS.

    Input and output filepaths are handled in config.py
    """

    print(f"\n---\nDATASETS - CTD^2:\n"
          f"Processing CTD^2 dataset data...\n---\n")

    # Define paths from config
    input_csv = config.CTD2_INPUT_CSV
    output_csv = config.CTD2_INTERMED_CSV

    # Load curated CTD^2 CSV
    df = pd.read_csv(input_csv)
    print(f"Loaded CTD^2 dataset data from {input_csv}")

    # Add blank required fields
    df['dataset_source_url'] = ''

    # Add reproducible uuid5
    uuid_fields = ['dataset_source_repo', 'dataset_title', 'description']
    df = get_composite_uuid5(df, uuid_fields, uuid_col='dataset_uuid')

    # Make output directories if they don't exist
    os.makedirs(os.path.dirname(output_csv), exist_ok=True)

    # Save output
    df.to_csv(output_csv, index=False)

    print(f"\n\nSuccess! CTD^2 datasets saved to {output_csv}.\n"
          f"Total CTD^2 dataset records:    {len(df)}")


# Run module as a standalone script when called directly
if __name__ == "__main__":

    print(f"Running {os.path.basename(__file__)} as standalone module...")

    # Run main function for data validation file generation
    gather_ctd2_data()