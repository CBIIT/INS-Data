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



def gather_ctd2_datasets():
    """
    Process a CSV of CTD^2 dataset metadata into an intermediate CSV for INS.
    Input and output filepaths are handled in config.py
    """

    print(f"\n---\nDATASETS - CTD^2:\n"
          f"Processing CTD^2 dataset data...\n---\n")

    # Define paths from config
    input_csv = config.CTD2_DATASET_INPUT_CSV
    output_csv = config.CTD2_DATASET_INTERMED_CSV

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
    
    return df



def gather_ctd2_filedata():
    """
    Process a CSV of file metadata into an intermediate CSV for INS dataloading.
    Input and output filepaths are handled in config.py
    """

    print(f"\n---\nFILES - CTD^2:\n"
        f"Processing CTD^2 file metadata...\n---\n")

    # Define paths from config
    input_csv = config.CTD2_FILE_INPUT_CSV

    # Load curated CTD^2 CSV
    df = pd.read_csv(input_csv)
    print(f"Loaded CTD^2 file metadata from {input_csv}")

    # Add reproducible uuid5
    uuid_fields = ['file_name', 'file_type', 'access_level']
    df = get_composite_uuid5(df, uuid_fields, uuid_col='file_id')

    return df



def gather_ctd2_data():
    """
    Main function for the CTD^2 data preparation process.

    Processes the dataset and file metadata inputs to generate aligned
    intermediate CSVs for the CTD^2 datasets and files. Maps files to their
    corresponding datasets using curated download file links, validates the
    mapping, and saves the processed file metadata.

    Steps performed:
        1. Process the CTD^2 datasets and file input CSVs into intermediates
        2. Map the files to datasets using curated download_file_link fields
        3. Validate the mappings

    Returns:
        pd.DataFrame: DataFrame containing the processed CTD^2 file metadata,
        including mapped dataset UUIDs.

    Raises:
        ValueError: If duplicate UUIDs are generated in the process.

    """
    # Process datasets input file, save csv, and return df
    datasets_df = gather_ctd2_datasets()

    # Process filedata input file and return df
    filedata_df = gather_ctd2_filedata()

    # FILE-DATASET MAPPING:
    # Map each file to a dataset uuid using the curated download_file_links

    print(f"\nMapping CTD^2 files to datasets...\n")

    # Split semicolon-separated download_file_links into new rows
    exploded_df = (datasets_df
        .assign(download_file_links=datasets_df['download_file_links']
                .str.split(';'))
        .explode('download_file_links')
        .dropna(subset=['download_file_links']))
    
    # Strip whitespace for reliable matching
    exploded_df['download_file_links'] = exploded_df['download_file_links'
                                                     ].str.strip()
    filedata_df['file_name'] = filedata_df['file_name'].str.strip()

    # Map file_name (from filedata_df) to download_file_links (from exploded_df)
    # Assumes that file_name matches exactly one of the download_file_links values

    mapping = exploded_df[['download_file_links', 'dataset_uuid']
                            ].drop_duplicates()

    filedata_df = filedata_df.merge(mapping,
                                    left_on='file_name',
                                    right_on='download_file_links',
                                    how='left')
    
    # Rename linking field in filedata
    filedata_df.rename(columns={'dataset_uuid': 'dataset.dataset_uuid'},
                       inplace=True)


    # MAPPING VALIDATION

    # 1. Mapped file count
    mapped_files = filedata_df[filedata_df['dataset.dataset_uuid'].notna()]
    print(f"Mapped file count: {len(mapped_files)}")

    # 2. Unmapped file count
    unmapped_files = filedata_df[filedata_df['dataset.dataset_uuid'].isna()]
    print(f"\nUnmapped file count: {len(unmapped_files)}")

    if not unmapped_files.empty:
        print("Unmapped file names:")
        print(unmapped_files['file_name'].to_list())

    # 3. Datasets with no mapped file counts
    mapped_dataset_uuids = set(mapped_files['dataset.dataset_uuid'].dropna())
    all_dataset_uuids = set(datasets_df['dataset_uuid'])
    datasets_with_no_files = all_dataset_uuids - mapped_dataset_uuids

    print(f"\nDatasets with no mapped files: {len(datasets_with_no_files)}")

    if datasets_with_no_files:
        print("Dataset titles with no mapped files:")
        titles = datasets_df[datasets_df['dataset_uuid'].isin(
            datasets_with_no_files)]['dataset_title'].tolist()
        for title in titles:
            print(title)

    # 4. Files mapped to more than 1 dataset
    file_dataset_counts = (filedata_df
                           .groupby('file_name')['dataset.dataset_uuid']
                           .nunique())
    files_with_multiple_datasets = file_dataset_counts[file_dataset_counts > 1]

    print(f"\nFiles mapped to more than one dataset: "
          f"{len(files_with_multiple_datasets)}")
    
    if not files_with_multiple_datasets.empty:
        print("Files mapped to multiple datasets:")
        print(files_with_multiple_datasets)


    # Save the updated filedata_df
    output_csv = config.CTD2_FILE_INTERMED_CSV
    os.makedirs(os.path.dirname(output_csv), exist_ok=True)
    filedata_df.to_csv(output_csv, index=False)

    print(f"\n\nSuccess! CTD^2 filedata saved to {output_csv}.\n"
          f"Total CTD^2 file records:    {len(filedata_df)}")

    return filedata_df


# Run module as a standalone script when called directly
if __name__ == "__main__":

    print(f"Running {os.path.basename(__file__)} as standalone module...")

    # Run main function for data validation file generation
    gather_ctd2_data()