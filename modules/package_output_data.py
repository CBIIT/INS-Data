"""
package_output_data.py
2023-12-21 ZD

This script defines primary function package_output_data that applies final 
formatting, filtering, etc. to the intermediate data gathering outputs in order
to prepare them for INS ingestion. These outputs are saved as TSVs. 
"""

import os
import sys
import unicodedata

import pandas as pd
import re

# Append the project's root directory to the Python path
# This allows for importing config when running as part of main.py or alone
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import config

# On hold until INS-807
# def remove_publications_before_projects(df_publications, df_projects):
#     """Remove publications published before the associated project date"""
#     # Placeholder for code
#     df_publications_filtered = df_publications # do stuff
#     return df_publications_filtered



def add_type_column(df, datatype):
    """Add a column for datatype and fill with specified datatype."""

    df_with_type = df.copy()
    df_with_type['type'] = datatype

    return df_with_type



def reorder_columns(df, column_configs, datatype):
    """Reorder and validate columns of df using dict from config.py. 
    Report any dropped columns and throw error if defined column not found.
    """
    # Check if the datatype exists in column_configs
    if datatype not in column_configs:
        raise ValueError(f"Invalid datatype: {datatype}")

    # Extract relevant configuration for the datatype
    config = column_configs[datatype]

    # Extract relevant information from the configuration
    keep_and_rename = config.get('keep_and_rename', {})

    # Check if all desired columns are present before renaming
    if not set(keep_and_rename.keys()).issubset(set(df.columns)):
        missing_columns = set(keep_and_rename.keys()) - set(df.columns)
        raise ValueError(f"Fields missing from intermediate {datatype} data: "
                         f"{', '.join(missing_columns)}")

    # Rename columns using keep_and_rename
    df.rename(columns=keep_and_rename, inplace=True)

    # Define list of columns to keep
    columns_to_keep = list(keep_and_rename.values())

    # Drop any columns not part of keep_and_rename
    dropped_columns = set(df.columns) - set(columns_to_keep)
    if dropped_columns:
        print(f"Columns dropped from {datatype} output: "
              f"{', '.join(dropped_columns)}")

    # Keep and reorder columns
    df_reordered = df[columns_to_keep]

    return df_reordered



def validate_first_columns(df, column_configs, datatype):
    """Validate that type, node_id, and link_id (opt) are first columns."""

    # Check if the datatype exists in column_configs and then get it
    if datatype not in column_configs:
        raise ValueError(f"Invalid datatype: {datatype}")
    config = column_configs.get(datatype)

    # Get values from column configuration
    node_id = config.get('node_id')
    link_id = config.get('link_id', None)

    # If link_id exists, add it to the list of expected columns 
    expected_cols = ['type', node_id]
    if link_id:
        expected_cols.append(link_id)

    # Check the order of columns in df
    if list(df.columns[:len(expected_cols)]) != expected_cols:
        raise ValueError(f"First columns in {datatype} output are not in the "
                         f"expected order. Reorder columns in config.py to fix.\n"
                         f"Expected order: {', '.join(expected_cols)}.\n"
                         f"Actual order:   "
                         f"{', '.join(df.columns[:len(expected_cols)])}")
    
    return None



def replace_defined_characters(text):
    """Use a mappping to replace specific non-standard characters in text."""

    # Return NaN and non-string values as-is
    if pd.isna(text) or not isinstance(text, str):
        return text

    # Define translation table to replace characters
    translation_table = str.maketrans({
        '\u2019': "'",  # Apostrophe with diacritic - replace with apostrophe
        '\u2028': ' ',  # Line Separator (LS) - replace with space
        '\u2029': ' ',   # Paragraph Separator (PS) - replace with space
        # Add mappings here in the future as needed
        })

    # Use translate only for strings
    replaced_text = (text.translate(translation_table) 
                     if isinstance(text, str) else text)
    
    # Use translate only for strings
    return replaced_text



def normalize_encoding(text):
    """Normalizes text encoding for consistency and compatibility."""

    # Return NaN and non-string values as-is
    if pd.isna(text) or not isinstance(text, str):
        return text

    # Encode and decode text into ASCII to normalize
    try:
        formatted_text = (unicodedata.normalize('NFKC', text)
                          .encode('ascii', 'ignore').decode('ascii'))
        return formatted_text
    
    # Second catch to handle non-string values as-is
    except TypeError:
        return text



def process_special_characters(df):
    """Cleans a DataFrame by normalizing special characters and encoding."""

    # Replace specific non-standard characters
    df = df.map(replace_defined_characters)
    # General normalization to ascii encoding
    df_cleaned = df.map(normalize_encoding)

    return df_cleaned



def format_list_like_columns(df, column_configs):
    """Use a list of column names from config.py to structure list-like strings
    into semicolon-separated strings. Remove whitespace and/or delimiters like 
    commas."""
    # Placeholder for code
    df_structured_lists = df # do stuff

    return df_structured_lists



def validate_and_clean_unique_nodes(df, column_configs, datatype):
    """
    Validate that all node_ids are either unique values or valid duplicates.

    Args:
        df (pd.DataFrame): Input DataFrame
        column_configs (dict): Configuration settings for columns, including 
            node_id and link_id
        datatype (str): Type of data in the DataFrame. Indicates which info
            will be accessed from column_configs

    Returns:
        pd.DataFrame: A new DataFrame with invalid duplicate rows removed.

    Raises:
        ValueError: If the specified node_id column is not found in the DataFrame.

    Notes:
        This function checks for valid uniqueness or duplicates in node_ids 
            based on the provided configuration.

        - For many-to-many nodes, it ensures that any duplicates only differ 
            in the link_id column.
        - True duplicates (identical values in all columns) are kept in the 
            first row for production.
    """

    # Get node_id and link_id from column configurations
    node_id = column_configs[datatype].get('node_id')
    link_id = column_configs[datatype].get('link_id', None)

    # Check if the specified node_id column exists
    if node_id not in df.columns:
        raise ValueError(f"Node ID column '{node_id}' not found in DataFrame.")

    # Get only rows with duplicated node_ids 
    dup_nodes = df[df.duplicated(subset=node_id, keep=False)]
    
    # Get true duplicates where values in all columns are identical
    true_duplicates = dup_nodes[dup_nodes.duplicated(keep=False)].assign(
        error_reason=f"True duplicate in all columns. "
                     f"Kept first row for production.")
    
    # Empty dataframe to store mismatch rows
    mismatch_df = pd.DataFrame()

    # Get list of value columns to check for mismatch
    value_cols = list(set(df.columns) - set([node_id, link_id]))

    # Group rows by node_id for shared value checking
    for node, node_df in dup_nodes.groupby(node_id):
        
        # Check for mismatched values in other columns within the same node
        value_mismatches = node_df[~node_df.duplicated(
            subset=value_cols, keep=False)
            ].assign(error_reason=f"Mismatched values in columns with same "
                                  f"node_id. Both dropped for production.")

        # Add to running list of rows with mismatched values
        mismatch_df = pd.concat([mismatch_df, value_mismatches])

    # Combine list of annotated invalid duplicates for export to csv
    invalid_duplicates_df = pd.concat([true_duplicates, mismatch_df])

    # Save any invalid duplicate rows to a report csv for troubleshooting
    if len(invalid_duplicates_df) > 0:
        invalid_duplicates_path = f"{config.REMOVED_DUPLICATES}{datatype}.csv"
        os.makedirs(os.path.dirname(invalid_duplicates_path), exist_ok=True)
        invalid_duplicates_df.to_csv(invalid_duplicates_path, index=False)
        print(f"Invalid duplicate rows cleaned from output and recorded in "
              f"{invalid_duplicates_path}.")

    # Return a cleaned df with only unique node_ids or valid duplicates
    # Remove all mismatched node_ids and any true duplicates after the first
    valid_df = df[~df[node_id].isin(mismatch_df[node_id]) 
                  & ~df.duplicated(subset=node_id, keep='first')
                  ].reset_index(drop=True)

    return valid_df



def standardize_data(df, column_configs, datatype):
    """Group standardization functions common to all data types"""

    # Edit data to standardize
    df = add_type_column(df, datatype)
    df = reorder_columns(df, column_configs, datatype)
    df = process_special_characters(df)
    df = format_list_like_columns(df, column_configs)
    df = validate_and_clean_unique_nodes(df, column_configs, datatype)

    # Validate that data meets loading standards
    validate_first_columns(df, column_configs, datatype)
    
    return df



# def package_programs(df_programs, column_configs):
#     """Placeholder
#     """
#     # Placeholder for code
#     print(f"Finalizing TSV for project data...")

#     #df_programs_output = standardize_data(df_programs, column_configs, datatype='program')

#     # During dev:
#     df_programs = add_type_column(df_programs, datatype='program')
#     df_programs = process_special_characters(df_programs)

#     df_programs_output = df_programs

#     # Export as TSV
#     output_filepath = config.PROGRAMS_OUTPUT_PATH
#     os.makedirs(os.path.dirname(output_filepath), exist_ok=True)
#     df_programs_output.to_csv(output_filepath, sep='\t', index=False, encoding='utf-8')

#    print(f"Done! Program data exported to TSV.")

#     return df_programs_output



def package_grants(df_grants, column_configs):
    """Placeholder
    """
    print(f"---\nFinalizing TSV for grant data...")

    # Placeholder for code
    df_grants_output = standardize_data(df_grants, column_configs, datatype='grant')

    # Export as TSV
    output_filepath = config.PROJECTS_OUTPUT_PATH
    os.makedirs(os.path.dirname(output_filepath), exist_ok=True)
    df_grants_output.to_csv(output_filepath, sep='\t', index=False, encoding='utf-8')

    print(f"Done! Final grant data saved as {output_filepath}.")

    return df_grants_output



# def package_projects(df_projects, column_configs):
#     """Placeholder
#     """
#     print(f"Finalizing TSV for project data...")

#     # Placeholder for code
#     df_projects_output = standardize_data(df_projects, column_configs, datatype='project')

#    print(f"Done! Project data exported to TSV.")

#     return df_projects_output # Also export as TSV



# def package_publications(df_publications, column_configs):
#     """Placeholder
#      """
#     print(f"Finalizing TSV for publication data...")

#     df_publications_output = standardize_data(df_publications, column_configs, datatype='publication')

#    print(f"Done! Publication data exported to TSV.")

#     return df_publications_output # Also export as TSV



def package_output_data():
    """Run all data packaging steps for all data types."""

    # Single data model dict defined in config
    column_configs = config.COLUMN_CONFIGS 

    # Load all files as dfs
    # df_programs = pd.read_csv(config.CLEANED_KEY_PROGRAMS_CSV)
    df_grants = pd.read_csv(config.PROJECTS_INTERMED_PATH)
    # df_projects = pd.read_csv(path_to_CSV)
    # df_publications = pd.read_csv(config.PUBLICATIONS_INTERMED_PATH)
    print(f"Loaded files from {config.GATHERED_DIR}.")

    # Special handling
    # df_publications = remove_publications_before_projects(df_publications, 
    #                                                       df_projects)

    # Final packaging
    # df_programs_output = package_programs(df_programs, column_configs)
    df_grants_output = package_grants(df_grants, column_configs)
    # df_projects_output = package_projects(df_projects, column_configs)
    # df_publications_output = package_publications(df_publications, column_configs)

    # Return a dictionary of DataFrames
    return {
        # 'programs_output': df_programs_output,
        'grants_output': df_grants_output,
        # 'projects_output': df_projects_output,
        # 'publications_output': df_publications_output
    }



# Run module as a standalone script when called directly
if __name__ == "__main__":

    print(f"Running {os.path.basename(__file__)} as standalone module...")

    # Placeholder for code
    outputs = package_output_data()

    # Access DataFrames using keys
    # df_programs_output = outputs['programs_output']
    df_grants_output = outputs['grants_output']
    # df_projects_output = outputs['projects_output']
    # df_publications_output = outputs['publications_output']