"""
package_output_data.py
2023-12-21 ZD

This script defines primary function package_output_data that applies final 
formatting, filtering, etc. to the intermediate data gathering outputs in order
to prepare them for INS ingestion. These outputs are saved as TSVs. 
"""

import os
import sys

import pandas as pd
import re

# Append the project's root directory to the Python path
# This allows for importing config when running as part of main.py or alone
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import config


def remove_publications_before_projects(df_publications, df_projects):
    """Remove publications published before the associated project date"""
    # Placeholder for code
    df_publications_filtered = df_publications # do stuff
    return df_publications_filtered



def reorder_columns(df, column_configs):
    """Reorder and validate columns of df using dict from config.py. 
    Report any dropped columns and throw error if defined column not found.
    """
    # Placeholder for code
    df_reordered = df # do stuff
    return df_reordered



def remove_special_characters(df):
    """Remove any special characters, unusuaul breakpoints, etc. from all values"""
    # Placeholder for code
    df_standard_characters = df # do stuff
    return df_standard_characters



def format_list_like_columns(df, column_configs):
    """Use a list of column names from config.py to structure list-like strings
    into semicolon-separated strings. Remove whitespace and/or delimiters like 
    commas."""
    # Placeholder for code
    df_structured_lists = df # do stuff

    return df_structured_lists



def validate_unique_node_id(df, column_configs):
    """Validate that each row with a matching node_id has identical values for 
    all other columns except for link_id"""
    # Placeholder for code
    return df # If all good, throw error if not



def standardize_data(df, column_configs):
    """Group standardization functions common to all data types"""
    df = reorder_columns(df, column_configs)
    df = remove_special_characters(df)
    df = format_list_like_columns(df, column_configs)
    df = validate_unique_node_id(df, column_configs)

    return df


def package_programs(df_programs, column_configs):
    """Placeholder
    """
    # Placeholder for code
    df_programs_output = standardize_data(df_programs)

    return df_programs_output # Also export as TSV



def package_grants(df_grants, column_configs):
    """Placeholder
    """
    # Placeholder for code
    df_grants_output = standardize_data(df_grants)

    return df_grants_output # Also export as TSV



def package_projects(df_projects, column_configs):
    """Placeholder
    """
    # Placeholder for code
    df_projects_output = standardize_data(df_projects)

    return df_projects_output # Also export as TSV



def package_publications(df_publications, column_configs):
    """Placeholder
     """
    
    df_publications_output = standardize_data(df_publications)

    return df_publications_output # Also export as TSV



def package_output_data():
    # Placeholder to run all data packaging steps

    column_configs = None # Single data model dict defined in config
    path_to_CSV = None # Defined in config for each file

    # Load all files as dfs
    df_programs = pd.read_csv(path_to_CSV)
    df_grants = pd.read_csv(path_to_CSV)
    df_projects = pd.read_csv(path_to_CSV)
    df_publications = pd.read_csv(path_to_CSV)

    # Special handling
    df_publications = remove_publications_before_projects(df_publications, 
                                                          df_projects)

    # Final packaging
    df_programs_output = package_programs(df_programs, column_configs)
    df_grants_output = package_grants(df_grants, column_configs)
    df_projects_output = package_projects(df_projects, column_configs)
    df_publications_output = package_publications(df_publications, column_configs)

    # Return a dictionary of DataFrames
    return {
        'programs_output': df_programs_output,
        'grants_output': df_grants_output,
        'projects_output': df_projects_output,
        'publications_output': df_publications_output
    }



# Run module as a standalone script when called directly
if __name__ == "__main__":

    print(f"Running {os.path.basename(__file__)} as standalone module...")

    # Placeholder for code
    outputs = package_output_data()

    # Access DataFrames using keys
    df_programs_output = outputs['programs_output']
    df_grants_output = outputs['grants_output']
    df_projects_output = outputs['projects_output']
    df_publications_output = outputs['publications_output']