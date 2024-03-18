"""
build_validation_files.py
2024-03-18 ZD

This module defines primary function build_validation_files that generates files
useful for summarizing the data gathered in this pipeline. The files are intended
to act as the expected source of truth when QA testing INS data functionality.
"""

import os
import sys

import pandas as pd


# Append the project's root directory to the Python path
# This allows for importing config when running as part of main.py or alone
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import config


def get_single_node_counts(df: pd.DataFrame):
    """Get counts of unique ids and total ids within dataframe."""

    # Pull dataframe label from the 'type column
    node_type = df['type'][0]

    # Pull the node_id from the second position column
    id_field = df.columns[1]

    #print(F"Counting {node_type} values in '{id_field}' field...")

    # Get simple counts
    unique_ids = df[id_field].nunique()
    total_ids = df[id_field].count()

    # Get the number of ids present in more than one row
    # These all have different linking nodes
    id_occurrence = df[id_field].value_counts()
    ids_in_more_than_one_row = id_occurrence[id_occurrence > 1].count()

    return node_type, unique_ids, total_ids, ids_in_more_than_one_row



def get_all_node_counts(df_list: pd.DataFrame) -> pd.DataFrame:
    """
    This function iterates through a list of node dataframes to get id counts.

    Args:
        df_list (list): A list containing the DataFrames to analyze.

    Returns:
        pd.DataFrame: A new DataFrame containing id counts for each 
            DataFrame in df_list.
    """

    # Build empty dataframe with defined column headers
    columns = ['node_type', 
               'unique_ids', 
               'total_ids', 
               'ids_in_more_than_one_row']
    df_node_counts = pd.DataFrame(columns=columns)

    # Iterate through each df
    for df in df_list:

        # Get node counts using the existing function
        (node_type, 
         unique_ids, 
         total_ids, 
         ids_in_more_than_one_row) = get_single_node_counts(df)

        # Create a dictionary with results
        result_row = {'node_type': node_type,
                    'unique_ids': unique_ids,
                    'total_ids': total_ids,
                    'ids_in_more_than_one_row': ids_in_more_than_one_row}
        
        # Convert to df
        result_row_df = pd.DataFrame(result_row, index=[0])

        # Append the dictionary as a new row to the result DataFrame
        df_node_counts = pd.concat([df_node_counts, result_row_df], 
                                   ignore_index=True)

    return df_node_counts



def build_validation_files(): 
    """Generate summary files useful for data validation."""

    # Load data and make df_list

    # Get unique and total id summary counts

    # Test Program filtering
        # Get program id(s) from program properties (DOC, Focus Area)
        # Get project(s) from program id(s)
        # Get grant(s) from project(s)
        # Get publication(s) from project(s)
    
    # Format output Excel for QA team
