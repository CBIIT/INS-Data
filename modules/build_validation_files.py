"""
build_validation_files.py
2024-03-18 ZD

This module defines primary function build_validation_files that generates files
useful for summarizing the data gathered in this pipeline. The files are intended
to act as the expected source of truth when QA testing INS data functionality.

Data output TSVs will not be changed by this module. It will only read them in 
order to build validation files. 
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



def get_downstream_node_records(df_downstream: pd.DataFrame,
                                link_id: str,
                                link_values: str | list):
    """
    Generic function to pull all associated records from a single downstream 
    node with link values. Handles either string or list for link_values.
    
    Args:
        df_downstream (pd.DataFrame): Dataframe containing a link_id that can
            associate records with the upstream Dataframes.
        link_id (str): String name of column containing linking ids
        link_values (list | str): ID or list of IDs to use for filteirng records

    Returns: 
        linked_records (pd.Dataframe): Dataframe of all downstream records with
            matching linking values. 

    Example: 
        To records for all projects associated with a program, use args:
            df_downstream = df_projects,
            link_id = 'program.program_id' # Column within df_projects
            link_values = 'ccdi' # Program ID within program.program_id column
    """

    # If a single link_value is provided:
    if isinstance(link_values, str):  

        # Filter for exact string match
        linked_records = df_downstream[df_downstream[link_id] == link_values]

    # If a list of linking values is provided:
    else:

        # Filter for all records containing any link values
        linked_records = df_downstream[df_downstream[link_id].isin(link_values)]

    return linked_records




def build_single_program_output_counts(df_programs: pd.DataFrame,
                                        df_projects: pd.DataFrame,
                                        df_grants: pd.DataFrame,
                                        df_publications: pd.DataFrame):
    """
    Generates a dataframe summarizing counts of outputs (projects, grants,
        and publications) for each program. 

    Args:
        df_programs (pd.DataFrame): Dataframe from program.tsv
        df_projects (pd.DataFrame): Dataframe from project.tsv
        df_grants (pd.DataFrame): Dataframe from grant.tsv
        df_publications (pd.DataFrame): Dataframe from publication.tsv
    """
    
    # Build empty dataframe to fill with results
    df_program_results = pd.DataFrame()

    # Get list of all unique program ids
    program_id_list = df_programs['program_id'].unique().tolist()

    # Iterate through each program
    for program_id in program_id_list:

        # Get single program record from programs df
        df_program_filtered = (df_programs[
                                df_programs['program_id'] == program_id]
                                ).reset_index()

        # Get program name for output
        program_name = df_program_filtered['program_name']
        
        # Get program acronym for output
        program_acronym = df_program_filtered['program_acronym']

        # Get projects linked to program
        linked_projects = get_downstream_node_records(df_projects, 
                                                      'program.program_id', 
                                                       program_id)
        project_id_list = linked_projects['project_id'].unique() 

        # Get grants linked to project (indirectly to program)
        linked_grants = get_downstream_node_records(df_grants,
                                                    'project.project_id',
                                                     project_id_list)
        grant_id_list = linked_grants['grant_id'].unique()

        # Get publications linked to project (indirectly to program)
        linked_publications = get_downstream_node_records(df_publications, 
                                                          'project.project_id', 
                                                           project_id_list)
        publication_id_list = linked_publications['pmid'].unique()

        # Create a dictionary with results
        program_result_row = {
            'program_id': program_id,
            'program_name': program_name,
            'program_acronym': program_acronym,
            'projects': len(project_id_list),
            'grants': len(grant_id_list),
            'publications': len(publication_id_list)
            }
        
        # Convert to df
        result_row_df = pd.DataFrame(program_result_row, index=[0])

        # Add to running list of full results
        df_program_results = pd.concat([df_program_results, result_row_df], 
                                       ignore_index=True)
    
    return df_program_results



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
