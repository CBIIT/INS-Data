"""
build_validation_file.py
2024-03-18 ZD

This module defines primary function build_validation_file that generates a file
useful for summarizing the data gathered in this pipeline. The file is intended
to act as the expected source of truth when QA testing INS data functionality.

Data output TSVs will not be changed by this module. It will only read them in 
order to build the validation file. 
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
    Iterates through a list of node dataframes to get summary counts of ids for 
        each node. Each node/dataframe is a specific data type (e.g. programs).

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
    id_count_summary_df = pd.DataFrame(columns=columns)

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
        id_count_summary_df = pd.concat([id_count_summary_df, result_row_df], 
                                   ignore_index=True)

    return id_count_summary_df



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
        To get records for all projects associated with a program, use args:
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

        # Get total of all downstream records
        total_records = sum([len(project_id_list),
                            len(grant_id_list),
                            len(publication_id_list),
                            ])

        # Get detail page urls            
        dev_url = get_detail_page_url(program_id, 'program', '-dev')
        qa_url = get_detail_page_url(program_id, 'program', '-qa')
        stage_url = get_detail_page_url(program_id, 'program', '-stage')
        prod_url = get_detail_page_url(program_id, 'program', '')

        # Create a dictionary with results
        program_result_row = {
            'program_id': program_id,
            'program_name': program_name,
            'program_acronym': program_acronym,
            'projects': len(project_id_list),
            'grants': len(grant_id_list),
            'publications': len(publication_id_list),
            'total_records': total_records,
            'dev_url': dev_url,
            'qa_url': qa_url,
            'stage_url': stage_url,
            'prod_url': prod_url,
        }
        
        # Convert to df
        result_row_df = pd.DataFrame(program_result_row, index=[0])

        # Add to running list of full results
        df_program_results = pd.concat([df_program_results, result_row_df], 
                                       ignore_index=True)
        
    # Sort output df by total records
    df_program_results.sort_values(by='total_records', ascending=False, 
                                   inplace=True)
    
    return df_program_results



def build_single_project_output_counts(df_projects: pd.DataFrame,
                                       df_grants: pd.DataFrame,
                                       df_publications: pd.DataFrame):
    """
    Generates a dataframe summarizing counts of outputs (grants and publications)
    for each project.

    Args:
        df_projects (pd.DataFrame): Dataframe from project.tsv
        df_grants (pd.DataFrame): Dataframe from grant.tsv
        df_publications (pd.DataFrame): Dataframe from publication.tsv
    """

    # Build empty dataframe to store results
    df_project_results = pd.DataFrame()

    # Get list of all unique project ids
    project_id_list = df_projects['project_id'].unique().tolist()

    # Iterate through each project id
    for project_id in project_id_list:

        # Get single project record from projects df
        df_project_filtered = (df_projects[
                                df_projects['project_id'] == project_id]
                                ).reset_index()

        # Get project title for output
        project_title = df_project_filtered['project_title']

        # Get grants linked to project
        linked_grants = get_downstream_node_records(df_grants,
                                                    'project.project_id',
                                                    project_id)
        grant_id_list = linked_grants['grant_id'].unique()

        # Get publications linked to project
        linked_publications = get_downstream_node_records(df_publications,
                                                          'project.project_id',
                                                          project_id)
        publication_id_list = linked_publications['pmid'].unique()

        # Get total of all downstream records
        total_records = sum([len(grant_id_list),
                            len(publication_id_list),
                            ])
        
        # Get detail page urls            
        dev_url = get_detail_page_url(project_id, 'project', '-dev')
        qa_url = get_detail_page_url(project_id, 'project', '-qa')
        stage_url = get_detail_page_url(project_id, 'project', '-stage')
        prod_url = get_detail_page_url(project_id, 'project', '')

        # Create a dictionary with project results
        project_result_row = {
            'project_id': project_id,
            'project_title': project_title,
            'grants': len(grant_id_list),
            'publications': len(publication_id_list),
            'total_records': total_records,
            'dev_url': dev_url,
            'qa_url': qa_url,
            'stage_url': stage_url,
            'prod_url': prod_url,
        }

        # Convert to dataframe format
        result_row_df = pd.DataFrame(project_result_row, index=[0])

        # Add to running list of project results
        df_project_results = pd.concat([df_project_results, result_row_df],
                                       ignore_index=True)
    # Sort output df by total records
    df_project_results.sort_values(by='total_records', ascending=False, 
                                   inplace=True)

    return df_project_results



def get_detail_page_url(node_id, node_type, tier):
    """Build the URL to the relevant program or project detail page on INS for
    easier QA testing.

    Args:
        node_id (str): ID of the program or project
        node_type (str): 'program' or 'project' describing the node
        tier (str): Site tier for url ('-dev', '-qa', '-stage', or None (prod))
    """

    if tier == '-prod' or '' or None:
        tier = ''

    url = f"https://studycatalog{tier}.cancer.gov/#/{node_type}/{node_id}"

    return url



def get_downstream_node_ids_from_list(df_downstream: pd.DataFrame, 
                                        link_id: str, 
                                        id_list: list, 
                                        target_id_col: str) -> int:
    """Input a list of upstream node ids and get the combined output ids. If 
    downstream outputs are duplicated, only report unique values. 

    Args:
        df_downstream (pd.DataFrame): Dataframe containing a link_id that can
            associate records with the upstream Dataframes
        link_id (str): String name of column containing linking ids
        id_list (list): List of ids to gather downstream data from 
        target_id_col (str): The name of the column containing the target IDs 
        whose unique values we want to find

    Returns:
        unique_id_list (list): List of all associated downstream node ids

    Example:
        To get counts for all projects associated with a list of programs, use:
            df_downstream = df_projects 
            link_id = 'program.program_id'      # Column within df_projects
            id_list = ['ccdi', 'pivot']         # list of program ids
            target_id_col (str): 'project_id'   # ID column within df_projects    
    """

    # Build empty dataframe to fill with results
    combined_records = pd.DataFrame()

    # Iterate through each id within the input list
    for node_id in id_list:
        
        # Get downstream records for individual node id
        linked_records = get_downstream_node_records(df_downstream, 
                                                        link_id, 
                                                        node_id)
        # Combine with running list of combined outputs
        combined_records = pd.concat([combined_records, linked_records], 
                                     ignore_index=True)
    
    # Keep only unique output ids
    unique_id_list = combined_records[target_id_col].unique().tolist()

    return unique_id_list



def build_multi_program_output_counts(program_list: list,
                                      df_projects: pd.DataFrame,
                                      df_grants: pd.DataFrame,
                                      df_publications: pd.DataFrame):
    
    """Generates a dataframe summarizing counts of programs and outputs 
    (projects, grants, and publications) for a list of programs.
    """

    # Get unique ids for all associated projects (from programs)
    project_id_list = get_downstream_node_ids_from_list(df_projects,
                                        'program.program_id',
                                        program_list,
                                        'project_id')

    # Get unique ids for all associated grants (from projects)
    grant_id_list = get_downstream_node_ids_from_list(df_grants,
                                        'project.project_id',
                                        project_id_list,
                                        'grant_id')

    # Get unique ids for all associated publications (from projects)
    publication_id_list = get_downstream_node_ids_from_list(df_publications,
                                        'project.project_id',
                                        project_id_list,
                                        'pmid')

    # Create a dictionary with results
    result_row = {
        'programs': len(program_list),
        'projects': len(project_id_list),
        'grants': len(grant_id_list),
        'publications': len(publication_id_list),
    }

    return result_row



def get_single_filter_value_results(filter_field: str, 
                                    filter_value: str,
                                    node_df: pd.DataFrame,
                                    node_id_col: str):
    
    """Retrieve a list of all node ids matching a single filter field value.
    Note: This will need to be reworked for non-string filter values.

    Args:
        filter_field (str): Filter field to search for value within
        filter_value (str): Filter value to match within filter_field
        node_df (pd.DataFrame): Node dataframe
        node_id (str): Column of node id

    Returns:
        id_list (list): List of all node ids matching the filter value

    Example:
        To get all program ids associated with Pancreas Cancer cancer_type, 
            use parameters:
            filter_field = 'cancer_type',
            filter_value = 'Pancreas Cancer',
            node_df = df_programs,
            node_id_col = 'program_id'
    """

    # Find matching filter strings
    df_filtered = node_df[node_df[filter_field].str.contains(filter_value)]

    # Get list of unique associated program ids
    id_list = df_filtered[node_id_col].unique().tolist()

    return id_list



def build_program_filter_output_counts(filter_field_list: list, 
                                       df_programs: pd.DataFrame, 
                                       df_projects: pd.DataFrame,
                                       df_grants: pd.DataFrame,
                                       df_publications: pd.DataFrame,
                                       max_values: int=None):
    
    """Build a dataframe of program and output counts associated with each 
    single facet filter selection.

    Args:
        df_programs (pd.DataFrame): Dataframe of programs
        filter_field_list (list): List of all columns to include as filters
        max_values (int, optional): Maximum number of values to include for 
            fields with many values. Default None to include all.
    
    Returns:
        filter_result_df (pd.DataFrame): Dataframe with columns for filter 
            labels and expected program, project, grant, and publication counts. 
    """

    # Build empty dataframe to fill with results
    filter_result_df = pd.DataFrame()

    # Iterate through list of columns to consider filters
    for filter_field in filter_field_list:
        
        # Separate list-like strings and get unique values
        filter_value_list = (df_programs[filter_field].str.split(';')
                                                .explode().unique().tolist())
        
        # Sort field alphabetically for consistency. Ignore case
        filter_value_list = sorted(filter_value_list, key=str.casefold)

        # Limit list size, if specified. Will select first alphabetically.
        if max_values is not None:
            filter_value_list = filter_value_list[0:max_values]
        
        # Iterate through each filter value option
        for filter_value in filter_value_list:
            
            # Get all programs associated with single filter
            program_list = get_single_filter_value_results(filter_field, 
                                                           filter_value, 
                                                           df_programs,
                                                           'program_id')

            # Get dictionary of output counts for program list
            result_row = build_multi_program_output_counts(program_list, 
                                                           df_projects, 
                                                           df_grants, 
                                                           df_publications)

            # Add labels for filter field and value to the output
            result_row['filter_field'] = filter_field
            result_row['filter_value'] = filter_value

            # Convert to df
            result_row_df = pd.DataFrame(result_row, index=[0])

            # Append as a new row in result DataFrame
            filter_result_df = pd.concat([filter_result_df, result_row_df], 
                                         ignore_index=True)
            
    # Rearrange columns in output dataframe for testing clarity
    filter_result_df = filter_result_df[['filter_field','filter_value',
                                         'programs','projects',
                                         'grants','publications']]

    return filter_result_df



def save_dataframes_to_excel(df_dict: dict, output_file: str):
    """
    Saves a dictionary of DataFrames as tabs within an Excel file. 
    including a report_info tab with information from a provided dictionary.

    Args:
        df_dict (dict): Dictionary containing DataFrames to save as tabs.
            Keys are sheet names, values are the DataFrames.
        output_file (str): Path to the output Excel file.
    """

    # Initialize excel writer
    writer = pd.ExcelWriter(output_file, engine='xlsxwriter')

    # Save each DataFrame to a separate sheet
    for sheet_name, df in df_dict.items():
        df.to_excel(writer, sheet_name=sheet_name, index=False)

    # Format width of (descriptive) first sheet columns for readability
    writer.sheets[list(df_dict.keys())[0]].set_column(0, 0, 40)
    writer.sheets[list(df_dict.keys())[0]].set_column(1, 1, 20)

    writer.close()



def build_validation_file(): 
    """Generate summary Excel file useful for data validation."""

    print(f"\n---\nDATA VALIDATION:\n"
        f"Generating file for data validation...\n---\n")

    # Load data from output TSVs
    df_programs = pd.read_csv(config.PROGRAMS_OUTPUT_PATH, sep='\t')
    df_projects = pd.read_csv(config.PROJECTS_OUTPUT_PATH, sep='\t')
    df_grants = pd.read_csv(config.GRANTS_OUTPUT_PATH, sep='\t')
    df_publications = pd.read_csv(config.PUBLICATIONS_OUTPUT_PATH, sep='\t')
    
    # Define list of dataframes for iteration
    df_list = [df_programs, df_projects, df_grants, df_publications]

    # Get unique and total ID counts for each dataframe
    id_count_summary_df = get_all_node_counts(df_list)

    # Get associated project, grant, and publication count for each program
    df_single_program_results = build_single_program_output_counts(
                                                    df_programs, 
                                                    df_projects, 
                                                    df_grants, 
                                                    df_publications)
    
    # Get associated grant and publication count for each project
    df_single_project_results = build_single_project_output_counts(
                                                    df_projects, 
                                                    df_grants, 
                                                    df_publications)

    # Define version and descriptive info to include on first tab of Excel
    # There's definitely a better way to format this, but it works...
    report_info = pd.DataFrame({
        'Qualtrics Version': config.QUALTRICS_VERSION,
        'Data Gathering Date': config.TIMESTAMP,
        'iCite Version': config.ICITE_VERSION,
        'Qualtrics Type': config.QUALTRICS_TYPE,
        '---': '---',
        'This report generated': pd.Timestamp.now(),
        }.items(),
        columns=['INS Data Validation Report', ''])

    # Define Excel tab names and dataframes for each tab
    df_dict = {
        'report_info': report_info,
        'total_counts': id_count_summary_df,
        'single_program_counts': df_single_program_results,
        'single_project_counts': df_single_project_results,
        }

    # Export file
    save_dataframes_to_excel(df_dict, config.DATA_VALIDATION_EXCEL)
    print(f"Done! Data validation file saved to {config.DATA_VALIDATION_EXCEL}.")
    


# Run module as a standalone script when called directly
if __name__ == "__main__":

    print(f"Running {os.path.basename(__file__)} as standalone module...")

    # Run main function for data validation file generation
    build_validation_file()