# summary_statistics.py
# 2023-09-08 ZD

# This script defines functions that will generate summary statistics relevant
# to data gathered for the INS project. The goal of these statistics are not 
# to ingest into the site, but rather for use in testing, validation, and 
# general reporting. 
# Summary statistics will be output to the reports/ directory with a versioning
# structure identical to the data/ directory. 

import pandas as pd
import os
from itertools import combinations
import sys
# Append the project's root directory to the Python path
# This allows for importing config when running as part of main.py or alone
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import config

def get_summary_statistics(all_grants_data:pd.DataFrame):
    """Create reports with summary statistics of high-level grants info"""

    # Define directory to store reports. Create if doesn't already exist
    reports_dir = config.REPORTS_DIR
    if not os.path.exists(reports_dir):
        os.makedirs(reports_dir)  

    # Summary 1: Group by Program and get aggregated grant stats
    get_grant_stats_by_program(all_grants_data)

    # Summary 2: Get pairs of Key Programs that share projects
    get_shared_projects_by_program_pair(all_grants_data)

    

def get_grant_stats_by_program(all_grants_data:pd.DataFrame):
    """Group grants by program and get aggregated summaries"""

    # Summarize columns of interest using aggregations defined in config
    grant_stats_by_program = (all_grants_data.groupby(
                                    config.PROGRAM_ID_FIELDNAME)
                                    .agg(config.STAT_AGG_FUNCS_BY_COL)
                                    .reset_index())
    
    # Rename fiscal year stat column for clarity
    fiscal_year_renamed = 'earliest_'+ config.STAT_FISCALYEAR_COL
    grant_stats_by_program.rename(columns={
                            config.STAT_FISCALYEAR_COL:fiscal_year_renamed},
                            inplace=True)

    # Export to reports
    grant_stats_by_program.to_csv(config.STAT_GRANTS_BY_PROGRAM_FILENAME, 
                                  index=False)
    

def get_shared_projects_by_program_pair(all_grants_data:pd.DataFrame):
    """Get pairs of Key Programs and the count of projects they share"""

    program_field = config.PROGRAM_ID_FIELDNAME
    project_field = config.STAT_CORE_PROJECT_COL

    # Group by core project and combine programs into unordered set
    grouped = all_grants_data.groupby(project_field)[program_field].apply(set)

    # Get combos of programs for each core project with itertools combinations
    df_shared_programs = grouped.apply(lambda x: list(combinations(x, 2)))

    # Flatten the program combos and count occurrences
    df_shared_programs = (df_shared_programs.explode()
                          .value_counts().reset_index())

    # Split program combo column into two separate columns
    # NOTE this will need to be reworked if a project has more than 2 programs
    df_shared_programs[['program_1','program_2']] = (df_shared_programs['index']
                                                     .apply(pd.Series))

    # Reorder and rename
    df_shared_programs.rename(columns={program_field:'shared_project_count'}, 
                              inplace=True)
    df_shared_programs = df_shared_programs[[
                                 'program_1',
                                 'program_2',
                                 'shared_project_count']]

    # Export as report
    shared_programs_filename = config.STAT_SHARED_PROJECT_PROGRAM_PAIRS_FILENAME
    df_shared_programs.to_csv(shared_programs_filename, index=False)


    # Run module as a standalone script when called directly
if __name__ == "__main__":

    print(f"Running {os.path.basename(__file__)} as standalone module...")

    # Load projects.tsv as a dataframe
    project_filename = os.path.join(config.PROCESSED_DIR, 'project.tsv')
    print(f"Generating report statistics on {project_filename}...")
    all_cleaned_grants = pd.read_csv(project_filename, sep='\t')

    # Run stats module
    get_summary_statistics(all_cleaned_grants)
    print(f"Summary reports for grants data successfully generated.\n"
          f"Results can be found in {config.REPORTS_DIR}.\n---")
