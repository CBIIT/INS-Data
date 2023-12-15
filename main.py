# main.py
# 2023-07-26 ZD
#
# Main function for the INS-Data project. This will accept a Qualtrics CSV of 
# Key Programs as input and generate the following outputs:
#
#   - data/cleaned/: 
#       - Clean `key_programs_{date}.csv`
#   - data/processed/{qualtrics version}/{api_gathering_date}:
#       - For loading into INS: 
#               - project.tsv' containing data for grants associated with 
#                       key programs
#               - 'publication.tsv' containing data for publication PMIDs 
#                       associated with core projects
#       - Intermediates for troubleshooting:
#               - 'projectPMIDs.csv' CSV of associated PMIDs for each project
#               - 'temp-publication-data/' CSV directory of chunk-sized PubMed  
#                       data gathered for PMIDs of interest
#               - 'icitePMIDData.csv' CSV of iCite data for PMIDs of interest
#               - 'mergedPMIDData.csv' CSV of merged iCite and PubMed data for 
#                       PMIDs of interest
#   - reports/:
#       - Summary statistic report csvs with high-level grant data
#       - Report of publications removed from publication.tsv with reasons

import os
import sys
import pandas as pd
import config

from modules.data_preparation import load_and_clean_programs
from modules.nih_reporter_api import get_nih_reporter_grants
from modules.clean_grants_data import clean_grants_data
from modules.summary_statistics import get_summary_statistics
import modules.gather_publication_data as gpub

def main():


    # STEP 1: Data Prep - Load and clean Key Programs data
    print(f"Loading and processing {config.QUALTRICS_CSV_PATH}...")

    # Load and clean Key Programs
    continue_bool, key_programs_df = load_and_clean_programs(
                                    csv_filepath = config.QUALTRICS_CSV_PATH, 
                                    col_dict = config.QUALTRICS_COLS)

    # If issues are found, user is prompted to stop or continue
    if continue_bool == True:
        # Export cleaned Key Programs file
        key_programs_df.to_csv(config.CLEANED_KEY_PROGRAMS_CSV, index=False)
        print(f"Success! Saved {config.CLEANED_KEY_PROGRAMS_CSV}.")

    else:
        sys.exit("Process ended by user. Cleaned Key Programs csv not saved. "
            "Make manual edits to input CSV and retry.")



    # STEP 2: NIH RePORTER API - Get grants info for each Key Program
    print(f"---\nGathering, cleaning, and saving grants data from NIH "
          f"Reporter API for each Key Program...")

    # Create empty df to fill with grants
    all_cleaned_grants = pd.DataFrame()

    # Initialize a counter for total records
    total_records_count = 0

    # Iterate through each Key Program to get name, nofos, and awards
    for index, row in key_programs_df.iterrows():
        program_name = row['program_name']
        nofo_str = row['nofo']
        award_str = row['award']

        print(f"---\n{program_name}")
        
        # Convert semicolon-separated award string into listed values
        if isinstance(award_str, str): 
            award_list = [x.strip() for x in award_str.split(';')]
        else: award_list = []
        # Convert semicolon-separated NOFO string into listed values
        if isinstance(nofo_str, str): 
            nofo_list = [x.strip() for x in nofo_str.split(';')]
        else: nofo_list = []        

        # If neither NOFO nor Awards provided
        if pd.isna(nofo_str) and pd.isna(award_str):
            print(f"No NOFOs or Awards defined. Skipping.")
            continue

        # Gather grants data from NIH RePORTER
        award_grants_data = get_nih_reporter_grants(award_list, 'award', 
                                                    print_meta=True)
        nofo_grants_data = get_nih_reporter_grants(nofo_list, 'nofo', 
                                                   print_meta=True)
        
        # Combine and report
        grants_data = award_grants_data + nofo_grants_data
        total_records_count = total_records_count + len(grants_data)
        print(f"Complete. {len(grants_data)} grants found.")
        if award_grants_data and nofo_grants_data:
            print(f"Both Awards and NOFOs provided. Results may include "
                  f"duplicate grants gathered by both NOFO and Award "
                  f"searches.")
            continue

        # Skip remaining steps if no grants data gathered
        if len(grants_data) == 0:
            print(f"NOTICE: No grants found for {program_name}")
            continue

        # STEP 3: Data Cleaning - Clean and format grants data

        # Run cleaning steps
        cleaned_grants_data = clean_grants_data(grants_data)

        # Add program id column from alphanumeric program name
        program_id = ''.join(filter(str.isalnum, program_name))
        cleaned_grants_data[config.PROGRAM_ID_FIELDNAME] = program_id

    # STEP 4: Combine and save cleaned grants data into single tsv

        all_cleaned_grants = pd.concat([all_cleaned_grants,
                                        cleaned_grants_data])   


    # Define versioned output directory using config.py
    clean_grants_directory = config.PROCESSED_DIR
    if not os.path.exists(clean_grants_directory):
        os.makedirs(clean_grants_directory)
    project_filename = os.path.join(clean_grants_directory,
                                    config.PROJECTS_OUTPUT_FILENAME)

    # Sort by program and project for consistency
    all_cleaned_grants.sort_values(by=[config.PROGRAM_ID_FIELDNAME,
                                       config.PROJECT_ID_FIELDNAME], 
                                       inplace=True, ignore_index=True)

    # Export to csv
    all_cleaned_grants.to_csv(project_filename, sep = '\t', index=False)

    
    print(f"---\nSuccess! NIH RePORTER API data gathered, cleaned, and saved.\n"
          f"{total_records_count} grants gathered across all Awards and NOFOs.\n"
          f"{len(all_cleaned_grants)} NCI-funded grants retained for INS. \n"
          f"Results can be found in {clean_grants_directory}.\n---") 

        # TODO: Save all the printed console output to versioned txt 


    # STEP 5: Stats Report - Generate and output summary stats
    
    # Run summary statistic module that processes and exports reports as csv
    get_summary_statistics(all_cleaned_grants)
    print(f"Summary statistic reports for grants data successfully generated. \n"
          f"Results can be found in {config.API_REPORTS_DIR}.\n---")


    # STEP 6: Publications

    # Use existing PMID list if present
    if os.path.exists(config.PROJECT_PMIDS):
        print(f"---\nReusing PMIDs already gathered in {config.PROJECT_PMIDS}.")
        df_pmids = pd.read_csv(config.PROJECT_PMIDS)
    
    else:
        # Gather PMIDs from projects using NIH RePORTER API
        print(f"---\nGathering associated PMIDs for each project from "
              f"NIH RePORTER API...")
        df_pmids = gpub.get_pmids_from_projects(all_cleaned_grants)

        # Store a checkpoint file of all gathered PMIDs and projects
        df_pmids.to_csv(config.PROJECT_PMIDS, index=False)
        print(f"Gathered PMIDs saved to {config.PROJECT_PMIDS}.")

    # Gather publication data and store in chunks of defined size 
    print(f"---\nGathering PubMed data for all PMIDs...\n"
          f"---NOTE: THIS STEP CAN TAKE 8+ HOURS.---\n"
          f"Saving partial files to {config.TEMP_PUBLICATION_DIR} throughout.")
    gpub.build_pmid_info_data_chunks(df_pmids, 
                                chunk_size=config.PUB_DATA_CHUNK_SIZE, 
                                output_folder=config.TEMP_PUBLICATION_DIR)
    
    # Run again to catch any timeouts or problematic PMIDs
    print(f"First pass successful!\n"
          f"Double-checking for any missing publication info...")
    gpub.build_pmid_info_data_chunks(df_pmids, 
                            chunk_size=config.PUB_DATA_CHUNK_SIZE, 
                            output_folder=config.TEMP_PUBLICATION_DIR)
    print(f"Success! PubMed Publication data gathered for all PMIDs.")

    # Gather iCite data for all PMIDs from projects
    print(f"---\nGathering iCite data for all PMIDs...")

    # Use existing iCite-enriched PMIDs if present
    if os.path.exists(config.ICITE_PMID_DATA):
        print(f"---\nReusing iCite data already gathered in {config.ICITE_PMID_DATA}.")
        df_icite = pd.read_csv(config.ICITE_PMID_DATA)
    else:
        df_icite = gpub.get_icite_data_for_pmids(df_pmids,
                                config.ICITE_FILENAME, 
                                config.ICITE_COLUMNS_TO_PULL,
                                chunk_size=250000, chunk_count_est=146)
        # Export iCite data as checkpoint
        df_icite.to_csv(config.ICITE_PMID_DATA, index=False)
        print(f"---\niCite data for PMIDS saved to {config.ICITE_PMID_DATA}.")

    # Load all partial PubMed files back into a single df
    print(f"---\nCombining and merging all publication data...")
    df_pubmed = gpub.load_all_directory_files_to_df(config.TEMP_PUBLICATION_DIR)
    # Reset string-ed dates back to datetime
    df_pubmed['publication_date'] = pd.to_datetime(df_pubmed['publication_date'])

    # Fill in missing PubMed data and columns with iCite
    df_pub_info = gpub.merge_pubmed_icite_pmid_data(df_pubmed, df_icite)
    # Export for troubleshooting
    df_pub_info.to_csv(config.MERGED_PMID_DATA, index=False)
    print(f"Merged PMID data saved to {config.MERGED_PMID_DATA}.")

    # Build final publication df
    print(f"---\nMerging PMID data back to Core Projects...")
    df_publications, df_removed_publications = gpub.merge_and_clean_project_pmid_info(
                                                df_pmids, 
                                                df_pub_info)
    
    # Export final publication data
    df_publications.to_csv(config.PUBLICATIONS_OUTPUT, sep='\t', index=False)
    df_removed_publications.to_csv(config.REMOVED_PUBLICATIONS, index=False)

    print(f"Success! Publication data saved to {config.PUBLICATIONS_OUTPUT}.\n"
          f"Removed publications saved to {config.REMOVED_PUBLICATIONS}")
    print(f"---\n")
    print(f"Total unique Publications saved:  {df_publications['pmid'].nunique():>8}")
    print(f"Total removed Publications:       {len(df_removed_publications):>8}")
    print(f"Project-Publication associations: {len(df_publications):>8}")

if __name__ == "__main__":
    main()
