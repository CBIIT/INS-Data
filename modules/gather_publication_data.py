# gather_publication_data.py
# 2023-11-20 ZD
#
# This script defines functions that will accept the project.tsv and use it 
# along with the NIH RePORTER API, PubMed API, and iCite bulk download to 
# gather associated publications and descriptive data for those publications.
# 
# The output `publications_df` and publication.tsv contains columns:
# 
# coreproject
# pmid
# title
# authors
# publication_date
# citation_count
# relative_citation_ratio
#

import os
import sys
from datetime import datetime
from time import sleep  # for retrying API calls
import re

import pandas as pd
import requests
from tqdm import tqdm  # for progress bars
from math import ceil  # for pagination logging
from Bio import Entrez  # for PubMed API

# Append the project's root directory to the Python path
# This allows for importing config when running as part of main.py or alone
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import config


# This is similar but separate from the `get_nih_reporter_grants` function
# within the nih_reporter_api.py module
def get_pmids_from_nih_reporter_api(project_id,
                                    print_meta=False):
    """Get PMIDs associated with a single provided Project ID. 
    
    :param project_id: String Core Project ID (e.g. 'R01CA263500')
    :param print_meta: boolean indicator. If True, print API gathering 
                        process results to console.
    :return pmid_data: JSON API response with 'coreproject', 'pmid' and 'applId'
    """

    base_url = "https://api.reporter.nih.gov/v2/publications/search"
    pmid_data = []

    # Set default values for params not likley to change
    LIMIT = 500
    MAX_ATTEMPTS = 5
    RETRY_TIME = 2

    # Set starting value for counters
    offset = 0
    page = 0
    attempts = 0

    # RePORTER API sets a max limit of 500 records per call.
    # Keep looping each call in "pages" until the number of records 
    # gathered reaches the total number of records available. 

    # Set a cap on the number of attempts at a failed call before 
    # moving on to the next award.
    while attempts < MAX_ATTEMPTS:
        # Set parameters for API call
        params = {
            "criteria": {
                "core_project_nums": [project_id]
            },
            "offset": offset,
            "limit": LIMIT,
            "sort_field":"core_project_nums",
            "sort_order":"desc"
        }

        try: 
            # Define response details
            response = requests.post(base_url, 
                                    json=params, 
                                    headers={
                                        "accept": "application/json", 
                                        "Content-Type": "application/json"})

            # If response is good, get results
            if response.status_code == 200:
                pmids = response.json()

                # Add grants to running list
                pmid_data.extend(pmids['results'])

                # Increase offset by limit to get next "page"
                total_records = pmids['meta']['total']
                offset = offset + LIMIT
                page = page + 1

                # Print paginated partial optional metadata
                # Consider replacing this with proper logging
                if print_meta == True:
                    total_pages = max(ceil(total_records/LIMIT),1)
                    print(f"{project_id}: "
                          f"({page}/{total_pages}): {pmids['meta']}")

                # Stop looping if offset has reached total record count
                if offset >= total_records:
                    break

            # Handle 500 errors by retrying after 2 second delay
            elif response.status_code == 500:
                attempts = attempts + 1
                print(f"Received a 500 error for "
                        f"{project_id}'. "
                        f"Retrying after {RETRY_TIME} seconds. "
                        f"Attempt {attempts}/{MAX_ATTEMPTS}")
                sleep(RETRY_TIME)
            else:
                print(f"Error occurred while fetching pmids for "
                        f"{project_id}': "
                        f"{response.status_code}")
                break

        except requests.exceptions.RequestException as e:
            print(f"An error occurred while making the API call for "
                    f"{project_id}': {e}")
            break

    return pmid_data


def get_pmids_from_projects(projects_df, print_meta=False):
    """Iterate through project data and use NIH RePORTER API to get all PMIDs.

    :param pd.Dataframe projects_df: DataFrame of all cleaned projects from the
                grants workflow
    :param bool print_meta: boolean indicator. If True, print status to console
    """ 

    # Get all unique project IDs from projects.tsv (loaded earlier in notebook)
    project_id_list = projects_df['queried_project_id'].unique().tolist()

    print(f"{len(project_id_list)} total unique project IDs.")

    # Create empty df
    all_pmids = []

    # Iterate through projects and add to running list of results
    for id in tqdm(project_id_list, ncols=80):
        results = get_pmids_from_nih_reporter_api(id, print_meta=print_meta)
        all_pmids.extend(results)

    # Reformat list of results as dataframe
    df_pmids = pd.DataFrame(all_pmids)

    # Drop irrelevant applid column and any duplicated rows
    df_pmids = df_pmids.drop('applid', axis=1).drop_duplicates()

    # Sort by PMIDs for consistent downstream handling
    df_pmids = df_pmids.sort_values(by='pmid')

    total_pubs = len(df_pmids)
    unique_pmids = df_pmids['pmid'].nunique()
    unique_pubs = df_pmids['coreproject'].nunique()

    # Print results
    print(f"---\nComplete! PMIDs successfully gathered. \n"
        f"{total_pubs} total Publication results. \n"
        f"{unique_pmids} unique PMIDs found across \n"
        f"{unique_pubs} unique Project IDs.")

    return df_pmids



def get_icite_data_for_pmids(df_pmid, icite_filepath, cols, 
                             chunk_size=250000, chunk_count_est=None):
    """
    Get iCite data for unique PMIDs from a DataFrame using chunks.

    :param df_pmid (pd.DataFrame): DataFrame containing PMID-related data.
    :param icite_filepath (str): Path to the zipped iCite CSV file.
    :param cols (list): Columns to pull from iCite and include in output df
    :param chunk_size (int): Number of rows to read per chunk.
    :param chunk_count (int): Estimated number of chunks. If None, ignored

    :return: pd.DataFrame: DataFrame with specified columns
    """
    # Get unique PMIDs from df_pmid
    unique_pmids = df_pmid['pmid'].unique()

    # Initialize an empty list to store the enriched data chunks
    df_enriched_chunks = []

    # Create a tqdm wrapper around the generator to track progress
    chunks = tqdm(pd.read_csv(icite_filepath, compression='zip', 
                              chunksize=chunk_size),
                  unit="chunk", total=chunk_count_est, ncols=80)

    # Iterate through chunks of the iCite DataFrame
    for chunk in chunks:
        # Filter the chunk to include only rows with PMIDs in unique_pmids
        chunk_filtered = chunk[chunk['pmid'].isin(unique_pmids)]

        # Append the filtered DataFrame to the list with selected columns
        df_enriched_chunks.append(chunk_filtered[cols])

    # Concatenate all the chunks into the final enriched DataFrame
    df_enriched = pd.concat(df_enriched_chunks, ignore_index=True)

    # Include unique PMIDs not found in iCite with NaN in iCite columns
    df_unique_pmids = pd.DataFrame(df_pmid['pmid'].unique(), columns=['pmid'])
    df_pmids_with_icite = pd.merge(df_unique_pmids, df_enriched, 
                                   how='left', on='pmid')

    # Sort by PMID for consistent handling
    df_pmids_with_icite.sort_values(by='pmid', inplace=True, ignore_index=True)

    return df_pmids_with_icite



def get_full_publication_record(pmid):
    """
    Get the full record for a given PMID without subselection. Helper function 
    useful for browsing the schema but not used in production workflow.

    :param pmid: PubMed ID (str)
    :return: Full record (dictionary)
    """
    # Get user email from hidden local env file. Use default if not defined
    Entrez.email = os.environ.get('NCBI_EMAIL', 'your-email@example.com')
    Entrez.api_key = os.environ.get('NCBI_API_KEY', '')

    try:
        handle = Entrez.efetch(db="pubmed", id=pmid, retmode="xml")
        records = Entrez.read(handle)
        handle.close()

        # Return the full record
        return records

    except Exception as e:
        print(f"Error fetching information for PMID {pmid}: {e}")
        return None



def format_authors(author_list):
    """
    Format author names as 'FirstName LastName'.

    :param author_list: List of authors
    :return: Formatted author names
    """
    formatted_authors = []
    for author in author_list:
        last_name = author.get('LastName', '')
        first_name = author.get('ForeName', '')
        formatted_author = f"{first_name} {last_name}".strip()
        formatted_authors.append(formatted_author)

    return ', '.join(formatted_authors)



def extract_medline_date_components(date_str):
    """Extracts date components from a Medline date string with various formats

    :param date_str: A string containing date information in the Medline format
    :return: Publication date as datetime or None if an error occurs
    """

    try:
        # Look for the first 4-digit number and store as the year
        year_match = re.search(r'\b\d{4}\b', date_str)
        if year_match:
            year = int(year_match.group())
        else:
            raise ValueError("No year found")

        # Find first month occurrence (Mmm) and store as the numerical month
        month_match = re.search(
                    r'\b(?:Jan|Feb|Mar|Apr|May|Jun|Jul|Aug|Sep|Oct|Nov|Dec)\b',
                    date_str, re.I)
        if month_match:
            month = datetime.strptime(month_match.group(), "%b").month
        else:
            # Look for season string (Winter, Spring, Summer, Fall) and store 
            # as the numerical month (12, 3, 6, 9).
            season_match = re.search(r'\b(?:Winter|Spring|Summer|Fall)\b',
                                      date_str, re.I)
            if season_match:
                season_to_month = {'Winter': 12, 'Spring': 3, 
                                   'Summer': 6, 'Fall': 9}
                month = season_to_month[season_match.group()]
            else:
                month = 1  # If no month or season, assign 1.

        # Look for first occurrence of 1 or 2 digits after space and store as day
        day_match_space = re.search(r'\b\d{1,2}\b', date_str)

        if day_match_space:
            day = int(day_match_space.group())
        else:
            day = 1  # If not identified, assign 1.

        return datetime(year, month, day) #.strftime('%Y-%m-%d')

    except Exception as e:
        print(f"Error processing date string '{date_str}': {e}")
        return None



def format_publication_date(pub_date):
    """Format publication date API response as a standardized datetime. 
    Replace any missing months or dates with 1. (i.e. January or the 1st)

    :param dict pub_date: JSON record of publication date from the PubMed API
    :return datetime: Publication date as datetime
    """
    # Return None if pub_date is empty
    if not pub_date:
        return None

    # Get publication year (default to None if not present)
    year = pub_date.get('Year', None)

    # If year is present, parse expected format into datetime
    if year:
        # Get month. If not present, default to January (01)
        month = pub_date.get('Month', 1)

        # If month is provided as a string, convert to numerical representation
        if isinstance(month, str):
            month_dict = {'Jan': 1, 'Feb': 2, 'Mar': 3, 'Apr': 4, 
                          'May': 5, 'Jun': 6, 'Jul': 7, 'Aug': 8, 
                          'Sep': 9, 'Oct': 10, 'Nov': 11, 'Dec': 12}
            month = month_dict.get(month, 1)

        # Get day. If not present, efault to 1st (01)
        day = pub_date.get('Day', 1)

        # Combine components into datetime
        date_out = datetime(int(year), int(month), int(day))

    # If year is not present, use Medline date parsing
    else: 
        date_out = extract_medline_date_components(pub_date['MedlineDate'])

    return date_out



def get_pubmed_info_from_pmid(pmid):
    """
    Get publication information for a given PMID.

    :param pmid: PubMed ID (str)
    :return: Dictionary containing publication information
    """
    # Get user email from hidden local env file. Use default if not defined
    Entrez.email = os.environ.get('NCBI_EMAIL', 'your-email@example.com')
    Entrez.api_key = os.environ.get('NCBI_API_KEY', '')
    # Reduce delay between failed attempts from default 15
    Entrez.sleep_between_tries = 2

    try:
        handle = Entrez.efetch(db="pubmed", id=pmid, retmode="xml")
        records = Entrez.read(handle)
        handle.close()

        # Access article details from full returned record
        record = records['PubmedArticle'][0]['MedlineCitation']['Article']

        # Isolate and format author names
        authors = record.get('AuthorList', [])
        formatted_authors = format_authors(authors)

        # Get publication date information
        # Use 'ArticleDate' if available
        article_date = record.get('ArticleDate', [])
        if article_date:
            publication_date = article_date[0]
        # Try using Journal PubDate field next
        elif record['Journal']['JournalIssue'].get('PubDate', {}):
            publication_date = record['Journal']['JournalIssue'].get('PubDate', {})
        # Use Medline date if no others avaialble
        else:
            publication_date = (record['Journal']['JournalIssue']
                                        ['PubDate'].get('MedlineDate', ''))

        # Use the format_publication_date function
        publication_date = format_publication_date(publication_date)

        publication_info = {
            'publication_id': pmid,
            'title': record['ArticleTitle'],
            'authors': formatted_authors,
            'publication_date': publication_date,
        }

        return publication_info

    except Exception as e:
        # Use tqdm.write() instead of print() for long processes
        tqdm.write(f"Error fetching information for PMID {pmid}: {e}")
        return None



def get_max_chunk_number(output_folder):
    """Get the maximum chunk number from existing filenames."""
    
    existing_chunk_files = [file for file in os.listdir(output_folder) 
                            if file.startswith('publicationDetails')]
    if not existing_chunk_files:
        return 0
    else:
        # Extract numbers from filenames and return the maximum
        return max([int(re.search(r'\d+', file)
                        .group()) for file in existing_chunk_files])



def load_all_directory_files_to_df(directory):
    """Load all identically-structured files within a folder into a single df."""
    
    filenames = [file for file in os.listdir(directory)]
    filepaths = [directory + '/' + file for file in filenames]
    df = pd.concat(map(pd.read_csv, filepaths), ignore_index=True)
    
    return df



def build_pmid_info_data_chunks(df_pmid, output_folder, chunk_size):
    """
    Get publication information for each PMID in the input DataFrame and export
    as separate csv files. 

    :param df_pmid: DataFrame containing 'coreproject' and 'pmid' columns.
    :param output_folder: Folder to store the output files.
    :param chunk_size: Number of records to process in each batch.
    :return: None. Files are stored for loading
    """
    # Create the output folder if it doesn't exist
    os.makedirs(output_folder, exist_ok=True)
    
    # Create empty df for iteration
    df_pmid_info = pd.DataFrame()
    
    # Check for pre-existing file chunks with publication details
    chunk_number = get_max_chunk_number(output_folder)
        
    # Track previously processed pmids
    processed_pmids = set()

    # If previously gathered publication data exist, reuse and reset counts
    if chunk_number > 0: 
        existing_data = load_all_directory_files_to_df(output_folder)
        processed_pmids.update(existing_data['pmid'].unique())
        print(f"---\nReusing existing Publication details for "
              f"{len(processed_pmids)} PMIDs found in {output_folder}.\n---")

    # Iterate through each unique PMID with tqdm progress bar
    for pmid in tqdm(df_pmid['pmid'].unique(), ncols=80):
        # Check if PMID is already processed and skip if so
        if pmid in processed_pmids:
            continue

        try:
            # Use PubMed API to get publication data
            publication_info = get_pubmed_info_from_pmid(pmid)

            if publication_info:
                # Combine the information with the original DataFrame
                df_current = pd.DataFrame({
                    'pmid': pmid,
                    'title': publication_info['title'],
                    'authors': publication_info['authors'],
                    'publication_date': publication_info['publication_date']
                }, index=[0])

                # Add the current DataFrame to df_pmid_info
                df_pmid_info = pd.concat([df_pmid_info, df_current], ignore_index=True)

                # Check if df should be saved to file
                if len(df_pmid_info) >= chunk_size:
                    # Save with standard filename and numbered suffix
                    chunk_number += 1
                    output_file = f"{output_folder}/publicationDetails_{chunk_number:03}.csv"
                    df_pmid_info.to_csv(output_file, index=False)

                    # Clear df_pmid_info for the next chunk
                    df_pmid_info = pd.DataFrame()

                    # Add processed PMID to the set
                    processed_pmids.add(pmid)

        except Exception as e:
            print(f"Error processing PMID {pmid}: {e}")
            # Fill in fields with NaN if not available
            df_current = pd.DataFrame({
                'pmid': pmid,
                'title': pd.NA,
                'authors': pd.NA,
                'publication_date': pd.NA
            }, index=[0])

            # Add the current DataFrame to df_pmid_info
            df_pmid_info = pd.concat([df_pmid_info, df_current], ignore_index=True)

            # Add processed PMID to the set
            processed_pmids.add(pmid)

    # Save the final chunk if present
    if len(df_pmid_info) >= 1:
        # Save with standard filename and numbered suffix
        chunk_number += 1
        output_file = f"{output_folder}/publicationDetails_{chunk_number:03}.csv"
        df_pmid_info.to_csv(output_file, index=False)

    return None



def merge_pubmed_icite_pmid_data(df_pubmed, df_icite):
    """Fill in missing PubMed data for PMIDs with iCite data for PMIDs.

    :param pd.DataFrame df_pubmed: pandas DataFrame of PubMed data with fields 
                    'pmid', 'title', 'authors', 'publication_date'
    :param pd.DataFrame df_icite: pandas DataFrame of iCite data with fields 
                    'pmid', 'title', 'authors', 'year', 
                    'citation_count', 'relative_citation_ratio'
    :return pd.DataFrame df_pub_info: pandas DataFrame of combined data with fields
                    'pmid', 'title', 'authors', 'publication_date', 
                    'citation_count', 'relative_citation_ratio'
    """

    # Format iCite year column as datetime
    df_icite['publication_year'] = pd.to_datetime(df_icite['year'])
    df_icite.drop(columns='year', inplace=True)

    # Fill in missing PubMed data and columns with iCite
    df_pub_info = df_pubmed.combine_first(df_icite)

    # Reorder columns and sort for consistency
    df_pub_info = df_pub_info[['pmid','title','authors','publication_date',
                               'citation_count','relative_citation_ratio']]
    df_pub_info.sort_values(by='pmid', inplace=True, ignore_index=True)

    return df_pub_info



def merge_and_clean_project_pmid_info(df_pmids, df_pub_info):
    """Merge publication information with the dataframe of projects and pmids. 
       Also perform some cleaning functions and store removed publications.

    :param df_pmid: DataFrame containing 'coreproject' and 'pmid' columns.
    :param df_pub_info: DataFrame containing 'pmid', 'title', 'authors', and 
                    'publication_date'.
    :return: df_merged: Clean DataFrame with projects, pmids, and pub info
    :return: df_removed_publications: DataFrame with errored publication info
    """
    # Merge the two DataFrames on the 'pmid' column
    df_merged = pd.merge(df_pmids, df_pub_info, on='pmid', how='left')

    # Remove rows where all values are NaN except for 'pmid'
    df_removed_publications = df_merged.loc[df_merged.drop('pmid', axis=1)
                                            .isnull().all(axis=1)].copy()
    df_merged = df_merged.dropna(subset=df_merged.columns.drop('pmid'), how='all')

    # Remove rows where 'publication_date' is below the cutoff in config
    cutoff_year = config.PUBLICATION_YEAR_CUTOFF
    df_removed_early_pubdate = df_merged[df_merged['publication_date']
                                         < datetime(cutoff_year, 1,1)].copy()
    df_merged = df_merged[df_merged['publication_date'] 
                                        >= datetime(cutoff_year, 1, 1)]

    # Add reasons for removal to 'df_removed_publications'
    df_removed_publications = pd.concat([df_removed_publications, 
                                         df_removed_early_pubdate], 
                                         ignore_index=True)
    df_removed_publications['reason'] = ''

    # Add reasons for removal based on conditions
    df_removed_publications.loc[df_removed_publications.drop('pmid', axis=1)
                                .isnull().all(axis=1), 'reason',] = 'No publication info'
    df_removed_publications.loc[df_removed_publications['publication_date']
                                < datetime(cutoff_year, 1, 1), 'reason',] = 'Published before 2000'

    return df_merged, df_removed_publications



# Run module as a standalone script when called directly
if __name__ == "__main__":

    print(f"Running {os.path.basename(__file__)} as standalone module...")

    # Load Grants data
    project_filepath = os.path.join(config.PROCESSED_DIR, 
                               config.PROJECTS_OUTPUT_FILENAME)
    all_cleaned_grants = pd.read_csv(project_filepath, sep='\t')
    print(f"Projects data loaded from {project_filepath}.")

    # Use existing PMID list if present
    if os.path.exists(config.PROJECT_PMIDS):
        print(f"---\nReusing PMIDs already gathered in {config.PROJECT_PMIDS}.")
        df_pmids = pd.read_csv(config.PROJECT_PMIDS)
    
    else:
        # Gather PMIDs from projects using NIH RePORTER API
        print(f"---\nGathering associated PMIDs for each project from "
              f"NIH RePORTER API...")
        df_pmids = get_pmids_from_projects(all_cleaned_grants)

        # Store a checkpoint file of all gathered PMIDs and projects
        df_pmids.to_csv(config.PROJECT_PMIDS, index=False)
        print(f"Gathered PMIDs saved to {config.PROJECT_PMIDS}.")

    # Gather publication data and store in chunks of defined size 
    print(f"---\nGathering PubMed data for all PMIDs...\n"
          f"---NOTE: THIS STEP CAN TAKE 8+ HOURS.---\n"
          f"Saving partial files to {config.TEMP_PUBLICATION_DIR} throughout.")
    build_pmid_info_data_chunks(df_pmids, 
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
        df_icite = get_icite_data_for_pmids(df_pmids,
                                config.ICITE_FILENAME, 
                                config.ICITE_COLUMNS_TO_PULL,
                                chunk_size=250000, chunk_count_est=146)
        # Export iCite data as checkpoint
        df_icite.to_csv(config.ICITE_PMID_DATA, index=False)
        print(f"---\niCite data for PMIDS saved to {config.ICITE_PMID_DATA}.")

    # Load all partial PubMed files back into a single df
    print(f"---\nCombining and merging all publication data...")
    df_pubmed = load_all_directory_files_to_df(config.TEMP_PUBLICATION_DIR)
    # Reset string-ed dates back to datetime
    df_pubmed['publication_date'] = pd.to_datetime(df_pubmed['publication_date'])

    # Fill in missing PubMed data and columns with iCite
    df_pub_info = merge_pubmed_icite_pmid_data(df_pubmed, df_icite)
    # Export for troubleshooting
    df_pub_info.to_csv(config.MERGED_PMID_DATA, index=False)
    print(f"Merged PMID data saved to {config.MERGED_PMID_DATA}.")

    # Build final publication df
    print(f"---\nMerging PMID data back to Core Projects...")
    df_publications, df_removed_publications = merge_and_clean_project_pmid_info(
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