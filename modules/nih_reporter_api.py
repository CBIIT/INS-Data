# nih_reporter_api.py
# 2023-08-04 ZD
# 
# This script defines the primary function get_nih_reporter_grants that calls
# the public NIH RePORTER API to gather grants data for the provided NOFO 
# (Notice of Funding Opportunity) and/or awards (i.e. grants, supplements, or
# parent projects) for each Key Program in the cleaned Key Programs csv. 
#
# NIH RePORTER Endpoint: https://api.reporter.nih.gov/v2/projects/search 
# NIH RePORTER Docs: https://api.reporter.nih.gov/
#
# NOTE: From a banner at the top of https://api.reporter.nih.gov/ on 8/4/23: 
#   Effective August 31, 2023, there will be a change in the way FOA Numbers 
#   are handled through the API Service. Moving forward, only the long FOA 
#   Number (Full_FOA) will be available, and it will be referred to as the 
#   "Opportunity Number." Additionally, the short FOA Number will no longer be 
#   supported and will be removed from the API."


import requests
import config
from datetime import datetime
from time import sleep # for retrying API calls
from math import ceil # for pagination logging

def get_nih_reporter_grants(search_values:str, 
                            search_type:str, 
                            print_meta=False):
    """Get grants info for either NOFOs or Awards as specified.
    
    :param search_values: string of values to search, separated by semicolons

    :param search_type: type of search (e.g., 'nofo' or 'award')

    :param print_meta: boolean indicator. If True, print API gathering 
                        process results to console.
    """

    # Set the search field based on the search_type
    if search_type == 'nofo':
        search_field = "opportunity_numbers"
    elif search_type == 'award':
        search_field = "project_nums"
    else:
        raise ValueError("Invalid search type.")

    base_url = "https://api.reporter.nih.gov/v2/projects/search"
    grants_data = []

    for search_value in search_values:
        # Check for blank search values
        if not search_value:
            print(f"Blank {search_type} value encountered. Skipping.")
            continue

        # Set default values for params not likley to change
        LIMIT = 500
        MAX_ATTEMPTS = 5
        RETRY_TIME = 2
        EARLIEST_FISCAL_YEAR = config.API_EARLIEST_FISCAL_YEAR

        # Get list of fiscal years to query
        current_year = datetime.now().year
        fiscal_year_list = [str(year) for year in range(
            EARLIEST_FISCAL_YEAR, current_year + 1)]

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
                    search_field: [search_value],
                    "exclude_subprojects": True,
                    "agencies": ["NCI"],
                    "is_agency_funding": True,
                    "fiscal_years": fiscal_year_list
                },
                "limit": LIMIT,
                "offset": offset,
                "sort_field": "FiscalYear",
                "sort_order": "desc",
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
                    grants = response.json()
                    # Add API source indicator
                    for grant in grants['results']:
                        grant['api_source_search'] = f"{search_type}_{search_value}"
                    # Add grants to running list
                    grants_data.extend(grants['results'])

                    # Increase offset by limit to get next "page"
                    total_records = grants['meta']['total']
                    offset = offset + LIMIT
                    page = page + 1

                    # Print paginated partial optional metadata
                    # Consider replacing this with proper logging
                    if print_meta == True:
                        total_pages = max(ceil(total_records/LIMIT),1)
                        print(f"{search_type}: {search_value} "
                              f"({page}/{total_pages}): {grants['meta']}")

                    # Stop looping if offset has reached total record count
                    if offset >= total_records:
                        break

                # Handle 500 errors by retrying after 2 second delay
                elif response.status_code == 500:
                    attempts = attempts + 1
                    print(f"Received a 500 error for "
                          f"{search_type} '{search_value}'. "
                          f"Retrying after {RETRY_TIME} seconds. "
                          f"Attempt {attempts}/{MAX_ATTEMPTS}")
                    sleep(RETRY_TIME)
                else:
                    print(f"Error occurred while fetching grants for "
                          f"{search_type} '{search_value}': "
                          f"{response.status_code}")
                    break

            except requests.exceptions.RequestException as e:
                print(f"An error occurred while making the API call for "
                      f"{search_type} '{search_value}': {e}")
                break

    return grants_data