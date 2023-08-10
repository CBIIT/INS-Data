# config.py for INS Data Processing
# 2023-07-26 ZD
# Make changes here to affect variables throughout repo


# FILEPATH CONFIGURATION

# Edit version and type below for each new QualtricS file received
# Inputs and outputs will use this versioning
# Version must match suffix in input filename

QUALTRICS_VERSION = "2023-07-19" # <-- CHANGE VERSION HERE
QUALTRICS_TYPE = "manual_fix" # <-- Define "raw" or "manual_fix" type of the input csv




# --- DO NOT EDIT BELOW FOR ROUTINE DATA GATHERING ---

# FILEPATH (CONTINUED)
QUALTRICS_CSV_PATH = "data/raw/qualtrics_output_" + QUALTRICS_VERSION +"_"+ QUALTRICS_TYPE + ".csv"
CLEANED_KEY_PROGRAMS_CSV = "data/cleaned/key_programs_" + QUALTRICS_VERSION + ".csv"
PROCESSED_DIR = "data/processed/" + QUALTRICS_VERSION
REPORTS_DIR = "reports/" + QUALTRICS_VERSION

# ---
# DATA LOADING CONFIGURATION

# Dictionary of old:new column names to keep from Qualtrics
QUALTRICS_COLS = {
    "Name of Key Program": "program_name",
    "Acronym for key program": "program_acronym",
    "Focus Area (select all that apply)": "focus_area",
    "DOC": "doc",
    "Primary Contact (PI)": "contact_pi",
    "Primary Contact (PI) email": "contact_pi_email",
    "NIH Contact (Program Officer/Program Director)": "contact_nih",
    "NIH Contact (Program Officer/Program Director) email": "contact_nih_email",
    "NOFO number (eg. format as \"RFA-CA-00-000\") (If more than one, separate with ; semicolon)": "nofo",
    "Grant/Award number {parent award FORMAT #LL#CA######, eg. 5UG3CA260607} (If more than one, separate with ; semicolon)": "award",
    "Link to program website": "program_link",
    "Link to data or DCC if available": "data_link",
    "What type of cancer is the primary focus of the program? (Check all that\napply)": "cancer_type",
    "Login ID": "login_id"
}

# ---
# GRANTS CLEANING CONFIGURATION

# Earliest fiscal year (int) to use when gathering grants
API_EARLIEST_FISCAL_YEAR = 2000

# Grant/project fields to keep from NIH RePORTER results
API_FIELDS = [
    'project_num',
    'core_project_num',
    'appl_id',
    'fiscal_year',
    'project_title',
    'abstract_text',
    'pref_terms',
    'org_name',
    'org_city',
    'org_state',
    'org_country',
    'principal_investigators',
    'program_officers',
    'award_amount',
    'agency_ic_fundings',
    'award_notice_date',
    'project_start_date',
    'project_end_date',
    'full_foa',
    'api_source_search' # Not from RePORTER. Added during data gathering step.
]

# List of name fields received from the API that require reformatting
API_NAME_FIELDS = [
    'principal_investigators',
    'program_officers',
]

# Field name containing nested JSON with NCI funding
API_AGENCY_FUNDING_FIELD = 'agency_ic_fundings'

# Dictionary of old:new column names. Rename API fields to match INS terms
# Any terms not included will remain as retrieved from API
API_FIELD_RENAMER = {
    "project_num": "project_id",
    "core_project_num": "queried_project_id",
    "appl_id": "application_id",
    "pref_terms": "keywords",
    "agency_ic_fundings": "nci_funded_amount"
}