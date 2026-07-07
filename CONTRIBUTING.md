# Contributing to INS-Data

Thank you for your interest in contributing to the INS-Data repository! This guide covers the essentials for getting started. If you have additional questions or interest, please reach out to the NCI Office of Data Sharing at [NCIOfficeofDataSharing@mail.nih.gov](NCIOfficeofDataSharing@mail.nih.gov).

## Setup

```bash
git clone https://github.com/CBIIT/INS-Data.git
cd INS-Data
uv venv .venv
.venv\Scripts\activate  # Windows. On macOS/Linux: source .venv/bin/activate
uv pip install -r requirements.txt
```

You will also need a `.env` file in the project root with your NCBI credentials:

```bash
NCBI_EMAIL=your-email@example.com
NCBI_API_KEY=your-api-key-here
```

Get an API key from [NCBI](https://ncbiinsights.ncbi.nlm.nih.gov/2017/11/02/new-api-keys-for-the-e-utilities/). See the [README](README.md#how-to-use-this-repository) for full prerequisites including the Qualtrics CSV and iCite database download.

## Running Tests

Run the full test suite with coverage:

```bash
pytest --cov=modules -v
```

The suite includes live API smoke tests that call external services (NIH RePORTER, NCBI). These run by default but can be skipped when offline:

```bash
pytest -m "not live_api"
```

All tests should pass before submitting a PR.

## Branch Naming

Use descriptive branch names that reference the relevant ticket or feature:

- `ins-1234-add-new-dataset-source`
- `fix-publication-date-parsing`
- `update-dbgap-curation-2026`

## Pull Requests

- Branch from `dev`
- Keep PRs focused on a single change or related set of changes
- Fill out the PR template with a summary, list of changes, and testing steps
- Ensure `pytest --cov=modules -v` passes before requesting review

## Project Structure

```bash
main.py              # Pipeline orchestrator
config.py            # All configuration, paths, and field mappings
modules/             # One module per pipeline step
tests/               # One test file per module (test_<module_name>.py)
data/00_input/       # Raw inputs (Qualtrics CSV, iCite, dbGaP, CEDCD)
data/01_intermediate/# Intermediate CSVs produced during the pipeline
data/02_output/      # Final TSVs ready for INS ingestion
reports/             # Validation reports and statistics
```

See the [README](README.md#data-gathering-workflow) for detailed documentation of each pipeline step.
