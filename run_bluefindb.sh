#!/bin/bash

# =======================================================
# BlueFinDB v2.0 - Publication-Ready Reference DB Builder
# User Interface Script (Wrapper)
# =======================================================

# --- Configuration & Defaults ---
DEFAULT_SLEEP=15
DB_NAME="BlueFinDB_12S_Fish"
MIFISH_U_F="GTCGGTAAAACTCGTGCCAGC"
MIFISH_U_R="CATAGTGGGGTATCTAATCCCAGTTTG"
MIFISH_E_F="GTCGGTAAAACTCGTGCCAGC"
MIFISH_E_R="CAGGTAGTGGGGTATCTAATCCCAGTTTG"

# Variables for argument passing
SPECIES_LIST=""
EMAIL=""
API_KEY=""
SLEEP_TIME=$DEFAULT_SLEEP
F_PRIMER=""
R_PRIMER=""
ALT_F_PRIMER=""
ALT_R_PRIMER=""

# --- 1. Argument Parsing ---
while [[ $# -gt 0 ]]; do
    key="$1"
    case $key in
        --email)
        EMAIL="$2"
        shift # past argument
        shift # past value
        ;;
        --api-key)
        API_KEY="$2"
        shift # past argument
        shift # past value
        ;;
        --sleep)
        SLEEP_TIME="$2"
        shift # past argument
        shift # past value
        ;;
        --trim-primers)
        F_PRIMER="$2"
        R_PRIMER="$3"
        shift 
        shift 
        shift
        ;;
        --mifish-preset)
        F_PRIMER="$MIFISH_U_F"
        R_PRIMER="$MIFISH_U_R"
        ALT_F_PRIMER="$MIFISH_E_F"
        ALT_R_PRIMER="$MIFISH_E_R"
        shift
        ;;
        *)
        if [[ -z "$SPECIES_LIST" ]]; then
            SPECIES_LIST="$1"
        fi
        shift
        ;;
    esac
done

# --- 2. Initial Checks and Validation ---
echo "--- BlueFinDB v2.0: Starting Setup ---"

# Check mandatory inputs
if [ -z "$SPECIES_LIST" ]; then
    echo "FATAL ERROR: Missing required species list file."
    echo "Usage: ./run_bluefindb.sh <species_list.txt> --email <your@email.com> [--mifish-preset]"
    exit 1
fi
if [ ! -f "$SPECIES_LIST" ]; then
    echo "FATAL ERROR: Species list file not found: $SPECIES_LIST"
    exit 1
fi
if [ -z "$EMAIL" ]; then
    echo "FATAL ERROR: NCBI requires an email address. Please provide using --email user@domain.com"
    exit 1
fi

# Check required commands (Dependencies)
missing=()
for cmd in python3 esearch efetch makeblastdb; do
    if ! command -v "$cmd" >/dev/null 2>&1; then
        missing+=("$cmd")
    fi
done

# Check cutadapt only if trimming is requested
if [ -n "$F_PRIMER" ]; then
    if ! command -v cutadapt >/dev/null 2>&1; then
        missing+=("cutadapt")
    fi
    echo "NOTICE: Trimming requested. Checking for cutadapt... OK."
fi

if [ "${#missing[@]}" -ne 0 ]; then
    echo "FATAL ERROR: Missing required commands: ${missing[*]}"
    echo "Install NCBI EDirect, python3 (with Biopython), BLAST+, and possibly cutadapt."
    exit 2
fi

# --- 3. Environment Setup and Info ---
export ENTREZ_EMAIL="$EMAIL"
if [[ -n "$API_KEY" ]]; then
    export NCBI_API_KEY="$API_KEY"
    echo "NOTICE: Using provided API Key for higher NCBI limits."
fi

# Create directories
RAW_DIR="fish_12s_raw"
FILTERED_DIR="fish_12s_filtered"
mkdir -p "$RAW_DIR" "$FILTERED_DIR"

# Estimate Time
SPECIES_COUNT=$(wc -l < "$SPECIES_LIST")
# Remove comment lines and empty lines from count
SPECIES_COUNT=$(grep -vE '^\s*#|^\s*$' "$SPECIES_LIST" | wc -l)

EST_SECONDS=$((SPECIES_COUNT * SLEEP_TIME))
EST_MINUTES=$((EST_SECONDS / 60))

echo "Species to process: $SPECIES_COUNT"
echo "Sleep time: ${SLEEP_TIME}s. Estimated Run Time: ~${EST_MINUTES} minutes."
echo "Output DB Name: $DB_NAME"
echo "-------------------------------------"

# --- 4. Launch Core Python Script (Stages 2, 3, 4, 5) ---
# Build the Python arguments dynamically
PYTHON_ARGS=(
    --list "$SPECIES_LIST"
    --sleep "$SLEEP_TIME"
)

if [[ -n "$F_PRIMER" ]]; then
    PYTHON_ARGS+=( --f-primer "$F_PRIMER" --r-primer "$R_PRIMER" )
fi

if [[ -n "$ALT_F_PRIMER" ]]; then
    PYTHON_ARGS+=( --alt-f-primer "$ALT_F_PRIMER" --alt-r-primer "$ALT_R_PRIMER" )
fi

# Execute the core Python script
python3 bluefindb_core.py "${PYTHON_ARGS[@]}"

# Check if Python script failed
if [ $? -ne 0 ]; then
    echo "FATAL ERROR: Core Python script failed. Check messages above."
    exit 3
fi

# --- 5. Automatic BLAST Indexing (Stage 6) ---
echo -e "\n[STEP 6/6] Automatically indexing final FASTA using makeblastdb..."
FINAL_FASTA="fish_12s_BLAST.fasta"

if [ ! -f "$FINAL_FASTA" ]; then
    echo "FATAL ERROR: Final FASTA file ($FINAL_FASTA) not found. Indexing aborted."
    exit 4
fi

makeblastdb -in "$FINAL_FASTA" -dbtype nucl -out "$DB_NAME" 

if [ -f "${DB_NAME}.nsq" ]; then 
    echo -e " SUCCESS! BLAST database is fully indexed and ready to use."
    echo -e "\n--- BlueFinDB COMPLETE ---"
    echo -e "Your searchable database is ready. Use the name \033[1m$DB_NAME\033[0m for your blastn queries."
    echo "Example: blastn -query my_otus.fasta -db $DB_NAME -out results.tsv -outfmt 6"
else
    echo "FATAL ERROR: makeblastdb failed to create the database index files."
    exit 5
fi