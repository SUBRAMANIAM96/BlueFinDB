#!/bin/bash

# ==============================================================================
# BlueFinDB v3.0 - Universal Reference Database Builder
# [PUBLICATION RELEASE]
#
# Author: SUBRAMANIAM VIJAYAKUMAR
#
# DESCRIPTION:
# Orchestrates the curation pipeline:
# 1. Sets marker-specific parameters (12S, 16S, 18S, COI).
# 2. Validates environment and dependencies.
# 3. Executes Core Engine (Data Retrieval & Trimming).
# 4. Triggers Rescue System for missing taxa.
# 5. Performs Final Quality Control (QC).
# ==============================================================================

# Exit on error
set -e

# --- VISUAL CONFIGURATION ---
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[1;34m'
CYAN='\033[0;36m'
BOLD='\033[1m'
NC='\033[0m' # No Color

# --- CONFIGURATION DEFAULTS ---
VERSION="v3.0"
DEFAULT_SLEEP=1 
DB_NAME_BASE="BlueFinDB"
DEFAULT_INPUT_FILE="Data/species_list_16S.txt"

# ==============================================================================
# 1. PRIMER & LENGTH DEFINITIONS
# ==============================================================================

# --- 12S (Ichthyans) ---
# MiFish-U and MiFish-E
P1_12S_F="GTCGGTAAAACTCGTGCCAGC"
P1_12S_R="CATAGTGGGGTATCTAATCCCAGTTTG"
P2_12S_F="GTCGGTAAAACTCGTGCCAGC"
P2_12S_R="CATAGTGGGGTATCTAATCCTAGTTTG"
LEN_12S_MIN=150
LEN_12S_MAX=550

# --- COI (Invertebrates) ---
# Leray-XT
P1_COI_F="GGWACWGGWTGAACWGTWTAYCCYCC"
P1_COI_R="TAIACYTCIGGRTGICCRAARAAYCA"
LEN_COI_MIN=100
LEN_COI_MAX=800

# --- 16S & 18S (UNIVERSAL SUITES) ---
# Primers are handled internally by the Python core for multi-primer sweeping.
# Length constraints are defined here.

# 16S: Wide range to accommodate short Mammal (140bp) and long Vertebrate (600bp) amplicons
LEN_16S_MIN=100
LEN_16S_MAX=600

# 18S: Wide range for Fungi (variable) and Nematodes
LEN_18S_MIN=100
LEN_18S_MAX=1000

# ==============================================================================
# 2. ARGUMENT PARSING
# ==============================================================================
MARKER="12S" # Default
NON_INTERACTIVE=false
SPECIES_LIST=""

while [[ "$#" -gt 0 ]]; do
    case $1 in
        --COI) MARKER="COI"; shift ;;
        --16S) MARKER="16S"; shift ;;
        --18S) MARKER="18S"; shift ;;
        --12S) MARKER="12S"; shift ;;
        --non-interactive) NON_INTERACTIVE=true; shift ;;
        *) 
            if [ -z "$SPECIES_LIST" ]; then
                SPECIES_LIST="$1"
            fi
            shift 
            ;;
    esac
done

if [ -z "$SPECIES_LIST" ]; then
    SPECIES_LIST="$DEFAULT_INPUT_FILE"
fi

# --- DISPLAY HEADER ---
clear
echo -e "${BLUE}"
echo "  ____  _            _____ _       _____  ____  "
echo " |  _ \| |          |  ___(_)     |  __ \|  _ \ "
echo " | |_) | |_   _  ___| |_   _ _ __ | |  | | |_) |"
echo " |  _ <| | | | |/ _ \  _| | | '_ \| |  | |  _ < "
echo " | |_) | | |_| |  __/ |   | | | | | |__| | |_) |"
echo " |____/|_|\__,_|\___|_|   |_|_| |_|_____/|____/ "
echo -e "${NC}"
echo -e "  ${CYAN}:: Universal Reference Curator :: $VERSION ::${NC}"
echo "===================================================="

# ==============================================================================
# 3. PRE-FLIGHT CHECKS & CREDENTIALS
# ==============================================================================

echo -e "${BOLD}1. System Diagnostics${NC}"

# A. Check for Required Commands
REQUIRED_TOOLS="python3 cutadapt vsearch makeblastdb realpath"
MISSING_TOOLS=0
for cmd in $REQUIRED_TOOLS; do
    if ! command -v "$cmd" >/dev/null 2>&1; then
        echo -e "   [${RED}FAIL${NC}] Command '$cmd' not found."
        MISSING_TOOLS=1
    fi
done
if [ $MISSING_TOOLS -eq 1 ]; then exit 2; fi
echo -e "   [${GREEN}OK${NC}] External dependencies verified."

# B. Check Scripts Existence
# Updated to point to the new core script name if you renamed it, otherwise keeps original
CORE_SCRIPT="core/Bluefindb_core.py" 
RESCUE_SCRIPT="core/rescue_core.py"
QC_SCRIPT="./Bluefindb_qc.sh"

if [ ! -f "$CORE_SCRIPT" ] || [ ! -f "$RESCUE_SCRIPT" ]; then
    echo -e "   [${RED}FAIL${NC}] Core scripts missing in 'core/' directory."
    exit 10
fi
if [ ! -f "$QC_SCRIPT" ]; then
    echo -e "   [${RED}FAIL${NC}] QC script not found."
    exit 11
fi
chmod +x "$QC_SCRIPT"

# C. Setup Credentials
EMAIL="${BFD_EMAIL}"
if [ -z "$EMAIL" ]; then
    echo -e "   [${RED}FAIL${NC}] BFD_EMAIL environment variable is not set."
    exit 1
fi
export ENTREZ_EMAIL="$EMAIL" 

if [ -n "${NCBI_API_KEY}" ]; then
    export NCBI_API_KEY="${NCBI_API_KEY}"
    API_MSG="(API Key Detected)"
else
    API_MSG="(No API Key - Rate limits will apply)"
fi

# D. NCBI Server Status Check
echo -ne "   [....] Pinging NCBI Database..."
python3 "$CORE_SCRIPT" --check-status
STATUS_CHECK=$?

if [ $STATUS_CHECK -eq 0 ]; then
    echo -e "\r   [${GREEN}OK${NC}] NCBI Connection Established $API_MSG"
elif [ $STATUS_CHECK -eq 5 ]; then
    echo -e "\r   [${RED}FAIL${NC}] Email configuration error in Python."
    exit 1
else
    echo -e "\r   [${RED}FAIL${NC}] NCBI is unreachable (Code: $STATUS_CHECK)."
    exit 99
fi

# ==============================================================================
# 4. ACTIVE CONFIGURATION & ETA
# ==============================================================================

# Set Active Primers & Lengths based on Mode
if [[ "$MARKER" == "12S" ]]; then
    FWD_1="$P1_12S_F"; REV_1="$P1_12S_R"
    FWD_2="$P2_12S_F"; REV_2="$P2_12S_R"
    MIN_L=$LEN_12S_MIN; MAX_L=$LEN_12S_MAX
    
elif [[ "$MARKER" == "COI" ]]; then
    FWD_1="$P1_COI_F"; REV_1="$P1_COI_R"
    FWD_2=""; REV_2=""
    MIN_L=$LEN_COI_MIN; MAX_L=$LEN_COI_MAX
    
elif [[ "$MARKER" == "16S" ]]; then
    # Universal 16S Mode: Pass NONE to trigger Python internal suites
    FWD_1="NONE"; REV_1="NONE"
    FWD_2=""; REV_2=""
    MIN_L=$LEN_16S_MIN; MAX_L=$LEN_16S_MAX
    echo -e "   • Mode:            ${GREEN}Universal 16S Suite Active${NC}"
    
elif [[ "$MARKER" == "18S" ]]; then
    # Universal 18S Mode: Pass NONE to trigger Python internal suites
    FWD_1="NONE"; REV_1="NONE"
    FWD_2=""; REV_2=""
    MIN_L=$LEN_18S_MIN; MAX_L=$LEN_18S_MAX
    echo -e "   • Mode:            ${GREEN}Universal 18S Suite Active${NC}"
fi

# Paths
SPECIES_LIST_PATH=$(realpath "$SPECIES_LIST")
DB_NAME="${DB_NAME_BASE}_${MARKER}"
RAW_DIR="output/raw_${MARKER}"
FILTERED_DIR="output/curated_${MARKER}"
LOG_DIR="output/logs"

# Calculate ETA
TOTAL_SPECIES=$(grep -cve '^\s*$' "$SPECIES_LIST_PATH" || echo "0")
EST_SECONDS=$(( TOTAL_SPECIES * (DEFAULT_SLEEP + 1) ))
EST_MINUTES=$(( EST_SECONDS / 60 ))
if [ "$EST_MINUTES" -eq 0 ]; then EST_MINUTES="<1"; fi

# Display Dashboard
echo -e "${BOLD}Mission Profile${NC}"
echo -e "   • Target Marker:   ${CYAN}$MARKER${NC} ($MIN_L - $MAX_L bp)"
echo -e "   • Species Count:   ${YELLOW}$TOTAL_SPECIES${NC}"
echo -e "   • Estimated Time:  ${YELLOW}~$EST_MINUTES minutes${NC}"
echo -e "   • Input List:      $(basename "$SPECIES_LIST_PATH")"
echo -e "   • Output DB:       $DB_NAME"
echo "----------------------------------------------------"

# Clean Setup
mkdir -p "$RAW_DIR" "$FILTERED_DIR" "$LOG_DIR"

# ==============================================================================
# 5. EXECUTION STAGE 1: CORE DOWNLOAD
# ==============================================================================
echo -e "\n${BOLD}3. Initializing Core Engine${NC}"
echo -e "   🚀 Launching retrieval and curation..."

# Construct Python Arguments
PYTHON_ARGS=(
    --list "$SPECIES_LIST_PATH"
    --sleep "$DEFAULT_SLEEP"
    --raw-dir "$RAW_DIR"
    --filtered-dir "$FILTERED_DIR"
    --f-primer "$FWD_1" --r-primer "$REV_1"
    --blast-name "$DB_NAME"
    --marker-name "$MARKER"
    --min-len "$MIN_L" 
    --max-len "$MAX_L"
)

if [ -n "$FWD_2" ]; then
    PYTHON_ARGS+=(--alt-f-primer "$FWD_2" --alt-r-primer "$REV_2")
fi

if [ "$NON_INTERACTIVE" = true ]; then
    PYTHON_ARGS+=(--non-interactive)
fi

# Run Core
python3 -W ignore "$CORE_SCRIPT" "${PYTHON_ARGS[@]}" 
CORE_EXIT_CODE=$?

if [ $CORE_EXIT_CODE -ne 0 ]; then
    echo -e "   [${RED}FATAL${NC}] Core Python script failed. Check $LOG_DIR/bluefindb_run.log."
    exit 3
fi

# ==============================================================================
# 6. EXECUTION STAGE 2: AUTOMATED RESCUE
# ==============================================================================
echo -e "\n${BOLD}4. Rescue & Gap Analysis${NC}"

ACCESSIONS_FILE="output/${MARKER}_accessions.tsv"
PATCH_SCRIPT="fix_taxonomy.sh"
FINAL_FASTA="output/${DB_NAME}.fasta"

if grep -E "FAILED|NA\tSUCCESS" "$ACCESSIONS_FILE" > /dev/null; then
    echo -e "   ${YELLOW}🔬 Gaps detected. Initializing Rescue Protocols...${NC}"
    
    RESCUE_ARGS=(
        --accessions-file "$ACCESSIONS_FILE"
        --blast-name "$DB_NAME"
        --raw-dir "$RAW_DIR"
        --f-primer "$FWD_1" --r-primer "$REV_1"
        --marker-name "$MARKER"
        --min-len "$MIN_L"
        --max-len "$MAX_L"
    )
    
    if [ -n "$FWD_2" ]; then
        RESCUE_ARGS+=(--alt-f-primer "$FWD_2" --alt-r-primer "$REV_2")
    else
        RESCUE_ARGS+=(--alt-f-primer "NONE" --alt-r-primer "NONE")
    fi

    python3 -W ignore "$RESCUE_SCRIPT" "${RESCUE_ARGS[@]}"
    
    if [ $? -ne 0 ]; then
        echo -e "   [${RED}FATAL${NC}] Rescue core failed."
        exit 4
    fi
    
    # --- PATCH EXECUTION LOGIC ---
    echo -e "   ⚠️  Patch Script Generated: ${YELLOW}$PATCH_SCRIPT${NC}"
    
    if command -v dos2unix >/dev/null 2>&1; then
        dos2unix "$PATCH_SCRIPT" >/dev/null 2>&1
    fi
    chmod +x "$PATCH_SCRIPT"

    if [ "$NON_INTERACTIVE" = true ]; then
        echo "   [AUTO] Executing patch script automatically..."
        ./"$PATCH_SCRIPT" --non-interactive
        if [ $? -ne 0 ]; then
             echo "   [ERROR] Auto-patch failed."
             exit 6
        fi
    else
        echo -e "   ${CYAN}ACTION REQUIRED:${NC} Please open '$PATCH_SCRIPT' and add missing Accession IDs."
        echo -e "   (Press [Enter] to continue once you have saved your changes...)"
        read -p ""
        
        if command -v dos2unix >/dev/null 2>&1; then
            dos2unix "$PATCH_SCRIPT" >/dev/null 2>&1
        fi
        
        ./"$PATCH_SCRIPT"
        if [ $? -ne 0 ]; then
            echo -e "   [${RED}FATAL${NC}] Patch script failed execution."
            exit 6
        fi
    fi

else
    # ==============================================================================
    # 7. EXECUTION STAGE 3: QC (If no rescue needed)
    # ==============================================================================
    echo -e "   [${GREEN}OK${NC}] Database integrity nominal. Proceeding to QC."
    QC_LOG_FILE="$LOG_DIR/final_qc_and_build.log"
    
    ./"$QC_SCRIPT" "$DB_NAME" "$FINAL_FASTA" "$QC_LOG_FILE" "$FWD_1" "$REV_1" "$VERSION"
    
    if [ $? -ne 0 ]; then
        echo -e "   [${RED}FATAL${NC}] QC script failed."
        exit 5
    fi
fi

# --- FINAL SUMMARY ---
echo -e "\n${BOLD}5. Smart Engine Complete${NC}"
echo -e "   [${GREEN}SUCCESS${NC}] Pipeline Finalized."
echo "   • Final FASTA: $FINAL_FASTA"
echo "   • BLAST DB:    ${DB_NAME}.n*"
echo "   • Detailed Log: $LOG_DIR/bluefindb_run.log"
echo ""
echo "   ----------------------------------------"
echo "   Thank you for using BlueFinDB."
echo "   Author: SUBRAMANIAM VIJAYAKUMAR"
echo "   ----------------------------------------"
echo ""