#!/bin/bash
# =====================================================================
# BlueFinDB QC Engine v3.0
# [PUBLICATION RELEASE]
#
# Author: SUBRAMANIAM VIJAYAKUMAR
#
# DESCRIPTION:
# Performs final Quality Control (QC) on the curated dataset:
# 1. Filters ambiguous bases (Ns).
# 2. Dereplicates sequences (100% identity).
# 3. Removes chimeric sequences (UCHIME).
# 4. Builds final BLAST database and Manifest.
# =====================================================================

# Exit immediately if a command exits with a non-zero status.
set -e

# --- 1. RECEIVE VARIABLES ---
# Check if required variables are provided
if [ $# -ne 6 ]; then
    echo "FATAL ERROR: Not all required arguments were passed to QC script."
    echo "Usage: ./Bluefindb_qc.sh <BLAST_DB_NAME> <FASTA_FILE> <QC_LOG_FILE> <FWD_PRIMER> <REV_PRIMER> <VERSION>"
    exit 1
fi

BLAST_DB_NAME="$1"
FASTA_FILE="$2"
QC_LOG_FILE="$3"
FWD_PRIMER_SEQ="$4"
REV_PRIMER_SEQ="$5"
SCRIPT_VERSION="$6"

PATCH_LOG="output/logs/taxonomy_patch.log" 

# --- 2. START LOGGING ---
echo -e "\n[STAGE 3] Starting final database QC (VSEARCH) and re-index..." | tee -a "$PATCH_LOG"
echo "--- Starting Final QC and Build (BlueFinDB $SCRIPT_VERSION) ---" > "$QC_LOG_FILE"
echo "Run Date: $(date)" >> "$QC_LOG_FILE"

TEMP_N_FILTER="temp_n_filtered.fasta"
TEMP_DEREPEATED="temp_derep.fasta"
TEMP_CHIMERA_FREE="temp_chimera_free.fasta"

# --- 3. QC PIPELINE ---
PRE_QC_COUNT=$(grep -c ">" "$FASTA_FILE")
echo "[QC STEP 1] Initial sequence count (automated + rescued): $PRE_QC_COUNT" | tee -a "$QC_LOG_FILE" "$PATCH_LOG"

# --- QC STEP: Filter Ns ---
echo "[QC STEP 2] Filtering sequences (Max Ns: 1)..." | tee -a "$QC_LOG_FILE" "$PATCH_LOG"
vsearch --fastx_filter "$FASTA_FILE" \
    --fastq_maxns 1 \
    --fastaout "$TEMP_N_FILTER" >> "$QC_LOG_FILE" 2>&1

if [ ! -s "$TEMP_N_FILTER" ] && [ "$PRE_QC_COUNT" -gt 0 ]; then
    echo "   [QC WARNING] VSEARCH N-filter produced an empty file. Retaining original sequences." | tee -a "$QC_LOG_FILE" "$PATCH_LOG"
    cp "$FASTA_FILE" "$TEMP_N_FILTER"
fi

POST_N_COUNT=$(grep -c ">" "$TEMP_N_FILTER")
REMOVED_N_COUNT=$(( PRE_QC_COUNT - POST_N_COUNT ))
echo "   [QC NOTE] Removed $REMOVED_N_COUNT sequences (failed N-filter)." | tee -a "$QC_LOG_FILE" "$PATCH_LOG"

# --- QC STEP: Dereplication ---
echo "[QC STEP 3] Removing 100% identical duplicates..." | tee -a "$QC_LOG_FILE" "$PATCH_LOG"
vsearch --derep_fulllength "$TEMP_N_FILTER" --output "$TEMP_DEREPEATED" --sizeout >> "$QC_LOG_FILE" 2>&1
POST_DEREP_COUNT=$(grep -c ">" "$TEMP_DEREPEATED")
REMOVED_DUP_COUNT=$(( POST_N_COUNT - POST_DEREP_COUNT ))
echo "   [QC NOTE] Removed $REMOVED_DUP_COUNT duplicate sequences." | tee -a "$QC_LOG_FILE" "$PATCH_LOG"

# --- QC STEP: Chimera Check ---
echo "[QC STEP 4] Checking for and removing chimeras..." | tee -a "$QC_LOG_FILE" "$PATCH_LOG"
vsearch --uchime_denovo "$TEMP_DEREPEATED" --nonchimeras "$TEMP_CHIMERA_FREE" >> "$QC_LOG_FILE" 2>&1
FINAL_COUNT=$(grep -c ">" "$TEMP_CHIMERA_FREE")
REMOVED_CHIM_COUNT=$(( POST_DEREP_COUNT - FINAL_COUNT ))
echo "   [QC NOTE] Removed $REMOVED_CHIM_COUNT chimeric sequences." | tee -a "$QC_LOG_FILE" "$PATCH_LOG"

# --- Cleanup and move ---
mv "$TEMP_CHIMERA_FREE" "$FASTA_FILE"
rm -f "$TEMP_N_FILTER" "$TEMP_DEREPEATED"

# --- 4. BUILD AND MANIFEST ---

echo "[QC STEP 5] Building final BLAST database..." | tee -a "$QC_LOG_FILE" "$PATCH_LOG"

# Derive output directory from the FASTA file path
DB_DIR=$(dirname "$FASTA_FILE")
FULL_DB_PATH="${DB_DIR}/${BLAST_DB_NAME}"

# UPDATE: Set Manifest path to be inside the output directory
MANIFEST_FILE="${DB_DIR}/BlueFinDB_Manifest.txt"


find "$DB_DIR" -maxdepth 1 -name "${BLAST_DB_NAME}.n*" -delete

# Build database specifically in the output directory
makeblastdb -in "$FASTA_FILE" -dbtype nucl -out "$FULL_DB_PATH" >> "$QC_LOG_FILE" 2>&1

echo -e "\n✅ QC & Build Complete! Final unique, QC-passed sequences: $FINAL_COUNT" | tee -a "$QC_LOG_FILE" "$PATCH_LOG"

echo "--- Writing Reproducibility Manifest ---" | tee -a "$PATCH_LOG"
echo "### BlueFinDB $SCRIPT_VERSION Reproducibility Manifest ###" > "$MANIFEST_FILE"
echo "" >> "$MANIFEST_FILE"
echo "Date of Final Build: $(date)" >> "$MANIFEST_FILE"
echo "Core Engine Log: output/logs/bluefindb_run.log" >> "$MANIFEST_FILE"
# Only list patch logs if they exist
if [ -f "$PATCH_LOG" ]; then
    echo "Rescue Patch Log: $PATCH_LOG" >> "$MANIFEST_FILE"
fi
echo "Final QC Log: $QC_LOG_FILE" >> "$MANIFEST_FILE"
echo "" >> "$MANIFEST_FILE"
echo "--- Filters Applied ---" >> "$MANIFEST_FILE"
echo "Primers (Primary): $FWD_PRIMER_SEQ / $REV_PRIMER_SEQ" >> "$MANIFEST_FILE"
echo "Length Filter: (Dynamic based on Marker Type)" >> "$MANIFEST_FILE"
echo "Max Ambiguous Ns: 1" >> "$MANIFEST_FILE"
echo "" >> "$MANIFEST_FILE"
echo "--- Curation Summary ---" >> "$MANIFEST_FILE"
echo "Final Unique Sequences: $FINAL_COUNT" >> "$MANIFEST_FILE"

echo "QC process finished successfully."