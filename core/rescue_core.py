#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# ==============================================================================
# BlueFinDB Rescue Engine v3.0
# [PUBLICATION RELEASE]
#
# Author: SUBRAMANIAM VIJAYAKUMAR
#
# DESCRIPTION:
# Handles advanced sequence retrieval and metadata correction for species
# that failed the primary curation pipeline. Generates a 2-Stage Patch Script
# (Sequence Injection + Taxonomy Correction) compatible with BlueFinDB standards.
# ==============================================================================

import sys
import os
import argparse
import time
import datetime
import subprocess
from Bio import SeqIO, Entrez

# --- CONFIGURATION CONSTANTS ---
VERSION = "v3.0" 
REQUIRED_RANKS = ["class", "order"] 
PATCH_LOG = "output/logs/taxonomy_patch.log" 
QC_SCRIPT_NAME = "./Bluefindb_qc.sh"
QC_LOG_FILE = "output/logs/final_qc_and_build.log"

# --- UNIVERSAL PRIMER SUITES (The "Dragnet") ---
# Contains ALL primers from BlueFinDB Core.
PRIMER_SUITES = {
    "12S": [
        ("GTCGGTAAAACTCGTGCCAGC", "CATAGTGGGGTATCTAATCCCAGTTTG", "MiFish_U"),
        ("GTCGGTAAAACTCGTGCCAGC", "CATAGTGGGGTATCTAATCCTAGTTTG", "MiFish_E")
    ],
    "COI": [
        ("GGWACWGGWTGAACWGTWTAYCCYCC", "TAIACYTCIGGRTGICCRAARAAYCA", "Leray_XT")
    ],
    "16S": [
        ("CGCCTGTTTATCAAAAACAT", "CCGGTCTGAACTCAGATCACGT", "16S_Vert_Std"),       # Frogs/Reptiles
        ("CGCCTGTTTACCAAAAACAT", "CCGGTCTGAACTCAGATCACGT", "16S_Mammal_Mod"),     # Whales/Mammals
        ("CGGTTGGGGTGACCTCGGA", "GCTGTTATCCCTAGGGTAACT", "16S_Mammal_Short"),     # Degraded DNA
        ("CGTTTTTTGGTCGACAGCC", "CCGGTCTGAACTCAGATCACGT", "16S_Cetacean_MarVer3"),# Cetaceans
        ("CCTACGGGNGGCWGCAG", "GACTACHVGGGTATCTAATCC", "16S_Bacteria_V3V4"),      # General Bacteria
        ("AACMGGATTAGATACCCKG", "ACGTCATCCCCACCTTCC", "16S_Bacteria_NoChloro")    # Water/Plants
    ],
    "18S": [
        ("GTACACACCGCCCGTC", "TGATCCTTCTGCAGGTTCACCTAC", "18S_Euk_General"),      # Plankton/General
        ("CGATAACGAACGAGACCT", "ANCCATTCAATCGGTANT", "18S_Fungi_Best"),           # Fungi
        ("GGCAAGTCTGGTGCCAG", "TCCGTCAATTYCTTTAAGT", "18S_Nema_Coverage"),        # Nematodes (Worms)
        ("GGTTAAAAMGYTCGTAGTTG", "TGGTGGTGCCCTTCCGTCA", "18S_Nema_Eco"),          # Nematodes (Specific)
        ("GGGGAAGTATGGTTGCAAA", "TACAAAGGGCAGGGACGTAAT", "18S_Nema_Classic"),     # Nematodes (Classic)
        ("GGGRAACTTACCAGGTCC", "GTGATRWGRTTTACTTRT", "18S_Ciliate_Set2"),         # Rumen Protozoa
        ("GCTTTCGWTGGTAGTGTATT", "CAACTGTCTCTATKAAYCG", "18S_Ciliate_Set1")       # Ciliate Diversity
    ]
}

# --- DEFAULT PRIMER FALLBACKS (Safety Net) ---
FALLBACK_PRIMERS = {
    "12S": ("GTCGGTAAAACTCGTGCCAGC", "CATAGTGGGGTATCTAATCCCAGTTTG"), 
    "16S": ("CGCCTGTTTATCAAAAACAT", "CCGGTCTGAACTCAGATCACGT"),     
    "18S": ("GTACACACCGCCCGTC", "TGATCCTTCTGCAGGTTCACCTAC"),       
    "COI": ("GGWACWGGWTGAACWGTWTAYCCYCC", "TAIACYTCIGGRTGICCRAARAAYCA") 
}

# --- HELPER FUNCTIONS ---

def rescue_log(message, level="INFO"):
    timestamp = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    print(f"[{timestamp}][RESCUE-{level}] {message}")

def setup_entrez():
    user_email = os.environ.get("ENTREZ_EMAIL")
    if not user_email:
        print("FATAL: ENTREZ_EMAIL environment variable is not set.")
        sys.exit(1)
    Entrez.email = user_email
    Entrez.api_key = os.environ.get("NCBI_API_KEY")
    rescue_log(f"Configuring Entrez connection...", "SETUP")

def get_taxid_from_accession(accession_id):
    """Retrieves TaxID from NCBI using the Accession number."""
    simple_acc = accession_id.split('.')[0]
    try:
        handle = Entrez.esummary(db="nucleotide", id=simple_acc, report="docsum")
        record = Entrez.read(handle)
        return int(record[0]['TaxId'])
    except Exception:
        return None

def get_full_lineage_rescue(taxid):
    """Retrieves full lineage information from NCBI Taxonomy."""
    tax_info = {rank.upper(): "NA" for rank in REQUIRED_RANKS}
    try:
        handle = Entrez.efetch(db="taxonomy", id=taxid, retmode="xml")
        records = Entrez.read(handle)
        if not records: return tax_info
        
        for rank_data in records[0].get("LineageEx", []):
            r_name = rank_data.get("Rank")
            if r_name in REQUIRED_RANKS:
                tax_info[r_name.upper()] = rank_data.get("ScientificName")
        return tax_info
    except Exception:
        return tax_info

# --- BASH SCRIPT GENERATOR (2-STAGE SYSTEM) ---

def generate_two_stage_patch_script(stage1_fixes, stage2_fixes, db_name, raw_dir, 
                                    active_suite):
    """Generates the structured bash script with Stage 1 (Seq) and Stage 2 (Tax) sections."""
    
    script_name = "fix_taxonomy.sh"
    fasta_file = f"{db_name}.fasta"
    
    # Construct Bash Arrays for Primers (Allows unlimited primers)
    fwd_list_str = " ".join([f'"{p[0]}"' for p in active_suite])
    rev_list_str = " ".join([f'"{p[1]}"' for p in active_suite])
    lbl_list_str = " ".join([f'"{p[2]}"' for p in active_suite])

    # Define the Bash Function with LOOP capability
    bash_func = r"""
# --- BASH FUNCTION: INJECT AND TRIM (MULTI-PRIMER) ---
inject_and_trim() {
    local SP_KEY="$1"; local INPUT_ACC="$2"; local TEMP_LOG="$3";
    local FASTA_FILE="$4"; local RAW_DIR="$5";
    local MIN_LEN=100; local MAX_LEN=1000;

    if [[ "$INPUT_ACC" == *"NA_FILL_ME"* || "$INPUT_ACC" == "SKIP" || "$INPUT_ACC" == "" ]]; then
        echo "  [SKIP] Skipping $SP_KEY (No Accession provided)." >> "$TEMP_LOG"
        return 0
    fi

    echo "[STEP] Processing $SP_KEY with Accession $INPUT_ACC..." | tee -a "$TEMP_LOG"
    local RAW_FASTA="${RAW_DIR}/${SP_KEY}_rescue.fasta";
    local TRIMMED_FASTA="${RAW_DIR}/${SP_KEY}_best.fasta";
    
    # 1. Download Sequence
    efetch -db nucleotide -id "$INPUT_ACC" -format fasta > "$RAW_FASTA" 2>> "$TEMP_LOG"
    if [ ! -s "$RAW_FASTA" ]; then
        echo "  [ERROR] Download failed for $INPUT_ACC. Skipping." | tee -a "$TEMP_LOG"
        return 0
    fi
    # Remove existing entry
    sed -i "/^>${SP_KEY}|/d" "$FASTA_FILE" 2>> "$TEMP_LOG"

    # 2. Iterate Through ALL Primers in Global Array
    BEST_LEN=0
    CHOSEN_FILE=""
    TRIM_STATUS="FAILED"

    # Get array length
    local LEN=${#PRIMER_SUITE_F[@]}
    
    for (( i=0; i<LEN; i++ )); do
        P_F="${PRIMER_SUITE_F[$i]}"
        P_R="${PRIMER_SUITE_R[$i]}"
        LABEL="${PRIMER_SUITE_LBL[$i]}"
        
        TEMP_TRIM="${RAW_DIR}/temp_trim_${SP_KEY}_${i}.fasta"
        
        # Reverse Complement via Python
        R_RC=$(python3 -c "from Bio.Seq import Seq; print(str(Seq(\"$P_R\").reverse_complement()))")
        
        cutadapt -g "^$P_F" -b "$R_RC" -e 0.20 --discard-untrimmed -o "$TEMP_TRIM" "$RAW_FASTA" >> "$TEMP_LOG" 2>&1
        
        CURR_LEN=0
        if [ -s "$TEMP_TRIM" ]; then
             CURR_LEN=$(grep -v ">" "$TEMP_TRIM" | tr -d "\n" | wc -c)
        fi
        
        # Logic: If valid length and longer than current best (Max Length Strategy for 18S/16S)
        if (( CURR_LEN >= MIN_LEN && CURR_LEN <= MAX_LEN )); then
             if (( CURR_LEN > BEST_LEN )); then
                 BEST_LEN=$CURR_LEN
                 CHOSEN_FILE="$TEMP_TRIM"
                 TRIM_STATUS="TRIMMED_${LABEL}"
             else
                 rm -f "$TEMP_TRIM" # Cleanup inferior
             fi
        else
             rm -f "$TEMP_TRIM" # Cleanup invalid
        fi
    done

    # 3. Fallback: Check Original Sequence Length
    if [[ "$CHOSEN_FILE" == "" ]]; then
        ORIG_LEN=$(grep -v ">" "$RAW_FASTA" | tr -d "\n" | wc -c)
        if (( ORIG_LEN >= MIN_LEN && ORIG_LEN <= MAX_LEN )); then
            CHOSEN_FILE="$RAW_FASTA"; TRIM_STATUS="PRETRIMMED_RESCUE"
            echo "  [OK] Primers failed, but original length ($ORIG_LEN) is valid. Keeping." >> "$TEMP_LOG"
        fi
    fi

    if [[ "$CHOSEN_FILE" != "" ]]; then
        mv "$CHOSEN_FILE" "$TRIMMED_FASTA"
        # Inject header structure compatible with Stage 2 sed replacement
        INJECT_HEADER=">${SP_KEY}|ACCESSION:${INPUT_ACC}|TAXID:NA|CLASS:NA|ORDER:NA|FAMILY:NA|COMMON_NAME:NA|CURATION:${TRIM_STATUS}"
        echo "$INJECT_HEADER" >> "$FASTA_FILE"
        awk '!/^>/{print}' "$TRIMMED_FASTA" >> "$FASTA_FILE"
        echo "  [SUCCESS] Injected $SP_KEY (Acc: $INPUT_ACC)." | tee -a "$TEMP_LOG"
        # Clean up other temp files if any exist
        rm -f "${RAW_DIR}/temp_trim_${SP_KEY}_"*.fasta
    else
        echo "  [FAIL] $INPUT_ACC failed all primers. Skipping." | tee -a "$TEMP_LOG"
    fi
    rm -f "$RAW_FASTA"
}
"""

    with open(script_name, "w") as f:
        # --- HEADER ---
        f.write("#!/bin/bash\n")
        f.write(f"# Auto-generated 2-Stage Rescue Script by BlueFinDB {VERSION}\n\n")
        
        # --- DEFINE ARRAYS ---
        f.write(f"declare -a PRIMER_SUITE_F=({fwd_list_str})\n")
        f.write(f"declare -a PRIMER_SUITE_R=({rev_list_str})\n")
        f.write(f"declare -a PRIMER_SUITE_LBL=({lbl_list_str})\n\n")
        
        f.write(f"BLAST_DB_NAME=\"{db_name}\"\n")
        f.write(f"FASTA_FILE=\"{fasta_file}\"\n")
        f.write(f"RAW_DIR=\"{raw_dir}\"\n")
        f.write(f"PATCH_LOG=\"{PATCH_LOG}\"\n")
        f.write(f"QC_LOG_FILE=\"{QC_LOG_FILE}\"\n")
        f.write(f"QC_SCRIPT=\"{QC_SCRIPT_NAME}\"\n\n")
        
        f.write(bash_func)
        f.write("\n")
        f.write("echo \"--- Starting Patch Process: $(date) ---\" > \"$PATCH_LOG\"\n\n")
        
        # --- STAGE 1: MISSING SEQUENCES ---
        f.write("# ===================================================================\n")
        f.write("# STAGE 1: MISSING SEQUENCES (Injection)\n")
        f.write("# ===================================================================\n\n")
        
        if not stage1_fixes:
            f.write("# No sequence failures detected.\n")
        
        for species, data in stage1_fixes.items():
            acc = data['accession'] if data['accession'] != "NA" else f"NA_FILL_ME__{species}_ACC"
            sp_var = species.upper().replace(".", "").replace(" ", "_").replace("-", "_")
            
            f.write(f"# 🐟 FIX REQUIRED: {species} (Reason: {data['status']})\n")
            f.write(f"CORRECT_ACCESSION_{sp_var}=\"{acc}\"\n")
            f.write(f"if [[ \"$CORRECT_ACCESSION_{sp_var}\" != *\"NA_FILL_ME\"* ]]; then\n")
            # NOTE: Removed P1_F/P1_R arguments because they are now global arrays
            f.write(f"  inject_and_trim \"{species}\" \"$CORRECT_ACCESSION_{sp_var}\" \"$PATCH_LOG\" \"$FASTA_FILE\" \"$RAW_DIR\"\n")
            f.write(f"fi\n\n")

        # --- STAGE 2: MISSING TAXONOMY ---
        f.write("# ===================================================================\n")
        f.write("# STAGE 2: MISSING TAXONOMY (Patching)\n")
        f.write("# ===================================================================\n\n")
        
        if not stage2_fixes:
            f.write("# No taxonomy gaps detected.\n")

        for species, data in stage2_fixes.items():
            sp_var = species.upper().replace(".", "").replace(" ", "_").replace("-", "_")
            acc = data['accession']
            tid = data['taxid'] if data['taxid'] != "NA" else f"NA_FILL_ME__{species}_TAXID"
            cls = data['class'] if data['class'] != "NA" else f"NA_FILL_ME__{species}_CLASS"
            ordr = data['order'] if data['order'] != "NA" else f"NA_FILL_ME__{species}_ORDER"
            
            f.write(f"# Taxonomy for: {species}\n")
            f.write(f"CORRECT_ACCESSION_{sp_var}=\"{acc}\"\n")
            f.write(f"CORRECT_TAXID_{sp_var}=\"{tid}\"\n")
            f.write(f"CORRECT_CLASS_{sp_var}=\"{cls}\"\n")
            f.write(f"CORRECT_ORDER_{sp_var}=\"{ordr}\"\n")
            f.write(f"if [[ \"$CORRECT_TAXID_{sp_var}\" != *\"NA_FILL_ME\"* ]]; then\n")
            f.write(f"  echo \"[INFO] Patching {species} header with verified taxonomy...\" | tee -a \"$PATCH_LOG\"\n")
            f.write(f"  sed -i \"/^>{species}/s/ACCESSION:[^|]*|TAXID:[^|]*|CLASS:[^|]*|ORDER:[^|]*/ACCESSION:$CORRECT_ACCESSION_{sp_var}|TAXID:$CORRECT_TAXID_{sp_var}|CLASS:$CORRECT_CLASS_{sp_var}|ORDER:$CORRECT_ORDER_{sp_var}/\" \"$FASTA_FILE\" 2>> \"$PATCH_LOG\"\n")
            f.write(f"  echo \"[SUCCESS] {species} metadata patched!\" | tee -a \"$PATCH_LOG\"\n")
            f.write(f"fi\n\n")

        # --- STAGE 3: QC ---
        f.write("# ===================================================================\n")
        f.write("# STAGE 3: FINAL QC\n")
        f.write("# ===================================================================\n")
        f.write("if [ -f \"$QC_SCRIPT\" ]; then\n")
        f.write("    chmod +x \"$QC_SCRIPT\"\n")
        f.write(f"    \"$QC_SCRIPT\" \"$BLAST_DB_NAME\" \"$FASTA_FILE\" \"$QC_LOG_FILE\" \"MULTI-PRIMER\" \"MULTI-PRIMER\" \"v3.0\"\n")
        f.write("else\n")
        f.write("    echo \"[ERROR] QC Script not found!\" | tee -a \"$PATCH_LOG\"\n")
        f.write("fi\n")
        f.write("echo -e \"\\n✅ Patching complete. Check logs for details.\"\n")

    # Automatic dos2unix conversion
    try:
        subprocess.run(["dos2unix", script_name], check=False, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        rescue_log(f"Applied dos2unix to {script_name}", "INFO")
    except Exception:
        pass

    rescue_log(f"Patch script generated: {script_name}", "SUCCESS")


def main():
    if not os.environ.get("ENTREZ_EMAIL"):
        print("FATAL: ENTREZ_EMAIL environment variable is not set.")
        sys.exit(1)
        
    parser = argparse.ArgumentParser()
    parser.add_argument("--accessions-file", required=True)
    parser.add_argument("--blast-name", required=True)
    parser.add_argument("--raw-dir", required=True)
    parser.add_argument("--f-primer", default="NONE")
    parser.add_argument("--r-primer", default="NONE")
    parser.add_argument("--alt-f-primer", default="NONE")
    parser.add_argument("--alt-r-primer", default="NONE")
    parser.add_argument("--marker-name", required=True)
    parser.add_argument("--min-len", type=int, default=100)
    parser.add_argument("--max-len", type=int, default=1000)
    
    args = parser.parse_args()
    setup_entrez()
    
    # --- SUITE SELECTION (Universal Mode Logic) ---
    active_suite = []
    
    # Ensure case-insensitive match for marker name
    target_marker = args.marker_name.upper()

    # Manual override has priority IF provided explicitly
    if args.f_primer != "NONE" and args.r_primer != "NONE":
        active_suite = [(args.f_primer, args.r_primer, "MANUAL_PRIMARY")]
        if args.alt_f_primer != "NONE" and args.alt_r_primer != "NONE":
            active_suite.append((args.alt_f_primer, args.alt_r_primer, "MANUAL_SECONDARY"))
        rescue_log(f"Manual Mode: Using provided primers.", "INFO")

    # Check the standardized marker against known suites
    elif target_marker in PRIMER_SUITES:
        active_suite = PRIMER_SUITES[target_marker]
        rescue_log(f"Universal Mode: Using auto-suite for {target_marker} ({len(active_suite)} pairs).", "INFO")
        
    else:
        # Fallback Mode
        rescue_log(f"Warning: Unknown marker '{target_marker}' with no manual primers. Using default fallback.", "WARN")
        if target_marker in FALLBACK_PRIMERS:
            f, r = FALLBACK_PRIMERS[target_marker]
            active_suite = [(f, r, "FALLBACK_DEFAULT")]

    # Analyze Accessions File to split failures
    seq_failures = []
    tax_failures = []
    
    if not os.path.exists(args.accessions_file):
        rescue_log("Accessions file not found.", "CRIT")
        sys.exit(1)

    with open(args.accessions_file, 'r') as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) < 4: continue
            species, acc, taxid, status = parts
            
            # Force species key to have underscores instead of spaces
            species = species.replace(" ", "_")

            if "FAILED" in status:
                if "TAXONOMY" in status:
                    tax_failures.append((species, acc))
                else:
                    seq_failures.append((species, status))

    if not seq_failures and not tax_failures:
        rescue_log("No failed species found to rescue.", "INFO")
        return 0

    rescue_log(f"Found {len(seq_failures)} sequence failures and {len(tax_failures)} taxonomy gaps.", "INFO")
    
    # --- STAGE 1 PREP: Find Accessions for missing sequences ---
    stage1_fixes = {}
    for sp_key, reason in seq_failures:
        rescue_log(f"Searching for missing sequence: {sp_key}...", "SEARCH")
        # Ensure search uses the standardized marker name
        search_term = f"{sp_key.replace('_', ' ')}[ORGN] AND ({target_marker}[Title] OR mitochondrion[Title] OR rRNA[Title])"
        
        try:
            handle = Entrez.esearch(db="nucleotide", term=search_term, retmax="1", sort="relevance")
            record = Entrez.read(handle)
            if record["IdList"]:
                current_accession = record["IdList"][0]
                h_sum = Entrez.esummary(db="nucleotide", id=current_accession)
                s_sum = Entrez.read(h_sum)
                real_acc = s_sum[0]['Caption']
                rescue_log(f"  Found: {real_acc}", "SUCCESS")
            else:
                real_acc = "NA"
                rescue_log(f"  No match found.", "WARN")
        except:
            real_acc = "NA"
            
        stage1_fixes[sp_key] = {'accession': real_acc, 'status': reason}

    # --- STAGE 2 PREP: Find Taxonomy for existing accessions ---
    stage2_fixes = {}
    for sp_key, acc in tax_failures:
        rescue_log(f"Fixing taxonomy for: {sp_key} ({acc})...", "TAX-FIX")
        
        tid = get_taxid_from_accession(acc)
        final_taxid = str(tid) if tid else "NA"
        
        class_val, order_val = "NA", "NA"
        if final_taxid != "NA":
            lin = get_full_lineage_rescue(final_taxid)
            class_val = lin.get('CLASS', 'NA')
            order_val = lin.get('ORDER', 'NA')
            
        stage2_fixes[sp_key] = {
            'accession': acc, 
            'taxid': final_taxid, 
            'class': class_val, 
            'order': order_val
        }

    # Generate Script
    generate_two_stage_patch_script(
        stage1_fixes, stage2_fixes, 
        args.blast_name, args.raw_dir,
        active_suite
    )
    
    return 0

if __name__ == "__main__":
    main()