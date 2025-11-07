#!/usr/bin/env python3

import sys
import os
import time
import argparse
import subprocess
import datetime 
from Bio import SeqIO, Entrez
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# --- CONFIG ---
BATCH_SIZE = 10 
MIFISH_U_F = "GTCGGTAAAACTCGTGCCAGC"
MIFISH_U_R = "CATAGTGGGGTATCTAATCCCAGTTTG"
MIFISH_E_F = "GTTGGTAAATCTCGTGCCAGC"
MIFISH_E_R = "CATAGTGGGGTATCTAATCCTAGTTTG" 
REQUIRED_RANKS = ["class", "order", "family", "genus"]
CUTADAPT_LOG = "cutadapt_details.log"
# --- END CONFIG ---

# --- LOGGING SETUP ---
LOG_FILE = "bluefindb_run.log"
def log(message, level="INFO"):
    """Appends a timestamped message to the log file."""
    timestamp = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    log_message = f"[{timestamp}][{level}] {message}"
    print(log_message) 
    with open(LOG_FILE, 'a') as f:
        f.write(log_message + "\n")
# --- END LOGGING SETUP ---


def setup_entrez():
    """Sets up Entrez email and API key from environment variables."""
    user_email = os.environ.get("ENTREZ_EMAIL")
    if not user_email or user_email == "your@email.com" or "@" not in user_email:
        log("FATAL ERROR: NCBI requires a valid email address to be provided via --email flag.", "CRIT")
        sys.exit(5)
    
    Entrez.email = user_email
    log(f"Configuring Entrez with email: {user_email}", "SETUP")
    
    api_key = os.environ.get("NCBI_API_KEY")
    if api_key:
        Entrez.api_key = api_key
        log("NOTICE: Using provided NCBI API Key for higher throughput.", "SETUP")


def get_full_lineage(accession_id):
    """
    Retrieves TaxID and full scientific lineage from NCBI.
    """
    tax_info = {rank.upper(): "NA" for rank in REQUIRED_RANKS}
    tax_info["TAXID"] = "NA"
    simple_acc = accession_id.split('.')[0] 

    try:
        # Step 1: Search Nucleotide DB with Accession to get the unique numeric ID
        handle_search = Entrez.esearch(db="nucleotide", term=simple_acc, retmax="1")
        record_search = Entrez.read(handle_search)
        handle_search.close()
        
        if not record_search["IdList"]:
             raise ValueError(f"Accession {simple_acc} not found in Nucleotide DB.")
        
        nucleotide_id = record_search["IdList"][0]

        # Step 2: Get TaxID from DocSum (using the numeric ID)
        handle_summary = Entrez.esummary(db="nucleotide", id=nucleotide_id)
        record_summary = Entrez.read(handle_summary)
        handle_summary.close()

        docsum = record_summary[0]
        raw_taxid = docsum.get('TaxId', 0)
        
        if raw_taxid == 0 or raw_taxid == 'NA':
            raise ValueError("TaxId field missing from Entrez DocSum.")
        
        tax_info["TAXID"] = raw_taxid 
        
        # Step 3: Fetch the full taxonomic record and parse lineage
        handle_fetch = Entrez.efetch(db="taxonomy", id=raw_taxid, retmode="xml")
        record_fetch = Entrez.read(handle_fetch)
        handle_fetch.close()

        lineage_ex = record_fetch[0].get("LineageEx", [])
        
        for rank_data in lineage_ex:
            rank = rank_data.get("Rank")
            name = rank_data.get("ScientificName")
            if rank in REQUIRED_RANKS:
                tax_info[rank.upper()] = name.replace(" ", "_")
        
        taxid_val = tax_info['TAXID'].real if hasattr(tax_info['TAXID'], 'real') else tax_info['TAXID']
        log(f"Taxonomy SUCCESS. TaxID: {taxid_val}. Order: {tax_info['ORDER']}", "DEBUG")
        return tax_info

    except Exception as e:
        error_msg = str(e).splitlines()[0] if e.args else "Unknown API Error"
        log(f"Taxonomy FAILURE for Accession {simple_acc}. Error: {error_msg}", "DEBUG")
        return tax_info 

def in_silico_pcr_trim(record, fwd_primer, rev_primer):
    """
    Performs mismatch-tolerant trimming using external cutadapt tool.
    """
    if not fwd_primer or not rev_primer:
        return 0, record

    temp_fasta = "temp_trim_in.fasta"
    temp_trimmed = "temp_trim_out.fasta"

    SeqIO.write(record, temp_fasta, "fasta")
    
    try:
        # -e "0.20" allows 20% mismatch tolerance.
        command = [
            "cutadapt",
            "-g", fwd_primer,
            "-a", rev_primer,
            "-e", "0.20", 
            "--discard-untrimmed",
            "-o", temp_trimmed,
            temp_fasta
        ]
        
        result = subprocess.run(command, capture_output=True, text=True, check=False)

        # --- AUDIT LOGGING ---
        with open(CUTADAPT_LOG, 'a') as f:
            # Safely check length of trimmed file
            trimmed_len = 0
            if os.path.exists(temp_trimmed):
                trimmed_records = list(SeqIO.parse(temp_trimmed, 'fasta'))
                if trimmed_records:
                    trimmed_len = len(trimmed_records[0].seq)
            
            f.write(f"\n--- Amplicon Analysis for {record.id} ---\n")
            f.write(f"Source Length: {len(record.seq)}. Primers: {fwd_primer[:5]}... / {rev_primer[:5]}...\n")
            f.write(f"Final Amplicon Length: {trimmed_len}\n")
            f.write(f"\nCutadapt STDOUT (Alignment Details):\n{result.stdout}\n")
            f.write(f"Cutadapt STDERR (Error Messages):\n{result.stderr}\n")
        # ---------------------
        
        trimmed_records = list(SeqIO.parse(temp_trimmed, "fasta"))
        
        if trimmed_records:
            amplicon_length = len(trimmed_records[0].seq)
            log(f"Trimming SUCCESS. Length: {amplicon_length}. Primers: {fwd_primer[:5]}.../{rev_primer[:5]}...", "DEBUG")
            return amplicon_length, trimmed_records[0]
        
        return 0, record

    except Exception as e:
        log(f"Trimming failed (Non-Cutadapt error): {str(e)}", "ERROR")
        return 0, record
    finally:
        if os.path.exists(temp_fasta): os.remove(temp_fasta)
        if os.path.exists(temp_trimmed): os.remove(temp_trimmed)

def run_bluefindb_core(args):
    """The main execution function for the core logic."""
    
    if os.path.exists(LOG_FILE): os.remove(LOG_FILE)
    if os.path.exists(CUTADAPT_LOG): os.remove(CUTADAPT_LOG)
    log("BlueFinDB Core Engine Starting.")

    setup_entrez()
    
    # --- Initialize Variables and Paths ---
    raw_dir = "fish_12s_raw"
    filtered_dir = "fish_12s_filtered"
    merged_fasta = "fish_12s_merged.fasta"
    blast_fasta = "fish_12s_BLAST.fasta"
    accessions_tsv = "fish_12s_accessions.tsv"
    
    species_info = {}
    total_processed = 0
    success_count = 0
    failed_species = []
    
    try:
        with open(args.list, 'r') as f:
            species_lines = [line.strip() for line in f if line.strip() and not line.strip().startswith('#')]
    except FileNotFoundError:
        log(f"FATAL ERROR: Species list file not found: {args.list}", "CRIT")
        sys.exit(7)

    log(f"Starting download and curation of {len(species_lines)} species.", "INFO")

    with open(accessions_tsv, 'w') as accessions_file:
        accessions_file.write("Species\tAccession\tTaxID\tStatus\n")
        
        for i, line in enumerate(species_lines):
            species, family, common = "NA", "NA", "NA"
            current_accession = "NA"
            
            try:
                # Parse fields
                parts = [p.strip() for p in line.split("|")]
                species, family, common = parts[0], parts[1], parts[2] if len(parts) >= 3 else "NA"
                
                sp_key = species.replace(" ", "_")
                log(f"--- Processing {species} ({i+1}/{len(species_lines)}) ---", "INFO")
                
                # --- A. Targeted Download (FINAL OPTIMIZED QUERY) ---
                raw_out = f"{raw_dir}/{sp_key}.fasta"
                # FINAL OPTIMIZED QUERY: Targets 12S fragments (100-1500 bp) to ensure trimmability.
                query = f"{species}[ORGN] AND (12S[Title] OR 12S rRNA[Title]) AND mitochondrion[Filter] NOT uncultured[Title] NOT unverified[Title] AND 100:1500[SLEN]"
                log(f"NCBI Query: {query}", "DEBUG")
                
                esearch_cmd = ["esearch", "-db", "nucleotide", "-query", query]
                efetch_cmd = ["efetch", "-format", "fasta"]
                
                esearch_proc = subprocess.run(esearch_cmd, capture_output=True, text=True, check=True)
                efetch_proc = subprocess.run(efetch_cmd, input=esearch_proc.stdout, capture_output=True, text=True, check=True)
                
                if not efetch_proc.stdout.strip():
                    raise ValueError("No sequences found for this query.")
                
                with open(raw_out, "w") as f:
                    f.write(efetch_proc.stdout)
                
                # --- B. Core QC: Longest Sequence Selection ---
                records = list(SeqIO.parse(raw_out, "fasta"))
                best_record = max(records, key=lambda r: len(r.seq))
                
                current_accession = best_record.id.split(" ")[0]
                log(f"Longest sequence selected: {current_accession} (Length: {len(best_record.seq)})", "INFO")
                
                best_record.id = sp_key
                best_record.description = ""
                
                # --- C. Taxonomic Validation (Stage 3) ---
                tax_data = get_full_lineage(current_accession)

                # --- D. Amplicon Trimming (Stage 4 - Dual Trim Logic) ---
                final_record = best_record
                final_status = "FULL_LENGTH"

                if args.f_primer and args.r_primer:
                    primer_sets = [
                        (args.f_primer, args.r_primer, "PRIMARY"), 
                    ]
                    
                    if args.alt_f_primer and args.alt_r_primer:
                        primer_sets.append((args.alt_f_primer, args.alt_r_primer, "SECONDARY"))
                    
                    best_amplicon_len = 0
                    
                    for fwd, rev, label in primer_sets:
                        log(f"Attempting trim with {label} primer set...", "DEBUG")
                        amplicon_len, trimmed_record = in_silico_pcr_trim(best_record, fwd, rev)
                        
                        if amplicon_len > best_amplicon_len:
                            best_amplicon_len = amplicon_len
                            final_record = trimmed_record
                            final_status = f"TRIMMED_{label}"

                    if best_amplicon_len == 0:
                        log("Trimming failed for all primer sets. Keeping original sequence.", "WARN")
                        final_record = best_record
                        final_status = "FULL_LENGTH_TRIM_FAILED"
                    else:
                         log(f"Best amplicon selected (Length: {best_amplicon_len}) via {final_status}.", "INFO")
                
                # --- E. Final Record Storage and Logging ---
                filtered_out = f"{filtered_dir}/{sp_key}_best.fasta"
                SeqIO.write(final_record, filtered_out, "fasta")
                
                species_info[sp_key] = {
                    "family_user": family, 
                    "common_name": common, 
                    "accession": current_accession,
                    "length": len(final_record.seq),
                    **tax_data 
                }
                
                # Accessions file logging (using final fixed TaxID conversion)
                raw_taxid_obj = tax_data['TAXID']
                taxid_log = str(raw_taxid_obj.real) if hasattr(raw_taxid_obj, 'real') else str(raw_taxid_obj)
                accessions_file.write(f"{species}\t{current_accession}\t{taxid_log}\tSUCCESS\n")
                success_count += 1
                
            except Exception as e:
                error_msg = str(e).splitlines()[0]
                failed_species.append((species, error_msg))
                accessions_file.write(f"{species}\t{current_accession}\tNA\tFAILED: {error_msg.replace(' ', '_')}\n")
                log(f"Failed to process {species}. Reason: {error_msg}", "ERROR")
                
            finally:
                total_processed += 1
                
                # --- BATCH PROCESSING LOGIC ---
                if args.sleep > 0 and (total_processed % BATCH_SIZE == 0) and (total_processed < len(species_lines)):
                    log(f"Batch of {BATCH_SIZE} complete. Sleeping for {args.sleep} seconds...", "API_SLEEP")
                    time.sleep(args.sleep)

    # --- Stage 5: Merging and Formatting ---
    log("\n[STEP 5/6] Merging and creating final BLAST-ready FASTA...", "INFO")
    
    fasta_files = [f"{filtered_dir}/{f}" for f in os.listdir(filtered_dir) if f.endswith(".fasta")]
    if not fasta_files:
        log("NOTICE: No filtered FASTA files found. Creating empty final files.", "WARN")
        with open(merged_fasta, "w") as f: f.write("")
        with open(blast_fasta, "w") as f: f.write("")
    else:
        merged_fasta = "fish_12s_merged.fasta"
        blast_fasta = "fish_12s_BLAST.fasta"
        
        with open(merged_fasta, "w") as merge_out:
            for f in fasta_files:
                with open(f, 'r') as infile:
                    merge_out.write(infile.read())
        
        with open(merged_fasta) as infile, open(blast_fasta, "w") as outfile:
            for record in SeqIO.parse(infile, "fasta"):
                sp_key = record.id
                data = species_info.get(sp_key, {})

                acc = data.get("accession", "NA")
                
                # --- FINAL TAXID FIX: Extract the integer value ---
                raw_taxid = data.get("TAXID", "NA") 
                if raw_taxid != 'NA' and hasattr(raw_taxid, 'real'):
                    taxid = str(raw_taxid.real) 
                else:
                    taxid = str(raw_taxid)
                # --- END FIX ---
                
                clss = data.get("CLASS", "NA")
                order = data.get("ORDER", "NA")
                family_ncbi = data.get("FAMILY", "NA")
                family_user = data.get("family_user", "NA")
                common = data.get("common_name", "NA").replace(" ", "_")
                
                final_family = family_ncbi if family_ncbi != "NA" else family_user

                header = f">{sp_key}|ACCESSION:{acc}|TAXID:{taxid}|CLASS:{clss}|ORDER:{order}|FAMILY:{final_family}|COMMON_NAME:{common}"
                outfile.write(header + "\n" + str(record.seq) + "\n")

    # --- Final Summary ---
    log("-" * 50)
    log(f"Final FASTA created: {blast_fasta}", "INFO")
    log(f"Accessions list created: {accessions_tsv}", "INFO")
    log(f"Trimming audit log created: {CUTADAPT_LOG}", "INFO")
    log(f"Run complete. Processed {total_processed} species. Success: {success_count}. Failed: {len(failed_species)}.", "INFO")
    if failed_species:
        log("\n⚠️ FAILED SPECIES SUMMARY:", "WARN")
        for sp, reason in failed_species:
            log(f"- {sp} (Reason: {reason})", "WARN")
    log("[END OF CORE ENGINE]", "INFO")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="BlueFinDB Core Sequence Curation Engine.")
    parser.add_argument("--list", required=True, help="Path to the species list file.")
    parser.add_argument("--sleep", type=int, default=15, help="Seconds to wait between NCBI queries.")
    parser.add_argument("--f-primer", default=None, help="Core Forward primer sequence for trimming.")
    parser.add_argument("--r-primer", default=None, help="Core Reverse primer sequence for trimming.")
    parser.add_argument("--alt-f-primer", default=None, help="Alternative Core Forward primer sequence (e.g., MiFish-E).")
    parser.add_argument("--alt-r-primer", default=None, help="Alternative Core Reverse primer sequence (e.g., MiFish-E).")
    
    run_bluefindb_core(parser.parse_args())