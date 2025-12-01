#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# ==============================================================================
# BlueFinDB Core Engine v3.0 - Universal Reference Database Builder
# [PUBLICATION RELEASE - CASCADE EDITION]
#
# Author: SUBRAMANIAM VIJAYAKUMAR
#
# DESCRIPTION:
# The computational heart of BlueFinDB. Handles:
# 1. Hybrid Data Mining:
#    - 18S: Checks Local SILVA DB first (Priority 1) -> Falls back to NCBI (Priority 2).
#    - Others: Standard NCBI Tier 1/Tier 2.
# 2. In-Silico PCR (Primer Trimming with Cutadapt).
# 3. Universal Primer Sweeping (16S/18S Dragnets).
# 4. Taxonomy Validation & Data Caching.
# ==============================================================================

import sys
import os
import time
import argparse
import subprocess
import datetime 
import re
import shutil
from Bio import SeqIO, Entrez
from Bio.Seq import Seq
from io import StringIO

# --- PROFESSIONAL UI IMPORTS ---
try:
    from rich.console import Console
    from rich.panel import Panel
    from rich.text import Text
    from rich.live import Live
    from rich.spinner import Spinner
    from rich.table import Table
    from rich import box
    RICH_AVAILABLE = True
    console = Console()
except ImportError:
    RICH_AVAILABLE = False
    # Fallback if rich is not installed

# --- CONFIGURATION CONSTANTS ---
VERSION = "v3.0-Cascade" 
REQUIRED_RANKS = ["class", "order"] 
CUTADAPT_LOG = "output/logs/cutadapt_details.log"
LOG_FILE = "output/logs/bluefindb_run.log"
BATCH_SIZE = 5

# --- DATABASE CONSTANTS (18S Optimization) ---
SILVA_URL = "https://www.arb-silva.de/fileadmin/silva_databases/release_138_1/Exports/SILVA_138.1_SSURef_NR99_tax_silva.fasta.gz"
LOCAL_DB_DIR = "Data/databases"
SILVA_FILENAME = "silva_138.1_ssu.fasta"
SILVA_PATH = os.path.join(LOCAL_DB_DIR, SILVA_FILENAME)

# --- PRIMER PRESETS (Fallbacks) ---
PRIMER_PRESETS = {
    "12S": { "type": "mito", "range": (150, 550), "desc": "Teleost Fish (MiFish)" },
    "COI": { "type": "mito", "range": (100, 800), "desc": "Invertebrates (Leray-XT)" },
    "16S": { "type": "mito", "range": (150, 600), "desc": "Vertebrates (Vences)" },
    "18S": { "type": "nuc", "range": (100, 600), "desc": "Eukaryotes (Stoeck)" }
}

# ==============================================================================
# UNIVERSAL PRIMER SUITES (The "Dragnet")
# ==============================================================================

# --- 16S SUITE: Vertebrates + Whales + Bacteria ---
SUITE_16S = [
    ("CGCCTGTTTATCAAAAACAT", "CCGGTCTGAACTCAGATCACGT", "16S_Vert_Std"),       # Frogs/Reptiles
    ("CGCCTGTTTACCAAAAACAT", "CCGGTCTGAACTCAGATCACGT", "16S_Mammal_Mod"),     # Whales/Mammals
    ("CGGTTGGGGTGACCTCGGA", "GCTGTTATCCCTAGGGTAACT", "16S_Mammal_Short"),     # Degraded DNA
    ("CGTTTTTTGGTCGACAGCC", "CCGGTCTGAACTCAGATCACGT", "16S_Cetacean_MarVer3"),# Cetaceans
    ("CCTACGGGNGGCWGCAG", "GACTACHVGGGTATCTAATCC", "16S_Bacteria_V3V4"),      # General Bacteria
    ("AACMGGATTAGATACCCKG", "ACGTCATCCCCACCTTCC", "16S_Bacteria_NoChloro")    # Water/Plants
]

# --- 18S SUITE: Eukaryotes + Fungi + Nematodes + Protozoa ---
SUITE_18S = [
    ("GTACACACCGCCCGTC", "TGATCCTTCTGCAGGTTCACCTAC", "18S_Euk_General"),      # Plankton/General
    ("CGATAACGAACGAGACCT", "ANCCATTCAATCGGTANT", "18S_Fungi_Best"),           # Fungi
    ("GGCAAGTCTGGTGCCAG", "TCCGTCAATTYCTTTAAGT", "18S_Nema_Coverage"),        # Nematodes (Worms)
    ("GGTTAAAAMGYTCGTAGTTG", "TGGTGGTGCCCTTCCGTCA", "18S_Nema_Eco"),          # Nematodes (Specific)
    ("GGGGAAGTATGGTTGCAAA", "TACAAAGGGCAGGGACGTAAT", "18S_Nema_Classic"),     # Nematodes (Classic)
    ("GGGRAACTTACCAGGTCC", "GTGATRWGRTTTACTTRT", "18S_Ciliate_Set2"),         # Rumen Protozoa
    ("GCTTTCGWTGGTAGTGTATT", "CAACTGTCTCTATKAAYCG", "18S_Ciliate_Set1")       # Ciliate Diversity
]

# ---------------------------------------------------------------
# 1. LOGGING AND UI UTILITIES
# ---------------------------------------------------------------

def log(message, level="INFO"):
    """Writes to log file and prints to console."""
    timestamp = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    log_message = f"[{timestamp}][{level}] {message}"
    
    # File Logging
    log_dir = os.path.dirname(LOG_FILE)
    if log_dir and not os.path.exists(log_dir):
        os.makedirs(log_dir)
    with open(LOG_FILE, 'a') as f:
        f.write(log_message + "\n")

    # Console Output (Rich or Plain)
    if level in ["CRIT", "ERROR"]:
        if RICH_AVAILABLE:
            console.print(f"[bold red]{message}[/bold red]")
        else:
            print(f"ERROR: {message}")
    elif level == "WARN":
        if RICH_AVAILABLE:
            console.print(f"[yellow]{message}[/yellow]")
        else:
            print(f"WARN: {message}")
    elif level == "SUCCESS":
        if RICH_AVAILABLE:
            console.print(f"[bold green]{message}[/bold green]")
        else:
            print(f"SUCCESS: {message}")
    elif level == "INFO":
        
        pass 

def print_header(version, marker, fwd, rev, alt_f, alt_r, desc):
    """Prints a publication-quality startup header."""
    if alt_f and alt_r and alt_f not in ["NONE", ""]:
        primer_text = f"P1: {fwd}\nP2: {alt_f}"
    else:
        primer_text = f"P1: {fwd}"

    if RICH_AVAILABLE:
        grid = Panel.fit(
            f"[bold white]BlueFinDB Reference Curator[/bold white] [cyan]{version}[/cyan]\n"
            f"[dim]Target Marker:[/dim] [bold yellow]{marker}[/bold yellow] ({desc})",
            title="🔬 Scientific Telemetry", border_style="blue"
        )
        console.print(grid)
    else:
        print(f"--- BlueFinDB {version} ({marker}) ---")
        print(f"Primers: {primer_text}")

def timed_sleep(seconds):
    """Professional spinner for API rate limiting."""
    if RICH_AVAILABLE:
        with Live(Spinner("dots", text="Respecting NCBI API rate limits..."), refresh_per_second=10, transient=True) as live:
            time.sleep(seconds)
    else:
        time.sleep(seconds)

# ---------------------------------------------------------------
# 2. DATABASE MANAGEMENT (SILVA AUTO-INSTALL)
# ---------------------------------------------------------------

def ensure_silva_db():
    """Checks if SILVA DB exists. If not, downloads and unzips it."""
    if os.path.exists(SILVA_PATH):
        return True
    
    if RICH_AVAILABLE: console.print(f"[bold yellow]⚠ SILVA Database not found. Downloading (One-time setup)...[/bold yellow]")
    else: print("Downloading SILVA Database...")

    try:
        os.makedirs(LOCAL_DB_DIR, exist_ok=True)
        zip_path = SILVA_PATH + ".gz"
        
        # Download
        subprocess.run(f"wget -O {zip_path} {SILVA_URL}", shell=True, check=True)
        
        # Unzip
        if RICH_AVAILABLE: console.print(f"[dim]Unzipping SILVA database...[/dim]")
        subprocess.run(f"gunzip {zip_path}", shell=True, check=True)
        
        if RICH_AVAILABLE: console.print(f"[bold green]✔ SILVA Database installed successfully![/bold green]")
        return True
    except Exception as e:
        log(f"Failed to download SILVA: {e}", "WARN")
        return False

def search_local_silva(species_name, db_path):
    """
    Fast grep-based search for species in local FASTA using awk for multiline records.
    """
    # Using awk to find the header
    cmd = f"awk '/{species_name}/{{n++; if(n>5) exit; p=1; print; next}} /^>/{{p=0}} p' {db_path}"
    
    try:
        result = subprocess.check_output(cmd, shell=True, text=True)
        if not result: return []
        
        local_records = []
        # Parse the string result as a FASTA file
        fasta_io = StringIO(result)
        
        for record in SeqIO.parse(fasta_io, "fasta"):
            # CLEANING: 
            # 1. Remove gaps 
            # 2. Convert RNA (U) to DNA (T) so primers bind
            clean_seq = str(record.seq).replace(".", "").replace("-", "").upper().replace("U", "T")
            
            # Update record with clean sequence
            record.seq = Seq(clean_seq)
            local_records.append(record)
                
        return local_records
    except subprocess.CalledProcessError:
        return [] 

# ---------------------------------------------------------------
# 3. NCBI CONNECTIVITY & ROBUSTNESS
# ---------------------------------------------------------------

def robust_esearch(db, term, retmax="1", retmode="xml", max_retries=3):
    """Wraps Entrez.esearch with retry logic to handle NCBI hiccups."""
    for attempt in range(max_retries):
        try:
            handle = Entrez.esearch(db=db, term=term, retmax=retmax, retmode=retmode)
            record = Entrez.read(handle)
            handle.close()
            return record
        except Exception as e:
            err_msg = str(e)
            if "Search Backend failed" in err_msg or "500" in err_msg or "502" in err_msg or "Network" in err_msg:
                wait_time = (attempt + 1) * 3
                log(f"NCBI Connection Hiccup: '{err_msg}'. Retrying in {wait_time}s...", "WARN")
                time.sleep(wait_time)
            else:
                raise e
    raise Exception(f"Failed to connect to NCBI after {max_retries} attempts.")

def check_entrez_status(email, api_key=None):
    """Pings NCBI to check if the service is reachable. Used by Wrapper."""
    try:
        Entrez.email = email
        if api_key: Entrez.api_key = api_key
        
        # Simple ping
        record = robust_esearch(db="nucleotide", term="Gadus morhua[ORGN]", retmax="1")
        if "IdList" in record:
            return 0 # Success
        else:
            return 10 # Connection okay, but empty result (weird)
            
    except Exception as e:
        return 99 # General Failure

def setup_entrez():
    """Configures global Entrez settings."""
    user_email = os.environ.get("ENTREZ_EMAIL")
    if not user_email: 
        log("FATAL ERROR: Entrez Email (BFD_EMAIL) not found.", "CRIT")
        sys.exit(5)
    
    Entrez.email = user_email
    api_key = os.environ.get("NCBI_API_KEY")
    if api_key:
        Entrez.api_key = api_key

# ---------------------------------------------------------------
# 4. CORE LOGIC: TAXONOMY & TRIMMING
# ---------------------------------------------------------------

def get_full_lineage(accession_id):
    """Retrieves TaxID and full scientific lineage from NCBI."""
    tax_info = {rank.upper(): "NA" for rank in REQUIRED_RANKS}
    tax_info["TAXID"] = "NA"; tax_info["FAMILY"] = "NA" 
    simple_acc = accession_id.split('.')[0] 

    try:
        record_search = robust_esearch(db="nucleotide", term=simple_acc, retmax="1")
        if not record_search["IdList"]: return tax_info
        
        # Get TaxID
        handle_summary = Entrez.esummary(db="nucleotide", id=record_search["IdList"][0])
        docsum = Entrez.read(handle_summary)[0]
        handle_summary.close()
        
        raw_taxid = docsum.get('TaxId', 0)
        if raw_taxid == 0 or raw_taxid == 'NA': return tax_info
        
        tax_info["TAXID"] = int(raw_taxid) 
        
        # Get Lineage
        handle_fetch = Entrez.efetch(db="taxonomy", id=tax_info["TAXID"], retmode="xml")
        record_fetch = Entrez.read(handle_fetch)
        handle_fetch.close()

        for rank_data in record_fetch[0].get("LineageEx", []):
            r_name = rank_data.get("Rank")
            if r_name in REQUIRED_RANKS + ["family"]:
                tax_info[r_name.upper()] = rank_data.get("ScientificName").replace(" ", "_")
        
        return tax_info

    except Exception as e:
        log(f"Taxonomy FAILURE for Accession {simple_acc}. Error: {e}", "DEBUG")
        return tax_info 

def in_silico_pcr_trim(record, fwd_primer, rev_primer, label, min_len, max_len):
    """
    Performs primer trimming using cutadapt.
    """
    if not fwd_primer or not rev_primer:
        return 0, record, "NO_PRIMERS"

    pid = os.getpid()
    temp_fasta = f"temp_trim_{pid}_in.fasta"
    temp_trimmed = f"temp_trim_{pid}_out.fasta"

    SeqIO.write(record, temp_fasta, "fasta")
    rev_comp_primer = str(Seq(rev_primer).reverse_complement()).upper() 
    
    try:
        # Cutadapt: -g (front), -b (back/both), -e (error rate)
        cmd = (f"cutadapt -g ^{fwd_primer} -b {rev_comp_primer} -e 0.20 " 
               f"--discard-untrimmed -o {temp_trimmed} {temp_fasta} >> {CUTADAPT_LOG} 2>&1")
        
        subprocess.run(cmd, shell=True, check=False)
        
        trimmed_records = list(SeqIO.parse(temp_trimmed, 'fasta'))
        
        if trimmed_records:
            amplicon_len = len(trimmed_records[0].seq)
            
            # Dynamic range check
            if not (min_len <= amplicon_len <= max_len):
                return 0, record, "FILTER_FAIL_LENGTH" 
            
            trimmed_records[0].id = record.id
            trimmed_records[0].description = f"TRIMMED_AMPLICON_LEN:{amplicon_len}"
            return amplicon_len, trimmed_records[0], f"TRIMMED_{label}"
        
        return 0, record, "PRIMER_NOT_FOUND"

    except Exception as e:
        log(f"Trim Error: {e}", "ERROR")
        return 0, record, "ERROR"
    finally:
        if os.path.exists(temp_fasta): os.remove(temp_fasta)
        if os.path.exists(temp_trimmed): os.remove(temp_trimmed)

def build_fasta_header(sp_key, data):
    """Constructs the BlueFinDB standardized header."""
    acc = data.get("accession", "NA")
    taxid = str(data.get("TAXID", "NA")) 
    family = data.get("FAMILY", "NA") if data.get("FAMILY", "NA") != "NA" else data.get("family_user", "NA")
    common = data.get("common_name", "NA").replace(" ", "_")
    return f">{sp_key}|ACCESSION:{acc}|TAXID:{taxid}|CLASS:{data.get('CLASS','NA')}|ORDER:{data.get('ORDER','NA')}|FAMILY:{family}|COMMON_NAME:{common}|CURATION:{data.get('status', 'NA')}"

def print_final_curation_report(species_info):
    """Generates the colorful table and summary stats at the end."""
    
    # 1. File Logging 
    log("\n" + "=" * 50, "SUMMARY")
    for sp, data in species_info.items():
        log(f"{sp} | {data['status']} | {data['length']}bp | {data['accession']}", "SUMMARY")

    # 2. Visual Table (Rich)
    if not RICH_AVAILABLE: return

    table = Table(title="Core Curation Results", box=box.ROUNDED, show_header=True, header_style="bold magenta")
    table.add_column("Species", style="cyan", no_wrap=True)
    table.add_column("Status", justify="center")
    table.add_column("Length", justify="right")
    table.add_column("Accession", style="dim")

    passed = 0
    failed = 0

    sorted_species = sorted(species_info.items(), key=lambda x: x[0])

    for sp_key, data in sorted_species:
        status = data.get('status', 'NA')
        length = data.get('length', 0)
        acc = data.get('accession', 'NA')
        
        if "FAILED" in status:
            status_str = f"[red]✘ {status}[/red]"
            len_str = "[red]-[/red]"
            failed += 1
        else:
            status_str = "[green]✔ PASS[/green]"
            len_str = f"[green]{length} bp[/green]"
            passed += 1
        
        table.add_row(sp_key, status_str, len_str, acc)

    console.print("\n")
    console.print(table)
    
    # Summary Panel
    summary_text = f"[bold green]Passed: {passed}[/bold green]   [bold red]Failed: {failed}[/bold red]   [bold]Total: {len(species_info)}[/bold]"
    console.print(Panel(summary_text, title="Mission Summary", border_style="yellow", expand=False))
    console.print("\n")

def run_tier_2_rescue(species, sp_key, primer_sets, raw_dir, marker_type, min_len, max_len, marker_name):
    """
    Tier 2: Downloads full genome/mitochondrion and attempts to slice the amplicon out.
    Optimized for 18S speed (Hybrid Approach).
    """
    raw_out = f"{raw_dir}/{sp_key}_GENOME.fasta"
    records = []

    # Check Cache
    if os.path.exists(raw_out) and os.path.getsize(raw_out) > 0:
         try: records = list(SeqIO.parse(raw_out, "fasta"))
         except: records = []
    
    # Download if not cached
    if not records:
        try:
            if marker_type == "mito":
                # 12S/COI/16S: Mitochondria are small, so downloading "complete genome" is fast/safe.
                query = f"{species}[ORGN] AND (complete mitochondrion[Title] OR complete genome[Title])"
            
            elif marker_name == "18S":
                
                query = f"{species}[ORGN] AND (18S[Title] OR small subunit[Title] OR SSU[Title] OR ribosomal RNA[Title]) AND 1000:50000[SLEN]"
            
            else:
                # Fallback for other nuclear markers
                query = f"{species}[ORGN] AND (complete genome[Title] OR whole genome shotgun[Title])"
            
            record_search = robust_esearch(db="nucleotide", term=query, retmax="3")
            
            if record_search["IdList"]:
                best_id = record_search["IdList"][0] # Prioritize first hit
                handle_fetch = Entrez.efetch(db="nucleotide", id=best_id, rettype="fasta", retmode="text")
                with open(raw_out, "w") as f: f.write(handle_fetch.read())
                records = list(SeqIO.parse(raw_out, "fasta"))
        except Exception:
            return None, None
            
    best_record_found = None
    
    # LOGIC FOR BEST HIT SELECTION
    if marker_name in ["16S", "18S"]:
        current_best_len = 0 # Start low, look for high
    else:
        current_best_len = float('inf') # Start high, look for low
    
    # Trim the Genome
    for potential_record in records: 
        for fwd, rev, label in primer_sets:
            # Use dynamic lengths
            amp_len, trimmed_rec, status = in_silico_pcr_trim(potential_record, fwd, rev, label, min_len, max_len)
            
            if amp_len > 0:
                # 16S/18S LOGIC: MAXIMIZE LENGTH
                if marker_name in ["16S", "18S"]:
                    if amp_len > current_best_len:
                        current_best_len = amp_len
                        data = {
                            "accession": potential_record.id.split(" ")[0],
                            "length": amp_len,
                            "status": "TRIMMED_GENOME_RESCUE"
                        }
                        best_record_found = ("TEMP", str(trimmed_rec.seq), data)
                
                # 12S/COI LOGIC: MINIMIZE LENGTH (Standard)
                else:
                    if amp_len < current_best_len:
                        current_best_len = amp_len
                        data = {
                            "accession": potential_record.id.split(" ")[0],
                            "length": amp_len,
                            "status": "TRIMMED_GENOME_RESCUE"
                        }
                        best_record_found = ("TEMP", str(trimmed_rec.seq), data)
                        
                        # Dynamic optimization break for standard logic
                        if min_len <= amp_len <= (min_len + 150):
                            break
                            
        if marker_name not in ["16S", "18S"] and min_len <= current_best_len <= (min_len + 150):
            break

    if best_record_found:
        _, seq, data = best_record_found
        tax_data = get_full_lineage(data["accession"])
        return best_record_found, tax_data
    
    return None, None

# ---------------------------------------------------------------
# 5. MAIN EXECUTION LOOP
# ---------------------------------------------------------------

def run_bluefindb_core(args):
    # [FIX START] Sanitize inputs to remove invisible characters
    if args.marker_name:
        args.marker_name = args.marker_name.strip().replace('\r', '').replace('\n', '').strip()
    if args.blast_name:
        args.blast_name = args.blast_name.strip().replace('\r', '').replace('\n', '').strip()
    if args.f_primer:
        args.f_primer = args.f_primer.strip().replace('\r', '').replace('\n', '').strip()
    if args.r_primer:
        args.r_primer = args.r_primer.strip().replace('\r', '').replace('\n', '').strip()
    # [FIX END]

    os.makedirs(os.path.dirname(LOG_FILE), exist_ok=True)
    if os.path.exists(LOG_FILE): os.remove(LOG_FILE)
    setup_entrez()
    
    # --- 18S SPECIAL: Ensure SILVA is ready ---
    if args.marker_name == "18S":
        ensure_silva_db()

    # --- DYNAMIC SETTINGS ---
    if args.min_len and args.max_len:
        min_l = args.min_len; max_l = args.max_len
        m_type = PRIMER_PRESETS[args.marker_name]["type"] if args.marker_name in PRIMER_PRESETS else "mito"
        desc = PRIMER_PRESETS[args.marker_name]["desc"] if args.marker_name in PRIMER_PRESETS else "Custom"
    else:
        # Fallback to presets if not provided
        preset = PRIMER_PRESETS.get(args.marker_name, {"range": (100, 1000), "type": "mito", "desc": "Custom"})
        min_l, max_l = preset["range"]
        m_type = preset["type"]; desc = preset["desc"]

    # Header call
    print_header(VERSION, args.marker_name, args.f_primer, args.r_primer, args.alt_f_primer, args.alt_r_primer, desc)
    
    # --- AUTOMATED PRIMER SELECTION ---
    primer_sets = []
    
    # MODE A: 16S AUTO-SUITE (Whales + Bacteria)
    if args.marker_name == "16S":
        if RICH_AVAILABLE: console.print(f"[bold cyan]ℹ 16S Mode: Activating Universal Suite (Vertebrate + Bacteria)[/bold cyan]")
        primer_sets = SUITE_16S
        
    # MODE B: 18S AUTO-SUITE (Fungi + Nematodes + Protozoa)
    elif args.marker_name == "18S":
        if RICH_AVAILABLE: console.print(f"[bold green]ℹ 18S Mode: Activating Robust Suite (Fungi/Nema/Ciliate)[/bold green]")
        primer_sets = SUITE_18S
        
    # MODE C: MANUAL / 12S / COI (Standard)
    else:
        if args.f_primer and args.r_primer:
            primer_sets.append((args.f_primer, args.r_primer, "PRIMARY"))
        if args.alt_f_primer and args.alt_r_primer:
            if args.alt_f_primer not in ["NONE", ""] and args.alt_r_primer not in ["NONE", ""]:
                primer_sets.append((args.alt_f_primer, args.alt_r_primer, "SECONDARY"))

    if not primer_sets:
        log("CRITICAL: No primers found. Exiting.", "CRIT")
        sys.exit(12)

    # Prepare Output Files
    output_dir = "output"
    os.makedirs(output_dir, exist_ok=True)
    
    # Use os.path.join for safe cross-platform path construction
    blast_fasta_path = os.path.join(output_dir, f"{args.blast_name}.fasta")
    if os.path.exists(blast_fasta_path): os.remove(blast_fasta_path)
    
    # Use the now-cleaned args.marker_name
    accessions_tsv = os.path.join(output_dir, f"{args.marker_name}_accessions.tsv")
    species_info = {} 
    success_count = 0
    
    with open(accessions_tsv, 'w') as accessions_file, \
         open(blast_fasta_path, 'a') as outfile_fasta:
        
        accessions_file.write("Species\tAccession\tTaxID\tStatus\n")
        
        try:
            with open(args.list, 'r') as f:
                species_lines = [l.strip() for l in f if l.strip() and not l.strip().startswith('#')]
        except FileNotFoundError:
            log(f"List file not found: {args.list}", "CRIT")
            sys.exit(7)

        total_species = len(species_lines)

        # Main Loop
        for i, line in enumerate(species_lines):
            if success_count > 0 and success_count % BATCH_SIZE == 0:
                timed_sleep(args.sleep)
            
            parts = [p.strip() for p in line.split("|")]
            species = parts[0]
            family = parts[1] if len(parts) > 1 else "NA"
            common = parts[2] if len(parts) > 2 else "NA"
            sp_key = species.replace(" ", "_")
            
            # Colorful Progress
            if RICH_AVAILABLE:
                console.print(f"[blue][{i+1}/{total_species}][/blue] Processing: {species}...")
            else:
                log(f"Processing {species}...", "INFO")
            
            # --- TIER 1: HYBRID RETRIEVAL (CASCADE) ---
            best_record_found = None
            raw_file = f"{args.raw_dir}/{sp_key}.fasta"
            source_used = "NONE"
            
            # 1. CHECK LOCAL SILVA (Priority 1)
            if args.marker_name == "18S" and os.path.exists(SILVA_PATH):
                local_records = search_local_silva(species, SILVA_PATH)
                if local_records:
                    if RICH_AVAILABLE: console.print(f"   [green]⚡ Found in Local SILVA DB ({len(local_records)} hits)[/green]")
                    
                    # --- TRIM LOCAL RECORDS ---
                    # Immediate processing of local data. If a valid barcode is generated, 
                    # the search terminates successfully. If trimming fails, the 
                    # process falls through to the NCBI fetch routine
                    
                    # For 16S/18S Suites, maximize length. For others, minimize.
                    if args.marker_name in ["16S", "18S"]:
                        current_best_len = 0 
                    else:
                        current_best_len = float('inf')

                    for potential_record in local_records:
                        for fwd, rev, label in primer_sets:
                            amp_len, trimmed_rec, trim_status = in_silico_pcr_trim(potential_record, fwd, rev, label, min_l, max_l)
                            
                            if amp_len > 0:
                                record_data = {
                                    "family_user": family, "common_name": common,
                                    "accession": potential_record.id.split(" ")[0],
                                    "length": amp_len, "status": trim_status
                                }
                                
                                if args.marker_name in ["16S", "18S"]:
                                    if amp_len > current_best_len:
                                        current_best_len = amp_len
                                        best_record_found = ("TEMP", str(trimmed_rec.seq), record_data)
                                        source_used = "LOCAL"
                                else:
                                    if amp_len < current_best_len:
                                        current_best_len = amp_len
                                        best_record_found = ("TEMP", str(trimmed_rec.seq), record_data)
                                        source_used = "LOCAL"

            # 2. CHECK CACHE (If Local failed or wasn't used)
            if not best_record_found:
                if os.path.exists(raw_file) and os.path.getsize(raw_file) > 0:
                     try: 
                        records = list(SeqIO.parse(raw_file, "fasta"))
                        
                        # --- TRIM CACHED RECORDS ---
                        if args.marker_name in ["16S", "18S"]: current_best_len = 0 
                        else: current_best_len = float('inf')

                        for potential_record in records:
                            for fwd, rev, label in primer_sets:
                                amp_len, trimmed_rec, trim_status = in_silico_pcr_trim(potential_record, fwd, rev, label, min_l, max_l)
                                if amp_len > 0:
                                    record_data = {
                                        "family_user": family, "common_name": common,
                                        "accession": potential_record.id.split(" ")[0],
                                        "length": amp_len, "status": trim_status
                                    }
                                    if args.marker_name in ["16S", "18S"]:
                                        if amp_len > current_best_len:
                                            current_best_len = amp_len
                                            best_record_found = ("TEMP", str(trimmed_rec.seq), record_data)
                                            source_used = "CACHE"
                                    else:
                                        if amp_len < current_best_len:
                                            current_best_len = amp_len
                                            best_record_found = ("TEMP", str(trimmed_rec.seq), record_data)
                                            source_used = "CACHE"
                     except: pass

            # --- TIER 1.5: FALLBACK TO NCBI (Priority 2) ---
            # Scenario: A record exists in the Local DB, but primer trimming failed.
            # Since the local sequence is valid but unusable in its current state,
            # the script attempts to fetch a fresh record from NCBI.
            if not best_record_found:
                if source_used == "LOCAL" and RICH_AVAILABLE: 
                     console.print(f"   [yellow]⚠ Local DB sequences failed trimming. Falling back to NCBI...[/yellow]")

                # 3. DOWNLOAD FROM NCBI (If Local/Cache failed)
                try:
                    if m_type == "mito":
                        query = f"{species}[ORGN] AND ({args.marker_name}[Title] OR mitochondrion[Title]) AND 100:20000[SLEN]"
                    else:
                        query = f"{species}[ORGN] AND ({args.marker_name}[Title] OR rRNA[Title]) AND 100:20000[SLEN]"
                    
                    rs = robust_esearch("nucleotide", query, retmax="500") 
                    if rs["IdList"]:
                        handle = Entrez.efetch(db="nucleotide", id=rs["IdList"], rettype="fasta", retmode="text")
                        with open(raw_file, "w") as f: f.write(handle.read())
                        records = list(SeqIO.parse(raw_file, "fasta"))
                        
                        # --- TRIM NCBI RECORDS ---
                        if args.marker_name in ["16S", "18S"]: current_best_len = 0 
                        else: current_best_len = float('inf')

                        for potential_record in records:
                            for fwd, rev, label in primer_sets:
                                amp_len, trimmed_rec, trim_status = in_silico_pcr_trim(potential_record, fwd, rev, label, min_l, max_l)
                                if amp_len > 0:
                                    record_data = {
                                        "family_user": family, "common_name": common,
                                        "accession": potential_record.id.split(" ")[0],
                                        "length": amp_len, "status": trim_status
                                    }
                                    if args.marker_name in ["16S", "18S"]:
                                        if amp_len > current_best_len:
                                            current_best_len = amp_len
                                            best_record_found = ("TEMP", str(trimmed_rec.seq), record_data)
                                    else:
                                        if amp_len < current_best_len:
                                            current_best_len = amp_len
                                            best_record_found = ("TEMP", str(trimmed_rec.seq), record_data)

                except Exception as e:
                    log(f"Download error {species}: {e}", "WARN")

            # --- Tier 2: Rescue (Priority 3) ---
            if best_record_found:
                h, s, d = best_record_found
                tax = get_full_lineage(d["accession"])
                d.update(tax)
                if RICH_AVAILABLE: console.print(f"   [green]✔ Direct Match ({d['length']}bp)[/green]")
            else:
                if RICH_AVAILABLE: console.print(f"   [yellow]⚠ Gene missing. Scanning genome...[/yellow]")
                
                best_record_found, rescue_tax = run_tier_2_rescue(species, sp_key, primer_sets, args.raw_dir, m_type, min_l, max_l, args.marker_name)
                if best_record_found:
                    _, seq, data = best_record_found
                    data.update({"family_user": family, "common_name": common})
                    if rescue_tax: data.update(rescue_tax)
                    if RICH_AVAILABLE: console.print(f"   [blue]✚ Rescue Successful! ({data['length']}bp obtained)[/blue]")
            
            # --- Final Writing ---
            if best_record_found:
                _, seq, data = best_record_found
                
                if "TAXID" not in data:
                    tax_data = get_full_lineage(data["accession"])
                    data.update(tax_data)
                
                has_tax_gap = (data.get("TAXID") == "NA")
                if not has_tax_gap:
                     for r in REQUIRED_RANKS: 
                         if data.get(r.upper()) == "NA": has_tax_gap = True
                
                if has_tax_gap:
                    data["status"] = "FAILED_TAXONOMY_GAP"
                    if RICH_AVAILABLE: console.print(f"   [yellow]⚠ Taxonomy Gap[/yellow]")
                    accessions_file.write(f"{species}\t{data['accession']}\t{data.get('TAXID','NA')}\tFAILED_TAXONOMY_GAP\n")
                else:
                    success_count += 1
                    accessions_file.write(f"{species}\t{data['accession']}\t{data.get('TAXID','NA')}\tSUCCESS\n")
                
                header = build_fasta_header(sp_key, data)
                outfile_fasta.write(f"{header}\n{seq}\n")
                species_info[sp_key] = data # [ADDED] Save for table
            else:
                accessions_file.write(f"{species}\tNA\tNA\tFAILED_ALL_PRIMERS\n")
                if RICH_AVAILABLE: console.print(f"   [red]✘ Failed[/red]")
                species_info[sp_key] = {"status": "FAILED", "length": 0, "accession": "NA"}

    log(f"Run Complete. Success: {success_count}", "INFO")
    
    # Visual Report
    print_final_curation_report(species_info)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="BlueFinDB Core Engine")
    
    # Mode Flags
    parser.add_argument("--check-status", action="store_true", help="Check NCBI connection and exit.")
    parser.add_argument("--non-interactive", action="store_true")
    
    # Inputs/Outputs
    parser.add_argument("--list", help="Species list file")
    parser.add_argument("--raw-dir", help="Directory for raw downloads")
    parser.add_argument("--filtered-dir", help="Directory for trimmed sequences")
    parser.add_argument("--blast-name", help="Output BLAST DB name prefix")
    parser.add_argument("--local-db", help="Path to local FASTA database (e.g., SILVA) for rapid search")
    
    # Primers
    parser.add_argument("--f-primer")
    parser.add_argument("--r-primer")
    parser.add_argument("--alt-f-primer")
    parser.add_argument("--alt-r-primer")
    
    # Configuration
    parser.add_argument("--marker-name", default="12S")
    parser.add_argument("--sleep", type=int, default=2)
    
    # Dynamic Length Args
    parser.add_argument("--min-len", type=int, help="Minimum amplicon length")
    parser.add_argument("--max-len", type=int, help="Maximum amplicon length")
    
    args = parser.parse_args()
    
    # --- 1. HANDLE STATUS CHECK ---
    if args.check_status:
        u_email = os.environ.get("ENTREZ_EMAIL")
        if not u_email:
            sys.exit(5) # Bad Config
        
        status_code = check_entrez_status(u_email, os.environ.get("NCBI_API_KEY"))
        sys.exit(status_code)

    # --- 2. EXECUTE CORE ---
    if not args.list or not args.raw_dir:
        # If running manually without args, show help
        parser.print_help()
        sys.exit(1)
        
    run_bluefindb_core(args)