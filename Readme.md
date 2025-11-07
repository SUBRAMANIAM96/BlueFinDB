# 🐟 BlueFinDB v2.0.0: High-Confidence 12S Reference Database Creator

BlueFinDB is a robust, single-command utility designed for ecologists and bioinformaticians to build publication-ready, high-fidelity 12S rRNA fish reference databases (libraries) for eDNA metabarcoding and taxonomic assignment.

It automates critical molecular and taxonomic validation steps that are often manual in simple bioinformatics workflows.

## Key Features (v2.0.0 — The Publication Version)

BlueFinDB integrates the best practices of complex pipelines (like in silico PCR and taxonomic API tracing) into a single, accessible Bash/Python workflow:
**100% Taxonomically Validated**: Automatically performs a 1:1 lookup against NCBI using the Accession ID to retrieve and integrate the official TaxID, Class, Order, and Family into the FASTA header, guaranteeing data traceability and credibility.

**Molecular Fidelity (Dual-Trimming)**: Employs mismatch-tolerant in silico PCR (cutadapt -e 0.20) with a competitive dual-trimming logic (MiFish-U and MiFish-E) to maximize the isolation of the optimal gene fragment for mixed fauna (Teleosts and Elasmobranchs).

**Targeted Curation**: Uses a sequence length filter (AND 100:1500[SLEN]) in the NCBI query to deliberately exclude massive mitochondrial genomes and prioritize trimmable 12S-specific gene fragments for reliable amplicon preparation.

**Zero-Step Deployment**: Automatically executes makeblastdb (Stage 6) after curation, providing a fully indexed, searchable database file that is immediately ready for blastn queries.

**Full Traceability**: Generates detailed audit logs (bluefindb_run.log and cutadapt_details.log) documenting every API call, trimming attempt, and final accession used for reproducible science.

## 📥 Requirements & Setup
OS: Linux, macOS, or WSL (Windows Subsystem for Linux).

Core Tools: python3 (with Biopython), NCBI EDirect (esearch, efetch), and BLAST+ (makeblastdb, blastn).

Trimming Tool: cutadapt (required if using the --mifish-preset or --trim-primers flags).
---

## 📂 Project Structure

BlueFinDB/
├── run_bluefindb.sh          # The main launcher script (user interface)
├── bluefindb_core.py         # The Python engine (all logic and validation)
├── species_list.txt          # ⬅️ User Input File
├── fish_12s_raw/             # Downloaded raw sequences
├── fish_12s_filtered/        # Best sequence per species (pre-merge)
├── BlueFinDB_12S_Fish.n* # Final BLAST index files
├── bluefindb_run.log         # Complete execution history and API calls
└── cutadapt_details.log      # Molecular audit log for trimming alignment
---

## 🔹 How to Run BlueFinDB (The One-Command Workflow)

The entire database curation process is executed using a single command that handles everything from downloading sequences to indexing the final database.

### 1. Create the Input List

Create a plain text file (e.g., my_fish_list.txt) with one species per line, using the pipe (|) separator:

Format: Genus species|Family|Common_Name

Gadus morhua|Gadidae|Atlantic cod
Scophthalmus maximus|Scophthalmidae|Turbot
Sardina pilchardus|Clupeidae|European sardine

### 2. The Execution Command

You must provide your species list file and your email address (mandatory for NCBI API access). You then choose your trimming strategy.

Recommended Command (using the robust MiFish Preset):

IMPORTANT: Replace user@email.com with your actual email address.

./run_bluefindb.sh my_list.txt \
    --email user@email.com \
    --mifish-preset \
    [--api-key ABCDEFGHIJKLMNOPQRSTUVWXYZ]

### 3. Command Line Arguments Breakdown

| Argument | Status | Description |
| :--- | :--- | :--- |
| `my_list.txt` | **Mandatory** | Path to your input file containing the list of target species. |
| **`--email`** | **Mandatory** | Your email address (required by NCBI for programmatic access). |
| **`--api-key`** | *Optional* | Your personal NCBI API key. Highly recommended for large lists to prevent rate-limiting and speed up the batch process. |
| **`--mifish-preset`** | *Optional* | Activates the Dual-Trimming Logic. Automatically loads the optimal core sequences for both MiFish-U (Teleost) and MiFish-E (Elasmobranch) and performs a competitive trim. |
| **`--trim-primers`** | *Optional* | Use this instead of `--mifish-preset` if you have your own single set of primers. Must be followed by both sequences in quotes. |
| **`--sleep`** | *Optional* (Default: 15) | Seconds to pause after every batch of 10 species to ensure stable API connection. |
3. Running a Local BLAST Query

### 3. Execution Flow and Output

The script executes the full 6-stage pipeline sequentially:

| Stage | Action | Verification |
| :--- | :--- | :--- |
| **Stage 1: Setup** | Checks system dependencies (Python, EDirect, BLAST+, cutadapt). | Console reports `NOTICE: Trimming requested. Checking for cutadapt... OK.` |
| **Stage 2-4: Curation** | Downloads data, selects longest sequence, validates taxonomy, and performs the Dual-Trimming on the gene fragment. | `bluefindb_run.log` shows `Taxonomy SUCCESS` and `Trimming SUCCESS` for each species. |
| **Stage 5: Merging** | Combines all curated sequences into `fish_12s_BLAST.fasta` with fully enriched headers. | FASTA file created. |
| **Stage 6: Deployment** | Automatically executes `makeblastdb` on the final FASTA file. | Console reports: `SUCCESS! BLAST database is fully indexed and ready to use.` |
## Final Output Files

### I. Core Products (The Final Database)

| File Name | Format/Type | Information Provided |
| :--- | :--- | :--- |
| **`fish_12s_BLAST.fasta`** | **FASTA** | The **final curated database file**. Contains optimized 12S gene fragments with fully enriched, validated taxonomic headers (TaxID, Class, Order, Family, etc.). |
| **`BlueFinDB_12S_Fish.n*`** | **BLAST Index Files** | The indexed files (.nhr, .nin, .nsq, etc.) created automatically by **`makeblastdb`**. These files make your database instantly searchable. |
| `fish_12s_merged.fasta` | **FASTA** | An intermediate file containing all filtered, curated sequences combined just before the final header enrichment step. |

### II. Audit and Reproducibility Files

| File Name | Content | Scientific Value |
| :--- | :--- | :--- |
| **`fish_12s_accessions.tsv`** | **TSV** (Tab-Separated) | **Provenance Log.** Lists the final chosen **Accession ID** and **TaxID** for every species, serving as the reproducibility manifest. |
| **`bluefindb_run.log`** | **Plain Text Log** | **Stability Audit.** Provides a chronological history of the execution, including NCBI queries, batch timing, and overall success/failure status. |
| **`cutadapt_details.log`** | **Verbose Log** | **Molecular Audit Trail.** Contains the raw output of **`cutadapt`** for every single trimming attempt. This proves in silico PCR occurred, showing alignment statistics and confirming the resulting amplicon length for validation. |
---



## 🔹 Tips & Notes

- Always run **run_bluefindb.sh** from inside the BlueFinDB folder to ensure all relative paths are found correctly.
- The system includes **Batch Processing** and a minimum sleep time to ensure responsible API usage and prevent your access from being blocked by NCBI.
- Check the **bluefindb_run.log** file first if your run is interrupted or if you encounter unexpected errors.
- A `test/` folder is included so new users can quickly verify functionality.
- The goal of the final process is to produce a 12S gene fragment (not necessarily a ∼170 bp amplicon) that is 100% validated for high-confidence identification.
---

## 👥 Contributors

**Subramaniam Vijayakumar**  
📧 [subramanyamvkumar@gmail.com]
🔗 [GitHub: SUBRAMANIAM96]



