# 🐟 BlueFinDB v3.0: High-Confidence Reference Database Creator

## Overview

BlueFinDB is an automated bioinformatics pipeline designed to construct taxonomically curated reference databases for environmental DNA (eDNA) metabarcoding. By accepting a target species list, it generates a validated, locally searchable BLAST database optimized for specific marker genes—COI, 12S, 16S, and 18S. 

The reproducibility of eDNA metabarcoding depends entirely on transparent, replicable workflows for creating complete, curated reference databases. Conventional curation tools (such as RESCRIPT, CRUX, and OBITools) often operate as rigid, inflexible systems. They frequently discard valid taxa due to minor primer mismatches or fail completely when specific marker sequences are absent from public repositories.

BlueFinDB addresses this gap with an open-source, multi-marker pipeline engineered for both flexibility and exceptional recovery rates. It constructs highly accurate, project-specific reference libraries while retaining valid data that other tools discard. By leveraging standard, publication-validated primer sets (MiFish, Leray-XT, Vences V4, and Stoeck V9), BlueFinDB ensures your database reflects your research goals rather than repository limitations.

---

## 🌟 Key Features of BlueFinDB v3.0

### 1. Enterprise-Grade NCBI Retrieval Architecture

**Fault-tolerant querying** integrates a custom `robust_esearch` wrapper with exponential backoff logic. The system automatically detects and recovers from HTTP 500/502 errors and NCBI backend failures, ensuring stable high-throughput retrieval even during network instability.

**Smart rate limiting** incorporates NCBI API key management to dynamically control request rates and prevent IP blocking, allowing seamless integration into production environments.

### 2. Tier-2 Genomic Rescue & Mining System

When standard marker sequences prove elusive, BlueFinDB doesn't simply report failure—it escalates automatically. The **Tier-2 protocol activates** when direct marker entries are unavailable, enabling comprehensive mining of full mitochondrial genomes and Whole Genome Shotgun (WGS) projects.

**Algorithmic extraction** applies precise coordinate mapping to locate and extract target amplicons from large-scale genomic records. This computational approach recovers valid taxa that conventional sequence scrapers routinely miss, effectively rescuing taxonomic information from the genomic frontier.

### 3. Advanced In-Silico PCR Simulation

Rather than relying on single primer pairs, BlueFinDB implements a sophisticated **"Dragnet" approach** that iteratively tests broad primer suites (e.g., MiFish-U vs. MiFish-E, or multiple 16S mammal/bacteria sets simultaneously). This multiplicative strategy maximizes taxonomic coverage by accounting for the biological variability inherent in natural populations.

**Degeneracy & mismatch handling** integrates Cutadapt with advanced logic to resolve IUPAC nucleotide ambiguity codes (W, Y, R, and others) and reverse complements. The system simulates laboratory PCR thermodynamics, ensuring that primer binding reflects realistic molecular biology rather than oversimplified sequence matching.

### 4. Automated Curation & Quality Control Pipeline

**Phylogenetic lineage verification** cross-references each sequence against the NCBI Taxonomy database to confirm the presence of required ranks (Class, Order, Family). Incomplete or taxonomically "orphan" entries are removed, ensuring database integrity.

**VSEARCH-based quality assurance** executes a rigorous, multi-stage QC workflow: global dereplication at 100% identity, strict ambiguous base filtering (≤1 N per sequence), and de novo chimera detection. The result is a publication-ready reference library free from PCR artifacts.

### 5. Self-Healing "Auto-Patch" Architecture

When gaps persist after Tier-1 and Tier-2 mining, the **Rescue System generates a custom Bash script** (`fix_taxonomy.sh`) with precise instructions to inject missing sequences or patch taxonomy. This allows one-step database remediation rather than manual reconstruction.

**Taxonomic metadata injection** automatically retrieves complete lineage strings and fills missing fields (Class, Order, Family), restoring database integrity without human intervention.

### 6. Specialized Feature: 18S Hybrid-Retrieval Engine

The 18S module uniquely implements a **priority-based search algorithm** that queries high-quality, locally curated repositories (SILVA 138.1 SSURef NR99) before falling back to NCBI. This hierarchical approach guarantees superior taxonomic resolution for eukaryotic markers.

**Latency & load optimization** dramatically reduces API overhead by resolving eukaryotic sequences locally—delivering ~50× faster performance compared to pure NCBI querying. Remote NCBI searches activate only when local curated matches are exhausted.

### 7. Real-Time Telemetry & Visual Analytics

**Rich console interface** via the `rich` library provides a professional, high-contrast CLI with real-time telemetry, status updates, color-coded indicators, and summary tables. Monitoring large-scale runs becomes intuitive rather than burdensome.

**Smart resume capability** leverages intelligent file checking to detect previously downloaded or trimmed sequences. If a job is interrupted, BlueFinDB resumes exactly where it left off without re-downloading gigabytes of data.

### 8. Strict Provenance & Reproducibility

**Automated manifest logging** generates a `BlueFinDB_Manifest.txt` file capturing timestamps, software version, active primer sets, and all filter parameters. This complete audit trail meets publication-grade reproducibility standards and supports transparent peer review.

---

## 🛠️ Dependencies

### Core Bioinformatics Tools

- **NCBI BLAST+** (`makeblastdb`, `blastn`): Builds and queries the final local reference database
- **VSEARCH** (`vsearch`): Performs chimera detection, dereplication, and strict quality filtering
- **Cutadapt** (`cutadapt`): Executes precise in-silico primer trimming with degeneracy handling
- **NCBI Entrez Utilities** (`esearch`, `efetch`): Powers high-throughput sequence retrieval and metadata access

### Python Requirements

- **Python 3.8+**
- **Biopython**: Handles FASTA I/O and Entrez API integration
- **Rich**: Delivers the professional CLI dashboard and progress visualization

### System Utilities

- **wget**: Required to download the SILVA database automatically
- **awk, sed, grep, gunzip**: Used extensively for efficient stream processing

### Operating System

- Linux, macOS, or WSL (Windows Subsystem for Linux)

---

## 📂 Repository Structure

The pipeline follows a modular architecture designed for clarity and maintainability:

```
BlueFinDB/
├── 📜 Bluefindb_run.sh          # 🚀 MASTER SCRIPT: Pipeline orchestrator
├── 📜 Bluefindb_qc.sh           # Standalone Quality Control (QC) module
├── 🧠 core/                     # Intelligence Engine
│   ├── Bluefindb_core.py        # Hybrid-mining & retrieval logic
│   └── rescue_core.py           # Tier-2 Rescue & Auto-Patch generator
├── 📁 Data/                     # Input directory (species lists)
├── 📂 output/                   # Results directory (databases & logs)
├── 🧪 test/                     # Validation Sandbox
│   ├── 📜 Bluefindb_test_run.sh # Test pipeline launcher
│   ├── 📜 Bluefindb_qc.sh       # QC script for testing
│   ├── 📁 Data/                 # Sample datasets
│   └── 📂 output/               # Test run results
├── ⚖️ LICENSE
└── 📖 README.md
```

---

## ⚙️ Script Architecture & Components

### 1. Bluefindb_run.sh (The Orchestrator)

This master entry point manages workflow logic without performing bioinformatic analysis directly:

- **Environment Validation**: Checks for all dependencies and API credentials before execution
- **Dynamic Configuration**: Loads active primer sequences based on user flags (--12S, --COI, --16S, --18S)
- **Flow Control**: Sequences Core Engine execution, monitors exit codes, triggers the Rescue System if needed, and launches the QC module

### 2. Bluefindb_core.py (The Mining Engine)

The computational heart of the pipeline—responsible for retrieval, in-silico PCR, and hybrid mining:

- **Hybrid Retrieval**: Implements local database checking (SILVA), result caching, and NCBI fallback queries
- **Smart Trimming**: Executes Cutadapt with configurable mismatch tolerance (default: 20%) to simulate realistic PCR conditions
- **Primer Dragnets**: For 16S and 18S markers, sweeps through multiple primer suites to maximize recovery
- **Genome Mining**: Downloads full genomes or mitochondria when direct markers are unavailable, then mathematically extracts target regions

### 3. rescue_core.py (The Self-Healing System)

Runs after core mining to fix retrieval failures without manual intervention:

- **Failure Forensics**: Analyzes the `accessions.tsv` file to identify specific failure modes (Missing Sequence vs. Missing Taxonomy)
- **Auto-Patching**: Dynamically generates a custom Bash script (`fix_taxonomy.sh`) with precise remediation instructions
- **Taxonomy Injection**: Retrieves complete lineage strings for orphan sequences, ensuring taxonomic completeness

### 4. Bluefindb_qc.sh (The Quality Assurance Module)

Ensures publication-ready output through rigorous filtering:

- **Ambiguity Filtering**: Removes sequences containing more than 1 ambiguous base (N)
- **Dereplication**: Merges 100% identical sequences to reduce redundancy
- **Chimera Detection**: Identifies and removes PCR artifacts using `uchime_denovo`
- **Database Compilation**: Generates the final BLAST-indexed reference library

---

## 💻 Usage & Execution

### Step 1: Populate Your Species List

Navigate to the `Data/` folder, where you'll find four pre-configured files:

- `species_list_12s.txt` — for 12S rRNA (Fish)
- `species_list_16S.txt` — for 16S rRNA (Vertebrates/Bacteria)
- `species_list_18S.txt` — for 18S rRNA (Eukaryotes/Fungi)
- `species_list_COI.txt` — for COI (Invertebrates)

Add your target species using the pipe-delimited format:

```
Scientific_Name|Family|Common_Name
```

**Example:**
```
Merluccius merluccius|Merlucciidae|European_hake
Gadus morhua|Gadidae|Atlantic_Cod
Dicentrarchus labrax|Moronidae|European_seabass
```

### Step 2: Configure the Master Script

Edit `Bluefindb_run.sh` and locate the **CONFIGURATION DEFAULTS** section (lines 20–25). Update the `DEFAULT_INPUT_FILE` variable to match your prepared species list:

```bash
# --- CONFIGURATION DEFAULTS ---
VERSION="v3.0"
DEFAULT_SLEEP=1
DB_NAME_BASE="BlueFinDB"

# CHANGE THIS LINE to your specific file
DEFAULT_INPUT_FILE="Data/species_list_16S.txt"
```

### Step 3: Set Environment Variables

Before running the pipeline, configure your NCBI credentials in your terminal:

```bash
export BFD_EMAIL="your.email@university.edu"          # REQUIRED
export NCBI_API_KEY="abcdef1234567890"                # RECOMMENDED
```

### Step 4: Verify Your Setup (Optional)

Use the built-in diagnostic tool to validate your environment:

```bash
bash Bluefindb_run.sh --check-status
```

### Step 5: Run the Pipeline

Execute the pipeline using the appropriate marker flag:

**🐟 For 12S (Fish)** — Targets mitochondrial 12S rRNA
```bash
bash Bluefindb_run.sh --12S
```

**🦀 For COI (Invertebrates)** — Targets mitochondrial COI
```bash
bash Bluefindb_run.sh --COI
```

**🐋 For 16S (Vertebrates/Mammals)** — Activates the Universal 16S Suite
```bash
bash Bluefindb_run.sh --16S
```

**🍄 For 18S (Eukaryotes/Fungi)** — Activates the Hybrid Mining Engine
```bash
bash Bluefindb_run.sh --18S
```

---

## 🔬 Targeted Primer Sets

BlueFinDB deploys publication-validated primer sets optimized for each marker. The system distinguishes between **Static Markers** (12S/COI with fixed primer pairs) and **Dynamic Markers** (16S/18S with sweeping primer suites).

### 12S rRNA (Fish) — Industry-Standard MiFish Primers

| Primer Set | Direction | Sequence (5' → 3') | Target Organism |
|:-----------|:----------|:------------------|:----------------|
| **MiFish-U** | Forward | `GTCGGTAAAACTCGTGCCAGC` | Bony Fish (Teleosts) |
| **MiFish-U** | Reverse | `CATAGTGGGGTATCTAATCCCAGTTTG` | Bony Fish (Teleosts) |
| **MiFish-E** | Forward | `GTCGGTAAAACTCGTGCCAGC` | Sharks & Rays (Elasmobranchs) |
| **MiFish-E** | Reverse | `CATAGTGGGGTATCTAATCCTAGTTTG` | Sharks & Rays (Elasmobranchs) |

### COI (Invertebrates) — Highly Degenerate Leray-XT Primers

| Primer Set | Direction | Sequence (5' → 3') | Target Organism |
|:-----------|:----------|:------------------|:----------------|
| **Leray-XT** | Forward | `GGWACWGGWTGAACWGTWTAYCCYCC` | Metazoan Invertebrates |
| **Leray-XT** | Reverse | `TAIACYTCIGGRTGICCRAARAAYCA` | Metazoan Invertebrates |

*Note: W, Y, R, I represent wobble bases to accommodate genetic variability in marine invertebrates.*

### 16S rRNA (Vertebrates & Bacteria) — Universal Suite (6 Primer Pairs)

| Suite Name | Forward Primer (5' → 3') | Reverse Primer (5' → 3') | Specific Target |
|:-----------|:------------------------|:------------------------|:----------------|
| **16S_Vert_Std** | `CGCCTGTTTATCAAAAACAT` | `CCGGTCTGAACTCAGATCACGT` | Standard Vertebrates (Frogs/Reptiles) |
| **16S_Mammal_Mod** | `CGCCTGTTTACCAAAAACAT` | `CCGGTCTGAACTCAGATCACGT` | Mammals & Cetaceans |
| **16S_Mammal_Short** | `CGGTTGGGGTGACCTCGGA` | `GCTGTTATCCCTAGGGTAACT` | Degraded DNA (Short fragments) |
| **16S_Bacteria** | `CCTACGGGNGGCWGCAG` | `GACTACHVGGGTATCTAATCC` | General Bacteria (V3-V4 Region) |
| **16S_NoChloro** | `AACMGGATTAGATACCCKG` | `ACGTCATCCCCACCTTCC` | Bacteria (Excludes Chloroplasts) |

### 18S rRNA (Eukaryotes) — Specialized Suite for Diverse Taxa

| Suite Name | Forward Primer (5' → 3') | Reverse Primer (5' → 3') | Specific Target |
|:-----------|:------------------------|:------------------------|:----------------|
| **18S_Euk_General** | `GTACACACCGCCCGTC` | `TGATCCTTCTGCAGGTTCACCTAC` | Plankton & General Eukaryotes |
| **18S_Fungi_Best** | `CGATAACGAACGAGACCT` | `ANCCATTCAATCGGTANT` | Fungi Kingdom |
| **18S_Nema_Cov** | `GGCAAGTCTGGTGCCAG` | `TCCGTCAATTYCTTTAAGT` | Nematodes (Coverage) |
| **18S_Nema_Eco** | `GGTTAAAAMGYTCGTAGTTG` | `TGGTGGTGCCCTTCCGTCA` | Nematodes (Ecological) |
| **18S_Ciliate** | `GGGRAACTTACCAGGTCC` | `GTGATRWGRTTTACTTRT` | Protozoa (Rumen Ciliates) |

---

## 🔄 Workflow Mechanics: A Tiered Approach to Sequence Recovery

BlueFinDB executes a linear, state-aware pipeline that transitions automatically between processing tiers. Data flows from direct mining (Tier 1) through genomic rescue (Tier 2) to quality assurance and database compilation.

### Tier 1: Hybrid Mining & In-Silico PCR

**Smart Caching:** The system checks the local `output/raw_[MARKER]` directory first. Previously downloaded sequences are loaded instantly, skipping network calls. Novel sequences are retrieved from NCBI.

**In-Silico PCR Simulation:** The core engine simulates laboratory PCR using Cutadapt with configurable parameters:
- `-e 0.20`: Allows 20% error rate (mismatches) to account for biological variation
- `--discard-untrimmed`: Removes sequences where primers fail to bind

**Critical Length Validation:** After trimming, amplicons must fall within marker-specific biological ranges:

| Marker | Allowed Range | Rationale |
|:-------|:--------------|:----------|
| **12S** | 150–550 bp | Excludes short fragments; fits standard MiFish amplicons |
| **COI** | 100–800 bp | Accommodates variable invertebrate barcodes |
| **16S** | 100–600 bp | Captures both short mammal fragments and longer vertebrate sequences |
| **18S** | 100–1000 bp | Allows for variable regions in diverse eukaryotes |

Sequences that fail length validation are discarded. If NO sequence passes this filter, the species is flagged as missing and escalates to Tier 2.

### Tier 2: Deep Genomic Mining (The Fallback)

**Activation:** Triggered immediately when Tier 1 fails to retrieve a valid, trimmable sequence.

**Strategy Shift:** Rather than searching for isolated marker entries, the system now targets larger genomic assemblies:
- For 12S/16S/COI: Downloads Complete Mitochondrial genomes
- For 18S: Queries Whole Genome Sequences (WGS)

**Computational Extraction:** The engine performs in-silico PCR on these assemblies, mathematically "slicing out" the target amplicon. This approach recovers data even when specific marker genes were never deposited to NCBI as standalone entries.

### The Rescue System (Interactive & Automated)

**Trigger:** Activates after Tier 1 and Tier 2 complete. If gaps remain, rescue begins automatically.

**Interactive Mode (Default):**
```bash
bash Bluefindb_run.sh --12S
```
The pipeline **PAUSES** and prompts you to intervene manually. You can edit `fix_taxonomy.sh` to provide missing accession IDs or taxonomy, then press **[Enter]** to resume.

**Non-Interactive Mode (Servers/Batch Jobs):**
```bash
bash Bluefindb_run.sh --12S --non-interactive
```
The pipeline executes the rescue script automatically without pausing, ideal for high-throughput computing environments.

---

## 🍄 Special Feature: The 18S Hybrid Engine

When invoked with `--18S`, BlueFinDB activates a specialized dual-source architecture for eukaryotic complexity. Unlike other markers, 18S implements a strict **Local-First** priority system:

**Step 1: Database Validation**
- The script checks for SILVA 138.1 SSURef NR99 in `Data/databases/`
- If absent, BlueFinDB automatically downloads and decompresses the database

**Step 2: Priority 1 — Local SILVA Mining**
- Uses optimized `awk` commands to scan SILVA headers without loading the 2GB+ file into RAM
- Sanitizes matches by removing alignment gaps and converting RNA bases (U) to DNA bases (T)

**Step 3: Priority 2 — NCBI Fallback**
- Queries NCBI Nucleotide database only if the species is not found locally

**Step 4: "Dragnet" Primer Sweeping**
- Tests all 18S primers against the retrieved sequence
- Selects the amplicon yielding the longest valid product within 100–1000 bp

**Performance Advantage:** Local resolution is ~50× faster than NCBI API querying, dramatically reducing overhead while maintaining taxonomic accuracy through SILVA's professional curation.

---

## ✅ Quality Control & Final Build (Bluefindb_qc.sh)

The final stage aggregates successful sequences from all tiers:

1. **Filter:** Removes sequences with >1 ambiguous base (N) using `vsearch --fastx_filter`
2. **Dereplicate:** Merges 100% identical sequences using `vsearch --derep_fulllength`
3. **Chimera Check:** Removes PCR artifacts using `vsearch --uchime_denovo`
4. **Database Build:** Generates the final local BLAST database via `makeblastdb`

---

## 🧪 Validation & Testing (Sandbox)

BlueFinDB includes a self-contained testing sandbox to verify your installation:

### Prerequisites

Set your credentials before testing:

```bash
export BFD_EMAIL="your.email@university.edu"
export NCBI_API_KEY="your_api_key_here"
```

### Running Tests

Navigate to the test directory and execute with the appropriate marker flag:

```bash
cd test
bash Bluefindb_test_run.sh --12S
```

The `test/Data/` folder contains pre-configured species lists for validation across all four markers.

---

## 🖼️ Workflow Diagram

BlueFinDB executes through a series of automated stages, each designed to maximize sequence recovery while maintaining quality standards.

![BlueFinDB Workflow Architecture](images/workflow_diagram.png)

## 📦 Output Files Reference

| File Name | Format | Purpose |
|:----------|:-------|:--------|
| `BlueFinDB_[MARKER].fasta` | FASTA | Final curated database with enriched taxonomic headers (TaxID, Class, Order, Family) |
| `BlueFinDB_[MARKER].n*` | BLAST Index | Binary index set (.nhr, .nin, .nsq) for immediate BLAST searchability |
| `BlueFinDB_Manifest.txt` | Text | Audit trail: run date, version, primers used, sequence counts, parameters |
| `[MARKER]_accessions.tsv` | TSV | Status report: species, Accession ID, TaxID, Pass/Fail status |
| `fix_taxonomy.sh` | Shell Script | Auto-generated rescue module (created only if needed) |
| `output/raw_[MARKER]/` | Directory | Cache: raw, untrimmed sequence downloads |
| `output/curated_[MARKER]/` | Directory | Trimmed amplicons: individual species FASTA files |
| `output/logs/` | Directory | Execution logs: `bluefindb_run.log`, `cutadapt_details.log`, `final_qc_and_build.log` |

---

## 📜 Citation

If you use BlueFinDB in your research, please cite:

> Vijayakumar, S. (2025). BlueFinDB: A Universal, Hybrid-Mining Reference Database Builder for eDNA Metabarcoding (Version 3.0) [Computer software]. GitHub. https://github.com/SUBRAMANIAM96/BlueFinDB

---

## 📄 License

BlueFinDB is open-source and available under the MIT License. You are free to copy, modify, and distribute this software, provided the original author is credited. See the LICENSE file for details.

---

## 👨‍💻 Author & Contact

**Subramaniam Vijayakumar**  
Lead Developer & Bioinformatics Architect  
Tool Version: v3.0  
Contact: [vijayakumar.subraman@mail.ucv.es](mailto:vijayakumar.subraman@mail.ucv.es)

---

**Last Updated:** December 2025  
**Repository:** [GitHub - BlueFinDB](https://github.com/SUBRAMANIAM96/BlueFinDB)