# 🐟 BlueFinDB v3.0: High-Confidence Reference Database Creator

## Overview

BlueFinDB is an automated bioinformatics pipeline designed to construct taxonomically curated reference databases for environmental DNA (eDNA) metabarcoding. By accepting a target species list, it generates a validated, locally searchable BLAST database optimized for specific marker genes—COI, 12S, 16S, and 18S.

1. BlueFinDB is an automated bioinformatics pipeline designed to construct taxonomically curated reference databases for environmental DNA (eDNA) metabarcoding. It accepts a target species list and generates a validated, local BLAST database for specific marker genes, including COI, 12S, 16S, and 18S.

2. The reproducibility of environmental DNA (eDNA) metabarcoding depends entirely on transparent and replicable workflows for creating complete, curated reference databases. Conventional curation tools (such as RESCRIPT, CRUX, and OBITools) often operate as rigid, non-flexible workflows. They frequently discard valid taxa due to minor primer mismatches or fail completely when specific marker sequences are absent in public repositories.

3. BlueFinDB addresses this gap with an open-source, multi-marker pipeline designed for flexibility and high recovery rates. It builds highly accurate, project-specific reference libraries for 12S, COI, 16S, and 18S markers. By utilizing standard primer sets (e.g., MiFish, Leray-XT, Vences V4, and Stoeck V9) according to researcher preference, it ensures that valid data is retained and curated rather than discarded.

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
├── 📁 images/                   # Documentation assets
│   └── workflow_diagram.png     # Pipeline architecture diagram
├── ⚖️ LICENSE
└── 📖 README.md
```

---

## ⚙️ Script Architecture & Components

### 1. Bluefindb_run.sh (The Orchestrator)

**Function:** Mission Control & Pipeline Management

This master entry point manages workflow logic without performing bioinformatic analysis directly:

- **Environment Validation**: Checks for all dependencies (cutadapt, vsearch, makeblastdb) and API credentials before execution
- **Dynamic Configuration**: Loads active primer sequences based on user flags (--12S, --COI, --16S, --18S)
- **Flow Control**: Sequences Core Engine execution, monitors exit codes, triggers the Rescue System if gaps are detected, and finally launches the QC module

### 2. Bluefindb_core.py (The Mining Engine)

**Function:** Retrieval, In-Silico PCR, & Hybrid Mining

The computational heart of the pipeline—responsible for retrieval, in-silico PCR, and hybrid mining:

- **Hybrid Retrieval**: Implements local database checking (SILVA), result caching, and NCBI fallback queries
- **Smart Trimming**: Executes Cutadapt with configurable mismatch tolerance (default: 20%) to simulate realistic PCR conditions
- **Primer Dragnets**: For 16S and 18S markers, sweeps through multiple primer suites to maximize recovery
- **Genome Mining**: Downloads full genomes or mitochondria when direct markers are unavailable, then mathematically extracts target regions

### 3. rescue_core.py (The Self-Healing System)

**Function:** Gap Analysis & Patch Generation

Runs after core mining to fix retrieval failures without manual intervention:

- **Failure Forensics**: Analyzes the `accessions.tsv` file to identify specific failure modes (Missing Sequence vs. Missing Taxonomy)
- **Auto-Patching**: Dynamically generates a custom Bash script (`fix_taxonomy.sh`) with precise remediation instructions
- **Taxonomy Injection**: Retrieves complete lineage strings for orphan sequences, ensuring taxonomic completeness

### 4. Bluefindb_qc.sh (The Quality Assurance Module)

**Function:** Filtering, Dereplication, & Database Building

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

This configuration is mandatory for NCBI audit trails and recommended for faster API access.

### Step 4: Verify Your Setup (Optional)

Use the built-in diagnostic tool to validate your environment:

```bash
bash Bluefindb_run.sh --check-status
```

BlueFinDB has a built-in diagnostic tool that automatically validates your connection and credentials before the start of the pipeline. Use this command to check your setup before starting.

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
| **16S_Cetacean** | `CGTTTTTTGGTCGACAGCC` | `CCGGTCTGAACTCAGATCACGT` | Cetaceans (MarVer3) |
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

**Interactive Mode (Default & Recommended):**
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

## ⚙️ Workflow Mechanics: Under the Hood

BlueFinDB executes a linear, state-aware pipeline. It transitions automatically between scripts, passing data from the Mining Phase (Tier 1) to the Rescue Phase (Tier 2) and finally to the Build Phase.

### 1. How to Configure

Run these commands in your terminal before starting the pipeline:

```bash
export BFD_EMAIL="your.email@university.edu"          # REQUIRED: Your email address for NCBI audit trails
export NCBI_API_KEY="abcdef1234567890"                # RECOMMENDED: Your NCBI API Key for faster mining
```

### 2. Verification System (--check-status)

BlueFinDB has a built-in diagnostic tool that automatically validates your connection and credentials before the start of the pipeline. If you wish to check it before starting the pipeline, you can run it using this command:

```bash
bash Bluefindb_run.sh --check-status
```

### 3. Activation & Configuration (Bluefindb_run.sh)

This shell script acts as the pipeline's "Switchboard," handling environment setup before any bioinformatics processing begins.

**General Syntax:**
```bash
bash Bluefindb_run.sh [MARKER_FLAG]
```

**Example: Run 12S pipeline**
```bash
bash Bluefindb_run.sh --12S
```

### 4. Targeted Primer Sets

BlueFinDB does not rely on generic searching. It utilizes specific, publication-validated primer sets to perform precise in-silico PCR. The system distinguishes between Static Markers (12S/COI) which use fixed pairs, and Dynamic Markers (16S/18S) which use sweeping "Suites."

**🐟 1. 12S rRNA (Fish):**

The pipeline deploys the industry-standard MiFish primers.

| Name | Direction | Sequence (5' → 3') | Target |
|:---|:---|:---|:---|
| **MiFish-U** | Forward | `GTCGGTAAAACTCGTGCCAGC` | Bony Fish (Teleosts) |
| **MiFish-U** | Reverse | `CATAGTGGGGTATCTAATCCCAGTTTG` | Bony Fish (Teleosts) |
| **MiFish-E** | Forward | `GTCGGTAAAACTCGTGCCAGC` | Sharks & Rays (Elasmobranchs) |
| **MiFish-E** | Reverse | `CATAGTGGGGTATCTAATCCTAGTTTG` | Sharks & Rays (Elasmobranchs) |

**🦀 2. COI (Invertebrates):**

The pipeline uses the highly degenerate Leray-XT primers. These sequences contain "Wobble Bases" (e.g., W, Y, R, I) to account for the high genetic variability found in marine invertebrates.

| Name | Direction | Sequence (5' → 3') | Target |
|:---|:---|:---|:---|
| **Leray-XT** | Forward | `GGWACWGGWTGAACWGTWTAYCCYCC` | Metazoan Invertebrates |
| **Leray-XT** | Reverse | `TAIACYTCIGGRTGICCRAARAAYCA` | Metazoan Invertebrates |

**🐋 3. 16S rRNA (Vertebrates & Bacteria):**

The pipeline activates a Universal Suite that iteratively tests 6 different primer pairs against every downloaded sequence to maximize recovery.

| Suite Name | Forward Primer (5' → 3') | Reverse Primer (5' → 3') | Specific Target |
|:---|:---|:---|:---|
| **16S_Vert_Std** | `CGCCTGTTTATCAAAAACAT` | `CCGGTCTGAACTCAGATCACGT` | Standard Vertebrates (Frogs/Reptiles) |
| **16S_Mammal_Mod** | `CGCCTGTTTACCAAAAACAT` | `CCGGTCTGAACTCAGATCACGT` | Mammals & Whales |
| **16S_Mammal_Short** | `CGGTTGGGGTGACCTCGGA` | `GCTGTTATCCCTAGGGTAACT` | Degraded DNA (Short fragments) |
| **16S_Cetacean** | `CGTTTTTTGGTCGACAGCC` | `CCGGTCTGAACTCAGATCACGT` | Cetaceans (MarVer3) |
| **16S_Bacteria** | `CCTACGGGNGGCWGCAG` | `GACTACHVGGGTATCTAATCC` | General Bacteria (V3-V4 Region) |
| **16S_NoChloro** | `AACMGGATTAGATACCCKG` | `ACGTCATCCCCACCTTCC` | Bacteria (Excludes Chloroplasts) |

**🍄 4. 18S rRNA (Eukaryotes):**

The 18S module sweeps through a specialized suite designed to capture distinct eukaryotic groups, ensuring that elusive taxa like Fungi and Nematodes are not missed during the in-silico PCR step.

| Suite Name | Forward Primer (5' → 3') | Reverse Primer (5' → 3') | Specific Target |
|:---|:---|:---|:---|
| **18S_Euk_General** | `GTACACACCGCCCGTC` | `TGATCCTTCTGCAGGTTCACCTAC` | Plankton / General Eukaryotes |
| **18S_Fungi_Best** | `CGATAACGAACGAGACCT` | `ANCCATTCAATCGGTANT` | Fungi Kingdom |
| **18S_Nema_Cov** | `GGCAAGTCTGGTGCCAG` | `TCCGTCAATTYCTTTAAGT` | Nematodes (Worms) - Coverage |
| **18S_Nema_Eco** | `GGTTAAAAMGYTCGTAGTTG` | `TGGTGGTGCCCTTCCGTCA` | Nematodes - Ecological |
| **18S_Ciliate** | `GGGRAACTTACCAGGTCC` | `GTGATRWGRTTTACTTRT` | Protozoa (Rumen Ciliates) |

### 5. Tier 1: Hybrid Mining & In-Silico PCR (Bluefindb_core.py)

Control is passed to the Python Core. It executes the Tier-1 Retrieval Protocol.

**Step A: Smart Caching:** It checks the `output/raw_[MARKER]` directory. If a sequence was downloaded in a previous run, it loads the local copy instantly, skipping the network call. If not, it downloads available user-defined species-specific marker sequences from NCBI.

**Step B: In-Silico PCR (Trimming):** The core engine simulates a laboratory PCR reaction using Cutadapt:
- `-e 0.20`: Allows a 20% error rate (mismatches) to account for biological variation
- `--discard-untrimmed`: Strictly removes sequences where primers do not bind

**⛔ Critical Step: Length Validation**

After trimming, the engine measures the amplicon. It must fall within the strict biological ranges defined below. If a sequence is too short (noise) or too long (read-through errors), it is discarded immediately.

| Marker | Allowed Range | Why? |
|:---|:---|:---|
| **12S** | 150–550 bp | Excludes short fragments; fits standard MiFish amplicons |
| **COI** | 100–800 bp | Wide range to accommodate variable invertebrate barcodes |
| **16S** | 100–600 bp | Captures both short mammal (140bp) and long vertebrate fragments |
| **18S** | 100–1000 bp | Allows for variable regions in diverse eukaryotes/fungi |

**Logic:** The engine loops through all available sequences. It discards invalid ones until a suitable amplicon is found. If NO sequence passes this filter, the species is flagged as "Missing," and the system automatically escalates to Tier 2.

### 6. Tier 2: Deep Genomic Mining (The Fallback)

**Trigger:** Activated immediately if Tier 1 fails to retrieve a valid, trimmable sequence.

**Strategy:** The engine shifts strategy to download the organism's larger genetic assemblies:
- For 12S/16S/COI: It downloads the Complete Mitochondrion
- For 18S: It searches for Whole Genome Sequences (WGS)

**Action:** It performs in-silico PCR on these large assemblies to mathematically "slice out" the target amplicon. This allows BlueFinDB to recover data even when the specific marker gene was never uploaded to NCBI as a standalone entry.

### 7. The Rescue System (`rescue_core.py`)

**Trigger:**

This phase begins **after** the Core Engine (Tier 1 & 2) has finished processing the entire species list. If gaps still remain (i.e., both Tier 1 and Tier 2 failed), the pipeline automatically launches the Rescue Script.

#### The Critical Decision: Interactive vs. Non-Interactive

Once rescue targets are identified (or if manual intervention is needed for missing accessions), the pipeline behavior depends entirely on your execution flag:

**A. Interactive Mode (Default & Recommended)**

The pipeline **PAUSES** and alerts you:
> `⚠️ Action Required: Please open fix_taxonomy.sh`

You can open the generated script, manually paste a specific Genome/Mitochondrion Accession ID to "help" the rescue engine, and press **[Enter]** to resume.

```bash
bash Bluefindb_run.sh --12S
```

**B. Non-Interactive Mode (--non-interactive):** 

The pipeline NEVER STOPS. It immediately executes the rescue script to attempt the automated genomic mining. Use this only for high-throughput server jobs.

```bash
bash Bluefindb_run.sh --12S --non-interactive
```

### 8. How to Edit the Rescue Script (fix_taxonomy.sh)

This script is the final safety net of the BlueFinDB pipeline. It is generated automatically by rescue_core.py only when specific species have failed both Tier 1 (Direct Mining) and Tier 2 (Genome Rescue).

**A. Why Am I Seeing This? (The "Last Resort" Principle)**

- If a species appears in this script, it means BlueFinDB has already scoured NCBI for direct marker entries, complete mitochondria, and whole genomes, but found nothing that matches your quality criteria.

- You are now manually intervening to point the tool toward a sequence it might have missed (e.g., a new submission, a synonym, or a raw unannotated file).

**B. How to Edit the File**

The script is custom-built based on the specific failures of your current run. You simply need to "fill in the blanks" for the variables marked `NA_FILL_ME`.

**Step A: Open the File**

When the pipeline pauses in Interactive Mode, open `fix_taxonomy.sh` in any text editor.

**Step B: Find a Valid Accession**

Search NCBI Nucleotide for your species and find a valid Accession ID (e.g., AF518210 or NC_001234) that contains the target gene.

**Step C: Fill in the Blank**

Paste the ID into the variable:

```bash
# --- BEFORE EDITING (What the script generates) ---
# 🐟 FIX REQUIRED: Sebastes_marinus (Reason: FAILED_ALL_PRIMERS)
CORRECT_ACCESSION_SEBASTES_MARINUS="NA_FILL_ME__SEBASTES_MARINUS_ACC"

# --- AFTER EDITING (What you save) ---
# 🐟 FIX REQUIRED: Sebastes_marinus
CORRECT_ACCESSION_SEBASTES_MARINUS="AF518210"
```

**C. Strategic Warning: "Amplicons > Genomes"**

Do not simply provide another Whole Genome or Complete Mitochondrion accession.

**Why?** The automated Tier-2 system has likely already attempted to mine these assemblies and failed due to assembly gaps or primer mismatches in that specific region.

**The Winning Strategy:** To guarantee a successful rescue, search NCBI for a specific, high-quality amplicon (e.g., a direct PCR submission or partial CDS) that strictly adheres to the Marker-Specific Target Range defined above (e.g., 150-550bp for 12S). This removes alignment ambiguity and maximizes the success rate of the in-silico trimming.

**D. Patching Missing Taxonomy**

Sometimes a sequence is valid, but NCBI returns incomplete lineage (e.g., missing Class or Order).

**Your Task:** Locate the specific rank marked `NA_FILL_ME` and type the correct scientific name.

```bash
# --- BEFORE EDITING ---
# Taxonomy for: Dicentrarchus_labrax
CORRECT_TAXID_DICENTRARCHUS_LABRAX="13489"
CORRECT_ORDER_DICENTRARCHUS_LABRAX="NA_FILL_ME__Dicentrarchus_labrax_ORDER"

# --- AFTER EDITING ---
# Taxonomy for: Dicentrarchus_labrax
CORRECT_TAXID_DICENTRARCHUS_LABRAX="13489"
CORRECT_ORDER_DICENTRARCHUS_LABRAX="Moroniformes"
```

**E. Execution**

Once you have saved your changes:

1. Close the text editor
2. Return to your terminal
3. Press **[Enter]**
4. BlueFinDB will execute your patched script, download the specific accessions, trim them, and inject them into the final database

### 9. 🍄 Special Feature: The 18S Hybrid Engine

When running with the `--18S` flag, BlueFinDB activates a specialized Dual-Source Architecture designed specifically for the complexities of eukaryotic taxonomy. Unlike the other markers (which rely primarily on NCBI), the 18S module implements a strict priority system to ensure speed and taxonomic fidelity.

**A. How It Works: The Hybrid Protocol**

The pipeline executes a "Local-First" search strategy for every species in your list.

**Step 1: Database Validation**

- **Logic:** Before processing, the script checks for the local presence of the SILVA 138.1 SSURef NR99 database in `Data/databases/`.

- **Auto-Install:** If the database is missing, BlueFinDB automatically downloads (wget) and decompresses (gunzip) the file.

**B. Step 2: Priority 1 - Local SILVA Mining**

- **Command:** awk (Fast Grep) Logic - The engine first searches the local SILVA database.

- It uses optimized awk commands to scan headers for your species name without loading the massive 2GB+ file into RAM.

- **Data Cleaning:** If a match is found, the script automatically sanitizes the sequence:
  - **Gap Removal:** Removes alignment gaps (`.` or `-`)
  - **DNA Conversion:** Converts RNA bases (`U`) to DNA bases (`T`) so primers can bind

**C. Step 3: Priority 2 - NCBI Fallback**

- **Logic:** If—and only if—the species is not found in SILVA, the engine "fails over" to the NCBI Nucleotide database.

**D. Step 4: "Dragnet" Primer Sweeping**

- **Logic:** Once a raw sequence is retrieved (from either SILVA or NCBI), the engine runs In-Silico PCR. Because eukaryotes are highly diverse, a single primer pair is insufficient. BlueFinDB sweeps through the 18S Primer Suite:

- **Selection Algorithm:** The engine tests all primers against the sequence. It selects the result that yields the longest valid amplicon within the target range (100–1000 bp).

**E. Why This Matters?**

- **Speed:** Resolving species locally via SILVA is ~50x faster than querying the NCBI API.

- **Accuracy:** SILVA is professionally curated. Prioritizing it ensures your database uses high-quality reference taxonomy before falling back to the "noisy" public data of NCBI.

### 10. QC & Final Build (Bluefindb_qc.sh)

The final stage aggregates all successful sequences from Tier 1, Tier 2, and the Rescue Phase:

- **Filter:** Removes sequences with >1 'N' (vsearch --fastx_filter)
- **Dereplicate:** Merges 100% identical duplicates (vsearch --derep_fulllength)
- **Chimera Check:** Removes PCR artifacts (vsearch --uchime_denovo)
- **Build:** Creates the final local BLAST database

---

## 🧪 Validation & Testing (Sandbox)

BlueFinDB includes a self-contained sandbox in the test/ directory to verify your installation before running real jobs.

### Prerequisites

Ensure your email and API key are set in your terminal, otherwise the test will fail:

```bash
export BFD_EMAIL="your.email@university.edu"
export NCBI_API_KEY="your_api_key_here"
```

### Verify Test Data

The `test/Data/` folder is already populated with sample species lists for verification:

- `species_list_12s.txt`
- `species_list_16S.txt`
- `species_list_18S.txt`
- `species_list_COI.txt`

### Configure the Test Runner

Just like the main pipeline, you must tell the script which file to process:

1. Open `test/Bluefindb_test_run.sh` in a text editor
2. Locate the `DEFAULT_INPUT_FILE` variable
3. Update it to match the marker you want to test

```bash
# Inside Bluefindb_test_run.sh
DEFAULT_INPUT_FILE="Data/species_list_12s.txt"
```

### Run the Test

Navigate to the test folder and execute the script using the corresponding marker flag:

```bash
cd test
bash Bluefindb_test_run.sh --12S
```

---

## ✅ Quality Control & Final Build (Bluefindb_qc.sh)

The final stage aggregates successful sequences from all tiers:

1. **Filter:** Removes sequences with >1 ambiguous base (N) using `vsearch --fastx_filter`
2. **Dereplicate:** Merges 100% identical sequences using `vsearch --derep_fulllength`
3. **Chimera Check:** Removes PCR artifacts using `vsearch --uchime_denovo`
4. **Database Build:** Generates the final local BLAST database via `makeblastdb`

---

## 📦 Output Files Reference

| File Name | Format / Type | Information Provided |
|:---|:---|:---|
| `BlueFinDB_[MARKER].fasta` | FASTA | The **final curated database file**. Contains optimized gene fragments (e.g., 12S, COI) with fully enriched, validated taxonomic headers (TaxID, Class, Order, Family, etc.). |
| `BlueFinDB_[MARKER].n*` | BLAST Index Files | The BLAST binary index set (including `.nhr`, `.nin`, `.nsq`). These files are generated by `makeblastdb` and make your FASTA file instantly searchable for local BLAST alignment. |
| `BlueFinDB_Manifest.txt` | Text File | Audit trail including run date, BlueFinDB version, active primer sequences, total sequence counts, and processing metadata. |
| `[MARKER]_accessions.tsv` | TSV | Status report listing each target species, retrieved Accession ID, NCBI TaxID, and the final Pass/Fail retrieval status. |
| `fix_taxonomy.sh` | Shell Script | Auto-generated rescue module created only if missing taxonomy or failed species are detected. Contains logic to mine genomes or inject missing metadata. |
| `output/raw_[MARKER]/` | Directory | Cache storage containing raw, untrimmed sequence downloads (NCBI or SILVA). Enables fast re-runs without redownloading. |
| `output/curated_[MARKER]/` | Directory | Trimmed amplicons. Holds individual FASTA files for each species, cut precisely to the target marker length. |
| `output/logs/` | Directory | Detailed execution logs, including `bluefindb_run.log` (master log), `cutadapt_details.log` (trimming stats), and `final_qc_and_build.log` (VSEARCH summary). |

---

## 🖼️ Workflow Diagram

BlueFinDB executes through a series of automated stages, each designed to maximize sequence recovery while maintaining rigorous quality standards. The diagram below illustrates the complete pipeline architecture, including multi-tier recovery mechanisms and quality assurance checkpoints:

![BlueFinDB Workflow Architecture](images/workflow_diagram.png)

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
