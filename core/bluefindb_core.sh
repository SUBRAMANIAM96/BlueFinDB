#!/bin/bash

# ================================
# CONFIG (minimal changes to accept CLI args)
# ================================
# Usage: ./bluefindb_core.sh [species_list.txt] [sleep_seconds]
SPECIES_LIST="${1:-fish_species_europe.txt}"  # species list file (or first CLI arg)
RAW_DIR="fish_12s_raw"
FILTERED_DIR="fish_12s_filtered"
MERGED_FASTA="fish_12s_merged.fasta"
BLAST_FASTA="fish_12s_BLAST.fasta"
SLEEP_TIME="${2:-15}"  # seconds between NCBI requests to avoid blocking (can override with 2nd arg)

# Create directories
mkdir -p "$RAW_DIR" "$FILTERED_DIR"

# ----------------------------
# Basic checks
# ----------------------------
if [ ! -f "$SPECIES_LIST" ]; then
    echo "ERROR: species list file not found: $SPECIES_LIST"
    echo "Usage: $0 [species_list.txt] [sleep_seconds]"
    exit 1
fi

# Check required commands
missing=()
for cmd in esearch efetch python3; do
    if ! command -v "$cmd" >/dev/null 2>&1; then
        missing+=("$cmd")
    fi
done

if [ "${#missing[@]}" -ne 0 ]; then
    echo "ERROR: Missing required commands: ${missing[*]}"
    echo "Install NCBI EDirect (esearch/efetch) and python3."
    exit 2
fi

# Check Biopython availability
if ! python3 -c "import Bio" >/dev/null 2>&1; then
    echo "ERROR: Biopython not installed for python3."
    echo "Install with: python3 -m pip install --user biopython"
    exit 3
fi

echo "Using species list: $SPECIES_LIST"
echo "Sleeping $SLEEP_TIME seconds between queries."

# ================================
# Download 12S sequences from NCBI
# ================================
# Read the whole list into an array to avoid stdin conflicts inside the loop
mapfile -t SPECIES_LINES < "$SPECIES_LIST"

for line in "${SPECIES_LINES[@]}"; do
    # Parse fields (expected: species|family|common)
    IFS="|" read -r species family common <<< "$line"

    # Remove possible carriage returns and trim spaces
    species=$(echo "${species:-}" | tr -d '\r' | xargs)
    family=$(echo "${family:-}" | tr -d '\r' | xargs)
    common=$(echo "${common:-}" | tr -d '\r' | xargs)

    # Skip empty or commented lines
    [[ -z "$species" ]] && continue
    case "$species" in \#*) continue ;; esac

    echo "Downloading 12S sequences for $species ..."
    raw_out="$RAW_DIR/${species// /_}.fasta"
    filtered_out="$FILTERED_DIR/${species// /_}_best.fasta"

    # Download all available 12S sequences for species
    # (keep same query you used; use || true so a failed fetch doesn't abort the whole run)
    esearch -db nucleotide -query "$species[ORGN] AND (12S[Title] OR 12S rRNA[Title]) AND mitochondrion[Filter] NOT uncultured[Title] NOT unverified[Title]" \
    | efetch -format fasta > "$raw_out" || true

    # Keep only the longest sequence (best representative)
    python3 << EOF
from Bio import SeqIO
sp_key = "${species}".replace(" ", "_")
records = list(SeqIO.parse("$raw_out", "fasta"))
if records:
    best = max(records, key=lambda r: len(r.seq))
    # Set the FASTA ID to species key so it matches your metadata later
    best.id = sp_key
    best.description = ""
    SeqIO.write(best, "$filtered_out", "fasta")
EOF

    if [[ -s "$filtered_out" ]]; then
        echo "Selected best (longest) 12S sequence for $species"
    else
        echo "⚠️ No valid 12S sequences found for $species"
    fi

    echo "Sleeping for $SLEEP_TIME seconds..."
    sleep "$SLEEP_TIME"
done

# ================================
# Merge all filtered FASTA
# ================================
# Only attempt to cat if there are files to merge
shopt -s nullglob
fasta_files=("$FILTERED_DIR"/*.fasta)
if [ ${#fasta_files[@]} -gt 0 ]; then
    # Overwrite merged fasta (don't append accidentally)
    cat "${fasta_files[@]}" > "$MERGED_FASTA"
    echo "Merged all best sequences into $MERGED_FASTA"
else
    echo "No filtered FASTA files found in $FILTERED_DIR. Nothing to merge."
    echo "✅ Exiting — no BLAST FASTA created."
    exit 0
fi
shopt -u nullglob

# ================================
# Format BLAST-ready headers
# ================================
python3 << EOF
from Bio import SeqIO

species_info = {}
with open("$SPECIES_LIST") as f:
    for line in f:
        parts = [p.strip() for p in line.strip().split("|")]
        if len(parts) != 3:
            continue
        sp, fam, com = parts
        species_info[sp.replace(" ", "_")] = (fam, com)

with open("$MERGED_FASTA") as infile, open("$BLAST_FASTA", "w") as outfile:
    for record in SeqIO.parse(infile, "fasta"):
        sp_key = record.id  # we set this to Genus_species above
        genus = sp_key.split("_")[0]
        fam, com = species_info.get(sp_key, ("NA", "NA"))
        header = f">{sp_key}|SPECIES:{sp_key}|GENUS:{genus}|FAMILY:{fam}|COMMON_NAME:{com}"
        outfile.write(header + "\n" + str(record.seq) + "\n")

print("BLAST-ready FASTA created:", "$BLAST_FASTA")
EOF

echo "✅ Done! Your 12S reference database is ready: $BLAST_FASTA"