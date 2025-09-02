#!/bin/bash
# Simple wrapper for BlueFinDB

if [ $# -ne 1 ]; then
    echo "=========================================="
    echo "   BlueFinDB: Custom 12S Reference DB Tool "
    echo "   Developed by Subramnaiam Vijayakumar        "
    echo "=========================================="
    echo "Usage: ./run_bluefindb.sh species_list.txt"
    exit 1
fi

INPUT=$1
OUTPUT="output_db"

echo "🚀 Running BlueFinDB..."
./core/bluefindb "$INPUT" 

echo "✅ Done! Your BLAST database is in: $OUTPUT"
