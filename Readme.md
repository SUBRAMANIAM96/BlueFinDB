## BlueFinDB – Custom 12S Reference Database Tool

# 🐟 BlueFinDB  

**BlueFinDB** is a lightweight and user-friendly tool for building custom BLAST databases from fish 12S rRNA sequences. Designed with simplicity in mind, it allows researchers, students, and conservationists to generate their own searchable database from a list of species names. With BlueFinDB, you can create a specialized reference for diet analysis, metabarcoding, and biodiversity studies without needing to manually download or curate sequences. Just provide a plain text species list, and BlueFinDB handles the rest.

- Downloads sequences from NCBI
- Selects the **best (longest) sequence per species**
- Merges sequences
- Formats **BLAST-ready headers**
- Produces a database ready for BLAST analyses

---

## 📂 Folder Structure

BlueFinBase/
├── core/
│ └── bluefindb # compiled binary (hidden core script)
├── run_bluefindb.sh # user wrapper
├── test/
│ └── species_list_test.txt # small test species list
├── species_list.txt # example species list
└── README.md # this file


---

## 🔹 Requirements

- Linux, WSL (Windows Subsystem for Linux), or macOS
- **BLAST+** installed (required for local BLAST)
- **EDirect** (esearch/efetch) for NCBI sequence download
- **Python 3** with **Biopython**:

```bash
pip install biopython

## 🔹 How to Make a Species List File

1. Create a plain text file (`.txt`) with each species on a separate line.  

2. include pipe-separated metadata in this order:  

    Genus species|Family|Common_Name

Example (species_list.txt):

Gadus morhua|Gadidae|Atlantic_cod
Merluccius merluccius|Merlucciidae|European_hake
Scomber scombrus|Scombridae|Atlantic_mackerel
Sardina pilchardus|Clupeidae|European_pilchard

## 🔹How to Run BlueFinDB

# Make the wrapper executable

chmod +x run_bluefindb.sh

Make sure your species list text file (e.g., `species_list.txt`) is inside the **BlueFinDB** folder.  
   - Example: `BlueFinDB/species_list.txt`  

Open a terminal and go to the **BlueFinDB** folder:  
cd BlueFinDB

# Run the tool with your species list (from inside BlueFinBase/):

./run_bluefindb.sh species_list.txt

# Outputs will be created inside the same folder:

fish_12s_raw/ → raw sequences downloaded from NCBI

fish_12s_filtered/ → best (longest) sequence per species

fish_12s_merged.fasta → merged FASTA of best sequences

fish_12s_BLAST.fasta → BLAST-ready database

⚡ Use fish_12s_BLAST.fasta for local BLAST or parsing.

##🔹 Using the Test Dataset

# Run the small test species list included in test/:
./run_bluefindb.sh test/species_list_test.txt

##🔹 Tips & Notes

Always run run_bluefindb.sh from the BlueFinBase folder for simplicity.

You can provide any species list file anywhere by giving the full path.

The tool automatically sleeps between NCBI requests to avoid blocking (default 15 seconds).

Including a test/ folder is recommended so new users can quickly verify functionality.

The tool selects one best (longest) 12S sequence per species automatically.

## 👥 Contributors
Subramaniam Vijayakumar 📧 Email: subramanyamvkumar@gmail.com 🔗 GitHub: SUBRAMANIAM96