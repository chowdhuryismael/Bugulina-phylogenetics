# Bugulina-phylogenetics
rom Bio import Entrez
Entrez.email = "iac63@humboldt.edu"  # Required by NCBI

# Search for Bugulidae COI sequences
handle = Entrez.esearch(db="nucleotide", term="txid10210[Organism] AND COI[Gene]", retmax=1000)
record = Entrez.read(handle)
ids = record["IdList"]
print(f"Found {len(ids)} sequences.")

# Download in GenBank format (with metadata)
handle = Entrez.efetch(db="nucleotide", id=ids, rettype="gb", retmode="text")
with open("bugulidae_ncbi.gb", "w") as f:
    f.write(handle.read())
print("Saved to bugulidae_ncbi.gb")


#!/bin/bash



# Usage: ./clean_one.sh GenusName



GENUS=$1

INPUT="${GENUS}/${GENUS}_sequence.fasta"

OUTPUT="${GENUS}/${GENUS}_sequence_cleaned.fasta"



echo "Cleaning $INPUT..."



# Extract sequences with 'COI' or 'ribosomal' in header, keep accession + species name, replace ambiguous bases

awk '/^>/ {header=$0; getline seq; if (header ~ /COI|ribosomal/) print header "\n" seq}' "$INPUT" | \

sed -E '/^>/ {

  s/^>([^ ]+) ([A-Za-z]+ [a-z]+).*/>\1_\2/

  b

}

s/[RYKWSMBDHV]/N/g

' > "$OUTPUT"

echo "Cleaned COI and ribosomal sequences saved to $OUTPUT"


!/bin/bash



# Check if the correct number of arguments is provided

if [ "$#" -ne 1 ]; then

    echo "Usage: $0 name.fa"

    exit 1

fi



input_fasta=$1

base_name=$(basename "$input_fasta" .fa)



mkdir -p outputs



# MAFFT alignment

if [ -f "outputs/mafft_${base_name}.fa" ]; then

    echo "MAFFT alignment already run. Results found in outputs/mafft_${base_name}.fa. Skipping Analysis"

else

    echo "Running MAFFT alignment for $input_fasta"

    mafft --auto "$input_fasta" > "outputs/mafft_${base_name}.fa"

fi



# MUSCLE alignment

if [ -f "outputs/muscle_${base_name}.fa" ]; then

    echo "MUSCLE alignment already run. Results found in outputs/muscle_${base_name}.fa. Skipping Analysis"

else

    echo "Running MUSCLE alignment for $input_fasta"

    muscle -align "$input_fasta" -output "outputs/muscle_${base_name}.fa"

fi
# Clustal Omega alignment

if [ -f "outputs/clustalo_${base_name}.fa" ]; then

    echo "CLUSTAL-omega alignment already run. Results found in outputs/clustalo_${base_name}.fa. Skipping Analysis"

else

    echo "Running CLUSTAL-omega alignment for $input_fasta"

    clustalo --force -i "$input_fasta" -o "outputs/clustalo_${base_name}.fa"

fi



# Trimming with trimal

if [ -f "outputs/trimal_mafft_${base_name}.fa" ] && \

   [ -f "outputs/trimal_muscle_${base_name}.fa" ] && \

   [ -f "outputs/trimal_clustalo_${base_name}.fa" ]; then

    echo "Trimming of sequences already done. Found in outputs/"

else

    echo "Running trimal for all alignments"

    trimal -in "outputs/mafft_${base_name}.fa" -out "outputs/trimal_mafft_${base_name}.fa" -gt 0.8 -st 0.01

    trimal -in "outputs/muscle_${base_name}.fa" -out "outputs/trimal_muscle_${base_name}.fa" -gt 0.8 -st 0.01

    trimal -in "outputs/clustalo_${base_name}.fa" -out "outputs/trimal_clustalo_${base_name}.fa" -gt 0.8 -st 0.01

fi



# Get the length of the first sequence in each trimmed file

get_first_seq_length() {

    awk '/^>/ {if (seq) {print length(seq); exit} else {seq=""}} !/^>/ {seq=seq$0} END {print length(seq)}' "$1"

}
mafft_length=$(get_first_seq_length "outputs/trimal_mafft_${base_name}.fa")

muscle_length=$(get_first_seq_length "outputs/trimal_muscle_${base_name}.fa")

clustalo_length=$(get_first_seq_length "outputs/trimal_clustalo_${base_name}.fa")



# Determine the file with the longest first sequence

longest_file="outputs/trimal_mafft_${base_name}.fa"

longest_length=$mafft_length



if [ "$muscle_length" -gt "$longest_length" ]; then

    longest_file="outputs/trimal_muscle_${base_name}.fa"

    longest_length=$muscle_length

fi



if [ "$clustalo_length" -gt "$longest_length" ]; then

    longest_file="outputs/trimal_clustalo_${base_name}.fa"

    longest_length=$clustalo_length

fi



echo "The longest sequence is in $longest_file with length $longest_length"

long_basename=$(basename "$longest_file" .fa)



# Run IQ-TREE on the longest alignment

if [ -f "outputs/iqtree_${long_basename}.contree" ]; then

    echo "Phylogenetic analysis already done. Found in output folder with 'iqtree_' prefix. Skipping."
else

    echo "Performing phylogenetic analysis"

    iqtree -s "$longest_file" -m GTR+G -B 1000 --bnni -alrt 1000 -T 16 --prefix "outputs/iqtree_${long_basename}"

fi
