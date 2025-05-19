#!/bin/bash

# this is just to create a quiet option build into the script
if [[ "$1" == "--quiet" ]]; then
    exec > /dev/null 2>&1
fi

# Define the root directory containing the species subdirectories
root_directory="/path/to/MitoFinder_output"

# Define the output root directory
output_root_directory="/path/to/output_extracted_mt_genes"

# Define the list of genes to extract
genes_to_extract=(
    "ATP8" "ATP6" "COX1" "COX2" "COX3" "CYTB" "ND1" "ND2" "ND3" "ND4"
    "ND4L" "ND5" "ND6" "rrnS" "rrnL"
)

# Define taxa
taxa="chaets"

# If you want to filter by taxon afterwards
taxonset="/path/to/species_list.conf"
filtered_output="/path/to/filtered_taxa_output"

# Function to extract gene sequences from a single fasta file
extract_gene_sequences() {
    local fasta_file="$1"
    # Extract species name using basename and dirname (more robust)
    local species_dir=$(dirname "$fasta_file")
    local species_name=$(basename "$species_dir")
    local current_gene=""
    local current_header=""
    local current_sequence=""

    # Check if the file is a final_genes_NT.fasta file
    if [[ "$fasta_file" != *"_final_genes_NT.fasta" ]]; then
        echo "Skipping non-final_genes_NT.fasta file: $fasta_file"
        return
    fi

    echo "Processing: $fasta_file" >> ${taxa}_extract_mt_genes.log
    echo "Species name: $species_name" >> ${taxa}_extract_mt_genes.log

    # Read the fasta file line by line
    while IFS= read -r line; do
        if [[ "$line" == ">"* ]]; then
            # Extract gene name from header
            current_gene=$(echo "$line" | sed 's/^.*@//')
            current_header="$line"
            current_sequence=""
            echo "Current gene: $current_gene"
            echo "Current header: $current_header"
        else
            # Append sequence data
            current_sequence+="$line"
        fi

        # Check if we have a complete entry (header + sequence)
        if [[ -n "$current_header" && -n "$current_sequence" ]]; then
            echo "Current sequence: $current_sequence"
            # Check if the current gene is in the list of genes to extract
            for gene in "${genes_to_extract[@]}"; do
                if [[ "$current_gene" == "$gene" ]]; then
                    # Create output directory if it doesn't exist
                    output_dir="$output_root_directory/${taxa}_${gene}"
                    mkdir -p "$output_dir"

                    # Construct the new output file name
                    # Extract the first two parts of the species name separated by _
                    species_prefix=$(echo "$species_name" | cut -d'_' -f1,2)
                    
                    # Create output file name
                    output_filename="${species_prefix}_${gene}.fasta"
                    output_filepath="$output_dir/$output_filename"

                    # Write the header and sequence to the output file
                    echo "$current_header" > "$output_filepath"
                    echo "$current_sequence" >> "$output_filepath"

                    echo "Extracted $gene sequence for $species_name and saved to $output_filename in $output_dir" >> ${taxa}_extract_mt_genes.log
                    

                    # Append to the all files
                    all_files="$output_root_directory/all_${taxa}_${gene}.fasta"
                    echo "$current_header" >> "$all_files"
                    echo "$current_sequence" >> "$all_files"
                    echo "Appended $species_name $gene to $all_files"

                    break
                fi
            done
            # Reset variables for the next entry
            current_header=""
            current_sequence=""
        fi
    done < "$fasta_file"
}

# Create the gene folders
for gene in "${genes_to_extract[@]}"; do
    mkdir -p "$output_root_directory/${taxa}_${gene}"
done

# Main loop to process all FASTA files in the directory
find "$root_directory" -type f -name "*_final_genes_NT.fasta" -print0 | while IFS= read -r -d $'\0' fasta_file; do
    extract_gene_sequences "$fasta_file"
done

echo "------------- Finished processing all FASTA files. -------------"

########################################## FILTER TAXON SPECIFIC ################################################
# Make a directory to store your filtered fasta files
# mkdir -p "$output_root_directory"/filtered

# filter taxa - here I'm using taxa from the phyluce taxonset configuration file so it matches 
# for i in "$output_root_directory"/${taxa}_*;
#   do
#   gene=$(basename "$i" | cut -d_ -f2)
#   grep -h -A1 --no-group-separator -F -f $taxonset "$output_root_directory"/${taxa}_"$gene"/*.fasta >> "$output_root_directory"/"$filtered_output"_"$gene".fasta
#   cat "$output_root_directory"/outgroups/outgroups_"$gene".fasta >> "$output_root_directory"/"$filtered_output"_"$gene".fasta
# done
