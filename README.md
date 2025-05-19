## Description

This repository describes a workflow to harvest and process partial and complete mitochondrial genomes from off-target reads of UCE captured data using MitoFinder (Allio et al. 2020) and MAFFT (Katoh & Standley, 2013). Associated data and publication are outlined in the Publication & Data section at the end of this repository.
<br>

## Getting started

### Softwares and Data Requirements

* MitoFinder #to extract mitochondrial genes from cleaned reads - see details on https://github.com/RemiAllio/MitoFinder
* MAFFT (v7.505) #to generate alignments for e.g. phylogenetic inference
* IQTREE #to conduct phylogenetic inferences
<br><br>
* reference_file.gb #containing mitochondrial genomes of reference - can be compiled from NCBI
* left_reads.fastq.gz #containing the left reads of paired-end sequencing
* right_reads.fastq.gz #containing the right reads of paired-end sequencing
* or single-end reads

## Running

Before running MitoFinder, you want to make sure that the nomenclature of mitochondrial genes is consistent across reference genome annotations. For example, the control region gene can also be named d-loop or the COXII / COX2. Pick one and standardize them across using the sed command.

```
sed -i 's/control\ region/D-loop/g' /path/to/mtDNA_reference_genomes.gb
# always run sed without the infile flag (-i) first to make sure that it looks good
```

### Extract mitochondrial data with MitoFinder

Here UCE captured data was used and the default value for --min-contig-size being 1,000 bp is too high for the shorter off-target sequences so it needs to be adapted according to contig coverage. Same applies for --blast-size, which was reduced to e-value = 1e-06 (default 0.00001, 30%) following Allio et al (2019). This will allow you to get more hits and detect more genes. You can do multiple test runs with different values for these options to check how the output assemblies look like. Generally, the lower the contig size the more hits you might get but they are likely to yield more fragmented contigs so there is bit of a trade-off here depending on what data you are going for.

```
#run MitoFinder - there is also a pbs script associated with this repository (run_MitoFinder.pbs)
to run MF as an array job on an HPC cluster and speed things up
mitofinder --metaspades \ #you can also use megahit which is faster but I found that MetaSPAdes gives better assemblies
    -j sample_name \
    -1 paired_end_R1 \
    -2 paired_end_R2 \
    -r mtDNA_reference_genomes.gb \
    -o 2 \ #organism type to specify genetic code, here for vertebrates
    -p 20 \
    -t mitfi \ # for tRNA annotation
    -e 0.000001 \ #this parameter and the following --min-contig-size is what worked best for my UCE captured data
    --min-contig-size 500 \
    --max-memory 100 \
    --new-genes
```
If the mitochondrial genome of your target organism contains introns and/or you want to look into nuclear mitochondrial DNA (NUMT), and have reference genomes from closely related species you can do a second of MitoFinder on the assemblies (--assembly) from the first run searching for introns (--allow-intron), specifying size (--intron-size) and preventing merging of exons (--cds-merge).

<br>

### Tidy up and plot away

When screening lots of samples for mitochondrial contigs, it might be useful to do some tidying up and take advantage of this to generate a summary log file for the results

```
# Make a new directory to store mt contigs
mkdir -p path/to/mt_contigs_MFRUN1

# go to MitoFinder output result folder - it should be folders named after target species 
cd /path/to/MitoFinder/output
for dir in */*_MitoFinder_mitfi_Final_Results; do
    sp=$(echo $dir | cut -d/ -f1) #store species name
# this loop first checks if any mt contigs were found
    if [ -z "$(ls -A "$dir")" ]; then
        echo "$sp: MitoFinder didn't find any mt contigs." >> ../MF_results.log
    else
        # if there are mt contigs for the species then let's create a directory for it
        mkdir -p path/to/mt_contigs_MFRUN1/$sp
        outdir="path/to/mt_contigs_MFRUN1/$sp"
        # Now we want to know if the retreived assemblies are fragmented or not. Up to you to keep one or the other or both
        if [ -f $dir/${sp}_mtDNA_contig.fasta ]; then
                echo "$sp: potentially complete mitchondrial genome found by MitoFinder." >> ../MF_results.log
                cat $sp/${sp}_MitoFinder_mitfi_Final_Results/${sp}_mtDNA_contig.fasta >> $outdir/${sp}_MFRUN1_mtDNA_contigs.fasta
            else
                echo "$sp: fragmented or partial mt genome found by MitoFinder." >> ../MF_results.log
                # concatenate all the mt contigs together. You can comment this part of if you are not interested in fragmented assemblies.
                cat $sp/${sp}_MitoFinder_mitfi_Final_Results/${sp}_mtDNA_contig_[1-9].fasta >> $outdir/${sp}_MFRUN2_mtDNA_contigs.fasta
        fi
    fi
done

```

With prior knowledge on the expected length of your species mitochondrial genomes, you can estimate approximately how many potential whole genomes were retrieved by roughly estimating the length of the assembly from the new output files we just created

```
for contigs in mt_contigs_MFRUN1/*/*.fasta; do
    sp=$(basename "$contigs" _MFRUN1_mtDNA_contigs.fasta) # store species names
    length=$(awk '/^>/ {next} {letters += length} END {print letters}' "$contigs") # calculate approx. length by counting number of base pairs
    echo -e "$sp\t$length"
done | sort -k2,2nr # sort by decreasing order of length
```

Maybe you want to plot the results of your MitoFinder run(s) so let's generate a csv file to input in R.

```
MF_output="/path/to/MitoFinder/output
csv_output="total_MF_genes.csv"

#Create header for the output csv file
echo -e "Species\tGenes\tCount" > "$csv_output"

for fasta in "$MF_output"/*/*_MitoFinder*_Final_Results/*_final_genes_NT.fasta; do
    sp=$(basename "${fasta%_final_genes_NT.fasta}") #store species name
# this loop checks if the genes listed in the text file gene_list.txt are present in the final genes found by MitoFinder. It gives it a "1" for yes and a "0" for no to create a presence/absence matrix
    while read -r line; do
        if grep -qwF "$line" "$fasta"; then
            present=1
        else
            present=0
        fi
        echo -e "$sp\t$line\t$present" >> "$csv_output"
    done < "genes_list.txt"
done
```

Next step is if you want to generate alignments to infer a phylogenetic tree for example. To do so, we need to extract the different mt genes for each species and store them fasta files that we can then align. The script I used is pretty old and to be honest, quite horrendous, but it works fine and eventually I will spend time on it to clean it up. Until then, it probably is best to run it as a bash script (see extract_mt_genes.sh in repository material).

A few things you'll need to change in the extract_mt_gene.sh script to fit your directory setup

```
# Define the root directory containing the species subdirectories (line 9)
root_directory="/path/to/MitoFinder/output"

# Define the output root directory (line 12)
output_root_directory="/path/to/extracted_mt_genes"

# Define taxa name - something that makes sense to your data e.g. if you are doing multiple species from one family you can use the family name (line 21)
taxa="dragons"

# If you want to subset for specific taxon afterwards, uncomment the lines 106-114 to allow the taxa filtering to run (it is turned off by default). You also need to specify a species list and an output directory (lines 24-25)
taxonset="/path/to/species_list.txt"
filtered_output="/path/to/filtered_taxa_output"

# Once you've adapted the script to your setup, make sure it is executable before running it
chmod +x extract_mt_gene.sh # add execute permission to script
./extract_mt_genes.sh --quiet # run script

```

Output should contain
  * 15 directories named ${taxa}_${gene} containing subdirectories with sequences for each species
  * 15 fasta files named all_${taxa}_${gene}.fasta containing sequences of all species for each gene separately
  * 1 log file named ${taxa}_extract_mt_genes.log

<br>

### Align mitochondrial genes with MAFFT

Now that we have concatenated all sequences into fasta files for individual mt genes, we can finally align them! This is an example of an easy loop that will generate alignments for each genes using the MAFFT aligner with default parameters. Best is probably to submit it as a job on a HPC cluster as MAFFT can take a while to run.

```
# load mafft aligner
module load mafft

# create output directory
mafft_out="/path/to/mafft_out"
mkdir -p $mafft_out

for i in /path/to/all_${taxa}_${gene}.fasta; do
    gene=$(basename "$i" .fasta | cut -d_ -f3) # store gene name
    # align mt genes using mafft aligner
    mafft --auto --thread 10 "$i" > ${mafft_out}/${taxa}_${gene}_mafft_alignment.fasta
    # clean loci name from alignments
    sed -i 's/@.*$//' ${mafft_out}/${taxa}_${gene}_mafft_alignment.fasta
done
```

That's it! After trimming your alignments with programs like trimAl or clipKIT you can now use them in downstream analyses e.g. phylogenetic inference
Best fishes and happy coding!
<br><br>
## Help

If you have any issues running the command lines or scripts in this repo feel free to reach out! - lauriane.baraf@my.jcu.edu.au

<br>

## Acknowledgements & Citations

Scripts were written by Lauriane Baraf, please cite this repo if using them.

If using MitoFinder, please cite:<br>
Allio, R., Schomaker‐Bastos, A., Romiguier, J., Prosdocimi, F., Nabholz, B., & Delsuc, F. (2020). MitoFinder: Efficient automated large‐scale extraction of mitogenomic data in target enrichment phylogenomics. Molecular Ecology Resources, 20(4), 892-905.

If using MAFFT v7, please cite:<br>
Katoh, K., & Standley, D. M. (2013). MAFFT multiple sequence alignment software version 7: improvements in performance and usability. Molecular biology and evolution, 30(4), 772-780.

<br>

## Publication & Data

This repository contains the code initially written to extract mitochondrial data from off-target reads for the reef-associated fish family Pomacanthidae (marine angelfishes). Results are presented in the following publication:<br>
Baraf et al. (2024). Comparative mitogenomics of marine angelfishes (F: Pomacanthidae). Ecology and Evolution, 14(8), e70127.

Associated complete mitogenomes for pomacanthid sepcimens are available on GenBank under the accession numbers PP316124-PP316129.


