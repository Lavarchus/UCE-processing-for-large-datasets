## Description

This repository outlines a workflow for processing UCE from targeted capture data with Phyluce (Faircloth, 2016) for large genomic datasets. Associated data and publication are outlined in the Publication & Data section at the end of this repository.
<br>

## Adapted steps from the original Phyluce pipeline for large genomic datasets

Some steps can take a while to run within Phyluce, especially for large datasets, so I find it's best to go in and out of the original pipeline to process this type of data.
Others steps can be facillitated by running easy shell scripts. This repository only goes through the steps that were adapted to run more efficiently on large datasets and does not go through the whole Phyluce pipeline. Please refer to the original Phyluce tutorial in parallel of running this workflow.

* Assembling contigs
* Structure your scripts and outputs
* Generate configuration file to extract UCE loci
* Prepare data to align UCE loci with MAFFT
* Align UCE loci with MAFFT outside of Phyluce
* Generate final alignment matrices

## Getting started

### Softwares and Data Requirements

* Phyluce v1.7.3 - to process UCE from targeted capture sequencing. See pipeline details on https://phyluce.readthedocs.io/en/latest/tutorials/tutorial-1.html
* SPAdes v3.15 - to assemble contig
* MAFFT v7.505 - to generate alignments for e.g. phylogenetic inference
* IQTREE v2.2 - to conduct phylogenetic inferences
<br><br>
* trimmed reads (left_reads.fastq.gz) - containing the left reads of paired-end sequencing
* trimmed reads (right_reads.fastq.gz) - containing the right reads of paired-end sequencing
* uce probe set (e.g. fish-uce-1k-probes.fasta)
<br>

### Assemble contigs with SPAdes

Before assembling your contigs, reads should have been trimmed from adapters and quality controlled (e.g. using FastQC).
This step will take some time to run, so let's run it as an array job to basically duplicate the job and process multiple samples in parallel.
To create an array job in you need to add the PBS directive -J followed by the number of samples or indices in the script PBS header. It will look something like this:

```
#!/bin/bash #specifies that the script should be executed using the Bash shell
#PBS -c s #if your environment supports job restart from checkpoints
#PBS -j oe #join standard output and standard error into a single file. This is useful for consolidating job logging
#PBS -m bae  #mail event - send me an email when job b(egins), a(lts/crashes), e(xits/finishes)
#PBS -M lauriane.baraf@my.jcu.edu.au #email recipient
#PBS -N name_job
#PBS -l select=1:ncpus=10:mem=20gb #amount of resources requested
#PBS -l walltime=08:00:00 #time allocated to run the job
#PBS -J 0-99
set -e #exit the script immediately if any command returns a non-zero (error) exit code. Useful Bash settings that helps catch failures early

shopt -s expand_aliases #allows use of aliases in the script
source /etc/profile.d/modules.sh #makes environment module system discoverable
```
Note that a -J directive needs to match both the number of samples to run but also how the job scheduler is set up.
The above example is for a PBS script to be submitted on a shared cluster, which in my case, is set to run a maximum of 100 duplicated jobs so a 100 samples. The -J directive here has a range of indices from 0-99 and that's because in the system I am using 0 counts as 1. Different HPC clusters might have different array job parameters.
To check your cluster-specific limits on array jobs you can use command-line tools:

```
# For PBS and TORQUE job schedulers
qmgr -c 'p s' | grep array

# For SLURM
scontrol show config | grep -i array
```

Now that the PBS job header is all set up properly, we can work on the rest of the script to run SPAdes. It will spit out a lot of heavy intermediate files while running so if you can, I recommend using a temporary directory to avoid overloading your storage. Then redirect only the final output to your output directory.

```
# load SPAdes
module load spades/3.15.5

# Get job index from the PBS job array (e.g., if -J 0-99 is used, this will be 0 through 99). 
id=${PBS_ARRAY_INDEX}

# Create an array of all folder names (each corresponding to a sample) in the input directory.
# Each sample directory should contain the trimmed paired-end FASTQ files for that sample.
READS_DIRS=($(ls /path/to/trimmed/reads/directory))

# Get the sample name based on the array index, basically linking the two together.
SAMPLE_NAME=${READS_DIRS[$id]}

# Define paths to paired-end trimmed forward (1) and reverse (2) read files
R1=$(echo "/path/to/trimmed/reads/directory/${SAMPLE_NAME}/${SAMPLE_NAME}-READ1.fastq.gz")
R2=$(echo "/path/to/trimmed/reads/directory/${SAMPLE_NAME}/${SAMPLE_NAME}-READ2.fastq.gz")

# Define directories
OUT_DIR="/path/to/where/to/store/spades_assemblies" #to store final output
TMP_DIR="/path/to/temporary/directory #to store heavy intermediate files

# Run SPAdes
spades -1 ${R1} -2 ${R2} --careful --threads 10 --cov-cutoff 2 --tmp-dir ${TMP_DIR} -o ${OUT_DIR}/${SAMPLE_NAME}
```

Once the job has finished, let's make sure that SPAdes actually successfully ran for all samples. Even if you get an Exit=0 job, it's good to check that there are not some cryptic run errors, especially when you have many samples. Thankfully, SPAdes is a very polite tool that will append the line "Thank you for using SPAdes!" to the log file when it has finished running successfully.
Let's use this to our advantage:

```
# Initialize counter to keep track of the number of files containing the string. Optional but a fun feature to use.
count=0

# Loop through all matching log files
for file in path/to/logs/jobname_SPADES.*; do
    # Check if the file exists to avoid errors when no matches are found
    [[ -e "$file" ]] || continue

    # Check if the success message is present
    if grep -q 'Thank you for using SPAdes!' "$file"; then
        ((count++))  # add +1 to count if successful count
    else
        # Print a message to rerun SPAdes for the sample
        echo "Error found in SPAdes run for sample $(basename "$file")"
    fi
done

# Check that the count matches the number of samples
echo "Total successful SPAdes runs: $count"
```

If all your assemblies ran successfully, we can create symbolic links and store them in a new directory with sample names.

```
# Loop through all the SPAdes output directories - they should be named with the respective sample names (*), the assembled contigs are stored in the file contigs.fasta
for fasta in /path/to/where/to/stored/spades_assemblies/*/contigs.fasta; do
    # get sample name (cut -d/ -f9) and replace underscore in sample names to hyphens (sed 's/_/-/g') because Phyluce is a bit finicky when it comes to names
        sample=$(echo ${fasta} | cut -d/ -f9 | sed 's/_/-/g') # adapt -f to number of / in your file path
    # always check that it looks good before running anything else
        #echo ${sample}
    # create symbolic link to assembled contigs into new contig directory
        ln -s ${fasta} /path/to/new/contig/directory/${sample}.contigs.fasta
done
```

This new directory containing assembled contig fasta files for all samples can now be inputted in Phyluce to find UCE loci by matching the contigs to a UCE probe set using the command _phyluce_assembly_match_contigs_to_probes_.

<br>

### Match contigs to UCE probe and locate UCE loci

This step is similar to the one outlined in the Phyluce pipeline (see https://phyluce.readthedocs.io/en/latest/tutorials/tutorial-1.html#finding-uce-loci).
I only want to share an example of how you can structure and organise your scripts and outputs to avoid file-margeddon in your directories.

```
# load Phyluce
module load phyluce

# Storing paths to the new contig directory and UCE probe will be useful to keep at the top of your scripts because you will reuse them down the Phyluce piepline.
CONTIGS="/path/to/new/contig/directory"
UCE_PROBE="/path/to/fish-uce-1k-probes.fasta"

# Create a directory to store all the log files that Phyluce will produce. It keeps things tidy and accessible.
mkdir -p /path/to/phyluce_output/logs
LOGS="/path/to/phyluce_output/logs"

# Go to the directory where you want the Phyluce outputs to be stored
cd /path/to/phyluce_output/

# Find UCE loci
phyluce_assembly_match_contigs_to_probes \
     --contigs ${CONTIGS} \
     --probes ${UCE_PROBE} \
     --output 3_uce_search_results \
     --log-path ${LOGS}
```

Note that my output directory containing UCE loci matches is labelled as the 3rd step in my workflow. This is optional and one way of organising your outputs.
My UCE processing directory looks something like this:

* 1_trimmed_reads #directory containing all trimmed reads and FastQC results
* 2_assemblies #directory containing all assemblies
* 3_UCE #directory containg UCE search results
* 4_mapping #directory containing mapping results (optional)
* 5_phasing #directory containing phased UCE loci (optional)
* 6_correction #directory containing corrected contigs (optional)
* 7_alignments #directory containing final alignments

Speaking of output directory, make sure that when you are running a Phyluce command, the output directory does not exist, especially if rerunning a command.
Example: if you want to have your UCE hits stored in a separate directory within a UCE search directory such as 3_uce_search_results/dragonidae - you should only create the /3_uce_search_results directory and Phyluce will create the /dragonidae directory when you specify --output 3_uce_search_results/dragonidae

<br>

### Extract UCE loci

To run this step, Phyluce requires a configuration file that lists names of the taxa/samples you want to extract UCE loci for. Let's run a small loop to quickly create it using the assembled contigs directories conviniently named after sample's names.

```
# Create a configuration file and add required header
echo "[all]" > 3_uce_search_results/dragonidae_taxonset.conf
# To run it for all Dragonidae samples let's loop through the assembled contig directories to get names
for dir in 2_assemblies/*; do
    # don't forget that we have the directory containing all the assembled contigs too, so let's skip this one
    if [ "$dir" != "all_contigs" ]; then
        sample=$(basename "$dir")
        echo $sample >> 3_uce_search_results/dragonidae_taxonset.conf
    fi
done
```

You can now run the Phyluce command _phyluce_assembly_get_match_counts_ using this configuration file and keep following the pipeline. Let's meet after you have generated exploded fasta files for each of your samples.<br>
If you are going down the Mapping, Phasing and Correction route - there are some helpful scripts in the **/ScriptCraft repository**.

<br>

### Prepare data to align UCE loci with MAFFT

Phyluce implements mafft and muscle aligners but it takes a while to run if you have large datsets so let's run it with MAFFT directly instead.
First, we need to explode the monolithic fasta files by loci rather than taxon, which is the default output format of the command _phyluce_assembly_explode_get_fastas_file_

```
phyluce_assembly_explode_get_fastas_file \
    --input /3_uce_search_results/all-taxa-incomplete.fasta \
    --output /3_uce_search_results/exploded_fastas_by_loci
```

If at this point you want to subset your dataset to keep only specific taxa in your final alignments, you can parse your fasta directories exploded by taxon how you'd like them to be. Best is to create new directories with symbolink links to target exploded files. Then run the script ./by_loci_from_exploded_taxa.sh from the **/ScriptCraft repository**.

Now that the UCE loci are exploded into separate fasta files, let's align the sequences within each files with MAFFT.

<br>

### Align UCE with MAFFT outside of Phyluce

You can reuse a similar PBS directive header as described in step 1 (assemble contigs). Add the following code to your PBS script, it looks more complicated than it is, there are a couple of extra control/checkpoints that are useful when you are working with lots of UCE loci.

```
# load MAFFT aligner
module load mafft

## Variables to adapt to your run ##
# Resources - must match the PBS directive
THREADS="10"
# if you have multiple set of taxa you can run the following script and just change the BATCH variable
BATCH="dragonidae"
# Create output directories and where to store MAFFT log file
OUTDIR="/path/to/7_alignments"
mkdir -p ${OUTDIR}/${BATCH}/{out_mafft,logs}
LOG=${OUTDIR}/${BATCH}/logs/mafft.log
# Set path to directory containing fasta file exploded by loci
IN="/path/to/3_uce_search_results/${BATCH}/exploded_fastas_by_loci/


echo "----- RUNNING MAFFT ALIGNMENT -----" >> ${LOG}
echo "----- RUNNING MAFFT ALIGNMENT -----"

for fasta in ${IN}/*.unaligned.fasta; do
  uce=$(basename "$fasta" .unaligned.fasta)

  # Count sequences in fasta files - mafft will stop running if there is only 1 sequence
  nseqs=$(grep -c "^>" ${fasta})

  if [ "$nseqs" -le 1 ]; then
    echo "${uce}: was not be aligned because contained only 1 sequence." >> ${LOG}
    continue #this makes sure that fasta files with only sequence will be skipped
  fi

  # Run MAFFT
  output_file=${OUTDIR}/${BATCH}/out_mafft/${uce}.fasta
  mafft --auto --thread ${THREAD} ${fasta} > ${output_file}

  # Check success
  if [ $? -eq 0 ]; then
    echo "${uce}: successfully aligned." >> "$LOG"
  else
    rm -f "$output_file"
    echo "Error: MAFFT failed to align $uce. File was removed from output directory." >> "$LOG"
  fi
done

echo "----- MAFFT ALIGNMENT COMPLETED -----" >> ${LOG}
echo "----- MAFFT ALIGNMENT COMPLETED -----"
```

Output directory should contain uce-number.fasta files with aligned sequences that can be inputted back into the Phyluce pipeline.

<br>

### Generate final alignment matrices
The commands below can be run back to back as a PBS script (make sure to have the set -e Bash setting in) pretty quickly and generally follow the Phyluce pipeline with some modifications, mainly for trimming-algorithm and checking for alignments issues.

```
# Load Phyluce
module load phyluce

# Give path to MAFFT alignment directory
ALN="${OUTDIR}/${BATCH}/out_mafft"
# Create and set path to log file
mkdir -p ${OUTDIR}/${BATCH}/logs
LOGS="{OUTDIR}/${BATCH}/logs"
# Get number of taxa for the final alignments from fasta files exploded by taxon
NTAX=$(ls /path/to/3_uce_search_results/exploded_fastas_by_taxon | wc -l)

# Go to output directory
cd ${OUTDIR}

# Clean  mafft alignments (remove locus names)
echo "----- Running: Removing locus name from alignments -----"
phyluce_align_remove_locus_name_from_files \
    --alignments ${ALN} \
    --output ${OUTDIR}/${BATCH}/1_mafft_clean \
    --input-format fasta \
    --output-format nexus \
    --${LOGS} \
    --cores 10

# Screen alignments for problems before any trimming
echo "----- Running: Screening Alignments for problems -----"
phyluce_align_screen_alignments_for_problems \
    --alignments ${OUTDIR}/${BATCH}/1_mafft_clean \
    --output ${OUTDIR}/${BATCH}/2_mafft_screened \
    --input-format nexus \
    --log-path ${LOGS} \
    --cores 10

# Remove taxa with no data from the screened alignments
echo "----- Running: Removing Empty Taxa -----"
phyluce_align_remove_empty_taxa \
    --alignments ${OUTDIR}/${BATCH}/2_mafft_screened \
    --output ${OUTDIR}/${BATCH}/3_mafft_filtered \
    --input-format nexus \
    --output-format nexus \
    --log-path ${LOGS} \
    --cores 10

# Edge trim filtered - for taxa that diverged more recently (< 30-50 mya)
echo "----- Running: Edge trimming on filtered alignments -----"
phyluce_align_get_trimmed_alignments_from_untrimmed \
    --alignments ${OUTDIR}/${BATCH}/3_mafft_filtered \
    --output ${OUTDIR}/${BATCH}/4_mafft_edge_trimmed \
    --input-format nexus \
    --log-path ${LOGS} \
    --cores 10

# Internal trimming with Gblocks is not supported anymore. Use TrimAl instead, which is an in between edge and internal trimming.
#echo "----- Running: Trimming filtered alignments with trimal -----"
#phyluce_align_get_trimal_trimmed_alignments_from_untrimmed \
#  --alignments ${OUTDIR}/${BATCH}/3_mafft_filtered \
#  --output ${OUTDIR}/${BATCH}/4_mafft_trimal_trimmed \
#  --input-format nexus \
#  --log-path ${LOGS} \
#  --cores 10

# Generate summary data csv file to get stats for trimmed alignments
phyluce_align_get_align_summary_data \
    --alignments ${OUTDIR}/${BATCH}/4_mafft_edge_trimmed \
    --cores 10 \
    --output-stats ${OUTDIR}/${BATCH}/mafft_edge_summary.csv

# Generate 75% complete alignment matrix from trimmed alignments
phyluce_align_get_only_loci_with_min_taxa \
    --alignments ${OUTDIR}/${BATCH}/4_mafft_edge_trimmed \
    --taxa $NTAX \
    --percent 0.75 \
    --output ${OUTDIR}/${BATCH}/${BATCH}_e75 \
    --log-path ${LOGS} \
    --cores 10
```

You've made it!! Final alignments can now be used in your downstream bioinformatic endeavours like infering phylogenomic trees.<br>
**Best fishes and happy coding!**
<br><br>

## Help

If you have any issues running the command lines or scripts in this repo feel free to reach out! - lauriane.baraf@my.jcu.edu.au

<br>

## Acknowledgements & Citations

Scripts were written by Lauriane Baraf, please cite this repo if using them.

If using the Phyluce processing pipeline, please cite:<br>
BC Faircloth, McCormack JE, Crawford NG, Harvey MG, Brumfield RT, Glenn TC. 2012. Ultraconserved elements anchor thousands of genetic markers spanning multiple evolutionary timescales. Systematic Biology 61: 717–726. doi:10.1093/sysbio/SYS004.

If using SPAdes please cite:<br>
Prjibelski, A., Antipov, D., Meleshko, D., Lapidus, A., & Korobeynikov, A. (2020). Using SPAdes de novo assembler. Current protocols in bioinformatics, 70(1), e102.

If using MAFFT v7, please cite:<br>
Katoh, K., & Standley, D. M. (2013). MAFFT multiple sequence alignment software version 7: improvements in performance and usability. Molecular biology and evolution, 30(4), 772-780.

<br>

## Publication & Data

This repository contains some of the code initially written to process UCE data harvested from targeted capture sequencing for the reef-associated fish family Pomacanthidae (marine angelfishes). Results are presented in the following publication:<br>
Baraf, L. M., Hung, J. Y., & Cowman, P. F. (2025). Phylogenomics of marine angelfishes: diagnosing sources of systematic discordance for an iconic reef fish family (F: Pomacanthidae). Systematic Biology, syaf016.

Associated raw sequencing data is available on the NCBI GenBank Database under the BioProject PRJNA1101094
Alignments and tree files can be accessed on the Dryad repository 10.5061/dryad.r4xgxd2n4


