#!/bin/bash

### USAGE ###
# THIS IS A TEMPORARY PART 3 SCRIPT FOR RUNNING A DIFFERENT ALGORITHM (UPARSE) FOR ASV PICKING.
# CODE FOLLOWS HERE #

set -e

# Custom error handler function
error_handler() {
    local error_message=$1
    echo "Error on line $error_message" | tee /dev/tty
}

# Trap errors and call the error handler
trap 'error_handler "$BASH_COMMAND"' ERR

# ARGUMENTS
while getopts ":d:j:l:o:m:n:u:" opt; do
  case $opt in
    d) DIR="$OPTARG"
    ;;
    j) IFS=',' read -ra JBS <<< "$OPTARG"
    ;;
    l) LEN="$OPTARG"
    ;;
    o) OUTDIR="$OPTARG"
    ;;    
    m) mmatchnum="$OPTARG"
    ;;
    n) MIN="$OPTARG"
    ;;   
    u) UPER="$OPTARG"
    ;;   
    \?) echo "Invalid option -$OPTARG" >&2
    ;;
  esac
done

# Check for mandatory arguments
if [ -z "$DIR" ] || [ -z "${JBS[*]}" ] || [ -z "$LEN" ] || [ -z "$OUTDIR" ]; then
    echo "Usage: $0 -d <directory_path> -j <comma-separated output directories> -l <trim length> -o <desired name of output dir> -u <uparse % similarity for clustering ASVs> [-m <mismatch number> -n <min size for asv calling>]"
    exit 1
fi

# Check if mmatchnum is not defined or is set to an empty string, set it to "0"
if [ -z "${mmatchnum+x}" ] || [ -z "$mmatchnum" ]; then
    mmatchnum="0"
fi

# Remove "_output" from each element in the array
for ((i=0; i<${#JBS[@]}; i++)); do
  JBS[$i]=${JBS[$i]%%_output}
done

# Check uparse value
if [[ ! $UPER =~ ^[0-9]+(\.[0-9]+)?$ ]]; then
  echo "Error: -u must be a valid floating-point number."
  exit 1
fi

echo " - -- --- ---- ---- --- -- -"
echo "Checking for input files"
echo " - -- --- ---- ---- --- -- -"

# show your fastq files and map
for JB in "${JBS[@]}"; do
    if [ ! -e "${DIR}/${JB}_output/${JB}_A1P1.M${mmatchnum}.fq" ]; then
        echo "File ${DIR}/${JB}_output/${JB}_A1P1.M${mmatchnum}.fq not found!"
        exit 1
    fi

    if [ ! -e "${DIR}/${JB}_output/${JB}_A1P2.M${mmatchnum}.fq" ]; then
        echo "File ${DIR}/${JB}_output/${JB}_A1P2.M${mmatchnum}.fq not found!"
        exit 1
    fi
done

# initiate log
timestamp="$(date +"%Y%m%d_%H:%M:%S")"
output_file="${DIR}/part3_${timestamp}.log"
exec > "$output_file" 2>&1

# log header
echo " - -- --- ---- ---- --- -- -"
echo "Log file for Part 3 of the Microbiome Pipeline. Processing the following arguments:
Working Directory: ${DIR}
Data IDs: ${JBS[@]}
Trim Length: ${LEN}
Output Directory: ${OUTDIR}
Mismatches if specified: ${mmatchnum}
 - -- --- ---- ---- --- -- -"

echo
echo " - -- --- ---- ---- --- -- -"
echo "Trimming and Filtering Reads"
echo " - -- --- ---- ---- --- -- -"

# set up and go to output directory
output_dir="${DIR}/${OUTDIR}"
mkdir -vp "${DIR}/${OUTDIR}"
echo "This folder contains the results of ${JBS[@]} trimmed at ${LEN}." > "${output_dir}/summary.txt"

# make tax directory
mkdir -vp "${DIR}/${OUTDIR}/tax_dir"
TAXDIR="${DIR}/${OUTDIR}/tax_dir"

#truncate reads at LEN
for JB in ${JBS[@]}; do
  usearch -fastx_truncate "${DIR}/${JB}_output/${JB}_A1P1.M${mmatchnum}.fq" -quiet -trunclen ${LEN} -fastqout "${DIR}/${JB}_output/${JB}_A1P1_${LEN}bp.fq"; # &
done

rm -f "${output_dir}/combined.fq"
if [ ${#JBS[@]} -gt 1 ]; then
  #EITHER combine if you have multiple files
  for JB in ${JBS[@]}; do
  cat "${DIR}/${JB}_output/${JB}_A1P1_${LEN}bp.fq" >> "${output_dir}/combined.fq"
  done
else
  #OR create a symlink if you have only one file
  ln -s "${DIR}/${JB}_output/${JB}_A1P1_${LEN}bp.fq" "${output_dir}/combined.fq"
fi

#maxee quality filtering of demultiplexed/truncated fq files (*** keep THREADS=1 for repeatability ***)
for JB in ${JBS[@]}; do
  usearch -threads 1 -fastq_filter "${DIR}/${JB}_output/${JB}_A1P1_${LEN}bp.fq" -quiet -fastq_maxee 1.0 -fastaout "${DIR}/${JB}_output/${JB}.filtered.fa"
done
echo
echo " - -- --- ---- ---- --- -- -"
echo "Pooling Samples and Creating ASVs"
echo " - -- --- ---- ---- --- -- -"

#sample pooling (https://www.drive5.com/usearch/manual/pool_samples.html)
rm -f "${output_dir}/filtered.fa"
if [ ${#JBS[@]} -gt 1 ]; then
  #EITHER combine if you have multiple files
  for JB in ${JBS[@]}; do
  cat "${DIR}/${JB}_output/${JB}.filtered.fa" >> "${output_dir}/filtered.fa"
  done
else
  #OR create a symlink if you have only one file
  ln -s "${DIR}/${JB}_output/${JB}.filtered.fa" "${output_dir}/filtered.fa"
fi

#find unique sequences
usearch -fastx_uniques "${output_dir}/filtered.fa" -quiet -fastaout "${output_dir}/uniques.fa" -sizeout -relabel Uniq

#make a subdirectory for the asvs
mkdir -vp "${output_dir}/asvs"

echo usearch -cluster_smallmem "${output_dir}/uniques.fa" -id ${UPER} -centroids "${output_dir}/asvs/asvs.fa"

# MODIFIED LINE - CLUSTERING at 100% IDENTITY
usearch -cluster_smallmem "${output_dir}/uniques.fa" -id ${UPER} -centroids "${output_dir}/asvs/asvs.fa"

# Convert headers
sed 's/>Uniq\([0-9]*\);.*$/>Asv\1/' "${output_dir}/asvs/asvs.fa" > "${output_dir}/asvs/z.fa"

# Check if the replacement was successful before overwriting
if grep -q '>Asv' "${output_dir}/asvs/z.fa"; then
    # Overwrite the original file if '>Asv' is found
    mv -v "${output_dir}/asvs/z.fa" "${output_dir}/asvs/asvs.fa"
else
    echo "Header replacement failed. Original file not overwritten."
fi
echo
echo " - -- --- ---- ---- --- -- -"
echo "Creating Initial ASV Table"
echo " - -- --- ---- ---- --- -- -"

#create an ASV table ("Input should be reads before quality filtering and before discarding low-abundance unique sequences, e.g. singletons")
usearch --otutab "${output_dir}/combined.fq" -quiet -zotus "${output_dir}/asvs/asvs.fa" -otutabout "${output_dir}/asvs/asv_table_00.txt"
sed -i 's/#OTU/#ASV/g' "${output_dir}/asvs/asv_table_00.txt"

# Run python script to sort qiime table
qiime_table_sorter.py "${output_dir}/asvs/asv_table_00.txt" "${output_dir}/asvs/asv_table_01.txt"

#source for bash helper functions
source "qiime_shell_helper_functions.sh"

#convert to biom
OTBL=asv_table_01
txt2biom_notax "${output_dir}/asvs/${OTBL}.txt" "${output_dir}/asvs/${OTBL}.biom"

#add counts to ASV file
otblfp="${output_dir}/asvs/asv_table_01.txt"
fastafp="${output_dir}/asvs/asvs.fa"
outfp="${output_dir}/asvs/seqs_chimera_filtered_ASVs.fasta"
Rscript -e "source('/helper_functions/pipeline_helper_functions.R'); add_counts_to_fasta_sequences('$otblfp', '$fastafp', '$outfp')"

mkdir -vp "${output_dir}/asvs/rep_set"
mv -v "${output_dir}/asvs/seqs_chimera_filtered_ASVs.fasta" "${output_dir}/asvs/rep_set"
echo
echo " - -- --- ---- ---- --- -- -"
echo "It's time to assign taxonomy! You will either do this by BLAST-ing your sequences against a 
local BLAST database or by using Qiime2 to classify them." | tee /dev/tty








