#!/bin/bash
# See mbio_tutorial.md for further guidance!

### USAGE ###
#-d: a working directory, which contains one folder for each of your fastq files named by ID
#-o: the name of your output directory
#-b: the path to your blast script file (only if you are running blast)
#-e: email of the user for NCBI purposes

#Optional arguments:
#-t This is a filter file - it will have four tab-delimited columns. 
#Include any taxonomic groups that you want to preferentially keep or reject as well as their taxonomic ID. 
#ALL taxonomies included under these taxonomic IDs will be treated accordingly so check NCBI
#and make sure. If you are doing a universal assay, do not include the -t flag and DO include the -u flag.
#-u: universal assay - causes final ASV tables to be split into taxonomic groups prior to normalizing
#-s: skip the blast - skips the blast portion - useful for troubleshooting or re-running taxonomy assignment 
    #steps etc. Note that if -s is enabled, -b is not required.
#-j: this flag creates a specialized excel summary output that Dr. Borneman specifically requested. 
    #Runtime will increase, as it requires an analysis examining the top 10 blast hits for each ASV.

# CODE FOLLOWS HERE #

#!/bin/bash

set -e

# Custom error handler function
error_handler() {
    local error_message=$1
    echo "Error on line $error_message" | tee /dev/tty
}

# Trap errors and call the error handler
trap 'error_handler "$BASH_COMMAND"' ERR

# ARGUMENTS
split_asv_table=false
skip_blast=false
james_sum_file_gen=false

while getopts ":d:o:b:e:t:usj" opt; do
  case $opt in
    d) DIR="$OPTARG"
    ;;
    o) OUTDIR="$OPTARG"
    ;;    
    b) blast_file="$OPTARG"
    ;;
    e) EMAIL="$OPTARG"
    ;;
    t) FILTERFILE="$OPTARG"
    ;;    
    u) split_asv_table=true
    ;;
    s) skip_blast=true
    ;;
    j) james_sum_file_gen=true
    ;;
    \?) echo "Invalid option -$OPTARG" >&2
    ;;
  esac
done


shift $((OPTIND -1))

# Check for mandatory arguments
if [ -z "$DIR" ] || [ -z "$OUTDIR" ] || [ -z "$EMAIL" ]; then
    echo "Usage: $0 -d <directory_path> -o <desired name of output dir> -e <email@email.com> [-b <blast parameter file> -t <filtertax_file> -u -s -j]"
    exit 1
fi

# If blast is not skipped, check for the blast_file arguments
if [ "$skip_blast" = false ]; then
    if [ -z "$blast_file" ]; then
        echo "Usage: $0 -d <directory_path> -o <desired name of output dir> -b <blast parameter file> -e <email@email.com> [-t <filtertax_file> -u -s -j]"
        exit 1
    fi
fi

# Check if email is in correct format
if ! [[ "${EMAIL}" =~ ^[a-zA-Z0-9._%+-]+@[a-zA-Z0-9.-]+\.[a-zA-Z]{2,}$ ]]; then
    echo "Email is not valid."
    exit 1
fi

# Define initial paths and run_type variable.
output_dir="${DIR}/${OUTDIR}"
TAXDIR="${DIR}/${OUTDIR}/tax_dir"
mkdir -vp $TAXDIR
run_type=local

echo " - -- --- ---- ---- --- -- -"
echo "Checking for input files"
echo " - -- --- ---- ---- --- -- -"

if [ ! -e "${output_dir}/asvs/rep_set/seqs_chimera_filtered_ASVs.fasta" ]; then
    echo "${output_dir}/asvs/rep_set/seqs_chimera_filtered_ASVs.fasta not found!"
    exit 1
fi

if [ ! -e "${output_dir}/asvs/asv_table_01.biom" ]; then
    echo "${output_dir}/asvs/asv_table_01.biom not found!"
    exit 1
fi

# initiate log
timestamp="$(date +"%Y%m%d_%H:%M:%S")"
output_file="${DIR}/part4_${timestamp}.log"
exec > "$output_file" 2>&1

# log header
echo " - -- --- ---- ---- --- -- -"
echo "Log file for Part 4 of the Microbiome Pipeline. Processed the following arguments:
Working directory: ${DIR}
Output directory: ${OUTDIR}
Email of user: ${EMAIL}"

if [ -z "$FILTERFILE" ]; then
    echo "No filter file provided."
else
    echo "Filter file used: ${FILTERFILE}"
fi

if [ "$skip_blast" = true ]; then
    echo "BLAST was skipped."
else
    echo "Blast run file: ${blast_file}"
fi

if [ "$split_asv_table" = true ]; then
    echo "Final ASV tables were split into three domains of life (for universal assay data)"
fi

if [ "$james_sum_file_gen" = true ]; then
    echo "A James Summary File was generated."
fi

echo " - -- --- ---- ---- --- -- -"


if [ "$skip_blast" = false ]; then
    echo
    echo " - -- --- ---- ---- --- -- -"
    echo "Running BLAST on ASVs"
    echo " - -- --- ---- ---- --- -- -"
    
    # Run BLAST script in the background
    blast_iterator.sh "${output_dir}" "${blast_file}" ${run_type} &
    
    # Get the process ID of the last background command
    blast_pid=$!
    
    # Wait for the process to finish
    while ps -p $blast_pid > /dev/null; do
        sleep 1
    done
    
    # Check the exit status of the BLAST script
    blast_exit_status=$?
    if [ $blast_exit_status -eq 0 ]; then
        echo "BLAST successfully completed."
    else
        echo "BLAST failed with exit status $blast_exit_status."
        exit 1
    fi
else
    echo "Skipping BLAST as per user instructions."
fi

echo
echo " - -- --- ---- ---- --- -- -"
echo "Determining Likely Taxonomy of ASVs Using Filters"
echo " - -- --- ---- ---- --- -- -"

# Create filter files in the taxonomy directory if needed. There are two sections - one for if the user didn't specify
# a filter file and another for if they did.
# Read the tab-delimited file, skipping the header

used_taxa=()

add_file_names() {
    local fand="$1"
    local fnot="$2"

    # Extract just the filenames from the paths
    fand_filename=$(basename "$fand")
    fnot_filename=$(basename "$fnot")

    used_taxa+=("$fand_filename" "$fnot_filename")
}

# Process based on FILTERFILE presence
if [ -z "$FILTERFILE" ]; then
    while IFS=$'\t' read -r Name ID Rank Action; do
        # Generate filenames
        FAND="${TAXDIR}/${Name}__${Rank}_txid${ID}_AND_Environmental_Samples.txt"
        FNOT="${TAXDIR}/${Name}__${Rank}_txid${ID}_NOT_Environmental_Samples.txt"
        # Check the action column and create files accordingly unless they already exist
        if [ "$Action" == "Keep" ]; then
            touch "$FAND"
            touch "$FNOT"
            echo "Creating or locating files ${Name}__${Rank}_txid${ID}_AND_Environmental_Samples.txt and ${Name}__${Rank}_txid${ID}_NOT_Environmental_Samples.txt"
            add_file_names "$FAND" "$FNOT"
        elif [ "$Action" == "Reject" ]; then
            touch "$FNOT"
            echo "Creating or locating file ${Name}__${Rank}_txid${ID}_NOT_Environmental_Samples.txt."
            add_file_names "" "$FNOT"
        fi
    done < <(tail -n +2 <(echo -e "Name\tID\tRank\tAction\nEukaryota\t2759\tk\tKeep\nBacteria\t2\tk\tKeep\nArchaea\t2157\tk\tKeep\nPlaceholder\t0\tk\tReject"))
else
    while IFS=$'\t' read -r Name ID Rank Action; do
        # Generate filenames
        FAND="${TAXDIR}/${Name}__${Rank}_txid${ID}_AND_Environmental_Samples.txt"
        FNOT="${TAXDIR}/${Name}__${Rank}_txid${ID}_NOT_Environmental_Samples.txt"
        # Check the action column and create files accordingly unless they already exist
        if [ "$Action" == "Keep" ]; then
            touch "$FAND"
            touch "$FNOT"
            echo "Creating or locating files ${Name}__${Rank}_txid${ID}_AND_Environmental_Samples.txt and ${Name}__${Rank}_txid${ID}_NOT_Environmental_Samples.txt"
            add_file_names "$FAND" "$FNOT"
        elif [ "$Action" == "Reject" ]; then
            touch "$FNOT"
            echo "Creating or locating file ${Name}__${Rank}_txid${ID}_NOT_Environmental_Samples.txt."
            add_file_names "" "$FNOT"
        fi
    done < <(tail -n +2 "$FILTERFILE")
fi

echo

### Update files 
TAXONS=()
TAXIDS=()

for file in "${used_taxa[@]}"; do
    if [[ $file == *NOT* ]]; then
        taxon=$(echo "$file" | perl -ne '@A=split/_NOT_/;@B=split/__\w_/,$A[0];print"$A[0]\n"')
        taxid=$(echo "$file" | perl -ne '@A=split/_NOT_/;@B=split/__\w_/,$A[0];print"$B[1]\n"')
        TAXONS+=("$taxon")
        TAXIDS+=("$taxid")
    fi
done

# create AND/NOT search-terms and filenames for esearch/efetch
# Function to retrieve taxonomy IDs
retrieve_taxonomy() {
    local query="$1"
    local output_file="$2"

    local max_attempts=3
    local attempt=1
    local success=false

    while [ $attempt -le $max_attempts ]; do
        echo "Attempt $attempt: Retrieving taxonomy for $output_file"
        esearch -db taxonomy -query "$query" | efetch -format uid > "$output_file"
        if [ $? -eq 0 ]; then
            success=true
            break
        else
            echo "Attempt $attempt failed. Retrying in 5 seconds..."
            sleep 5
            ((attempt++))
        fi
    done

    if ! $success; then
        echo "Failed to retrieve taxonomy for $output_file after $max_attempts attempts. This may be because you're requesting too many in quick succession. Run part 4 again with the s flag to skip the BLAST."
        exit 1
    fi
}

N=${#TAXONS[@]}

# Main loop
for ((i = 0; i < N; i++)); do
    FAND="${TAXDIR}/${TAXONS[i]}_AND_Environmental_Samples.txt"
    FNOT="${TAXDIR}/${TAXONS[i]}_NOT_Environmental_Samples.txt"

    # Check if the file is not empty and was updated within the last 24 hours
    if [[ -s "$FAND" && $(find "$FAND" -mtime -1 2>/dev/null) ]] && [[ -s "$FNOT" && $(find "$FNOT" -mtime -1 2>/dev/null) ]]; then
        echo "Skipping ${TAXONS[i]} - already up-to-date."
        continue
    fi

    echo "Retrieving taxonomic IDs for ${TAXONS[i]}."
    retrieve_taxonomy "${TAXIDS[i]}[subtree] AND \"Environmental Samples\"[subtree]" "$FAND"
    retrieve_taxonomy "${TAXIDS[i]}[subtree] NOT \"Environmental Samples\"[subtree]" "$FNOT"
done

echo

# download and unzip updated merged.dmp file if it was not already downloaded today
if [[ -f "${TAXDIR}/merged.dmp" ]]; then
    if [[ $(find "${TAXDIR}/merged.dmp" -mtime -1) ]]; then
        echo "Updated merged.dmp file present."
    else
        echo "File ${TAXDIR}/merged.dmp exists but was not updated within the last day."
        download_file=true
    fi
else
    echo "File ${TAXDIR}/merged.dmp does not exist."
    download_file=true
fi

if [ "$download_file" = true ]; then
    wget -q ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdmp.zip -P "${TAXDIR}"
    if [ $? -ne 0 ]; then
        echo "Error: Failed to download taxdmp.zip."
        exit 1
    fi
    unzip -o "${TAXDIR}/taxdmp.zip" -d "${TAXDIR}"
    if [ $? -ne 0 ]; then
        echo "Error: Failed to unzip taxdmp.zip."
        exit 1
    fi
    rm -rf "${TAXDIR}/taxdmp.zip"
fi


touch "${TAXDIR}/Placeholder__k_txid0_NOT_Environmental_Samples.txt"

echo

# find the user's usearch path, which is where the AccnsWithDubiousTaxAssigns.txt file will be stored.
#userpath=$(which usearch)
#userdir=$(dirname "$userpath")

#Mario: Within the container, usearch will ALWAYS exist within the bind directory. These commands write to the binded directory and anything that has been written to this directory will persist on the host's binded directory, even after the container is closed
tax_files_dir="/bind/mbio_taxa_fz"
mkdir -vp "$tax_files_dir" 
touch "${tax_files_dir}/AccnsWithDubiousTaxAssigns.txt" 
chmod 777 "${tax_files_dir}/AccnsWithDubiousTaxAssigns.txt" 

# This section will create the ASVs2filter.log, which will be used to assign taxonomy. 
# Again, there are two sections - one for if the user didn't specify a filter file and another for if they did.
if [ -z "$FILTERFILE" ]; then
  blast_taxa_categorizer.py \
    -i "${output_dir}/asvs/rep_set/final.blastout" \
    -k $(awk -F'\t' '$4=="Keep"{print "'${TAXDIR}'/"$1"__"$3"_txid"$2"_NOT_Environmental_Samples.txt"}' <(echo -e "Name\tID\tRank\tAction\nEukaryota\t2759\tk\tKeep\nBacteria\t2\tk\tKeep\nArchaea\t2157\tk\tKeep\nPlaceholder\t0\tk\tReject") | paste -sd, -) \
    -e $(awk -F'\t' '$4=="Keep"{print "'${TAXDIR}'/"$1"__"$3"_txid"$2"_AND_Environmental_Samples.txt"}' <(echo -e "Name\tID\tRank\tAction\nEukaryota\t2759\tk\tKeep\nBacteria\t2\tk\tKeep\nArchaea\t2157\tk\tKeep\nPlaceholder\t0\tk\tReject") | paste -sd, -) \
    -r $(awk -F'\t' '$4=="Reject"{print "'${TAXDIR}'/"$1"__"$3"_txid"$2"_NOT_Environmental_Samples.txt"}' <(echo -e "Name\tID\tRank\tAction\nEukaryota\t2759\tk\tKeep\nBacteria\t2\tk\tKeep\nArchaea\t2157\tk\tKeep\nPlaceholder\t0\tk\tReject") | paste -sd, -) \
    -t "$tax_files_dir" \
    -m "${TAXDIR}/merged.dmp"
else
  blast_taxa_categorizer.py \
      -i "${output_dir}/asvs/rep_set/final.blastout" \
      -k $(awk -F'\t' '$4=="Keep"{print "'${TAXDIR}'/"$1"__"$3"_txid"$2"_NOT_Environmental_Samples.txt"}' "$FILTERFILE" |paste -sd, -) \
      -e $(awk -F'\t' '$4=="Keep"{print "'${TAXDIR}'/"$1"__"$3"_txid"$2"_AND_Environmental_Samples.txt"}' "$FILTERFILE" |paste -sd, -) \
      -r $(awk -F'\t' '$4=="Reject"{print "'${TAXDIR}'/"$1"__"$3"_txid"$2"_NOT_Environmental_Samples.txt"}' "$FILTERFILE" |paste -sd, -)  \
      -t "$tax_files_dir" \
      -m "${TAXDIR}/merged.dmp" #-f
fi

#Mario: blast_taxa_categorizer.py generates bad_accns.txt in the current directory from which the script was executed.
#Before when the script changed directory into ${output_dir} the file(bad_accns) would be generated in the correct place; now it is generated in the current directory.
#This is because the python script defaults to generating output files in the current directory. Having the file not generated in the correct location would cause 
#issuses later on as the file which the variable bad_accans references, doesn't exist in that directory.
#To solve this, I simply mv the generated file into the output_dir; similar to how the ASV2 files are moved.

#I did notice that blast_taxa_categorizer.py has an output flag (-o) that outputs files into the passed in directory path, but it's commented out
#and the code to append the directory path to the files isn't present. I assume this was done for a reason so, I didn't change anything within blast_taxa_categorizer.py

#mv bad_accns.txt generated into output_dir
mv bad_accns.txt ${output_dir}

bad_accns="${output_dir}/bad_accns.txt"
dubious_accns="${tax_files_dir}/AccnsWithDubiousTaxAssigns.txt"

# Function to check if a number exists in the dubious_accns file
check_existence() {
    local number="$1"
    if grep -q "^$number$" "$dubious_accns"; then
        return 0 # Number exists in the file
    else
        return 1 # Number does not exist in the file
    fi
}

# Filter out lines starting with # or empty lines
# Check if $bad_accns is non-empty and contains lines that are not empty and do not start with #
if [[ -s "$bad_accns" && $(grep -cve '^#|^$' "$bad_accns") -gt 0 ]]; then
    if filtered_lines=$(grep -vE '^#|^$' "$bad_accns" 2>/dev/null); then
        echo "New bad accns: $filtered_lines"
    else
        echo "No new bad accns."
    fi
fi

while IFS= read -r line; do
    number=$(echo "$line" | cut -f1)
    
    if [ -z "$number" ]; then
        continue
    fi

    if ! check_existence "$number"; then
        echo "$number" >> "$dubious_accns" # Append the number if it doesn't exist
        echo "Adding $number to AccnsWithDubiousTaxAssigns.txt"
        new_addition=true
    fi
done <<< "$filtered_lines"

# If new additions were made, re-run the loop
if [ "$new_addition" = true ]; then
    # Run the blast_taxa_categorizer.py script
    if [ -z "$FILTERFILE" ]; then
        blast_taxa_categorizer.py \
            -i "${output_dir}/asvs/rep_set/final.blastout" \
            -k $(awk -F'\t' '$4=="Keep"{print "'${TAXDIR}'/"$1"__"$3"_txid"$2"_NOT_Environmental_Samples.txt"}' <(echo -e "Name\tID\tRank\tAction\nEukaryota\t2759\tk\tKeep\nBacteria\t2\tk\tKeep\nArchaea\t2157\tk\tKeep\nPlaceholder\t0\tk\tReject") | paste -sd, -) \
            -e $(awk -F'\t' '$4=="Keep"{print "'${TAXDIR}'/"$1"__"$3"_txid"$2"_AND_Environmental_Samples.txt"}' <(echo -e "Name\tID\tRank\tAction\nEukaryota\t2759\tk\tKeep\nBacteria\t2\tk\tKeep\nArchaea\t2157\tk\tKeep\nPlaceholder\t0\tk\tReject") | paste -sd, -) \
            -r $(awk -F'\t' '$4=="Reject"{print "'${TAXDIR}'/"$1"__"$3"_txid"$2"_NOT_Environmental_Samples.txt"}' <(echo -e "Name\tID\tRank\tAction\nEukaryota\t2759\tk\tKeep\nBacteria\t2\tk\tKeep\nArchaea\t2157\tk\tKeep\nPlaceholder\t0\tk\tReject") | paste -sd, -) \
            -t "$tax_files_dir" \
            -m "${TAXDIR}/merged.dmp"
    else
        blast_taxa_categorizer.py \
            -i "${output_dir}/asvs/rep_set/final.blastout" \
            -k $(awk -F'\t' '$4=="Keep"{print "'${TAXDIR}'/"$1"__"$3"_txid"$2"_NOT_Environmental_Samples.txt"}' "$FILTERFILE" |paste -sd, -) \
            -e $(awk -F'\t' '$4=="Keep"{print "'${TAXDIR}'/"$1"__"$3"_txid"$2"_AND_Environmental_Samples.txt"}' "$FILTERFILE" |paste -sd, -) \
            -r $(awk -F'\t' '$4=="Reject"{print "'${TAXDIR}'/"$1"__"$3"_txid"$2"_NOT_Environmental_Samples.txt"}' "$FILTERFILE" |paste -sd, -)  \
            -t "$tax_files_dir" \
            -m "${TAXDIR}/merged.dmp" #-f
    fi
fi

new_addition=false

# Move the generated ASVs files to the rep_set folder
if ls ASVs2* 1> /dev/null 2>&1; then
    mv ASVs2* "${output_dir}/asvs/rep_set"
fi

mkdir -vp "${output_dir}/asvs/rep_set/assgntax"

#Mario: Created helper script that activates python virtual environment containing the necessary pip modules
#for the next steps.
source pymods.sh || { echo "Error: Unable to activate python-pip-modules environment"; exit 1; }

echo
echo " - -- --- ---- ---- --- -- -"
echo "Assigning Taxonomy With Filters"
echo " - -- --- ---- ---- --- -- -"
rm -rf "${output_dir}/asvs/rep_set/assgntax/taxonomyDB.json"
rm -f *.xml
blast_assign_taxonomy.py -i "${output_dir}/asvs/rep_set/ASVs2filter.log" \
    --db "${output_dir}/asvs/rep_set/assgntax/taxonomyDB.json" --assign_all --add_sizes \
    -m "${EMAIL}"\
    -o "${output_dir}/asvs/rep_set/assgntax/seqs_chimera_filtered_tax_assignments.txt"

if [ ! -f "${output_dir}/asvs/rep_set/assgntax/seqs_chimera_filtered_tax_assignments.txt" ]; then
    echo "Error: Output file not found."
    exit 1
else
    echo "Blast assign taxonomy completed successfully."
fi

rm *.xml

echo
echo " - -- --- ---- ---- --- -- -"
echo "Printing Taxa Levels With Filters"
echo " - -- --- ---- ---- --- -- -"

# Source Qiime shell helper functions
source qiime_shell_helper_functions.sh || { echo "Error: Unable to source Qiime shell helper functions"; exit 1; }

# Run count_taxa_levels command
count_taxa_levels ${output_dir}/asvs/rep_set/assgntax/seqs_chimera_filtered_tax_assignments.txt > ${output_dir}/asvs/rep_set/assgntax/taxa_levels.txt || { echo "Error: count_taxa_levels command failed"; }

echo
echo " - -- --- ---- ---- --- -- -"
echo "Determining Likely Taxonomy of ASVs Without Filters"
echo " - -- --- ---- ---- --- -- -"
blast_taxa_categorizer.py \
    -i "${output_dir}/asvs/rep_set/final.blastout" \
    -k $(awk -F'\t' '$4=="Reject"{print "'${TAXDIR}'/"$1"__"$3"_txid"$2"_NOT_Environmental_Samples.txt"}' <(echo -e "Name\tID\tRank\tAction\nPlaceholder\t0\tk\tReject") | paste -sd, -) \
    -e $(awk -F'\t' '$4=="Reject"{print "'${TAXDIR}'/"$1"__"$3"_txid"$2"_NOT_Environmental_Samples.txt"}' <(echo -e "Name\tID\tRank\tAction\nPlaceholder\t0\tk\tReject") | paste -sd, -) \
    -r $(awk -F'\t' '$4=="Reject"{print "'${TAXDIR}'/"$1"__"$3"_txid"$2"_NOT_Environmental_Samples.txt"}' <(echo -e "Name\tID\tRank\tAction\nPlaceholder\t0\tk\tReject") | paste -sd, -) \
    -t "$tax_files_dir" \
    -m "${TAXDIR}/merged.dmp" #-f

# Rename the generated ASVs files and move to the rep_set folder
mv ./ASVs2filter.log ${output_dir}/asvs/rep_set/nf_ASVs2filter.log
mv ./ASVs2summary.txt ${output_dir}/asvs/rep_set/nf_ASVs2summary.txt
mv ./ASVs2reject.txt ${output_dir}/asvs/rep_set/nf_ASVs2reject.txt
mv ./ASVs2keep.txt ${output_dir}/asvs/rep_set/nf_ASVs2keep.txt
#Mario: Another bad_accns is generated, which overwrites previous one. Unknown as to what to do with this file so
#I just moved it and renamed it to nf_bad_accns.txt
mv ./bad_accns.txt ${output_dir}/nf_bad_accns.txt

echo
echo " - -- --- ---- ---- --- -- -"
echo "Assigning Taxonomy Without Filters"
echo " - -- --- ---- ---- --- -- -"
rm -rf ${output_dir}/asvs/rep_set/assgntax/nf_taxonomyDB.json
rm -f *.xml
blast_assign_taxonomy.py -i ${output_dir}/asvs/rep_set/nf_ASVs2filter.log \
  --db ${output_dir}/asvs/rep_set/assgntax/nf_taxonomyDB.json --assign_all --add_sizes \
  -m "${EMAIL}" \
  -o ${output_dir}/asvs/rep_set/assgntax/nf_seqs_chimera_filtered_tax_assignments.txt

  if [ ! -f "${output_dir}/asvs/rep_set/assgntax/nf_seqs_chimera_filtered_tax_assignments.txt" ]; then
    echo "Error: Output file not found."
    exit 1
else
    echo "Blast assign taxonomy completed successfully."
fi

rm *.xml

deactivate

# adding extra code to check for and replace non-ascii characters in taxonomies, as this can throw an error later on with the biom functions.
nonascii_finder.py "${output_dir}/asvs/rep_set/assgntax/*seqs_chimera_filtered_tax_assignments.txt"

echo
echo " - -- --- ---- ---- --- -- -"
echo "Printing Taxa Levels Without Filters"
echo " - -- --- ---- ---- --- -- -"
# count taxa levels
# Source Qiime shell helper functions
source qiime_shell_helper_functions.sh || { echo "Error: Unable to source Qiime shell helper functions"; exit 1; }

# Run count_taxa_levels command
count_taxa_levels ${output_dir}/asvs/rep_set/assgntax/seqs_chimera_filtered_tax_assignments.txt > ${output_dir}/asvs/rep_set/assgntax/taxa_levels.txt || { echo "Error: count_taxa_levels command failed"; }

echo
echo " - -- --- ---- ---- --- -- -"
echo "Adding Taxa and Sequences to ASV Tables"
echo " - -- --- ---- ---- --- -- -"

#add taxa to ASV table
OTBL=asv_table_01
biomAddObservations ${output_dir}/asvs/${OTBL}.biom ${output_dir}/asvs/asv_table_02_add_taxa.biom ${output_dir}/asvs/rep_set/assgntax/seqs_chimera_filtered_tax_assignments.txt

# create three additional taxonomic levels of ASV tables
OTBL="asv_table_02_add_taxa"
 
#Mario: Added -o flag which specifies the output directory location
#http://qiime.org/scripts/summarize_taxa.html
summarize_taxa.py -i "${output_dir}/asvs/${OTBL}.biom" -a -L 2,6,7 -o "${output_dir}/asvs"

# will need to chnage if we add more levels
to_process=(
    "${output_dir}/asvs/asv_table_02_add_taxa_L2.biom"
    "${output_dir}/asvs/asv_table_02_add_taxa_L6.biom"
    "${output_dir}/asvs/asv_table_02_add_taxa_L7.biom"
    "${output_dir}/asvs/asv_table_02_add_taxa.biom"
)

for F in "${to_process[@]}"; do
    if [ ! -f "$F" ]; then
        echo "Error: File $F not found."
        exit 1
    fi

    FNAME=$(basename "$F" | sed 's|^asv_table_02_add_taxa||' | sed 's/.biom//')
    ID=$(basename "$F" | sed 's/asv_table_02_add_taxa//' | sed 's/.biom//')

    biom2txt "$F" "${output_dir}/asvs/${OTBL}${ID}.txt"

    if [[ "$split_asv_table" == true ]]; then
        KDOMS=("k__Archaea" "k__Bacteria" "k__Eukaryota")
        for K in "${KDOMS[@]}"; do
            grep -P "(#|$K)" "${output_dir}/asvs/${OTBL}${ID}.txt" > "${output_dir}/asvs/${OTBL}${ID}.${K}.txt"
            NEW_OTBL="${OTBL}${ID}.${K}"
            if ! grep -q "$K" "${output_dir}/asvs/${NEW_OTBL}.txt"; then
                rm "${output_dir}/asvs/${NEW_OTBL}.txt"
            else
                if [[ "${NEW_OTBL}" == *"_L"* ]]; then
                    txt2biom_notax -f "${output_dir}/asvs/${NEW_OTBL}.txt" "${output_dir}/asvs/${NEW_OTBL}.biom"
                else 
                    txt2biom -f "${output_dir}/asvs/${NEW_OTBL}.txt" "${output_dir}/asvs/${NEW_OTBL}.biom"
                fi
                biom_table_math_ops.py -i "${output_dir}/asvs/${NEW_OTBL}.biom" -o "${output_dir}/asvs/${NEW_OTBL}_norm.biom" --normalize2unity
                biom2txt "${output_dir}/asvs/${NEW_OTBL}_norm.biom" "${output_dir}/asvs/${NEW_OTBL}_norm.txt"
            fi
        done
    fi

    biom_table_math_ops.py -i "$F" -o "${output_dir}/asvs/${OTBL}${ID}_norm.biom" --normalize2unity
    biom2txt "${output_dir}/asvs/${OTBL}${ID}_norm.biom" "${output_dir}/asvs/${OTBL}${ID}_norm.txt"
done

# add seqs to L8 
otblfp="${output_dir}/asvs/asv_table_02_add_taxa.txt"
outfp="${output_dir}/asvs/asv_table_03_add_seqs.txt"

Rscript -e "source('${HDIR}/pipeline_helper_functions.R'); add_sequences_to_asv_table('$otblfp', '${output_dir}/asvs/rep_set/seqs_chimera_filtered_ASVs.fasta', '$outfp')"

otblfp="${output_dir}/asvs/asv_table_02_add_taxa_norm.txt"
outfp="${output_dir}/asvs/asv_table_03_add_seqs_norm.txt"

Rscript -e "source('${HDIR}/pipeline_helper_functions.R'); add_sequences_to_asv_table('$otblfp', '${output_dir}/asvs/rep_set/seqs_chimera_filtered_ASVs.fasta', '$outfp')"

to_process2=($(find "${output_dir}/asvs" -maxdepth 1 -type f -name "*taxa.k*txt"))

for F in "${to_process2[@]}"; do
    if [ ! -f "$F" ]; then
        echo "Error: File $F not found."
        exit 1
    fi

    FNAME=$(basename "$F" | sed 's|^./asv_table_02_add_taxa||')
    otblfp="${F}"
    outfp="${output_dir}/asvs/asv_table_03_add_seqs${FNAME}"
    Rscript -e "source('${HDIR}/pipeline_helper_functions.R'); add_sequences_to_asv_table('$otblfp', '${output_dir}/asvs/rep_set/seqs_chimera_filtered_ASVs.fasta', '$outfp')"
done
echo

if [ "$james_sum_file_gen" = true ]; then
    echo " - -- --- ---- ---- --- -- -"
    echo "Creating Summary File"
    echo " - -- --- ---- ---- --- -- -"

    source pymods.sh || { echo "Error: Unable to activate python-pip-modules environment"; exit 1; }

    # Error check for mixed_family_checker.py
    mixed_family_checker.py "${output_dir}/asvs/rep_set/final.blastout" --email "${EMAIL}" || { echo "Error: mixed_family_checker.py failed"; exit 1; }

    # Move output file with error check
    mv mixed_family_checker_out.txt "${output_dir}/asvs/rep_set/" || { echo "Error: Unable to move mixed_family_checker_out.txt"; exit 1; }

    # Error check for Rscript
    Rscript -e "source('${HDIR}/pipeline_helper_functions.R'); process_data_and_write_excel('${output_dir}/asvs/asv_table_03_add_seqs_norm.txt', '${output_dir}/asvs/rep_set/assgntax/nf_seqs_chimera_filtered_tax_assignments.txt', '${output_dir}/asvs/rep_set/assgntax/seqs_chimera_filtered_tax_assignments.txt', '${output_dir}/asvs/asv_table_03_add_seqs.txt', '${output_dir}/asvs/rep_set/mixed_family_checker_out.txt')" || { echo "Error: Rscript failed"; exit 1; }
   
    #Move generated ASVSumamry file
    mv ASV_summary* ${output_dir}/asvs


else
    echo "No Borneman summary file generated, as requested."
fi

echo
echo "All ASV tables have been generated." | tee /dev/tty
echo " - -- --- ---- ---- --- -- -"  | tee /dev/tty
echo "Final Recommendations"  | tee /dev/tty
echo " - -- --- ---- ---- --- -- -"  | tee /dev/tty
echo "Sometimes ASVs can be primarily associated with your controls (likely a source of contamination,
which may be from other samples or even from the individual building the library in the first place). Make sure 
to check your ASV tables for these ASVs - just look at your control columns and see if any of the ASVs are 
relatively high in them and not present in the regular samples. These ASVs should be removed before further 
analysis.

Additionally, there may be ambiguous taxonomic assignments in your ASV table. This may not be problematic if the 
ambigulously assigned ASV isn't abundant, but if it is then you should manually blast the sequence and determine
taxonomy for yourself on NCBI prior to further analyses. If you make changes to taxonomy, note that your L2, L6, 
and L7 ASV tables will likely change as well and you may want to regenerate them. 

Here are some examples of ambiguous taxa:

k__Bacteria;Other
k__Fungi;Other
k__Eukaryota;Other
k__Bacteria;p__unclassified_Bacteria
k__Fungi;p__unclassified_Fungi
k__Eukaryota;p__unclassified_Eukaryota
k__Bacteria_OR_k__unclassified_;Other
k__Fungi_OR_k__unclassified_;Other
k__Eukaryota_OR_k__unclassified_;Other
k__Unassigned;Other
" | tee /dev/tty
