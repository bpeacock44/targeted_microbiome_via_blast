#!/bin/bash 
HDIR=/sw/paul_helper_scripts

# Set strict mode
set -euo pipefail

# Define initial values
reblast_iteration="rb0"
maxseqs=5000

DIR=$1
BLAST_FILE=$2
RUN_TYPE=$3

timestamp="$(date +"%Y%m%d_%H:%M:%S")"
output_file="${DIR}/blast.${timestamp}.log"
exec > "$output_file" 2>&1

first_run=true  # Add a flag for the first run

criteria_met() {
    if [ "$first_run" == true ]; then
    	echo "Skipping criteria check for the first run."
    	return 0
	fi

    echo "Checking criteria..."
    reblast_file="${DIR}/zotus/rep_set/${reblast_iteration}.fasta"

    # Check if the reblast file has any lines
    if [ ! -s "$reblast_file" ] || [ "$(wc -l < "$reblast_file")" -eq 0 ]; then
        echo "${reblast_file} is empty. Ending the loop."
        return 1 # End if the reblast file is empty
    else
        echo "${reblast_file} has content. Continuing..."
    fi

    return 0 # Continue if reblast file has content
}

while criteria_met; do
    job_ids=()
    echo "Starting first batch of jobs..."

    if "$first_run"; then
        input_file="${DIR}/zotus/rep_set/seqs_chimera_filtered_otus.fasta"
    else
        input_file="${DIR}/zotus/rep_set/${reblast_iteration}.fastq"
    fi

    if [[ "$RUN_TYPE" == local ]]; then
        # Run locally
        ${BLAST_FILE} "${input_file}" "${maxseqs}" > "${DIR}/zotus/rep_set/${maxseqs}.${reblast_iteration}.blastout"
        echo "BLAST completed locally for $input_file"
    else
        # Submit as sbatch job
        echo sbatch -o "${DIR}/zotus/rep_set/${maxseqs}.${reblast_iteration}.blastout" ${BLAST_FILE} "${input_file}" "${maxseqs}" 
        job_id=$(sbatch -o "${DIR}/zotus/rep_set/${maxseqs}.${reblast_iteration}.blastout" ${BLAST_FILE} "${input_file}" "${maxseqs}" | awk '{print $NF}')
        echo "Blast job submitted with ID: $job_id for $input_file"
        job_ids+=("$job_id")

        # Wait until the job starts running
        while true; do
            # Check the job status
            job_status=$(squeue -j "$job_id" -h -o "%t" 2>/dev/null)

            # Break the loop if the job is in the "running" state
            if [ "$job_status" == "R" ]; then
                echo "Job $job_id is now running."
                break
            fi

            # Sleep for a short interval before checking again
            sleep 10
        done

        # Wait for the submitted job to finish
        while squeue -j "$job_id" -h &>/dev/null; do
            sleep 10
        done

        echo "Job $job_id has completed."
    fi

    next_reblast_value() {
        local current="$1"
        case $current in
            "rb0") echo "rb1" ;;
            "rb1") echo "completed" ;;
            *) echo "error"; exit 1 ;;
        esac
    }

    next_value=$(next_reblast_value "$reblast_iteration")
    echo $next_value
    if [[ "$next_value" == "completed" ]]; then
        echo "All reblast iterations have been completed."
        exit 0
    elif [[ "$next_value" == "error" ]]; then
        echo "Unknown reblast iteration value: $reblast_iteration"
        exit 1
    fi

    bout="${DIR}/zotus/rep_set/${maxseqs}.${reblast_iteration}.blastout"
 	outfile="${DIR}/zotus/rep_set/${maxseqs}.${reblast_iteration}.blastout.not_enough_hits.txt"
	${HDIR}/reblast_check.pl ${bout} ${outfile}
    echo $bout 
    echo $outfile
	# Get the number of reads contributing to the biggest OTU
	total=$(awk 'NR==1{print $NF}' ${DIR}/zotus/rep_set/seqs_chimera_filtered_otus.fasta)
    echo "Determining which OTUs are worth blasting."
    echo ${total}

    not_enough_hits_file="${DIR}/zotus/rep_set/${maxseqs}.${reblast_iteration}.blastout.not_enough_hits.txt"
    big_OTUs_file="${DIR}/zotus/rep_set/${maxseqs}.${reblast_iteration}.blastout.not_enough_hits.big_OTUs.txt"
    reblast_seqs_file="${DIR}/zotus/rep_set/${next_value}.fasta"
    
    # Check if the not_enough_hits_file is empty
    if [ -s "$not_enough_hits_file" ]; then
        # Proceed with grep and subsequent commands
        grep -w -f "$not_enough_hits_file" ${DIR}/zotus/rep_set/seqs_chimera_filtered_otus.fasta | \
        awk -v total="$total" '{
            # Remove ">" from the first column
            gsub(">", "", $1)
    
            # Calculate the ratio and store the result after the white space
            ratio = $2 / total
    
            # Print lines where the ratio is 0.01 or less
            if (ratio <= 0.01) {
                print $1, $2, ratio
            }
        }' | \
        awk '$3 > 0.01' > "$big_OTUs_file"
    
        # Additional commands
        awk '/^>/{sub(/ .*/, "");}1' "${DIR}/zotus/rep_set/seqs_chimera_filtered_otus.fasta" > ${DIR}/zotus/rep_set/modified_seqs.fasta
    
        echo "Extracting reblast seqs."
        # Extract the fasta sequences needing a reblast
        seqkit grep -n -f "$big_OTUs_file" ${DIR}/zotus/rep_set/modified_seqs.fasta -o "$reblast_seqs_file"
    else
        # If the not_enough_hits_file is empty, create an empty file for reblast_seqs
        touch "$reblast_seqs_file"
    fi

    case $reblast_iteration in
        "rb0") maxseqs=30000; reblast_iteration="rb1" ;;
        "rb1") echo "Completed all reblast iterations. Any OTUs that may require further reBLASTing are stored in ${DIR}/zotus/rep_set/${next_value}.fasta"; exit 0 ;;
    esac

    first_run=false
    echo "First loop completed."
done

echo "Merging all blastout files."

rm -f ${DIR}/zotus/rep_set/final.blastout
cat ${DIR}/zotus/rep_set/5000.rb0.blastout | grep -v "# BLAST processed" >> ${DIR}/zotus/rep_set/final.blastout 
if [ -e "${DIR}/zotus/rep_set/30000.rb1.blastout" ]; then
    cat "${DIR}/zotus/rep_set/30000.rb1.blastout" | grep -v "# BLAST processed" >> "${DIR}/zotus/rep_set/final.blastout"
fi

echo "# BLAST processed" >> ${DIR}/zotus/rep_set/final.blastout 


