#!/usr/bin/env python3
# Usage: python mixed_family_checker.py <blastout.file> --email <email>

# This checks the first 10 blast results for each ASV that can return a family in it's taxonomic info
# and checks if the families in that top 10 group are identical or not. If they are not, the ASV is written
# to a new file called "mixed_family_checker_out.txt." This information may be helpful in assessing
# the taxonomic assignment's quality.

from collections import defaultdict
from Bio import Entrez
from urllib.error import HTTPError
import time
import argparse

# Function to parse BLAST output file
def parse_blast_output(file_path):
    queries = defaultdict(list)
    current_query = None
    with open(file_path, 'r') as file:
        for line in file:
            if line.startswith("# Query:"):
                if current_query:
                    yield current_query, queries[current_query]  # Yield the current query and its tax IDs
                    queries[current_query] = []  # Reset the tax IDs list
                current_query = line.strip().split(":")[-1].strip().split()[0]
            elif not line.startswith("#") and current_query:
                fields = line.strip().split("\t")
                tax_id = fields[7].split(";")[0]
                queries[current_query].append(tax_id)
    if current_query:
        yield current_query, queries[current_query]  # Yield the last query and its tax IDs


# Initialize tax_cache as an empty dictionary
tax_cache = {}
failed_requests = set()  # Track failed requests

def fetch_taxonomy(tax_id, email):
    try:
        # If taxonomy for this tax_id is already fetched, return from cache
        if tax_id in tax_cache:
            # If the cached value is not None, return it directly
            if tax_cache[tax_id] is not None:
                return tax_cache[tax_id]
        
        # Check if this tax_id has previously failed
        if tax_id in failed_requests:
            return None  # Return None directly if it previously failed

        Entrez.email = email  
        handle = Entrez.efetch(db="taxonomy", id=tax_id, retmode="xml")
        record = Entrez.read(handle)[0]
        lineage = record.get('LineageEx', [])
        for entry in lineage:
            if entry.get('Rank') == 'family':
                # Store fetched taxonomy in cache
                tax_cache[tax_id] = entry['ScientificName']
                return entry['ScientificName']
    except HTTPError as err:
        if 500 <= err.code <= 599 or err.code == 400:
            print("Received error from server for tax_id {}: {}".format(tax_id, err))
            failed_requests.add(tax_id)  # Add tax_id to failed requests
        else:
            print("Error from server for tax_id {}: {}".format(tax_id, err))
    except Exception as e:
        print("Error fetching taxonomy for tax_id {}: {}".format(tax_id, e))
    
    return None

def process_identifiers(identifiers, output_file, email):
    first_family = None
    family_count = 0
    
    with open(output_file, "a") as file:
        for identifier, tax_id in identifiers.items():
            for tax_id_sub in tax_id:
                family = fetch_taxonomy(tax_id_sub, email)
                if family:
                    if first_family is None:
                        first_family = family
                        family_count = 1
                    elif family != first_family:
                        file.write(identifier + "\n")  # Write directly to the file
                        break  # No need to continue processing this identifier
                    else:
                        family_count += 1
                    
                    if family_count >= 10:
                        break  # Exit the loop if the first family appears 10 times

# Main function
def main():
    parser = argparse.ArgumentParser(description='Process BLAST output file.')
    parser.add_argument('blastout_file', type=str, help='Path to BLAST output file')
    parser.add_argument('--email', type=str, help='Your email address for Entrez')
    args = parser.parse_args()
    output_file = "mixed_family_checker_out.txt"

    blast_output_file = args.blastout_file
    email = args.email
    for query, tax_ids in parse_blast_output(blast_output_file):
        process_identifiers({query: tax_ids}, output_file, email)

if __name__ == "__main__":
    main()
