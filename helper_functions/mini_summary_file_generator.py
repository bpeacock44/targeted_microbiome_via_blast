#!/usr/bin/env python3

import pandas as pd
import numpy as np
import argparse

def process_data_and_write_summary(norm, tax, tax_c, output_file):
    # Read the main data file without recognizing any comment characters
    df = pd.read_csv(norm, sep='\t', comment=None, header=1)
    df.rename(columns={df.columns[0]: 'ID'}, inplace=True)
        
    # Read primary taxonomy file
    tax_df = pd.read_csv(tax, sep='\t', comment=None, header=0)
    tax_df.columns = ['ID', 'taxonomy', 'bitscore', 'per_ID', 'per_qcov']
    df = pd.merge(df, tax_df[['ID', 'taxonomy', 'per_ID', 'per_qcov']], on='ID', how='left')
    
    # Check if tax_c is provided and process it if available
    if tax_c and tax_c.lower() != 'null':
        # Read classifier taxonomy file if given
        tax_df_C = pd.read_csv(tax_c, sep='\t', comment=None, header=0)
        tax_df_C.columns = ['ID', 'c_taxonomy', 'confidence']
        df = pd.merge(df, tax_df_C, on='ID', how='left')
        # Include c_taxonomy and confidence in excluded columns
        excluded_columns = ['ID', 'taxonomy', 'per_ID', 'per_qcov', 'c_taxonomy', 'confidence', 'sequence']
    else:
        # Exclude c_taxonomy and confidence if tax_c is not provided
        excluded_columns = ['ID', 'taxonomy', 'per_ID', 'per_qcov', 'sequence']
    
    # Calculate mean abundance using only numeric columns
    numeric_cols = df.select_dtypes(include=[np.number]).columns.difference(excluded_columns)
    df['avg_abun'] = df[numeric_cols].mean(axis=1, skipna=True)

    # Remove numeric columns after 'avg_abun' has been created
    df.drop(columns=numeric_cols, inplace=True)

    # Write dataframe to the specified output file
    df.to_csv(output_file, sep='\t', index=False, quoting=False)

def main():
    # Create an argument parser
    parser = argparse.ArgumentParser(description='Process data and write summary to file.')

    # Define command-line arguments
    parser.add_argument('norm', type=str, help='Path to the normalized data file')
    parser.add_argument('tax', type=str, help='Path to the primary taxonomy file')
    parser.add_argument('tax_c', type=str, nargs='?', default=None, help='Path to the classifier taxonomy file (optional, use "NULL" to skip)')
    parser.add_argument('output_file', type=str, help='Path to the output file')

    # Parse command-line arguments
    args = parser.parse_args()

    # Call the function with the arguments
    process_data_and_write_summary(args.norm, args.tax, args.tax_c, args.output_file)

if __name__ == '__main__':
    main()
