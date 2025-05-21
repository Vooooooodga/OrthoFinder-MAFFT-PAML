import csv
import os
import shutil
import argparse

def main():
    parser = argparse.ArgumentParser(description="Search for Orthogroups for candidate genes and copy MSA files.")
    parser.add_argument("--gene_list_file", default='DEG.csv', help="Path to the gene list CSV file.")
    parser.add_argument("--orthogroups_file", default='/home/yuhangjia/data/AlternativeSplicing/exon_expansion_test/orthofinder/output_msa_2/Results_Jul05/Orthogroups/Orthogroups.tsv', help="Path to the Orthogroups TSV file.")
    parser.add_argument("--msa_dir", default='/home/yuhangjia/data/AlternativeSplicing/exon_expansion_test/orthofinder/output_msa_2/Results_Jul05/MultipleSequenceAlignments', help="Directory containing MSA files.")
    parser.add_argument("--output_dir", default='selected_msa_files_deg', help="Directory to save selected MSA files.")
    parser.add_argument("--species", default='Apis_mellifera', help="Name of the species column in Orthogroups.tsv to search for genes.")
    args = parser.parse_args()

    # Create output directory if it doesn't exist
    if not os.path.exists(args.output_dir):
        os.makedirs(args.output_dir)

    # Read gene list
    with open(args.gene_list_file, newline='') as csvfile:
        gene_reader = csv.reader(csvfile)
        gene_list = {row[0] for row in gene_reader}

    # Read Orthogroups file and map gene names to Orthogroups
    gene_to_orthogroup = {}
    species_column_index = -1

    with open(args.orthogroups_file, newline='') as csvfile:
        tsv_reader = csv.reader(csvfile, delimiter='\t')
        header = next(tsv_reader)  # Read header row
        try:
            species_column_index = header.index(args.species)
        except ValueError:
            print(f"Error: Species '{args.species}' not found in the header of {args.orthogroups_file}.")
            print(f"Available columns are: {header}")
            exit(1)

        for row in tsv_reader:
            orthogroup = row[0]
            if species_column_index < len(row):
                species_specific_genes = row[species_column_index]
                if species_specific_genes:
                    genes = species_specific_genes.split(', ')
                    for gene in genes:
                        if gene in gene_list:
                            gene_to_orthogroup[gene] = orthogroup
            # else: # Handle cases where a row might be shorter than expected, though TSV should be consistent
            #     print(f"Warning: Row for orthogroup {orthogroup} is shorter than expected and does not contain species column {args.species}")

    # Confirm found Orthogroups
    if not gene_to_orthogroup:
        print(f"No genes from the list found for species '{args.species}' in {args.orthogroups_file}.")
        # Optionally, exit if no genes found, or proceed to see if any orthogroups were found (e.g. if multiple species were relevant)
    found_orthogroups = set(gene_to_orthogroup.values())
    if not found_orthogroups:
        print("No orthogroups found for the selected genes and species.")
        # exit(0) # Depending on desired behavior

    # Iterate through MSA file directory and find matching files
    copied_files_count = 0
    for msa_file in os.listdir(args.msa_dir):
        if msa_file.split('.')[0] in found_orthogroups:
            source_file = os.path.join(args.msa_dir, msa_file)
            destination_file = os.path.join(args.output_dir, msa_file)
            shutil.copy(source_file, destination_file)
            print(f'Copied {msa_file} to {args.output_dir}')
            copied_files_count += 1

    if copied_files_count > 0:
        print(f"MSA files selection completed. Copied {copied_files_count} files.")
    else:
        print("MSA files selection completed. No matching MSA files were found to copy.")

if __name__ == '__main__':
    main()