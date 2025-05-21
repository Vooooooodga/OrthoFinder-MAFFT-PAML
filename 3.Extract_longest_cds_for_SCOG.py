from Bio import SeqIO
import os
import argparse

# List of species with subspecies names
subspecies_included = ['Bombus_vancouverensis_nearcticus', 'Osmia_bicornis_bicornis'] # This might be better as a command-line argument if it changes often

# Custom parsing function to extract the most suitable gene identifier from the description
def parse_description(description):
    return description.strip('>')

# Create index, with special handling for parsing gene IDs
def build_index(nucleotide_db_dir):
    index = {}
    for nucleotide_file in os.listdir(nucleotide_db_dir):
        file_path = os.path.join(nucleotide_db_dir, nucleotide_file)
        species_id = os.path.splitext(nucleotide_file)[0]
        record_dict = {}
        for record in SeqIO.parse(file_path, "fasta"):
            key = parse_description(record.description)
            # Check and select the longest transcript
            if key in record_dict:
                if len(record.seq) > len(record_dict[key].seq):
                    record_dict[key] = record
            else:
                record_dict[key] = record
        index[species_id] = record_dict
    return index

def main():
    parser = argparse.ArgumentParser(description="Extract longest CDS for single-copy orthologous groups.")
    parser.add_argument("--protein_msa_dir", default='single_copy_selected_msa_files_deg', help="Directory containing protein MSA files.")
    parser.add_argument("--nucleotide_db_dir", default='../longest_cds', help="Directory containing nucleotide database files.")
    parser.add_argument("--output_dir", default='cds_sequences_deg', help="Directory to save extracted CDS sequences.")
    parser.add_argument("--subspecies_included", nargs='*', default=subspecies_included, help="List of species with subspecies names.") # Use the global or allow override
    args = parser.parse_args()

    if not os.path.exists(args.output_dir):
        os.makedirs(args.output_dir)

    index = build_index(args.nucleotide_db_dir)

    # Read protein sequence files and extract corresponding nucleotide sequences
    for msa_file in os.listdir(args.protein_msa_dir):
        if msa_file.endswith('.fa'):
            orthogroup_id = os.path.splitext(msa_file)[0]
            protein_msa_path = os.path.join(args.protein_msa_dir, msa_file)
            output_path = os.path.join(args.output_dir, f'{orthogroup_id}_cds.fa')

            with open(protein_msa_path, "r") as file, open(output_path, "w") as output_file:
                for line in file:
                    if line.startswith('>'):
                        header = line[1:].strip()
                        parts = header.split('_')
                        subspecies_found = False
                        for sub in args.subspecies_included:
                            if sub in header:
                                species_id = sub
                                gene_id = header.replace(species_id, '').strip('_')
                                subspecies_found = True
                                break
                        if not subspecies_found:
                            species_id = '_'.join(parts[:2])
                            gene_id = '_'.join(parts[2:])
                        identifier = f"{species_id}_{gene_id}"
                        # Check if the index contains this species and gene ID
                        if species_id in index and gene_id in index[species_id]:
                            seq_record = index[species_id][gene_id]
                            output_file.write(f">{identifier}\n{str(seq_record.seq)}\n")
                        else:
                            print(f"No valid identifier found for {msa_file} of {species_id}: {gene_id}.")
    print("CDS sequence extraction completed.")

if __name__ == '__main__':
    main()