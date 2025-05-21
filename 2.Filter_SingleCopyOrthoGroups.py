import os
import shutil
import argparse

def main():
    parser = argparse.ArgumentParser(description="Filter and copy single-copy orthogroup files.")
    parser.add_argument("--single_copy_orthogroups_file", default='/home/yuhangjia/data/AlternativeSplicing/exon_expansion_test/orthofinder/output_msa_2/Results_Jul05/Orthogroups/Orthogroups_SingleCopyOrthologues.txt', help="Path to the single-copy orthogroups file.")
    parser.add_argument("--selected_msa_files_dir", default='./selected_msa_files_deg', help="Directory containing selected MSA files.")
    parser.add_argument("--output_dir", default='single_copy_selected_msa_files_deg', help="Directory to save filtered single-copy MSA files.")
    args = parser.parse_args()

    # Create output directory if it doesn't exist
    if not os.path.exists(args.output_dir):
        os.makedirs(args.output_dir)

    # Read list of single-copy orthogroups
    with open(args.single_copy_orthogroups_file, 'r') as file:
        single_copy_orthogroups = {line.strip() for line in file}

    # Filter and copy single-copy orthogroup files
    for msa_file in os.listdir(args.selected_msa_files_dir):
        orthogroup_id = os.path.splitext(msa_file)[0]  # Remove file extension
        if orthogroup_id in single_copy_orthogroups:
            source_file = os.path.join(args.selected_msa_files_dir, msa_file)
            destination_file = os.path.join(args.output_dir, msa_file)
            shutil.copy(source_file, destination_file)
            print(f'Copied {msa_file} to {args.output_dir}')

    print("Single-copy orthogroups selection completed.")

if __name__ == '__main__':
    main()

