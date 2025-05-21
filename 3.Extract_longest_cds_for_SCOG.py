from Bio import SeqIO
import os
import argparse

# 定义文件路径

protein_msa_dir = 'single_copy_selected_msa_files_deg'

nucleotide_db_dir = '../longest_cds'

output_dir = 'cds_sequences_deg'

if not os.path.exists(output_dir):

    os.makedirs(output_dir)

# 物种具有亚种名的列表

subspecies_included = ['Bombus_vancouverensis_nearcticus', 'Osmia_bicornis_bicornis']

# 自定义解析函数，从描述中提取最合适的基因标识符

def parse_description(description):

    return description.strip('>')

# 创建索引，特别处理基因ID的解析

def build_index(nucleotide_db_dir):

    index = {}

    for nucleotide_file in os.listdir(nucleotide_db_dir):

        file_path = os.path.join(nucleotide_db_dir, nucleotide_file)

        species_id = os.path.splitext(nucleotide_file)[0]

        record_dict = {}

        for record in SeqIO.parse(file_path, "fasta"):

            key = parse_description(record.description)

            # 检查并选择最长的转录本

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
    parser.add_argument("--subspecies_included", nargs='*', default=['Bombus_vancouverensis_nearcticus', 'Osmia_bicornis_bicornis'], help="List of species with subspecies names.")
    args = parser.parse_args()

    if not os.path.exists(args.output_dir):
        os.makedirs(args.output_dir)

    index = build_index(args.nucleotide_db_dir)

    # 读取蛋白质序列文件并提取对应的核酸序列

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

                        # 检查索引是否包含此物种和基因ID

                        if species_id in index and gene_id in index[species_id]:

                            seq_record = index[species_id][gene_id]

                            output_file.write(f">{identifier}\n{str(seq_record.seq)}\n")

                        else:

                            print(f"No valid identifier found for {msa_file} of {species_id}: {gene_id}.")

    print("CDS sequence extraction completed.")

if __name__ == '__main__':
    main()