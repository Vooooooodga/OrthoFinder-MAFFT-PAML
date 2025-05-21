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
    args = parser.parse_args()

    # 创建输出目录（如果不存在）
    if not os.path.exists(args.output_dir):
        os.makedirs(args.output_dir)

    # 读取基因列表
    with open(args.gene_list_file, newline='') as csvfile:
        gene_reader = csv.reader(csvfile)
        gene_list = {row[0] for row in gene_reader}

    # 读取Orthogroups文件并建立基因名到Orthogroups的映射
    gene_to_orthogroup = {}
    with open(args.orthogroups_file, newline='') as csvfile:
        tsv_reader = csv.reader(csvfile, delimiter='\t')
        header = next(tsv_reader)  # 跳过表头
        for row in tsv_reader:
            orthogroup = row[0]
            # Assuming Apis_mellifera genes are in the 6th column (index 5)
            if len(row) > 5:
                apis_mellifera_genes = row[5]  # Apis_mellifera基因列
                if apis_mellifera_genes:
                    genes = apis_mellifera_genes.split(', ')
                    for gene in genes:
                        if gene in gene_list:
                            gene_to_orthogroup[gene] = orthogroup

    # 确认找到的Orthogroups
    found_orthogroups = set(gene_to_orthogroup.values())

    # 遍历MSA文件目录，找到匹配的文件
    for msa_file in os.listdir(args.msa_dir):
        if msa_file.split('.')[0] in found_orthogroups:
            source_file = os.path.join(args.msa_dir, msa_file)
            destination_file = os.path.join(args.output_dir, msa_file)
            shutil.copy(source_file, destination_file)
            print(f'Copied {msa_file} to {args.output_dir}')

    print("MSA files selection completed.")

if __name__ == '__main__':
    main()