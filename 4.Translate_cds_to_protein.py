import os
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import argparse

def translate_cds_to_protein(input_folder, output_folder, genetic_code=1):
    # 确保输出文件夹存在
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)
    # 遍历输入文件夹中的所有文件
    for filename in os.listdir(input_folder):
        if filename.endswith(".fa") or filename.endswith(".fasta"):
            input_path = os.path.join(input_folder, filename)
            output_filename = filename.replace(".fasta", "_protein.fasta").replace(".fa", "_protein.fa")
            output_path = os.path.join(output_folder, output_filename)
            # 读取CDS序列，翻译成蛋白质，并保存
            with open(output_path, 'w') as output_file:
                protein_records = []
                for record in SeqIO.parse(input_path, "fasta"):
                    protein_seq = record.seq.translate(table=genetic_code, to_stop=True)
                    protein_record = SeqRecord(protein_seq, id=record.id, description="")
                    protein_records.append(protein_record)
                # 将翻译后的蛋白质序列写入新文件
                SeqIO.write(protein_records, output_file, "fasta")

def main():
    parser = argparse.ArgumentParser(description="Translate CDS sequences to protein sequences.")
    parser.add_argument("--input_folder", default="cds_sequences_deg", help="Folder containing CDS FASTA files.")
    parser.add_argument("--output_folder", default="translated_proteins_deg", help="Folder to save translated protein FASTA files.")
    parser.add_argument("--genetic_code", type=int, default=1, help="Genetic code table number (default: 1 for standard code).")
    args = parser.parse_args()

    translate_cds_to_protein(args.input_folder, args.output_folder, args.genetic_code)
    print(f"Protein translation completed. Output saved to {args.output_folder}")

if __name__ == '__main__':
    main()