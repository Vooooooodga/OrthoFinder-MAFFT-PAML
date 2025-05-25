import os
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import argparse

def translate_cds_to_protein(input_folder, output_folder, genetic_code=1):
    # Ensure output folder exists
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)
    # Iterate through all files in the input folder
    for filename in os.listdir(input_folder):
        if filename.endswith(".fa") or filename.endswith(".fasta"): # Ensure FASTA files are processed
            input_path = os.path.join(input_folder, filename)
            output_filename = filename.replace(".fasta", "_protein.fasta").replace(".fa", "_protein.fa")
            output_path = os.path.join(output_folder, output_filename)
            # Read CDS sequences, translate to protein, and save
            with open(output_path, 'w') as output_file:
                protein_records = []
                for record in SeqIO.parse(input_path, "fasta"):
                    original_seq_len = len(record.seq)
                    protein_seq = record.seq.translate(table=genetic_code, stop_symbol='X')

                    if original_seq_len % 3 != 0:
                        print(f"详细警告: 文件 {filename} (Orthogroup) 中的序列 {record.id} 长度为 {original_seq_len}, 不是3的倍数。翻译结果末尾将添加一个 'X' 代表不完整或未知的氨基酸。")
                        protein_seq += Seq("X")
                    
                    protein_record = SeqRecord(protein_seq, id=record.id, description="")
                    protein_records.append(protein_record)
                # Write translated protein sequences to a new file
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