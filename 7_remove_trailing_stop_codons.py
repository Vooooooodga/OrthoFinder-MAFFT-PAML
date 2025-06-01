import os
import argparse
from Bio import SeqIO
from Bio.Seq import Seq

def remove_trailing_stop_codon(sequence_str):
    """
    Checks if the trailing codon is a stop codon and removes it.
    Standard stop codons: TAA, TAG, TGA.
    """
    stop_codons = ["TAA", "TAG", "TGA"]
    seq_len = len(sequence_str)

    if seq_len >= 3:
        last_codon = sequence_str[-3:].upper()
        if last_codon in stop_codons:
            return sequence_str[:-3]
    return sequence_str

def process_fasta_file(input_file_path, output_file_path):
    """
    Reads a FASTA file, removes trailing stop codons from each sequence,
    and writes the modified sequences to a new FASTA file.
    """
    records_to_write = []
    try:
        for record in SeqIO.parse(input_file_path, "fasta"):
            original_seq_str = str(record.seq)
            modified_seq_str = remove_trailing_stop_codon(original_seq_str)
            
            # Create a new SeqRecord with the modified sequence
            # Retain the original ID and description
            modified_record = record[:] # Create a shallow copy to modify sequence
            modified_record.seq = Seq(modified_seq_str)
            records_to_write.append(modified_record)
        
        SeqIO.write(records_to_write, output_file_path, "fasta")
        # print(f"Processed: {input_file_path} -> {output_file_path}")
        return True
    except Exception as e:
        print(f"Error processing file {input_file_path}: {e}")
        return False

def main():
    parser = argparse.ArgumentParser(
        description="Removes trailing stop codons (TAA, TAG, TGA) from sequences in FASTA files.",
        formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument(
        "--input_dir", 
        required=True, 
        help="Directory containing input FASTA files (codon alignments)."
    )
    parser.add_argument(
        "--output_dir", 
        required=True, 
        help="Directory to save FASTA files with trailing stop codons removed."
    )
    
    args = parser.parse_args()

    if not os.path.isdir(args.input_dir):
        print(f"Error: Input directory '{args.input_dir}' not found.")
        return

    if not os.path.exists(args.output_dir):
        os.makedirs(args.output_dir)
        print(f"Created output directory: {args.output_dir}")

    print(f"Starting to process files from: {args.input_dir}")
    print(f"Output will be saved to: {args.output_dir}")

    processed_files_count = 0
    failed_files_count = 0

    for filename in os.listdir(args.input_dir):
        if filename.endswith((".fa", ".fasta", ".aln")): # Adjust extensions as needed
            input_file_path = os.path.join(args.input_dir, filename)
            output_file_path = os.path.join(args.output_dir, filename) # Keep the same filename

            if os.path.isfile(input_file_path):
                print(f"Processing {filename}...")
                if process_fasta_file(input_file_path, output_file_path):
                    processed_files_count += 1
                else:
                    failed_files_count +=1
    
    print(f"\\n--- Summary ---")
    print(f"Total files processed: {processed_files_count}")
    if failed_files_count > 0:
        print(f"Files failed to process: {failed_files_count}")
    print("Script execution finished.")

if __name__ == '__main__':
    main() 