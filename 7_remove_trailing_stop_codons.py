import os
import argparse
from Bio import SeqIO
from Bio.Seq import Seq

def replace_internal_stop_codons(sequence_str):
    """
    Iterates through the sequence codon by codon.
    If an internal stop codon (TAA, TAG, TGA, case-insensitive) is found,
    it's replaced with "NNN". Trailing stop codons are NOT modified.
    """
    stop_codons = ["TAA", "TAG", "TGA"]
    seq_upper = sequence_str.upper()
    original_chars = list(sequence_str)
    seq_len = len(seq_upper)
    
    # Determine the end index for the loop to cover only full codons
    # This also helps identify the start index of the last potential full codon
    last_codon_loop_end = seq_len - (seq_len % 3)

    for i in range(0, last_codon_loop_end, 3):
        codon = seq_upper[i:i+3]
        if codon in stop_codons:
            # Check if this stop codon is the last complete codon in the sequence
            is_trailing_stop_codon = (i + 3 == last_codon_loop_end)
            
            if not is_trailing_stop_codon:
                # It's an internal stop codon, so replace it
                original_chars[i] = 'N'
                original_chars[i+1] = 'N'
                original_chars[i+2] = 'N'
            # If it is_trailing_stop_codon, do nothing
            
    return "".join(original_chars)

def process_fasta_file(input_file_path, output_file_path):
    """
    Reads a FASTA file, replaces internal stop codons (excluding trailing ones)
    with 'NNN' in each sequence, and writes the modified sequences to a new FASTA file.
    """
    records_to_write = []
    try:
        for record in SeqIO.parse(input_file_path, "fasta"):
            original_seq_str = str(record.seq)
            modified_seq_str = replace_internal_stop_codons(original_seq_str)
            
            modified_record = record[:] 
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
        description="Replaces internal stop codons (TAA, TAG, TGA) with 'NNN' in FASTA sequences, checking codon by codon. Trailing stop codons are preserved.",
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