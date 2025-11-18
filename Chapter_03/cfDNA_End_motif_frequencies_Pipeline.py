import argparse
import pysam
from collections import Counter
import pandas as pd

def is_valid_dna(sequence):
    return all(base in 'ATCG' for base in sequence)

def extract_5_prime_end_motifs(bam_file, motif_length=4):
    end_motifs = []
    with pysam.AlignmentFile(bam_file, "rb") as bam:
        for read in bam:
            # Get the 5' end sequence (motif) from the read
            end_motif_5 = read.query_sequence[:motif_length]

            # Check if the motif contains only valid DNA bases
            if is_valid_dna(end_motif_5):
                # Add the 5' end motif to the list
                end_motifs.append(end_motif_5)

    return end_motifs

def save_motif_counts_to_excel(bam_file_list, motif_length=4):
    all_motif_counts = {}
    total_motifs_per_sample = {}
    for bam_file in bam_file_list:
        sample_name = bam_file.split(".")[0]  # Extract sample name from BAM file name

        end_motifs_5_prime = extract_5_prime_end_motifs(bam_file, motif_length)
        motif_counts = Counter(end_motifs_5_prime)

        total_motifs = len(end_motifs_5_prime)
        unique_motifs = len(motif_counts)

        all_motif_counts[sample_name] = motif_counts
        total_motifs_per_sample[sample_name] = total_motifs

    # Convert motif_counts dictionary to a pandas DataFrame
    df = pd.DataFrame.from_dict(all_motif_counts, orient='index')
    df.fillna(0, inplace=True)  # Fill missing values (for motifs not present) with 0

    # Calculate the proportion for each motif in each sample
    for sample_name, total_motifs in total_motifs_per_sample.items():
        df.loc[sample_name] = df.loc[sample_name] / total_motifs * 100

    # Add total motif counts for each sample
    df['Total_Count'] = df.sum(axis=1)

    # Create an Excel writer
    excel_file = "motif_counts.xlsx"
    writer = pd.ExcelWriter(excel_file, engine='xlsxwriter')

    # Save the DataFrame to an Excel file in one sheet
    df.to_excel(writer, sheet_name='Motif_Counts', index_label='Sample')

    # Close the Excel writer
    writer.save()

    print(f"Results saved in {excel_file}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Generate motif counts from BAM files and save to Excel.')
    parser.add_argument('sample_list_file', type=argparse.FileType('r'), help='Text file containing a list of BAM files')
    parser.add_argument('-l', '--motif_length', type=int, default=4, help='Length of motifs to extract')

    args = parser.parse_args()
    bam_file_list = [line.strip() for line in args.sample_list_file]
    motif_length = args.motif_length

    save_motif_counts_to_excel(bam_file_list, motif_length)

