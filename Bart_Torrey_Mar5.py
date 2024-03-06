import argparse
import numpy as np
from Bio import motifs
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from scipy.stats import norm
import random
from collections import Counter
from concurrent.futures import ThreadPoolExecutor

# Function to calculate p-value
def calculate_p_value(input_scores, background_scores):
    # Normalize match scores
    mean_background = np.mean(background_scores)
    std_background = np.std(background_scores)
    normalized_scores = (input_scores - mean_background) / std_background
    # Calculate p-value
    p_value = norm.sf(input_scores, loc=mean_background, scale=std_background)
    return p_value

# Function to generate random background sequences
def generate_randomized_background(input_file, output_file):
    # Read sequences from FASTA file
    input_sequences = []
    with open(input_file, 'r') as fasta_file:
        for record in SeqIO.parse(fasta_file, 'fasta'):
            input_sequences.append(str(record.seq))

    # Concatenate sequences
    concatenated_sequence = ''.join(input_sequences)

    # Shuffle the concatenated sequence
    shuffled_sequence = list(concatenated_sequence)
    random.shuffle(shuffled_sequence)
    shuffled_sequence = ''.join(shuffled_sequence)

    # Determine the length of random sequences based on the input sequences
    sequence_length = len(input_sequences[0])

    # Random sequences
    num_sequences = len(input_sequences)
    random_sequences = []
    for _ in range(num_sequences):
        start_index = random.randint(0, len(shuffled_sequence) - sequence_length)
        random_sequences.append(shuffled_sequence[start_index:start_index + sequence_length])

    # Prepare SeqRecord objects for each sequence
    seq_records = [SeqRecord(Seq(seq), id=f"random_seq_{i}", description="") for i, seq in enumerate(random_sequences, 1)]

    # Write sequences to a FASTA file
    with open(output_file, 'w') as output_fasta:
        SeqIO.write(seq_records, output_fasta, "fasta")

# Function to search for motifs 
def process_jaspar_motif(jaspar_motif, sequences, scores_across_motifs):
    pwm = jaspar_motif.counts.normalize(pseudocounts=0.5)
    pssm = pwm.log_odds()

    position_SNP = 10
    all_scores = {}

    for sequence_id, sequence in sequences.items():
        sequence_str = str(sequence.seq).upper()

        for strand in [sequence_str, str(Seq(sequence_str).reverse_complement())]:
            for position, score in pssm.search(strand, threshold=-9999):
                if strand != sequence_str:
                    position = len(sequence_str) - position - 1

                motif_start = position
                motif_end = position + len(pwm) - 1

                if score < 0 or motif_end < position_SNP or motif_start > position_SNP:
                    continue

                pair_key = (jaspar_motif.name, sequence_id)
                if pair_key not in all_scores or score > all_scores[pair_key][2]:
                    all_scores[pair_key] = (jaspar_motif.name, sequence_id, position, score)

    # Update the specific dictionary (input or background)
    scores_across_motifs.update(all_scores)

# Main function
def main():
    parser = argparse.ArgumentParser(description="Motif Finder Script")
    parser.add_argument("-i", "--input", required=True, help="Input sequences file (FASTA format)")
    parser.add_argument("-o", "--output", required=True, help="Output results file")
    parser.add_argument("-m", "--motif", required=True, help="Input motif file (JASPAR format)")
    parser.add_argument("-t", "--score_threshold", type=float, default=0.05, help="Top scores threshold (default = 0.05)")
    parser.add_argument("-l", "--output_limit", type=int, default=500, help="Output limit for the number of hits (default = 500)")
    args = parser.parse_args()

    # Generate random background sequences
    generate_randomized_background(args.input, "random_sequences.fasta")

    # Read input and background sequences
    input_sequences = SeqIO.to_dict(SeqIO.parse(args.input, "fasta"))
    background_sequences = SeqIO.to_dict(SeqIO.parse("random_sequences.fasta", "fasta"))

    # Separate dictionaries for input and background scores
    input_scores_across_motifs = {}
    background_scores_across_motifs = {}

    with open(args.output, 'w', encoding='utf-8') as output_file:
        output_file.write("Motif_ID\tSequence_ID\tMatch_Start\tMatch_Score\tP_Value\n")

        with open(args.motif) as f:
            motif_list = list(motifs.parse(f, "jaspar"))

        with ThreadPoolExecutor() as executor:
            for jaspar_motif in motif_list:
                # Process input sequences
                executor.submit(process_jaspar_motif, jaspar_motif, input_sequences, input_scores_across_motifs)
                # Process background sequences
                executor.submit(process_jaspar_motif, jaspar_motif, background_sequences, background_scores_across_motifs)

    # Sort input scores by match score in descending order
    sorted_input_scores = sorted(input_scores_across_motifs.values(), key=lambda x: x[3], reverse=True)

    # Select the top x% of match scores
    score_threshold_index = int(len(sorted_input_scores) * args.score_threshold)
    top_input_scores = sorted_input_scores[:score_threshold_index] if sorted_input_scores else []

    with open(args.output, 'a', encoding='utf-8') as output_file:  # Use 'a' for append
        # Calculate p-value by normalizing to background sequence scores
        for motif_id, sequence_id, position, score in top_input_scores[:min(args.output_limit, len(top_input_scores))]:
            background_sequence_scores = [item[3] for item in background_scores_across_motifs.values()]
            background_scores = np.array(background_sequence_scores)
            p_value = calculate_p_value(score, background_scores)
            
            # Write output file
            output_file.write(f"{motif_id}\t{sequence_id}\t{position}\t{score:.2f}\t{p_value:.2e}\n")

if __name__ == "__main__":
    main()
