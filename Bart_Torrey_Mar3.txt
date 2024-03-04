import argparse
from Bio import motifs
from Bio import SeqIO
from scipy.stats import hypergeom
import random
from collections import Counter
from concurrent.futures import ThreadPoolExecutor
from decimal import Decimal
import argparse

# Function to calculate p-value (FIX THIS)
def calculate_p_value(matches, total_sequences, total_background_sequences):
    M = len(matches)  
    n = total_sequences
    N = total_background_sequences
    x = len(set(matches))
    p_value = hypergeom.cdf(x - 1, N, n - N, M)
    return p_value

# Function to generate random background sequences
def generate_randomized_background(input_sequences):
    concatenated_input = ''.join(str(seq) for seq in input_sequences)
    nucleotide_counts = Counter(concatenated_input)
    total_nucleotides = len(concatenated_input)
    nucleotide_probs = {nuc: count / total_nucleotides for nuc, count in nucleotide_counts.items()}
    random_sequence = ''.join(random.choices(list(nucleotide_probs.keys()), weights=list(nucleotide_probs.values()), k=total_nucleotides))
    background_sequences = [random_sequence[i:i+len(seq)] for i, seq in enumerate(input_sequences)]
    return background_sequences

# Function to search for motifs
def process_jaspar_motif(jaspar_motif, input_sequences, background_sequences, output_file):
    pwm = jaspar_motif.counts.normalize(pseudocounts=0.5)
    pssm = pwm.log_odds()

    # specify the position of interest (SNP = 11th nucleotide, adjust as needed)
    position_of_interest = 11 - 1

    # collect all scores with associated p-values
    all_scores = {}

    for sequence_id, sequence in input_sequences.items():
        sequence_str = str(sequence.seq).upper()

        for position, score in pssm.search(sequence_str, threshold=-9999):
            motif_start = position
            motif_end = position + len(pwm) - 1

            if score < 0 or motif_end < position_of_interest or motif_start > position_of_interest:
                continue

            p_value = calculate_p_value([position], len(input_sequences), len(background_sequences))

            # keep track of the highest score per sequence-motif pair
            pair_key = (sequence_id, jaspar_motif.name)
            if pair_key not in all_scores or score > all_scores[pair_key][2]:
                all_scores[pair_key] = (sequence_id, position, score, p_value)

    # sort all_scores by score in descending order
    sorted_scores = sorted(all_scores.values(), key=lambda x: x[2], reverse=True)

    # select the top % based on score (adjust as needed)
    score_threshold_index = int(len(sorted_scores) * 0.1)
    selected_scores = sorted_scores[:score_threshold_index] if sorted_scores else []

    # sort selected_scores by p-value in ascending order
    selected_scores.sort(key=lambda x: x[3])

    # output the top matches based on p-value (max = 100)
    for sequence_id, position, score, p_value in selected_scores[:min(100, len(selected_scores))]:
        formatted_p_value = Decimal(p_value).__format__('e')
        output_file.write(f"{jaspar_motif.name}\t{sequence_id}\t{position}\t{score}\t{formatted_p_value}\n")

def limit_output_lines(output_file_path, limit):
    with open(output_file_path, 'r', encoding='utf-8') as output_file:
        lines = output_file.readlines()

    if len(lines) > limit:
        with open(output_file_path, 'w', encoding='utf-8') as output_file:
            output_file.write("".join(lines[:limit]))

# Command line function            
def main():
    parser = argparse.ArgumentParser(description="Motif Finder Script")
    parser.add_argument("-i", "--input", required=True, help="Input sequences file (FASTA format)")
    parser.add_argument("-o", "--output", required=True, help="Output results file")
    parser.add_argument("-m", "--motif", required=True, help="Motif file (JASPAR format)")
    args = parser.parse_args()

    input_sequences = SeqIO.to_dict(SeqIO.parse(args.input, "fasta"))
    background_sequences = generate_randomized_background(input_sequences.values())

    with open(args.output, 'w', encoding='utf-8') as output_file:
        output_file.write("Motif_ID\tSequence_ID\tMatch_Start\tMatch_Score\tP_Value\n")

        with open(args.motif) as f:
            motif_list = list(motifs.parse(f, "jaspar"))

        with ThreadPoolExecutor() as executor:
            for jaspar_motif in motif_list:
                executor.submit(process_jaspar_motif, jaspar_motif, input_sequences, background_sequences, output_file)
    
    limit_output_lines(args.output, 100)

if __name__ == "__main__":
    main()
