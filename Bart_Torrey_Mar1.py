import argparse
from Bio import motifs
from Bio import SeqIO
from scipy.stats import hypergeom
import random
from collections import Counter
from concurrent.futures import ThreadPoolExecutor

def calculate_p_value(matches, total_sequences, total_background_sequences):
    M = len(matches)
    n = total_sequences
    N = total_background_sequences
    x = len(set(matches))
    p_value = hypergeom.sf(x - 1, N, M, n)
    return p_value

def generate_randomized_background(input_sequences):
    concatenated_input = ''.join(str(seq) for seq in input_sequences)
    nucleotide_counts = Counter(concatenated_input)
    total_nucleotides = len(concatenated_input)
    nucleotide_probs = {nuc: count / total_nucleotides for nuc, count in nucleotide_counts.items()}
    random_sequence = ''.join(random.choices(list(nucleotide_probs.keys()), weights=list(nucleotide_probs.values()), k=total_nucleotides))
    background_sequences = [random_sequence[i:i+len(seq)] for i, seq in enumerate(input_sequences)]
    return background_sequences

def process_jaspar_motif(jaspar_motif, input_sequences, background_sequences, output_file):
    pwm = jaspar_motif.counts.normalize(pseudocounts=0.5)
    pssm = pwm.log_odds()

    # Collect all scores
    all_scores = []
    for sequence_id, sequence in input_sequences.items():
        sequence_str = str(sequence.seq).upper()
        for position, score in pssm.search(sequence_str, threshold=-9999):
            if score < 0:
                continue
            all_scores.append((sequence_id, position, score))

    # Calculate the score threshold for the top 5%
    score_threshold = sorted(all_scores, key=lambda x: x[2], reverse=True)[:int(len(all_scores) * 0.05)][-1][2] if all_scores else 0

    # Process sequences with scores above the threshold
    for sequence_id, position, score in all_scores:
        if score > score_threshold:
            p_value = calculate_p_value([position], len(input_sequences), len(background_sequences))
            output_file.write(f"{jaspar_motif.name}\t{sequence_id}\t{position}\t{score}\t{p_value}\n")

def limit_output_lines(output_file_path, limit):
    with open(output_file_path, 'r', encoding='utf-8') as output_file:
        lines = output_file.readlines()

    if len(lines) > limit:
        print(f"Output file has more than {limit} lines. Truncating to the top {limit} lines.")
        with open(output_file_path, 'w', encoding='utf-8') as output_file:
            output_file.write("".join(lines[:limit]))

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

    # Limit the output to the top 500 lines
    limit_output_lines(args.output, 500)

if __name__ == "__main__":
    main()
