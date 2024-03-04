import argparse
from Bio import motifs
from Bio import SeqIO
from scipy.stats import hypergeom
import random
from collections import Counter
from concurrent.futures import ThreadPoolExecutor
from decimal import Decimal
from Bio.Seq import Seq

# Function to calculate p-value (FIX)
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
def process_jaspar_motif(jaspar_motif, input_sequences, background_sequences, all_scores_across_motifs):
    pwm = jaspar_motif.counts.normalize(pseudocounts=0.5)
    pssm = pwm.log_odds()

    # specify the position of the SNP
    position_SNP = 10

    # store all sequence-motif scores
    all_scores = {}

    for sequence_id, sequence in input_sequences.items():
        sequence_str = str(sequence.seq).upper()

        for strand in [sequence_str, str(Seq(sequence_str).reverse_complement())]:
            for position, score in pssm.search(strand, threshold=-9999):
                if strand != sequence_str: 
                    position = len(sequence_str) - position - 1

                motif_start = position
                motif_end = position + len(pwm) - 1

                if score < 0 or motif_end < position_SNP or motif_start > position_SNP:
                    continue

                p_value = calculate_p_value([position], len(input_sequences), len(background_sequences))

                # store the highest score per sequence-motif pair
                pair_key = (jaspar_motif.name, sequence_id)
                if pair_key not in all_scores or score > all_scores[pair_key][2]:
                    all_scores[pair_key] = (jaspar_motif.name, sequence_id, position, score, p_value)

    # update
    all_scores_across_motifs.update(all_scores)

# Main function
def main():
    parser = argparse.ArgumentParser(description="Motif Finder Script")
    parser.add_argument("-i", "--input", required=True, help="Input sequences file (FASTA format)")
    parser.add_argument("-o", "--output", required=True, help="Output results file")
    parser.add_argument("-m", "--motif", required=True, help="Motif file (JASPAR format)")
    parser.add_argument("-t", "--score_threshold", type=float, default=0.05, help="Top scores threshold")
    parser.add_argument("-l", "--output_limit", type=int, default=500, help="Output limit for the number of hits")
    args = parser.parse_args()

    input_sequences = SeqIO.to_dict(SeqIO.parse(args.input, "fasta"))
    background_sequences = generate_randomized_background(input_sequences.values())

    with open(args.output, 'w', encoding='utf-8') as output_file:
        output_file.write("Motif_ID\tSequence_ID\tMatch_Start\tMatch_Score\tP_Value\n")

        with open(args.motif) as f:
            motif_list = list(motifs.parse(f, "jaspar"))

        all_scores_across_motifs = {}
        with ThreadPoolExecutor() as executor:
            for jaspar_motif in motif_list:
                executor.submit(process_jaspar_motif, jaspar_motif, input_sequences, background_sequences, all_scores_across_motifs)

        # Sort all_scores_across_motifs by match score in descending order
        sorted_scores = sorted(all_scores_across_motifs.values(), key=lambda x: x[3], reverse=True)

        # Select the top x% of match scores
        score_threshold_index = int(len(sorted_scores) * args.score_threshold)
        top_scores = sorted_scores[:score_threshold_index] if sorted_scores else []

        # Sort top_scores by p-value in ascending order
        # top_scores.sort(key=lambda x: x[4])

        # Output the top hits ordered by match score (p-value if possible)
        for motif_id, sequence_id, position, score, p_value in top_scores[:min(args.output_limit, len(top_scores))]:
            formatted_p_value = Decimal(p_value).__format__('e')
            output_file.write(f"{motif_id}\t{sequence_id}\t{position}\t{score}\t{formatted_p_value}\n")

if __name__ == "__main__":
    main()
