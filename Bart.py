import os
import argparse
import numpy as np
from Bio import motifs, SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from scipy.stats import norm 
import random
from collections import Counter
from concurrent.futures import ThreadPoolExecutor

# Function to calculate p-value
def calculate_p_value(input_scores, background_scores):
    # normalize input scores to background score distribution
    mean_background = np.mean(background_scores)
    std_background = np.std(background_scores)
    normalized_scores = (input_scores - mean_background) / std_background
    
    # calculate p-value
    p_value = norm.sf(input_scores, loc=mean_background, scale=std_background)
    return p_value

# Function to generate random background sequences of the same length / composition as input sequences
def generate_randomized_background(input_file, output_file):
    # read input sequences from FASTA file
    input_sequences = []
    with open(input_file, 'r') as fasta_file:
        for record in SeqIO.parse(fasta_file, 'fasta'):
            input_sequences.append(str(record.seq))

    # determine length of input sequences
    sequence_length = len(input_sequences[0])

    # calculate nucleotide composition of input sequences
    composition = {
        'A': sum(seq.count('A') for seq in input_sequences),
        'T': sum(seq.count('T') for seq in input_sequences),
        'C': sum(seq.count('C') for seq in input_sequences),
        'G': sum(seq.count('G') for seq in input_sequences),
    }

    # set probabilities for random sequences
    total_nucleotides = sum(composition.values())
    probabilities = {base: count / total_nucleotides for base, count in composition.items()}

    # generate random sequences
    num_sequences = len(input_sequences)
    random_sequences = []
    for _ in range(num_sequences):
        random_sequence = ''.join(random.choices(list(probabilities.keys()), weights=list(probabilities.values()), k=sequence_length))
        random_sequences.append(random_sequence)

    # write the random sequences to the output file
    with open(output_file, 'w') as output_fasta:
        for i, seq in enumerate(random_sequences):
            output_fasta.write(f'>random_sequence_{i}\n{seq}\n')

# Function to search for motifs 
def process_jaspar_motif(jaspar_motif, sequences, scores_across_motifs):
    pwm = jaspar_motif.counts.normalize(pseudocounts=0.5)
    pssm = pwm.log_odds()

    # set position of SNP
    position_SNP = 10

    # collect all scores
    all_scores = {}

    for sequence_id, sequence in sequences.items():
        sequence_str = str(sequence.seq).upper()

        # reverse complement sequences
        for strand in [sequence_str, str(Seq(sequence_str).reverse_complement())]:
            for position, score in pssm.search(strand, threshold=-9999):
                if strand != sequence_str:
                    position = len(sequence_str) - position - 1

                # set positions
                motif_start = position
                motif_end = position + len(pwm) - 1

                # filter out negative scores and motifs that don't overlap the SNP
                if score < 0 or motif_end < position_SNP or motif_start > position_SNP:
                    continue

                # store the highest score for each motif-sequence pair
                pair_key = (jaspar_motif.name, sequence_id)
                if pair_key not in all_scores or score > all_scores[pair_key][2]:
                    all_scores[pair_key] = (jaspar_motif.name, sequence_id, position, score)

    # update the specific dictionary (input or background)
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

    # generate random background sequences
    generate_randomized_background(args.input, "random_sequences.fasta")

    # read input and background sequences
    input_sequences = SeqIO.to_dict(SeqIO.parse(args.input, "fasta"))
    background_sequences = SeqIO.to_dict(SeqIO.parse("random_sequences.fasta", "fasta"))

    # dictionaries for input and background scores
    input_scores_across_motifs = {}
    background_scores_across_motifs = {}

    with open(args.output, 'w', encoding='utf-8') as output_file:
        output_file.write("Motif_ID\tSequence_ID\tMatch_Start\tMatch_Score\tP_Value\n")

        with open(args.motif) as f:
            motif_list = list(motifs.parse(f, "jaspar"))

        with ThreadPoolExecutor() as executor:
            for jaspar_motif in motif_list:
                # process input sequences
                executor.submit(process_jaspar_motif, jaspar_motif, input_sequences, input_scores_across_motifs)
                
                # process background sequences
                executor.submit(process_jaspar_motif, jaspar_motif, background_sequences, background_scores_across_motifs)

    # sort input scores by match score in descending order
    sorted_input_scores = sorted(input_scores_across_motifs.values(), key=lambda x: x[3], reverse=True)

    # select the top x% of match scores
    score_threshold_index = int(len(sorted_input_scores) * args.score_threshold)
    top_input_scores = sorted_input_scores[:score_threshold_index] if sorted_input_scores else []

    with open(args.output, 'a', encoding='utf-8') as output_file:
        
        # calculate p-value of top match scores
        for motif_id, sequence_id, position, score in top_input_scores[:min(args.output_limit, len(top_input_scores))]:
            background_sequence_scores = [item[3] for item in background_scores_across_motifs.values()]
            background_scores = np.array(background_sequence_scores)
            p_value = calculate_p_value(score, background_scores)
            
            # write output file
            output_file.write(f"{motif_id}\t{sequence_id}\t{position}\t{score:.2f}\t{p_value:.2e}\n")

if __name__ == "__main__":
    main()

    # remove the random sequences file
    os.remove("random_sequences.fasta")

