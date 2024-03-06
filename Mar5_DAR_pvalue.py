import argparse
from Bio import motifs
from Bio import SeqIO
from scipy.stats import hypergeom
import random
from collections import Counter
from concurrent.futures import ThreadPoolExecutor
from decimal import Decimal
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import numpy as np
from scipy.stats import norm

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

    # Generate random sequences
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
def process_jaspar_motif(jaspar_motif, input_sequences, all_scores_across_motifs):
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
                    
                #p_value = calculate_p_value([position], len(input_sequences), len(background_sequences))


                # store the highest score per sequence-motif pair
                pair_key = (jaspar_motif.name, sequence_id)
                if pair_key not in all_scores or score > all_scores[pair_key][2]:
                    all_scores[pair_key] = (jaspar_motif.name, sequence_id, position, score)

    # update
    all_scores_across_motifs.update(all_scores)


    # function to search for motifs IN BACKGROUND SEQS
def process_jaspar_motif_BS(jaspar_motif, background_sequences, all_scores_across_motifs):
    pwm = jaspar_motif.counts.normalize(pseudocounts=0.5)
    pssm = pwm.log_odds()

    # specify the position of the SNP
    position_SNP = 10

    # store all sequence-motif scores
    all_scores = {}

    for sequence_id, sequence in background_sequences.items():
        sequence_str = str(sequence.seq).upper()
        #print(sequence_str)

        for strand in [sequence_str, str(Seq(sequence_str).reverse_complement())]:
            for position, score in pssm.search(strand, threshold=-9999):
                if strand != sequence_str: 
                    position = len(sequence_str) - position - 1

                motif_start = position
                motif_end = position + len(pwm) - 1

                if score < 0 or motif_end < position_SNP or motif_start > position_SNP:
                    continue               

                # store the highest score per sequence-motif pair
                pair_key = (jaspar_motif.name, sequence_id)
                if pair_key not in all_scores or score > all_scores[pair_key][2]:
                    all_scores[pair_key] = (jaspar_motif.name, sequence_id, position, score)

    # update
    all_scores_across_motifs.update(all_scores)    

    
# Main function
def main():
    parser = argparse.ArgumentParser(description="Motif Finder Script")
    parser.add_argument("-i", "--input", required=True, help="Input sequences file (FASTA format)")
    parser.add_argument("-o", "--output", required=True, help="Output results file")
    parser.add_argument("-m", "--motif", required=True, help="Input motif file (JASPAR format)")
    parser.add_argument("-t", "--score_threshold", type=float, default=0.05, help="Top scores threshold (default = 0.05)")
    parser.add_argument("-l", "--output_limit", type=int, default=500, help="Output limit for the number of hits (default = 500)")
    args = parser.parse_args()
    
    
    # Specify the path to the output file where you want to save the sequences
    generate_randomized_background(args.input, "random_sequences.fasta")

    input_sequences = SeqIO.to_dict(SeqIO.parse(args.input, "fasta"))
    background_sequences = SeqIO.to_dict(SeqIO.parse("random_sequences.fasta", "fasta"))

########background seqs
    with open(args.motif) as f:
        motif_list = list(motifs.parse(f, "jaspar"))        
        
    all_scores_across_motifs_BS = {}
    with ThreadPoolExecutor() as executor:
        for jaspar_motif in motif_list:
            executor.submit(process_jaspar_motif_BS, jaspar_motif, background_sequences, all_scores_across_motifs_BS)
                
        # Sort all_scores_across_motifs by match score in descending order
    sorted_scores = sorted(all_scores_across_motifs_BS.values(), key=lambda x: x[3], reverse=True)

        # Select the top x% of match scores
    score_threshold_index = int(len(sorted_scores) * args.score_threshold)
    top_scores = sorted_scores[:score_threshold_index] if sorted_scores else []
        
        # Output the top hits ordered by match score (p-value if possible)
    all_background_scores_list = []
    for motif_id, sequence_id, position, score in top_scores[:min(args.output_limit, len(top_scores))]:
        all_background_scores_list.append(score)
            
    mean_background = np.mean(all_background_scores_list)
    std_background = np.std(all_background_scores_list)
        

########input seqs       
    with open(args.motif) as f:
        motif_list = list(motifs.parse(f, "jaspar"))

    all_scores_across_motifs = {}
    with ThreadPoolExecutor() as executor:
        for jaspar_motif in motif_list:
            executor.submit(process_jaspar_motif, jaspar_motif, input_sequences, all_scores_across_motifs)
                
        # Sort all_scores_across_motifs by match score in descending order
    sorted_scores = sorted(all_scores_across_motifs.values(), key=lambda x: x[3], reverse=True)

        # Select the top x% of match scores
    score_threshold_index = int(len(sorted_scores) * args.score_threshold)
    top_scores = sorted_scores[:score_threshold_index] if sorted_scores else []

        # Sort top_scores by p-value in ascending order
        
    all_input_scores_list = []
    for motif_id, sequence_id, position, score in top_scores[:min(args.output_limit, len(top_scores))]:
        all_input_scores_list.append(score)
    
    normalized_scores = (all_input_scores_list - mean_background) / std_background
   
    with open(args.output, 'w', encoding='utf-8') as output_file:
        output_file.write("Motif_ID\tSequence_ID\tMatch_Start\tMatch_Score\tP_Value\n")

        for motif_id, sequence_id, position, score in top_scores[:min(args.output_limit, len(top_scores))]:
                # Calculate p-value for each score
            p_value = norm.sf(score, loc=mean_background, scale=std_background)
                # Write the output line including the calculated p-value
            output_file.write(f"{motif_id}\t{sequence_id}\t{position}\t{score}\t{p_value:.2e}\n")


if __name__ == "__main__":
    main()
