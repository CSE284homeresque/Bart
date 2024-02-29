import argparse
from Bio.Seq import Seq
from Bio import motifs
from Bio import SeqIO
from Bio.SeqUtils import nt_search
from scipy.stats import hypergeom
import random
from collections import Counter

# Function to calculate hypergeometric p-value
def calculate_p_value(matches, total_sequences, total_background_sequences):
    M = len(matches)
    n = total_sequences
    N = total_background_sequences
    x = len(set(matches))
    
    p_value = hypergeom.sf(x - 1, N, M, n)
    return p_value

# Function to get the reverse complement of a sequence
def reverse_complement(sequence):
    return str(Seq(sequence).reverse_complement())

# Function to generate randomized background sequences
def generate_randomized_background(input_sequences):
    # Concatenate all input sequences into one string
    concatenated_input = ''.join(str(seq) for seq in input_sequences)

    # Calculate nucleotide frequencies in the concatenated input
    nucleotide_counts = Counter(concatenated_input)

    # Calculate nucleotide probabilities
    total_nucleotides = len(concatenated_input)
    nucleotide_probs = {nuc: count / total_nucleotides for nuc, count in nucleotide_counts.items()}

    # Generate a random sequence with the same length and nucleotide composition
    random_sequence = ''.join(random.choices(list(nucleotide_probs.keys()), 
                                             weights=list(nucleotide_probs.values()), 
                                             k=total_nucleotides))

    # Split the random sequence back into individual sequences
    background_sequences = [Seq(random_sequence[i:i+len(seq)]) for i, seq in enumerate(input_sequences)]

    return background_sequences

def main():
    # Create argument parser
    parser = argparse.ArgumentParser(description="Motif Finder Script")
    
    # Add arguments
    parser.add_argument("-i", "--input", required=True, help="Input sequences file (FASTA format)")
    parser.add_argument("-o", "--output", required=True, help="Output results file")
    parser.add_argument("-m", "--motif", required=True, help="Motif file (JASPAR format)")
    parser.add_argument("-p", "--pvalue", type=float, default=0.05, help="P-value threshold (default: 0.05)")

    # Parse command-line arguments
    args = parser.parse_args()

    # Load input sequences
    input_sequences = SeqIO.to_dict(SeqIO.parse(args.input, "fasta"))

    # Generate randomized background sequences
    background_sequences = generate_randomized_background(input_sequences.values())

    # Open output file for writing results
    with open(args.output, 'w') as output_file:
        # Write header
        output_file.write("Motif_ID\tSequence_ID\tMatch_Positions\tP_Value\n")

        # Container for results
        results = []

        # Load JASPAR motif
        jaspar_motif = motifs.read(open(args.motif), "jaspar")

        # Loop over input sequences
        for sequence_id, sequence in input_sequences.items():
            # Perform motif scan on the original sequence
            matches = nt_search(str(jaspar_motif), str(sequence))

            # Perform motif scan on the reverse complement
            rev_comp_matches = nt_search(str(jaspar_motif), reverse_complement(str(sequence)))

            # Combine matches from both orientations
            all_matches = matches + [(pos, "rev_comp") for pos in rev_comp_matches]

            # Calculate hypergeometric p-value
            p_value = calculate_p_value(all_matches, len(input_sequences), len(input_sequences))

            # Append results to container
            results.append((jaspar_motif.name, sequence_id, str(all_matches), p_value))

        # Sort results by p-value (most significant first)
        results.sort(key=lambda x: x[3])

        # Write sorted results to the output file
        for result in results:
            output_file.write(f"{result[0]}\t{result[1]}\t{result[2]}\t{result[3]}\n")

if __name__ == "__main__":
    main()
