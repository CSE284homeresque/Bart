# Bart
HOMER Jr. Motif Finder  
Disclaimer: We are biologists trying our best!

**Overview**  
Our team is interested in identifying potential transcription factor (TF)-binding motifs that are disrupted by single nucleotide polymorphisms (SNPs) associated with disease. To this end, we developed a HOMER-esque motif finder named Bart. Bart searches genomic sequences surrounding SNP coordinates for high similarity to known TF-binding motifs, represented as position frequency matrices. Each sequence-motif alignment that includes the position of the SNP receives a 'match score.' Bart then outputs the top x% (command-line option) of positive scores, including the motif ID, sequence ID (rsID), and the starting position of the alignment.

**Usage**  
Bart.py [-h] -i INPUT -o OUTPUT -m MOTIF [-t SCORE_THRESHOLD] [-l OUTPUT_LIMIT]

-h, --help = help message  
-i INPUT, --input INPUT = input sequences file (FASTA format)  
-m MOTIF, --motif MOTIF = motifs file (JASPAR format)  
-o OUTPUT, --output OUTPUT = output results file  
-t SCORE_THRESHOLD, --score_threshold SCORE_THRESHOLD = top scores threshold (optional, default = 0.05)  
-l OUTPUT_LIMIT, --output_limit OUTPUT_LIMIT = max number of matches outputted (optional, default = 500)

**Installation instructions**  
1.	Open terminal on JupyterHub (8CPU) and navigate to the directory where you want to clone our GitHub repository.
2.	Run the following commands to clone the repository and install required packages:

  	git clone https://github.com/CSE284homeresque/Bart.git

  	pip install biopython

Note: If you encounter any issues running our script, please refer to the "pip_list.txt" file. This file contains a list of all the packages / version numbers that were used successfully.

**Run a test example**   
The example dataset used during development of this tool was SNPs associated with Type 2 Diabetes (T2D), downloaded from the NHGRI-EBI GWAS Catalog.   

Input sequences = T2D_SNP_demo.fasta   
We filtered out any duplicate SNPs and those with missing information, resulting in a final dataset of approximately 3,500 unique SNPs. Using these coordinates, we extracted genomic sequences extending +/- 10 nucleotides from the SNP positions (given that TF-binding motifs are typically 6-10 nucleotides). Thus, each input sequence is 21 nucleotides long, with the SNP located at the center position. The complete sequence file is named 'T2D_SNP_full.fasta', but for testing purposes, we also provided 'T2D_SNP_demo.fasta,' which contains 150 sequences.   

Motif file = motifs.txt   
The motif file was downloaded from the JASPAR CORE database (2020, Vertebra).

Instructions:

1.	Open terminal on JupyterHub (8CPU) and navigate to the repository directory.
2.	Run the following command (optionally, you can change the score threshold and output limit options from their default values):

    python Bart.py -i T2D_SNP_demo.fasta -m motifs.txt -o output_demo

This should take approximately 1 minute to run. When finished, the output_demo file will be in the repository directory.

**Future improvements**  
Our initial goal was to calculate a p-value for each match and then reorder the output file by p-values. However, we encountered challenges in calculating appropriate p-values and implementing the code. We are actively working to overcome this challenge, but in the meanwhile, the matches are simply ranked by score.


