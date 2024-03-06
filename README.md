# Bart
HOMER Jr. Motif Finder  
Disclaimer: We are biologists trying our best!

**Overview**  
Our team is interested in identifying potential transcription factor (TF)-binding motifs that are disrupted by single nucleotide polymorphisms (SNPs) associated with disease. To this end, we developed a HOMER-esque motif finder named Bart. Bart searches genomic sequences surrounding SNP coordinates for high similarity to known TF-binding motifs, represented as position frequency matrices. Each sequence-motif alignment that includes the position of the SNP receives a 'match score.' Bart then outputs the top x% (command-line option) of positive scores, including the motif ID, sequence ID (rsID), and the starting position of the alignment.

**Usage**  
Bart.py [-h] -i INPUT -o OUTPUT -m MOTIF [-t SCORE_THRESHOLD] [-l OUTPUT_LIMIT]

-h, --help = Help message  
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
1.	Open terminal on JupyterHub (8CPU) and navigate to the repository directory.
2.	Run the following command:

    python Bart.py -i T2D_SNP_demo.fasta -m motifs.txt -o output_demo

Optionally, you can change the score threshold (-t) and output limit (-l) options from their default values.

**Future improvements**  
Our initial goal was to calculate a p-value for each match and then reorder the output file by p-values. However, we encountered challenges in calculating appropriate p-values. We are actively working to overcome this challenge, but in the meanwhile, the matches are simply ranked by score.


