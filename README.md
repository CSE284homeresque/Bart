# Bart   
HOMER Jr. Motif Finder   

Group 9: Jackie Lanzalotto, Delisa Ramos, Torrey Rhyne    
Disclaimer: We are biologists trying our best!

**Overview**  
Our team is interested in identifying potential transcription factor (TF)-binding motifs that are disrupted by single nucleotide polymorphisms (SNPs) associated with disease. To this end, we developed Bart, a Python script motif finder inspired by HOMER's findMotifs.pl. Bart searches genomic sequences surrounding SNP coordinates for high similarity to known TF-binding motifs. Each sequence-motif alignment that includes the position of the SNP receives a 'match score.' The script then outputs the top percentage or number of positive scores, including the motif ID, sequence ID, alignment position, and p-value.

**Usage**  
Bart.py [-h] -i INPUT -o OUTPUT -m MOTIF [-t SCORE_THRESHOLD] [-l OUTPUT_LIMIT]

-h, --help = help message  
-i INPUT, --input INPUT = input sequences file (FASTA format)  
-m MOTIF, --motif MOTIF = motifs file (JASPAR PFM format)  
-o OUTPUT, --output OUTPUT = output results file name    
-t SCORE_THRESHOLD, --score_threshold SCORE_THRESHOLD = top scores threshold (optional, default = 0.05)  
-l OUTPUT_LIMIT, --output_limit OUTPUT_LIMIT = max number of hits outputted (optional, default = 500)

**Installation instructions**  
1.	Open terminal on JupyterHub (8CPU) and navigate to the directory where you want to clone our GitHub repository.
2.	Run the following commands to clone the repository and install required packages:
```
git clone https://github.com/CSE284homeresque/Bart.git
pip install Biopython
```
Bart also requires scipy and numpy, but these packages should already be installed. If you encounter any issues running our script, please refer to the "pip_list.txt" file, which contains a list of all the packages and version numbers that were used successfully.

**Run a test example**   
The example dataset used during the development of this tool was SNPs associated with Type 2 Diabetes (T2D), downloaded from the NHGRI-EBI GWAS Catalog.   

Input sequences = 'T2D_SNP_demo.fasta'   
We filtered out duplicate SNPs and those with missing information, resulting in a final dataset of ~3,500 unique SNPs. Using these coordinates, we extracted genomic sequences extending +/- 10 nucleotides from the SNP positions (given that TF-binding motifs are typically 6-10 nucleotides in length). Thus, each input sequence is 21 nucleotides long, with the SNP positioned at the center. The resulting sequences were saved in FASTA format, with the rsIDs serving as sequence identifiers. The sequence file is named 'T2D_SNP_full.fasta', but for testing purposes, we provided 'T2D_SNP_demo.fasta,' which contains 150 sequences.  

Motif file = 'motifs.txt'   
We downloaded the TF-binding motifs in PFM (position frequency matrix) format from the JASPAR CORE database (2020, Vertebra). 

Instructions:

1.	Open terminal on JupyterHub (8CPU) and navigate to the repository directory.
2.	Run the following command (optionally, you can change the score threshold and output limit options from their default values):
```
python Bart.py -i T2D_SNP_demo.fasta -m motifs.txt -o output_demo
```
This should take about 5 minutes to run. When finished, the 'output_demo' file will be in the repository directory.   
We also included the 'output_full' file so you can see the output from the full dataset. This took about 1.5 hours to run.    
Note: Bart processes the full dataset in about 15 minutes on my computer. These times are using the JupyterHub server (8CPU).

**Future improvements**   
1. There is variability in the p-values between runs because they are calculated using the score distribution of randomly-generated background sequences, with the number of background sequences matching the number of input sequences. As a result, the variation increases as sample size decreases. We considered generating a larger set of background sequences regardles of the number of input sequences, but we opted against it because the gain in precision did not justify the associated increase in runtime, particularly for demonstration purposes.
2. Further exploration of optimization and parallelization strategies may result in improved efficiency. Currently, Bart is about three times slower than HOMER at processing our dataset of ~3,500 sequences.
3. Improving the data output by incorporating a summary table would streamline result presentation and improve interpretation, especially when specific motifs recur frequently among the top hits. Additionally, including motif logos would be helpful for visualization purposes. These improvements would align Bart’s output more closely with that of HOMER’s.
