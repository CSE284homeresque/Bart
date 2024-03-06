# Bart
HOMER Jr. Motif Finder
Disclaimer: We are biologists trying our best!

**Overview**
Our team is interested in identifying potential transcription factor (TF)-binding motifs that are disrupted by single nucleotide polymorphisms (SNPs) associated with disease. To this end, we developed a HOMER-esque motif finder named Bart. Bart searches genomic sequences surrounding SNP coordinates for high similarity to known TF-binding motifs, represented as position frequency matrices. Every sequence-motif alignment that includes the position of the SNP receives a 'match score.' Bart then outputs the top x% (command-line option) of positive scores, including the motif ID, sequence ID (rsID), and the starting position of the alignment.



**Future improvements**
Our initial goal was to calculate a p-value for each match and then reorder the output file by p-values. However, we encountered challenges in calculating appropriate p-values. We are actively working to overcome this challenge, but in the meanwhile, the matches are simply ranked by score.







