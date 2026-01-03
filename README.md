match_islands.py is an in-house Python pipeline for comparative analysis and network construction of genomic islands (e.g. virulence, resistance, symbiotic or metabolic islands) predicted by GIPSy2. The pipeline performs protein-level, sequence-based comparisons between islands across multiple bacterial strains and generates outputs ready for downstream interpretation and visualization in Cytoscape.

The main goal of match-islands is to:
- Compare genomic islands across two or more strains using protein sequence similarity;
- Identify shared, partially related, and exclusive islands;
- Generate interpretable summaries and similarity networks for exploratory and comparative genomics analyses.

# How the pipeline works?
1. Input detection: The script must be placed inside a directory containing the folders of interest (in the same directory of LineageName_Islands_fisher/ folders). The script automatically searches for GIPSy2 amino acid outputs following the structure:
LineageName_Islands_fisher/
   └── Amino_acids/
       └── Virulence_Island_*.faa
All .faa files found are stored internally and compared in a pairwise manner.

2. Protein-level sequence alignment: For each unique pair of island files, all protein sequences will be read using Bio.SeqIO, by Biopython and then, each protein from "Island A" is aligned against all proteins from "Island B", using  Bio.Align.PairwiseAligner. In the global mode, match_score = 1; mismatch_score = 0. So, in this way, sequence identity is calculated as "identity = number of matches / length of the longest sequence". This alighnment is perform at protein level and this approach ensures true sequence-based comparison, independent of annotation inconsistencies.

3. Best-hit selection and identity thresholds: For each protein, only the best match in the target island is retained. Each one of those matches are classified using a default identity threshold such as "High similarity" as ≥ 90% and "Low similarity" as < 90%. Thresholds can be easily adjusted within the script.

4. Parallel processing: Pairwise comparisons are computationally expensive. In this sense, to improve performance, the pipeline uses "multiprocessing.Pool", that automatically detects the number of available CPU cores to use. This makes the pipeline suitable for datasets with dozens of strains and islands. Also, alignment is computationally expensive so I recommend run using screen or tmux on Linux.
   
5. Similarity network construction: From the pairwise comparison summaries, the pipeline will extract average similarity values and classify island relationships into three categories: High (≥ 90%) or Medium (50–89%) or Low (< 50%). This analyse will be export to a CSV file compatible with Cytoscape with node names are automatically simplified to improve network readability.

6. Global summary: A final report is generated containing number of strains analyzed, total number of islands, number of islands per strain and distribution of similarity levels in the network.

# Output files
After execution, the following files are produced:
- "SUMMARY_virulence_island_comparison_results.txt": Pairwise summaries including average similarity and unmatched genes;
- "cytoscape_virulence_network_filtered.csv": Similarity network file for Cytoscape;
- "cytoscape_virulence_network_filtered_renamed.csv": Network file with simplified node names.
- "SUMMARY_virulence_islands_global_analysis.txt": Global overview of the analysis.

# Dependences 
Python libraries:
-> Bio.SeqIO: Reading FASTA and .faa files. 
-> Bio.Align.PairwiseAligner: Pairwise global alignment.
-> os, pathlib: Directory and file handling.
-> collections: Auxiliary data structures.

