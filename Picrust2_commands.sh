# The below command was used for running Picrust2. Ensure you have Picrust2 installed and your input files are correctly specified.

picrust2_pipeline.py -s ../Pediatric_cohort/Qiime2/Trimmed/fasta_rep_seqs/dna-sequences.fasta \
    -i ../Pediatric_cohort/Qiime2/Files_For_Phyloseq/feature_table_w_taxonomy.biom \
    -o picrust2_pipeline_output
    
