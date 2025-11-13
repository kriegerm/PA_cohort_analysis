# The following commands are for processing 16S rRNA gene sequencing data using QIIME 2 (version 2022.2).

conda activate /opt/miniconda3/envs/qiime2-2022.2


qiime tools import \
  --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-path ./rawreads \
  --input-format CasavaOneEightSingleLanePerSampleDirFmt \
  --output-path ./Qiime2/demux-paired-end.qza



qiime demux summarize \
  --i-data ./Qiime2/demux-paired-end.qza \
  --o-visualization ./Qiime2/demux-paired-end.qzv


qiime tools view ./Qiime2/demux-paired-end.qzv

qiime dada2 denoise-paired \
  --i-demultiplexed-seqs ./Qiime2/demux-paired-end.qza \
  --p-trunc-len-f 250 \
  --p-trunc-len-r 230 \
  --p-trim-left-f 10 \
  --p-trim-left-r 10 \
  --p-n-threads 0 \
  --output-dir ./Qiime2/Trimmed \
  --verbose


#Generate Rep Seqs
qiime metadata tabulate \
--m-input-file ./Qiime2/Trimmed/representative_sequences.qza \
--o-visualization ./Qiime2/Trimmed/representative_sequences.qzv


#Add Metadata to Sequences
qiime feature-table summarize \
  --i-table ./Qiime2/Trimmed/table.qza \
  --o-visualization ./Qiime2/Trimmed/table.qzv \
  --m-sample-metadata-file ./metadata.tsv



qiime tools view ./Qiime2/Trimmed/table.qzv

#MAKING A TREE
qiime phylogeny align-to-tree-mafft-fasttree \
  --i-sequences ./Qiime2/Trimmed/representative_sequences.qza \
  --output-dir ./Qiime2/Tree \
  --verbose


#Export Tree (needed for beta diversity analysis in Phyloseq – We used the rooted tree which is in newick format)
qiime tools export \
  --input-path ./Qiime2/Tree/rooted_tree.qza \
  --output-path ./Qiime2/Tree/Unfiltered_Rooted_tree_for_phyloseq



#BLAST to Reference

qiime feature-classifier classify-consensus-blast \
  --i-query ./Qiime2/Trimmed/representative_sequences.qza  \
  --i-reference-reads ../../HOMD_classifier_MK/homd-ref-seqs_V3V4.qza \
  --i-reference-taxonomy ../../HOMD_classifier_MK/homd-ref-taxonomy.qza \
  --p-perc-identity 0.7 \
  --p-maxaccepts 1 \
  --output-dir ./Qiime2/blast
  --verbose

#Add Taxonomy to Rep Seqs
qiime metadata tabulate \
  --m-input-file ./Qiime2/Trimmed/representative_sequences.qza \
  --m-input-file ./Qiime2/blast/classification.qza \
  --o-visualization ./Qiime2/blast/classification.qzv

mkdir Qiime2/blast/Files_For_Phyloseq

#EXPORT THE (unfiltered) TABLE
qiime tools export \
--input-path ./Qiime2/Trimmed/table.qza \
--output-path ./Qiime2/blast/Files_For_Phyloseq

#Take the ./Qiime2/blast/classification.qza and open in text editor (BBedit) and access data/.tsv Change headers to [#OTUID taxonomy confidence] and save in the Files_for_phyloseq folder. use that as the observation-#metadata-fp input and use the mapping file for the –m input for the biom add-metadata plugin! 


#Add Taxonomy to .biom file
biom add-metadata \
  -i ./Qiime2/blast/Files_For_Phyloseq/feature-table.biom \
  -m ./metadata.tsv \
 --observation-metadata-fp ./Qiime2/blast/Files_For_Phyloseq/taxonomy.tsv \
  -o ./Qiime2/blast/Files_For_Phyloseq/feature_table_w_taxonomy.biom \
  --sc-separated taxonomy \
  --observation-header OTUID,taxonomy



#For exporting the representative seqs as Fasta files for Picrust2:
qiime tools export --input-path representative_sequences.qza --output-path fasta_rep_seqs     
