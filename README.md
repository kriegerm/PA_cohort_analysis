These are all the scripts to produce the results in the manuscript "Paired oral clinical specimens reveal the underlying ecology supporting the emergence of inflammophilic microbiome communities"

Relevant sequencing files can be found in SRA PRJNA1312071.


*A brief walkthrough of the scripts in this repository is below:*

**MK_Microbiome_Functions.R**: This R script contains custom functions used throughout the analysis.

**Qiime2_commands.sh**: This file contains all the commands used to process the sequencing data using Qiime2/DADA2. They can be run in terminal from a conda environment with Qiime2 version 2022.2 installed.

**01_Process_Qiime2.Rmd**: This R markdown file contains the code to import the Qiime2 data into R and process it using the phyloseq package. It includes steps for filtering and transforming the data.

**02_cohort_analysis.Rmd**: This R markdown file contains the code to perform most of the cohort-level analyses, including alpha and beta diversity analyses, taxonomic summaries, and visualizations.

**anvio_commands.txt**: This file contains the commands used to run anvi'o for inferred metagenomic profiling.

*Picrust2_commands.sh**: This file contains the commands used to run PICRUSt2 for functional prediction of the microbiome data.

**03_functional_analysis.Rmd**: This R markdown file contains the code to perform functional profiling of the microbiome data using PICRUSt2 and anvi'o data.




Full knit versions of these files that were used for publication are also included in this repository as HTML files:
**01_Process_Qiime2.html**
**02_cohort_analysis.html**
**03_functional_analysis.html**