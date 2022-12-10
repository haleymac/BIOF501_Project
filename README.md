BIOF501 Term Project: Snakemake Alignment QC of Chromosome 22 Workflow
===

By: Haley MacDonald 
---

___

Repository Contents
---

**Note that for this project the students in this class have been provided with a Jupyter Notebook Environment to host this workflow in. The Jupyter notebook can be found at: http://class.cidgoh.ca/ .**

**Due to the fact that this workflow can be run on any human .bam file, as well file size constraints on github and because the efficacy of this workflow will be assesed in the student Jupyter notebook, I have decided not to include the single cell whole genome sequencing .bam input files and the reference genome file in this github repository. THEY CAN BE FOUND ALREADY DOWNLOADED IN THE STUDENT JUPYTER NOTEBOOK. If you wish to use these files to run the workflow on your own computer, please reach out to me and I will provide them.**

### Directories
* input: In the student jupyter notebook this will contain .bam files that can be used as input in the pipeline. For more information on the .bam file data included in the notebook see the *data* section. For more information on this folder see the .txt file within it.
* output: is an empty directory that output files from the workflow can be directed into. For more information on this folder see the .txt file within it.
* reference_genome: In the student jupyter notebook this will contain the fasta file for the GRCH37-lite human reference genome. For more information on this genome and where to download it see the *data* section. For more information on this folder see the .txt file within it. 

### Files
* biof501env.yaml: This is a text file that contains the list of dependencies and channels to prioritize downloading them from that are relevant for this workflow. This file can be used to create an environment to run the workflow in. 
* config.yaml: This is a text file that contains user-specified parameters relevant to the workflow. Currently this file contains information necessary to run the workflow on the included data in the provided student environment. 
* snakefile: Snakemake source code used to run the workflow
* dag.pdf: Representation of the overall workflow in PDF format. Viewable in web browser.
* A90553C_metadata.csv: File containing information about the .bam files included in the input file of this repository. This file contains columns that give information about a single cell's condition prior to sequencing - the cell-call column indicates whether the cell is alive or not (C1 indicates a live cell), and the experimental condtion indicates which cell-cycle phase the cell was in prior to sequencing. 

___

Introduction
---

### Background

Alignment of whole genome sequencing (WGS) data to a reference genome is a critical step in many next generation sequencing experiments, and the quality of the resulting data can often impact conclusions drawn from a study. As such, performing quality control on alignment data is an important process to ensure that erroneous conclusions are not being drawn from bad data.1 

Alignment quality can be assessed through several metrics, with the most used being the average read depth, as well as the percentage of the genome covered at that depth (breadth). Though these two metrics are helpful, other metrics should also be examined as just examining average depth and breadth of sequencing can be misleading. Several factors can skew read depth, such as GC bias, duplicate reads and insert size distribution so these factors should be examined before using any alignment data in further downstream analysis.2,3 

In addition to ensuring good quality data is used in downstream analysis, alignment quality metrics can also be used to compare and optimize experimental conditions for later experimentation. If one were to compare the quality scores of WGS data one source but with differing experimental conditions, it would likely be possible to ascertain which of the methods used produced better quality alignment data. 

This workflow employs the use of several GATK tools to create files containing alignment quality metrics from aligned human WGS data, and attempts to improve the ability to compare quality metrics between two datasets by subsetting the data and reducing the intensity of the computations required to generate these files. 

### Purpose

This pipeline aims to compute and output several quality control metrics that indicate alignment quality. These metrics can be used to ensure alignment data is of high quality before using the data for further analysis. 

These metrics can also be used to identify the effect of different experimental conditions on alignment quality. For example, the bam files in this repository contain single cell WGS data from cells that were FACS-sorted by cell-cycle state prior to sequencing. This pipeline could be run on this data to determine if there are differences in quality of alignment between cells in different cell cycle states. Information on the particular states of the cells included in this repository is available in the accompanying metadata file. 

### Dataset 

Though this pipeline has been designed to run on any human whole genome sequencing aligned .bam file, a dataset has been included in this repository that the pipeline can be run on. The config file has already been set up to run this particular dataset in the provided student environment (see *config* section for more information on config file setup for alternative datasets or environments).

The sequence data included in this repository was generated by Laks et al. in 2019 for their paper ‘Clonal Decomposition and DNA Replication States Defined by Scaled Single-Cell Genome Sequencing.’ (DOI: 

[Laks et al. 2019](https://doi.org/10.1016/j.cell.2019.10.026). 

The cell line used for this project was the GM18507 lymphoblastoid cell line originally from the HapMap Project. For more information see: 

[GM18507 lymphoblastoid cell line](https://catalog.coriell.org/0/Sections/Search/Sample_Detail.aspx?Ref=GM18507).


Four .bam files containing single cell whole genome sequencing data that has been aligned to the GRCH37-lite version of the human reference genome have been included in the /input folder in this repository. These cells were FACS-sorted by cell-cycle state prior to sequencing, and imaged to ensure they were alive. The individual information for each cell can be accessed in the A90553C_metadata.csv file included in this repository. The four cells included in this repository are all alive. Two of the cells were in G1 prior to sequencing, and two were in S phase. The workflow can be run on these cells to generate alignment quality metrics for each cell, which can then be compared between the two conditions (see *results* section). 


The reference genome these cell genomes were originally aligned to, and therefore the reference included in this repository to be run with cells in the workflow is the GRCH37-lite version of the human reference genome. It can be downloaded here: 

[GRCH37-lite](https://www.bcgsc.ca/downloads/genomes/9606/hg19/1000genomes/bwa_ind/genome)


___
 
Workflow Overview
---
 
This workflow was built using Snakemake, a useful tool to create reproducible and scalable data analyses. The workflow has been designed to run on any previously aligned human whole genome sequencing .bam files (see *config* section for more information on how to specify specific files) and will output several useful quality control metrics for each file it runs on. 
 
##### The key steps in the workflow are:
1. Use samtools to index .bam files 
2. Create a subset of the .bam files that only includes chromosome 22 and save as a .sam file 
    *Note: Chromosome 22 was chosen arbitrarily as a means to reduce the computational intensity of the workflow
3. Compress the chromosome 22 .sam files back into .bam files for further processing 
4. Use GATK tools to produce several quality control metrics that can be used to assess the quality of the alignment. GATK tools include:
    ..* CollectWgsMetrics: Collects metrics about the fractions of reads that pass base- and mapping-quality filters as well as coverage (read-depth) levels for WGS analyses
    ..* CollectGCBiasMetrics: Collects information about the relative proportions of guanine (G) and cytosine (C) nucleotides in a sample - this can be used to determine if there is GC bias in an alignment 
    ..* CollectDuplicateMetrics: Identifies duplicate reads (reads that are mapped to the exact same location) - this can be used to ensure that duplicate reads are not falsely inflating the average sequencing depth 
    ..* CollectInsertSizeMetrics: Identifies insert size that was used in sequencing - this can show the insert size distribution and can give insight into how library preperation methods can affect sequencing and alignment
    
A figure illustrating these key steps:
*Note: the indexing step appears to be beside the gatk steps because the direct output of the indexing step (a .bai file) is not used as input for the next step, and as such snakemake doesn't register that it is occuring first in the workflow when it creates this workflow image. In actuality this step takes place first as it is a necessary precurser to subsetting the .bam file by chromosome. 
 
![alt text](link_to_image_on_gihub "Pipeline Overview")
 
 
 ___
 
 Installation and Running of the Workflow
 ---
 
 *Note: The course instructors have provided a Student Jupyter Notebook to host the workflow in. The installation/usage commands differ for this environment vrs. the  
 
 Installing this pipeline requires conda and git. Instructions for installing these software can be found at the following links:
 
 [I'm an inline-style link](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html)
 
 [I'm an inline-style link](https://git-scm.com/book/en/v2/Getting-Started-Installing-Git)
 
 
### Dependencies 
 
 The dependencies and channels they can be downloaded from that this workflow relies upon are listed in the biof501env.yaml file included in this repository. 
 
### Installation

##### On Your Own Computer:
To access this workflow on your own computer, open a terminal and navigate to where you would like to save the repository. Then clone the repository by typing: 

```
git clone LINK TO GIT REPOSITORY 
```

Following this, navigate to the new directory 

```
cd /path-to-the-workflow-directory
```

Once in the new directory, create a new environment from the biof501env.yaml file included in this repository by typing:

```
conda env create -f biof501env.yaml
```

Once the environment has been resolved, activate it with: 

```
conda activate biof501env
```

You're now ready to move on to the *config* step!

##### In the BIOF501 Student Jupyter Notebook

The environment containing all of the necessary dependencies for this workflow has already been created in the instructor-provided student jupyter notebook. 

To access the workflow in this environment, navigate to Haley MacDonald's notebook by typing the follow command into a terminal:
```
cd /home/jupyter-macdonaldhaley7
```

Then activate the environment by typing:
```
conda activate biof501env
```
You're now ready to move on to the *config* step!

### Config 

This workflow has been designed to run for any human whole genome sequencing .bam files. The parameters the work flow will run on are set by the user in a config.yaml file. This config.yaml file has several mandatory field that must be filled prior to running the workflow, including:

..* input_dir: user must specify path to an input directory containing .bam files they would like to run through the workflow 
..* output_dir: user must specify path to an output directory they would like their output files to be generated in
..* reference_genome: user must specify path to the downloaded reference genome .fa file that their .bam files were originally aligned to - **this MUST be the same reference genome version that was used for alignment for the workflow to work properly**
..* sample_ids: user must specify the file names (without .bam extention) that they would like to run through the workflow (ex. to run SA928-A90553C-R44-C19.bam through the workflow, SA928-A90553C-R44-C19 must be included in the config file in the sample_id section)

The config.yaml file included in this repository has already been configured to run the workflow on all .bam files included in this repository. If the user would like to run the workflow on these samples, they should skip this step. 

You're now ready to move on to the *usage* step!

### Usage 

Once you have activated the biof501env environment in the directory containing the necessary files, and your config file is set up to run the files you would like, you can run the workflow by typing:
```
snakemake --cores 1
```
*Note: If you would like to parallelize the workflow over several samples, specify this by changing the number of cores indicated in the above command.*

___

Outputs and Results
---

### Outputs

Each GATK step of this workflow produces it's own output of file(s) containing quality metrics. All of these files should appear in the output folder specified in the config.yaml file after the workflow has been run. 

Output Files for each step:

1. CollectWgsMetrics
    ..* wgs_metrics.txt: A text file containing whole genome sequencing metrics such as read depth and coverage breadth. These files will be in the format "output_directory_specified_in_config/sample_id_specified_in_config.wgs_metrics.txt"

2. CollectGCBiasMetrics
    ..* gc_bias_metrics.txt: A text file containing gc bias metrics. These files will be in the format "output_directory_specified_in_config/sample_id_specified_in_config.gc_bias_metrics.txt"
    ..* gc_bias_metrics_chart.pdf: A pdf file containing a GB bias plot of all reads. These files will be in the format "output_directory_specified_in_config/sample_id_specified_in_config.gc_bias_metrics_chart.pdf"
    ..* gc_bias_summary.txt:  A text file containing a summary of the GC bias metrics for all reads. These files will be in the format "output_directory_specified_in_config/sample_id_specified_in_config.gc_bias_summary.txt"
  
3. CollectDuplicateMetrics
    ..* duplicate_metrics.txt: A text file containing duplicate reads metrics. These files will be in the format "output_directory_specified_in_config/sample_id_specified_in_config.duplicate_metrics.txt"

4. CollectInsertSizeMetrics
    ..* insert_metrics.txt: A text file containing insert size metrics. These files will be in the format "output_directory_specified_in_config/sample_id_specified_in_config.insert_metrics.txt"
    ..* insert_metrics_histo.pdf: A pdf containing a histogram of the distribution of insert size sequenced. These files will be in the format "output_directory_specified_in_config/sample_id_specified_in_config.insert_metrics_histo.pdf"


### Results
*These Results were generated from running the workflow on the .bam files included in this workflow*

To illustrate the results without overwhelming a viewer, I have included the result files from only one of the included samples run through the workflow - those from cell SA928-A90553C-R44-C19 which was in G1 at time of sequencing. 

Output by GATK step for cell SA928-A90553C-R44-C19:
1. CollectWgsMetrics

**wgs_metrics.txt**

The following image is a screenshot of the 'METRICS CLASS' section of the wgs_metrics.txt file generated from running the workflow on this cell. 
![alt text](link "WGS METRICS CLASS")

This file contains several WGS quality metrics, including MEAN COVERAGE - otherwise called average depth. These metrics are useful for understanding the overall quality of alignment.  

2. CollectGCBiasMetrics

**gc_bias_summary.txt**

The following image is a screenshot of the 'METRICS CLASS' section of the gc_bias_summary.txt file generated from running the workflow on this cell: 
![alt text](link "GC BIAS METRICS CLASS")

This file contains several GC bias metrics such as GC_dropout that are useful for assessing the effect of GC bias on alignment. 


**gc_bias_metrics_chart.pdf**

The following image is the plot generated for the GC bias metrics step in the workflow:
![alt text](link "GC BIAS CHART")


3. CollectDuplicateMetrics

**duplicate_metrics.txt**

The following image is a screenshot of the 'METRICS CLASS' section of the duplicate_metrics.txt file generated from running the workflow on this cell: 
![alt text](link "DUPLICATE METRICS CLASS")

This file contains several duplicate reads metrics such as READ_PAIR_DUPLICATES that are useful for assessing the effect of duplicate reads on the average read depth. 

4. CollectInsertSizeMetrics

**insert_metrics.txt**

The following image is a screenshot of the 'METRICS CLASS' section of the insert_metrics.txt file generated from running the workflow on this cell: 
![alt text](link "INSERT METRICS CLASS")

This file contains several insert distribution metrics such as the MEDIAN_INSERT_SIZE which are useful for assessing the efficacy of library preperation methods on data quality. 

**insert_metrics_histo.pdf***

The following image is the histogram generated for the insert size metrics step of the workflow: 
![alt text](link "INSERT SIZE HISTOGRAM")

This plot does an excellent job of illustrating the insert size distribution that was sequenced. 

___

References
---

1. Hadfield, J. & Eldridge, M.D. Multi-genome alignment for quality control and contamination screening of next-generation sequencing data. Front. Genet. 20 (2014).
2. Chen, Y.C. et al.  Effects of GC Bias in Next-Generation-Sequencing Data on De Novo Genome Assembly. PLOS ONE 8, 4 (2013).
3. Krasnenko, A. et al. Effect of DNA insert length on whole-exome sequencing enrichment efficiency: an observational study. Advances in Genomics and Genetics 8, 13-15 (2018).
4. Laks, E. Clonal Decomposition and DNA Replication States Defined by Scaled Single-Cell Genome Sequencing Cell 179, 1207-1221 (2019).

