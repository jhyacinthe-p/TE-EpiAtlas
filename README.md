

# TE analysis

Usages

- Prerequisites
- Standalone
- With additional analysis
- Applied to directory
- Plots and figures
    - Merge Tables
    - R code
- Overall Pipeline Process


## Prerequisites

Bedtools >= 2.29  
python pandas  
python pybedtools

>pip install pandas  
>pip install pybedtools

or

>pip install -r requirements.txt

be sure that bedtools is loaded

>module load bedtools


------

Full data available at https://ihec-epigenomes.org/epiatlas/data/

Small sample subset for demo purposes distributed in `demo_data`.


## Standalone

Used as a single python script to obtain the amount of peaks (from input file) that overlapped TEs and the associated count from simulated samples (iteration count) as control.

Intersects a given peak file with TE and returns the overlap count in comparison to associated random simulated samples.  
returns a csv table, a csv of peak distribution relative to tss and a bed file of intersection with repeatmaker

>Usage:  
>analyze_peaks.py input_file -g genome -l length -n size -r iteration_count -o output_file -p processes

input_file: [str] standard bed file  
genome: default=38 [int] 19 or 38  
lenght: default=200 [int] lenght of region directly around peak to consider  
size: default=-1 [int] number of simulated peaks per trial. -1 sets to input peak count  
iteration_count: default=10 [int] random iteration count, recommended 1000 (Slow).  
output_file: default=analyze_peaks_output.csv [str] csv table output file  
processes: parralelize the work using multiple processes  

>[!IMPORTANT]  
>If an error occurs related to parallel_apply:
>  
>>parallel_apply  
>>        raise StopIteration  
>>StopIteration  
>  
>>it may be fixed by replacing **"raise StopIteration"** by **"return"** on line 2962 and 2982. (as of *v.0.11.0*) with parallel_apply script in myenv/lib/python3.11/site-packages/pybedtools/bedtool.py or wherever your pybedtools installation is.

## With additional Analysis

The launcher bash script, will do bigbed conversion if necessary, launch the standalone analysis and report GC content and peak counts within 10kbp windows of the whole genome

>Usage:  
>launcher_analyze_peaks.sh [input] [iteration count] [output directory] [optional name] [optional length]

## Applied to directory

The analyze directory bash script will use an input folder and use all narrowPeak.gz, Bigbed or bed files within and perform the analysis (with additional analysis) of all those files. This is what was used for the EpiAtlas analysis.

>Usage:  
>analyze_directory.sh [input directory] [trial count] [Output folder]

>[!NOTE]  
>analyze_directory uses Slurm (sbatch) to run launcher_analyze_peaks.sh on each sample in parrallel. To avoid inadvertedly flooding jobs, a maximum of 10 iteration is implemented. (can be changed or removed in code)\
>Do not forget to add your own Slurm attributes (or change the code to sh)

------

## Plot and Figures

The aformentionned analysis outputs a csv table containing the amount of peaks overlapping repeats from repeatmasker (and expected counts). It is the main input of all the downstream analysis and figures generated from our R markdown code.

### Merge Tables

The analysis outputs many csv tables and other files (such as genome windows peak count). We merge the tables and organize all the necessary data into one R object used as input for the plots.

>Usage:  
>Rscript TE_table_merger_script.R [input directory] [output file(.Rdata)] [start_path] [trial count] [pval tresh] 

### Plot Figures

Obtain a HTML page with the figures and data from publication (pre any post-processing)

>Usage:  
>R -e "rmarkdown::render('TE_figures.Rmd', params = list(input = 'core_data.Rdata', start_path = '../files'))"

input: core_data.Rdata file (Output of the table merge)  
start_path: working directory where the subfolders (lib, data, etc) are.

------


## Overall Pipeline Process

The many scripts should be run according to the folowing sequence:

- analyze_directory.sh
  - use as input a directory with all the bed/narrowpeak.gz files you want analyzed
  - In the output directory will be contain all files resulting from the analysis
- Rscript TE_table_merger_script.R
  - use as input the output directory of `analyze_directory.sh`
  - the output will be a `core_data.Rdata` file containing compiled tables from the files
- TE_figures.Rmd
  - use as input the `core_data.Rdata` from `TE_table_merger_script.R`
  - the output will be a html page `TE_figures.html` with plots, figures and data from varous analysis.
