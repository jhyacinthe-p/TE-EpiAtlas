

# TE analysis

Code and processed data for analyses and figures of the Hyacinthe et al (2024) preprint.

Code for measurement and analysis of Transposable Element (TE) overlap in ChIP-seq samples with random simulations as controls (`analyze_peaks_p2.py`). In addition, the downstream analysis adapted to apply the TE analysis to the EpiAtlas dataset is available.  


Usages

- Prerequisites & installation
- Usage Notes
- Standalone
- With additional analysis
- Applied to directory
- Plots and figures
    - Merge Tables
    - R code
- Overall Pipeline Process
- tested versions


# Prerequisites & installation

Bedtools >= 2.29  
python pandas  
python pybedtools

>pip install pandas  
>pip install pybedtools

or

>pip install -r requirements.txt

be sure that bedtools is loaded

>module load bedtools  

prerequisite installation should be brief and take less than 10 minutes.

------

If using github, library and other required files are stored seperately due their size:
They are available in the zenodo repository as "libfiles.zip". https://doi.org/10.5281/zenodo.18343281 

Simply unzip lib and files folder directly within the TE-EpiAtlas folder (this repository). 

(te_epiatlas_code.zip includes everything: the code, the lib and files and demo data but is much larger)


------

Full EpiAtlas data available at https://ihec-epigenomes.org/epiatlas/data/

Small sample subset for demo purposes distributed in `demo_data`.

# Usage Note

The standalone analysis `analyze_peaks_p2.py` can be used  on any bed file-like input and run on a normal computer. 10 iterations take about 5 minutes and runs linearly. Thus, it should take about 500 min (or 8h) for the 1000 iterations. If multiple process are available, it can run in parallel which should be much faster.  

The applied to directory `analyze_directory.sh` script was made specifically to run on the EpiAtlas full dataset. The full dataset being multiple terabytes and total running time of thousands of samples for a thousand of iteration being days of runtime, it is not expected to run on normal computers. It was run on Digital Research Alliance of Canada server's Slurm system.  

The content of TE_figures can be run to output the main figures. These scripts are also closely intertwined and dependent on the EpiAtlas dataset. It was run on Digital Research Alliance of Canada server's Slurm system. It takes about 5h on the full epiatlas dataset (`epiatlas_core_data.Rdata`) and about 30 min on the smaller test data (`core_data_test.Rdata`).

# Standalone

Used as a single python script to obtain the amount of peaks (from input file) that overlapped TEs and the associated count from simulated samples (iteration count) as control.

Intersects a given peak file with TE and returns the overlap count in comparison to associated random simulated samples.  
returns a csv table, a csv of peak distribution relative to tss and a bed file of intersection with repeatmaker

>Usage:  
>python analyze_peaks_p2.py input_file -g genome -l length -n size -r iteration_count -o output_file -p processes

input_file: [str] standard bed file  
genome: default=38 [int] 38  
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

# With additional Analysis

The launcher bash script, will do bigbed conversion if necessary, launch the standalone analysis and report GC content and peak counts within 10kbp windows of the whole genome

>Usage:  
>sh launcher_analyze_peaks.sh [input] [iteration count] [output directory] [optional name] [optional length]

# Applied to directory

The analyze directory bash script will use an input folder and use all narrowPeak.gz, Bigbed or bed files within and perform the analysis (with additional analysis) of all those files. This is what was used for the EpiAtlas analysis.

>Usage:  
>sh analyze_directory.sh [input directory] [trial count] [Output folder]

>[!NOTE]  
>analyze_directory uses Slurm (sbatch) to run launcher_analyze_peaks.sh on each sample in parrallel. To avoid inadvertedly flooding jobs, a maximum of 10 iteration is implemented. (can be changed or removed in code)\
>Do not forget to add your own Slurm attributes (or change the code to sh)

------

# Plot and Figures

The aformentionned analysis outputs a csv table containing the amount of peaks overlapping repeats from repeatmasker (and expected counts). It is the main input of all the downstream analysis and figures generated from our R markdown code.

## Merge Tables

The analysis outputs many csv tables and other files (such as genome windows peak count). We merge the tables and organize all the necessary data into one R object used as input for the plots.

>Usage:  
>Rscript TE_table_merger_script.R [input directory] [output file(.Rdata)] [start_path] [trial count] [pval tresh] 

## Plot Figures

Obtain a HTML page with the figures and data from publication (pre any post-processing)

>Usage:  
>R -e "rmarkdown::render('TE_figures.Rmd', params = list(input = 'core_data.Rdata', start_path = '../files'))"

input: core_data.Rdata file (Output of the table merge)  
start_path: working directory where the subfolders (lib, data, etc) are.

------


# Overall Pipeline Process

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
  
# Tested versions

program was tested with:  
python 3.8.10 and 3.11.4  
bedtools v2.27.1 and v2.31.0  
pybedtools v.0.11.0 and v.0.12.0  
pandas v2.0.3  

see requirements.txt and renv.lock for extended details.  

This tool features and includes data from UCSC genome/table browser tracks, genome data and tools.  
https://genome.ucsc.edu/

