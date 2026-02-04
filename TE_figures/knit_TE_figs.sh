#!/bin/bash
#SBATCH --time=16:00:00
#SBATCH --mem=240G
#SBATCH --account=rrg-user-ad
#
#for epiatlas data was 16h, 240G (conservative, more around 5h)
#for test (default) was 5h, 60G (conservative, handful of minutes)
module load r/4.4.0
#R -e "rmarkdown::render('TE_figures.Rmd', params = list(input = '../demo_data/epiatlas_core_data.Rdata'))"
R -e "rmarkdown::render('TE_figures.Rmd')"
