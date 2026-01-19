#!/bin/bash
#SBATCH --time=3:00:00
#SBATCH --cpus-per-task=10
#SBATCH --mem=20G
#SBATCH --account=rrg-user-ad
__usage="Usage: $0
launcher_analyze_peaks.sh <input> <iteration count> <output directory> <optional name> <optional length>
analyze_peaks launcher, will convert bigbed if necessary and then launch analyse_peaks.
	analyse_peaks: Intersects a given peak file with TE and returns the overlap count in comparison to associated random simulated samples.
returns a csv table, a csv of peak distribution relative to tss and a bed file of intersection with repeatmaker

In addition, the launcher will also report GC contents and peak counts within 10kbp windows of the whole genome.

input: [str] standard bed file
iteration_count: [int] random iteration count, recommended 1000 (Slow).
output directory: [str] csv table output file
name: default=input file name [str] output files base name
lenght: default=200 [int] lenght of region directly around peak to consider
---
Fixed arguments (can only be changed using analyse_peaks.py directly)
genome: 38 [int] 19 or 38
size: -1 [int] number of simulated peaks per trial. -1 sets to input peak count
"

if [[ $# < 3 || $# > 5 ]]; then
	echo "$__usage"
	exit 1
fi

module load bedtools
source myenv/bin/activate


INBED=$1
TRIALS=$2
FOLD=$3
NAME=$4 #For the odd cases where name isnt in the file_name.
LEN=${5:-200}
NAMEBED=${NAME:-${INBED##*/}}
timestamp_start=`date "+%Y-%m-%d_%T"`
echo $INBED
#In case of bigbed, convert
if [[ "$INBED" == *.bb ]]; then
	echo "convert_to_bed $INBED"
    CONVBED=analyze_peaks_results/${FOLD}/Beds/${NAMEBED%.bb}.bed
    mkdir -p "$(dirname $CONVBED)"
    lib/bigBedToBed $INBED $CONVBED
    INBED=$CONVBED
    NAMEBED=${NAMEBED%.bb}.bed
fi
#Set up directories
OUTBED=analyze_peaks_results/${FOLD}/${NAMEBED%.bed}_bed_results_${LEN}bp_${TRIALS}t.csv
OUTGC=analyze_peaks_results/${FOLD}/GC/${NAMEBED%.bed}_GC.bed
OUTWIND=analyze_peaks_results/${FOLD}/windows/${NAMEBED%.bed}_10kwin_overlap.bed
mkdir -p "$(dirname $OUTBED)"
mkdir -p "$(dirname $OUTGC)"
mkdir -p "$(dirname $OUTWIND)"
#Launch analysis
python3 analyze_peaks_p2.py $INBED -r $TRIALS -o $OUTBED -p 10
bedtools nuc -fi lib/hg38.fa -bed $INBED > $OUTGC
bedtools intersect -a lib/hg38_10k_window.bed -b $INBED -c > $OUTWIND
timestamp_end=`date "+%Y-%m-%d_%T"`
echo "$INBED,$timestamp_start,$timestamp_end" >> analyze_peaks_results/${FOLD}/${FOLD}_completed_files.txt
sleep 5
