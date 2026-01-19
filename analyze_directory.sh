#!/bin/bash
INDIR=$1
TRIALS=$2
FOLD=$3
STARTI=$4
#NAME=$4 #For the odd cases where name isnt in the file_name.
LEN=200
FOLDPEAKS=$INDIR/peak_call
FOLDMETRICS=$INDIR/ihec_metrics

module load bedtools

echo $FOLDPEAKS
echo $FOLDPEAKS/*.bb
i_from=0
i_to=`expr $i_from + 10`
i=0
for f in $INDIR/*K.narrowPeak.gz $INDIR/*/*K.narrowPeak.gz $FOLDPEAKS/*.bed $FOLDPEAKS/*.bb;
do
## Only go through a limited set of samples (indices "from" up "to")
#if [ $i -lt $i_from ]; then
#echo skipped $i $f
#i=`expr $i + 1`
#continue
#fi
#if [ $i -ge $i_to ]; then
#break
#fi
if [ -e "$f" ]; then
echo $i $f
echo ${f##*/}
sbatch launcher_analyze_peaks.sh $f $TRIALS $FOLD ${f##*/}
fi
sleep 2
i=`expr $i + 1`
done

OUTMETRICS=analyze_peaks_results/${FOLD}/metrics/metrics_${FOLD}_file_list.txt
OUTJOINTMETRICS=analyze_peaks_results/${FOLD}/metrics/metrics_${FOLD}_joint.csv
mkdir -p "$(dirname $OUTMETRICS)"
if [ -d "$FOLDMETRICS" ]; then
echo $FOLDMETRICS
> $OUTMETRICS
for f in $FOLDMETRICS/*.txt;
do
echo $f
echo "--"
if [ -e "$f" ]; then
echo $f >> $OUTMETRICS
fi
done
python lib/merge_metrics.py $OUTMETRICS $OUTJOINTMETRICS
else
python lib/merge_jsonqc_metrics.py $INDIR $OUTJOINTMETRICS
fi

echo "Done"