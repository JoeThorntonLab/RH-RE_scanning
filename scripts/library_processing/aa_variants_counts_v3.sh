#!/bin/sh
#SBATCH --job-name=AAvar_count
#SBATCH --output=AAvar_count_REP4.out
#SBATCH --error=AAvar_count_REP4.err
#SBATCH --time=02:00:00
#SBATCH --partition=caslake
#SBATCH --ntasks=15
#SBATCH --cpus-per-task=10
#SBATCH --mem-per-cpu=2000
#SBATCH --account=pi-joet1

start=`date +%s.%N`

module load parallel/latest
module load python

DIR=/project/joet1/santiagoherrera/DMS_RH_RE/221202_A00639_1438_AHLYYJDRX2-JT-JP-SHA-1s/demultiplexed_reads/merged_dmx_splits/REP4
cd $DIR

# create output directory
mkdir AA_variant_counts
OUTDIR=$PWD/AA_variant_counts

# Function to count AA variants using bash (faster than python and potentially uses less memory)
count_RH_vars () {
# '$1' is the first argument: a fasta file

# Get the file info (REBC, BBC and REP)
if echo "$1" | grep -q "GGKM\|rh"; then
        REBC=$(echo "$1" | cut -d'_' -f1-3)
        BBC=$(echo "$1" | cut -d'_' -f4)
        REP=$(echo "$1" | cut -d'_' -f5 | cut -d'.' -f1)

elif echo "$1" | grep -q "GFPneg"; then
        REBC=$(echo "$1" | cut -d'_' -f1-2)
        BBC=$(echo "$1" | cut -d'_' -f3 | cut -d'.' -f1)
        REP=$(echo "Debulk")

else
        REBC=$(echo "$1" | cut -d'_' -f1-2)
        BBC=$(echo "$1" | cut -d'_' -f3)
        REP=$(echo "$1" | cut -d'_' -f4 | cut -d'.' -f1)
fi

OUT=$(echo "$REBC"_"$BBC"_"$REP"_AA_var_count.csv) # name of output file

# Translate sequences with seqkit using frame 1 (default)
# Reformat fasta file (a single line per sequence)
# Print only the sequences (i.e. the line after the match ">")
# Extract the RH region from each sequence (positions 21-25, excluding the constant 'C' - pos23)
# Sort sequences (necessary for counting with 'uniq'. It counts adjacent lines that are equal)
# Count occurrence of unique RH sequences and replace spaces by tabs

seqkit translate "$1" | sed ':a;N;/>/!s/\n//;ta;P;D' | awk '/^>/{y=1;next}y' | cut -c 21-22,24-25 | sort | uniq -c | sed -e 's/^[[:space:]]*//' -e 's/[[:space:]]/\t/g' | tr -d ' ' > $OUT

# Append REBC, BBC and REP info to output
# Reorganize columns and transform to csv file
# Add header
awk -i inplace -v var1=$REBC -v var2=$BBC -v var3=$REP '{ print $0 "\t"var1"\t"var2"\t"var3}' $OUT
awk -i inplace -F $'\t' '{print $2","$1","$3","$4","$5}' $OUT
sed -i '1 i AA_var,Count,REBC,BinBC,REP' $OUT

# Report completion of file
echo file "$1"... completed!
}
export -f count_RH_vars


source activate /project/joet1/santiagoherrera/python_env # Python env to use seqkit
parallel --delay .2 -j $SLURM_NTASKS count_RH_vars {} ::: *.fasta
mv *.csv $OUTDIR 
conda deactivate #Leave the environment

end=`date +%s.%N`
runtime=$( echo "$end - $start" | bc -l )

echo $'·········\nAA counts DONE!\n·········'
echo "$runtime seconds total time"