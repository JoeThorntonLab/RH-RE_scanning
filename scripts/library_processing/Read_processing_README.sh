## Read processing of NovaSeq1 data: 221202_A00639_1438_AHLYYJDRX2-JT-JP-SHA-1s
## Last modified: SHA, 12/08/22

# Load modules in cluster
module load python
module load sickle/1.33

## Global variables:
MAIN_FOLDER=/project2/joet1/santiagoherrera/DMS_RH_RE/221202_A00639_1438_AHLYYJDRX2-JT-JP-SHA-1s/FastQ #modify this
R1_READS=JT-JP-SHA-1s_S1_R1_001.fastq.gz #modify this
R2_READS=JT-JP-SHA-1s_S1_R2_001.fastq.gz #modify this

cd $MAIN_FOLDER

## Step 1: Check quality statistics of the fastq files
# Create a table-histogram of the read-length distribution (awk is specially useful for big datasets: https://www.biostars.org/p/72433/)
zcat $R1_READS | awk 'NR%4 == 2 {lengths[length($0)]++} END {for (l in lengths) {print l, lengths[l]}}' > read_lengths_R1.txt
zcat $R2_READS | awk 'NR%4 == 2 {lengths[length($0)]++} END {for (l in lengths) {print l, lengths[l]}}' > read_lengths_R2.txt

## Step 2: Remove reads with average Phred score < Q30 and keep minimum length of 79
sickle pe -f $R1_READS -r $R2_READS -t sanger \
-o $(echo $R1_READS | rev | cut -d'/' -f1 | rev | cut -d'.' -f1)_trimmed.fastq.gz -p $(echo $R2_READS | rev | cut -d'/' -f1 | rev | cut -d'.' -f1)_trimmed.fastq.gz \
-s trimmed_singles_file.fastq -q 30 -l 79 -g

#Save trimmed reads to new variables
R1_TRIMMED=$(find $PWD -name "*R1_001_trimmed.fastq.gz" -print)
R2_TRIMMED=$(find $PWD -name "*R2_001_trimmed.fastq.gz" -print)

# Create a table-histogram of the trimmed-read-length distribution
zcat $R1_TRIMMED | awk 'NR%4 == 2 {lengths[length($0)]++} END {for (l in lengths) {print l, lengths[l]}}' > read_lengths_R1_trimmed.txt
zcat $R2_TRIMMED | awk 'NR%4 == 2 {lengths[length($0)]++} END {for (l in lengths) {print l, lengths[l]}}' > read_lengths_R2_trimmed.txt

## Step 3: Merge paired-end reads and check aseembled read length distribution
## NOTE: to use pear load a python/conda environment: source activate path/to/python_env (see: https://rcc-uchicago.github.io/user-guide/midway23/software/apps_and_envs/python/#managing-environments)
source activate /project/joet1/santiagoherrera/python_env
pear -f $R1_TRIMMED -r $R2_TRIMMED -o DMS_Rep1-3_S0_L001_4 -n 100 -y 800M -j 10
conda deactivate #Leave the environment

ASSEMBLED=$(find $PWD -name "*.assembled.fastq" -print)
gzip -v9 $ASSEMBLED

# Create a table-histogram of the assembled-read-length distribution
#python plot_histogram_read_lengths.py DMS_Rep1_S0_L001_4.assembled.fastq.gz "DMS_Rep1_assemled"
zcat $ASSEMBLED | awk 'NR%4 == 2 {lengths[length($0)]++} END {for (l in lengths) {print l, lengths[l]}}' > read_length_assembled_reads.txt

## Step 4: Compute base composition of assembled reads (based on a random sample of 1M reads)
python Calc_base_composition_stats_assembled_reads_v3.py $ASSEMBLED "DMS_Rep1-3_S0_L001_4_assembled" 1000000

## Step 5: Split main assembled file into smaller files for demultiplex
# The assembled file is divided in chunks of 100M reads (18 splits total). Each chunck is demuxed independently in 8 batches
zcat $ASSEMBLED | wc -l
zcat $ASSEMBLED | split -l 400000000 -d --filter='gzip > $FILE.gz' - DMS_Rep1_S0_L001_4.assembled.fastq_
mkdir splits
mv DMS_Rep1_S0_L001_4.assembled.fastq_* splits/

## Step 6: Demultiplexing reads by bin barcode and background: 
# Use 'demultiplex_REBC_binBC_v9.py' script to sort reads from each split into 512 files of GFP+ cells (32 libs, 4 bins, 4 reps)
sbatch dmx_split00.sbatch
sbatch dmx_split01.sbatch
sbatch dmx_split02.sbatch
sbatch dmx_split03.sbatch
sbatch dmx_split04.sbatch
sbatch dmx_split05.sbatch
sbatch dmx_split06.sbatch
sbatch dmx_split07.sbatch
sbatch dmx_split08.sbatch
sbatch dmx_split09.sbatch
sbatch dmx_split10.sbatch
sbatch dmx_split11.sbatch
sbatch dmx_split12.sbatch
sbatch dmx_split13.sbatch
sbatch dmx_split14.sbatch
sbatch dmx_split15.sbatch
sbatch dmx_split16.sbatch
sbatch dmx_split17.sbatch

# Use 'demultiplex_REBC_binBC_v9.py' script to sort reads from each split into 32 files of GFP- cells
sbatch dmx_split00_neg.sbatch
sbatch dmx_split01_neg.sbatch
sbatch dmx_split02_neg.sbatch
sbatch dmx_split03_neg.sbatch
sbatch dmx_split04_neg.sbatch
sbatch dmx_split05_neg.sbatch
sbatch dmx_split06_neg.sbatch
sbatch dmx_split07_neg.sbatch
sbatch dmx_split08_neg.sbatch
sbatch dmx_split09_neg.sbatch
sbatch dmx_split10_neg.sbatch
sbatch dmx_split11_neg.sbatch
sbatch dmx_split12_neg.sbatch
sbatch dmx_split13_neg.sbatch
sbatch dmx_split14_neg.sbatch
sbatch dmx_split15_neg.sbatch
sbatch dmx_split16_neg.sbatch
sbatch dmx_split17_neg.sbatch

# Use 'demultiplex_spikein_binBC_v9.py' script to demultiplex the spike-in sequences per bin (148 files: 37 controls, 4 bins)
sbatch dmx_split00_spikein.sbatch
sbatch dmx_split01_spikein.sbatch
sbatch dmx_split02_spikein.sbatch
sbatch dmx_split03_spikein.sbatch
sbatch dmx_split04_spikein.sbatch
sbatch dmx_split05_spikein.sbatch
sbatch dmx_split06_spikein.sbatch
sbatch dmx_split07_spikein.sbatch
sbatch dmx_split08_spikein.sbatch
sbatch dmx_split09_spikein.sbatch
sbatch dmx_split10_spikein.sbatch
sbatch dmx_split11_spikein.sbatch
sbatch dmx_split12_spikein.sbatch
sbatch dmx_split13_spikein.sbatch
sbatch dmx_split14_spikein.sbatch
sbatch dmx_split15_spikein.sbatch
sbatch dmx_split16_spikein.sbatch
sbatch dmx_split17_spikein.sbatch

## Step 6: Count number of reads per variant per demultiplexed file
# bash script (aa_variants_counts_v3.sh)
sbatch run_AAcounts_parallel_v2.sbatch

## Step 7: Estimate meanF for all avriants in all backgrounds per replicate. 
Rscript Calc_meanF_variants_REP1.R
Rscript Calc_meanF_variants_REP2.R
Rscript Calc_meanF_variants_REP3.R
Rscript Calc_meanF_variants_REP4.R





