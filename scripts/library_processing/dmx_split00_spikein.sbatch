#!/bin/bash
#SBATCH --job-name=dmx_s00_SI
#SBATCH --output=/project/joet1/santiagoherrera/DMS_RH_RE/221202_A00639_1438_AHLYYJDRX2-JT-JP-SHA-1s/demultiplexed_reads/dmx_split_00/demultiplex_spikein_log3.txt
#SBATCH --error=dmx_s00_SI.err
#SBATCH --time=10:00:00
#SBATCH --partition=caslake
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=14
#SBATCH --mem-per-cpu=2000
#SBATCH --account=pi-joet1


# This script demultiplxes split_00 into all 148 files of spike-in controls (37 controls, 4 bins)

module load python

start=`date +%s.%N`

DMXDIR=/project/joet1/santiagoherrera/DMS_RH_RE/221202_A00639_1438_AHLYYJDRX2-JT-JP-SHA-1s/demultiplexed_reads/dmx_split_00 #Modify this
SCRIPTS=/project/joet1/santiagoherrera/DMS_RH_RE/scripts
FASTQDIR=/project/joet1/santiagoherrera/DMS_RH_RE/221202_A00639_1438_AHLYYJDRX2-JT-JP-SHA-1s/FastQ/splits

cd $DMXDIR

FASTQ=$FASTQDIR/DMS_Rep1-3_S0_L001_4.assembled.fastq_00.gz #Modify this
BBC=../bin_barcodes_forward_R1.txt
REBC=../spike_in_refSeq_forward.txt
DMX=../spikein_dmx_files_4.txt

python $SCRIPTS/demultiplex_spikein_binBC_v9.py $FASTQ $BBC $REBC $DMX

end=`date +%s.%N`
runtime=$( echo "$end - $start" | bc -l )

echo $'·········\ndemultiplexing DONE!\n·········'
echo "$runtime seconds total time"
