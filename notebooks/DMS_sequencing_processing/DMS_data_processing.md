DMS sequencing data processing
================
Santiago Herrera
2022-08-22

This is the workflow for cleaning, assembling and demultiplexing DMS reads based on DBD background, REBC (the synonymous mutations identifying the RE strain), and BBC (FACS bin barcode: enrichment GFP- bin, and four GFP+ bins per experimental replicate).

The Python and R scripts used during the workflow are found in the `scripts/library_processing` folder. The script `DMS_data_processing.sh` is a bash script to run everything on a cluster. The same bash script is used to process sequencing data from NextSeq (Rep1) and NovaSeq (Reps 1-4).

    #!/bin/bash
    # Load modules in cluster
    module load python
    module load sickle/1.33
    module load R

    ## Global variables:

    # NovaSeq sequencing data
    SEQ_FOLDER=/project/joet1/RH-RE_scanning/data/sequencing/NovaSeq
    R1_READS=$SEQ_FOLDER/NovaSeq/JT-JP-SHA-1s_S1_R1_001.fastq.gz
    R2_READS=$SEQ_FOLDER/NovaSeq/JT-JP-SHA-1s_S1_R2_001.fastq.gz

    # NextSeq sequencing data
    #SEQ_FOLDER=/project/joet1/RH-RE_scanning/data/sequencing/NextSeq
    #R1_READS=$SEQ_FOLDER/NextSeq/DMS_Rep1_S0_L001_4_R1_001.fastq.gz
    #R2_READS=$SEQ_FOLDER/NextSeq/DMS_Rep1_S0_L001_4_R2_001.fastq.gz

    cd $SEQ_FOLDER

First, we will check how many sequencing reads we got. Replicate 1 binned sort libraries were sequenced on NextSeq High Output run. The remaining libraries were sequenced on a NovaSeq S1 2 x 100bp flow cell. We used standard read primers and 86 cycles for read 1 and 80 cycles for read 2. This enabled us to sequence through almost the entire scanning and RE barcode region (short one nucleotide for read 2) in both directions (forward: R1 and reverse: R2).

    ## Step 1: Check quality statistics of the fastq files
    # Create a table-histogram of the read-length distribution
    zcat $R1_READS | awk 'NR%4 == 2 {lengths[length($0)]++} END {for (l in lengths) {print l, lengths[l]}}' > read_lengths_R1.txt
    zcat $R2_READS | awk 'NR%4 == 2 {lengths[length($0)]++} END {for (l in lengths) {print l, lengths[l]}}' > read_lengths_R2.txt

``` r
# histogram
read_lengths_R1_nova <- read.table("./NovaSeq/read_lengths_R1.txt",header = F) %>% dplyr::rename(length = V1, freq = V2) %>% mutate(file="R1",seq="NovaSeq")
read_lengths_R2_nova <- read.table("./NovaSeq/read_lengths_R2.txt",header = F) %>% dplyr::rename(length = V1, freq = V2) %>% mutate(file="R2",seq="NovaSeq")

read_lengths_R1_next <- read.table("./NextSeq/read_lengths_R1.txt",header = F) %>% dplyr::rename(length = V1, freq = V2) %>% mutate(file="R1",seq="NextSeq")
read_lengths_R2_next <- read.table("./NextSeq/read_lengths_R2.txt",header = F) %>% dplyr::rename(length = V1, freq = V2) %>% mutate(file="R2",seq="NextSeq")


knitr::kable(rbind(read_lengths_R1_nova,read_lengths_R2_nova,read_lengths_R1_next,read_lengths_R2_next))
```

|  length|        freq| file | seq     |
|-------:|-----------:|:-----|:--------|
|     118|  2053763453| R1   | NovaSeq |
|     118|  2053763453| R2   | NovaSeq |
|      86|   127490257| R1   | NextSeq |
|      80|   127490257| R2   | NextSeq |

As expected, we got about 2 x 10^9 reads with a length of 118bp for the NovaSeq run, and 1.2 x 10^8 reads with a length of ~80 bp for the NextSeq run.

Now, let's trim the reads based on Phred quality and length. A per-base Phred score of 30 corresponds to an error probability of ~0.001, so we will retain reads that have an average Phred score of 30. We will also only retain reads with a minimum length of 79 bp, because this is the minimum length for a forward or reverse read to cover the DMS sites.

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

``` r
# plot histograms
trimmed_lengths_R1_nova <- read.table("./NovaSeq/read_lengths_R1_trimmed.txt",header = F) %>% dplyr::rename(length = V1, freq = V2) %>% mutate(file="R1",seq="NovaSeq")
trimmed_lengths_R2_nova <- read.table("./NovaSeq/read_lengths_R2_trimmed.txt",header = F) %>% dplyr::rename(length = V1, freq = V2) %>% mutate(file="R2",seq="NovaSeq")

trimmed_lengths_R1_next <- read.table("./NextSeq/read_lengths_R1_trimmed.txt",header = F) %>% dplyr::rename(length = V1, freq = V2) %>% mutate(file="R1",seq="NextSeq")
trimmed_lengths_R2_next <- read.table("./NextSeq/read_lengths_R2_trimmed.txt",header = F) %>% dplyr::rename(length = V1, freq = V2) %>% mutate(file="R2",seq="NextSeq")


knitr::kable(rbind(trimmed_lengths_R1_nova,trimmed_lengths_R2_nova,trimmed_lengths_R1_next,trimmed_lengths_R2_next) %>% 
               group_by(file,seq) %>% reframe(reads_kept = sum(freq)))
```

| file | seq     |  reads\_kept|
|:-----|:--------|------------:|
| R1   | NextSeq |    120710413|
| R1   | NovaSeq |   1849725285|
| R2   | NextSeq |    120710413|
| R2   | NovaSeq |   1849725285|

``` r
mean_lengths <- rbind(trimmed_lengths_R1_nova,trimmed_lengths_R2_nova,trimmed_lengths_R1_next,trimmed_lengths_R2_next) %>% 
  group_by(file,seq) %>% mutate(total = sum(freq), fr = freq/total) %>% reframe(mean_rl = sum(length*fr))

rbind(trimmed_lengths_R1_nova,trimmed_lengths_R2_nova,trimmed_lengths_R1_next,trimmed_lengths_R2_next) %>%
  ggplot(aes(x=length,y=freq)) + geom_bar(stat="identity") +
  labs(x="Trimmed read length",y="Frequency") +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
              labels = trans_format("log10", math_format(10^.x))) +
  geom_vline(data=mean_lengths, aes(xintercept = mean_rl),linetype=2,col="red") +
  theme_classic() + 
  theme(axis.title = element_text(size=13)) +
  facet_wrap(seq~file, scales = "free_y")
```

![](DMS_data_processing_files/figure-markdown_github/unnamed-chunk-2-1.png)

After the trimming process, we retained ~1.8 x 10^9 (90%) of the NovaSeq reads, and ~1.2 x 10^8 (94%) of the NextSeq reads.

Now, we will assemble the trimmed reads to get full length amplicon sequences, with a minimum assembly length of 100 nucleotides.

    ## Step 3: Merge paired-end reads and check aseembled read length distribution
    ## NOTE: to use pear load a python/conda environment: source activate path/to/python_env (see: https://rcc-uchicago.github.io/user-guide/midway23/software/apps_and_envs/python/#managing-environments)
    source activate /project/joet1/santiagoherrera/python_env
    pear -f $R1_TRIMMED -r $R2_TRIMMED -o DMS_Rep1-3_S0_L001_4 -n 100 -y 800M -j 10
    conda deactivate #Leave the environment

    ASSEMBLED=$(find $PWD -name "*.assembled.fastq" -print)
    gzip -v9 $ASSEMBLED

    # Create a table-histogram of the assembled-read-length distribution
    zcat $ASSEMBLED | awk 'NR%4 == 2 {lengths[length($0)]++} END {for (l in lengths) {print l, lengths[l]}}' > read_length_assembled_reads.txt

``` r
# plot histogram
assembled_lengths_nova <- read.table("./NovaSeq/read_length_assembled_reads.txt",header = F) %>% dplyr::rename(length = V1, freq = V2) %>% mutate(seq="NovaSeq")
assembled_lengths_next <- read.table("./NextSeq/read_length_assembled_reads.txt",header = F) %>% dplyr::rename(length = V1, freq = V2) %>% mutate(seq="NextSeq")

knitr::kable(rbind(assembled_lengths_nova,assembled_lengths_next) %>% group_by(seq) %>% reframe(file="assembled",total_assembled=sum(freq)))
```

| seq     | file      |  total\_assembled|
|:--------|:----------|-----------------:|
| NextSeq | assembled |         118878340|
| NovaSeq | assembled |        1749639397|

``` r
mean_lengths <- rbind(assembled_lengths_nova,assembled_lengths_next) %>% group_by(seq) %>% 
  mutate(total = sum(freq), fr = freq/total) %>% reframe(mean = sum(length*fr))

rbind(assembled_lengths_nova,assembled_lengths_next) %>%
  ggplot(aes(x=length,y=freq)) + geom_bar(stat="identity") +
  labs(x="Assembled read length",y="Frequency") +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
              labels = trans_format("log10", math_format(10^.x))) +
  geom_vline(data=mean_lengths,aes(xintercept = mean),linetype=2,col="red") +
  theme_classic() + 
  theme(axis.title = element_text(size=13)) +
  facet_wrap(.~seq)
```

![](DMS_data_processing_files/figure-markdown_github/unnamed-chunk-3-1.png)

We see that we assembled successfully ~1.2 x 10^8 (&gt;98%) and 1.75 x 10^9 (&gt;94%) of the trimmed reads for NextSeq and NovaSeq, respectively.

------------------------------------------------------------------------

Now, we will do a quick check of the quality of the sequenced libraries by mapping a sample of 1M reads from each run to the reference sequence of AncSR1 or AncSR2. This will allow us to check the base composition per site: there should only be amino acid variation at the background sites that differ between AncSR1 and AncSR2, the sites corresponding to the REBC, and at the DMS sites.

For the background differences, there should only be 2 nucleotide states; for the REBC sites there are some sites with only two states and other with four; and for the DMS sites, we used NNS codons, that is, the first and second position of the codon should contain every nucleotide in a frequency of ~0.25, and the third position of the codon should have only G/C at a frwquency of ~0.5.

    ## Step 4: Compute base composition of assembled reads (based on a sample of 1M reads)
    ASSEMBLED=$(find $PWD -name "*.assembled.fastq.gz" -print)
    python Calc_base_composition_stats_assembled_reads_v4.py $ASSEMBLED "DMS_Rep1-3_S0_L001_4_assembled" 

``` r
#plot base composition
bc_nova <- read.csv("./NovaSeq/DMS_Rep1-3_S0_L001_4_assembled_base_composition.csv",h=T)
rownames(bc_nova) <- paste0("s",1:117)
bc_next <- read.csv("./NextSeq/DMS_Rep1_S0_L001_4_assembled_base_composition.csv",h=T)
rownames(bc_next) <- paste0("s",1:117)

DMS <- c(55:60,64:69)
REBC <- c(42,45,46,48,49,50,51,54,63)
Background_subs <- c(4:7,13,14,19,21,32,47:49,68,82,85,88,96,97,100,106,108,112,114,116)

DMS_sites_df <- data.frame(site = 1:117) %>% mutate(type = case_when(site %in% DMS ~ "DMS",
                                                                     site %in% REBC ~ "REBC",
                                                                     site %in% Background_subs ~ "Bg subs",
                                                                     T ~ "fixed")) %>% select(type)

row.names(DMS_sites_df) <- rownames(bc_nova)
DMS_sites_colors <- list(type = c("fixed" = "gray60","DMS" = "yellow","REBC" = "green","Bg subs" = "red"))
  
pheatmap(t(bc_nova),cluster_rows = F,cluster_cols = F, border_color="white",
         fontsize_row=10,fontsize_col=10,cellheight=10,main = "NovaSeq",color=rev(viridis::inferno(20)),
         annotation_col = DMS_sites_df,annotation_colors = DMS_sites_colors,show_colnames=F)
```

![](DMS_data_processing_files/figure-markdown_github/unnamed-chunk-4-1.png)

``` r
# Show table per site with the frequencies of NNS at the DMS sites
print("NovaSeq dataset:")
```

    ## [1] "NovaSeq dataset:"

``` r
knitr::kable(bc_nova %>% mutate(site=row_number()) %>% filter(site %in% c(55:60,64:69)) %>% 
               mutate(class = rep(c("N","N","S"),4)))
```

|     |          A|          C|          G|          T|  site| class |
|:----|----------:|----------:|----------:|----------:|-----:|:------|
| s55 |  0.2306249|  0.1920567|  0.2581156|  0.3192028|    55| N     |
| s56 |  0.2391984|  0.1856869|  0.2615419|  0.3135728|    56| N     |
| s57 |  0.0044923|  0.4636597|  0.5255879|  0.0062602|    57| S     |
| s58 |  0.2298676|  0.1905384|  0.2565314|  0.3230626|    58| N     |
| s59 |  0.2263753|  0.1985171|  0.2583017|  0.3168059|    59| N     |
| s60 |  0.0002167|  0.4679917|  0.5213229|  0.0104686|    60| S     |
| s64 |  0.2051860|  0.1894537|  0.2683027|  0.3370576|    64| N     |
| s65 |  0.2028574|  0.1891922|  0.3006837|  0.3072667|    65| N     |
| s66 |  0.0105028|  0.4166228|  0.5727048|  0.0001696|    66| S     |
| s67 |  0.1977751|  0.1933770|  0.2943328|  0.3145151|    67| N     |
| s68 |  0.1955619|  0.1875868|  0.3077790|  0.3090723|    68| N     |
| s69 |  0.0001378|  0.2795345|  0.7098956|  0.0104321|    69| S     |

``` r
#plot base composition
pheatmap(t(bc_next),cluster_rows = F,cluster_cols = F, border_color="white",
         fontsize_row=10,fontsize_col=10,cellheight=10,main = "NextSeq",color=rev(viridis::inferno(20)),
         annotation_col = DMS_sites_df,annotation_colors = DMS_sites_colors,show_colnames=F)
```

![](DMS_data_processing_files/figure-markdown_github/unnamed-chunk-5-1.png)

``` r
# Show table per site with the frequencies of NNS at the DMS sites
print("NextSeq dataset:")
```

    ## [1] "NextSeq dataset:"

``` r
knitr::kable(bc_next %>% mutate(site=row_number()) %>% filter(site %in% c(55:60,64:69)) %>% 
               mutate(class = rep(c("N","N","S"),4)))
```

|     |          A|          C|          G|          T|  site| class |
|:----|----------:|----------:|----------:|----------:|-----:|:------|
| s55 |  0.2307965|  0.1950473|  0.2567112|  0.3174449|    55| N     |
| s56 |  0.2397776|  0.1871780|  0.2660186|  0.3070258|    56| N     |
| s57 |  0.0039661|  0.4681445|  0.5194862|  0.0084033|    57| S     |
| s58 |  0.2268906|  0.1919239|  0.2607838|  0.3204017|    58| N     |
| s59 |  0.2248092|  0.2077677|  0.2550504|  0.3123726|    59| N     |
| s60 |  0.0005531|  0.4752940|  0.5125606|  0.0115923|    60| S     |
| s64 |  0.2106303|  0.1932214|  0.2654668|  0.3306816|    64| N     |
| s65 |  0.2036623|  0.1934412|  0.2998380|  0.3030584|    65| N     |
| s66 |  0.0119856|  0.4211390|  0.5663933|  0.0004821|    66| S     |
| s67 |  0.1969826|  0.1961153|  0.2981978|  0.3087043|    67| N     |
| s68 |  0.1932869|  0.1901116|  0.3070559|  0.3095456|    68| N     |
| s69 |  0.0003592|  0.2823404|  0.7057683|  0.0115322|    69| S     |

This shows that the libraries seem to have a good quality based on the sample of reads analyzed.

Finally, we will do the demultiplexing of the reads into libraries by background, REBC and BBC.

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

After having demultiplxed the reads, we will now count the number of reads per DBD-RE complex present in each file.

    ## Step 6: Count number of reads per variant per demultiplexed file
    # bash script (aa_variants_counts_v3.sh)
    sbatch run_AAcounts_parallel_v2.sbatch

Finally, using the read counts per variant per bin per replicate per background, we can compute their mean fluorescence.

    ## Step 7: Estimate meanF for all avriants in all backgrounds per replicate. 
    Rscript Calc_meanF_variants_REP1.R
    Rscript Calc_meanF_variants_REP2.R
    Rscript Calc_meanF_variants_REP3.R
    Rscript Calc_meanF_variants_REP4.R