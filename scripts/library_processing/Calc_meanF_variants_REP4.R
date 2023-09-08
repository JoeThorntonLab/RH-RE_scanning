#!/usr/bin/env Rscript
# Estimate fluorescence for REP4
setwd("/project/joet1/santiagoherrera/DMS_RH_RE/221202_A00639_1438_AHLYYJDRX2-JT-JP-SHA-1s/demultiplexed_reads/merged_dmx_splits/REP4/AA_variant_counts")
library(tidyverse)

###################
# This script estimates mean fluorescence for each variant

# Example of usage:
# Rscript Calc_meanF_variants_v3.R
###################

#### FUNCTIONS ####
#This function simply joins and reorganizes the AA counts per bin per library into a single dataframe.
#It joins ALL AA counts across libraries into a dataframe for downstream meanF estimation.
bind_dfs <- function(names){
  df <- data.frame('AA_var'=NULL,'REBC'=NULL,'Count_b1'=NULL,'Count_b2'=NULL,'Count_b3'=NULL,'Count_b4'=NULL,'Count_total'=NULL)
  for(i in 1:length(names)){
    b1 <- paste(names[i],"_BB1_REP4_AA_var_count.csv",sep="")
    b2 <- paste(names[i],"_BB2_REP4_AA_var_count.csv",sep="")
    b3 <- paste(names[i],"_BB3_REP4_AA_var_count.csv",sep="")
    b4 <- paste(names[i],"_BB4_REP4_AA_var_count.csv",sep="")
    b1_df <- read.csv(b1,h=T)
    b2_df <- read.csv(b2,h=T)
    b3_df <- read.csv(b3,h=T)
    b4_df <- read.csv(b4,h=T)
    tmp <- full_join(b1_df,b2_df,by=c("AA_var","REBC")) %>% 
      full_join(b3_df,by=c("AA_var","REBC")) %>% 
      full_join(b4_df,by=c("AA_var","REBC")) %>%
      mutate(Count.x = replace(Count.x,is.na(Count.x),0)) %>% #Replace NAs counts with zero (AAvar found in one df but in the other))
      mutate(Count.y = replace(Count.y,is.na(Count.y),0)) %>%
      mutate(Count.x.x = replace(Count.x.x,is.na(Count.x.x),0)) %>%
      mutate(Count.y.y = replace(Count.y.y,is.na(Count.y.y),0)) %>%
      select(.,-c(BinBC.x,BinBC.y,BinBC.x.x,BinBC.y.y,REP.y,REP.x.x,REP.y.y)) %>%
      relocate(Count.x,.after = REBC) %>%
      rename(Count_b1 = Count.x) %>% rename(Count_b2 = Count.y) %>% rename(Count_b3 = Count.x.x) %>% rename(Count_b4 = Count.y.y) %>% rename(REP = REP.x) %>%
      mutate(Count_total = rowSums(select(.,c(Count_b1,Count_b2,Count_b3,Count_b4)))) %>%
      relocate(REP,.after = Count_total) %>% mutate(REP = replace(REP,is.na(REP),"REP4")) # replace NAs with REP#
    df <- rbind(df,tmp)
  }
  df
}

#This function stores the meanF per bin (as estimated from cytometry data during sorting using FlowJo. Mean of 4 measurements per hour)
bin_meanF <- function(){
  #log10(FITC-A/FSC-A^1.5)
  bin_1 <-  -4.68475
  bin_2 <-  -4.302706
  bin_3 <-  -3.548033
  bin_4 <-  -2.650779
  r <- c(bin_1,bin_2,bin_3,bin_4)
  names(r) <- c("Bin1","Bin2","Bin3","Bin4")
  r
}

#This function stores the total cell count per bin during the sorting experiment.
bin_Cellcount <- function(bin = "all"){
  if(bin == 1) c <- 16233786
  if(bin == 2) c <- 12616013
  if(bin == 3) c <- 822095
  if(bin == 4) c <- 579424
  if(bin == "all") c <- c(16233786,12616013,822095,579424)
  c
}

###############################
###############################
### Analsyses ####

#AA count files (per bin and library)
aa_files <- list.files(".",pattern = "AA_var_count")
aa_files_REBCs <- unique(gsub("_BB[0-9]_REP4_AA_var_count.csv","",aa_files))

aa_counts <- bind_dfs(aa_files_REBCs)

# Read count in bin B (across *all* backgrounds): total number of reads per bin
reads_b1 <- sum(aa_counts[,3]) 
reads_b2 <- sum(aa_counts[,4])
reads_b3 <- sum(aa_counts[,5])
reads_b4 <- sum(aa_counts[,6])

# Estimating number of cells per variant per background (library) for meanF
# Estimated fraction of cells with each AA variant in Bin B: Fraction of reads in bin B from variant i times the number of cells sorted in bin B.
aa_counts <- aa_counts %>% 
  rowwise() %>%
  mutate(cellCount_b1 = (Count_b1/reads_b1)*bin_Cellcount(1)) %>% #Estimated fraction of cells with each AA variant in Bin1
  mutate(cellCount_b2 = (Count_b2/reads_b2)*bin_Cellcount(2)) %>% #Estimated fraction of cells with each AA variant in Bin2
  mutate(cellCount_b3 = (Count_b3/reads_b3)*bin_Cellcount(3)) %>% #Estimated fraction of cells with each AA variant in Bin3
  mutate(cellCount_b4 = (Count_b4/reads_b4)*bin_Cellcount(4)) %>% #Estimated fraction of cells with each AA variant in Bin4
  mutate(meanF = (sum(c(cellCount_b1,cellCount_b2,cellCount_b3,cellCount_b4) * bin_meanF())/sum(c(cellCount_b1,cellCount_b2,cellCount_b3,cellCount_b4)))) # Estimate meanF for all variants (across bins)

#Export dataframe
outfile <- "DMS_meanF_rep4.csv"
write.csv(x = aa_counts,file = outfile,quote = F,row.names = F)
print(paste(outfile,"... complete!",sep=""))
