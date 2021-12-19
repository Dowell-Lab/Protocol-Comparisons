#!/usr/bin/env Rscript
## Author: Samuel Hunter <sahu0957@colorado.edu>
## Adapted from Analysis-Flow Pipeline from Zachary Maas
## Maintainer: Samuel Hunter <sahu0957@colorado.edu>
##
######################################################################
##
### Commentary:
##
## Plots Metagenes for multiple conditions and replicates around TSSs
##
######################################################################
##
### Code:
install.packages("magrittr")
R.version
suppressMessages(library("tidyverse"))
suppressMessages(library("ggplot2"))
suppressMessages(library("argparse"))
suppressMessages(library("reshape2"))
suppressMessages(library("digest"))
library(dplyr)

rm(list=ls())
workdir='/path/to/workdir/'
setwd(workdir)
sense <- 'sensecounts.txt'
antisense <- 'antisensecounts.txt'

num_bins <- 100
region_length<-10

outfile <- 'metagene.pdf'
countsSense <- read_delim(sense, delim="\t", col_names = TRUE,skip = 1)
countsAntiSense <- read_delim(antisense, delim="\t", col_names =TRUE,skip=1)

names(countsSense)

df_sense <- countsSense %>% separate(Geneid, into = c("geneid", "coord"), sep="/")
df_sense$coord <- as.numeric(df_sense$coord)
df_sense <- df_sense %>% arrange(geneid, coord)

names(df_sense)
# column names of samples to be tested from count files. Change these to match directorys and file names
samples_of_interest <- c("../remap_reanalysis/hg38_nomultimaps/mapped/bams/testing//GRO15_DMSO_062111_042215.sorted.bam",
                         "../remap_reanalysis/hg38_nomultimaps/mapped/bams/testing//GRO18_DMSO10_062111_042215.sorted.bam",
                         "../remap_reanalysis/hg38_nomultimaps/mapped/bams/testing//GRO20_DMSO100_062111_042215.sorted.bam",
                         "../remap_reanalysis/hg38_nomultimaps/mapped/bams/testing//i13_AGTCAA_L008_R1_001.sorted.bam",
                         "../remap_reanalysis/hg38_nomultimaps/mapped/bams/testing//i17_GTAGAG_L008_R1_001.sorted.bam",
                         "../remap_reanalysis/hg38_nomultimaps/mapped/bams/testing//SRR1105736.sorted.bam",
                         "../remap_reanalysis/hg38_nomultimaps/mapped/bams/testing//SRR1105737.sorted.bam",
                         "../remap_reanalysis/hg38_nomultimaps/mapped/bams/testing//SRR8429046.sorted.bam",
                         "../remap_reanalysis/hg38_nomultimaps/mapped/bams/testing//SRR8429047.sorted.bam",
                         "../remap_reanalysis/hg38_nomultimaps/mapped/bams/testing//SRR8429054.sorted.bam",
                         "../remap_reanalysis/hg38_nomultimaps/mapped/bams/testing//SRR8429055.sorted.bam"
                         )

# classification of samples to be tested from count files (should match up w/ samples_of_interest)
preps_of_interest <- c("GRO-RPR",
                       "GRO-RPR",
                       "GRO-RPR",
                       "GRO-LIG",
                       "GRO-LIG",
                       "GRO-CIRC",
                       "GRO-CIRC",
                       "PUBLIC-GRO-RPR",
                       "PUBLIC-GRO-RPR",
                       "PUBLIC-GRO-RPR",
                       "PUBLIC-GRO-RPR")


identity_df <- as.data.frame(cbind(samples_of_interest,preps_of_interest))
vital <- c("geneid",
           "coord",
           "Chr",                                                                                      
           "Start",                                                                                    
           "End",                                                                                      
           "Strand",                                                                                   
           "Length")

df_sense <- unique(df_sense[,colnames(df_sense) %in% identity_df$samples_of_interest | colnames(df_sense) %in% vital])
head(df_sense)

df_antisense <- countsAntiSense %>% separate(Geneid, into = c("geneid", "coord"), sep = "/")
df_antisense$coord <- as.numeric(df_antisense$coord)
df_antisense <- df_antisense %>% arrange(geneid, coord)

df_antisense <- unique(df_antisense[,colnames(df_antisense) %in% identity_df$samples_of_interest | colnames(df_antisense) %in% vital])


# melt dfs into signal by coord
meltify_frame <- function(frame) {
    new_frame <- frame %>%
        subset(select = -c(geneid, Chr, Start, End)) %>%
        melt(id = c("coord", "Strand", "Length")) %>%
        mutate(condition = variable,
               coord = ifelse(Strand == '-', (num_bins - 1) - coord, coord),
               reads = value) %>%
        subset(select = -c(variable, value)) %>%
        as_tibble()
    return(new_frame)
}

# same as above, but for the opposite strand (only necessary depending on 
# how you separate the strands to count the metagenes. For example, with Tfit, I need this function because
# bidirectionals always have the same shape, so I essentially count each region twice
# for each strand. Genes can be + or - strand though, so we only use meltify_frame to account for - strand genes)
antimeltify_frame <- function(frame) {
  new_frame <- frame %>%
    subset(select = -c(geneid, Chr, Start, End)) %>%
    melt(id = c("coord", "Strand", "Length")) %>%
    mutate(condition = variable,
           coord = ifelse(Strand == '+', (num_bins - 1) - coord, coord),
           reads = value) %>%
    subset(select = -c(variable, value)) %>%
    as_tibble()
  return(new_frame)
}

head(df_sense)
sense_melt <- meltify_frame(df_sense)
antisense_melt <- meltify_frame(df_antisense)
head(sense_melt)
head(antisense_melt)

## Normalize the counts for each region by TPM
preps <- unique(sense_melt$condition)
preps

# per sample TPM manual
tpm_normalize <- function(fst, snd, sample) {
    ## First, calculate RPK for the strand only
    RPK_sample <- (fst[[sample]] / (fst$Length / (10 ^ 3)))
    ## Then, calculate RPK for both strands (for scaling factor)
    RPK <- ((fst[[sample]] + snd[[sample]]) /
            (fst$Length / (10 ^ 3)))
    ## Then, calculate the scaling factor
    scale <- sum(RPK) / 1000000
    ## Divide RPK values by scaling factor
    out <- RPK_sample / scale
    return(out)
}
save_sense <- sense_melt
save_antisense <- antisense_melt

sense_tpm_melt <- data.frame()
antisense_tpm_melt <- data.frame()

sense_melt <- save_sense
antisense_melt <- save_antisense
head(sense_melt)

# Automated TPM by condition (if they're named the same, use manual method above for replicates, 
# and average them later!). Note that the region length is always 10, so only the depth adjustment actually matters
for (i in preps){
  tmp_sense_df <- sense_melt[sense_melt$condition==i,]
  tmp_antisense_df <- antisense_melt[antisense_melt$condition==i,]
  rpk <- ((tmp_sense_df$reads + tmp_antisense_df$reads))/(region_length/1000)
  scale <- sum(rpk) / 1000000 
  
  tmp_sense_df$tpm_reads <- (tmp_sense_df$reads)/(region_length/1000)/scale
  tmp_antisense_df$tpm_reads <- -(tmp_antisense_df$reads)/(region_length/1000)/scale
  
  sense_tpm_melt <- rbind(sense_tpm_melt,tmp_sense_df)
  antisense_tpm_melt <- rbind(antisense_tpm_melt,tmp_antisense_df)
}

# After we tpmify everything, we should rename all of the samples to be reflective of the preps instead


sense_melt <- merge(sense_tpm_melt, identity_df, by.x=4, by.y=1, all.x=TRUE)
antisense_melt <- merge(antisense_tpm_melt, identity_df, by.x=4, by.y=1, all.x=TRUE)

# There are some towers near the ends of the metagene that are really driving signal. 
# Remove regions by a simple filter of TPM > 100 if they're near either end of the metagene

sense_melt <- sense_melt[!(sense_melt$tpm_reads > 100 & sense_melt$coord > 75),]
sense_melt <- sense_melt[!(sense_melt$tpm_reads > 100 & sense_melt$coord < 25),]

antisense_melt <- antisense_melt[!(antisense_melt$tpm_reads < -100),]

quantile(sense_melt$tpm_reads)

# Generqte summary stats, based on replicates, read signal, coord
summary_stats <- function(frame) {
    new_frame <- group_by(frame, coord, preps_of_interest) %>%
        summarize(mu = mean(tpm_reads),
                  var = sd(tpm_reads) / sqrt(length(tpm_reads))) %>%
        mutate(min = mu - var, max = mu + var)
    return(new_frame)
}

sense_melt_final <- summary_stats(sense_melt)
antisense_melt_final <- summary_stats(antisense_melt)

hex_to_int <- function(h) {
    xx = strsplit(tolower(h), "")[[1L]]
    pos = match(xx, c(0L:9L, letters[1L:6L]))
    sum((pos - 1L) * 16^(rev(seq_along(xx) - 1)))
}
sample_color <- function(name) {
    hue <- (hex_to_int(digest(name, algo='xxhash32')) %% 1009) / 1009
    new_color <- hsv(hue, 1, 0.75)
    return(new_color)
}
#install.packages("ggthemes")
library('ggthemes')

head(sense_melt_final)


my_comparisons <- c("GRO-LIG","GRO-RPR","GRO-CIRC")

sense_melt_final$coord <- (sense_melt_final$coord*region_length)-((num_bins*region_length)/2)
antisense_melt_final$coord <- (antisense_melt_final$coord*region_length)-((num_bins*region_length)/2)


test_sense <- sense_melt_final[sense_melt_final$preps_of_interest %in% my_comparisons,]
test_antisense <- antisense_melt_final[antisense_melt_final$preps_of_interest %in% my_comparisons,]

sense_melt_final
ggplot() + theme_tufte(base_size = 17) +
  scale_color_manual(values =  c(   "GRO-CIRC"="#1CBF00",
                                     "GRO-LIG"="#0024BF",
                                     "GRO-RPR"="#bf2419"
  )) +
  scale_fill_manual(values =  c(    "GRO-CIRC"="#1CBF00",
                                    "GRO-LIG"="#0024BF",
                                    "GRO-RPR"="#bf2419"
  )) +
 # geom_line(data = sense_melt_final, aes(x = coord, y = mu, color = preps_of_interest)) +
#  geom_line(data = antisense_melt_final, aes(x = coord, y = mu, color = preps_of_interest)) +
  stat_smooth(data = test_sense,method="gam", aes(x = coord, y = mu,color=preps_of_interest,fill=preps_of_interest)) +
  stat_smooth(data = test_antisense,method="gam", aes(x = coord, y = mu,color=preps_of_interest,fill=preps_of_interest)) +
    geom_hline(yintercept = 0) +
    labs(title = paste0("Comparing Preps:5' End"),
         x = "Relative Position (bp)", y = "Normalized Read Depth",
         color = "Condition", fill = "SD of Mean") +
  geom_vline(xintercept = 0) 


my_comparisons <- c("GRO-RPR","PUBLIC-GRO-RPR")

test_sense <- sense_melt_final[sense_melt_final$preps_of_interest %in% my_comparisons,]
test_antisense <- antisense_melt_final[antisense_melt_final$preps_of_interest %in% my_comparisons,]

ggplot() + theme_tufte(base_size = 17) +
  # geom_line(data = sense_melt_final, aes(x = coord, y = mu, color = preps_of_interest)) +
  #  geom_line(data = antisense_melt_final, aes(x = coord, y = mu, color = preps_of_interest)) +
  stat_smooth(data = test_sense,method="gam",color="#bf2419",fill="#bf2419", aes(x = coord, y = mu,linetype=preps_of_interest)) +
  stat_smooth(data = test_antisense,method="gam", color="#bf2419",fill="#bf2419", aes(x = coord, y = mu,linetype=preps_of_interest)) +
  geom_hline(yintercept = 0) +
  labs(title = paste0("Comparing Preps:5' End"),
       x = "Relative Position (bp)", y = "Normalized Read Depth",
       color = "Condition", fill = "SD of Mean") +
  geom_vline(xintercept = 0) 

ggsave("metagene.png",device="png")
ggsave("metagene.svg",device="svg")
## Output Data Frame DataFrame
