library(tidyverse)
library(dplyr)

all_snps <- read.table("chr3R_indels.txt", header = TRUE)

##=====================================nAchr beta 2 RA===================================

FGbn0004118_RA <- read.table("FGbn0004118.gtf", header = FALSE) %>%
  filter(V3 == "exon" & V19 == "nAChRbeta2-RA") %>%
  dplyr::select(V1, V2, V3, V4, V5, V13) %>%
  unique()

FGbn0004118_RB <- read.table("FGbn0004118.gtf", header = FALSE) %>%
  filter(V3 == "exon" & V19 == "nAChRbeta2-RB") %>%
  dplyr::select(V1, V2, V3, V4, V5, V13) %>%
  unique()

FGbn0004118_RA_SNPs <- tibble(chrom = character(), pos = numeric(), 
                         ref = character(), alt = character(), codon = integer(), founder = character())

for (i in 1:nrow(FGbn0004118_RA)) {
  
  exon <- FGbn0004118_RA[i,]
  start <- as.numeric(exon[4])
  end <- as.numeric(exon[5])
  index <- data.frame(pos = as.numeric(start:end))
  index$codon <- rep_len(1:3, nrow(index))
  exon_snps <- subset(all_snps, POS > start & POS < end)
  exon_snps$codon <- index$codon[match(exon_snps$POS, index$pos)]
  
  FGbn0004118_RA_SNPs <- rbind(FGbn0004118_RA_SNPs, exon_snps)
  
}

FGbn0004118_RA_count <- FGbn0004118_RA_SNPs %>%
  group_by(CHROM, POS) %>%
  summarise(n = n()) %>%
  filter(n == '1')

FGbn0004118_RA_unique <- subset(FGbn0004118_RA_SNPs, FGbn0004118_RA_SNPs$POS %in% FGbn0004118_RA_count$POS)
FGbn0004118_RA_B7 <- subset(FGbn000411_RA_SNPs, FOUNDER == "B7")

##=====================================nAchr beta 2 RB===================================

#24513667 A -> G, M -> V [A3]
#24513656 AT -> A, frameshift [A3]
#24509435 C -> A, V -> F [A3]

#24513637 A -> T, F -> I [A2, A3, A5, A6, A7, B1, B4, B5, B6]


check <- subset(FGbn0004118_RB_SNPs, POS == "24513637")

FGbn0004118_RB_SNPs <- tibble(chrom = character(), pos = numeric(), 
                              ref = character(), alt = character(), codon = integer(), founder = character())

for (i in 1:nrow(FGbn0004118_RB)) {
  
  exon <- FGbn0004118_RB[i,]
  start <- as.numeric(exon[4])
  end <- as.numeric(exon[5])
  index <- data.frame(pos = as.numeric(start:end))
  index$codon <- rep_len(1:3, nrow(index))
  exon_snps <- subset(all_snps, POS > start & POS < end)
  exon_snps$codon <- index$codon[match(exon_snps$POS, index$pos)]
  
  FGbn0004118_RB_SNPs <- rbind(FGbn0004118_RB_SNPs, exon_snps)
  
}

FGbn0004118_RB_count <- FGbn0004118_RB_SNPs %>%
  group_by(CHROM, POS) %>%
  summarise(n = n()) %>%
  filter(n == '1')

FGbn0004118_RB_unique <- subset(FGbn0004118_RB_SNPs, FGbn0004118_RB_SNPs$POS %in% FGbn0004118_RB_count$POS)
FGbn0004118_RB_B7 <- subset(FGbn000411_RB_SNPs, FOUNDER == "B7")


#=========================================nAchr alpha 1 RA ===============================================

FBgn0000036_RA <- read.table("FBgn0000036.gtf", header = FALSE) %>%
  filter(V3 == "exon" & V19 == "nAChRalpha1-RA") %>%
  dplyr::select(V1, V2, V3, V4, V5, V13) %>%
  unique()

FBgn0000036_RA_SNPs <- tibble(chrom = character(), pos = numeric(), 
                           ref = character(), alt = character(), codon = integer(), founder = character())

for (i in 1:nrow(FBgn0000036_RA)) {
  
  exon <- FBgn0000036_RA[i,]
  start <- as.numeric(exon[4])
  end <- as.numeric(exon[5])
  index <- data.frame(pos = as.numeric(start:end))
  index$codon <- rep_len(1:3, nrow(index))
  exon_snps <- subset(all_snps, POS > start & POS < end)
  exon_snps$codon <- index$codon[match(exon_snps$POS, index$pos)]
  
  FBgn0000036_RA_SNPs <- rbind(FBgn0000036_RA_SNPs, exon_snps)
  
}

FBgn0000036_RA_count <- FBgn0000036_RA_SNPs %>%
  group_by(CHROM, POS) %>%
  summarise(n = n()) %>%
  filter(n == '2')

FBgn0000036_RA_unique <- subset(FBgn0000036_RA_SNPs, FBgn0000036_RA_SNPs$POS %in% FBgn0000036_RA_count$POS)

#=========================================nAchr alpha 1 RB ===============================================

FBgn0000036_RB <- read.table("FBgn0000036.gtf", header = FALSE) %>%
  filter(V3 == "exon" & V19 == "nAChRalpha1-RB") %>%
  dplyr::select(V1, V2, V3, V4, V5, V13) %>%
  unique()

FBgn0000036_RB_SNPs <- tibble(chrom = character(), pos = numeric(), 
                              ref = character(), alt = character(), codon = integer(), founder = character())

for (i in 1:nrow(FBgn0000036_RB)) {
  
  exon <- FBgn0000036_RB[i,]
  start <- as.numeric(exon[4])
  end <- as.numeric(exon[5])
  index <- data.frame(pos = as.numeric(start:end))
  index$codon <- rep_len(1:3, nrow(index))
  exon_snps <- subset(all_snps, POS > start & POS < end)
  exon_snps$codon <- index$codon[match(exon_snps$POS, index$pos)]
  
  FBgn0000036_RB_SNPs <- rbind(FBgn0000036_RB_SNPs, exon_snps)
  
#=========================================nAchr alpha 1 RC ===============================================
  
  FBgn0000036_RC <- read.table("FBgn0000036.gtf", header = FALSE) %>%
    filter(V3 == "exon" & V19 == "nAChRalpha1-RC") %>%
    dplyr::select(V1, V2, V3, V4, V5, V13) %>%
    unique()
  
  FBgn0000036_RC_SNPs <- tibble(chrom = character(), pos = numeric(), 
                                ref = character(), alt = character(), codon = integer(), founder = character())
  
  for (i in 1:nrow(FBgn0000036_RC)) {
    
    exon <- FBgn0000036_RC[i,]
    start <- as.numeric(exon[4])
    end <- as.numeric(exon[5])
    index <- data.frame(pos = as.numeric(start:end))
    index$codon <- rep_len(1:3, nrow(index))
    exon_snps <- subset(all_snps, POS > start & POS < end)
    exon_snps$codon <- index$codon[match(exon_snps$POS, index$pos)]
    
    FBgn0000036_RC_SNPs <- rbind(FBgn0000036_RC_SNPs, exon_snps)
  
}

#==================================================nAchr alpha 2===============================
  
FBgn0000039_RA <- read.table("FBgn0000039.gtf", header = FALSE) %>%
  filter(V3 == "exon" & V19 == "nAChRalpha2-RA") %>%
  dplyr::select(V1, V2, V3, V4, V5, V13) %>%
  unique()

FBgn0000039_RA_SNPs <- tibble(chrom = character(), pos = numeric(), 
                           ref = character(), alt = character(), codon = integer(), founder = character())

for (i in 1:nrow(FBgn0000039_RA)) {
  
  exon <- FBgn0000039_RA[i,]
  start <- as.numeric(exon[4])
  end <- as.numeric(exon[5])
  index <- data.frame(pos = as.numeric(start:end))
  index$codon <- rep_len(1:3, nrow(index))
  exon_snps <- subset(all_snps, POS > start & POS < end)
  exon_snps$codon <- index$codon[match(exon_snps$POS, index$pos)]
  
  FBgn0000039_RA_SNPs <- rbind(FBgn0000039_RA_SNPs, exon_snps)
  
}

FBgn0000039_RA_count <- FBgn0000039_RA_SNPs %>%
  group_by(CHROM, POS) %>%
  summarise(n = n()) %>%
  filter(n == '1')

FBgn0000039_RA_unique <- subset(FBgn0000039__RASNPs, FBgn0000039__RASNPs$POS %in% FBgn0000039__RAcount$POS)

FBgn0000039_RA_B7 <- subset(FBgn0000039_RA_SNPs, FOUNDER == "B7")

#====================================================Cftr======================================================


FBgn0039207_RA <- read.table("FBgn0039207.gtf", header = FALSE) %>%
  filter(V3 == "exon" & V19 == "Cftr-RA") %>%
  dplyr::select(V1, V2, V3, V4, V5, V13) %>%
  unique()

FBgn0039207_RA_SNPs <- tibble(chrom = character(), pos = numeric(), 
                           ref = character(), alt = character(), codon = integer(), founder = character())

for (i in 1:nrow(FBgn0039207_RA)) {
  
  exon <- FBgn0039207_RA[i,]
  start <- as.numeric(exon[4])
  end <- as.numeric(exon[5])
  index <- data.frame(pos = as.numeric(start:end))
  index$codon <- rep_len(1:3, nrow(index))
  exon_snps <- subset(all_snps, POS > start & POS < end)
  exon_snps$codon <- index$codon[match(exon_snps$POS, index$pos)]
  
  FBgn0039207_RA_SNPs <- rbind(FBgn0039207_RA_SNPs, exon_snps)
  
}

FBgn0039207_RA_count <- FBgn0039207_RA_SNPs %>%
  group_by(CHROM, POS) %>%
  summarise(n = n()) %>%
  filter(n == '2')

FBgn0039207_RA_unique <- subset(FBgn0039207_RA_SNPs, FBgn0039207_RA_SNPs$POS %in% FBgn0039207_RA_count$POS)

FBgn0039207_B7 <- subset(FBgn0039207_SNPs, FOUNDER == "B7")
check <- subset(FBgn0039207_SNPs, POS == "24540374")

#24541500 first pos, A3, A7, B2, B7
#24539888 A2-A7, B2, B4, B7
