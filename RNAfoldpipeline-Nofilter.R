#!/usr/bin/env Rscript
############################################################
print("Loading the libraries")
############################################################
lapply(c("Biostrings", "tidyverse", "stringr", "dplyr"), library, character.only=T)
############################################################
print("Reading the input files")
############################################################
args = commandArgs(trailingOnly = T)
fastadirectory <- args[1]
temp <- readDNAStringSet(fastadirectory)
seqname <- names(temp)
sequence <- paste(temp)
fasta <- data.frame(seqname, sequence)

inputdir <- args[2]
myfiles <- list.files(inputdir, pattern="\\.csv$")
for (i in 1:length(myfiles)) assign(myfiles[i], read.csv(myfiles[i]))

top <- TopStrandSigPeaks_EnhancedCompare.csv
top <- top %>% add_column(Sequence..50nt.upstream = "", Utrack = "") %>%
  add_column(strand = "+", .after = "PeakCoord")

comp <- CompStrandSigPeaks_EnhancedCompare.csv
comp <- comp %>% add_column(Sequence..50nt.upstream = "", Utrack = "") %>%
  add_column(strand = "-", .after = "PeakCoord")
############################################################
print("Calling 50 nt Upstream of every TTS")
############################################################
for(i in 1:nrow(top)){
  for(j in 1:nrow(comp)){
    top$Sequence..50nt.upstream[i] <- substr(fasta$sequence, top$PeakCoord[i] - 50, top$PeakCoord[i])
    comp$Sequence..50nt.upstream[j] <- substr(fasta$sequence, comp$PeakCoord[j], comp$PeakCoord[j] + 50)
  }
}
comp$Sequence..50nt.upstream <- chartr("ATGC","TACG",comp$Sequence..50nt.upstream)
strReverse <- function(x)
  sapply(lapply(strsplit(x, NULL), rev), paste, collapse="")
comp$Sequence..50nt.upstream <-  strReverse(comp$Sequence..50nt.upstream)
############################################################
print("Running RNA fold")
############################################################
RNAfold <- rbind(top, comp)
RNAfold <- RNAfold %>% arrange(RNAfold$PeakCoord)
RNAfoldinput <- as.data.frame(RNAfold$Sequence..50nt.upstream) 
setwd(args[2])
write_csv(RNAfoldinput, file = "RNAfoldinput.csv")

system("RNAfold RNAfoldinput.csv > RNAfoldoutput.csv")

RNAfoldoutput <- read.csv("RNAfoldoutput.csv")

RNAfold <- RNAfold %>% add_column(FoldingPattern = "",RNAfolding = "", FreefoldingEnergy = "", .before = "Utrack")
RNAfoldoutput <- as.data.frame(RNAfoldoutput[c(-1),])
colnames(RNAfoldoutput) <- "sequence"
Foldingpattern <- as.data.frame(RNAfoldoutput[substr(RNAfoldoutput$sequence,1,1) %in% "A" | substr(RNAfoldoutput$sequence,1,1) %in% "U" | substr(RNAfoldoutput$sequence,1,1) %in% "G" | substr(RNAfoldoutput$sequence,1,1) %in% "C",])
colnames(Foldingpattern) <- "sequence"
RNAfold$FoldingPattern <- Foldingpattern$sequence
RNAfoldoutput <- as.data.frame(RNAfoldoutput[!substr(RNAfoldoutput$sequence,1,1) %in% "A" &!substr(RNAfoldoutput$sequence,1,1) %in% "U" & !substr(RNAfoldoutput$sequence,1,1) %in% "G" & !substr(RNAfoldoutput$sequence,1,1) %in% "C",])
colnames(RNAfoldoutput) <- "folding"
peakcoord <- as.data.frame(RNAfold[,1])
RNAfoldoutput <-  cbind(RNAfoldoutput, peakcoord)
colnames(RNAfoldoutput) <- c("sequence", "PeakCoord")

for(i in 1:nrow(RNAfold)){
  for(j in 1:nrow(RNAfoldoutput)){
    if(RNAfold$PeakCoord[i] == RNAfoldoutput$PeakCoord[j]){
      RNAfold$RNAfolding[i] <- RNAfoldoutput$sequence[j]
    }
  }
  RNAfold$FreefoldingEnergy[i] <- substr(RNAfold$RNAfolding[i], (nchar(RNAfold$RNAfolding[i])-8), (nchar(RNAfold$RNAfolding[i])))
  RNAfold$RNAfolding[i] <- strsplit(RNAfold$RNAfolding[i], split = " ", fixed = T)[[1]][1]
  RNAfold$FreefoldingEnergy[i] <- substr(RNAfold$FreefoldingEnergy[i], (nchar(RNAfold$FreefoldingEnergy[i])-6),(nchar(RNAfold$FreefoldingEnergy[i])-1))
}
############################################################
print("Calculating the U-track rich regions")
############################################################
for(x in 1:nrow(RNAfold)){
  breakCond <- 0
  
  while(breakCond == 0){
    
    for(y in 1:(nchar(RNAfold$FoldingPattern[x])-6)){
      
      if(substr(RNAfold$FoldingPattern[x],y,y) %in% "U"){
        
        tempString <- substr(RNAfold$FoldingPattern[x],y,(y+6))
        
        if(str_count(tempString,"U") >= 5){
          RNAfold$Utrack[x] <- "yes"
          breakCond <- 1
        }
      }
    }
    
    if(breakCond == 0){
      RNAfold$Utrack[x] <- "no"
      breakCond <- 1
    }
  }
}
RNAfold <- RNAfold %>% add_column(Case = "")
top <- RNAfold %>% filter(strand == "+")
comp <- RNAfold %>% filter(strand== "-")
############################################################
#  Calculating Cases
############################################################
case1.top <-  top %>% filter(Sig_Control == 1 & Sig_Cond1 == 1 & Sig_Cond2 == 1)
case2.top <- top %>% filter(Sig_Control == 1 & Sig_Cond1 == 0 & Sig_Cond2 == 0)
case3.top <- top %>% filter((Sig_Control == 0 & Sig_Cond1 == 0 & Sig_Cond2 == 1)| 
  (Sig_Control == 0 & Sig_Cond1 == 1 & Sig_Cond2 == 1))
case1.comp <-  comp %>% filter(Sig_Control == 1 & Sig_Cond1 == 1 & Sig_Cond2 == 1)
case2.comp <- comp %>% filter(Sig_Control == 1 & Sig_Cond1 == 0 & Sig_Cond2 == 0)
case3.comp <- comp %>% filter((Sig_Control == 0 & Sig_Cond1 == 0 & Sig_Cond2 == 1)| 
  (Sig_Control == 0 & Sig_Cond1 == 1 & Sig_Cond2 == 1))
for(i in 1:length(top$PeakCoord)){
  for(j in 1:length(case1.top$PeakCoord)){
    if(top$PeakCoord[i] == case1.top$PeakCoord[j]){
      top$Case[i] <- "case1"  
    }
  }
  for(k in 1:length(case2.top$PeakCoord)){
    if(top$PeakCoord[i] == case2.top$PeakCoord[k]){
      top$Case[i] <- "case2"  
    }
  }
  for(l in 1:length(case3.top$PeakCoord)){
    if(top$PeakCoord[i] == case3.top$PeakCoord[l]){
      top$Case[i] <- "case3"  
    }
  }
}
for(i in 1:length(comp$PeakCoord)){
  for(j in 1:length(case1.comp$PeakCoord)){
    if(comp$PeakCoord[i] == case1.comp$PeakCoord[j]){
      comp$Case[i] <- "case1"  
    }
  }
  for(k in 1:length(case2.comp$PeakCoord)){
    if(comp$PeakCoord[i] == case2.comp$PeakCoord[k]){
      comp$Case[i] <- "case2"  
    }
  }
  for(l in 1:length(case3.comp$PeakCoord)){
    if(comp$PeakCoord[i] == case3.comp$PeakCoord[l]){
      comp$Case[i] <- "case3"  
    }
  }
}
RNAfold <- rbind(top, comp)
RNAfold <- RNAfold %>% arrange(PeakCoord)
write_csv(RNAfold, file = "RNAfold_RawFile.csv")
############################################################
print("Annotating the genome")
############################################################
top <- RNAfold %>% filter(strand == "+")
comp <- RNAfold %>% filter(strand == "-")
genes <- list.files(inputdir,pattern="\\.genes$")
genes <- read.delim(genes[1], header = T, sep = "")
colnames(genes) <- c("genome", "from", "to", "strand", "name", "old.name", "length", "protein.Name", "symbol")

genes.top <- genes %>% filter(strand == "+")
genes.comp <- genes %>% filter(strand == "-")

# List of columns to convert to numeric
cols_to_numeric <- c("from", "to")
# Loop through the selected columns and convert them to numeric
genes.comp[cols_to_numeric] <- lapply(genes.comp[cols_to_numeric], as.numeric)
genes.top[cols_to_numeric] <- lapply(genes.top[cols_to_numeric], as.numeric)
comp$PeakCoord <- as.numeric(comp$PeakCoord)
top$PeakCoord <- as.numeric(top$PeakCoord)
top <- top %>% 
  add_column(TSSUTR ="", Locus ="", termUTR ="", .after = "PeakCoord")
comp <- comp %>% 
  add_column(TSSUTR ="", Locus ="", termUTR ="",.after = "PeakCoord")
for (i in 1:nrow(top)) {
  for (k in 1:nrow(genes.top)) {
    if (!is.na(genes.top$to[k]) && !is.na(genes.top$from[k + 1])) {
      #first annotate 3'UTR
      if (
        top$PeakCoord[i] > genes.top$to[k] &&
        top$PeakCoord[i] <= genes.top$from[k + 1] + 50 &&
        top$PeakCoord[i] <= genes.top$to[k] + 300
      ) {
        top$termUTR[i] <- genes.top$old.name[k]
      }
    }
    if (!is.na(genes.top$from[k]) && top$termUTR[i] == "") {
      #second annotate genic
      if (
        top$PeakCoord[i] >= genes.top$from[k] + 50 &&
        top$PeakCoord[i] <= genes.top$to[k]
      ) {
        top$Locus[i] <- genes.top$old.name[k]
      }
    }
    
    if (!is.na(genes.top$to[k]) && !is.na(genes.top$from[k + 1]) && top$termUTR[i] =="") {
      #third annotate 5'UTR
      if (
        top$PeakCoord[i] > genes.top$to[k] + 300 &&
        top$PeakCoord[i] < genes.top$from[k + 1] &&
        top$PeakCoord[i] > genes.top$from[k + 1] - 300
      ) {
        top$TSSUTR[i] <- genes.top$old.name[k + 1]
      }
    }
  }
}

for (i in 1:nrow(comp)) {
  for (k in 1:nrow(genes.comp)) {
    if (!is.na(genes.comp$to[k]) && !is.na(genes.comp$from[k + 1])) {
      #first annotate 3'UTR
      if (
        comp$PeakCoord[i] >= genes.comp$from[k+1] - 300 &&
        comp$PeakCoord[i] >= genes.comp$to[k] - 50 &&
        comp$PeakCoord[i] < genes.comp$from[k+1]
      ) {
        comp$termUTR[i] <- genes.comp$old.name[k+1]
      }
    }
    if (!is.na(genes.comp$from[k]) && comp$termUTR[i] == "") {
      #second annotate genic
      if (
        comp$PeakCoord[i] >= genes.comp$from[k] + 50 &&
        comp$PeakCoord[i] <= genes.comp$to[k]
      ) {
        comp$Locus[i] <- genes.comp$old.name[k]
      }
    }
    
    if (!is.na(genes.comp$to[k]) && !is.na(genes.comp$from[k + 1]) && comp$termUTR[i] == "") {
      #third annotate 5'UTR
      if (
        comp$PeakCoord[i] < genes.comp$from[k+1] - 300 &&
        comp$PeakCoord[i] > genes.comp$to[k] &&
        comp$PeakCoord[i] <= genes.comp$to[k] + 300
      ) {
        comp$TSSUTR[i] <- genes.comp$old.name[k]
      }
    }
  }
}

############################################################
print("Analysing the transcriptomic distribution")
############################################################
############### Gather Cases - TOP STRAND ##################
############################################################ 
case1.top <-  top %>% filter(Sig_Control == 1 & Sig_Cond1 == 1 & Sig_Cond2 == 1) 
case2.top <- top %>% filter(Sig_Control == 1 & Sig_Cond1 == 0 & Sig_Cond2 == 0)
case3.top <- top %>% filter((Sig_Control == 0 & Sig_Cond1 == 0 & Sig_Cond2 == 1) | 
                              (Sig_Control == 0 & Sig_Cond1 == 1 & Sig_Cond2 == 1))
#############################################################
############### Gather Cases - COMP STRAND ##################
############################################################# 
case1.comp <-  comp %>% filter(Sig_Control == 1 & Sig_Cond1 == 1 & Sig_Cond2 == 1)
case2.comp <- comp %>% filter(Sig_Control == 1 & Sig_Cond1 == 0 & Sig_Cond2 == 0)
case3.comp <- comp %>% filter((Sig_Control == 0 & Sig_Cond1 == 0 & Sig_Cond2 == 1) |
                                (Sig_Control == 0 & Sig_Cond1 == 1 & Sig_Cond2 == 1))
#########################################################
##########  Case 1 - Transcriptome Distribution #########
#########################################################
############  Top
genic.case1.top <- case1.top[((!case1.top$Locus == "") & 
                                (case1.top$TSSUTR == "") & 
                                (case1.top$termUTR == "")),]
TssUTR.case1.top <- case1.top[((!case1.top$TSSUTR == "") & 
                                 (case1.top$Locus == "") & 
                                 (case1.top$termUTR == "")),]
orphan.case1.top <- case1.top[((case1.top$Locus == "") & 
                                 (case1.top$TSSUTR == "") & 
                                 (case1.top$termUTR == "")),]
termUTR.case1.top <- case1.top[((case1.top$Locus == "") & 
                                  (case1.top$TSSUTR == "") &
                                  (!case1.top$termUTR == "")),]
other1 <- case1.top[((!case1.top$Locus == "") & (!case1.top$TSSUTR == "") & (!case1.top$termUTR == "")),]
other2 <- case1.top[((case1.top$Locus == "") & 
                       (!case1.top$TSSUTR == "") & 
                       (!case1.top$termUTR == "")),]
other3 <- case1.top[((!case1.top$Locus == "") & 
                       (case1.top$TSSUTR == "") & 
                       (!case1.top$termUTR == "")),]
other4 <- case1.top[((!case1.top$Locus == "") & 
                       (!case1.top$TSSUTR == "") & 
                       (case1.top$termUTR == "")),]
other.case1.top <- rbind(other1, other2, other3, other4)
############  Comp
genic.case1.comp <- case1.comp[((!case1.comp$Locus == "") & 
                                  (case1.comp$TSSUTR == "") & 
                                  (case1.comp$termUTR == "")),]
TssUTR.case1.comp <- case1.comp[((!case1.comp$TSSUTR == "") & 
                                   (case1.comp$Locus == "") & 
                                   (case1.comp$termUTR == "")),]
orphan.case1.comp <- case1.comp[((case1.comp$Locus == "") & 
                                   (case1.comp$TSSUTR == "") & 
                                   (case1.comp$termUTR == "")),]
termUTR.case1.comp <- case1.comp[((case1.comp$Locus == "") & 
                                    (case1.comp$TSSUTR == "") & 
                                    (!case1.comp$termUTR == "")),]
other1 <- case1.comp[((!case1.comp$Locus == "") & 
                        (!case1.comp$TSSUTR == "") & 
                        (!case1.comp$termUTR == "")),]
other2 <- case1.comp[((case1.comp$Locus == "") & 
                        (!case1.comp$TSSUTR == "") & 
                        (!case1.comp$termUTR == "")),]
other3 <- case1.comp[((!case1.comp$Locus == "") & 
                        (case1.comp$TSSUTR == "") & 
                        (!case1.comp$termUTR == "")),]
other4 <- case1.comp[((!case1.comp$Locus == "") &
                        (!case1.comp$TSSUTR == "") &
                        (case1.comp$termUTR == "")),]
other.case1.comp <- rbind(other1, other2, other3, other4)
#########################################################
##########  Case 2 - Transcriptome Distribution #########
#########################################################
############  Top
genic.case2.top <- case2.top[((!case2.top$Locus == "") & 
                                (case2.top$TSSUTR == "") & 
                                (case2.top$termUTR == "")),]
TssUTR.case2.top <- case2.top[((!case2.top$TSSUTR == "") & 
                                 (case2.top$Locus == "") & 
                                 (case2.top$termUTR == "")),]
orphan.case2.top <- case2.top[((case2.top$Locus == "") & 
                                 (case2.top$TSSUTR == "") & 
                                 (case2.top$termUTR == "")),]
termUTR.case2.top <- case2.top[((case2.top$Locus == "") & 
                                  (case2.top$TSSUTR == "") & 
                                  (!case2.top$termUTR == "")),]
other1 <- case2.top[((!case2.top$Locus == "") & 
                       (!case2.top$TSSUTR == "") & 
                       (!case2.top$termUTR == "")),]
other2 <- case2.top[((case2.top$Locus == "") & 
                       (!case2.top$TSSUTR == "") & 
                       (!case2.top$termUTR == "")),]
other3 <- case2.top[((!case2.top$Locus == "") & 
                       (case2.top$TSSUTR == "") &
                       (!case2.top$termUTR == "")),]
other4 <- case2.top[((!case2.top$Locus == "") & 
                       (!case2.top$TSSUTR == "") & 
                       (case2.top$termUTR == "")),]
other.case2.top <- rbind(other1, other2, other3, other4)
############  Comp
genic.case2.comp <- case2.comp[((!case2.comp$Locus == "") & 
                                  (case2.comp$TSSUTR == "") & 
                                  (case2.comp$termUTR == "")),]
TssUTR.case2.comp <- case2.comp[((!case2.comp$TSSUTR == "") & 
                                   (case2.comp$Locus == "") & 
                                   (case2.comp$termUTR == "")),]
orphan.case2.comp <- case2.comp[((case2.comp$Locus == "") &
                                   (case2.comp$TSSUTR == "") & 
                                   (case2.comp$termUTR == "")),]
termUTR.case2.comp <- case2.comp[((case2.comp$Locus == "") &
                                    (case2.comp$TSSUTR == "") & 
                                    (!case2.comp$termUTR == "")),]
other1 <- case2.comp[((!case2.comp$Locus == "") & 
                        (!case2.comp$TSSUTR == "") & 
                        (!case2.comp$termUTR == "")),]
other2 <- case2.comp[((case2.comp$Locus == "") & 
                        (!case2.comp$TSSUTR == "") & 
                        (!case2.comp$termUTR == "")),]
other3 <- case2.comp[((!case2.comp$Locus == "") &
                        (case2.comp$TSSUTR == "") & 
                        (!case2.comp$termUTR == "")),]
other4 <- case2.comp[((!case2.comp$Locus == "") & 
                        (!case2.comp$TSSUTR == "") & 
                        (case2.comp$termUTR == "")),]
other.case2.comp <- rbind(other1, other2, other3, other4)
#########################################################
##########  Case 3 - Transcriptome Distribution #########
#########################################################
############  Top
genic.case3.top <- case3.top[((!case3.top$Locus == "") & 
                                (case3.top$TSSUTR == "") & 
                                (case3.top$termUTR == "")),]
TssUTR.case3.top <- case3.top[((!case3.top$TSSUTR == "") &
                                 (case3.top$Locus == "") &
                                 (case3.top$termUTR == "")),]
orphan.case3.top <- case3.top[((case3.top$Locus == "") & 
                                 (case3.top$TSSUTR == "") & 
                                 (case3.top$termUTR == "")),]
termUTR.case3.top <- case3.top[((case3.top$Locus == "") &
                                  (case3.top$TSSUTR == "") & 
                                  (!case3.top$termUTR == "")),]
other1 <- case3.top[((!case3.top$Locus == "") & 
                       (!case3.top$TSSUTR == "") & 
                       (!case3.top$termUTR == "")),]
other2 <- case3.top[((case3.top$Locus == "") &
                       (!case3.top$TSSUTR == "") & 
                       (!case3.top$termUTR == "")),]
other3 <- case3.top[((!case3.top$Locus == "") & 
                       (case3.top$TSSUTR == "") &
                       (!case3.top$termUTR == "")),]
other4 <- case3.top[((!case3.top$Locus == "") & 
                       (!case3.top$TSSUTR == "") & 
                       (case3.top$termUTR == "")),]
other.case3.top <- rbind(other1, other2, other3, other4)
############  Comp
genic.case3.comp <- case3.comp[((!case3.comp$Locus == "") & 
                                  (case3.comp$TSSUTR == "") & 
                                  (case3.comp$termUTR == "")),]
TssUTR.case3.comp <- case3.comp[((!case3.comp$TSSUTR == "") & 
                                   (case3.comp$Locus == "") & 
                                   (case3.comp$termUTR == "")),]
orphan.case3.comp <- case3.comp[((case3.comp$Locus == "") &
                                   (case3.comp$TSSUTR == "") & 
                                   (case3.comp$termUTR == "")),]
termUTR.case3.comp <- case3.comp[((case3.comp$Locus == "") & 
                                    (case3.comp$TSSUTR == "") & 
                                    (!case3.comp$termUTR == "")),]
other1 <- case3.comp[((!case3.comp$Locus == "") & 
                        (!case3.comp$TSSUTR == "") & 
                        (!case3.comp$termUTR == "")),]
other2 <- case3.comp[((case3.comp$Locus == "") & 
                        (!case3.comp$TSSUTR == "") & 
                        (!case3.comp$termUTR == "")),]
other3 <- case3.comp[((!case3.comp$Locus == "") & 
                        (case3.comp$TSSUTR == "") &
                        (!case3.comp$termUTR == "")),]
other4 <- case3.comp[((!case3.comp$Locus == "") &
                        (!case3.comp$TSSUTR == "") & 
                        (case3.comp$termUTR == "")),]
other.case3.comp <- rbind(other1, other2, other3, other4)
# Define the cases
cases <- c("case1.top", "case2.top", "case3.top", 
           "case1.comp", "case2.comp", "case3.comp")
# Create an empty data frame
EnrichedOutput <- as.data.frame(matrix(0, nrow = length(cases), ncol = 13))
# Set column names
colnames(EnrichedOutput) <- c("cases", "totalnumber", "genic.number", "orphan.number", "UTR.number", "termUTR.number", "other","genicpercentage", "orphan.percentage", "UTR.percentage", "termUTR.percentage", "other.percentage", "total.percentage")
# Populate the 'cases' column
EnrichedOutput$cases <- cases
for(i in 1:length(cases)) {
  case <- cases[i]
  if (exists(case)) {
    EnrichedOutput$totalnumber[i] <- nrow(get(case))
    genic_var <- paste("genic.", case, sep = "")
    orphan_var <- paste("orphan.", case, sep = "")
    UTR_var <- paste("TssUTR.", case, sep = "")
    termUTR_var <- paste("termUTR.", case, sep = "")
    other_var <- paste("other.", case, sep = "")
    EnrichedOutput$genic.number[i] <- if(exists(genic_var)) nrow(get(genic_var)) else 0
    EnrichedOutput$orphan.number[i] <- if(exists(orphan_var)) nrow(get(orphan_var)) else 0
    EnrichedOutput$UTR.number[i] <- if(exists(UTR_var)) nrow(get(UTR_var)) else 0
    EnrichedOutput$termUTR.number[i] <- if(exists(termUTR_var)) nrow(get(termUTR_var)) else 0
    EnrichedOutput$other[i] <- if(exists(other_var)) nrow(get(other_var)) else 0
  } else {
    # Handle the case where the variable doesn't exist
    warning(paste("Variable", case, "not found. Setting values to 0."))
  }
}
# Calculate percentages
EnrichedOutput[, 8:13] <- 100 * EnrichedOutput[, 3:8] / EnrichedOutput$totalnumber
EnrichedOutput$total.percentage <- rowSums(EnrichedOutput[, 8:13])
EnrichedOutput[7:9,1] <- c("case1.combined", "case2.combined", "case3.combined")
for (i in 2:7) {
  for (j in 7:9) {
    EnrichedOutput[j, i] <- EnrichedOutput[j-6, i] + EnrichedOutput[j-3, i]
  }
}
for(i in 8:13){
  for(j in 7:9){
    EnrichedOutput[j,i] <- (EnrichedOutput[j-6, i] + EnrichedOutput[j-3,i])/2  
  }
}
write.csv(EnrichedOutput, file = "Termseq-EnrichedOutput.csv")
RNAfold <- rbind(top, comp)
RNAfold <- RNAfold %>% arrange(PeakCoord)
outputfilename <- args[3]
write_csv(RNAfold, outputfilename, col_names = T)

