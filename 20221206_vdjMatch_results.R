library(tidyverse)
library(dplyr)
library(reshape2)

setwd("/Volumes/PRIVATE/Individual folders")
setwd("/Volumes/PRIVATE/Individual folders/Vladimir/R_analysis/VDJmatch/TRB_scRNAseq/")


temp = list.files(pattern="*.txt")
myfiles = lapply(temp, read.delim)
names <- gsub(".txt", "", temp)
names <- gsub("20221205_", "", names)

names(myfiles) <- names
# C1
C1 <- myfiles$C1_TRB[!is.na(myfiles$C1_TRB$v), ]
C1 <- C1[!duplicated(C1$cdr3aa), ]
C1 <- C1[!C1$vdjdb.score == 0, ]
C1 <- as.data.frame(table(C1$antigen.species))
colnames(C1)[2] <- "C1"

# C2
C2 <- myfiles$C2_TRB[!is.na(myfiles$C2_TRB$v), ]
C2 <- C2[!duplicated(C2$cdr3aa), ]
C2 <- C2[!C2$vdjdb.score == 0, ]
C2 <- as.data.frame(table(C2$antigen.species))
colnames(C2)[2] <- "C2"

# C3
C3 <- myfiles$C3_TRB[!is.na(myfiles$C3_TRB$v), ]
C3 <- C3[!duplicated(C3$cdr3aa), ]
C3 <- C3[!C3$vdjdb.score == 0, ]
C3 <- as.data.frame(table(C3$antigen.species))
colnames(C3)[2] <- "C3"

# C6
C6 <- myfiles$C6_TRB[!is.na(myfiles$C6_TRB$v), ]
C6 <- C6[!duplicated(C6$cdr3aa), ]
C6 <- C6[!C6$vdjdb.score == 0, ]
C6 <- as.data.frame(table(C6$antigen.species))
colnames(C6)[2] <- "C6"

# C7
C7 <- myfiles$C7_TRB[!is.na(myfiles$C7_TRB$v), ]
C7 <- C7[!duplicated(C7$cdr3aa), ]
C7 <- C7[!C7$vdjdb.score == 0, ]
C7 <- as.data.frame(table(C7$antigen.species))
colnames(C7)[2] <- "C7"

# C14
C14 <- myfiles$C14_TRB[!is.na(myfiles$C14_TRB$v), ]
C14 <- C14[!duplicated(C14$cdr3aa), ]
C14 <- C14[!C14$vdjdb.score == 0, ]
C14 <- as.data.frame(table(C14$antigen.species))
colnames(C14)[2] <- "C14"

# C16
C16 <- myfiles$C16_TRB[!is.na(myfiles$C16_TRB$v), ]
C16 <- C16[!duplicated(C16$cdr3aa), ]
C16 <- C16[!C16$vdjdb.score == 0, ]
C16 <- as.data.frame(table(C16$antigen.species))
colnames(C16)[2] <- "C16"

TCRB_sum <- merge(C1, C2, by = "Var1", all = TRUE)
TCRB_sum <- merge(TCRB_sum, C6, by = "Var1", all = TRUE)
TCRB_sum <- merge(TCRB_sum, C7, by = "Var1", all = TRUE)
TCRB_sum$C3 <- c(NA)
TCRB_sum$C14 <- c(NA)
TCRB_sum$C16 <- c(NA)
rownames(TCRB_sum) <- TCRB_sum$Var1
#TCRB_sum <- TCRB_sum[-1]
TCRB_sum[is.na(TCRB_sum)] <- 0
TCRB_sum_melt <- melt(TCRB_sum)
TCRB_sum_melt$variable <- factor(TCRB_sum_melt$variable, levels =  c("C1", "C2","C3", "C6", "C7", "C14", "C16"))

ggplot(TCRB_sum_melt, aes(fill=Var1, y=value, x=variable)) + 
  geom_bar(position="stack", stat="identity") +
  ylab("Number of TCR clonotypes perfectly
       matched with published TCRs") +
  xlab(NULL) +
  scale_fill_manual(values = gray.colors(4)) +
  theme_bw()
