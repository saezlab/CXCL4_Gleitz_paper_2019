#This is a script designed to clean and perform differential analysis on a transcriptomic dataset.
#Copyright (C) 2019  Aurelien Dugourd

#This program is free software: you can redistribute it and/or modify
#it under the terms of the GNU General Public License as published by
#the Free Software Foundation, either version 3 of the License, or
#(at your option) any later version.

#This program is distributed in the hope that it will be useful,
#but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#GNU General Public License for more details.

#You should have received a copy of the GNU General Public License
#along with this program.  If not, see <http://www.gnu.org/licenses/>.


#Set the working directory
setwd("~/Dropbox/CXCL4_gleitz_paper/") #There you should put the path to the main directory of the analysis "CXCL4_gleitz_paper"

#Package loading
source("scripts/limmaWrapper.R")
library(readxl)
library(dplyr)
library(vsn)
library(limma)
library(readr)
library(ggplot2)

##############################################################
#######           DATA LOADING AND CLEANING          #########
##############################################################


### LOAD THE EARLY FIBROSIS DATASET
early_fibro <- as.data.frame(read_excel("data/early/combined (Rebekka Schneiders in Konflikt stehende Kopie 2017-01-18).xls"))

### CLEANING
names(early_fibro) <- paste(early_fibro[1,], names(early_fibro),sep = "")

early_fibro <- early_fibro[,c(1,grep("[.][.]",names(early_fibro), invert = T))]

early_fibro <- early_fibro[-1,]
early_fibro <- early_fibro[,-c(14:17)]

for(i in 2:length(early_fibro[1,]))
{
  early_fibro[,i] <- as.numeric(early_fibro[,i])
}

names(early_fibro)[1] <- "symbole"
early_fibro$symbole <- toupper(early_fibro$symbole)

## Sum duplicated transcripts
early_fibro <- early_fibro %>% group_by(symbole) %>% summarise_each(funs(sum(., na.rm = TRUE)))
early_fibro <- as.data.frame(early_fibro)

row.names(early_fibro) <- early_fibro[,1]
early_fibro <- early_fibro[,-1]

### LOAD THE LATE FIBROSIS DATASET
late_fibro <- as.data.frame(read_excel("data/late/expression.xls"))

### CLEANING
names(late_fibro) <- paste(late_fibro[1,], names(late_fibro),sep = "")

late_fibro <- late_fibro[,c(1,grep("[.][.]",names(late_fibro), invert = T))]

late_fibro <- late_fibro[-1,]
late_fibro <- late_fibro[,-11]

for(i in 2:length(late_fibro[1,]))
{
  late_fibro[,i] <- as.numeric(late_fibro[,i])
}

names(late_fibro)[1] <- "symbole"
late_fibro$symbole <- toupper(late_fibro$symbole)

## Sum duplicated transcripts
late_fibro <- late_fibro %>% group_by(symbole) %>% summarise_each(funs(sum(., na.rm = TRUE)))
late_fibro <- as.data.frame(late_fibro)

row.names(late_fibro) <- late_fibro[,1]
late_fibro <- late_fibro[,-1]

names(late_fibro) <- gsub("FPKM","late_Stroma_",names(late_fibro))
names(late_fibro) <- gsub("control","EV",names(late_fibro))
names(early_fibro) <- gsub("FPKM","early_",names(early_fibro))

### CREATE DATAFRAMES TO SUMMARISE THE EXPERIMENTAL DESIGN
targets_late <- as.data.frame(matrix(NA,length(late_fibro[1,]),2))
names(targets_late) <- c("sample","condition")
targets_late$sample <- names(late_fibro)
targets_late$condition <- gsub("[0-9]","",targets_late$sample)

targets_early <- as.data.frame(matrix(NA,length(early_fibro[1,]),2))
names(targets_early) <- c("sample","condition")
targets_early$sample <- names(early_fibro)
targets_early$condition <- gsub("[0-9]","",targets_early$sample)

### REMOVE LOW EXPRESSION GENES FOR LATE FIBROSIS
late_fibro[late_fibro < 0.001] <- NA
late_fibro <- late_fibro[complete.cases(late_fibro),]

### NORMALISATION WITH VSN
fit <- vsnMatrix(as.matrix(late_fibro))
meanSdPlot(fit)
late_fibro <- as.data.frame(predict(fit,as.matrix(late_fibro)))

### REMOVE LOW EXPRESSION GENES FOR EARLY FIBROSIS
early_fibro[early_fibro < 0.001] <- NA
early_fibro <- early_fibro[complete.cases(early_fibro),]

### NORMALISATION WITH VSN
fit <- vsnMatrix(as.matrix(early_fibro))
meanSdPlot(fit)
early_fibro <- as.data.frame(predict(fit,as.matrix(early_fibro)))

##############################################################
#######             DIFFERENTIAL ANALYSIS            #########
##############################################################

### RUN LIMMA
unique(targets_late$condition)
limma_late <- runLimma(late_fibro, targets_late, comparisons = list(c(2,-1)))

### FORMAT LIMMA OUTPUT
ttop_late_stroma_TPO_vs_stroma_EV <- ttopFormatter(topTable(limma_late[[1]], coef = 1, number = length(late_fibro[,1]), adjust.method = "fdr"))

unique(targets_early$condition)
limma_early <- runLimma(early_fibro, targets_early, comparisons = list(c(3,-1),c(4,-2)))

ttop_early_HSPC_TPO_vs_HSPC_EV <- ttopFormatter(topTable(limma_early[[1]], coef = 1, number = length(early_fibro[,1]), adjust.method = "fdr"))
ttop_early_stroma_TPO_vs_stroma_EV <- ttopFormatter(topTable(limma_early[[1]], coef = 2, number = length(early_fibro[,1]), adjust.method = "fdr"))

### WRITE RESULTS
write_csv(ttop_late_stroma_TPO_vs_stroma_EV, "results/ttop_late_stroma_TPO_vs_stroma_EV.csv")
write_csv(ttop_early_stroma_TPO_vs_stroma_EV, "results/ttop_early_stroma_TPO_vs_stroma_EV.csv")
write_csv(ttop_early_HSPC_TPO_vs_HSPC_EV, "results/ttop_early_HSPC_TPO_vs_HSPC_EV.csv")
write_csv(late_fibro,"data/late_fibro_vsv.csv")
write_csv(early_fibro, "data/early_fibro_vsn.csv")
write_csv(targets_late,"support/targets_late.csv")
write_csv(targets_early,"support/targets_early.csv")

### MAKE SOME FIGURES
ggsave("visualisation/late/vsn/late_stroma_TPO_vs_stroma_EV.pdf",volcano_nice(ttop_late_stroma_TPO_vs_stroma_EV, FCIndex = 2, IDIndex = 1, pValIndex = 5, label = T, nlabels = 30, manual_labels = c("PF4","S100A8","S100A9")))
ggsave("visualisation/early/vsn/early_stroma_TPO_vs_stroma_EV.pdf",volcano_nice(ttop_early_stroma_TPO_vs_stroma_EV, FCIndex = 2, IDIndex = 1, pValIndex = 5, label = T, nlabels = 30, manual_labels = c("PF4","S100A8","S100A9")))
ggsave("visualisation/early/vsn/early_HSPC_TPO_vs_HSPC_EV.pdf",volcano_nice(ttop_early_HSPC_TPO_vs_HSPC_EV, FCIndex = 2, IDIndex = 1, pValIndex = 5, label = T, nlabels = 30, manual_labels = c("PF4","S100A8","S100A9")))
#

       