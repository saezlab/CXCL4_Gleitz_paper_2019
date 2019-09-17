#This is a script designed to perform Enrichment analysis with PIANO.
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
library(readr)

source("scripts/support_functions.R")
source("scripts/limmaWrapper.R")
source("scripts/heatmap_support.R")

### IMPORT DIFFERENTIAL ANALYSIS RESULTS
ttop_late_stroma_TPO_vs_stroma_EV <- as.data.frame(read_csv("results/ttop_late_stroma_TPO_vs_stroma_EV.csv"))
ttop_early_HSPC_TPO_vs_HSPC_EV <- as.data.frame(read_csv("results/ttop_early_HSPC_TPO_vs_HSPC_EV.csv"))
ttop_early_stroma_TPO_vs_stroma_EV <- as.data.frame(read_csv("results/ttop_early_stroma_TPO_vs_stroma_EV.csv"))

### IMPORT NABA GENESET
NABA <- gmt_to_csv("support/NABA_matrisone.gmt")
NABA_dummy <- NABA
NABA_dummy$term <- "dummy" #There are technical issues when running piano with only one pathway, so we just duplicate the pathway
NABA <- as.data.frame(rbind(NABA,NABA_dummy))

### RUN PIANO
NABA_late_stroma_TPO_vs_stroma_EV <- runPIANO(ttop_late_stroma_TPO_vs_stroma_EV, NABA)
NABA_late_stroma_TPO_vs_stroma_EV <- as.data.frame(NABA_late_stroma_TPO_vs_stroma_EV[[1]])

NABA_early_HSPC_TPO_vs_HSPC_EV <- runPIANO(ttop_early_HSPC_TPO_vs_HSPC_EV, NABA)
NABA_early_HSPC_TPO_vs_HSPC_EV <- as.data.frame(NABA_early_HSPC_TPO_vs_HSPC_EV[[1]])

NABA_early_stroma_TPO_vs_stroma_EV <- runPIANO(ttop_early_stroma_TPO_vs_stroma_EV, NABA)
NABA_early_stroma_TPO_vs_stroma_EV <- as.data.frame(NABA_early_stroma_TPO_vs_stroma_EV[[1]])

### MAKE VOLCANO PLOT FOR NABA SPECIFIC GENES
library(ggplot2)
volcano_nice(ttop_early_HSPC_TPO_vs_HSPC_EV[ttop_early_HSPC_TPO_vs_HSPC_EV$ID %in% NABA$gene,], FCIndex = 2, pValIndex = 5, IDIndex = 1, label = T, nlabels = 30, manual_labels = c("PF4", "S100A8", "S100A9")) + ylab("-log(p-value)")
volcano_nice(ttop_early_stroma_TPO_vs_stroma_EV[ttop_early_stroma_TPO_vs_stroma_EV$ID %in% NABA$gene,], FCIndex = 2, pValIndex = 5, IDIndex = 1, label = T, nlabels = 30, manual_labels = c("PF4", "S100A8", "S100A9")) + ylab("-log(p-value)")
volcano_nice(ttop_late_stroma_TPO_vs_stroma_EV[ttop_late_stroma_TPO_vs_stroma_EV$ID %in% NABA$gene,], FCIndex = 2, pValIndex = 5, IDIndex = 1, label = T, nlabels = 30, manual_labels = c("PF4", "S100A8", "S100A9")) + ylab("-log(p-value)")
