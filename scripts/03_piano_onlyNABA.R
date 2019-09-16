setwd("~/Dropbox/CXCL4_gleitz_paper/")

source("scripts/heatmap_support.R")

library(readr)

source("scripts/support_functions.R")
source("scripts/limmaWrapper.R")


ttop_late_stroma_TPO_vs_stroma_EV <- as.data.frame(read_csv("results/ttop_late_stroma_TPO_vs_stroma_EV.csv"))
ttop_early_HSPC_TPO_vs_HSPC_EV <- as.data.frame(read_csv("results/ttop_early_HSPC_TPO_vs_HSPC_EV.csv"))
ttop_early_stroma_TPO_vs_stroma_EV <- as.data.frame(read_csv("results/ttop_early_stroma_TPO_vs_stroma_EV.csv"))

NABA <- gmt_to_csv("support/NABA_matrisone.gmt")
NABA_dummy <- NABA
NABA_dummy$term <- "dummy"
NABA <- as.data.frame(rbind(NABA,NABA_dummy))

NABA_late_stroma_TPO_vs_stroma_EV <- runPIANO(ttop_late_stroma_TPO_vs_stroma_EV, NABA)
NABA_late_stroma_TPO_vs_stroma_EV <- as.data.frame(NABA_late_stroma_TPO_vs_stroma_EV[[1]])

NABA_early_HSPC_TPO_vs_HSPC_EV <- runPIANO(ttop_early_HSPC_TPO_vs_HSPC_EV, NABA)
NABA_early_HSPC_TPO_vs_HSPC_EV <- as.data.frame(NABA_early_HSPC_TPO_vs_HSPC_EV[[1]])

NABA_early_stroma_TPO_vs_stroma_EV <- runPIANO(ttop_early_stroma_TPO_vs_stroma_EV, NABA)
NABA_early_stroma_TPO_vs_stroma_EV <- as.data.frame(NABA_early_stroma_TPO_vs_stroma_EV[[1]])

library(ggplot2)
volcano_nice(ttop_early_HSPC_TPO_vs_HSPC_EV[ttop_early_HSPC_TPO_vs_HSPC_EV$ID %in% NABA$gene,], FCIndex = 2, pValIndex = 5, IDIndex = 1, label = T, nlabels = 30, manual_labels = c("PF4", "S100A8", "S100A9")) + ylab("-log(p-value)")
volcano_nice(ttop_early_stroma_TPO_vs_stroma_EV[ttop_early_stroma_TPO_vs_stroma_EV$ID %in% NABA$gene,], FCIndex = 2, pValIndex = 5, IDIndex = 1, label = T, nlabels = 30, manual_labels = c("PF4", "S100A8", "S100A9")) + ylab("-log(p-value)")
volcano_nice(ttop_late_stroma_TPO_vs_stroma_EV[ttop_late_stroma_TPO_vs_stroma_EV$ID %in% NABA$gene,], FCIndex = 2, pValIndex = 5, IDIndex = 1, label = T, nlabels = 30, manual_labels = c("PF4", "S100A8", "S100A9")) + ylab("-log(p-value)")
