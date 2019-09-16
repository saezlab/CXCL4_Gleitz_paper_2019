setwd("~/Dropbox/CXCL4_gleitz_paper/")

source("scripts/limmaWrapper.R")
library(readxl)
library(dplyr)
library(vsn)
library(limma)
library(readr)
library(ggplot2)

early_fibro <- as.data.frame(read_excel("data/early/combined (Rebekka Schneiders in Konflikt stehende Kopie 2017-01-18).xls"))

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

early_fibro <- early_fibro %>% group_by(symbole) %>% summarise_each(funs(sum(., na.rm = TRUE)))
early_fibro <- as.data.frame(early_fibro)

row.names(early_fibro) <- early_fibro[,1]
early_fibro <- early_fibro[,-1]

late_fibro <- as.data.frame(read_excel("data/late/expression.xls"))

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

late_fibro <- late_fibro %>% group_by(symbole) %>% summarise_each(funs(sum(., na.rm = TRUE)))
late_fibro <- as.data.frame(late_fibro)

row.names(late_fibro) <- late_fibro[,1]
late_fibro <- late_fibro[,-1]

names(late_fibro) <- gsub("FPKM","late_Stroma_",names(late_fibro))
names(late_fibro) <- gsub("control","EV",names(late_fibro))
names(early_fibro) <- gsub("FPKM","early_",names(early_fibro))


targets_late <- as.data.frame(matrix(NA,length(late_fibro[1,]),2))
names(targets_late) <- c("sample","condition")
targets_late$sample <- names(late_fibro)
targets_late$condition <- gsub("[0-9]","",targets_late$sample)

targets_early <- as.data.frame(matrix(NA,length(early_fibro[1,]),2))
names(targets_early) <- c("sample","condition")
targets_early$sample <- names(early_fibro)
targets_early$condition <- gsub("[0-9]","",targets_early$sample)

###########
###########

late_fibro[late_fibro < 0.001] <- NA
late_fibro <- late_fibro[complete.cases(late_fibro),]

fit <- vsnMatrix(as.matrix(late_fibro))
meanSdPlot(fit)
late_fibro <- as.data.frame(predict(fit,as.matrix(late_fibro)))

early_fibro[early_fibro < 0.001] <- NA
early_fibro <- early_fibro[complete.cases(early_fibro),]

fit <- vsnMatrix(as.matrix(early_fibro))
meanSdPlot(fit)
early_fibro <- as.data.frame(predict(fit,as.matrix(early_fibro)))

############
############

unique(targets_late$condition)
limma_late <- runLimma(late_fibro, targets_late, comparisons = list(c(2,-1)))

ttop_late_stroma_TPO_vs_stroma_EV <- ttopFormatter(topTable(limma_late[[1]], coef = 1, number = length(late_fibro[,1]), adjust.method = "fdr"))

unique(targets_early$condition)
limma_early <- runLimma(early_fibro, targets_early, comparisons = list(c(3,-1),c(4,-2)))

ttop_early_HSPC_TPO_vs_HSPC_EV <- ttopFormatter(topTable(limma_early[[1]], coef = 1, number = length(early_fibro[,1]), adjust.method = "fdr"))
ttop_early_stroma_TPO_vs_stroma_EV <- ttopFormatter(topTable(limma_early[[1]], coef = 2, number = length(early_fibro[,1]), adjust.method = "fdr"))
# 
write_csv(ttop_late_stroma_TPO_vs_stroma_EV, "results/ttop_late_stroma_TPO_vs_stroma_EV.csv")
write_csv(ttop_early_stroma_TPO_vs_stroma_EV, "results/ttop_early_stroma_TPO_vs_stroma_EV.csv")
write_csv(ttop_early_HSPC_TPO_vs_HSPC_EV, "results/ttop_early_HSPC_TPO_vs_HSPC_EV.csv")
write_csv(late_fibro,"data/late_fibro_vsv.csv")
write_csv(early_fibro, "data/early_fibro_vsn.csv")
write_csv(targets_late,"support/targets_late.csv")
write_csv(targets_early,"support/targets_early.csv")
#
ggsave("visualisation/late/vsn/late_stroma_TPO_vs_stroma_EV.pdf",volcano_nice(ttop_late_stroma_TPO_vs_stroma_EV, FCIndex = 2, IDIndex = 1, pValIndex = 5, label = T, nlabels = 30, manual_labels = c("PF4","S100A8","S100A9")))
ggsave("visualisation/early/vsn/early_stroma_TPO_vs_stroma_EV.pdf",volcano_nice(ttop_early_stroma_TPO_vs_stroma_EV, FCIndex = 2, IDIndex = 1, pValIndex = 5, label = T, nlabels = 30, manual_labels = c("PF4","S100A8","S100A9")))
ggsave("visualisation/early/vsn/early_HSPC_TPO_vs_HSPC_EV.pdf",volcano_nice(ttop_early_HSPC_TPO_vs_HSPC_EV, FCIndex = 2, IDIndex = 1, pValIndex = 5, label = T, nlabels = 30, manual_labels = c("PF4","S100A8","S100A9")))
#

       