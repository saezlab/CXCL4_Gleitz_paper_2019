#This is a script designed to perform PROGENY et DOROTHEA analysis.
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
library(pheatmap)
library(viper)
library(grid)
library(ggplot2)
library(ggrepel)
library(grid)
library(gridExtra)

source("scripts/heatmap_support.R")
source("scripts/support_functions.R")

##############################################################
#######               PROGENY ANALYSIS               #########
##############################################################

### LOAD THE DIFFERENTIAL ANALYSIS RESULTS
ttop_late_stroma_TPO_vs_stroma_EV <- as.data.frame(read_csv("results/ttop_late_stroma_TPO_vs_stroma_EV.csv"))
ttop_early_HSPC_TPO_vs_HSPC_EV <- as.data.frame(read_csv("results/ttop_early_HSPC_TPO_vs_HSPC_EV.csv"))
ttop_early_stroma_TPO_vs_stroma_EV <- as.data.frame(read_csv("results/ttop_early_stroma_TPO_vs_stroma_EV.csv"))

### LOAD THE PROGENY MODEL MATRIX
progeny_matrix_mouse_v1 <- as.data.frame(read_csv("support/progeny_matrix_mouse_v1.csv"))
progeny_matrix_mouse_v1$X1 <- toupper(progeny_matrix_mouse_v1$X1)

### FORMAT THE DIFFERENTIAL ANALYSIS RESUTS AS INPUTS FOR PROGENY
t_table <- merge(ttop_late_stroma_TPO_vs_stroma_EV[,c(1,4)], ttop_early_stroma_TPO_vs_stroma_EV[,c(1,4)], by = "ID")
t_table <- merge(t_table, ttop_early_HSPC_TPO_vs_HSPC_EV[,c(1,4)], by = "ID")

names(t_table) <- c("ID","late_stroma","early_stroma","early_HSPC")
t_table <- t_table[,c(1,4,3,2)]

### MAKE SCATTER PLOTS TO REPRESENT THE PATHWAY RESPONSIVE GENES WITH THEIR T VALUES AND WEIGHTS
scats <- progenyScatter(t_table, progeny_matrix_mouse_v1, statName = "t_value")
names(scats) <- names(t_table)[c(2,3,4)]

saveProgenyPlots(scats, names(scats), "visualisation/progeny_scatters/")

### RUN PROGENY
progeny_res <- runProgenyFast(df = t_table, weight_matrix = progeny_matrix_mouse_v1, k = 10000, z_score = T)

### FORMAT PROGENY RESULTS
progeny_res_df <- progeny_res
progeny_res_df$pathway <- row.names(progeny_res_df)
progeny_res_df <- progeny_res_df[,c(4,1,2,3)]
names(progeny_res_df) <- names(t_table)

### GENERATE P-VALUES FROM PROGENY Z-SCORES
progeny_pvals <- as.data.frame(apply(progeny_res, 2, pnorm))
progeny_pvals[] <- lapply(progeny_pvals, function(x){ifelse(x > 0.5, 1 - x, x)})
progeny_pvals <- progeny_pvals[c(6,2,9,7,11,5,10,4,8,1,13,14,3,12),]

## MAKE HEATMAP FOR PVALUE
t <- as.vector(t(progeny_pvals))
palette1 <- createLinearColors(t[t > 0],withZero = F , maximum = 99, my_colors =  c("red","royalblue3","white"))
# palette2 <- createLinearColors(t[t < 0],withZero = F , maximum = 1)
palette <- c(palette1)

pheatmap(progeny_pvals, display_numbers = T, cluster_cols = F, cluster_rows = F, color = palette)

## MAKE HEATMAP FOR PROGENY Z-SCORES
t <- as.vector(t(progeny_res_df[,-1]))
palette1 <- createLinearColors(t[t < 0],withZero = F , maximum = 51)
palette2 <- createLinearColors(t[t > 0],withZero = F , maximum = 49)
palette <- c(palette1,palette2)

pheatmap(progeny_res_df[,-1], display_numbers = T, cluster_cols = F, color = palette, cluster_rows = F)

##############################################################
#######              DOROTHEA ANALYSIS               #########
##############################################################

### IMPORT DOROTHEA REGULONS
dorothea_regulon_mouse_v1 <- as.data.frame(read_csv("support/dorothea_regulon_mouse_v1.csv"))

### FORMAT DOROTHEA REGULON FOR VIPER 
dorothea_regulon_mouse_v1$tf <- toupper(dorothea_regulon_mouse_v1$tf)
dorothea_regulon_mouse_v1$target <- toupper(dorothea_regulon_mouse_v1$target)
dorothea_regulon_mouse_v1 <- dorothea_regulon_mouse_v1[dorothea_regulon_mouse_v1$confidence %in% c("A","B","C"),]
dorothea_regulon_mouse_v1 <- dorothea_regulon_mouse_v1[,c(3,1,4)]

dorothea_regulon_viper <- df_to_viper_regulon(dorothea_regulon_mouse_v1)

row.names(t_table) <- t_table$ID

### RUN VIPER WITH DOROTHEA REGULONS
viperRes <- as.data.frame(t(viper(eset = t_table[,-1], regulon = dorothea_regulon_viper, minsize = 25, adaptive.size = F, eset.filter = F, pleiotropy = T)))

### FORMAT RESULTS AND CREATE HEATMAP
k <- c()
for(i in 1:length(viperRes[1,]))
{
  if(max(abs(viperRes[,i])) > 4)
  {
    k <- c(k,i)
  }
}

viperRes <- viperRes[,k]

t <- as.vector(viperRes)
palette1 <- createLinearColors(t[t < 0],withZero = F , maximum = 64)
palette2 <- createLinearColors(t[t > 0],withZero = F , maximum = 36)
palette <- c(palette1,palette2)

pheatmap(viperRes, cluster_cols = F, display_numbers = T, color = palette)

### WRITE RESULTS
write_csv(progeny_res_df, "results/progeny_early_late.csv")
viperRes$condition <- names(t_table[,-1])
viperRes <- viperRes[,c(length(viperRes),1:(length(viperRes[,-1])))]
write_csv(viperRes,"results/TF_activities.csv")

