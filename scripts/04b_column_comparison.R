#This is a script designed to analysis an RNA counts dataset.
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
library(vsn)
library(limma)
library(viper)
library(pheatmap)
library(ggplot2)
library(grid)
library(gridExtra)

source("scripts/support_functions.R")
source("scripts/heatmap_support.R")
source("scripts/limmaWrapper.R")

##############################################################
#######           DATA LOADING AND CLEANING          #########
##############################################################

### CREATE DATAFRAMES TO SUMMARISE THE EXPERIMENTAL DESIGN
full_targets <- as.data.frame(read_delim("support/full_targets.csv", 
                                         ";", escape_double = FALSE, trim_ws = TRUE))

full_targets$`Mouse group` <- gsub(" ","_",full_targets$`Mouse group`)
full_targets$cellTypeorMarker <- gsub("CD41","MegaK",full_targets$cellTypeorMarker)
full_targets$cellTypeorMarker <- gsub("TdTomato_Gli1","Stromal",full_targets$cellTypeorMarker)
full_targets$`Mouse group` <- paste(full_targets$cellTypeorMarker,full_targets$`Mouse group`, sep = "_")

full_targets <- full_targets[order(full_targets$`Mouse group`),]
full_targets$replicates <- rep(c(1,2,3),8)
full_targets$sample <- paste(full_targets$`Mouse group`,full_targets$replicates, sep = "_")

targets <- full_targets[,c(12,10)]
names(targets) <- c("sample","condition")

### IMPORT THE COUNT FILES AND COMBINE THEM IN A SINGLE DATAFRAME
count_file_list <- list.files("data/count_files/", full.names = T)
count_file_names <- gsub("[.]counts","",list.files("data/count_files/", full.names = F))

count_df_list <- list()
for(count_file in count_file_list)
{
  count_df <- as.data.frame(read_delim(count_file, 
                                       "\t", escape_double = FALSE, col_names = FALSE, 
                                       trim_ws = TRUE))
  count_df_list[[gsub("[.]counts","",count_file)]] <- count_df
}

i <- 1
for(count_df in count_df_list)
{
  if(i == 1)
  {
    batches <- count_df
  }
  else
  {
    batches <- merge(batches,count_df, by = "X1", all = T)
  }
  i <- i+1
}

names(batches) <- c("ID", count_file_names)

batches <- batches[-c(1:5),]

batches <- batches[,c("ID",full_targets$`Sample Name`)]
names(batches) <- c("ID",targets$sample)

row.names(batches) <- batches$ID
batches <- batches[,-1]

### REMOVE LOW COUNTS
SDs <- apply(batches,1,sd)
means <- apply(batches,1,mean)

batches <- batches[means > 50,]

batches[batches == 0] <- 0.1

### NORMALISATION WITH VSN

fit <- vsnMatrix(as.matrix(batches))
meanSdPlot(fit)
batches <- as.data.frame(predict(fit,as.matrix(batches)))


##############################################################
#######             DIFFERENTIAL ANALYSIS            #########
##############################################################

### DEFINE THE CONTRASTS OF THE DIFFERENTIAL ANALYSIS
comparisons_2 <- list("MegaK_EV" = c(1,-3), "MegaK_TPO" = c(2,-4), "Stromal_EV" = c(5,-7), "Stromal_TPO" = c(6,-8))

### RUN LIMMA FOR DIFFERENTIAL ANALYSIS
limmaRes_2 <- runLimma(batches,targets,comparisons_2)

### FORMAT LIMMA RESULTS
ttop_list <- list()
for(i in 1:length(comparisons_2))
{
  ttop_list[[i]] <- ttopFormatter(topTable(limmaRes_2[[1]], coef = i, number = length(batches[,1]), adjust.method = "fdr"))
}
names(ttop_list) <- names(comparisons_2)
View(ttop_list[[2]])
View(ttop_list[[4]])
View(ttop_list[[1]])
View(ttop_list[[3]])

library(QQperm)

### CHECK THE PVALUE DISTRIBUTION WITH QQPLOTS
expected <- pnorm(rnorm(length(ttop_list[[1]][,1])))

qqplot(P.perm = sort(expected), P.observed = sort(ttop_list[[2]][,5]), xlim = c(0,6), ylim = c(0,6))
qqplot(P.perm = sort(expected), P.observed = sort(ttop_list[[4]][,5]), xlim = c(0,6), ylim = c(0,6))

### FORMAT THE DIFFERENTIAL ANALYSIS RESULTS
for(i in 1:length(comparisons_2))
{
  if( i == 1)
  {
    t_table <- ttop_list[[i]][,c(1,4)]
  }
  else
  {
    t_table <- merge(t_table,ttop_list[[i]][,c(1,4)], by = "ID")
  }
}
names(t_table) <- c("ID",names(comparisons_2))

##############################################################
#######               PROGENY ANALYSIS               #########
##############################################################

### LOAD THE PROGENY MODEL MATRIX
model_matrix <- as.data.frame(read_csv("support/progeny_matrix_mouse_v1.csv"))

### RUN PROGENY
progeny_res <- runProgenyFast(t_table, model_matrix)

### FORMAT PROGENY RESULTS AND MAKE HEATMAP
progeny_df <- progeny_res

t <- as.vector(t(progeny_df))
palette1 <- createLinearColors(t[t < 0],withZero = F , maximum = 39)
palette2 <- createLinearColors(t[t > 0],withZero = F , maximum = 61)
palette <- c(palette1,palette2)

pheatmap(progeny_df, display_numbers = T, cluster_cols = F, color = palette)

### COMPUTE PVALUE FOR PROGENY Z-SCORES AND MAKE HEATMAP
progeny_pvals <- as.data.frame(apply(progeny_df, 2, pnorm))
progeny_pvals[] <- lapply(progeny_pvals, function(x){ifelse(x > 0.5, 1 - x, x)})

t <- as.vector(t(progeny_pvals))
palette1 <- createLinearColors(t[t > 0],withZero = F , maximum = 99, my_colors =  c("red","royalblue3","white"))
# palette2 <- createLinearColors(t[t < 0],withZero = F , maximum = 1)
palette <- c(palette1)

pheatmap(progeny_pvals, display_numbers = T, cluster_cols = F, cluster_rows = F, color = palette)

##############################################################
#######              DOROTHEA ANALYSIS               #########
##############################################################

source("scripts/viperWrapper.R")

### IMPORT DOROTHEA REGULONS
dorothea_regulon_mouse_v1 <- as.data.frame(read_csv("support/dorothea_regulon_mouse_v1.csv"))

dorothea_regulon_mouse_viper_AB <- df_to_viper_regulon(dorothea_regulon_mouse_v1[dorothea_regulon_mouse_v1$confidence %in% c("A","B"),c(3,1,4)])

### RUN VIPER
TF_activity_viper <- runViper(ttop_list, regulon = dorothea_regulon_mouse_viper_AB)

TF_activity_viper_df <- makeViperResDf(viperResList = TF_activity_viper)
row.names(TF_activity_viper_df) <- TF_activity_viper_df$ID
TF_activity_viper_df <- TF_activity_viper_df[,-1]

###MAKE HEATMAP OF TF ACTIVITIES
t <- as.vector(t(TF_activity_viper_df))
palette1 <- createLinearColors(t[t < 0],withZero = F , maximum = 60)
palette2 <- createLinearColors(t[t > 0],withZero = F , maximum = 40)
palette <- c(palette1,palette2)

pheatmap(t(TF_activity_viper_df[rowSums(abs(TF_activity_viper_df)) > 5,]), cluster_cols = T, display_numbers = T, color = palette, height = 3, width = 12)

