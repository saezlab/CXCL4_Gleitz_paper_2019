library(ggplot2)

progenyScatter <- function(df,weight_matrix,dfID = 1, weightID = 1, statName = "gene stats")
{
  plot_list_contrasts <- list(0)
  for (i in 2:length(df[1,]))
  {
    plot_list_pathways <- list(0)
    for (j in 2:length(weight_matrix[1,]))
    {
      sub_df <- df[,c(dfID,i)]
      
      
      
      pathway_weights <- weight_matrix[,c(weightID,j)]
      names(sub_df) <- c("ID","stat")
      
      minstat <- min(sub_df$stat)
      maxstat <- max(sub_df$stat)
      histo <- ggplot(sub_df, aes(x = stat, fill = "blue")) + geom_density() + coord_flip() + scale_fill_manual( values = c("#00c5ff")) + xlim(minstat, maxstat) + theme_minimal() + theme(legend.position = "none", axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank())
      
      names(pathway_weights) <- c("ID","weight")
      pathway_weights <- pathway_weights[pathway_weights$weight != 0,]
      
      percentile <- ecdf(sub_df$stat)
      
      sub_df <- merge(sub_df,pathway_weights,by = "ID")
      
      sub_df$color <- "3"
      sub_df[(sub_df$weight > 0 & sub_df$stat > 0),"color"] <- "1"
      sub_df[(sub_df$weight > 0 & sub_df$stat < 0),"color"] <- "2"
      sub_df[(sub_df$weight < 0 & sub_df$stat > 0),"color"] <- "2"
      sub_df[(sub_df$weight < 0 & sub_df$stat < 0),"color"] <- "1"
      
      sub_df[(percentile(sub_df$stat) < .95 & percentile(sub_df$stat) > .05),1] <- NA
      
      print(paste("weights of ",names(weight_matrix)[j], sep = ""))
      
      title <- paste("weights of ",names(weight_matrix)[j], sep = "")
      
      scatterplot <- ggplot(sub_df, aes(x = weight, y = stat, color = color)) + geom_point() +
        # scale_colour_manual(values = c("#15ff00","#ff0000","#c9c9c9")) + #green and red
        scale_colour_manual(values = c("red","royalblue3","grey")) +
        geom_label_repel(aes(label = ID)) +
        ylim(minstat, maxstat) + theme_minimal() + theme(legend.position = "none") + geom_vline(xintercept = 0, linetype = 'dotted') + geom_hline(yintercept = 0, linetype = 'dotted') + labs(x = title, y = statName)
      
      lay <- t(as.matrix(c(1,1,1,1,2)))
      gg <- arrangeGrob(scatterplot, histo, nrow = 1, ncol = 2, layout_matrix = lay)
      
      #grid.arrange(gg)
      plot_list_pathways[[j-1]] <- gg
    }
    names(plot_list_pathways) <- names(weight_matrix[,-weightID])
    plot_list_contrasts[[i-1]] <- plot_list_pathways
  }
  return(plot_list_contrasts)
}

makeViperResDf <- function(viperResList)
{
  i <- 1
  for (viperRes in viperResList)
  {
    viperScores <- viperRes
    
    viperScores <- data.frame(viperScores)
    print(viperScores)
    viperScores$ID <- row.names(viperScores)
    
    names(viperScores)[1] <- i
    
    if (i == 1)
    {
      viperResDf <- viperScores
    }
    else
    {
      viperResDf <- merge(viperResDf, viperScores, by = "ID", all = T)
    }
    
    i <- i+1
  }
  if(i > 2)
  {
    names(viperResDf)[2:length(viperResDf[1,])] <- names(viperResList)
  }
  else
  {
    viperResDf <- viperResDf[,c(2,1)]
    names(viperResDf) <- c("ID", names(viperResList))
  }
  
  return(viperResDf)
}

heatProgeny <- function(t_table, model_mat, pathway, out_dir, height = 5, width = 5)
{
  row.names(t_table) <- t_table[,1]
  model_mat <- model_mat[order(model_mat[,pathway], decreasing = T),]
  model_mat <- model_mat[model_mat[,1] %in% row.names(t_table),]
  
  responsive_gene_up <- model_mat[model_mat[,pathway] > 0,1]
  if(length(responsive_gene_up) > 0)
  {
    responsive_gene_up <- responsive_gene_up[1:10]
  }
  
  model_mat <- model_mat[order(model_mat[,pathway], decreasing = F),]
  responsive_gene_down <- model_mat[model_mat[,pathway] < 0,1]
  if(length(responsive_gene_down) > 0)
  {
    responsive_gene_down <- responsive_gene_down[1:10]
  }
  
  t_table_up <- t_table[t_table[,1] %in% responsive_gene_up,-1]
  t_table_down <- t_table[t_table[,1] %in% responsive_gene_down,-1]
  t <- as.vector(t(t_table_up))
  
  if(min(t) < 0 & max(t) > 0)
  {
    row.names(t_table_up) <- toupper(row.names(t_table_up))
    spread <- max(t) - min(t)
    lower_part <- abs(round((min(t)/spread) * 100))
    upper_part <- round((1 - abs((min(t)/spread))) * 100)
    palette1 <- createLinearColors(t[t < 0],withZero = F , maximum = lower_part)
    palette2 <- createLinearColors(t[t > 0],withZero = F , maximum = upper_part)
    palette <- c(palette1,palette2)
    file_name_up <- paste(out_dir, paste("progeny_", paste(pathway, "_up.pdf", sep= ""), sep = ""), sep = "/")
    pheatmap(t_table_up, display_numbers = T, cluster_cols = F, color = palette, filename = file_name_up, height = height, width = width, cluster_rows = F)
  }
  
  t <- as.vector(t(t_table_down))
  
  if(min(t) < 0 & max(t) > 0)
  {
    row.names(t_table_down) <- toupper(row.names(t_table_down))
    spread <- max(t) - min(t)
    lower_part <- abs(round((min(t)/spread) * 100))
    upper_part <- round((1 - abs((min(t)/spread))) * 100)
    palette1 <- createLinearColors(t[t < 0],withZero = F , maximum = lower_part)
    palette2 <- createLinearColors(t[t > 0],withZero = F , maximum = upper_part)
    palette <- c(palette1,palette2)
    file_name_down <- paste(out_dir, paste("progeny_", paste(pathway, "_down.pdf", sep= ""), sep = ""), sep = "/")
    pheatmap(t_table_down,display_numbers = T, cluster_cols = F, color = palette, filename = file_name_down, height = height, width = width, cluster_rows = F)
  }
}

createLinearColors <- function(numbers, withZero = T, maximum = 100, my_colors = c("royalblue3","white","red"))
{
  if (min(numbers, na.rm = T) > 0)
  {
    if(withZero)
    {
      numbers <- c(0,numbers)
    }
    myPalette <- colorRampPalette(my_colors[c(2,3)])
    myColors <- myPalette(maximum)
  }
  else
  {
    if (max(numbers, na.rm = T) < 0)
    {
      if(withZero)
      {
        numbers <- c(0,numbers)
      }
      myPalette <- colorRampPalette(my_colors[c(1,2)])
      myColors <- myPalette(maximum)
    }
    else
    {
      myPalette_pos <- colorRampPalette(c("white","red"))
      myPalette_neg <- colorRampPalette(c("royalblue3","white"))
      npos <- length(numbers[numbers >= 0]) + 1
      nneg <- length(numbers[numbers <= 0]) + 1
      
      myColors_pos <- myPalette_pos(npos)
      myColors_neg <- myPalette_neg(nneg)
      
      #print(myColors_neg)
      #print(myColors_pos)
      
      myColors <- c(myColors_neg[-(nneg)], myColors_pos[-1])
    }
  }
  return(myColors)
}

saveProgenyPlots <- function(plots, contrast_names, dirpath)
{
  i <- 1
  for (condition in plots)
  {
    dirname <- paste(dirpath,contrast_names[i], sep = "")
    dir.create(dirname, recursive = T, showWarnings = F)
    j <- 1
    for (pathway in condition)
    {
      filename <- paste(dirname,names(condition)[j],sep = "/")
      filename <- paste(filename,".pdf",sep = "")
      print(filename)
      ggsave(filename, pathway,device = "pdf", dpi = 300)
      j <- j+1
    }
    i <- i+1
  }
}

runProgenyFast <- function(df,weight_matrix,k = 10000, z_scores = T, get_nulldist = F)
{
  resList <- list()
  if(get_nulldist)
  {
    nullDist_list <- list()
  }

  for(i in 2:length(df[1,]))
  {
    current_df <- df[,c(1,i)]
    current_df <- current_df[complete.cases(current_df),]
    t_values <- current_df[,2]

    current_weights <- weight_matrix

    names(current_df)[1] <- "ID"
    names(current_weights)[1] <- "ID"

    common_ids <- merge(current_df, current_weights, by = "ID")
    common_ids <- common_ids$ID
    common_ids <- as.character(common_ids)

    row.names(current_df) <- current_df$ID
    current_df <- as.data.frame(current_df[common_ids,-1])

    row.names(current_weights) <- current_weights$ID
    current_weights <- as.data.frame(current_weights[common_ids,-1])

    current_mat <- as.matrix(current_df)
    current_weights <- t(current_weights)

    scores <- as.data.frame(current_weights %*% current_mat)

    null_dist_t <- replicate(k, sample(t_values,length(current_mat[,1]), replace = F))

    null_dist_scores <- current_weights %*% null_dist_t

    if(get_nulldist)
    {
      nullDist_list[[i-1]] <- null_dist_scores
    }

    if(z_scores)
    {
      scores$mean <- apply(null_dist_scores,1,mean)
      scores$sd <- apply(null_dist_scores,1,sd)
      resListCurrent <- (scores[,1]-scores[,2])/scores[,3]
      names(resListCurrent) <- names(weight_matrix[,-1])
      resList[[i-1]] <- resListCurrent
    }
    else
    {
      for(j in 1:length(weight_matrix[,-1]))
      {
        ecdf_function <- ecdf(null_dist_scores[j,])
        scores[j,1] <- ecdf_function(scores[j,1])
      }
      score_probas <- scores*2-1

      resListCurrent <- score_probas[,1]
      names(resListCurrent) <- names(weight_matrix[,-1])
      resList[[i-1]] <- resListCurrent
    }
  }
  names(resList) <- names(df[,-1])
  resDf <- as.data.frame(resList)
  if(get_nulldist)
  {
    names(nullDist_list) <- names(df[,-1])
    return(list(resDf, nullDist_list))
  }
  else
  {
    return(resDf)
  }

}

df_to_viper_regulon <- function(df)
{
  names(df) <- c("feature","pathway","sign")
  df <- df[complete.cases(df),]

  pathway_regulon <- list(0)
  i <- 1
  for(pathway in unique(df$pathway))
  {
    pathway_feature_list <- list(0)
    features <- df[df$pathway == pathway, 3]
    names(features) <- df[df$pathway == pathway, 1]
    pathway_feature_list[[1]] <- features
    pathway_feature_list[[2]] <- rep(1,length(features))
    names(pathway_feature_list) <- c("tfmode","likelihood")

    pathway_regulon[[i]] <- pathway_feature_list
    i <- i+1
  }
  names(pathway_regulon) <- unique(df$pathway)
  return(pathway_regulon)
}

library(piano)
library(parallel)
library(GSEABase)
library(snowfall)
library(readr)

runPIANO <- function(topTable, gene_to_term, nCores = 1000, IDIndex = 1, FCIndex = 2, PvalIndex = 5, TvalIndex = 4, nPerm = 10000)
{

  nCores <- min(c(nCores,detectCores()-1))

  print(paste(nCores, " cpu(s) will be used."), sep = "")

  topTable <- as.data.frame(topTable)

  gene_to_term <- as.data.frame(gene_to_term)
  gene_to_term[,1] <- toupper(gene_to_term[,1])
  gene_to_term[,2] <- toupper(gene_to_term[,2])
  #print(head(gene_to_term))
  geneSet <- loadGSC(gene_to_term)

  myFC <- topTable[,FCIndex]
  names(myFC) <- toupper(topTable[,IDIndex])

  myPval <- topTable[,PvalIndex]
  names(myPval) <- toupper(topTable[,IDIndex])

  myTval <- topTable[,TvalIndex]
  names(myTval) <- toupper(topTable[,IDIndex])

  nPerm <- nPerm-nPerm%%nCores

  print(paste(as.character(nPerm), " permutations will be made (so that there is an integer x such that x*nCores=nPerm)", sep = ""))

  ###Run the GSA
  gsaRes1 <- runGSA(myTval, gsc=geneSet, adjMethod = "fdr", geneSetStat = "mean", ncpus = nCores, nPerm = nPerm)
  #gsaRes1 <- runGSA(myTval, gsc=geneSet, adjMethod = "fdr", geneSetStat = "mean")
  gsaRes2 <- runGSA(myTval, gsc=geneSet, adjMethod = "fdr", geneSetStat = "median", ncpus = nCores, nPerm = nPerm)
  gsaRes3 <- runGSA(myTval, gsc=geneSet, adjMethod = "fdr", geneSetStat = "sum", ncpus = nCores, nPerm = nPerm)
  gsaRes4 <- runGSA(myTval, gsc=geneSet, adjMethod = "fdr", geneSetStat = "maxmean", ncpus = nCores, nPerm = nPerm)
  gsaRes5 <- runGSA(geneLevelStats = myPval, directions = myFC, gsc=geneSet, adjMethod = "fdr", geneSetStat = "fisher", ncpus = nCores, nPerm = nPerm)
  gsaRes6 <- runGSA(geneLevelStats = myPval, directions = myFC, gsc=geneSet, adjMethod = "fdr", geneSetStat = "stouffer", ncpus = nCores, nPerm = nPerm)
  gsaRes7 <- runGSA(geneLevelStats = myPval, directions = myFC, gsc=geneSet, adjMethod = "fdr", geneSetStat = "tailStrength", ncpus = nCores, nPerm = nPerm)
  gsaRes8 <- runGSA(geneLevelStats = myPval, directions = myFC, gsc=geneSet, adjMethod = "fdr", geneSetStat = "wilcoxon", ncpus = nCores, nPerm = nPerm)
  gsaRes9 <- runGSA(myTval, gsc=geneSet, adjMethod = "fdr", geneSetStat = "page", ncpus = nCores, nPerm = nPerm)
  gsaRes10 <- runGSA(geneLevelStats = myPval, directions = myFC, gsc=geneSet, adjMethod = "fdr", geneSetStat = "reporter", ncpus = nCores, nPerm = nPerm)
  gsaRes11 <- runGSA(myTval, gsc=geneSet, adjMethod = "fdr", geneSetStat = "fgsea", ncpus = nCores, nPerm = nPerm)

  resList <- list(gsaRes1,gsaRes2,gsaRes3,gsaRes4,gsaRes5,gsaRes6,gsaRes7,gsaRes8,gsaRes9,gsaRes10,gsaRes11)
  names(resList) <- c("mean","median","sum","maxmean","fisher", "stouffer","tailStrength","wilcoxon","page", "reporter","fgsea")

  ch <- consensusHeatmap(resList,cutoff=50000,method="median", ncharLabel = 1000, cellnote = "medianPvalue", cex = 0.2, plot = FALSE) ##The results are strange

  consensus <- ch$pMat

  return(list(consensus,ch,resList))
}

gmt_to_csv <- function(gmtfile, fast = T)
{
  if(fast)
  {
    genesets = GSEABase::getGmt(con = gmtfile)
    genesets = unlist(genesets)

    gene_to_term =plyr::ldply(genesets,function(geneset){
      temp <- geneIds(geneset)
      temp2 <- setName(geneset)
      temp3 <- as.data.frame(cbind(temp,rep(temp2,length(temp))))

    },.progress = plyr::progress_text())
    names(gene_to_term) <- c("gene","term")
    return(gene_to_term[complete.cases(gene_to_term),])
  }
  else
  {
    genesets = getGmt(con = gmtfile)
    genesets = unlist(genesets)

    gene_to_term <- data.frame(NA,NA)
    names(gene_to_term) <- c("gene","term")
    for (geneset in genesets)
    {
      temp <- geneIds(geneset)
      temp2 <- setName(geneset)
      temp3 <- as.data.frame(cbind(temp,rep(temp2,length(temp))))
      names(temp3) <- c("gene","term")
      gene_to_term <- rbind(gene_to_term,temp3)
    }

    return(gene_to_term[complete.cases(gene_to_term),])
  }
}