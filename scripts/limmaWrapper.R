#This is a script designed to provide usefull functions to run LIMMA.
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

library(limma)

checkInputs <- function(measurments, targets)
{
  if(class(measurments) != "data.frame")
  {
    error_message <- paste("The measurments argument should be a data.frame. It's currently a", paste(class(measurments), ".",sep = ""))
    return(list(FALSE, error_message))
  }
  else
  {
    if(dim(measurments)[1] == 0)
    {
      error_message <- "The measurments dataframe doesn't seem to contain any measurments..."
      return(list(FALSE, error_message))
    }
    else
    {
      if(dim(measurments)[2] == 0)
      {
        error_message <- "The measurments dataframe doesn't seem to contain any samples..."
        return(list(FALSE, error_message))
      }
      else
      {
        if(class(as.matrix(measurments)[,1]) != "numeric")
        {
          return(list(FALSE, "The measurments dataframe should contain only numerical values (or NAs)."))
        }
        else
        {
          if(class(targets) != "data.frame")
          {
            error_message <- paste("The targets argument should be a data.frame. It's currently a", paste(class(targets), ".",sep = ""))
            return(list(FALSE, error_message))
          }
          else
          {
            if(dim(targets)[2] < 2)
            {
              return(list(FALSE,"The targets dataframe should have at least two columns, sample names and conditions."))
            }
            else
            {
              if(dim(targets)[1] != dim(measurments)[2])
              {
                error_message <- paste("The targets dataframe should have as many samples (targets rows) as the measurements (measurments columns). Currently, the targets dataframe has", paste(dim(targets)[1], "samples and the measurements have", paste(dim(measurments)[2],"samples.")))
                return(list(FALSE, error_message))
              }
              else
              {

              }
            }
          }
        }
      }
    }
  }
  return(list(TRUE, "All seems to be in order..."))
}

normalise <- function(measurments, targets, batch = NULL, method = 3)
{
  input_check <- checkInputs(measurments, targets)
  if (input_check[[1]])
  {

  }
  else
  {
    print(input_check[[2]])
    return(input_check[[1]])
  }
}

makeContrastsAlt <- function(targets, comparisons)
{
  cont.matrix <- matrix(0,nrow = length(unique(targets$condition)), ncol = length(comparisons))
  i <- 1
  for (comparison in comparisons)
  {
    for (j in 1:length(comparison))
    {
      cont.matrix[abs(comparison[j]),i] <- cont.matrix[abs(comparison[j]),i]+(comparison[j]/abs(comparison[j]))
    }
    i <- i + 1
  }
  return(cont.matrix)
}

poolContrasts <- function(cont.matrix, pool)
{
  max_size <- 1
  for (pooler in pool)
  {
    max_size <- ifelse(max_size > length(pooler), max_size, length(pooler))
  }

  max_iterations <- (max_size-1)*length(pooler)

  for (k in 1:max_iterations)
  {
    for (pooler in pool)
    {
      cont.matrix <- apply(cont.matrix, 2, function(x)
      {
        if (x[pooler[1]] == 1)
        {
          x[pooler[2]] <- -1
          return(x)
        }
        if (x[pooler[1]] == -1)
        {
          x[pooler[2]] <- 1
          return(x)
        }
        return(x)
      })
    }
  }
  return(cont.matrix)
}

runLimma <- function(measurements, targets, comparisons = NULL, pool = NULL, regress_out = NULL)
{
  input_check <- checkInputs(measurements, targets)
  if (input_check[[1]]) #input has correct format
  {
    if (!is.null(comparisons))
    {
      if (!is.null(regress_out))
      {
        for (regressor in regress_out)
        {
          measurements <- removeBatchEffect(measurements, targets[,regressor])
        }
      }

      cont.matrix <- makeContrastsAlt(targets, comparisons)

      if (!is.null(pool))
      {
        cont.matrix <- poolContrasts(cont.matrix, pool)
      }

      cont.matrix <- as.data.frame(cont.matrix)
      row.names(cont.matrix) <- unique(targets$condition)
      cont.matrix <- as.matrix(cont.matrix)

      fcond <- factor(targets$condition, levels = unique(targets$condition))

      design <- model.matrix(~0+fcond)
      design <- as.data.frame(design)
      names(design) <- unique(targets$condition)
      design <- as.matrix(design)

      print(cont.matrix)

      fit <- lmFit(measurements, design)
      fit2 <- contrasts.fit(fit, cont.matrix)
      fit2 <- eBayes(fit2)

      return(list(fit2, cont.matrix, fit))
    }
  }
  else
  {
    print(input_check[[2]])
    return(input_check[[1]])
  }
}


ttopFormatter <- function(ttop)
{
  ttop$ID <- row.names(ttop)
  ttop <- ttop[,c(7,1,2,3,4,5,6)]
  ttop <- ttop[complete.cases(ttop),]
  return(ttop)
}

library(ggplot2)
library(ggrepel)
#df = result from limma differnecial expression (with column X as identifier column)
#gtt = gene to term table
#consensus = consensus table of piano analisys (with column X as term column)
#outpath = path to the directory were the figures will be stored
#mt = mapping table between the identifiers of limma analisys and gene to term table

#'\code{volcano_nice}
#'
#'This function is designed to generate a stylish and practical volcano plot to visual the result of a differential analisys, with simple required inputs.
#'The plot features two bi-symptotic curves to give a visual support for p-value and fold change threshold.
#'
#'@param df a dataframe of n*m dimension, where n is the number of omic features (genes,protein,etc...) and m is at least 3 (identifiers, foldchanges and p-values).
#'@param hAss the p-value threshold (example : 0.05)
#'@param FCIndex the column number corresponding to the foldchanges
#'@param pValIndex the column number corresponding to the p-values
#'@param IDIndex the column number corresponding to the identifiers
#'@param vAss the foldchange threshold (example : 0.5)
#'@param label a boolean argument to indicate if the labels of the significant genes should be displayed or not (maximum 30 labels)
#'@param straight a boolean argument to indicate if the plot should feature bi-symptotic curves or straight lines to visualise the thresholds.
#'@param nlabels number of labels to display
#'@param manual_labels a vector of labels that will be kept even if they are not in the top "nlabels"
#'
#'@return a ggplot object of the volcano plot
volcano_nice <- function (df, hAss = 0.05, FCIndex, pValIndex, IDIndex, vAss = NULL,
                          label = FALSE, straight = FALSE, nlabels, manual_labels = NA)
{
  df <- df[complete.cases(df), ]
  names(df)[1] <- "X1"
  hAssOri <- hAss
  hAss <- -log(hAss)
  names(df) <- gsub("adj.P.Val", "FDR", names(df))
  names(df)[FCIndex] <- "logFC"
  names(df)[pValIndex] <- "adj.P.Val"
  if (max(abs(df[, FCIndex])) >= 1) {
    xlimAbs <- ceiling(max(abs(df[, FCIndex])))
    ylimAbs <- ceiling(max(abs(-log(df[, pValIndex]))))
  }
  else {
    xlimAbs <- max(abs(df[, FCIndex]))
    ylimAbs <- max(abs(-log(df[, pValIndex])))
  }
  if (is.null(vAss)) {
    vAss <- xlimAbs/10
  }
  xneg <- function(x) abs(hAss - 1 + x/(x + vAss))
  xpos <- function(x) abs(hAss - 1 + x/(x - vAss))
  test <- function(x, y, vAss) {
    if (x < -vAss) {
      if (xneg(x) < -log(y)) {
        return("1")
      }
      else {
        return("0")
      }
    }
    else {
      if (x > vAss) {
        if (xpos(x) < -log(y)) {
          return("1")
        }
        else {
          return("0")
        }
      }
      else {
        return("0")
      }
    }
  }
  if (straight) {
    df$couleur <- ifelse(abs(df$logFC) >= vAss & df$adj.P.Val <=
                           hAssOri, "1", "0")
  }
  else {
    df$couleur <- "0"
    df$couleur <- apply(df, 1, FUN = function(x) test(as.numeric(x[FCIndex]),
                                                      as.numeric(x[pValIndex]), vAss))
  }
  df <- df[order(df$adj.P.Val, decreasing = F), ]
  df$condLabel <- df[, IDIndex]
  df[df$couleur == "0", "condLabel"] <- NA
  labels_to_keep <- c(df[c(1:nlabels), "condLabel"],manual_labels)
  df[!(df$condLabel %in% labels_to_keep), "condLabel"] <- NA
  df$couleur <- ifelse(df$couleur == "1" & df$logFC < 0, "2",
                       df$couleur)
  if (label) {
    a <- ggplot(df, aes(x = logFC, y = -log(adj.P.Val), color = couleur)) +
      geom_point(alpha = 0.5) + geom_label_repel(aes(label = condLabel)) +
      stat_function(fun = xneg, xlim = c(-xlimAbs, -vAss),
                    color = "black", alpha = 0.7) + ylim(c(0, ylimAbs)) +
      xlim(c(-xlimAbs, xlimAbs)) + stat_function(fun = xpos,
                                                 xlim = c(vAss, xlimAbs), color = "black", alpha = 0.7) +
      scale_colour_manual(values = c("grey30", "red",
                                     "royalblue3")) + theme_minimal() + theme(legend.position = "none")
  }
  else {
    if (straight) {
      a <- ggplot(df, aes(x = logFC, y = -log(adj.P.Val),
                          color = couleur)) + geom_point(alpha = 0.5) +
        geom_vline(xintercept = -vAss, color = "blue") +
        geom_vline(xintercept = vAss, color = "blue") +
        ylim(c(0, ylimAbs)) + xlim(c(-xlimAbs, xlimAbs)) +
        geom_hline(yintercept = hAss, color = "red") +
        scale_colour_manual(values = c("grey30", "red",
                                       "royalblue3")) + theme_minimal() + theme(legend.position = "none")
    }
    else {
      a <- ggplot(df, aes(x = logFC, y = -log(adj.P.Val),
                          color = couleur)) + geom_point(alpha = 0.5) +
        stat_function(fun = xneg, xlim = c(-xlimAbs,
                                           -vAss), color = "black", alpha = 0.7) + ylim(c(0,
                                                                                          ylimAbs)) + xlim(c(-xlimAbs, xlimAbs)) + stat_function(fun = xpos,
                                                                                                                                                 xlim = c(vAss, xlimAbs), color = "black", alpha = 0.7) +
        scale_colour_manual(values = c("grey30", "red",
                                       "royalblue3")) + theme_minimal() + theme(legend.position = "none")
    }
  }
  return(a)
}

