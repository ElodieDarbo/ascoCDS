#!/usr/bin/env Rscript 

# Copyright: Patricia Thébaut (patricia.thebault@u-bordeaux.fr) and Elodie Darbo (elodie.darbo@u-bordeaux.fr), Université de Bordeaux, 2024

# Script to run classification tasks using various methods, sampling strategies, and classification levels

# Parse command line arguments
args = commandArgs(trailingOnly=TRUE)

k <- as.numeric(args[1]) # nb fold: 1 to 4
type_sampling <- args[2] #  V1 / V3 / smote
level <- args[3] # subphylum class order superfamily genus species
method.classif <- args[4] # svm rapide knn RF
shap <- as.boolean(args[5]) # TRUE ou FALSE

# Set main directory and results directory
main.dir <- "path to project"
res.dir <- file.path(main.dir,"results",type_sampling,level)
dir.create(file.path(res.dir,method.classif),recursive=TRUE)
setwd(main.dir)

# Load required libraries for classification and data manipulation
library(ggplot2)
library(caret)
library(glmnet)
library(data.table)
library(e1071)
library(class)
library(randomForest)

# Source custom functions for classification tasks
source(file.path(main.dir,"scripts","functions_for_classif.R"))

# Define the features to be used in the classification
features <- c("GCsur", "GCcds", "CpG", "GC1", "GC2", "GC3", "G3skew", "A3skew", # nucleic sequence
                      "Lcod", "F", "Nc", "Psyn", "SCUO", # codon composition
                      "CAI", "CBI", "Fop", # optimal codons
                      "FPC", "FAC", "BOC", "IPC", "IAC", "BIC", # codon context
                      "Gravy", "Aromo", "pI", "II", "class") # physico-chemical properties


# Log the current configuration being run
print("####################")
message("Sampling ", type_sampling ,"; Classif : ", method.classif ,"; Running : ", level, " ; fold : ", k)

# Run the classification task for the given fold, sampling strategy, and classification level
res.classif <- run_all_classifications(level, k, features, res.dir, type_sampling, method.classif,shap)

# Add fold information to the results
res.classif$overall$fold <- k
res.classif$byClass <- res.classif$byClass[,c(1,7,12)]
res.classif$byClass$fold <- k

# Define prefix for saving the results
prefix <- file.path(res.dir,method.classif,paste0("fold",k,"_",level))

# Save overall and class-specific results to tab-separated files
write.table(res.classif$overall,paste0(prefix,"_overall.tab"),sep="\t",quote=F)
write.table(res.classif$byClass,paste0(prefix,"_byClass.tab"),sep="\t",quote=F)

# Save the entire result object for later use
save(res.classif,file=paste0(prefix,".RData"))
