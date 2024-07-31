#!/usr/bin/env Rscript 

args = commandArgs(trailingOnly=TRUE)

k <- as.numeric(args[1]) # nb fold: 1 to 4
type_sampling <- args[2] #  V1 / V3 / smote
level <- args[3] # subphylum class order superfamily genus species
method.classif <- args[4] # svm rapide knn RF
shap <- as.boolean(args[5]) # TRUE ou FALSE

main.dir <- "path to project"
res.dir <- file.path(main.dir,"results",type_sampling,level)
dir.create(file.path(res.dir,method.classif),recursive=TRUE)
setwd(main.dir)

library(ggplot2)
library(caret)
library(glmnet)
library(data.table)
library(e1071)
library(class)
library(randomForest)

source(file.path(main.dir,"scripts","functions_for_classif.R"))

features <- c("GCsur", "GCcds", "CpG", "GC1", "GC2", "GC3", "G3skew", "A3skew", # nucleic sequence
                      "Lcod", "F", "Nc", "Psyn", "SCUO", # codon composition
                      "CAI", "CBI", "Fop", # optimal codons
                      "FPC", "FAC", "BOC", "IPC", "IAC", "BIC", # codon context
                      "Gravy", "Aromo", "pI", "II", "class") # physico-chemical properties


print("####################")
message("Sampling ", type_sampling ,"; Classif : ", method.classif ,"; Running : ", level, " ; fold : ", k)

res.classif <- run_all_classifications(level, k, features, res.dir, type_sampling, method.classif,shap)
res.classif$overall$fold <- k
res.classif$byClass <- res.classif$byClass[,c(1,7,12)]
res.classif$byClass$fold <- k

prefix <- file.path(res.dir,method.classif,paste0("fold",k,"_",level))

write.table(res.classif$overall,paste0(prefix,"_overall.tab"),sep="\t",quote=F)
write.table(res.classif$byClass,paste0(prefix,"_byClass.tab"),sep="\t",quote=F)

save(res.classif,file=paste0(prefix,".RData"))
