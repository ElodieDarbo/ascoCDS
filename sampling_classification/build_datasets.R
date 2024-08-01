#!/usr/bin/env Rscript 

# Copyright: Patricia Thébaut (patricia.thebault@u-bordeaux.fr) and Elodie Darbo (elodie.darbo@u-bordeaux.fr), Université de Bordeaux, 2024

# Script to run dataset building following the three sampling strategies 
# for all combinations of level - sublevels stratification


## Define the main directory and set working directory
main.dir <- "path to project"
fold.dir <- file.path(main.dir,"FOLD_DATA")
setwd(main.dir)

# Source necessary functions for classification
source(file.path(main.dir,"scripts","functions_for_classif.R"))

# Load required libraries
library(themis)
library(groupdata2)
library(data.table)
library(sail)

# Load the main dataset containing CDS metrics and taxonomy
Data_phylum_Asco <- fread(file.path(main.dir,"CDS_metrics_taxonomy.tab"), data.table = F)

# Define features and scale them for consistency
features <- colnames(Data_phylum_Asco[,1:26])

anno <- Data_phylum_Asco[,-(1:26)]
featab <- Data_phylum_Asco[,features]
featab <- as.data.frame(scale(featab))

# Replace original features with the scaled features
Data_phylum_Asco[,features] <- featab

# Create and save training and test folds if they do not already exist
if (!dir.exists(fold.dir)){
	message("Creating folds")
	dir.create(fold.dir)
	folds <- createFolds(Data_phylum_Asco$species, k =4)
	for (k in 1:4){
	  Liste_data=NULL
	  print(paste0("FOLD :: ",k))
	  filename=paste(c(REG_FOLD,k,"TRAIN.txt"),sep="_",collapse="_") 
	  writeLines(Data_phylum_Asco$codID[-folds[[k]]],file.path(fold.dir,paste("FOLD", k, "TRAIN.txt" , sep = "_")),quote = F, row.names = F)
	  writeLines(Data_phylum_Asco$codID[folds[[k]]],file.path(fold.dir,paste("FOLD", k, "TEST.txt" , sep = "_")),quote = F, row.names = F)
	}

}

# Define the sampling strategies and classification levels
sampling <- c("V1","V3","smote")
levels <- c("subphylum","class", "order", "superfamily", "genus","species")

# Iterate over each fold to build datasets based on different sampling strategies and classification levels
for (k in 1:4){
	# Load training and test IDs for the current fold
	IDs.train <- readLines(file.path(fold.dir,paste("FOLD", k, "TRAIN.txt", sep="_")))[-1]
	Liste_data <- list()
	Liste_data$train <- Data_phylum_Asco[row.names(Data_phylum_Asco)%in%IDs.train,]
	IDs.test <- readLines(file.path(fold.dir,paste("FOLD", k, "TEST.txt" , sep="_")))[-1]
	Liste_data$test <- Data_phylum_Asco[row.names(Data_phylum_Asco)%in%IDs.test,]
	# Iterate over each classification level and sampling strategy to build datasets
	for (for_col in levels){
		for (select_sampling in sampling){
			# Define the directory to save the results
			res.dir <- file.path(main.dir,"results",select_sampling,for_col,"datasets")
			suppressWarnings(dir.create(res.dir,recursive=TRUE))
			message("building dataset for fold",k,", ",select_sampling,", ",for_col)
			# Generate file names and build datasets
			name_files <- file.path(res.dir,paste0("fold",k,"_",select_sampling,"_",for_col))
			build_datasets(Liste_data,for_col,features,name_files,select_sampling,over=1)
		}
	}
}







