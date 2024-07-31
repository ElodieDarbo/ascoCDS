#!/usr/bin/env Rscript 

## define main directory
main.dir <- "path to project"
fold.dir <- file.path(main.dir,"FOLD_DATA")
setwd(main.dir)

source(file.path(main.dir,"scripts","functions_for_classif.R"))

library(themis)
library(groupdata2)
library(data.table)
library(sail)

Data_phylum_Asco <- fread(file.path(main.dir,"CDS_metrics_taxonomy.tab"), data.table = F)

features <- colnames(Data_phylum_Asco[,1:26])

anno <- Data_phylum_Asco[,-(1:26)]
featab <- Data_phylum_Asco[,features]
featab <- as.data.frame(scale(featab))

Data_phylum_Asco[,features] <- featab

if (!dir.exists(fold.dir)){
	message("Creating folds")
	dir.create(fold.dir)
	folds <- createFolds(Data_phylum_Asco$species, k =4)
	for (k in 1:4){
	  Liste_data=NULL
	  print(paste0("FOLD :: ",k))
	  filename=paste(c(REG_FOLD,k,"TRAIN.txt"),sep="_",collapse="_") 
	  ### Pour TRAIN et à chaque FOLD, on sauvegarde les trois autres FOLD
	  writeLines(Data_phylum_Asco$codID[-folds[[k]]],file.path(fold.dir,paste("FOLD", k, "TRAIN.txt" , sep = "_")),quote = F, row.names = F)
	  writeLines(Data_phylum_Asco$codID[folds[[k]]],file.path(fold.dir,paste("FOLD", k, "TEST.txt" , sep = "_")),quote = F, row.names = F)
	}

}

sampling <- c("V1","V3","smote")
levels <- c("subphylum","class", "order", "superfamily", "genus","species")
for (k in 1:4){
	IDs.train <- readLines(file.path(fold.dir,paste("FOLD", k, "TRAIN.txt", sep="_")))[-1]
	Liste_data <- list()
	Liste_data$train <- Data_phylum_Asco[row.names(Data_phylum_Asco)%in%IDs.train,]
	IDs.test <- readLines(file.path(fold.dir,paste("FOLD", k, "TEST.txt" , sep="_")))[-1]
	Liste_data$test <- Data_phylum_Asco[row.names(Data_phylum_Asco)%in%IDs.test,]
	for (for_col in levels){
		for (select_sampling in sampling){
			res.dir <- file.path(main.dir,"results",select_sampling,for_col,"datasets")
			suppressWarnings(dir.create(res.dir,recursive=TRUE))
			message("building dataset for fold",k,", ",select_sampling,", ",for_col)
			name_files <- file.path(res.dir,paste0("fold",k,"_",select_sampling,"_",for_col))
			if (select){
			  name_files <- paste0(name_files,"_no_Eurotiomycetes")
			}
			build_datasets(Liste_data,for_col,features,name_files,select_sampling,over=1)
		}
	}
}







