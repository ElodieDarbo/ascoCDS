# Copyright: Patricia Thébaut (patricia.thebault@u-bordeaux.fr) and Elodie Darbo (elodie.darbo@u-bordeaux.fr), Université de Bordeaux, 2024

# This script enable to to prepare datasets depending on sampling strategies and run all classification tasks

####################################################################################################
####################################################################################################
###################################### TRAINING and TEST SETS PREPARATION ##########################
####################################################################################################
####################################################################################################
# Main function
build_datasets <- function(Liste_data,for_col,features,name_files,select_sampling,over=1){
  strates <- rev(c("subphylum","class", "order", "superfamily", "genus","species"))
  s <- which(strates==for_col)
  set.seed(2023)
  Liste_data$test <- downsample(Liste_data$test,cat_col = for_col, id_method = "n_rows_c")
    print(colnames(Liste_data$test))
    res <- lapply(1:s,function(x,Liste_data,for_col,features,strates,name_files,select_sampling){
    cl <- strates[x]
    name_files <- paste0(name_files,"_by_col_",cl)
    if (file.exists(paste0(name_files,".RData"))){
      message("File ",paste0(name_files,".RData")," already exists")
    } else {
      message("Generating : ",paste0(name_files,".RData"))
      liste_by <- sample_data(Liste_data,cl,for_col,select_sampling,features,annotation,over)
      message("by col: ",cl,", nb cds: ",nrow(liste_by$train))
      liste_by$test=Liste_data$test
      message("Saving : ",paste0(name_files,".RData"))
      save(liste_by,file=paste0(name_files,".RData"))
    }
    },Liste_data,for_col,features,strates,name_files,select_sampling)
}


#######################
# downsample categories by keeping balanced subcategories
sample_data_by_stratification <- function(Data,by_col,for_col,min_f_classes) {
  
  f_classes= names(table(Data[for_col]))
  one_class=NULL
  tmp_train=NULL
  for (cl in 1:length(f_classes)){
    one_class=Data[Data[for_col]==f_classes[cl],]
    nb_classes= dim(table(droplevels(one_class[by_col])))
    set.seed(116)
    df=balance(one_class,size=as.integer(min_f_classes/nb_classes),cat_col=by_col)
    tmp_train=rbind(tmp_train,df)
  }
  return(tmp_train)
  
}

#######################
# find.thr.dataset
find.thr.dataset <- function(train,for_col,opt,features=NULL,over=NULL){
  all.levels <- c("subphylum","class", "order", "superfamily", "genus","species")
  i <- which(all.levels==for_col)
  switch(opt,
         V3={
            mins <- unlist(lapply(i:length(all.levels),function(x,train,all.levels,for_col){
             by_col <- all.levels[x]
             test <- as.data.frame.matrix(table(train[,for_col],train[,by_col]))
             test[test==0] <- NA
             min.byCol <- min(test,na.rm=T)
             min.forCol <- min(apply(test,1,function(x) sum(!is.na(x))))
             max.element <- min.byCol*min.forCol
             return(max.element)
           },train,all.levels,for_col))
           res <- min(mins)
         },
         smote={
           mins <- unlist(lapply(i:length(all.levels),function(x,train,all.levels,for_col){
             by_col <- all.levels[x]
             test <- as.data.frame.matrix(table(train[,for_col],train[,by_col])) 
             test.over <- as.data.frame.matrix(table(train[,for_col],train[,by_col]))* over
             test[test==0] <- NA
             max.byCol <- max(test.over,na.rm=T)
             test[test<max.byCol] <- max.byCol
             #min.forCol <- min(apply(test,1,function(x) sum(!is.na(x))))
             #min.byCol <- min(test,na.rm=T)
             max.forCol <- min(apply(test,1,function(x) sum(!is.na(x))))
             max.byCol <- max(test,na.rm=T)
             
             #max.element <- min.byCol*min.forCol
             max.element <- max.byCol*max.forCol
             return(max.element)
           },train,all.levels,for_col))
           res <- min(mins)
         }
  )
  return(res)
}


#######################
# apply sampling according to stategies
sample_data <- function(Liste_data, by_col, for_col, select_sampling, features,annotation=NULL,over) {
    if (select_sampling == "V1") {
      set.seed(20231)
      df2 = balance(Liste_data$train, size = "min", cat_col = by_col)
      df = balance(df2, size = "min", cat_col = for_col)
    }
    else if (select_sampling == "V3") {
      set.seed(20233)
      min.per.levels <- find.thr.dataset(Liste_data$train, for_col, opt = "V3",over=over)
      df = sample_data_by_stratification(Liste_data$train, by_col, for_col, min.per.levels)
    }
    else  if (select_sampling == "smote") {
      Data_tmp = Liste_data$train[, features]
      Data_tmp[, by_col] = as.factor(Liste_data$train[, by_col])
      set.seed(20231122)
      Data_tmp = smote(Data_tmp, var = by_col, over_ratio = over)
      annotation <- unique(Liste_data$train[,c(for_col,cl)])
      Data_tmp[,for_col] <- annotation[match(Data_tmp[,cl], annotation[,cl]),for_col]  
      min.per.levels <- find.thr.dataset(Liste_data$train, for_col, opt = "smote", features,over)
      print(min.per.levels)
      df = sample_data_by_stratification(Data_tmp, by_col, for_col, min.per.levels)
    }
    newList <- list("train" = df,
                      "test" = Liste_data$test,
                      "Features" = features)    
    return(newList)
}


####################################################################################################
####################################################################################################
#################### CLASSIFICATION, EVALUATION, EXPLAINABILITY ####################################
####################################################################################################
####################################################################################################
# Main function to launch classification tasks, model evaluations and explainability for all stratification sublevels 
run_all_classifications <- function(for_col,k,features,res.dir,select_sampling,met_classif,shap=F){
  strates <- rev(c("subphylum","class", "order", "superfamily", "genus","species"))
  s <- which(strates==for_col)
  results_liste <- mclapply(1:s,function(x,k,for_col,features,strates,res.dir,select_sampling,met_classif,n_neighbors){
    cl <- strates[x]
    i <- which(c("subphylum","class", "order", "superfamily", "genus","species")==cl)
    data_file <- file.path(res.dir,"datasets",paste0("fold",k,"_",select_sampling,"_",for_col,"_by_col_",cl,".RData"))
    res_file <- file.path(res.dir,met_classif,paste0("fold",k,"_",for_col,"_byCol_",cl))
    message("loading ",data_file)
    load(data_file)   
    message("Running ",met_classif)
    results_liste <- run_one_classif(liste_by,features,for_col,cl,res.dir,met_classif)
    
    if (shap){
      message("Running Shapley")
      shapleys <- shapleys_for_train(met_classif,results_liste$fit,liste_by$train,cl)
      save(shapleys,file=file.path(res.dir,paste0("SHAPs/shaps_",for_col,"_",met_classif,"_",cl,".RData")))
    }

    df <- as.data.frame(results_liste$byClass)
    df$Sampling <- paste0(i,": by ",cl)
    
    message("Writing ",paste0(res_file,"_byClass.txt"))
    write.table(df,file=paste0(res_file,"_byClass.txt"),sep=";")
    df_overall=as.data.frame(results_liste$overall)
    df_overall$by_class <- cl
    
    message("Writing ",paste0(res_file,"_overall.txt"))
    write.table(df_overall, file=paste0(res_file,"_overall.txt") ,sep=";")
    
    message("Saving ",paste0(res_file,"_classif.RData"))
    save(results_liste,df,df_overall,file=paste0(res_file,"_classif.RData"))
    return(list(df,df_overall))
  },k,for_col,features,strates,res.dir,select_sampling,met_classif,n_neighbors,mc.cores=s)
  
  overall <- do.call("rbind",lapply(results_liste,function(x) {
    x[[2]]
  }))
  byClass <- do.call("rbind",lapply(results_liste,function(x) {
    x[[1]]
  }))
  return (list(overall=overall, byClass = byClass))
}

####################
# Classification task and model evaluation
run_one_classif <- function(liste_sample, features, for_col, by_col, name_file, met_classif) {
    corr <- cor(liste_sample$train[, features])
    F_corr = caret::findCorrelation(corr, cutoff = 0.9)
    x_train = liste_sample$train[, features]
    y_train = as.factor(liste_sample$train[, for_col])
    x_test = liste_sample$test[, features]
    y_test = as.factor(liste_sample$test[, for_col])
    
    if (met_classif == "rapide") {
      multinom_model <-glmnet(as.data.frame(x_train),
               y_train,
               family = "multinomial",
               alpha = 1)#### OCT_23
      last_lam = multinom_model$lambda[length(multinom_model$lambda)]
      y_predict_test = predict(
        multinom_model,
        newx = as.matrix(x_test),
        type = "response",
        s = last_lam,
        exact = FALSE
      )
      fit <- multinom_model
      pred_label_test = levels(as.factor(y_test))[max.col(y_predict_test[, , 1])]
    }
    else if (met_classif == "svm")
    {
      model <- svm(
          x_train,
          y_train,
          method = "C-classification",
          kernal = "radial",
          gamma = 0.1,
          cost = 10,
          probability=T
        )
      pred_label_test <- predict(model, x_test)
      fit <- model

    }
    else if (met_classif == "RF")
    {
      model.rf = randomForest(x_train, y_train, ntree = 500)
      pred_label_test <- predict(model.rf, newdata = x_test)
      fit <- model.rf
    }
    
    else if (met_classif == "knn")
    {
      pred_label_test <- knn(
        train = x_train,
        test = x_test,
        cl = y_train,
        k = 10
      )
    }
    conf = caret::confusionMatrix(table(pred_label_test, y_test))
    conf.mat <- as.data.frame.matrix(table(pred_label_test, y_test))
    x_test <- as.data.frame(x_test)
    x_test$pred <- pred_label_test
    x_test$true <- y_test
    x_test$codID <- liste_sample$test$codID
    x_test$species <- liste_sample$test$species
    return (list(
      "overall" = conf$overall,
      "byClass" = conf$byClass,
      conf.mat = conf.mat,
      prediction = x_test,
      fit=fit
    ))
}

###################################################################
# Explainability task: compute shapley values
shapleys_for_train <- function(met_classif,fit,x_train,for_col,by_col){
  categories <- unique(x_train[,for_col])
  all.shaps <- rbindlist(lapply(categories,function(lev,met_classif,fit,x_train){
    message(lev)
    switch(met_classif,
           rapide={
             pfun <- function(object, newdata) {
               last_lam = object$lambda[length(object$lambda)]
               test <- predict(object, newx = as.matrix(newdata),type = "link",
                               s = last_lam,
                               exact = FALSE)
               pred <- test[,,1]
               pred <- pred[,lev]
               unname(pred)
             }
           },
           svm={
             pfun <- function(object, newdata) {
               pred <- predict(object, newdata = as.matrix(newdata),probability=T)
               unname(attr(pred, "probabilities")[,lev])
             }
           },
           RF={
             pfun <- function(object, newdata) {
               pred <- predict(object, newdata = as.matrix(newdata),type="prob")
               unname(pred[,lev])
             }
               }
    )
    
    
    shp <- explain(fit, 
                       X = x_train, 
                       pred_wrapper = pfun,
                       newdata=x_train[y_train==lev,]
    )
    shp <- as.data.frame(shp)
    shp$class <- lev
    shp
  },met_classif,fit,x_train))
  all.shaps.melt <- melt(all.shaps)
  metrics.tab <- rbind(data.frame(metrics=c("GCsur", "GCcds", "CpG", "GC1", "GC2", "GC3", "G3skew", "A3skew"),group="nucleic sequence"),
                     data.frame(metrics=c("Lcod", "F","Nc", "Psyn", "SCUO"),group="codon composition"),
                     data.frame(metrics=c("CAI", "CBI", "Fop"),group="optimal codons"),
                     data.frame(metrics=c("FPC", "FAC", "BOC", "IPC", "IAC", "BIC"),group="codon context"),
                     data.frame(metrics=c("Gravy", "Aromo", "pI", "II"),group="physico-chemical properties"))
  all.shaps.melt$metrics.group <- metrics.tab$group[match(all.shaps.melt$variable,metrics.tab$metrics)]
  all.shaps.melt$variable <- factor(as.vector(all.shaps.melt$variable),levels=metrics.tab$metrics)
  all.shaps.melt$sublevel <- by_col
  return(all.shaps.melt)
}
