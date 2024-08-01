# Copyright: Patricia Thébaut (patricia.thebault@u-bordeaux.fr) and Elodie Darbo (elodie.darbo@u-bordeaux.fr), Université de Bordeaux, 2024

# This script give the functions that enable to produce the plots present in the paper doi:xxx (in preparation)

####################
# Function to retrieve and aggregate accuracy metrics
# This function loads the classification results across different folds, sampling strategies, methods, and taxonomic levels,
# and compiles them into a single data table for further analysis.

get_accuracy <- function(){
  # Iterate through sampling strategies, classification methods, and taxonomic levels
  predictions <- rbindlist(lapply(c("V1","V3","smote"),function(s){
    rbindlist(lapply(c("rapide","knn","RF","svm"),function(m,s){
      rbindlist(lapply(c("subphylum","class","order","superfamily","genus","species"),function(l,s,m){
        rbindlist(lapply(1:4,function(k,l,s,m){
          cls <- c("subphylum","class","order","superfamily","genus","species")
          subl <- cls[which(cls==l):6]
          rbindlist(lapply(subl,function(sl,k,l,s,m){
            message("Fold: ",k,", Sampling: ",s,", method: ",m,", objective: ",l,", byCol: ",sl)
            if (file.exists(file.path(s,l,m,paste0("fold",k,"_",l,"_byCol_",sl,"_classif.RData")))){
              # Load the classification results and dataset for the current fold and stratification sublevel
              load(file.path(s,l,m,paste0("fold",k,"_",l,"_byCol_",sl,"_classif.RData")))
              load(file.path(s,l,"datasets",paste0("fold",k,"_",s,"_",l,"_by_col_",sl,".RData")))
              res <- data.table(fold=k,method=m,sampling=s,objective=l,by_sublevel=sl,accuracy=results_liste$overall["Accuracy"],trainsize=nrow(liste_by$train))
            }
            else {
              res <- data.table(fold=k,method=m,sampling=s,objective=l,by_sublevel=sl,accuracy=NA,trainsize=NA)
            }
            return(res)        
          },k,l,s,m))
        },l=l,s=s,m=m))    
      },s=s,m=m))
    },s=s))
  }))
  # Label the sampling strategy and set factors for objective and sublevel
  predictions$sampling.name <- ifelse(as.vector(predictions$sampling)=="V1","BD",ifelse(as.vector(predictions$sampling)=="V3","hBD","hBA"))
  predictions$objective <- factor(as.vector(predictions$objective),levels=c("subphylum","class","order","superfamily","genus","species"))
  predictions$by_sublevel <- factor(as.vector(predictions$by_sublevel),levels=c("subphylum","class","order","superfamily","genus","species"))
  predictions$method <- toupper(predictions$method)
  predictions$method[predictions$method=="RAPIDE"] <- "MLR"
  return(predictions)
}


####################
# Function to create a barplot of accuracy metrics
# This function generates barplots to visualize the accuracy of different classification methods 
# across sampling strategies and taxonomic levels, stratified by sublevels.

barplot_accuracy <- function(ALL_Acc,for_col, Titre,nb){
  cols=c("subphylum"=colours()[403],"class"=colours()[135],"order"=colours()[550],"family"=colours()[144], "genus"=colours()[129],"species"=colours()[256])
  # Subset data for the given classification level
  ALL_Acc_sub <- ALL_Acc[ALL_Acc$objective==for_col,]
  # Create the ggplot barplot with accuracy metrics stratified by sublevel
  pl_sub <- ggplot(ALL_Acc_sub, aes(x = method, y = accuracy, fill = by_sublevel)) +
    geom_col(width = nb-0.1, position = position_dodge(nb-0.1)) + 
    facet_grid(cols = vars(sampling.name), scales="free", space = "free")+
    theme_bw() +
    theme(panel.grid.major = element_line(colour = "grey"),plot.title = element_text(hjust = 0.5))+
    theme(strip.text.x = element_text(size = 10), plot.title = element_text(hjust = 0.5,size = 12), axis.text.x=element_text(color = "black", size=10, hjust=0.5),legend.position="bottom", legend.direction = "horizontal")+
    scale_colour_manual(
      values = cols,
      aesthetics = c("colour", "fill"),
      name="Sublevel of stratified samplings")+
    labs(title = Titre, y="Accuracy")+
    coord_cartesian(ylim = c(0.2,1)) 

  g <- pl_sub + theme(axis.title.x = element_blank())
  return(g)
}

####################
# Function to retrieve F1 scores for classification results
# This function loads and compiles F1 scores from classification tasks across different methods and sampling strategies,
# organizing the results for visualization and comparison.

get_F1 <- function(classif_level="class",method_sampling=c("V3","smote")){
  # Iterate through sampling strategies and classification methods to load results
  all_F1 <- rbindlist(lapply(c(method_sampling),function(s){
    rbindlist(lapply(c("rapide","knn","RF","svm"),function(m,s){
      load(file.path(s,classif_level,m,paste0("fold",k,"_",classif_level,".RData")))
      F1.tmp <- res.classif$byClass
      F1.tmp$M_classif <- m
      F1.tmp$M_samplings <- s
      F1.tmp$Level <- gsub("[1-6]: ","",F1.tmp$Sampling)
      F1.tmp$Class <- gsub("Class: ","",row.names(F1.tmp))
      return(F1.tmp)
    },s=s))
  }))
  # Label the sampling strategy and classification method
  all_F1$M_samplings <- ifelse(as.vector(all_F1$M_samplings)=="V3","hBD sampling","hBA sampling")
  all_F1$M_classif <- toupper(all_F1$M_classif)
  all_F1$M_classif[all_F1$M_classif=="RAPIDE"] <- "MLR"
  all_F1$Level <- factor(as.vector(all_F1$Level),levels=paste("by",c("class","order","superfamily","genus","species")))
  return(all_F1)
}

####################
# Function to retrieve classification predictions
# This function loads the prediction results for each fold, method, and sampling strategy,
# combining them into a comprehensive dataset for analysis.

get_predictions <- function(classif_level="class"){
  # Iterate through folds, sampling strategies, and methods to load prediction results
  predictions.elements <- rbindlist(lapply(c("V1","V3","smote"),function(s){
    rbindlist(lapply(c("rapide","knn","RF","svm"),function(m,s){
      rbindlist(lapply(c("class","order","superfamily","genus","species"),function(l,s,m){
        rbindlist(lapply(1:4,function(k,l,s,m){
          message("Fold: ",k,", Sampling: ",s,", method: ",m,", byCol: ",l)
          load(file.path(s,classif_level,m,paste0("fold",k,"_",classif_level,"_byCol_",l,"_classif.RData")))
          res <- results_liste$prediction
          res$fold <- k
          res$method <- m
          res$sampling <- s
          res$byCol <- l
          return(data.table(res))
        },l=l,s=s,m=m))
        
      },s=s,m=m))
    },s=s))
  }))
}

####################
# Function to reorder data for heatmap plotting
# This function applies hierarchical clustering to reorder rows and columns in a data matrix,
# ensuring that related groups are visually grouped together in the heatmap.

order_for_heatmap <- function(tmp.tab){
  # Perform hierarchical clustering on the data matrix
  clust <- pheatmap(tmp.tab, 
                    cluster_cols = F,
                    clustering_method = "ward.D",
                    silent = T)
  # Reorder data by genus, superfamily, order, and class to enhance visualization
  tmp.tab <- tmp.tab[clust$tree_row$order,]
  anno.species.tmp <- anno.species.tmp[clust$tree_row$order,]

  tmp.tab <- tmp.tab[order(anno.species.tmp$genus),]
  anno.species.tmp <- anno.species.tmp[order(anno.species.tmp$genus),]


  tmp.tab <- tmp.tab[order(anno.species.tmp$superfamily),]
  anno.species.tmp <- anno.species.tmp[order(anno.species.tmp$superfamily),]


  tmp.tab <- tmp.tab[order(anno.species.tmp$order),]
  anno.species.tmp <- anno.species.tmp[order(anno.species.tmp$order),]

  tmp.tab <- tmp.tab[order(anno.species.tmp$class),]
  anno.species.tmp <- anno.species.tmp[order(anno.species.tmp$class),]

  return(list(mat=tmp.tab,anno=anno.species.tmp))
}

####################
# Function to plot a heatmap
# This function generates a heatmap with hierarchical clustering, applying custom color schemes
# to differentiate between classes, species, and metric groups.

plot_heatmap <- function(tmp.tab,colors,annot_colors,anno.species.tmp){
  # Generate the heatmap with custom annotations and color schemes
  print(pheatmap(tmp.tab,cluster_rows = F,cluster_cols = F, 
               show_rownames = T,
               show_colnames = F,
               color = colors,
               annotation_row = anno.species.tmp,
               annotation_col = anno.class[,"class",drop=F],
               clustering_distance_rows = "correlation",
               clustering_method = "ward.D",
               breaks=seq(0,1,0.01),
               border_color = NA,
               fontsize=8,
               annotation_colors = annot_colors
               )
  )
}

####################
# Function to plot SHAP values
# This function visualizes SHAP values for model explainability, highlighting the impact of different metrics 
# on predictions across taxonomic classes, with custom color schemes for metric groups and classes.

plot_SHAP_values <- function(shp.dir,res.dir,k,select_sampling,met_classif,for_col,by_col){
  # Load the training dataset and SHAP values for the current classification model
  load(file.path(res.dir,"datasets",paste0("fold",k,"_",select_sampling,"_",for_col,"_by_col_",by_col,".RData"))) # liste_by
  x_train = liste_by$train
  features <- colnames(x_train)[1:26]
  x_train <- x_train[,c(features,for_col)]
  y_train = as.factor(liste_by$train[, for_col])

  load(file.path(shp.dir,paste0("shaps_class_",met_classif,"_",by_col,".RData"))) # all.shaps.melt
  # Prepare the SHAP values and metric groups for plotting
  color.metrics <- c("codon composition"=colours()[132],"codon context"=colours()[76],"nucleic sequence"="darkgrey","optimal codons"=colours()[35],"physico-chemical"=colours()[617])
  metrics.tab <- rbind(data.frame(metrics=c("GCsur", "GCcds", "CpG", "GC1", "GC2", "GC3", "G3skew", "A3skew"),group="nucleic sequence"),
                       data.frame(metrics=c("Lcod", "F","Nc", "Psyn", "SCUO"),group="codon composition"),
                       data.frame(metrics=c("CAI", "CBI", "Fop"),group="optimal codons"),
                       data.frame(metrics=c("FPC", "FAC", "BOC", "IPC", "IAC", "BIC"),group="codon context"),
                       data.frame(metrics=c("Gravy", "Aromo", "pI", "II"),group="physico-chemical"))

  # Adjust metrics grouping to ensure consistency
  all.shaps.melt$metrics.group[all.shaps.melt$metrics.group=="physico-chemical properties"] <- "physico-chemical"
  # Melt the training data and SHAP values for plotting
  train.melt <- melt(x_train,id.vars="class")
  train.melt <- train.melt[order(train.melt$class),]
  all.shaps.melt <- all.shaps.melt[order(as.vector(all.shaps.melt$class))]

  # Ensure that the variable order matches between the melted data frames
  print(all.equal(as.vector(train.melt$variable),as.vector(all.shaps.melt$variable)))

  # Combine SHAP values with the corresponding feature values for plotting
  all.shaps.melt$measure <- train.melt$value

   # Define colors for each class in the plot 
  fills <- c(Saccharomycetes=colours()[33],Pneumocystidomycetes=colours()[614],Schizosaccharomycetes=colours()[613],Taphrinomycetes=colours()[104],Dothideomycetes=colours()[132],Sordariomycetes=colours()[131],Xylonomycetes=colours()[124],Eurotiomycetes=colours()[430])

  # Cap extreme values for better visualization
  all.shaps.melt$measure[all.shaps.melt$measure>2 | all.shaps.melt$measure < -2] <- sign(all.shaps.melt$measure[all.shaps.melt$measure>2 | all.shaps.melt$measure < -2]) * 2

   # Create the ggplot object for SHAP value visualization 
  set.seed(20240708)
  g <- ggplot(all.shaps.melt[sample(1:nrow(all.shaps.melt),round(nrow(all.shaps.melt)/2)),],aes(x=value,y=variable)) + theme_bw() +
    geom_jitter(aes(color=measure),size=0.1,alpha=0.2) + geom_vline(xintercept = 0,linetype="dashed",linewidth=0.2) +
    facet_grid(metrics.group~class,scales = "free_y", switch = 'y') + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1),legend.position = "top",strip.text = element_text(color="white")) +
    labs(y="Metrics") +
    scale_colour_gradient(low = colours()[142],high = "purple") 

  # Adjust the strip background colors for classes and metric groups in the plot
  g <- ggplot_gtable(ggplot_build(g))
  stripr <- which(grepl('strip-t', g$layout$name))
  fills <- c(Saccharomycetes=colours()[33],Pneumocystidomycetes=colours()[614],Schizosaccharomycetes=colours()[613],Taphrinomycetes=colours()[104],Dothideomycetes=colours()[132],Sordariomycetes=colours()[131],Xylonomycetes=colours()[124],Eurotiomycetes=colours()[430])
  k <- 1
  for (i in stripr) {
    j <- which(grepl('rect', g$grobs[[i]]$grobs[[1]]$childrenOrder))
    g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
    k <- k+1
  }

  stripr <- which(grepl('strip-l', g$layout$name))
  fills <- color.metrics[names(color.metrics)%in%all.shaps.melt$metrics.group]
  k <- 1
  for (i in stripr) {
    j <- which(grepl('rect', g$grobs[[i]]$grobs[[1]]$childrenOrder))
    g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
    k <- k+1
  }

  # Adjust the heights of different sections of the plot for better readability
  g$heights[10] <- g$heights[10]*0.7
  g$heights[c(11,13,15,17)] <- g$heights[10]*0
  g$heights[12] <- g$heights[12]*0.9

  g$heights[16] <- g$heights[16]*0.55
  g$heights[18] <- g$heights[18]*0.65
  
  # Return the final ggplot object
  return(g)
}



