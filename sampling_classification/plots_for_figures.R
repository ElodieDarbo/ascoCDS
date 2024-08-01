# Copyright: Patricia Thébaut (patricia.thebault@u-bordeaux.fr) and Elodie Darbo (elodie.darbo@u-bordeaux.fr), Université de Bordeaux, 2024

# This script enable to produce the plots present in the paper doi:xxx (in preparation)

## define main directory
main.dir <- "path to project"
data.dir <- file.path(main.dir,"data")
res.dir <- file.path(main.dir,"results","plots")

source(file.path(main.dir,"scripts","myUtils.R"))
source(file.path(main.dir,"scripts","functions_for_figure.R"))

library(data.table)
library(factoextra)

library(ggplot2)
library(pheatmap)
library(cowplot)
library(FactoMineR)

### get data
# raw CDS metrics
CDS_raw <- fread(file.path(data.dir,"CDS_metrics_taxonomy.tab"))
# get accuracy values
accuracies <- get_accuracy()
# get F1 values at class level for sample V3 and smote
all_F1 <- get_F1(classif_level="class",method_sampling=c("V3","smote"))
# get predictions at class level
predictions.elements <- get_predictions(classif_level="class")

########################################### Figure 1 ###########################################
# Panel A: CDS number per species per subphylum
count.per.species <- CDS_raw[,list(nb_CDS=length(codID)),by=c("subphylum","species")]

ggplot(count.per.species, aes(x=subphylum, y=nb_CDS, fill=subphylum))+
  scale_fill_manual(values=c("Saccharomycotina"=colours()[137],"Taphrinomycotina"=colours()[139],"Pezizomycotina"=colours()[132]))+
  scale_colour_manual(values=c("Saccharomycotina"=colours()[137],"Taphrinomycotina"=colours()[139],"Pezizomycotina"=colours()[132]))+
  stat_halfeye(adjust = 0.5, # permet de régler le lissage
               width=0.5, # permet de gérer la hauteur des curves
               justification = -0.2, # permet de déplacer les courbes sur la droite
               # permet de supprimer l'affichage d'un intervalle présent par défaut
               .width = 0, point_colour = NA) + 
  labs(y = "Number of CDS by species ",x="")+ 
  theme_bw()+ 
  theme(axis.title.y =element_text( size=12), axis.text.x = element_text(size = 12,face="bold" ) )+
  geom_boxplot(
  width = 0.12, # gère la largeur des boites
  outlier.color ="black", # passe les outliers en orange
  outlier.shape = 20,
  outlier.size = 0.6,
  # ajoute une transparence dans les boites (elles sont plus claires que l'aire des courbes de densités)
  alpha = 0.2) + 
  geom_point(aes(colour=subphylum), # permet d'ajouter une couleur sur le spoints
             size = 1.3, # gère la taille des points
             alpha = .3, #ajoute de la transparence sur les points
             position = position_jitter( # permet d'obtenir des points décalés
               seed = 1, # permet de toujours obtenir la même représnetation aléatoire des points 
               width = .09 # permet de gérer la largeur du décalage
             ))

export.plot(file.path(res.dir,"Fig1A"),width=10,height=5,export.formats = c("jpeg"))


########################################### Figure 3 ###########################################
# compute median CDS values per species

ordered.metrics <- c("GCsur", "GCcds", "CpG", "GC1", "GC2", "GC3", "G3skew", "A3skew", # nucleic sequence
					"Lcod", "F", "Nc", "Psyn", "SCUO", # codon composition
					"CAI", "CBI", "Fop", # optimal codons
                    "FPC", "FAC", "BOC", "IPC", "IAC", "BIC", # codon context
                    "Gravy", "Aromo", "pI", "II", "class") # physico-chemical properties

level <- "subphylum"

CDS_raw <- CDS_raw[,c(ordered.metrics,"species","subphylum"),with=F]
CDS_median <- CDS_raw[,lapply(.SD, median),by=c("subphylum","species")]
CDS_median <- as.data.frame(CDS_median)

CDS <- CDS_median[,c(ordered.metrics,level)]

# compute MFA
groups <- c(8,5, 3, 6, 4,1)
t_group <- c("s","s","s","s","s","n")

res.MFA <- MFA(CDS,groups, type = t_group, ind.sup = NULL, 
                 name.group = c("nucleic sequence","codon composition","optimal codons","codon context","physico-chemical properties","Taxonomy"), 
                 num.group.sup = NULL, 
                 graph = FALSE)


##### Panel A 
contrib1 <- fviz_contrib(res.MFA, choice = "quanti.var", axes = 1, top = 20, palette = "jco", title = "Contribution to Dim 1",xtickslab.rt = 90) + 
  theme(legend.text = element_text(size = 8)) + theme(plot.title = element_text(hjust = 0.5))
contrib2 <-  fviz_contrib(res.MFA, choice = "quanti.var", axes = 2, top = 20, palette = "jco", title = "Contribution to Dim 2",xtickslab.rt = 90) + 
  theme(legend.text = element_text(size = 8)) + theme(plot.title = element_text(hjust = 0.5)) 
contrib3 <-  fviz_contrib(res.MFA, choice = "quanti.var", axes = 3, top = 20, palette = "jco", title = "Contribution to Dim 3",xtickslab.rt = 90) + 
  theme(legend.text = element_text(size = 8)) + theme(plot.title = element_text(hjust = 0.5)) 

contrib <- plot_grid(contrib1,contrib2,contrib3,align = 'h', ncol = 3)


##### Panel B
Var1 <- fviz_mfa_var(res.MFA, "quanti.var", col.var = "contrib", 
                    col.var.sup = "violet", repel = TRUE,
                    geom = c("point", "text"))  + theme(plot.title = element_text(hjust = 0.5))

Ind_12 <- fviz_mfa_ind(res.MFA,
                      habillage  = level, # colorer par groupes
                      addEllipses = TRUE, # Ellipses de concentrations
                      ellipse.type = "convex",
                      legend.title = level,
                      repel = TRUE,
                      label = "none"
) + guides(fill=guide_legend(ncol =1))  + theme(plot.title = element_text(hjust = 0.5))

MFA12 <- plot_grid(Var1,Ind_12)


##### Panel C
Ind_13 <- fviz_mfa_ind(res.MFA,
                      axes = c(1,3),
                      habillage  = level, # colorer par groupes
                      addEllipses = TRUE, # Ellipses de concentrations
                      ellipse.type = "convex",
                      repel = TRUE,
                      label = "none") + 
	guides(fill=guide_legend(ncol =1))  + theme(plot.title = element_text(hjust = 0.5))


Var2 <- fviz_mfa_var(res.MFA, "quanti.var", col.var = "contrib",
                    col.var.sup = "violet", repel = TRUE, axes = c(1,3),
                    geom = c("point", "text"))  + theme(plot.title = element_text(hjust = 0.5))

MFA13 <- plot_grid(Var2,Ind_13)


# combine panels
plot_grid(contrib, MFA12,MFA13, nrow = 3)
ggsave(filename = file.path(res.dir,"Fig3.jpeg",sep=""),width = 1600, height = 1230,dpi = 300, units = "px", bg = 'white')

########################################### Figure 4 ###########################################
# compute median accuracy accross folds
med.fold <- accuracies[,list(accuracy=median(accuracy)),by=colnames(accuracies)[-1]]

# su
wi <- 1
p1 <- barplot_accuracy(med.fold,"subphylum", "SUBPHYLUM level",1*wi)
p2 <- barplot_accuracy(med.fold,"class", "CLASS level",5/6*wi)
p3 <- barplot_accuracy(med.fold,"order", "ORDER level",4/6*wi)
p4 <- barplot_accuracy(med.fold,"family", "FAMILY level",3/6*wi)
p5 <- barplot_accuracy(med.fold,"genus", "GENUS level",2/6*wi)
p6 <- barplot_accuracy(med.fold,"species", "SPECIES level",1/6*wi+0.05)

ggarrange(p1,p2,p3,p4,p5,p6,
                    labels = c("A", "B", "C","D","E","F"),
                    ncol = 2, nrow = 3,
                    common.legend = TRUE, legend = "bottom")

export.plot(file.path(res.dir,"Fig4"),width=14,height=17,export.formats = c("jpeg"))


########################################### Figure 5 ###########################################
# panel A: compare by stratification
predictions$params <- paste(predictions$method,predictions$objective,predictions$sampling.name,sep="_")
predictions$params.med <- paste(predictions$method,predictions$objective,sep="_")

predictions.med <- predictions[,list(mean.acc=mean(accuracy,na.rm=T),med.acc=median(accuracy,na.rm=T)),by=params.med]
predictions$relative.accuracy <- predictions$accuracy-predictions.med$med.acc[match(predictions$params.med,predictions.med$params.med)]

g <- ggplot(predictions[fold==1 & objective!="species"],aes(x=by_sublevel,y=relative.accuracy,group=params)) + theme_bw() +
  geom_point(aes(color=sampling.name)) + 
  geom_line(aes(color=sampling.name)) + 
  #geom_linerange(aes(ymin = Q25, ymax = Q75,color=method)) +
  facet_grid(method~objective,scales="free_x") + labs(x="Stratification sublevel",y="relative accuracy") +
  scale_colour_manual(values=colours()[c(24,124,613)],labels=c("BD","hBD","hBA")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

g <- ggplot_gtable(ggplot_build(g))

g$widths[9] <- g$widths[7]*0.9
g$widths[11] <- g$widths[9]*0.9
g$widths[13] <- g$widths[11]*0.85
g$widths[15] <- g$widths[13]*0.8

grid::grid.draw(g)

export.plot(file.path(res.dir,"stratif_comparison_inv"),width=8,height=5,export.formats = c("jpeg"))

# panel B: compare by size / BD sampling
V1 <- predictions[sampling=="BD"]
V1.med <- V1[,list(accuracy=mean(accuracy),trainsize=mean(trainsize)),by=c("method","objective","by_sublevel")]
V1.med <- V1.med[,list(med.acc=median(accuracy)),by=c("method","objective")]

V1 <- merge(V1,V1.med,by=c("method","objective"))
V1$acc.rel <- V1$accuracy-V1$med.acc
V1$objective <- factor(as.vector(V1$objective),levels=c("subphylum","class","order","superfamily","genus","species"))

ggplot(V1[objective!="species",list(accuracy=mean(acc.rel),trainsize=mean(trainsize)/1000),by=c("method","objective","by_sublevel")],aes(x=trainsize,y=accuracy)) + theme_bw() +
  geom_point(aes(color=method)) + labs(x="training set size (K)",y="relative accuracy") +
  geom_smooth(method="lm",fill=NA,aes(color=method))+ 
  facet_wrap(~objective,scales="free_x",ncol = 5) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ylim(c(-0.1,0.1)) + 
  scale_colour_manual(values=colours()[c(24,124,613,53)],labels=c("KNN","MLR","RF","SVM"))

export.plot(file.path(res.dir,"V1_size_comparison"),width=8,height=2,export.formats = c("jpeg"))

########################################### Figure 6 ###########################################
sublevel_cols=c("by subphylum"=colours()[403],"by class"=colours()[144],"by order"=colours()[135],"by superfamily"=colours()[550], "by genus"=colours()[256],"by species"=colours()[129]    )
ordered.class <- c("Saccharomycetes","Pneumocystidomycetes","Schizosaccharomycetes","Taphrinomycetes","Dothideomycetes","Sordariomycetes","Xylonomycetes","Eurotiomycetes")

all_F1$Class <- factor(as.vector(all_F1$Class),levels=ordered.class)

ggplot(all_F1, aes(x = Class, y = F1, color = Level)) + 
      geom_boxplot() + 
      theme_bw() + 
      facet_grid(cols = vars(M_samplings), row = vars(M_classif))+
      theme_bw() +
      theme(panel.grid.major = element_line(colour = "lightgrey"),
            plot.title = element_text(hjust = 0.5),
            strip.text.x = element_text(size = 12), 
            plot.title = element_text(hjust = 0.5), 
            axis.text.x=element_text(colour=c(rep("darkred",1), rep("darkolivegreen",3), rep("darkblue",4)), 
                                      angle=45, 
                                      size=12, 
                                      hjust=0.9),
            legend.direction = "vertical",
            legend.title = element_text(size=12),legend.text = element_text(size=10)) +
      labs(fill=labels_sampling)  +
      scale_color_manual("Sub-levels of stratified samplings",values=sublevel_cols)

export.plot(file.path(res.dir,"Fig6"),width=10,height=8,export.formats = c("jpeg"))

########################################### Figure 7 ###########################################
predictions.elements$class <- factor(as.vector(predictions.elements$class),levels=ordered.class)
predictions.elements$order <- CDS_raw[match(predictions.elements$codID,CDS_raw$codID),]$order
predictions.elements$genus <- CDS_raw[match(predictions.elements$codID,CDS_raw$codID),]$genus
predictions.elements$superfamily <- CDS_raw[match(predictions.elements$codID,CDS_raw$codID),]$superfamily


# count CDS per species
species.counts <- predictions.elements[,list(count=length(codID)),by=c("fold","method","sampling","byCol","species")]
# average accross folds
species.counts <- species.counts[,list(sum.count=sum(count),mean.count=mean(count)),by=c("method","sampling","byCol","species")]

predictions.elements <- predictions.elements
predictions.elements <- predictions.elements[,list(method,sampling,fold,byCol,species,true,pred)]
# count good and bad predicted CDS per species
predictions.elements$count <- 1
predictions.elements <- predictions.elements[,list(count=sum(count)),by=c("fold","method","sampling","byCol","species","true","pred")]
predictions.elements <- predictions.elements[,list(sum.count=sum(count),mean.count=mean(count),min.count=min(count),max.count=max(count)),by=c("method","sampling","byCol","species","true","pred")]
predictions.elements$params <- paste(predictions.elements$sampling,predictions.elements$method,predictions.elements$byCol,sep="_")

predictions.elements$total <- species.counts[match(predictions.elements$species,species)]$mean.count
predictions.elements$mean.freq <- predictions.elements$mean.count / predictions.elements$total 

all.sp <- unique(CDS_raw$species)

# make aggregated counts: species ~ prediction
pred.per.species <- lapply(unique(predictions.elements$params),function(p,pred.class,ordered.class,all.sp){
  print(p)
  tmp <- pred.class[params==p]
  tmp <- as.data.frame.matrix(xtabs(mean.freq~species+pred,tmp))
  if (sum(!all.sp%in%row.names(tmp))>0){
    sp.0 <- as.data.frame(matrix(0,ncol=8,nrow=sum(!all.sp%in%row.names(tmp)),dimnames = list(all.sp[!all.sp%in%row.names(tmp)],c(colnames(tmp)))))
    tmp <- rbind(tmp,sp.0)
    tmp <- tmp[order(row.names(tmp)),]
  }
  return(tmp[,ordered.class])
},pred.class=predictions.elements,ordered.class=ordered.class,all.sp=all.sp)

names(pred.per.species) <- unique(predictions.elements$params)

# create annotation for rows
anno.species <- as.data.frame(unique(CDS_raw[,c("subphylum","class",'order',"superfamily","genus","species")]))
row.names(anno.species) <- anno.species$species
anno.species$class <- factor(anno.species$class,levels=ordered.class)
anno.species <- anno.species[order(anno.species$class),]

# create annotation for columns
anno.c <- as.data.frame(unique(CDS_raw[,c("class","subphylum")]))
anno.class <- rbind(anno.c,anno.c,anno.c,anno.c,anno.c)
row.names(anno.class) <- c(anno.c$class,paste(anno.c$class,rep(1:4,each=8),sep="."))
row.names(anno.class) <- anno.class$class

## Analysis of prediction from RF after hBA
ptmp <- "smote_RF"
tmp.tabs <- pred.per.species[grep(ptmp,id.tmp)]
colnames(tmp.tab) <- sub(paste0(ptmp,"_[a-z]+[.]"),"",colnames(tmp.tab))
tmp.tab <- na.omit(tmp.tab[anno.species$species,])

## order rows to plot heatmap
ordered.obj <- order_for_heatmap(tmp.tab)
tmp.tab <- ordered.obj$mat
anno.species.tmp <- ordered.obj$anno

## get colors for sublevel annotations
color.dev = grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)]
order <- sample(color.dev, length(unique(anno.species.tmp$order)))
names(order) <- unique(anno.species.tmp$order)

superfamily <- sample(color.dev, length(unique(anno.species.tmp$superfamily)))
names(superfamily) <- unique(anno.species.tmp$superfamily)

genus <- sample(color.dev, length(unique(anno.species.tmp$genus)))
names(genus) <- unique(anno.species.tmp$genus)

class.colors <- c(Saccharomycetes=colours()[33],Pneumocystidomycetes=colours()[614],Schizosaccharomycetes=colours()[613],Taphrinomycetes=colours()[104],Dothideomycetes=colours()[132],Sordariomycetes=colours()[131],Xylonomycetes=colours()[124],Eurotiomycetes=colours()[430])

## plot for sublevel class
class.tab <- tmp.tab[,1:8]
colors <- colorRampPalette(c("white",colours()[c(34,555,556)]))(100)
annot_colors <- list(class=class.colors, order=order)
anno.rows <- anno.species.tmp[,c("class","order")]
plot_heatmap(class.tab,colors,annot_colors,anno.rows)
export.plot(file.path(res.dir,"Fig7_byClass"),width=6,height=5)


## plot for sublevel order
class.tab <- tmp.tab[,9:16]
colors <- colorRampPalette(c("white",colours()[c(97,98,99)]))(100)
annot_colors <- list(class=class.colors, order=order)
anno.rows <- anno.species.tmp[,c("class","order")]
plot_heatmap(class.tab,colors,annot_colors,anno.rows)
export.plot(file.path(res.dir,"Fig7_byOrder"),width=6,height=5)

## plot for differences between sublevel class and order
class.tab <- tmp.tab[,1:8] - tmp.tab[,9:16]
class.tab <- class.tab[anno.species.tmp$class!="Saccharomycetes",]
colors <- colorRampPalette(c(colours()[c(99,98,97)],"white",colours()[c(34,555,556)]))(100)
annot_colors <- list(class=class.colors, order=order,superfamily=superfamily, genus=genus)
anno.rows <- anno.species.tmp[,c("class","order","superfamily","genus")]
plot_heatmap(class.tab,colors,annot_colors,anno.rows)
export.plot(file.path(res.dir,"Fig7_diff_byClassVsOrder_zoom"),width=6,height=5)

########################################### Figure 8 ###########################################
shp.dir <- file.path(res.dir,"SHAPs")
results.dir <- file.path(main.dir,"results",select_sampling,for_col)

k <- 1
select_sampling <- "smote"
for_col <- "class"
by_col <- "class"


##### MLR
met_classif <- "rapide"

g <- plot_SHAP_values(shp.dir,results.dir,k,select_sampling,met_classif,for_col,by_col)
grid::grid.draw(g)

export.plot(file.path(res.dir,"Fig8A"),height=10,width=13,export.formats=c("jpeg"))

##### RF
met_classif <- "RF"

g <- plot_SHAP_values(shp.dir,results.dir,k,select_sampling,met_classif,for_col,by_col)
grid::grid.draw(g)

export.plot(file.path(res.dir,"Fig8B"),height=10,width=13,export.formats=c("jpeg"))


##### SVM
met_classif <- "svm"

g <- plot_SHAP_values(shp.dir,results.dir,k,select_sampling,met_classif,for_col,by_col)
grid::grid.draw(g)

export.plot(file.path(res.dir,"Fig8C"),height=10,width=13,export.formats=c("jpeg"))



########################################### SupplementaryFigure 1 ###########################################
# panel A
load(file.path(results.dir,"datasets",paste0("fold",k,"_",select_sampling,"_",for_col,"_by_col_",by_col,".RData"))) # liste_by
x_train = liste_by$train
x_train <- x_train[,c(ordered.metrics,for_col)]

cor.feat <- cor(x_train[,ordered.metrics],method = "spearman")

color.metrics <- c("codon composition"=colours()[132],"codon context"=colours()[76],"nucleic sequence"="darkgrey","optimal codons"=colours()[35],"physico-chemical"=colours()[617])
metrics.tab <- rbind(data.frame(metrics=c("GCsur", "GCcds", "CpG", "GC1", "GC2", "GC3", "G3skew", "A3skew"),group="nucleic sequence"),
                     data.frame(metrics=c("Lcod", "F","Nc", "Psyn", "SCUO"),group="codon composition"),
                     data.frame(metrics=c("CAI", "CBI", "Fop"),group="optimal codons"),
                     data.frame(metrics=c("FPC", "FAC", "BOC", "IPC", "IAC", "BIC"),group="codon context"),
                     data.frame(metrics=c("Gravy", "Aromo", "pI", "II"),group="physico-chemical"))

row.names(metrics.tab) <- metrics.tab$metrics
metrics.tab$order <- match(metrics.tab$group,names(color.metrics))
metrics.tab <- metrics.tab[order(metrics.tab$order),]

clust <- pheatmap(cor.feat,
         clustering_distance_rows = "correlation",
         clustering_distance_cols = "correlation",
         clustering_method = "ward.D",
         silent=T)

cor.feat <- cor.feat[clust$tree_row$order,clust$tree_col$order]
metrics.tab <- metrics.tab[match(row.names(cor.feat),metrics.tab$metrics),]

cor.feat <- cor.feat[order(metrics.tab$order),order(metrics.tab$order)]

pheatmap(cor.feat,
         cluster_rows = F,
         cluster_cols = F,
         breaks=seq(-1,1,0.02),color = colorRampPalette(colours()[c(131,1,34)])(100),
         annotation_row = metrics.tab[,"group",drop=F],
         annotation_col = metrics.tab[,"group",drop=F],
         annotation_colors = list(group=color.metrics))

export.plot(file.path(res.dir,"SuppFig1A"),height=8,width=8,export.formats=c("jpeg"))

# panel B
train.melt <- melt(x_train,id.vars=for_col)
train.melt$metrics.group <- metrics.tab$group[match(train.melt$variable,metrics.tab$metrics)]
train.melt$metrics.group <- factor(as.vector(train.melt$metrics.group),levels=names(color.metrics))
train.melt <- train.melt[order(train.melt$metrics.group)]
train.melt$variable <- factor(as.vector(train.melt$variable),levels=unique(train.melt$variable))

set.seed(20240708)
g <- ggplot(train.melt[sample(1:nrow(train.melt),round(nrow(train.melt)/2)),],aes(x=class,y=measure)) + theme_bw() +
  geom_violin(aes(color=class,fill=class),alpha=0.2,draw_quantiles = c(0.25,0.5,0.75)) + geom_hline(yintercept = 0,linetype="dashed",linewidth=0.2) +
  facet_wrap(~variable,ncol=6) + 
  theme(axis.text.x = element_blank(),axis.ticks.x = element_blank(),legend.position = "top",strip.text = element_text(color="white")) +
  labs(y="Metrics") +
  scale_colour_manual(values=fills)  +
  scale_fill_manual(values=fills)  +
  coord_cartesian(ylim=c(-5,5)) 

g <- ggplot_gtable(ggplot_build(g))
stripr <- which(grepl('strip-t', g$layout$name))

fills.c <- rev(color.metrics[metrics.tab$group])
fills.c <- fills.c[c(1:2,rev(3:8),rev(9:14),rev(15:20),rev(21:26))]

m <- 1
for (i in stripr[-c(3:6)]) {
  j <- which(grepl('rect', g$grobs[[i]]$grobs[[1]]$childrenOrder))
  g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills.c[m]
  m <- m+1
}
grid::grid.draw(g)

export.plot(file.path(res.dir,"SuppFig1B"),height=8,width=10,export.formats=c("jpeg"))






