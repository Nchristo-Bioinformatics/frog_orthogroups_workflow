library(rjson)
library(tidyverse)
library(topGO)
library(RJSONIO)
#############RELAX RESULTS#############
########################################
jsonFiles <- list.files(path = './RELAX', pattern = "*.json", full.names = TRUE)
jsonFiles <- sort(jsonFiles, decreasing = TRUE)
jsonFiles <- jsonFiles[sapply(jsonFiles, file.size) > 0]
relaxJSONProcessing <- function(i) {
  # Read in the JSON file:
  OG0000698_relax  <- rjson::fromJSON(file=i)
  
  # Create a data frame with the information for each branch tested:
  OG0000698_relaxResult <- data.table::rbindlist(OG0000698_relax[["branch attributes"]][["0"]], fill = TRUE, idcol = TRUE)
  #OG0000698_absrelResult <- subset(OG0000698_absrelResult, select = -c(`Rate Distributions`))
  OG0000698_relaxResult <- OG0000698_relaxResult %>% drop_na('original name')
  OG0000698_relaxResult <- distinct (OG0000698_relaxResult)
  
  #Add a column to hold the orthogroup name:
  orthogroup <- sapply(strsplit(OG0000698_relax[["input"]][["file name"]],"/"), tail, 1)
  orthogroup <- sapply(strsplit(orthogroup, "[.]"), head, 1)
  OG0000698_relaxResult$Orthogroup <- orthogroup
  OG0000698_relaxResult$"original name" <- sub("Scutiger_cf_sikimmensis","Scutigercf_sikimmensis",OG0000698_relaxResult$"original name")
  #species <- paste(unlist(sapply(strsplit(OG0000698_absrelResult[["original name"]],"_"), head, 2)),collapse="_")
  species <- sapply(strsplit(OG0000698_relaxResult[["original name"]],"_"), function(x) paste(x[1:2], collapse = '_')) 
  #species <- as.factor(species)
  #species <- as.vector(species)
  #species <- paste(unlist(species),collapse="_")
  OG0000698_relaxResult$Species <- species
  #OG0000698_absrelResult$Species <- sub(",", "_", OG0000698_absrelResult$Species)
  # Add a column to label genes as having been part of the foreground or background
  
  
  return(OG0000698_relaxResult)
}
possiblyRelaxJSONProcessing <- possibly(relaxJSONProcessing, otherwise = "File empty.")

relaxResults_small <- map(jsonFiles, possiblyRelaxJSONProcessing)
relaxResults_small <- as.data.frame(do.call(rbind, relaxResults_small))
relaxResults_small$branch <- as.factor(relaxResults_small$.id)

rep_na_mean <- function(vec){
  m <- mean(vec, na.rm = TRUE) 
  vec[is.na(vec)] <- m 
  return(vec) 
}


relax_sp_kvalue <-relaxResults_small %>%
  dplyr::select(Species,Orthogroup, 'k (general descriptive)') %>%
  #mutate(omega_sig=ifelse(omega_sig==0,1e-18,omega_sig)) %>%
  spread(Orthogroup, 'k (general descriptive)') %>%
  mutate_at(vars(-Species),rep_na_mean)
log_use_0s_omega <- function(column){
  column[column == 0] <- 1
  return(log(column))
}
bsrel_sp_kvalue_logtrans <- relax_sp_kvalue  %>%
  mutate_at(vars(-Species),log_use_0s_omega)
bsrel_sp_kvalue_logtrans_forpca <- bsrel_sp_kvalue_logtrans %>% as.data.frame
write.csv(bsrel_sp_kvalue_logtrans_forpca, 'kvalue_transformed_values_pca.csv')
#bsrel_sp_kvalue_logtrans_forpca <- read.csv('omega_transformed_values_pca.csv')
rownames(bsrel_sp_kvalue_logtrans_forpca) <- bsrel_sp_kvalue_logtrans_forpca$Species
bsrel_sp_kvalue_logtrans_forpca <- subset(bsrel_sp_kvalue_logtrans_forpca, select = -c(Species))
bsrel_sp_kvalue_logtrans_forpca <- bsrel_sp_kvalue_logtrans_forpca[,which(colSums(bsrel_sp_kvalue_logtrans_forpca) != 0)]

pca_kvalue <- prcomp(bsrel_sp_kvalue_logtrans_forpca)



#Variance
summary(pca_kvalue)



#Save kvalue results
# write_csv(data.frame(pca_kvalue$rotation),"pca_sp_kvalue_loadings_pristi.csv")
# write_csv(data.frame(pca_kvalue$x),"pca_sp_kvalue_coordinates_pristi.csv")
# save(pca_kvalue,relax_sp_kvalue,bsrel_sp_kvalue_logtrans, bsrel_sp_kvalue_logtrans_forpca,file="pca_sp_kvalue_all_res_pristi.rDat")


save(pca_kvalue,relax_sp_kvalue,bsrel_sp_kvalue_logtrans, bsrel_sp_kvalue_logtrans_forpca,file="pca_sp_kvalue_all_res_dirdev.rDat")
write_csv(data.frame(pca_kvalue$rotation),"pca_sp_kvalue_loadings_dirdev.csv")
write_csv(data.frame(pca_kvalue$x),"pca_sp_kvalue_coordinates_dirdev.csv")

###START HERE FROM RSTUDIO###
library(ape)
library(phytools)
library(geiger)
library(ggplot2)

#setwd('~/Orthogroup_paper_edits_2022/TO_DO_FOR_PAPER/simple_results_for_table_processing_r/pristimantis')

#load("pca_sp_kvalue_all_res_pristi.rDat")

setwd('~/Orthogroup_paper_edits_2022/TO_DO_FOR_PAPER/simple_results_for_table_processing_r/direct_developers')

load("pca_sp_omega_all_res_dirdev.rDat")
pca_sp <- pca_kvalue
tree <- read.tree('SpeciesTree_rooted.txt')
plotBranchbyTrait<-function(tree,x,mode=c("edges","tips","nodes"),palette="rainbow",legend=TRUE,xlims=NULL,...){
  mode<-mode[1]
  if(mode=="tips"){
    x<-c(x[tree$tip.label],fastAnc(tree,x))
    names(x)[1:length(tree$tip.label)]<-1:length(tree$tip.label)
    XX<-matrix(x[tree$edge],nrow(tree$edge),2)
    x<-rowMeans(XX)
  } else if(mode=="nodes"){
    XX<-matrix(x[tree$edge],nrow(tree$edge),2)
    x<-rowMeans(XX)
  }
  # begin optional arguments
  if(hasArg(tol)) tol<-list(...)$tol
  else tol<-1e-6
  if(hasArg(prompt)) prompt<-list(...)$prompt
  else prompt<-FALSE
  if(hasArg(type)) type<-list(...)$type
  else type<-"phylogram"
  if(hasArg(show.tip.label)) show.tip.label<-list(...)$show.tip.label
  else show.tip.label<-TRUE
  if(hasArg(show.node.label)) show.node.label<-list(...)$show.node.label
  else show.node.label<-FALSE
  if(hasArg(edge.width)) edge.width<-list(...)$edge.width
  else edge.width<-4
  if(hasArg(edge.lty)) edge.lty<-list(...)$edge.lty
  else edge.lty<-1
  if(hasArg(font)) font<-list(...)$font
  else font<-3
  if(hasArg(cex)) cex<-list(...)$cex
  else cex<-par("cex")
  if(hasArg(adj)) adj<-list(...)$adj
  else adj<-NULL
  if(hasArg(srt)) srt<-list(...)$srt
  else srt<-0
  if(hasArg(no.margin)) no.margin<-list(...)$no.margin
  else no.margin<-TRUE
  if(hasArg(label.offset)) label.offset<-list(...)$label.offset
  else label.offset<-0.01*max(nodeHeights(tree))
  if(hasArg(underscore)) underscore<-list(...)$underscore
  else underscore<-FALSE
  if(hasArg(x.lim)) x.lim<-list(...)$x.lim
  else x.lim<-NULL
  if(hasArg(y.lim)) y.lim<-list(...)$y.lim
  else y.lim<-if(legend&&!prompt) c(1-0.06*length(tree$tip.label),length(tree$tip.label)) else NULL
  if(hasArg(direction)) direction<-list(...)$direction
  else direction<-"rightwards"
  if(hasArg(lab4ut)) lab4ut<-list(...)$lab4ut
  else lab4ut<-"horizontal"
  if(hasArg(tip.color)) tip.color<-list(...)$tip.color
  else tip.color<-"black"
  if(hasArg(plot)) plot<-list(...)$plot
  else plot<-TRUE
  if(hasArg(rotate.tree)) rotate.tree<-list(...)$rotate.tree
  else rotate.tree<-0
  if(hasArg(open.angle)) open.angle<-list(...)$open.angle
  else open.angle<-0
  if(hasArg(use.edge.length)) use.edge.length<-list(...)$use.edge.length
  else use.edge.length<-TRUE
  # end optional arguments
  if(palette=="heat.colors") cols<-heat.colors(n=1000)
  if(palette=="gray") cols<-gray(1000:1/1000)
  if(palette=="rainbow")	cols<-rainbow(1000,start=0.7,end=0) # blue->red
  if(palette=="custom") {colfunct<-colorRampPalette(c("#3D52A1","#77B7E5","#FFFAD2","#ED875E","#AE1C3E"))
  cols <- colfunct(1000)}
  if(is.null(xlims)) xlims<-range(x)+c(-tol,tol)
  breaks<-0:1000/1000*(xlims[2]-xlims[1])+xlims[1]
  whichColor<-function(p,cols,breaks){
    i<-1
    while(p>=breaks[i]&&p>breaks[i+1]) i<-i+1
    cols[i]
  }
  colors<-sapply(x,whichColor,cols=cols,breaks=breaks)
  par(lend=2)
  # now plot
  tree <- compute.brlen(tree)
  xx<-plot.phylo(tree,type=type,show.tip.label=show.tip.label,show.node.label=show.node.label,edge.color=colors,
                 edge.width=edge.width,edge.lty=edge.lty,font=font,cex=cex,adj=adj,srt=srt,no.margin=no.margin,root.edge=root.edge,
                 label.offset=label.offset,underscore=underscore,x.lim=x.lim,y.lim=y.lim,direction=direction,lab4ut=lab4ut,
                 tip.color=tip.color,plot=plot,rotate.tree=rotate.tree,open.angle=open.angle,lend=2,new=FALSE,use.edge.length=FALSE)
  if(legend==TRUE&&is.logical(legend)) legend<-round(15*max(nodeHeights(tree)),2)
  if(legend){
    if(hasArg(title)) title<-list(...)$title
    else title<-"trait value"
    if(hasArg(digits)) digits<-list(...)$digits
    else digits<-1
    add.color.bar(legend,cols,title,xlims,digits,prompt=prompt)
  }
  invisible(xx)
}


species_info_df <- read.csv('absrel_counts_both_hypotheses.csv', row.names=1, header=TRUE)
species_info_df$Species <- rownames(species_info_df)
sp_names_PC1 <- pca_sp$x[,"PC1"]
names(sp_names_PC1) <- rownames(species_info_df)
sp_names_PC2 <- pca_sp$x[,"PC2"]
names(sp_names_PC2) <- rownames(species_info_df)
sp_names_PC3 <- pca_sp$x[,"PC3"]
names(sp_names_PC3) <- rownames(species_info_df)
plotBranchbyTrait(tree,sp_names_PC1,mode="tips",palette = "custom",use.edge.length=FALSE)


selection <- species_info_df$selection_avg
counts <- species_info_df$value
grp <- as.factor(species_info_df$Clade)
grp_pristi <- as.factor(species_info_df$Pristi)
selection_percent <- as.vector(species_info_df$selection_percent_pristi)
selection_percent_sigonly <- as.vector(species_info_df$selection_percent_pristi_sigonly)
branch_percentages <- as.vector(species_info_df$Proportion_sig_branches)
sp_names_PC1_vector <- as.vector(sp_names_PC1)
names(grp) <- rownames(species_info_df)
names(counts) <- rownames(species_info_df)
names(selection) <- rownames(species_info_df)
names(grp_pristi) <- rownames(species_info_df)
names(selection_percent) <- rownames(species_info_df)
names(selection_percent_sigonly) <- rownames(species_info_df)	
names(branch_percentages) <- rownames(species_info_df)
names(sp_names_PC1_vector) <- rownames(species_info_df)

###### Are the k values different for foreground vs background species? ###PHYLANOVA

model_kvalues_pristi <- phylANOVA(tree, grp_pristi, sp_names_PC1, nsim=1000, posthoc=TRUE, p.adj="BH")
# ANOVA table: Phylogenetic ANOVA
# 
# Response: y
# Sum Sq   Mean Sq   F value Pr(>F)
# x        274.8488 274.84882 15.574404  0.069
# Residual 688.2513  17.64747                 
# 
# P-value based on simulation.
# ---------
#   
#   Pairwise posthoc test using method = "BH"
# 
# Pairwise t-values:
#   Other Pristimantis
# Other         0.000000     3.946442
# Pristimantis -3.946442     0.000000
# 
# Pairwise corrected P-values:
#   Other Pristimantis
# Other        1.000        0.069
# Pristimantis 0.069        1.000
model_kvalues_dirdev <- phylANOVA(tree, grp, sp_names_PC1, nsim=1000, posthoc=TRUE, p.adj="BH")

# Response: y
# Sum Sq   Mean Sq   F value Pr(>F)
# x        265.3368 265.33681 14.830437  0.095
# Residual 697.7633  17.89137                 
# 
# P-value based on simulation.
# ---------
#   
#   Pairwise posthoc test using method = "BH"
# 
# Pairwise t-values:
#   Direct_developer     Other
# Direct_developer         0.000000 -3.851031
# Other                    3.851031  0.000000
# 
# Pairwise corrected P-values:
#   Direct_developer Other
# Direct_developer            1.000 0.095
# Other                       0.095 1.000
# ---------

#########Correlations with PC axes and life history traits
species_info_df$Colour <- as.factor(species_info_df$Colour)
cols<-setNames(c("blue","red"),
               levels(species_info_df$Colour))
species_info_df$elevation <- ifelse(species_info_df$elevation == 0 , 1, species_info_df$elevation)
sp_coord <- pca_sp$x %>%
  as.data.frame %>%
  rownames_to_column("sp_abbr") %>%
  as.tibble
sp_coord$sp_abbr <- species_info_df$Species
sp_coord_anno <- sp_coord %>%
  left_join(species_info_df,by=c("sp_abbr" = "Species")) %>%
  mutate(log_elevation = elevation)

# sp_coord_anno %>%
#   ggplot(aes(log(elevation),PC1)) +
#   geom_point() +
#   xlab("log(elevation)")
# ggsave("PC1_sp_omega_log_elevation.pdf",width=6,height=4)


sp_coord_anno_df <- sp_coord_anno %>%
  as.data.frame()
rownames(sp_coord_anno_df) <- sp_coord_anno_df$sp_abbr
elevation_scaled <- setNames(sp_coord_anno_df[,"log_elevation"],rownames(sp_coord_anno_df))
PC1_scaled <- setNames(sp_coord_anno_df[,"PC1"],rownames(sp_coord_anno_df))

##############Plot phylomorphospace plot of elevation by % selected PC1
pdf("PC1_sp_kvalue_elevation_phylomorphospace.pdf",width=8,height=8)
phylomorphospace(tree,cbind(elevation_scaled ,PC1_scaled),node.size=c(0,1),label="off",xlab="Elevation", ylab="Weight PC1")
points(elevation_scaled,PC1_scaled,pch=16,col=species_info_df$Colour)
legend("topleft", legend=c("Pristimantis"), col=c("red"), pch=19)
title(main="Phylomorphospace of elevation and PC1 k-values")
dev.off()
########CORELATE % SELECTION AND ELEVATION########
brwn_weight_pc1 <- gls(PC1 ~ elevation_scaled,correlation=corBrownian(1,tree,form=~sp_abbr),data=sp_coord_anno_df)
summary(brwn_weight_pc1)
# Generalized least squares fit by REML
# Model: PC1 ~ elevation_scaled 
# Data: sp_coord_anno_df 
# AIC      BIC    logLik
# 259.7686 264.7593 -126.8843
# 
# Correlation Structure: corBrownian
# Formula: ~sp_abbr 
# Parameter estimate(s):
#   numeric(0)
# 
# Coefficients:
#   Value Std.Error   t-value p-value
# (Intercept)      0.04129179  3.639311 0.0113460  0.9910
# elevation_scaled 0.00057146  0.000588 0.9711564  0.3375
# 
# Correlation: 
#   (Intr)
# elevation_scaled -0.103
# 
# Standardized residuals:
#   Min           Q1          Med           Q3          Max 
# -1.646045293 -0.374659998 -0.006617844  0.287140417  0.919970964 
# 
# Residual standard error: 8.359464 
# Degrees of freedom: 41 total; 39 residual
