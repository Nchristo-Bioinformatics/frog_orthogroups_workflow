library(rjson)
library(tidyverse)
library(topGO)
library(RJSONIO)

# Write a function to read in the aBSREL JSON result file and transform it into a dataframe:
absrelJSONProcessing <- function(file, analysisType) {
  # Read in the JSON file:
  OG0000698_absrel <- RJSONIO::fromJSON(content = file)
  
  # Create a data frame with the information for each branch tested:
  OG0000698_absrelResult <- data.table::rbindlist(OG0000698_absrel[["branch attributes"]][["0"]], fill = TRUE, idcol = TRUE)
  #OG0000698_absrelResult <- subset(OG0000698_absrelResult, select = -c(`Rate Distributions`))
  OG0000698_absrelResult <- distinct (OG0000698_absrelResult)
  
  #Add a column to hold the orthogroup name:
  orthogroup <- sapply(strsplit(OG0000698_absrel[["input"]][["file name"]],"/"), tail, 1)
  orthogroup <- sapply(strsplit(orthogroup, "[.]"), head, 1)
  OG0000698_absrelResult$Orthogroup <- orthogroup
  OG0000698_absrelResult$"original name" <- sub("Scutiger_cf_sikimmensis","Scutigercf_sikimmensis",OG0000698_absrelResult$"original name")
  #species <- paste(unlist(sapply(strsplit(OG0000698_absrelResult[["original name"]],"_"), head, 2)),collapse="_")
  species <- sapply(strsplit(OG0000698_absrelResult[["original name"]],"_"), head, 2)
  #species <- as.factor(species)
  #species <- as.vector(species)
  #species <- paste(unlist(species),collapse="_")
  OG0000698_absrelResult$Species <- species
  #OG0000698_absrelResult$Species <- sub(",", "_", OG0000698_absrelResult$Species)
  # Add a column to label genes as having been part of the foreground or background
  OG0000698_absrelResult$ForegroundBackground <- analysisType
  
  return(OG0000698_absrelResult)
}
foregroundFiles <- list.files(path = "./foreground_absrel", pattern = "*.json", full.names = TRUE)
# Drop any files with file size of zero:
foregroundFiles <- foregroundFiles[sapply(foregroundFiles, file.size) > 0]

# Make a "Safe" version of the function that will return an error if there is a problem with a particular file:
possiblyabsrelJSONProcessing <- possibly(absrelJSONProcessing, otherwise = "File empty.")
# Map the function over all foreground files:
foregroundaBSRELResults <- map_dfr(foregroundFiles, analysisType = "foreground", possiblyabsrelJSONProcessing)
do <- as.data.frame(do.call(rbind, lapply(foregroundaBSRELResults$Species, as.vector)))
foregroundaBSRELResults <- cbind(Species_new=paste(do$V1,do$V2, sep='_'), foregroundaBSRELResults )
foregroundaBSRELResults <- foregroundaBSRELResults[,-c('Species')]
colnames(foregroundaBSRELResults)[1] <- "Species"
# Do the same for background files:
backgroundFiles <- list.files(path = "./background_absrel", pattern = "*.json", full.names = TRUE)
# Drop any files with file size of zero:
backgroundFiles <- backgroundFiles[sapply(backgroundFiles, file.size) > 0]

# Map the function over all background results:
backgroundaBSRELResults <- map_dfr(backgroundFiles, analysisType = "background", possiblyabsrelJSONProcessing)
do <- as.data.frame(do.call(rbind, lapply(backgroundaBSRELResults$Species, as.vector)))
backgroundaBSRELResults <- cbind(Species_new=paste(do$V1,do$V2, sep='_'), backgroundaBSRELResults )
backgroundaBSRELResults <- backgroundaBSRELResults[,-c('Species')]
colnames(backgroundaBSRELResults)[1] <- "Species"
allResults <- rbind(backgroundaBSRELResults, foregroundaBSRELResults,fill=TRUE)
allResults$branch <- as.factor(allResults$.id)
# Remove rows with NAs in the corrected p-value columns
allResults <- drop_na(allResults, `Corrected P-value`) ##DON"T CAREW ABOUT NODES###
allResults <- drop_na(allResults, "original name")
allResults <- drop_na(allResults, "Species")
allResults <- allResults %>%
  mutate(Rate = map_dbl(allResults$'Rate Distributions', 1),
         Percent = map_dbl(allResults$'Rate Distributions', 2))
allResults$Percent_transformed <- ifelse(allResults$'Rate classes' ==1, 0, allResults$'Percent')
allResults$Percent_transformed <- ifelse(allResults$'Corrected P-value' > 0.05, 0, allResults$Percent_transformed)
allResults$Significant <- ifelse(allResults$'Corrected P-value' > 0.05, 0, allResults$Percent_transformed)	
allResults$Significant <- ifelse(allResults$'Corrected P-value' <= 0.05, 1, allResults$Significant)	
allResults$omega_transformed <- ifelse(allResults$'Rate classes' ==1, 1, allResults$'Rate')
#allResults$omega_transformed <- ifelse(allResults$'Rate classes' ==1, 1, allResults$omega_transformed)
allResults$omega_transformed <- ifelse(allResults$'Corrected P-value' > 0.05, 1, allResults$omega_transformed)
allResults$omega_transformed <- ifelse(allResults$omega_transformed >= 10000, 10000, allResults$omega_transformed)
rep_na_mean <- function(vec){
  m <- mean(vec, na.rm = TRUE) 
  vec[is.na(vec)] <- m 
  return(vec) 
}
percentages <- allResults %>%
  group_by(branch) %>%
  filter(Rate %in% c(max(Rate)))
percentages <- ungroup(percentages)
percentages <- percentages[!duplicated(percentages[,c('branch')]),]
bsrel_sp_weight <- percentages %>%
  dplyr::select(Species,Orthogroup,Percent_transformed) %>%
  spread(Orthogroup,Percent_transformed) %>%
  mutate_at(vars(-Species),rep_na_mean) 

omgs <- allResults
omgs$avg_rate <- as.numeric(omgs$omega_transformed)*as.numeric(omgs$Percent_transformed)
omgs <- omgs %>% group_by(branch) %>% mutate(total_rate = sum(avg_rate)) %>% filter(Rate %in% c(max(Rate)))
omgs <- ungroup(omgs)
omgs <- omgs[!duplicated(omgs[,c('branch')]),]
bsrel_sp_omega_theirway <- omgs %>%
  dplyr::select(Species,Orthogroup,omega_transformed) %>%
  #mutate(omega_sig=ifelse(omega_sig==0,1e-18,omega_sig)) %>%
  spread(Orthogroup,omega_transformed) %>%
  mutate_at(vars(-Species),rep_na_mean)
log_use_0s_omega <- function(column){
  column[column == 0] <- 1
  return(log(column))
}
bsrel_sp_omega_theirway <- bsrel_sp_omega_theirway %>%
  mutate_at(vars(-Species),log_use_0s_omega)


bsrel_sp_omega_myway <- omgs %>%
  dplyr::select(Species,Orthogroup,total_rate) %>%
  #mutate(omega_sig=ifelse(omega_sig==0,1e-18,omega_sig)) %>%
  spread(Orthogroup,total_rate) %>%
  mutate_at(vars(-Species),rep_na_mean)
bsrel_sp_omega_myway <- bsrel_sp_omega_myway %>%
  mutate_at(vars(-Species),log_use_0s_omega)

###########PREPARE FOR PCA###############
OGs_to_keep <- colnames(bsrel_sp_omega_theirway)
bsrel_sp_weight_forpca <- bsrel_sp_weight %>%
  dplyr::select(Species,OGs_to_keep) %>%
  as.data.frame
rownames(bsrel_sp_weight_forpca) <- bsrel_sp_weight_forpca$Species
bsrel_sp_weight_forpca$Species <- NULL
#write.csv(bsrel_sp_weight_forpca, 'percentweight_transformed_values_pca.csv')
#bsrel_sp_weight_forpca <- read.csv('percentweight_transformed_values_pca.csv')
#rownames(bsrel_sp_weight_forpca) <- bsrel_sp_weight_forpca$X
#bsrel_sp_weight_forpca <- subset(bsrel_sp_weight_forpca, select = -c(`X`))
bsrel_sp_weight_forpca <- bsrel_sp_weight_forpca[,which(colSums(bsrel_sp_weight_forpca) != 0)]
bsrel_sp_omega_theirway_forpca <- bsrel_sp_omega_theirway %>%
  dplyr::select(Species,OGs_to_keep) %>%
  as.data.frame
#write.csv(bsrel_sp_omega_theirway_forpca, 'omega_transformed_values_pca.csv')
#bsrel_sp_omega_theirway_forpca <- read.csv('omega_transformed_values_pca.csv')
rownames(bsrel_sp_omega_theirway_forpca) <- bsrel_sp_omega_theirway_forpca$Species
bsrel_sp_omega_theirway_forpca$Species <- NULL
bsrel_sp_omega_theirway_forpca <- bsrel_sp_omega_theirway_forpca[,which(colSums(bsrel_sp_omega_theirway_forpca) != 0)]


bsrel_sp_omega_myway_forpca <- bsrel_sp_omega_myway %>%
  dplyr::select(Species,OGs_to_keep) %>%
  as.data.frame
#write.csv(bsrel_sp_omega_myway_forpca, 'omega_transformed_values_pca.csv')
#bsrel_sp_omega_myway_forpca <- read.csv('omega_transformed_values_pca.csv')
rownames(bsrel_sp_omega_myway_forpca) <- bsrel_sp_omega_myway_forpca$Species
bsrel_sp_omega_myway_forpca$Species <- NULL
bsrel_sp_omega_myway_forpca <- bsrel_sp_omega_myway_forpca[,which(colSums(bsrel_sp_omega_myway_forpca) != 0)]

pca_omega <- prcomp(bsrel_sp_omega_theirway_forpca)
pca_omega_myway <- prcomp(bsrel_sp_omega_myway_forpca)
pca_weights <- prcomp(bsrel_sp_weight_forpca)
#pca_kvalue <- prcomp(expr_ks)

summary(pca_omega)
summary(pca_weights)
summary(pca_omega_myway)

#######PRISTIMANTIS########
#Save omega results
write_csv(data.frame(pca_omega$rotation),"pca_sp_omega_loadings_pristi.csv")
write_csv(data.frame(pca_omega$x),"pca_sp_omega_coordinates_pristi.csv")

save(pca_omega,bsrel_sp_omega_theirway, bsrel_sp_omega_theirway_forpca,file="pca_sp_omega_all_res_pristi.rDat")
write_csv(data.frame(pca_omega_myway$rotation),"pca_sp_omega_loadings_pristi_myway.csv")
write_csv(data.frame(pca_omega_myway$x),"pca_sp_omega_coordinates_pristi_myway.csv")
# 
save(pca_omega_myway,bsrel_sp_omega_myway, bsrel_sp_omega_myway_forpca,file="pca_sp_omega_all_res_pristi_myway.rDat")
# 
#Save weights results
write_csv(data.frame(pca_weights$rotation),"pca_sp_weights_loadings_pristi.csv")
write_csv(data.frame(pca_weights$x),"pca_sp_weights_coordinates_pristi.csv")

save(pca_weights,bsrel_sp_weight,bsrel_sp_weight_forpca, file="pca_sp_weighs_all_res_pristi.rDat")

#####DIRECT DEVS######
#Save omega results
# write_csv(data.frame(pca_omega$rotation),"pca_sp_omega_loadings_dirdev.csv")
# write_csv(data.frame(pca_omega$x),"pca_sp_omega_coordinates_dirdev.csv")
# 
# save(pca_omega,bsrel_sp_omega_theirway, bsrel_sp_omega_theirway_forpca,file="pca_sp_omega_all_res_dirdev.rDat")
# write_csv(data.frame(pca_omega_myway$rotation),"pca_sp_omega_loadings_dirdev_myway.csv")
# write_csv(data.frame(pca_omega_myway$x),"pca_sp_omega_coordinates_dirdev_myway.csv")
# 
# save(pca_omega_myway,bsrel_sp_omega_myway, bsrel_sp_omega_myway_forpca,file="pca_sp_omega_all_res_dirdev_myway.rDat")
# 
# #Save weights results
# write_csv(data.frame(pca_weights$rotation),"pca_sp_weights_loadings_dirdev.csv")
# write_csv(data.frame(pca_weights$x),"pca_sp_weights_coordinates_dirdev.csv")

#save(pca_weights,bsrel_sp_weight,bsrel_sp_weight_forpca, file="pca_sp_weighs_all_res_dirdev.rDat")

#allResults2 <- apply(allResults,2,as.character)
#write.csv(allResults2, 'all_absrel_results_2ndtime.csv')
############################################START HERE FROM RSTUDIO AT HOME##########################
library(ape)
library(phytools)
library(geiger)
library(ggplot2)
setwd('~/Orthogroup_paper_edits_2022/TO_DO_FOR_PAPER/simple_results_for_table_processing_r/pristimantis')
load("pca_sp_omega_all_res_pristi.rDat")
load("pca_sp_weighs_all_res_pristi.rDat")
load("pca_sp_omega_all_res_pristi_myway.rDat")

# setwd('~/Orthogroup_paper_edits_2022/TO_DO_FOR_PAPER/simple_results_for_table_processing_r/direct_developers')
# load("pca_sp_omega_all_res_dirdev.rDat")
# load("pca_sp_weighs_all_res_dirdev.rDat")
# load("pca_sp_omega_all_res_dirdev_myway.rDat")

####################################################################################################
########DO STUFF WITH OMEGA, % SELECTED#############
####################################################################################################
#pca_sp <- pca_omega
pca_sp <- pca_weights
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

######ANOVAS, JUST USING PRISTIMANTIS VALUES SINCE THEY ARE ABOUT THE SAME##########

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

###### Are the omega values different for foreground vs background species? ###PHYLANOVA

model_omega_pristi <- phylANOVA(tree, grp_pristi, sp_names_PC1, nsim=1000, posthoc=TRUE, p.adj="BH")
model_omega_pristi
# ANOVA table: Phylogenetic ANOVA
# 
# Response: y
# Sum Sq    Mean Sq   F value Pr(>F)
# x        1300.460 1300.46036 14.149904  0.078
# Residual 3584.332   91.90595                 
# 
# P-value based on simulation.
# ---------
#   
#   Pairwise posthoc test using method = "BH"
# 
# Pairwise t-values:
#   Other Pristimantis
# Other        0.000000    -3.761636
# Pristimantis 3.761636     0.000000
# 
# Pairwise corrected P-values:
#   Other Pristimantis
# Other        1.000        0.078
# Pristimantis 0.078        1.000


model_omega_dirdev <- phylANOVA(tree, grp, sp_names_PC1, nsim=1000, posthoc=TRUE, p.adj="BH")
model_omega_dirdev
# ANOVA table: Phylogenetic ANOVA
# 
# Response: y
# Sum Sq   Mean Sq  F value Pr(>F)
# x         995.3636 995.36357 9.980689  0.144
# Residual 3889.4289  99.72894                
# 
# P-value based on simulation.
# ---------
#   
#   Pairwise posthoc test using method = "BH"
# 
# Pairwise t-values:
#   Direct_developer    Other
# Direct_developer         0.000000 3.159223
# Other                   -3.159223 0.000000
# 
# Pairwise corrected P-values:
#   Direct_developer Other
# Direct_developer            1.000 0.144
# Other                       0.144 1.000

###### Are the proportion of positivelys selected sites different for foreground vs background species? ###PHYLANOVA
model_weights_pristi <- phylANOVA(tree, grp_pristi, sp_names_PC1, nsim=1000, posthoc=TRUE, p.adj="BH")
model_weights_pristi
# ANOVA table: Phylogenetic ANOVA
# 
# Response: y
# Sum Sq  Mean Sq  F value Pr(>F)
# x        0.001112 0.001112 0.063839  0.924
# Residual 0.679229 0.017416                
# 
# P-value based on simulation.
# ---------
#   
#   Pairwise posthoc test using method = "BH"
# 
# Pairwise t-values:
#   Other Pristimantis
# Other        0.000000    -0.252664
# Pristimantis 0.252664     0.000000
# 
# Pairwise corrected P-values:
#   Other Pristimantis
# Other        1.000        0.924
# Pristimantis 0.924        1.000
model_weights_dirdev <- phylANOVA(tree, grp, sp_names_PC1, nsim=1000, posthoc=TRUE, p.adj="BH")
model_weights_dirdev
# ANOVA table: Phylogenetic ANOVA
# 
# Response: y
# Sum Sq  Mean Sq  F value Pr(>F)
# x        0.001325 0.001325 0.076098  0.928
# Residual 0.679016 0.017411                
# 
# P-value based on simulation.
# ---------
#   
#   Pairwise posthoc test using method = "BH"
# 
# Pairwise t-values:
#   Direct_developer    Other
# Direct_developer         0.000000 0.275858
# Other                   -0.275858 0.000000
# 
# Pairwise corrected P-values:
#   Direct_developer Other
# Direct_developer            1.000 0.928
# Other                       0.928 1.000
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
elevation_scaled <- setNames(sp_coord_anno_df[,"elevation"],rownames(sp_coord_anno_df))
PC1_scaled <- setNames(sp_coord_anno_df[,"PC1"],rownames(sp_coord_anno_df))

##############Plot phylomorphospace plot of elevation by % selected PC1
pdf("PC1_sp_weight_log_elevation_phylomorphospace.pdf",width=8,height=8)
phylomorphospace(tree,cbind(elevation_scaled ,PC1_scaled),node.size=c(0,1),label="off",xlab="Elevation", ylab="Weight PC1")
points(elevation_scaled,PC1_scaled,pch=16,col=species_info_df$Colour)
legend("topleft", legend=c("Pristimantis"), col=c("red"), pch=19)
title(main="Phylomorphospace of elevation and PC proportion values")
dev.off()
########CORELATE % SELECTION AND ELEVATION########

brwn_weight_pc1 <- gls(PC1 ~ elevation_scaled,correlation=corBrownian(1,tree,form=~sp_abbr),data=sp_coord_anno_df)
summary(brwn_weight_pc1)
# Generalized least squares fit by REML
# Model: PC1 ~ elevation_scaled 
# Data: sp_coord_anno_df 
# AIC       BIC  logLik
# -24.1132 -19.12252 15.0566
# 
# Correlation Structure: corBrownian
# Formula: ~sp_abbr 
# Parameter estimate(s):
#   numeric(0)
# 
# Coefficients:
#   Value  Std.Error    t-value p-value
# (Intercept)      -0.007705809 0.09558722 -0.0806155  0.9362
# elevation_scaled  0.000019208 0.00001546  1.2427818  0.2214
# 
# Correlation: 
#   (Intr)
# elevation_scaled -0.103
# 
# Standardized residuals:
#   Min          Q1         Med          Q3         Max 
# -2.49013300 -0.05913182  0.12917186  0.20088968  0.56069554 
# 
# Residual standard error: 0.219563 
# Degrees of freedom: 41 total; 39 residual

##############Plot phylomorphospace plot of elevation by omega PC1
pdf("PC1_sp_omega_elevation_phylomorphospace.pdf",width=8,height=8)
phylomorphospace(tree,cbind(elevation_scaled ,PC1_scaled),node.size=c(0,1),label="off",xlab="Elevation", ylab="Omega PC1")
points(elevation_scaled,PC1_scaled,pch=16,col=species_info_df$Colour)
legend("topleft", legend=c("Pristimantis"), col=c("red"), pch=19)
title(main="Phylomorphospace of elevation and PC1 omega values")
dev.off()
########CORELATE Omega AND ELEVATION########
brwn_weight_pc1 <- gls(PC1 ~ elevation_scaled,correlation=corBrownian(1,tree,form=~sp_abbr),data=sp_coord_anno_df)
summary(brwn_weight_pc1)
# Generalized least squares fit by REML
# Model: PC1 ~ elevation_scaled 
# Data: sp_coord_anno_df 
# AIC      BIC    logLik
# 332.1737 337.1644 -163.0869
# 
# Correlation Structure: corBrownian
# Formula: ~sp_abbr 
# Parameter estimate(s):
#   numeric(0)
# 
# Coefficients:
#   Value Std.Error    t-value p-value
# (Intercept)      -4.263785  9.207930 -0.4630557  0.6459
# elevation_scaled  0.000219  0.001489  0.1472029  0.8837
# 
# Correlation: 
#   (Intr)
# elevation_scaled -0.103
# 
# Standardized residuals:
#   Min          Q1         Med          Q3         Max 
# -0.29535184 -0.12505801 -0.02020213  0.23249494  1.98116218 
# 
# Residual standard error: 21.15053 
# Degrees of freedom: 41 total; 39 residual