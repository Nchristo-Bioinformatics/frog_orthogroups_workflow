###PREPARING YOUR ABSREL RESULTS###
library(RJSONIO)
library(tidyverse)
library(topGO)

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
foregroundFiles <- list.files(path = "./foreground", pattern = "*.json", full.names = TRUE)
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
backgroundFiles <- list.files(path = "./background", pattern = "*.json", full.names = TRUE)
# Drop any files with file size of zero:
backgroundFiles <- backgroundFiles[sapply(backgroundFiles, file.size) > 0]

# Map the function over all background results:
backgroundaBSRELResults <- map_dfr(backgroundFiles, analysisType = "background", possiblyabsrelJSONProcessing)
do <- as.data.frame(do.call(rbind, lapply(backgroundaBSRELResults$Species, as.vector)))
backgroundaBSRELResults <- cbind(Species_new=paste(do$V1,do$V2, sep='_'), backgroundaBSRELResults )
backgroundaBSRELResults <- backgroundaBSRELResults[,-c('Species')]
colnames(backgroundaBSRELResults)[1] <- "Species"


# Combine the foreground and background tables:
allResults <- rbind(backgroundaBSRELResults, foregroundaBSRELResults,fill=TRUE)
allResults$branch <- as.factor(allResults$.id)
# Remove rows with NAs in the corrected p-value columns
allResults <- drop_na(allResults, `Corrected P-value`) ##DON"T CAREW ABOUT NODES###
allResults <- drop_na(allResults, "original name")
allResults <- drop_na(allResults, "Species")

########TRANSFORMMMM##############
allResults <- allResults %>%
  mutate(Rate = map_dbl(allResults$'Rate Distributions', 1),
         Percent = map_dbl(allResults$'Rate Distributions', 2))
allResults$Percent_transformed <- ifelse(allResults$'Rate classes' ==1, 0, allResults$'Percent')
allResults$Percent_transformed <- ifelse(allResults$'Corrected P-value' >= 0.05, 0, allResults$Percent_transformed)
allResults$Significant <- ifelse(allResults$'Corrected P-value' >= 0.05, 0, allResults$Percent_transformed)	
allResults$Significant <- ifelse(allResults$'Corrected P-value' <= 0.05, 1, allResults$Significant)	
allResults$omega_transformed <- ifelse(allResults$'Rate classes' ==1, 0, allResults$'Rate')
#allResults$omega_transformed <- ifelse(allResults$'Rate classes' ==1, 1, allResults$omega_transformed)
allResults$omega_transformed <- ifelse(allResults$'Corrected P value' >= 0.05, 0, allResults$omega_transformed)
rep_na_mean <- function(vec){
  m <- mean(vec, na.rm = TRUE) 
  vec[is.na(vec)] <- m 
  return(vec) 
}
percentages <- allResults %>%
	group_by(branch) %>%
	filter(Rate %in% c(max(Rate)))
percentages <- ungroup(percentages)
bsrel_sp_weight <- percentages %>%
  dplyr::select(Species,Orthogroup,Percent_transformed) %>%
  spread(Orthogroup,Percent_transformed) %>%
  mutate_at(vars(-Species),rep_na_mean) 
  
bsrel_sp_omega <- percentages %>%
  dplyr::select(Species,Orthogroup,omega_transformed) %>%
  #mutate(omega_sig=ifelse(omega_sig==0,1e-18,omega_sig)) %>%
  spread(Orthogroup,omega_transformed) %>%
  mutate_at(vars(-Species),rep_na_mean)
log_use_0s_omega <- function(column){
  column[column == 0] <- 1
  return(log(column))
}
bsrel_sp_omega_logtrans <- bsrel_sp_omega %>%
  mutate_at(vars(-Species),log_use_0s_omega)
  
  ###########PREPARE FOR PCA###############
OGs_to_keep <- colnames(bsrel_sp_omega_logtrans)
bsrel_sp_weight_forpca <- bsrel_sp_weight %>%
  dplyr::select(Species,OGs_to_keep) %>%
  as.data.frame
rownames(bsrel_sp_weight_forpca) <- bsrel_sp_weight_forpca$Species
bsrel_sp_weight_forpca$Species <- NULL
write.csv(bsrel_sp_weight_forpca, 'percentweight_transformed_values_pca.csv')
#bsrel_sp_weight_forpca <- read.csv('percentweight_transformed_values_pca.csv')
rownames(bsrel_sp_weight_forpca) <- bsrel_sp_weight_forpca$X
bsrel_sp_weight_forpca <- subset(bsrel_sp_weight_forpca, select = -c(`X`))
bsrel_sp_weight_forpca <- bsrel_sp_weight_forpca[,which(colSums(bsrel_sp_weight_forpca) != 0)]
bsrel_sp_omega_logtrans_forpca <- bsrel_sp_omega_logtrans %>%
  dplyr::select(Species,OGs_to_keep) %>%
  as.data.frame
write.csv(bsrel_sp_omega_logtrans_forpca, 'omega_transformed_values_pca.csv')
#bsrel_sp_omega_logtrans_forpca <- read.csv('omega_transformed_values_pca.csv')
rownames(bsrel_sp_omega_logtrans_forpca) <- bsrel_sp_omega_logtrans_forpca$X
bsrel_sp_omega_logtrans_forpca <- subset(bsrel_sp_omega_logtrans_forpca, select = -c(`X`))
bsrel_sp_omega_logtrans_forpca <- bsrel_sp_omega_logtrans_forpca[,which(colSums(bsrel_sp_omega_logtrans_forpca) != 0)]
#bsrel_sp_omega_logtrans_forpca$Species <- NULL
pca_omega <- prcomp(bsrel_sp_omega_logtrans_forpca)
pca_weights <- prcomp(bsrel_sp_weight_forpca)
pca_kvalue <- prcomp(expr_ks)

loading_weights <- pca_weights$rotation %>%
  as.data.frame %>%
  rownames_to_column("Orthogroup") %>%
  as.tibble

loading_omega <- pca_omega$rotation %>%
  as.data.frame %>%
  rownames_to_column("Orthogroups") %>%
  as.tibble

#Scree plot

pdf("pca_sp_omega_scree_plot.pdf")
plot(pca_omega)
dev.off()

pdf("pca_sp_weights_scree_plot.pdf")
plot(pca_weights)
dev.off()

#Variance
summary(pca_pvals)
summary(pca_omega)
summary(pca_weights)


#Save omega results
write_csv(data.frame(pca_omega$rotation),"pca_sp_omega_loadings.csv")
write_csv(data.frame(pca_omega$x),"pca_sp_omega_coordinates.csv")

save(pca_omega,bsrel_sp_omega,bsrel_sp_omega_logtrans, bsrel_sp_omega_logtrans_forpca,file="pca_sp_omega_all_res.rDat")

#Save weights results
write_csv(data.frame(pca_weights$rotation),"pca_sp_weights_loadings.csv")
write_csv(data.frame(pca_weights$x),"pca_sp_weights_coordinates.csv")

save(pca_weights,bsrel_sp_weight,bsrel_sp_weight_forpca, file="pca_sp_weighs_all_res.rDat")


allResults2 <- apply(allResults,2,as.character)
write.csv(allResults2, 'all_absrel_results_2ndtime.csv')
############################################START HERE FROM RSTUDIO AT HOME##########################
load("pca_sp_omega_all_res.rDat")
load("pca_sp_weighs_all_res.rDat")
load("pca_sp_kvalue_all_res.rDat")
####################################################################################################
########DO STUFF WITH OMEGA#############
####################################################################################################
pca_sp <- pca_omega
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

pdf("PC_individual_PC1_sp_omega.pdf")
plotBranchbyTrait(tree,sp_names_PC1,mode="tips",palette = "custom",use.edge.length=FALSE)
dev.off()
pdf("PC_individual_PC2_sp_omega.pdf")
plotBranchbyTrait(tree,sp_names_PC2,mode="tips",palette = "custom",use.edge.length=FALSE)
dev.off()
pdf("PC_individual_PC3_sp_omega.pdf")
plotBranchbyTrait(tree,sp_names_PC3,mode="tips",palette = "custom",use.edge.length=FALSE)
dev.off()

#########Correlations with PC axes and life history traits
species_info_df$Colour <- as.factor(species_info_df$Colour)
cols<-setNames(c("blue","red"),
    levels(sig_regressions_df$Colour))
species_info_df$elevation <- ifelse(species_info_df$elevation == 0 , 1, species_info_df$elevation)
sp_coord <- pca_sp$x %>%
  as.data.frame %>%
  rownames_to_column("sp_abbr") %>%
  as.tibble
sp_coord$sp_abbr <- species_info_df$Species
sp_coord_anno <- sp_coord %>%
  left_join(species_info_df,by=c("sp_abbr" = "Species")) %>%
  mutate(log_elevation = log(elevation))

sp_coord_anno %>%
  ggplot(aes(log(elevation),PC1)) +
  geom_point() +
  xlab("log(elevation)")
ggsave("PC1_sp_omega_log_elevation.pdf",width=6,height=4)


sp_coord_anno_df <- sp_coord_anno %>%
  as.data.frame()
rownames(sp_coord_anno_df) <- sp_coord_anno_df$sp_abbr
elevation_scaled <- setNames(sp_coord_anno_df[,"log_elevation"],rownames(sp_coord_anno_df))
PC1_scaled <- setNames(sp_coord_anno_df[,"PC1"],rownames(sp_coord_anno_df))
#Plot phylomorphospace plot of body size by PC1
pdf("PC1_sp_omega_log_elevation_phylomorphospace.pdf",width=8,height=8)
phylomorphospace(tree,cbind(elevation_scaled ,PC1_scaled),node.size=c(0,1),label="off",xlab="log(elevation)", ylab="Omega PC1")
points(elevation_scaled,PC1_scaled,pch=16,col=species_info_df$Colour)
legend("topleft", legend=c("Pristimantis"), col=c("red"), pch=19)
title(main="Phylomorphospace of elevation and PC omega values")
dev.off()


########CORELATE OMEGA AND ELEVATION########
brwn_omega_pc1 <- gls(PC1 ~ elevation_scaled,correlation=corBrownian(1,tree,form=~sp_abbr),data=sp_coord_anno_df)
summary(brwn_omega_pc1)
# Generalized least squares fit by REML
  # Model: PC1 ~ elevation_scaled 
  # Data: sp_coord_anno_df 
       # AIC      BIC    logLik
  # 331.6771 336.6678 -162.8385

# Correlation Structure: corBrownian
 # Formula: ~sp_abbr 
 # Parameter estimate(s):
# numeric(0)

# Coefficients:
                     # Value Std.Error    t-value p-value
# (Intercept)      -9.936018 12.288313 -0.8085746  0.4237
# elevation_scaled  0.925157  1.047971  0.8828083  0.3827

 # Correlation: 
                 # (Intr)
# elevation_scaled -0.488

# Standardized residuals:
        # Min          Q1         Med          Q3         Max 
# -0.31951024 -0.14846401 -0.06780641  0.28843210  2.03274073 

####################################################################################################
########DO STUFF WITH WEIGHT (% selected) #############
####################################################################################################
pca_sp <- pca_weights
tree <- read.tree('SpeciesTree_rooted.txt')



species_info_df <- read.csv('absrel_counts_both_hypotheses.csv', row.names=1, header=TRUE)
species_info_df$Species <- rownames(species_info_df)
sp_names_PC1 <- pca_sp$x[,"PC1"]
names(sp_names_PC1) <- rownames(species_info_df)
sp_names_PC2 <- pca_sp$x[,"PC2"]
names(sp_names_PC2) <- rownames(species_info_df)
sp_names_PC3 <- pca_sp$x[,"PC3"]
names(sp_names_PC3) <- rownames(species_info_df)
plotBranchbyTrait(tree,sp_names_PC1,mode="tips",palette = "custom",use.edge.length=FALSE)

pdf("PC_individual_PC1_sp_weights.pdf")
plotBranchbyTrait(tree,sp_names_PC1,mode="tips",palette = "custom",use.edge.length=FALSE)
dev.off()
pdf("PC_individual_PC2_sp_weights.pdf")
plotBranchbyTrait(tree,sp_names_PC2,mode="tips",palette = "custom",use.edge.length=FALSE)
dev.off()
pdf("PC_individual_PC3_sp_weights.pdf")
plotBranchbyTrait(tree,sp_names_PC3,mode="tips",palette = "custom",use.edge.length=FALSE)
dev.off()

#########Correlations with PC axes and life history traits
species_info_df$Colour <- as.factor(species_info_df$Colour)
cols<-setNames(c("blue","red"),
    levels(sig_regressions_df$Colour))
species_info_df$elevation <- ifelse(species_info_df$elevation == 0 , 1, species_info_df$elevation)
sp_coord <- pca_sp$x %>%
  as.data.frame %>%
  rownames_to_column("sp_abbr") %>%
  as.tibble
sp_coord$sp_abbr <- species_info_df$Species
sp_coord_anno <- sp_coord %>%
  left_join(species_info_df,by=c("sp_abbr" = "Species")) %>%
  mutate(log_elevation = log(elevation))

sp_coord_anno %>%
  ggplot(aes(log(elevation),PC1)) +
  geom_point() +
  xlab("log(elevation)")
ggsave("PC1_sp_omega_log_elevation.pdf",width=6,height=4)


sp_coord_anno_df <- sp_coord_anno %>%
  as.data.frame()
rownames(sp_coord_anno_df) <- sp_coord_anno_df$sp_abbr
elevation_scaled <- setNames(sp_coord_anno_df[,"log_elevation"],rownames(sp_coord_anno_df))
PC1_scaled <- setNames(sp_coord_anno_df[,"PC1"],rownames(sp_coord_anno_df))
#Plot phylomorphospace plot of body size by PC1
pdf("PC1_sp_weight_log_elevation_phylomorphospace.pdf",width=8,height=8)
phylomorphospace(tree,cbind(elevation_scaled ,PC1_scaled),node.size=c(0,1),label="off",xlab="log(elevation)", ylab="Weight PC1")
points(elevation_scaled,PC1_scaled,pch=16,col=species_info_df$Colour)
legend("topleft", legend=c("Pristimantis"), col=c("red"), pch=19)
title(main="Phylomorphospace of elevation and PC proportion values")
dev.off()


########CORELATE % SELECTION AND ELEVATION########
brwn_weight_pc1 <- gls(PC1 ~ elevation_scaled,correlation=corBrownian(1,tree,form=~sp_abbr),data=sp_coord_anno_df)
summary(brwn_weight_pc1)
#Generalized least squares fit by REML
  Model: PC1 ~ elevation_scaled 
  Data: sp_coord_anno_df 
       AIC       BIC   logLik
  -39.2287 -34.23802 22.61435

Correlation Structure: corBrownian
 Formula: ~sp_abbr 
 Parameter estimate(s):
numeric(0)

Coefficients:
                       Value  Std.Error   t-value p-value
(Intercept)       0.11043842 0.10576306  1.044206  0.3028
elevation_scaled -0.01961975 0.00901968 -2.175217  0.0357

 Correlation: 
                 (Intr)
elevation_scaled -0.488

Standardized residuals:
       Min         Q1        Med         Q3        Max 
-0.8777550 -0.2209063 -0.1128560  0.1501799  2.2375832 

Residual standard error: 0.2131349 
Degrees of freedom: 41 total; 39 residual
####################################################################################################
########DO STUFF WITH K VALUES (RELAX) #############
####################################################################################################
pca_sp <- pca_kvalue
tree <- read.tree('SpeciesTree_rooted.txt')

species_info_df <- read.csv('absrel_counts_both_hypotheses.csv', row.names=1, header=TRUE)
species_info_df$Species <- rownames(species_info_df)
sp_names_PC1 <- pca_sp$x[,"PC1"]
names(sp_names_PC1) <- rownames(species_info_df)
sp_names_PC2 <- pca_sp$x[,"PC2"]
names(sp_names_PC2) <- rownames(species_info_df)
sp_names_PC3 <- pca_sp$x[,"PC3"]
names(sp_names_PC3) <- rownames(species_info_df)
plotBranchbyTrait(tree,sp_names_PC1,mode="tips",palette = "custom",use.edge.length=FALSE)

pdf("PC_individual_PC1_sp_kvalue.pdf")
plotBranchbyTrait(tree,sp_names_PC1,mode="tips",palette = "custom",use.edge.length=FALSE)
dev.off()
pdf("PC_individual_PC2_sp_kvalue.pdf")
plotBranchbyTrait(tree,sp_names_PC2,mode="tips",palette = "custom",use.edge.length=FALSE)
dev.off()
pdf("PC_individual_PC3_sp_kvalue.pdf")
plotBranchbyTrait(tree,sp_names_PC3,mode="tips",palette = "custom",use.edge.length=FALSE)
dev.off()

#########Correlations with PC axes and life history traits
species_info_df$Colour <- as.factor(species_info_df$Colour)
cols<-setNames(c("blue","red"),
    levels(sig_regressions_df$Colour))
species_info_df$elevation <- ifelse(species_info_df$elevation == 0 , 1, species_info_df$elevation)
sp_coord <- pca_sp$x %>%
  as.data.frame %>%
  rownames_to_column("sp_abbr") %>%
  as.tibble
sp_coord$sp_abbr <- species_info_df$Species
sp_coord_anno <- sp_coord %>%
  left_join(species_info_df,by=c("sp_abbr" = "Species")) %>%
  mutate(log_elevation = log(elevation))

sp_coord_anno %>%
  ggplot(aes(log(elevation),PC1)) +
  geom_point() +
  xlab("log(elevation)")
ggsave("PC1_sp_kvalue_log_elevation.pdf",width=6,height=4)


sp_coord_anno_df <- sp_coord_anno %>%
  as.data.frame()
rownames(sp_coord_anno_df) <- sp_coord_anno_df$sp_abbr
elevation_scaled <- setNames(sp_coord_anno_df[,"log_elevation"],rownames(sp_coord_anno_df))
PC1_scaled <- setNames(sp_coord_anno_df[,"PC1"],rownames(sp_coord_anno_df))
#Plot phylomorphospace plot of body size by PC1
pdf("PC1_sp_kvalue_log_elevation_phylomorphospace.pdf",width=8,height=8)
phylomorphospace(tree,cbind(elevation_scaled ,PC1_scaled),node.size=c(0,1),label="off",xlab="log(body mass)", ylab="k relaxation value PC1")
points(elevation_scaled,PC1_scaled,pch=16,col=species_info_df$Colour)
legend("topleft", legend=c("Pristimantis"), col=c("red"), pch=19)
title(main="Phylomorphospace of elevation and PC k relaxation values")
dev.off()


########CORELATE kvalue AND ELEVATION########
brwn_kvalue_pc1 <- gls(PC1 ~ elevation_scaled,correlation=corBrownian(1,tree,form=~sp_abbr),data=sp_coord_anno_df)
summary(brwn_kvalue_pc1)
# Generalized least squares fit by REML
  # Model: PC1 ~ elevation_scaled 
  # Data: sp_coord_anno_df 
       # AIC      BIC    logLik
  # 247.8218 252.8125 -120.9109

# Correlation Structure: corBrownian
 # Formula: ~sp_abbr 
 # Parameter estimate(s):
# numeric(0)

# Coefficients:
                      # Value Std.Error     t-value p-value
# (Intercept)       0.5466762  4.193691  0.13035679  0.8970
# elevation_scaled -0.0322688  0.357646 -0.09022563  0.9286

 # Correlation: 
                 # (Intr)
# elevation_scaled -0.488

# Standardized residuals:
        # Min          Q1         Med          Q3         Max 
# -1.43514174 -0.34183376  0.01969883  0.38688907  1.15970617 

# Residual standard error: 8.451174 
# Degrees of freedom: 41 total; 39 residual

####################MAKE DEM HEATMAPS############################
library(ComplexHeatmap)
library(circlize)
SAMPLE <- read.csv('absrel_counts_both_hypotheses.csv', row.names=1)
SAMPLE$Clade <- SAMPLE$Pristi
SAMPLE <- subset(SAMPLE, select = c(Clade))
SAMPLE
colorder=c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41)
colerz <- list(Clade= c("Other"="coral2","Pristimantis"="#66cc33"))
expr_1 <- t(bsrel_sp_weight_forpca)
p1 <- Heatmap(expr_1, name = "% Selected sites", col = colorRamp2(c(0,.15,.30,.5), c("darkblue","lightblue","green", "yellow")),
  show_row_dend = FALSE,top_annotation = HeatmapAnnotation(df = SAMPLE,show_annotation_name = TRUE, annotation_name_side = "right", col=colerz),column_order=colorder,
 show_row_names = FALSE, show_column_names = TRUE, heatmap_legend_param = list(legend_width = unit(3, "mm")), column_names_gp=gpar(fontsize=5), show_heatmap_legend = TRUE, width=unit(6.5, "cm"))



expr_2 <- t(bsrel_sp_omega_logtrans_forpca)
endo = read.csv("GO_endopeptidase.csv", row.names=1)
endo1 = rownames(expr_2) %in% rownames(endo)
atp <- read.csv("GO_ATP.csv", row.names=1)
atp1 <- rownames(expr_2) %in% rownames(atp)
heme <- read.csv("GO_heme.csv", row.names=1)
heme1 <- rownames(expr_2) %in% rownames(heme)
oxido <- read.csv("GO_oxidoreductase.csv", row.names=1)
oxido1 <- rownames(expr_2) %in% rownames(oxido)
prot <- read.csv("GO_protein.csv", row.names=1)
prot1 <- rownames(expr_2) %in% rownames(prot)
translate <- read.csv("GO_translation.csv", row.names=1)
translate1 <- rownames(expr_2) %in% rownames(translate)
p2 <- Heatmap(expr_2, name = "log(omega)", col = colorRamp2(c(-2, 0, 5, 15), c("dodgerblue2", "white", "gold", "coral")), 
show_row_dend = FALSE,top_annotation = HeatmapAnnotation(df = SAMPLE,show_annotation_name = FALSE, col=colerz), column_order=colorder, show_row_names = FALSE, show_column_names = TRUE, heatmap_legend_param = list(legend_width = unit(3, "mm")),
column_names_gp=gpar(fontsize=5), show_heatmap_legend = TRUE, width=unit(8, "cm")) + Heatmap(endo1 + 0, show_column_names = FALSE, name = "Endopeptidase inhibitor activity", col = c("0" = "white","1" = "blue"), show_heatmap_legend = TRUE, width = unit(2, "mm"))+Heatmap(atp1 + 0, show_column_names = FALSE, name = "ATP binding", col = c("0" = "white","1" = "goldenrod1"), show_heatmap_legend = TRUE, width = unit(2, "mm"))+Heatmap(heme1 + 0, show_column_names = FALSE, name = "Heme binding", col = c("0" = "white","1" = "red"), show_heatmap_legend = TRUE, width = unit(2, "mm"))+Heatmap(oxido1 + 0, show_column_names = FALSE, name = "Oxidoreductase activity", col = c("0" = "white","1" = "purple"), show_heatmap_legend = TRUE, width = unit(2, "mm"))+Heatmap(prot1 + 0, show_column_names = FALSE, name = "Protein binding", col = c("0" = "white","1" = "brown"), show_heatmap_legend = TRUE, width = unit(2, "mm")) +Heatmap(translate1 + 0, show_column_names = FALSE, name = "Translation", col = c("0" = "white","1" = "green"), show_heatmap_legend = TRUE, width = unit(2, "mm"))

expr_3 <- t(bsrel_sp_kvalue_logtrans_forpca)
p3 <- Heatmap(expr_3, name = "log (K selection value)", col = colorRamp2(c(-10,-5, 0, 5), c("dodgerblue2", "lightblue","white", "coral")), 
show_row_dend = FALSE,top_annotation = HeatmapAnnotation(df = SAMPLE,show_annotation_name = FALSE, col=colerz), column_order=colorder, show_row_names = FALSE, show_column_names = TRUE, heatmap_legend_param = list(legend_width = unit(3, "mm")),
column_names_gp=gpar(fontsize=5), show_heatmap_legend = TRUE, width=unit(6.5, "cm"))
htlist <- p3 + p2 + p1
draw(htlist, row_km =.1, cluster_rows=TRUE)

expr_3 <- expr_3[(row.names(expr_3) %in% row.names(expr_2)), ]
expr_3smol <- expr_2[!(row.names(expr_2) %in% row.names(expr_3)), ]
expr_3 <- rbind(expr_3, expr_3smol)













ForegroundSignificant <- filter(allResults, ForegroundBackground == "foreground" & Significant == 1)
ForegroundSignificant$Orthogroup <- as.factor(ForegroundSignificant$Orthogroup)
length(levels(ForegroundSignificant$Orthogroup))
#[1] 95
length(levels(as.factor(ForegroundSignificant$branch)))
#] 127
#ForegroundSignificant2 <- apply(ForegroundSignificant,2,as.character)
#write.csv(ForegroundSignificant2, "sig_foreground_results.csv")#

BackgroundSignificant <- filter(allResults, ForegroundBackground == "background" & Significant == 1)
BackgroundSignificant$Orthogroup <- as.factor(BackgroundSignificant$Orthogroup)
length(levels(BackgroundSignificant$Orthogroup))
#[1] 188
#BackgroundSignificant2 <- apply(BackgroundSignificant,2,as.character)
#write.csv(BackgroundSignificant2, "sig_background_results.csv")

#####How many orthogroups are selected upon in at least two of the pristimantis species?####
pristi_species_count <- ForegroundSignificant[,c('Species', 'Orthogroup')]
pristi_species_count <- pristi_species_count %>% group_by(Orthogroup) %>% tally()
morethan1 <- filter(pristi_species_count, n > 1)
morethan1
# # A tibble: 26 x 2
   # Orthogroup     n
   # <fct>      <int>
 # 1 OG0000008      3
 # 2 OG0000042      3
 # 3 OG0000068      2
 # 4 OG0000153      2
 # 5 OG0000287      3
 # 6 OG0000292      2
 # 7 OG0000382      2
 # 8 OG0000740      2
 # 9 OG0000844      2
# 10 OG0001050      2
# # ... with 16 more rows
dim(morethan1)
#[1] 26  2

###########AGAIN BUT RESULTS SECOND TIME##########

pristi_species_count <- ForegroundSignificant[,c('Species', 'Orthogroup')]
pristi_species_count <- pristi_species_count %>% group_by(Orthogroup) %>% tally()
morethan1 <- filter(pristi_species_count, n > 1)
morethan1
# A tibble: 92 x 2
   # Orthogroup     n
   # <fct>      <int>
 # 1 OG0000008      6
 # 2 OG0000027      2
 # 3 OG0000040      2
 # 4 OG0000041      2
 # 5 OG0000042      6
 # 6 OG0000047      2
 # 7 OG0000059      2
 # 8 OG0000067      2
 # 9 OG0000068      5
# 10 OG0000153      4






####Which orthogroups are selected upon in pristimantis but not others?####
onlyForeground <- ForegroundSignificant$Orthogroup[!(ForegroundSignificant$Orthogroup %in% BackgroundSignificant$Orthogroup)] 
write.csv(data.frame(unique(onlyForeground)),'only_foreground_orthogroups.csv')
length(unique(onlyForeground))
[1] 20
####Which orthogroups are selected upon in others but not pristimantis?####
onlyBackground <- BackgroundSignificant$Orthogroup[!(BackgroundSignificant$Orthogroup %in% ForegroundSignificant$Orthogroup)]
write.csv(data.frame(unique(onlyBackground)),'only_background_orthogroups.csv')
length(unique(onlyBackground))
[1] 113










###GETTING PERCENTAGES FOR REGRESSIONS####
allsignificant <- filter(allResults, `Corrected P-value` <= 0.05)
sig_percentages <- allsignificant %>%
	group_by(branch) %>%
	filter(Rate %in% c(max(Rate)))
sig_percentages_species <- data.frame(tapply(sig_percentages$Percent, sig_percentages$Species, mean))
write.csv(sig_percentages_species,  'percentages_species-significant.CORRECT.csv')
####GET HORIZONTAL transformed omega FOR HEATMAP for orthogroup by species#####	
allResults$omega_transformed <- ifelse(allResults$'Rate.classes' ==1, 1, allResults$'Baseline.MG94xREV.omega.ratio')
allResults$omega_transformed <- ifelse(allResults$'Corrected.P.value' >= 0.05, 1, allResults$omega_transformed)
allResults$omega_transformed <- ifelse(allResults$'Corrected.P.value' >= 0.05, 1, allResults$omega_transformed)
allResults$omega_transformed <- ifelse(allResults$omega_transformed >= 30, 15, allResults$omega_transformed)
omegas <- allResults %>%
	group_by(branch) %>%
	filter(Rate %in% c(max(Rate)))
omegas_transformed <- omegas[,c('Species', 'Orthogroup','omega_transformed')]
omegas_transformed <- omegas_transformed[-7543,]
omegas_transformed_final <- omegas_transformed %>% spread(Species, omega_transformed,convert=TRUE,fill=.0000001)
write.csv(omegas_transformed_final, 'heatmapvalues_omegas_transformed.csv')
####GET HORIZONTAL transformed percentages FOR HEATMAP for orthogroup by species#####
percentages <- allResults %>%
	group_by(branch) %>%
	filter(Rate %in% c(max(Rate)))
percentages_transformed <- percentages[,c('Species', 'Orthogroup','Percent_transformed')] 
percentages_transformed <- percentages_transformed[-7543,]
percentages_transformed_final <- percentages_transformed %>% spread(Species, Percent_transformed,convert=TRUE,fill=-600)
write.csv(percentages_transformed_final, 'heatmapvalues_percent_selected_transformed.csv')


#GET HORIZONTAL transformed p-values for heatmap#
#allsignificant <- filter(allResults, `Corrected.P.value` <= 0.05)
#allsignificant$Orthogroup <- as.factor(allsignificant$Orthogroup)
#length(levels(allsignificant$Orthogroup))
pvalues <- allResults %>%
	group_by(branch) %>%
	filter(Rate %in% c(max(Rate)))
pvalues <- pvalues[,c('Species', 'Orthogroup','Corrected.P.value')]
pvalues <- pvalues[-7543,]
pvalues_final <- pvalues %>% spread(Species, Corrected.P.value,convert=TRUE,fill=4)
write.csv(pvalues_final, 'heatmapvalues_pvalues_nottransformed.csv')


####GO ENRICHMENT FOR ABSREL STUFF###


###Make Orthogroup annotations into branch annotations for GO enrichment analyses###
cut -f1 cluster_domain_annotation.GO.txt | while read -r LINE; do grep $LINE all_results.csv | cut -f 17 -d, | while read -r PINE; do grep $LINE cluster_domain_annotation.GO.txt | sed "s/$LINE/$PINE/g" >> analyzed_branch_annotations.tsv;done ;done
cut -f1,3 analyzed_branch_annotations.tsv | sed 's/\t/,/g' | sed 's/"//g' > analyzed_branch_annotations_short.csv



###making shortened go term list with go terms distributed among branches for enrichment test##				
GOannotations2 <- read_delim("analyzed_branch_annotations_short.csv", delim=",")		
longAnnotations2 <- dplyr::select(GOannotations2, cluster_id, domain_id)
longAnnotations2 <- longAnnotations2 %>%
  dplyr::rename(branch = cluster_id)
longAnnotations2 <- longAnnotations2 %>% distinct()
wideListAnnotations2 <- tapply(longAnnotations2$domain_id, longAnnotations2$branch, function(x)x)
geneNames <- names(wideListAnnotations2)





####What are the enriched pathways/GO terms in Absrel positive branches branches unique to foreground-only orthogroups?#####
significanceInfo <- dplyr::select(allResults, branch, Orthogroup, holm_pvalue,ForegroundBackground )
significanceInfo <- significanceInfo %>%
	group_by(branch) %>%
	slice(which.min(holm_pvalue))
significanceInfo_pristionly <- ifelse(significanceInfo$Orthogroup %in% onlyForeground & significanceInfo$holm_pvalue < .05, 1, 0)


geneList <- significanceInfo_pristionly
names(geneList) <- significanceInfo$branch


###BIOOGICAL PROCESSES###
GOdataBP <- new("topGOdata",
              ontology = "BP",
              allGenes = geneList,
              geneSelectionFun = function(x)(x == 1),
              annot = annFUN.gene2GO, gene2GO = wideListAnnotations2)
resultFisherBP <- runTest(GOdataBP, algorithm = "elim", statistic = "fisher")
resultFisherBP
resultsFisherBPTable <- GenTable(GOdataBP, raw.p.value = resultFisherBP, topNodes = length(resultFisherBP@score),
numChar = 120)
write.table(resultsFisherBPTable, 'foreground-only-BP-enriched.txt')
"GO.ID" "Term" "Annotated" "Significant" "Expected" "raw.p.value"
"1" "GO:0018279" "protein N-linked glycosylation via asparagine" 40 2 0.12 "0.006"
"2" "GO:0006790" "sulfur compound metabolic process" 154 2 0.47 "0.076"
"3" "GO:0006807" "nitrogen compound metabolic process" 1805 8 5.51 "0.102"
"4" "GO:0050790" "regulation of catalytic activity" 38 1 0.12 "0.110"
"5" "GO:0065009" "regulation of molecular function" 38 1 0.12 "0.110"
"6" "GO:0000096" "sulfur amino acid metabolic process" 39 1 0.12 "0.113"
"7" "GO:0009069" "serine family amino acid metabolic process" 39 1 0.12 "0.113"
"8" "GO:0019346" "transsulfuration" 39 1 0.12 "0.113"
"9" "GO:0009092" "homoserine metabolic process" 39 1 0.12 "0.113"

###MOLECULAR FUNCTIONS###
GOdataMF <- new("topGOdata",
              ontology = "MF",
              allGenes = geneList,
              geneSelectionFun = function(x)(x == 1),
              annot = annFUN.gene2GO, gene2GO = wideListAnnotations2)
resultFisherMF <- runTest(GOdataMF, algorithm = "elim", statistic = "fisher")
resultFisherMF
resultsFisherMFTable <- GenTable(GOdataMF, raw.p.value = resultFisherMF, topNodes = length(resultFisherMF@score),
numChar = 120)
write.table(resultsFisherMFTable, 'foreground-only-MF-enriched.txt')
"GO.ID" "Term" "Annotated" "Significant" "Expected" "raw.p.value"
"1" "GO:0030170" "pyridoxal phosphate binding" 39 1 0.12 "0.11"
"2" "GO:0019842" "vitamin binding" 39 1 0.12 "0.11"
"3" "GO:0070279" "vitamin B6 binding" 39 1 0.12 "0.11"
"4" "GO:0008199" "ferric iron binding" 41 1 0.13 "0.12"
"5" "GO:0003700" "DNA binding transcription factor activity" 55 1 0.17 "0.16"
"6" "GO:0140110" "transcription regulator activity" 55 1 0.17 "0.16"
"7" "GO:0004869" "cysteine-type endopeptidase inhibitor activity" 60 1 0.19 "0.17"
"8" "GO:0008194" "UDP-glycosyltransferase activity" 63 1 0.2 "0.18"
"9" "GO:0004197" "cysteine-type endopeptidase activity" 79 1 0.25 "0.22"

####What are the enriched pathways/GO terms in Absrel positive branches branches unique to background-only orthogroups?#####
significanceInfo_background_only <- ifelse(significanceInfo$Orthogroup %in% onlyBackground & significanceInfo$holm_pvalue < .05, 1, 0)
significanceInfo_background_only
geneList <- significanceInfo_background_only
names(geneList) <- significanceInfo$branch

###BIOOGICAL PROCESSES###
GOdataBP <- new("topGOdata",
              ontology = "BP",
              allGenes = geneList,
              geneSelectionFun = function(x)(x == 1),
              annot = annFUN.gene2GO, gene2GO = wideListAnnotations2)
resultFisherBP <- runTest(GOdataBP, algorithm = "elim", statistic = "fisher")
resultFisherBP
resultsFisherBPTable <- GenTable(GOdataBP, raw.p.value = resultFisherBP, topNodes = length(resultFisherBP@score),
numChar = 120)
resultsFisherBPTable
write.table(resultsFisherBPTable, 'background-only-BP-enriched.txt')
"GO.ID" "Term" "Annotated" "Significant" "Expected" "raw.p.value"
"1" "GO:0006399" "tRNA metabolic process" 40 2 0.26 "0.027"
"2" "GO:0034660" "ncRNA metabolic process" 40 2 0.26 "0.027"
"3" "GO:0043038" "amino acid activation" 40 2 0.26 "0.027"
"4" "GO:0043039" "tRNA aminoacylation" 40 2 0.26 "0.027"
"5" "GO:0006418" "tRNA aminoacylation for protein translation" 40 2 0.26 "0.027"
"6" "GO:0016070" "RNA metabolic process" 136 3 0.88 "0.054"
"7" "GO:0090304" "nucleic acid metabolic process" 240 4 1.56 "0.063"
"8" "GO:0006520" "cellular amino acid metabolic process" 148 3 0.96 "0.067"
"9" "GO:0005975" "carbohydrate metabolic process" 72 2 0.47 "0.078"

###MOLECULAR FUNCTIONS###
GOdataBP <- new("topGOdata",
              ontology = "MF",
              allGenes = geneList,
              geneSelectionFun = function(x)(x == 1),
              annot = annFUN.gene2GO, gene2GO = wideListAnnotations2)
resultFisherMF <- runTest(GOdataBP, algorithm = "elim", statistic = "fisher")
resultFisherMF
resultsFisherMFTable <- GenTable(GOdataMF, raw.p.value = resultFisherMF, topNodes = length(resultFisherMF@score),
numChar = 120)
resultsFisherMFTable
write.table(resultsFisherMFTable, 'background-only-MF-enriched.txt')
"GO.ID" "Term" "Annotated" "Significant" "Expected" "raw.p.value"
"1" "GO:0140098" "catalytic activity, acting on RNA" 76 0 0.24 "0.024"
"2" "GO:0016874" "ligase activity" 40 0 0.12 "0.042"
"3" "GO:0016875" "ligase activity, forming carbon-oxygen bonds" 40 0 0.12 "0.042"
"4" "GO:0004812" "aminoacyl-tRNA ligase activity" 40 0 0.12 "0.042"
"5" "GO:0140101" "catalytic activity, acting on a tRNA" 40 0 0.12 "0.042"
"6" "GO:0005515" "protein binding" 889 1 2.76 "0.050"
"7" "GO:0005488" "binding" 3401 9 10.55 "0.077"
"8" "GO:0016788" "hydrolase activity, acting on ester bonds" 73 0 0.23 "0.120"
"9" "GO:0051213" "dioxygenase activity" 29 0 0.09 "0.213"



####Are the number of branches showing signs of selection significantly different for foreground vs background? ###CHI squared###
grouped_results <- allResults %>%
	group_by(branch) %>%
	slice(which.min(Corrected.P.value))
Pristimantis_branch_num_test <- table(grouped_results$ForegroundBackground,grouped_results$Significant)
>Pristimantis_branch_num_test
                0    1
  background 6831  414
  foreground 1170  127

> chisq.test(Pristimantis_branch_num_test)

        Pearson's Chi-squared test with Yates' continuity correction

data:  Pristimantis_branch_num_test
X-squared = 30.148, df = 1, p-value = 4.003e-08
>Prop.table(Pristimantis_branch_num_test,1)

                      0          1
  background 0.94285714 0.05714286
  foreground 0.90208173 0.09791827
-------- > SIGNIFICANTLY HIGHER PERCENTAGE FOR FOREGROUND <----------------

####Are the number of orthogroups showing signs of selection significantly different for foreground vs background? ###CHI squared###
grouped_results_by_orthogroup_background <- allResults %>%
	filter(ForegroundBackground == "background") %>%
	group_by(Orthogroup) %>%
	slice(which.min(Corrected.P.value))
grouped_results_by_orthogroup_foreground <- allResults %>%
	filter(ForegroundBackground == "foreground") %>%
	group_by(Orthogroup) %>%
	slice(which.min(Corrected.P.value))
grouped_results_orthogroup <- rbind(grouped_results_by_orthogroup_background, grouped_results_by_orthogroup_foreground)	
Pristimantis_orthogroup_num_test <- table(grouped_results_orthogroup$ForegroundBackground,grouped_results_orthogroup$Significant)

> Pristimantis_orthogroup_num_test

               0   1
  background  46 188
  foreground 139  95
chisq.test(Pristimantis_orthogroup_num_test)

        Pearson's Chi-squared test with Yates' continuity correction

data:  Pristimantis_orthogroup_num_test
X-squared = 75.659, df = 1, p-value < 2.2e-16
----------------------- > SIGNIFICANTLY HIGHER PERCENTAGE FOR BACKGROUND <------------------------
	
	
	
############### Are the proportion of branches different for foreground vs background species? ###PHYLANOVA
grouped_results <- allResults %>%
	group_by(branch) %>%
	slice(which.min(FDR_pvalue))
Pristimantis_branch_percentage_table <- table(grouped_results$Species,grouped_results$Significant)
write.table(Pristimantis_branch_percentage_table, "branch_percentages.tsv")
#put them in CSV file for local testing#	
setwd("~/GARD JSON FILES/Regressions/NO_CDHIT")
library(ape)
library(phytools)
library(geiger)
library(ggplot2)
tree <- read.tree('SpeciesTree_rooted.txt')
counts.df <- read.csv('absrel_counts_both_hypotheses.csv', row.names=1, header=TRUE)
selection <- counts.df$selection_avg
counts <- counts.df$value
grp <- as.factor(counts.df$Clade)
grp_pristi <- as.factor(counts.df$Pristi)
selection_percent <- as.vector(counts.df$selection_percent_pristi)
selection_percent_sigonly <- as.vector(counts.df$selection_percent_pristi_sigonly)
branch_percentages <- as.vector(counts.df$Proportion_sig_branches)
names(grp) <- rownames(counts.df)
names(counts) <- rownames(counts.df)
names(selection) <- rownames(counts.df)
names(grp_pristi) <- rownames(counts.df)
names(selection_percent) <- rownames(counts.df)
names(selection_percent_sigonly) <- rownames(counts.df)	
names(branch_percentages) <- rownames(counts.df)

> model <- phylANOVA(tree, grp_pristi, branch_percentages, nsim=1000, posthoc=TRUE, p.adj="bonferroni")
> model
ANOVA table: Phylogenetic ANOVA

Response: y
           Sum Sq  Mean Sq  F value Pr(>F)
x        0.013658 0.013658 3.363368  0.415
Residual 0.158375 0.004061                

P-value based on simulation.
---------

Pairwise posthoc test using method = "bonferroni"

Pairwise t-values:
                Other Pristimantis
Other        0.000000    -1.833949
Pristimantis 1.833949     0.000000

Pairwise corrected P-values:
             Other Pristimantis
Other        1.000        0.415			----------------------------------------------------------->	###NOT SIGNIFICANT####
Pristimantis 0.415        1.000

p = ggplot(counts.df, aes(x = Group.1, y = branch_percentages, group = Pristi, fill = Pristi)) +
  geom_bar(stat="identity",color='black', position = "dodge") +
  xlab('Species') + ylab('Rates of positively selected branches') + theme_bw()+
  theme(text = element_text(size=16),
        axis.text.x = element_text(angle=0, hjust=.5),
        plot.title = element_text(size=12,hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5))+
  ggtitle(" % aBSREL branches under selection")



	
############### Are the proportion of positivelys selected sites different for foreground vs background species? ###PHYLANOVA
model <- phylANOVA(tree, grp_pristi, selection_percent_sigonly, nsim=1000, posthoc=TRUE, p.adj="bonferroni")
> model
ANOVA table: Phylogenetic ANOVA

Response: y
           Sum Sq  Mean Sq  F value Pr(>F)
x        0.001224 0.001224 0.940875  0.656
Residual 0.050752 0.001301                

P-value based on simulation.
---------

Pairwise posthoc test using method = "bonferroni"

Pairwise t-values:
                Other Pristimantis
Other        0.000000    -0.969987
Pristimantis 0.969987     0.000000

Pairwise corrected P-values:
             Other Pristimantis
Other        1.000        0.656
Pristimantis 0.656        1.000
---------

> model <- phylANOVA(tree, grp_pristi, selection_percent, nsim=1000, posthoc=TRUE, p.adj="bonferroni")
> model
ANOVA table: Phylogenetic ANOVA

Response: y
           Sum Sq  Mean Sq   F value Pr(>F)
x        0.249993 0.249993 13.409459  0.085
Residual 0.727080 0.018643                 

P-value based on simulation.
---------

Pairwise posthoc test using method = "bonferroni"

Pairwise t-values:
                Other Pristimantis
Other        0.000000    -3.661893
Pristimantis 3.661893     0.000000

Pairwise corrected P-values:
             Other Pristimantis
Other        1.000        0.085
Pristimantis 0.085        1.000