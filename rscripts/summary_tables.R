

library(plyr)
library(tidyverse)
library(rjson)
#####ABSREL RESULTS(VERY BASIC)################

# Write a function to read in the aBSREL JSON result file and transform it into a dataframe:
foregroundabsrelJSONProcessing <- function(file) {
  # Read in the JSON file:
  absrel <- rjson::fromJSON(file = file)
  absrelResult <- data.table::rbindlist(absrel[["branch attributes"]][["0"]], fill = TRUE, idcol = TRUE)
    #Add a column to hold the orthogroup name:
  orthogroup <- sapply(strsplit(absrel[["input"]][["file name"]],"/"), tail, 1)
  orthogroup <- sapply(strsplit(orthogroup, "[.]"), head, 1)
  absrelResult$Orthogroup <- orthogroup
  absrelResult <- absrelResult[!duplicated(absrelResult[,c('Orthogroup')]),]
  # Create a data frame with the information for each branch tested:
  positive <- absrel[["test results"]][2]
  tested <- absrel[["test results"]][3]
  branch_ratio <- paste(positive,tested, sep='/')
  absrelResult$foregroundpositive <- branch_ratio
  absrelResult <- absrelResult %>% select(Orthogroup,foregroundpositive)
  return(absrelResult)
}
foregroundFiles <- list.files(path = "./foreground_absrel", pattern = "*.json", full.names = TRUE)
foregroundFiles <- foregroundFiles[sapply(foregroundFiles, file.size) > 0]
possiblyforegroundabsrelJSONProcessing <- possibly(foregroundabsrelJSONProcessing, otherwise = "File empty.")

foregroundaBSRELResults <- map_dfr(foregroundFiles, possiblyforegroundabsrelJSONProcessing)

################ Do the same for background files:##############
backgroundabsrelJSONProcessing <- function(file) {
  # Read in the JSON file:
  absrel <- rjson::fromJSON(file = file)
  absrelResult <- data.table::rbindlist(absrel[["branch attributes"]][["0"]], fill = TRUE, idcol = TRUE)
    #Add a column to hold the orthogroup name:
  orthogroup <- sapply(strsplit(absrel[["input"]][["file name"]],"/"), tail, 1)
  orthogroup <- sapply(strsplit(orthogroup, "[.]"), head, 1)
  absrelResult$Orthogroup <- orthogroup
  absrelResult <- absrelResult[!duplicated(absrelResult[,c('Orthogroup')]),]
  # Create a data frame with the information for each branch tested:
  positive <- absrel[["test results"]][2]
  tested <- absrel[["test results"]][3]
  branch_ratio <- paste(positive,tested, sep='/')
  absrelResult$backgroundpositive <- branch_ratio


  # Add a column to label genes as having been part of the foreground or background

  absrelResult <- absrelResult %>% select(Orthogroup,backgroundpositive)
  return(absrelResult)
}
backgroundFiles <- list.files(path = "./background_absrel", pattern = "*.json", full.names = TRUE)
backgroundFiles <- backgroundFiles[sapply(backgroundFiles, file.size) > 0]
possiblybackgroundabsrelJSONProcessing <- possibly(backgroundabsrelJSONProcessing, otherwise = "File empty.")
# Map the function over all foreground files:
backgroundaBSRELResults <- map_dfr(backgroundFiles, possiblybackgroundabsrelJSONProcessing)
allResults <- foregroundaBSRELResults %>% left_join(backgroundaBSRELResults,by=c("Orthogroup" = "Orthogroup"))


###############################################
####### BUSTED RESULTS #####
###############################################


# Construct a list of all of the json files:
bustedjsonFiles <- list.files(path = './BUSTED', pattern = "*.json", full.names = TRUE)
bustedjsonFiles <- sort(bustedjsonFiles, decreasing = TRUE)
# Write a function that will process each individual json file and extract the file name, orthogroup number, p-value, and return a text description of the p-value:
bustedJSONProcessing <- function(file) {
  busted <- rjson::fromJSON(file = file)
  bustedResult <- data.frame(cbind(busted[["test results"]][["p-value"]],
  busted[["fits"]][["Unconstrained model"]][["Rate Distributions"]][["Background"]][["0"]][["omega"]],
  busted[["fits"]][["Unconstrained model"]][["Rate Distributions"]][["Background"]][["0"]][["proportion"]],
  busted[["fits"]][["Unconstrained model"]][["Rate Distributions"]][["Background"]][["1"]][["omega"]],
  busted[["fits"]][["Unconstrained model"]][["Rate Distributions"]][["Background"]][["1"]][["proportion"]],
  busted[["fits"]][["Unconstrained model"]][["Rate Distributions"]][["Background"]][["2"]][["omega"]],
  busted[["fits"]][["Unconstrained model"]][["Rate Distributions"]][["Background"]][["2"]][["proportion"]],
  busted[["fits"]][["Unconstrained model"]][["Rate Distributions"]][["Test"]][["0"]][["omega"]],
  busted[["fits"]][["Unconstrained model"]][["Rate Distributions"]][["Test"]][["0"]][["proportion"]],
  busted[["fits"]][["Unconstrained model"]][["Rate Distributions"]][["Test"]][["1"]][["omega"]],
  busted[["fits"]][["Unconstrained model"]][["Rate Distributions"]][["Test"]][["1"]][["proportion"]],
  busted[["fits"]][["Unconstrained model"]][["Rate Distributions"]][["Test"]][["2"]][["omega"]],
  busted[["fits"]][["Unconstrained model"]][["Rate Distributions"]][["Test"]][["2"]][["proportion"]]))
    #Add a column to hold the orthogroup name:
  backrates1 <- paste(bustedResult$X2,bustedResult$X3, sep=':')
  backrates2 <- paste(bustedResult$X4,bustedResult$X5, sep=':')
  backrates3 <- paste(bustedResult$X6,bustedResult$X7, sep=':')
  forerates1 <- paste(bustedResult$X8,bustedResult$X9, sep=':')
  forerates2 <- paste(bustedResult$X10,bustedResult$X11, sep=':')
  forerates3 <- paste(bustedResult$X12,bustedResult$X13, sep=':')
  forerates <- paste(forerates1, forerates2, forerates3, collapse='|')
  backrates <- paste(backrates1, backrates2, backrates3, collapse='|')
  orthogroup <- sapply(strsplit(busted[["input"]][["file name"]],"/"), tail, 1)
  orthogroup <- sapply(strsplit(orthogroup, "[.]"), head, 1)
  bustedResult$Orthogroup <- orthogroup
  bustedResult$forerates <- forerates
  bustedResult$backrates <- backrates
  names(bustedResult)[1] <- "busted_pvalue"
  bustedResult <- bustedResult %>% select(Orthogroup,busted_pvalue,forerates,backrates) 

}
possiblyBustedJSONProcessing <- possibly(bustedJSONProcessing, otherwise = "File empty.")
bustedResults <- map(bustedjsonFiles, possiblyBustedJSONProcessing)

# Convert the results to a dataframe:
bustedResults <- as.data.frame(do.call(rbind, bustedResults))   
bustedResults$busted_adj_pvalue <- p.adjust(bustedResults$busted_pvalue, "BH")
bustedResults <- bustedResult %>% select(Orthogroup,busted_adj_pvalue,forerates,backrates)
allResults <- allResults %>% left_join(bustedResults,by=c("Orthogroup" = "Orthogroup"))


#############RELAX RESULTS#############
########################################
jsonFiles <- list.files(path = './RELAX', pattern = "*.json", full.names = TRUE)
jsonFiles <- sort(jsonFiles, decreasing = TRUE)
jsonFiles <- jsonFiles[sapply(jsonFiles, file.size) > 0]
#########################################################################
##################ORTHOGROUP LEVEL RESULTS#################
#########################################################################
# Write a function to read in the relax JSON result file and transform it into a dataframe:
relaxJSONProcessing <- function(file) {
  relaxResult <- rjson::fromJSON(file = file)
  if (relaxResult[["test results"]][["p-value"]] < 0.05) {
    if (relaxResult[["test results"]][["relaxation or intensification parameter"]] > 1) {
      data <- c(relaxResult[["input"]][["file name"]], relaxResult[["test results"]][["p-value"]], relaxResult[["test results"]][["relaxation or intensification parameter"]], "Intensification of selection along foreground branches")
      return(data)
    } else {
      data <- c(relaxResult[["input"]][["file name"]], relaxResult[["test results"]][["p-value"]], relaxResult[["test results"]][["relaxation or intensification parameter"]], "Relaxation of selection along foreground branches")
      return(data)
    }
  } else {
    if (relaxResult[["test results"]][["relaxation or intensification parameter"]] > 1) {
      data <- c(relaxResult[["input"]][["file name"]], relaxResult[["test results"]][["p-value"]], relaxResult[["test results"]][["relaxation or intensification parameter"]], "Nonsignificant intensification")
      return(data)
    } else {
      data <- c(relaxResult[["input"]][["file name"]], relaxResult[["test results"]][["p-value"]], relaxResult[["test results"]][["relaxation or intensification parameter"]], "Nonsignificant relaxation")
      return(data)
    }
  }
}
possiblyRelaxJSONProcessing <- possibly(relaxJSONProcessing, otherwise = "File empty.")

relaxResults <- map(jsonFiles, possiblyRelaxJSONProcessing)

# Convert the results to a dataframe:
relaxResults <- as.data.frame(do.call(rbind, relaxResults))   
colnames(relaxResults) <- c("fileName", "relax_pValue", "kValue", "shortDescription")
orthogroup <- sapply(strsplit(as.character(relaxResults$fileName),"/"), tail, 1)
orthogroup <- orthogroup %>% word(1, sep = fixed('.'))
relaxResults$Orthogroup <- orthogroup
#relaxResults$relax_adj_pvalue <- p.adjust(relaxResults$pValue, "BH")
relaxResults <- relaxResults %>% select(Orthogroup,relax_pValue,kValue,shortDescription)
allResults <- allResults %>% left_join(relaxResults,by=c("Orthogroup" = "Orthogroup"))
#write.csv(allResults, 'direct_dev_table.csv')
write.csv(allResults, 'pristimantis_table.csv')


#######Which are our candidate genes? Things that have at least one significant result###
library(ComplexHeatmap)
library(ComplexUpset)
allResults <- read.csv('dir_dev.csv')
allResults <- separate(allResults, foregroundpositive,  c("foregroundpositive", "foregroundtested"), sep="/")
Significant_atleast1 <- filter(allResults, foregroundpositive > 1 | busted_adj_pvalue <= .05 | relax_pValue <= .05 )
Orthogroups <- allResults$Orthogroup
relaxSignificant <- filter(allResults,  relax_pValue <= .05 )
bustedSignificant <- filter(allResults, busted_adj_pvalue <= .05)
absrelSignificant <- filter(allResults, foregroundpositive > 1)
is_relax <- data.frame(as.integer(Orthogroups %in% relaxSignificant$Orthogroup))
is_positive <- data.frame(as.integer(Orthogroups %in% bustedSignificant$Orthogroup))
is_diverse <- data.frame(as.integer(Orthogroups %in% absrelSignificant$Orthogroup))

d <- cbind(Orthogroups, is_relax, is_positive, is_diverse)
colnames(d) <- c("Orthogroups", "Intensified or Relaxed Selection (RELAX)", "Positive selection (BUSTED)", "Positive selection (aBSREL)")
tests <- names(d)
dcomb <- make_comb_mat(d, mode ="distinct")
upset(data = d, intersect = tests, name = "Selection test groupings by frequency. Total pool is 232 orthogroups.") + labs(title= "Overlap of selection type (orthogroups)")

###Is foreground different than background?####
frequencyTablepristi <- (c(172, 616, 1879, 13392),
                                           nrow = 2,
                                           byrow = TRUE,
                                           dimnames = list("category" = c("evidenceForSelection", "nonsignificantResult") ,
"selection" = c("foreground", "background")))
fishersExactTest_pristi <- fisher.test(frequencyTablepristi)
fishersExactTest_pristi[["p.value"]]

frequencyTabledd <- matrix(c(201, 604, 2336, 13095),
                                           nrow = 2,
                                           byrow = TRUE,
                                           dimnames = list("category" = c("evidenceForSelection", "nonsignificantResult") ,
"selection" = c("foreground", "background")))
fishersExactTest_dd <- fisher.test(frequencyTabledd)
fishersExactTest_dd[["p.value"]]




############PLOTTING DIFFERENCES IN AVERAGE RATES FROM BUSTED BETWEEN FOREGROUND AND BACKGROUND##########
###PRISTIMANTIS###
allResults_avgs_pristi <- allResults %>% separate(forerates, c("A","B","C"), sep = " ")
allResults_avgs_pristi <- allResults_avgs_pristi %>% separate(backrates, c("D","E","F"), sep = " ")
allResults_avgs_pristi <- allResults_avgs_pristi %>% separate(A, c("A1","A2"), sep = ":")
allResults_avgs_pristi$forerateA <- as.numeric(allResults_avgs_pristi$A1)*as.numeric(allResults_avgs_pristi$A2)
allResults_avgs_pristi <- allResults_avgs_pristi %>% separate(B, c("B1","B2"), sep = ":")
allResults_avgs_pristi$forerateB <- as.numeric(allResults_avgs_pristi$B1)*as.numeric(allResults_avgs_pristi$B2)
allResults_avgs_pristi <- allResults_avgs_pristi %>% separate(C, c("C1","C2"), sep = ":")
allResults_avgs_pristi$forerateC <- as.numeric(allResults_avgs_pristi$C1)*as.numeric(allResults_avgs_pristi$C2)
allResults_avgs_pristi <- allResults_avgs_pristi %>% separate(D, c("D1","D2"), sep = ":")
allResults_avgs_pristi$backrateD <- as.numeric(allResults_avgs_pristi$D1)*as.numeric(allResults_avgs_pristi$D2)
allResults_avgs_pristi <- allResults_avgs_pristi %>% separate(E, c("E1","E2"), sep = ":")
allResults_avgs_pristi$backrateE <- as.numeric(allResults_avgs_pristi$E1)*as.numeric(allResults_avgs_pristi$E2)
allResults_avgs_pristi <- allResults_avgs_pristi %>% separate(F, c("F1","F2"), sep = ":")
allResults_avgs_pristi$backrateF <- as.numeric(allResults_avgs_pristi$F1)*as.numeric(allResults_avgs_pristi$F2)
allResults_avgs_pristi$backrate_avg <- allResults_avgs_pristi$backrateF + allResults_avgs_pristi$backrateE + allResults_avgs_pristi$backrateD
allResults_avgs_pristi$forerate_avg <- allResults_avgs_pristi$forerateA + allResults_avgs_pristi$forerateB + allResults_avgs_pristi$forerateC
write.csv(allResults_avgs_pristi, 'pristiavgs.csv')
d <- read.csv('pristiavgs.csv')
ylim1 = boxplot.stats(d$rate_avg)$stats[c(1, 5)]
p <- ggplot(data = d, aes(x = foreground_background, y = rate_avg)) + geom_boxplot(outlier.shape=NA, alpha = 0.5) + coord_cartesian(ylim = ylim1*1.05) + geom_point(aes(fill = foreground_background), size = 2, shape = 21, position = position_jitterdodge()) +  theme(legend.position="none", plot.title = element_text(size = 11, hjust = .5)) + ggtitle("Pristimantis foreground dN/dS distribution") + ylab("average dN/dS")
###direct developers###
allResults_avgs_dirdev <- allResults2 %>% separate(forerates, c("A","B","C"), sep = " ")
allResults_avgs_dirdev <- allResults_avgs_dirdev %>% separate(backrates, c("D","E","F"), sep = " ")
allResults_avgs_dirdev <- allResults_avgs_dirdev %>% separate(A, c("A1","A2"), sep = ":")
allResults_avgs_dirdev$forerateA <- as.numeric(allResults_avgs_dirdev$A1)*as.numeric(allResults_avgs_dirdev$A2)
allResults_avgs_dirdev <- allResults_avgs_dirdev %>% separate(B, c("B1","B2"), sep = ":")
allResults_avgs_dirdev$forerateB <- as.numeric(allResults_avgs_dirdev$B1)*as.numeric(allResults_avgs_dirdev$B2)
allResults_avgs_dirdev <- allResults_avgs_dirdev %>% separate(C, c("C1","C2"), sep = ":")
allResults_avgs_dirdev$forerateC <- as.numeric(allResults_avgs_dirdev$C1)*as.numeric(allResults_avgs_dirdev$C2)
allResults_avgs_dirdev <- allResults_avgs_dirdev %>% separate(D, c("D1","D2"), sep = ":")
allResults_avgs_dirdev$backrateD <- as.numeric(allResults_avgs_dirdev$D1)*as.numeric(allResults_avgs_dirdev$D2)
allResults_avgs_dirdev <- allResults_avgs_dirdev %>% separate(E, c("E1","E2"), sep = ":")
allResults_avgs_dirdev$backrateE <- as.numeric(allResults_avgs_dirdev$E1)*as.numeric(allResults_avgs_dirdev$E2)
allResults_avgs_dirdev <- allResults_avgs_dirdev %>% separate(F, c("F1","F2"), sep = ":")
allResults_avgs_dirdev$backrateF <- as.numeric(allResults_avgs_dirdev$F1)*as.numeric(allResults_avgs_dirdev$F2)
allResults_avgs_dirdev$backrate_avg <- allResults_avgs_dirdev$backrateF + allResults_avgs_dirdev$backrateE + allResults_avgs_dirdev$backrateD
allResults_avgs_dirdev$forerate_avg <- allResults_avgs_dirdev$forerateA + allResults_avgs_dirdev$forerateB + allResults_avgs_dirdev$forerateC
write.csv(allResults_avgs_dirdev, 'dirdevavgs.csv')
dd <- read.csv('dirdevavgs.csv')
ylim1 = boxplot.stats(dd$rate_avg)$stats[c(1, 5)]
q <- ggplot(data = dd, aes(x = foreground_background, y = rate_avg)) + geom_boxplot(outlier.shape=NA, alpha = 0.5) + coord_cartesian(ylim = ylim1*1.05) + geom_point(aes(fill = foreground_background), size = 2, shape = 21, position = position_jitterdodge()) +  theme(legend.position="none", plot.title = element_text(size = 11, hjust = .5)) + ggtitle("Direct developer foreground dN/dS distribution") + ylab("average dN/dS")





######MAKING PLOTS OF RELAX RESULTS AND FISHER TESTS FOR RELAXED/INTENSIFIED ORTHOGROUPS#########
###############################################################################################

library(plyr)
library(tidyverse)
library(rjson)
library(ggplot2)
library(grid)
library(gridExtra)
relaxResults <- read.csv('dirdev_relax_results.csv')
#relaxResults <- read.csv('pristimantis_relax_results.csv')
relaxframe <- relaxResults
relaxframe <- subset(relaxframe, fileName != "File empty.")
relaxframe$qval <- NA
relaxframe$sig <- NA

relaxframe$qval <- p.adjust(relaxframe$relax_pValue, method = "fdr", n = nrow(relaxframe))
Bonfcutoff <- (0.05 / nrow(relaxframe))

# We want to count the number of genes that are under relaxed and intensified selection. This person chose to do this with a for loop, so first create values set at zero. Each time the for loop encounters a row that would fit into one of these categories, it will increase the counter by 1. 
relaxedWithNoCorrection <- 0
relaxedWithFDR <- 0
relaxedWithBonferroni <- 0
intensifiedWithNoCorrection <- 0
intensifiedWithFDR <- 0
intensifiedWithBonferroni <- 0
for (i in 1:nrow(relaxframe)) {
 if ((relaxframe$relax_pValue[i] < 0.05) && (relaxframe$kValue[i] < 1)) {
    relaxedWithNoCorrection <- relaxedWithNoCorrection + 1
    relaxframe$sig[i] <- "relaxedWithNoCorrection"
  } 
  if ((relaxframe$qval[i] < 0.10) && (relaxframe$kValue[i] < 1)) {
    relaxedWithFDR <- relaxedWithFDR + 1
    relaxframe$sig[i] <- "relaxedWithFDR"
  } 
  if ((relaxframe$relax_pValue[i] < Bonfcutoff) && (relaxframe$kValue[i] < 1)) {
    relaxedWithBonferroni <- relaxedWithBonferroni + 1
    relaxframe$sig[i] <- "relaxedWithBonferroni"
  } 
  if ((relaxframe$relax_pValue[i] < 0.05) && (relaxframe$kValue[i] > 1)) {
    intensifiedWithNoCorrection <- intensifiedWithNoCorrection + 1
    relaxframe$sig[i] <- "intensifiedWithNoCorrection"
  } 
  if ((relaxframe$qval[i] < 0.10) && (relaxframe$kValue[i] > 1)) {
    intensifiedWithFDR <- intensifiedWithFDR + 1
    relaxframe$sig[i] <- "intensifiedWithFDR"
  } 
  if ((relaxframe$relax_pValue[i] < Bonfcutoff) && (relaxframe$kValue[i] > 1)) {
    intensifiedWithBonferroni <- intensifiedWithBonferroni + 1
    relaxframe$sig[i] <- "intensifiedWithBonferroni"
  }
}

frequencyTableWithNoCorrection <- matrix(c(relaxedWithNoCorrection, intensifiedWithNoCorrection, 
                                           (nrow(relaxframe) - relaxedWithNoCorrection),
                                           (nrow(relaxframe) - intensifiedWithNoCorrection)),
                                         nrow = 2,
                                         byrow = TRUE,
                                         dimnames = list("category" = c("outliers", "non-outliers") ,
                                                         "selection" = c("relaxed", "intensified")))
frequencyTableWithFDR <- matrix(c(relaxedWithFDR, intensifiedWithFDR, 
                                  (nrow(relaxframe) - relaxedWithFDR),
                                  (nrow(relaxframe) - intensifiedWithFDR)),
                                nrow = 2,
                                byrow = TRUE,
                                dimnames = list("category" = c("outliers", "non-outliers") ,
                                                "selection" = c("relaxed", "intensified")))

frequencyTableWithBonferroni <- matrix(c(relaxedWithBonferroni, intensifiedWithBonferroni, 
                                         (nrow(relaxframe) - relaxedWithBonferroni),
                                         (nrow(relaxframe) - intensifiedWithBonferroni)),
                                       nrow = 2,
                                       byrow = TRUE,
                                       dimnames = list("category" = c("outliers", "non-outliers") ,
                                                       "selection" = c("relaxed", "intensified")))

below1 <- relaxframe$kValue[which(relaxframe$kValue < 1)]
above1 <- relaxframe$kValue[which(relaxframe$kValue > 1)]

freq_table_below_or_above_1 <- matrix(c(length(below1), length(above1), 
                                        1347.5, 1347.5),
                                      nrow = 2,
                                      byrow = TRUE,
                                      dimnames = list("category" = c("outliers", "non-outliers") ,
                                                      "selection" = c("relaxed", "intensified")))

fishersExactTest <- fisher.test(frequencyTableWithNoCorrection)
fisher.test(frequencyTableWithNoCorrection)
fisher.test(frequencyTableWithFDR)
fisher.test(frequencyTableWithBonferroni)
fisher.test(freq_table_below_or_above_1)

freq_table_selparam <- matrix(c((length(relaxframe$kValue[relaxframe$kValue != 1]) / 2), 
                                (length(relaxframe$kValue[relaxframe$kValue != 1]) / 2), 
                                length(relaxframe$kValue[relaxframe$kValue < 1]),
                                length(relaxframe$kValue[relaxframe$kValue > 1])),
                              nrow = 2,
                              byrow = TRUE,
                              dimnames = list("category" = c("expected", "observed") ,
"selection" = c("relaxed", "intensified")))

print(fishersExactTest[["p.value"]])

pValue <- if (fishersExactTest[["p.value"]] < 0.000001) {
  formatC(fishersExactTest[["p.value"]], format = "e", digits = 2)
} else {
  round(fishersExactTest[["p.value"]], digits = 4)
}
print(pValue)
pValueLabel <- paste("Fisher p-value = ", pValue, sep = "")


pValue2 <- if (fishersExactTest[["p.value"]] < 0.000001) {
  formatC(fishersExactTest[["p.value"]], format = "e", digits = 2)
} else {
  round(fishersExactTest[["p.value"]], digits = 4)
}
print(pValue2)
pValueLabel2 <- paste("Fisher p-value = ", pValue2, sep = "")




relaxResults$shortDescription <-
  factor(relaxResults$shortDescription,
         levels = c("Intensification of selection along foreground branches", 
                    "Relaxation of selection along foreground branches",
                    "Nonsignificant intensification",
                    "Nonsignificant relaxation",
                    "File empty."))
relaxResults$kValue <- as.numeric(as.character(relaxResults$kValue))
plot <- ggplot(data = filter(relaxResults, 
shortDescription != "File empty."))

                       

# Combine that base plot with a bar plot to show the number of genes intensified, relaxed, and NS:
p <- plot + 
  geom_bar(mapping = aes(x = shortDescription)) +
  labs(x = "Selective regime", 
       y = "Count of orthogroups", 
       title = "Distribution of RELAX selection regimes with Pristimantis foreground") +
  theme(axis.text.x = element_text(angle = 0,size = 1,
                                 hjust = 0.5,
                                 vjust = 1),axis.title.x = element_text(size=2),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank()) + 
  scale_x_discrete(labels=c("Nonsignificant relaxation"  = "Nonsignificant\nrelaxation" ,
                            "Intensification of selection along foreground branches" = "Intensification of\nselection along\nforeground branches",
                            "Nonsignificant intensification" = "Nonsignificant\nintensification",
                            "Relaxation of selection along foreground branches" = "Relaxation of\nselection along\nforeground branches",
                            "File empty." = "File empty.")) + theme(title = element_text(size = 3)) +
  geom_text(aes(x = 2, 
                y = (100),
                label = pValueLabel),
            stat = "unique", size=3) + 
  theme_bw() 
p <- p + theme(plot.title = element_text(size =11))
######## Create a plot that shows the distribution of k-values:

q <- ggplot(data = filter(relaxResults, 
                     shortDescription != "File empty.", 
                     kValue < 2.5)) + 
  geom_histogram(mapping = aes(x = kValue),
                 binwidth = 0.05,
                 center = 0,
                 color = "grey26",
                 alpha = 0) + 
  geom_vline(xintercept = 1,
             color = "darkred",
             size = 0.75) +
  geom_text(aes(x = 1.05, y = 28,
                label = "Neutral expectation, k = 1"),
            stat = "unique",
            hjust = 0, 
            color = "darkred", size=2.5) +
  theme_bw() +   ylim(0, 30) +
  labs(x = "Relaxation parameter, k", 
       y = "Count of orthogroups") +
  theme(axis.text.x = element_text(angle = 0,
                                   hjust = 0.5,
                                   vjust = 1),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank())
# Combine that base plot with a bar plot to show the number of genes intensified, relaxed, and NS:
r <- plot + 
  geom_bar(mapping = aes(x = shortDescription)) +
  labs(x = "Selective regime", 
       y = "Count of orthogroups", 
       title = "Distribution of RELAX selective regimes with direct developers foreground") +
  theme(axis.text.x = element_text(angle = 0, size=1,
                                 hjust = 0.5,
                                 vjust = 1),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank()) + theme(plot.title = element_text(size = 5)) +
  scale_x_discrete(labels=c("Nonsignificant relaxation"  = "Nonsignificant\nrelaxation" ,
                            "Intensification of selection along foreground branches" = "Intensification of\nselection along\nforeground branches",
                            "Nonsignificant intensification" = "Nonsignificant\nintensification",
                            "Relaxation of selection along foreground branches" = "Relaxation of\nselection along\nforeground branches",
                            "File empty." = "File empty.")) +
  geom_text(aes(x = 2, 
                y = (100),
                label = pValueLabel2),
            stat = "unique",size=3) + 
  theme_bw() 
r <- r + theme(plot.title = element_text(size = 11))
#r + theme(text = element_text(size=rel(3)))
######## Create a plot that shows the distribution of k-values:

s <- ggplot(data = filter(relaxResults, 
                     shortDescription != "File empty.", 
                     kValue < 2.5)) + 
  geom_histogram(mapping = aes(x = kValue),
                 binwidth = 0.05,
                 center = 0,
                 color = "grey26",
                 alpha = 0) + 
  geom_vline(xintercept = 1,
             color = "darkred",
             size = 0.75) +
  geom_text(aes(x = 1.05, y = 28,
                label = "Neutral expectation, k = 1"),
            stat = "unique",
            hjust = 0, 
            color = "darkred", size=2.5) +
  theme_bw() +   ylim(0, 30) +
  labs(x = "Relaxation parameter, k", 
       y = "Count of orthogroups") +
  theme(axis.text.x = element_text(angle = 0,
                                   hjust = 0.5,
                                   vjust = 1),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank())
#png("dirdev_pristi_distributions_relax.png")
#grid.arrange(p, q,r,s)
#dev.off()theme(axis.text.x  = element_text(angle = 90, size = 8))
#lay = rbind(c(1,2,2))
lay = rbind(c(1,2,2),
            c(3,4,4))
	    
png("pristi_direct_distributions_final.png",width = 9, height = 5, units = "in", res = 600)                        
print(grid.arrange(arrangeGrob(q, left = textGrob("a)", x = unit(1, "npc"), 
                               y = unit(.95, "npc"))), 
                   arrangeGrob(p, left =textGrob("b)", x = unit(1, "npc"), 
                               y = unit(.95, "npc"))),
		arrangeGrob(s, left = textGrob("c)", x = unit(1, "npc"), 
                               y = unit(1, "npc"))), 
                   arrangeGrob(r, left =textGrob("d)", x = unit(1, "npc"), 
                               y = unit(1, "npc"))),
                   layout_matrix = lay))
dev.off()
###########Idk#####
relaxResults$pValue <- as.numeric(relaxResults$pValue)
relaxResults$kValue <- as.numeric(relaxResults$kValue)

ggplot(data = relaxResults) +
  geom_hex(mapping = aes(x = pValue,
                         y = kValue),
           bins = 50) +
geom_vline(xintercept = 0.05)
