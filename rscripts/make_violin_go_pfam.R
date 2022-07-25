


####NOW RUN IT####library(tidyverse)
library(readr)
library(timecourseRnaseq)
library(dplyr)
setwd("~/Orthogroup_paper_edits_2022/TO_DO_FOR_PAPER")

###GOS###
GOannotations <- read_delim("cluster_domain_annotation.GO.txt", delim = "\t")
GOannotations
longGOAnnotations <- dplyr::select(GOannotations, `#cluster_id`, domain_id, domain_description)
longGOAnnotations

longGOAnnotations <- longGOAnnotations %>% group_by(`#cluster_id`)  %>% mutate(gos_by_foo = paste(domain_id, collapse = ";"))
longGOAnnotations <- longGOAnnotations %>% group_by(`#cluster_id`)  %>% mutate(desc_by_foo = paste(domain_description, collapse = ";"))
longGOAnnotations <- longGOAnnotations[!duplicated(longGOAnnotations$`#cluster_id`),]
shortGOAnnotations <- dplyr::select(longGOAnnotations, `#cluster_id`,gos_by_foo, desc_by_foo)

#####PFAMS####
pfamannotations <- read_delim("cluster_domain_annotation.Pfam.txt", delim = "\t")
pfamannotations
longpfamAnnotations <- dplyr::select(pfamannotations, `#cluster_id`, domain_id, domain_description)
longpfamAnnotations
longpfamAnnotations <- cbind(pfams_new=paste(longpfamAnnotations$domain_id,longpfamAnnotations$domain_description, sep='|'), longpfamAnnotations)


longpfamAnnotations <- longpfamAnnotations %>% group_by(`#cluster_id`)  %>% mutate(pfams_by_foo = paste(domain_id, collapse = ";"))
longpfamAnnotations <- longpfamAnnotations %>% group_by(`#cluster_id`)  %>% mutate(pdesc_by_foo = paste(domain_description, collapse = ";"))
longpfamAnnotations <- longpfamAnnotations[!duplicated(longpfamAnnotations$`#cluster_id`),]
shortpfamAnnotations <- dplyr::select(longpfamAnnotations, `#cluster_id`, pfams_by_foo,pdesc_by_foo)
total <- merge(shortpfamAnnotations,shortGOAnnotations,by="#cluster_id")
write.csv(total,'annotation_file_violin_plot.csv')

violin_annots <- read.csv('annotation_file_violin_plot.csv', header=T)
relaxframe <- read.csv('pristimantis_relax_results.csv', stringsAsFactors = F, header = T)
combframe <- merge(violin_annots, relaxframe, by = "Orthogroup")
# foo <- 
GO_vec <- c()
GO_code_vec <- c()
for (j in 1:nrow(combframe)) {
  if (!is.na(combframe$pfams_by_foo[j])) {
    GO_vec <- append(GO_vec, unlist(strsplit(as.character(combframe$pfams_by_foo[j]), split = ";", fixed = TRUE)))
    GO_code_vec <- append(GO_code_vec, unlist(strsplit(as.character(combframe$pdesc_by_foo[j]), split = ";", fixed = TRUE)))
  }
}
selparam_list <- list()
GO_frame <- data.frame("GO_category" = GO_code_vec, "GO_code" = GO_vec, "pval_median" = NA, "selparam_median" = NA,
                       "GO_frequency" = NA)
for (k in 1:length(GO_vec)) {
  pval_median_vec <- c()
  selparam_median_vec <- c()
  print(k)
  for (l in 1:nrow(combframe)) {
    invec <- FALSE
    invec <- grepl(GO_vec[k], combframe$pfams_by_foo[l], fixed = TRUE)
    if (invec) {
      pval_median_vec <- append(pval_median_vec, combframe$relax_pValue[l])
      selparam_median_vec <- append(selparam_median_vec, combframe$kValue[l])
    }
  }
  if (!is.null(pval_median_vec)) {
    pval_median <- median(pval_median_vec)
    GO_frame$pval_median[k] <- pval_median
  } else {
    GO_frame$pval_median[k] <- NA
  }
  if (!is.null(selparam_median_vec)) {
    selparam_median <- median(selparam_median_vec)
    GO_frame$selparam_median[k] <- selparam_median
    GO_frame$GO_frequency[k] <- length(selparam_median_vec)
  } else {
    GO_frame$selparam_median[k] <- NA
    GO_frame$GO_frequency[k] <- NA
  }
  selparam_list[[GO_code_vec[k]]] <- selparam_median_vec
}
GO_frame <- GO_frame[!duplicated(GO_frame), ]
#write.table(GO_frame, file = "RELAX_and_GO_analysis_out_table.csv", quote = FALSE, sep = ";",
            #row.names = FALSE)

#save(selparam_list, file = "RELAX_and_GO_analysis_session.RData")




GO_frame <- GO_frame[!duplicated(GO_frame[c(3,5)]),]
GO_frame %>% arrange(desc(GO_frequency))

wc.frame <- data.frame("GO_ID" = names(selparam_list), "psig" = NA, "qsig" = NA, "pval" = NA, "qval" = NA, "Vstat" = NA)
for (o in 1:length(selparam_list)) {
  pval <- NA
  qval <- NA
  psig <- NA
  qsig <- NA
  Vstat <- NA
  wc.test <- wilcox.test(selparam_list[[o]], mu = 1, alternative = "two.sided", na.action = na.omit)
  pval <- wc.test$p.value
  qval <- p.adjust(pval, method = "fdr", n = length(selparam_list))
  Vstat <- wc.test$statistic
  if (!is.na(pval)) {
    if (pval < 0.05) {
      psig <- TRUE
    } else {
      psig <- FALSE
    }
  }
  if (!is.na(qval)) {
    if (qval < 0.10) {
      qsig <- TRUE
    } else {
      qsig <- FALSE
    }
  }
  wc.frame$psig[o] <- psig
  wc.frame$qsig[o] <- qsig
  wc.frame$pval[o] <- pval
  wc.frame$qval[o] <- qval
  wc.frame$Vstat[o] <- Vstat
}
wc.frame<- wc.frame[!duplicated(wc.frame[c(4,5,6)]),]
violin_list <- selparam_list[which(wc.frame$pval < 1)]

k_vec <- c()
GO_term_vec <- c()
for (m in 1:length(violin_list)) {
  for (n in 1:length(violin_list[[m]])) {
    k_vec <- append(k_vec, violin_list[[m]][n])
    GO_term_vec <- append(GO_term_vec, names(violin_list)[m])
  }
}
violin_df <- data.frame("k" = k_vec, "GO" = factor(GO_term_vec))
violin_df <- violin_df[!duplicated(violin_df[c(1)]),]
vplot <- ggplot(violin_df, aes(GO, k, fill = GO))
vplot + geom_violin(draw_quantiles = c(0.25, 0.5, 0.75)) + geom_jitter(height = 0, width = 0.1) + guides(fill=FALSE) + geom_hline(yintercept=1, linetype = "dashed", 
                                                                                                                                  colour = "black", size=1) + coord_cartesian(ylim = c(0, 3)) +
  theme_grey(base_size = 18) + theme(axis.text.x = element_text(angle = 20, hjust = 1)) +
  scale_x_discrete(name = "GO term")
vplot + geom_violin(draw_quantiles = c(0.25, 0.5, 0.75)) + geom_jitter(height = 0, width = 0.1) + guides(fill=FALSE) + geom_hline(yintercept=1, linetype = "dashed", 
                                                                                                                                  colour = "black", size=1) + coord_cartesian(ylim = c(0, 3)) +
  theme_grey(base_size = 20)
#vplot + geom_boxplot() + geom_jitter(height = 0, width = 0.1) + guides(fill=FALSE) + geom_hline(yintercept=1, linetype = "dashed", 
wilcox.test(relaxframe$kValue, mu = 1, alternative = "two.sided", na.action = na.omit)

#Wilcoxon signed rank test with continuity correction

#data:  relaxframe$kValue
#V = 8772, p-value = 5.388e-06
#alternative hypothesis: true location is not equal to 1                                                                                                colour = "black", size=1) + coord_cartesian(ylim = c(0, 3))                                                                                              