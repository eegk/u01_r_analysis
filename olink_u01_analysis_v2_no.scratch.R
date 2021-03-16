#############################################################################################################################################
#############################################################################################################################################
libs<-c("corrplot","ROTS","Hmisc","tcR","tidyr","limma","NMF","tsne","Rtsne" ,"UpSetR",
        "AnnotationDbi","RColorBrewer","papmap","sva","GO.db","fitdistrplus","ff","plyranges",
        "annotables","Rsamtools","GenomicFeatures","ggarrange","pheatmap","rms",
        "dplyr","DESeq2","ggplot2","ggrepel","edgeR","doParallel","ggradar","colorspace",
        "reshape2","rmarkdown","org.Hs.eg.db","treemapify",
        "variancePartition","ggpubr","factoextra","Mfuzz", "universalmotif","ggbio","GenomicRanges",
        "scales","tibble","RColorBrewer","tidyr","ggplot2","reshape2","circlize",
        "colorspace","Vennerable","enrichR","cowplot","data.table",
        "ROCR","mlbench","caretEnsemble","gbm","rpart","UpSetR","SCENIC","AUCell","RcisTarget","plyr",
        "tidyverse","hrbrthemes","fmsb","colormap","viridis","survminer","survival","ggalluvial")
lapply(libs, require, character.only = TRUE) ; rm(libs)
#if (!requireNamespace("BiocManager", quietly = TRUE)) ; install.packages("BiocManager") ; BiocManager::install("survminer")
#############################################################################################################################################
setwd("/Users/gonzae34/Documents/projects_gnjatic/U01_Colitis_Project")
#############################################################################################################################################
### Olink data analysis for U01 Sacha Gnjatic project
### Initial analysis: Jan 04 2021
### Last review: Mar 15 2021
### Edgar E. G. Kozlova
#############################################################################################################################################
### flat corr function
flattenCorrMatrix <- function(cormat, pmat) { ut <- upper.tri(cormat)
  data.frame( row = rownames(cormat)[row(cormat)[ut]], column = rownames(cormat)[col(cormat)[ut]],cor  =(cormat)[ut], p = pmat[ut] ) }
#############################################################################################################################################
### Prepare the data
### WARNING: Data prep. is strictly dependent on the input file structure. Changes here may cause errors downstream
#############################################################################################################################################
### Load data from EXCEL
all_olink_data <- as.data.frame(readxl::read_excel("/Users/gonzae34/Documents/projects_gnjatic/U01_Colitis_Project/u01_olink_colitis_data/SAGN20_IO&Inflammation_All_Plates_NPX_Combined_COPY.xlsx"))
### Inspect excel file and separate
### First metadata
colitis_metadata <- all_olink_data[-c(191:193),c(185:206)]
### Then protein ids
protein_metadata <- all_olink_data[c(191:193),-c(185:206)]
### Finally only values
colitis_data_matrix <- all_olink_data[-c(191:193),-c(185:206)]
protein_metadata <- as.data.frame(t(protein_metadata))
protein_metadata$Marker <- rownames(protein_metadata)
### Transform to numeric
colitis_data_matrix <- lapply(colitis_data_matrix, function(x) { if(is.character(x)) as.numeric(as.character(x)) else x })
sapply(colitis_data_matrix, class) ; colitis_data_matrix <- do.call(cbind,colitis_data_matrix)
colitis_data_matrix <- as.data.frame(colitis_data_matrix)
### relabel metadata columns
colnames(colitis_metadata)[colnames(colitis_metadata) %in% c("QC Warning...190","Assay...195","Assay...198","Assay...199","Assay...200","Assay...202","Assay...203","Assay...204")] <- c("QC","Plate","Type","Patient","Timepoint","Fluid","Study","Colitis")
### modify few items
colitis_metadata$Colitis[colitis_metadata$Colitis %in% "Disease"] <- "Yes"
colitis_metadata$Colitis[colitis_metadata$Colitis %in% "Healthy"] <- "No"
### number the samples
colitis_metadata$sample_id <- paste("sample",1:nrow(colitis_metadata),sep="_")
### add the ids to each row
rownames(colitis_data_matrix) <- colitis_metadata$sample_id
### samples problematic, either pass or NA in one of the experiments
table(!is.na(colitis_data_matrix$TNFRSF9...2))
### remove NA samples
ix <- !is.na(colitis_data_matrix$TNFRSF9...2)
### data
colitis_metadata <- colitis_metadata[ix,]
colitis_data_matrix <- colitis_data_matrix[ix,]
### Observe data distribution and composition with a heatmap (small enough to visualize)
pdf(file="u01_olink_colitis_figures/u01_olink_colitis_raw_data_heatmap.pdf",width = 8,height = 8)
pheatmap(colitis_data_matrix,show_rownames = FALSE,show_colnames = FALSE,scale = "none",na_col = "white")
dev.off()
dev.off()
### transpose for tools
colitis_data_matrix <- as.data.frame(t(colitis_data_matrix))
###
table(colitis_metadata$Fluid,colitis_metadata$Study)
table(colitis_metadata$Colitis,colitis_metadata$Study)
###
colnames(protein_metadata)<-c("QC1","QC2","Experiment","Marker")
protein_metadata$Experiment_Short <- protein_metadata$Experiment
# table(protein_metadata$Experiment_Short)
protein_metadata$Experiment_Short[protein_metadata$Experiment_Short %in% "Olink Target 96 Immuno-Oncology(v.3112)"] <- "IO"
protein_metadata$Experiment_Short[protein_metadata$Experiment_Short %in% "Olink Target 96 Inflammation(v.3022)"] <- "I"
protein_metadata$Molecule_Experiment <- paste(protein_metadata$Experiment_Short,protein_metadata$Marker,sep="___")
###
colitis_metadata$Study_Type<- paste(colitis_metadata$Study,colitis_metadata$Type,sep="___")
colitis_metadata$Colitis_Timepoint <- paste(colitis_metadata$Colitis,colitis_metadata$Timepoint,sep="___")
### STORE CHECKPOINT DATA
### save.image(file="u01_RData_files/olink_experiment.RData")
#############################################################################################################################################
### variance analysis
#############################################################################################################################################
### Load checkpoint
load(file="u01_RData_files/olink_experiment.RData")
### Inspect data
dim(colitis_data_matrix) # 184 183
dim(colitis_metadata) # 183 25
dim(protein_metadata) # 184 6
### Select only immune panel (EXP 1)
ix <- which(protein_metadata$Experiment_Short %in% "I")
### NOTE1: Study and fluid are correlated
### NOTE2: Select only the first experiment
### SUBSET & Replace NAs with average value
colitis_data_matrix[ix,][is.na(colitis_data_matrix[ix,])] <- summary(c(unlist(colitis_data_matrix[ix,])))[3]
### 
form <- ~ (1|QC) + (1|Patient) + (1|Plate) + (1|Timepoint) + (1|Study) + (1|Colitis) 
variance_exprs_matrix <- fitExtractVarPartModel(colitis_data_matrix[ix,], form, colitis_metadata)
###  figure for variance
#############################################################################################################################################
pdf(file="u01_olink_colitis_figures/VP_olink_inflammation_experiment.pdf",width = 4, height = 4)
plotVarPart( sortCols( variance_exprs_matrix[,] , decreasing=FALSE) ) + 
  ggtitle('Olink Target 96 Inflammation(v.3022)') + 
  ylim(0,75) + ### Reduce the amount of white space
  theme(axis.text.x = element_text(angle = 45,hjust = 1)) + 
  coord_flip() + 
  ylab("Variance explained (%)") +
  theme(text = element_text(size=rel(4.8)),
        legend.text=element_text(size=rel(3.0)),
        axis.title.x = element_text(size=rel(3.0)),
        plot.title=element_text(size=rel(2)) ) 
dev.off()
#############################################################################################################################################
### Select EXP 2
ix <- which(protein_metadata$Experiment_Short %in% "IO")
### average replace
colitis_data_matrix[ix,][is.na(colitis_data_matrix[ix,])] <- summary(c(unlist(colitis_data_matrix[ix,])))[3]
### 
variance_exprs_matrix <- fitExtractVarPartModel(colitis_data_matrix[ix,], form, colitis_metadata)
###
#############################################################################################################################################
pdf(file="u01_olink_colitis_figures/VP_olink_IO_experiment.pdf",width = 4, height = 4)
plotVarPart( sortCols( variance_exprs_matrix , decreasing=FALSE) ) + 
  ggtitle('Olink Target 96 Immuno-Oncology(v.3112)') + ylim(0,100) + 
  theme(axis.text.x = element_text(angle = 45,hjust = 1)) + coord_flip() + 
  theme(text = element_text(size=rel(4.8)),
        legend.text=element_text(size=rel(3.0)),
        axis.title.x = element_text(size=rel(3.0)),
        plot.title=element_text(size=rel(2)) ) +
  ylab("Variance explained (%)")
dev.off()
#############################################################################################################################################
### Repeat including Type, which estimates the variance between the controls
### Select only immune panel (EXP 1)
ix <- which(protein_metadata$Experiment_Short %in% "I")
### NOTE1: Study and fluid are correlated
### NOTE2: Select only the first experiment
### Replace NAs with average value (few NAs, conversion is possible)
colitis_data_matrix[ix,][is.na(colitis_data_matrix[ix,])] <- summary(c(unlist(colitis_data_matrix[ix,])))[3]
### 
form <- ~ (1|QC) + (1|Type) + (1|Patient) + (1|Plate) + (1|Timepoint) + (1|Study) + (1|Colitis) 
variance_exprs_matrix <- fitExtractVarPartModel(colitis_data_matrix[ix,], form, colitis_metadata, )
###
#pdf(file="u01_olink_colitis_figures/VP_olink_inflammation_experiment_w_type.pdf",width = 4, height = 4)
plotVarPart( sortCols( variance_exprs_matrix[,] , decreasing=FALSE) ) + ggtitle('Olink Target 96 Inflammation(v.3022)') + ylim(0,75) + theme(axis.text.x = element_text(angle = 45,hjust = 1)) + coord_flip() + theme(plot.title = element_text(size=12)) + ylab("Variance explained (%)")
#dev.off()
###
### Select EXP 2
ix <- which(protein_metadata$Experiment_Short %in% "IO")
### average replace
colitis_data_matrix[ix,][is.na(colitis_data_matrix[ix,])] <- summary(c(unlist(colitis_data_matrix[ix,])))[3]
### 
variance_exprs_matrix <- fitExtractVarPartModel(colitis_data_matrix[ix,], form, colitis_metadata, )
###
#pdf(file="u01_olink_colitis_figures/VP_olink_IO_experiment_w_type.pdf",width = 4, height = 4)
plotVarPart( sortCols( variance_exprs_matrix , decreasing=FALSE) ) + ggtitle('Olink Target 96 Immuno-Oncology(v.3112)') + ylim(0,100) + theme(axis.text.x = element_text(angle = 45,hjust = 1)) + coord_flip() + theme(plot.title = element_text(size=12)) + ylab("Variance explained (%)")
#dev.off()
###
### clean up space
rm(variance_exprs_matrix,text,ix,form,test)
###
#############################################################################################################################################
### Profiling & Clustering
#############################################################################################################################################
### +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
### Select 1 experiment
### +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
ix <- which(protein_metadata$Experiment_Short %in% "I")
### subset
colitis_data_matrix_I <- colitis_data_matrix[ix,]
### select IPC
ix <- which(colitis_metadata$Patient %in% "IPC")
### loop correlations
my_matrix <- colitis_data_matrix_I[,ix]
the_rest <- colnames(colitis_data_matrix_I[,ix])
rm(my_storage) ; my_storage <- list() ; count=1
for (xxi in 1:length(the_rest) ){
  for( xxy in 1:length(the_rest)){
    my_test <- cor.test(my_matrix[,the_rest[xxi]],my_matrix[,the_rest[xxy]],method="spearman")    
    my_storage[[count]] <- data.frame( p.value = my_test$p.value , S = my_test$statistic , Rho = my_test$estimate ,
                                       Sample_A = the_rest[xxi] , Sample_B = the_rest[xxy] )
    count=count+1  } }
### data process
my_storage <- do.call(rbind,my_storage)
my_storage$nP.value <- -log10(my_storage$p.value)
my_storage$nP.value[my_storage$nP.value > 10]  <- 10
### plot
#pdf(file="u01_olink_colitis_figures/correlations_internal_controls_i.pdf",width = 4, height = 3.5)
ggplot(my_storage, aes(Sample_A, Sample_B, color= Rho, size=nP.value )) + 
  geom_point(shape=16)+ geom_text(aes(label=round(Rho,2)),size=2,color="black") +
  scale_colour_gradient2(low='firebrick', mid="white",high='steelblue', name='Rho') + #,breaks = c(-2.8, 0, 1, 5)) + #scale_color_viridis() +
  theme_bw() +rotate_x_text(angle=45) + 
  labs(x ='', y='', title='IPCs-I') 
#dev.off()
### +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
### Select 2 experiment
### +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
ix <- which(protein_metadata$Experiment_Short %in% "IO")
### subset
colitis_data_matrix_I <- colitis_data_matrix[ix,]
### select IPC
ix <- which(colitis_metadata$Patient %in% "IPC")
### loop correlations
my_matrix <- colitis_data_matrix_I[,ix]
the_rest <- colnames(colitis_data_matrix_I[,ix])
rm(my_storage) ; my_storage <- list() ; count=1
for (xxi in 1:length(the_rest) ){
  for( xxy in 1:length(the_rest)){
    my_test <- cor.test(my_matrix[,the_rest[xxi]],my_matrix[,the_rest[xxy]],method="spearman")    
    my_storage[[count]] <- data.frame( p.value = my_test$p.value , S = my_test$statistic , Rho = my_test$estimate ,
                                       Sample_A = the_rest[xxi] , Sample_B = the_rest[xxy] )
    count=count+1  } }
### data process
my_storage <- do.call(rbind,my_storage)
my_storage$nP.value <- -log10(my_storage$p.value)
my_storage$nP.value[my_storage$nP.value > 10]  <- 10
### plot
#pdf(file="u01_olink_colitis_figures/correlations_internal_controls_io.pdf",width = 4, height = 3.5)
ggplot(my_storage, aes(Sample_A, Sample_B, color= Rho, size=nP.value )) + 
  geom_point(shape=16)+ geom_text(aes(label=round(Rho,2)),size=2,color="black") +
  scale_colour_gradient2(low='firebrick', mid="white",high='steelblue', name='Rho') + #,breaks = c(-2.8, 0, 1, 5)) + #scale_color_viridis() +
  theme_bw() +rotate_x_text(angle=45) + 
  labs(x ='', y='', title='IPCs-IO') 
#dev.off()
### +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
### Select both experiment
### +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
### subset
colitis_data_matrix_I <- colitis_data_matrix
### select IPC
ix <- which(colitis_metadata$Patient %in% "IPC")
### loop correlations
my_matrix <- colitis_data_matrix_I[,ix]
the_rest <- colnames(colitis_data_matrix_I[,ix])
rm(my_storage) ; my_storage <- list() ; count=1
for (xxi in 1:length(the_rest) ){
  for( xxy in 1:length(the_rest)){
    my_test <- cor.test(my_matrix[,the_rest[xxi]],my_matrix[,the_rest[xxy]],method="spearman")    
    my_storage[[count]] <- data.frame( p.value = my_test$p.value , S = my_test$statistic , Rho = my_test$estimate ,
                                       Sample_A = the_rest[xxi] , Sample_B = the_rest[xxy] )
    count=count+1  } }
### data process
my_storage <- do.call(rbind,my_storage)
my_storage$nP.value <- -log10(my_storage$p.value)
my_storage$nP.value[my_storage$nP.value > 10]  <- 10
### plot
#pdf(file="u01_olink_colitis_figures/correlations_internal_controls_io+i.pdf",width = 4, height = 3.5)
ggplot(my_storage, aes(Sample_A, Sample_B, color= Rho, size=nP.value )) + 
  geom_point(shape=16)+ geom_text(aes(label=round(Rho,2)),size=2,color="black") +
  scale_colour_gradient2(low='firebrick', mid="white",high='steelblue', name='Rho') + #,breaks = c(-2.8, 0, 1, 5)) + #scale_color_viridis() +
  theme_bw() +rotate_x_text(angle=45) + 
  labs(x ='', y='', title='IPCs-I') 
#dev.off()
### +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
### All samples in experiment 1 (samples vs samples)
### +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
ix <- which(protein_metadata$Experiment_Short %in% "I")
### subset
colitis_data_matrix_I <- colitis_data_matrix[ix,]
### loop correlations
my_matrix <- colitis_data_matrix_I
the_rest <- colnames(colitis_data_matrix_I)
rm(my_storage) ; my_storage <- list() ; count=1
for (xxi in 1:length(the_rest) ){
  for( xxy in 1:length(the_rest)){
    my_test <- cor.test(my_matrix[,the_rest[xxi]],my_matrix[,the_rest[xxy]],method="spearman")    
    my_storage[[count]] <- data.frame( p.value = my_test$p.value , S = my_test$statistic , Rho = my_test$estimate ,
                                       Sample_A = the_rest[xxi] , Sample_B = the_rest[xxy] )
    count=count+1  } }
### data process
my_storage <- do.call(rbind,my_storage)
my_storage$nP.value <- -log10(my_storage$p.value)
my_storage$nP.value[my_storage$nP.value > 10]  <- 10
###
head(my_storage)
### arrange
rho.matrix <- dcast(data = my_storage,formula = Sample_A~Sample_B,value.var = "Rho")
rownames(rho.matrix) <- rho.matrix[,1]
rho.matrix <- rho.matrix[,-1]
###
pval.matrix <- dcast(data = my_storage,formula = Sample_A~Sample_B,value.var = "p.value")
rownames(pval.matrix) <- pval.matrix[,1]
pval.matrix <- pval.matrix[,-1]
###
ix <- match( colitis_metadata$sample_id,colnames(rho.matrix))
###
rho.matrix <- rho.matrix[ix,ix]
pval.matrix <- pval.matrix[ix,ix]
### verify
identical(colitis_metadata$sample_id,rownames(rho.matrix))
### subset metadata
ix <- which(colnames(colitis_metadata) %in% c("QC","Plate","Timepoint","Study","Type","Fluid","Colitis"))
annotation_data <- colitis_metadata[,ix]
rownames(annotation_data) <- colitis_metadata$sample_id
###
pdf(file="u01_olink_colitis_figures/correlations_panel_i_samples_x_samples.pdf",width = 12, height = 10)
pheatmap( as.matrix(rho.matrix),
          main="I-Panel Spearman's Rho",
          show_rownames = FALSE,
          show_colnames = FALSE,
          scale="none",
          annotation_col = annotation_data)
dev.off()
### store
my_correlations_results <- list( raw_I_rho = rho.matrix , raw_I_pval = pval.matrix, raw_I_flat = my_storage )
### +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
### All samples in experiment 2
### +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
ix <- which(protein_metadata$Experiment_Short %in% "IO")
### subset
colitis_data_matrix_I <- colitis_data_matrix[ix,]
### loop correlations
my_matrix <- colitis_data_matrix_I
the_rest <- colnames(colitis_data_matrix_I)
rm(my_storage) ; my_storage <- list() ; count=1
for (xxi in 1:length(the_rest) ){
  for( xxy in 1:length(the_rest)){
    my_test <- cor.test(my_matrix[,the_rest[xxi]],my_matrix[,the_rest[xxy]],method="spearman")    
    my_storage[[count]] <- data.frame( p.value = my_test$p.value , S = my_test$statistic , Rho = my_test$estimate ,
                                       Sample_A = the_rest[xxi] , Sample_B = the_rest[xxy] )
    count=count+1  } }
### data process
my_storage <- do.call(rbind,my_storage)
my_storage$nP.value <- -log10(my_storage$p.value)
my_storage$nP.value[my_storage$nP.value > 10]  <- 10
###
head(my_storage)
### arrange
rho.matrix <- dcast(data = my_storage,formula = Sample_A~Sample_B,value.var = "Rho")
rownames(rho.matrix) <- rho.matrix[,1]
rho.matrix <- rho.matrix[,-1]
###
pval.matrix <- dcast(data = my_storage,formula = Sample_A~Sample_B,value.var = "p.value")
rownames(pval.matrix) <- pval.matrix[,1]
pval.matrix <- pval.matrix[,-1]
###
ix <- match( colitis_metadata$sample_id,colnames(rho.matrix))
###
rho.matrix <- rho.matrix[ix,ix]
pval.matrix <- pval.matrix[ix,ix]
### verify
identical(colitis_metadata$sample_id,rownames(rho.matrix))
### subset metadata
ix <- which(colnames(colitis_metadata) %in% c("QC","Plate","Timepoint","Study","Type","Fluid","Colitis"))
annotation_data <- colitis_metadata[,ix]
rownames(annotation_data) <- colitis_metadata$sample_id
###
pdf(file="u01_olink_colitis_figures/correlations_panel_io_samples_x_samples.pdf",width = 12, height = 10)
pheatmap( as.matrix(rho.matrix),
          main="IO-Panel Spearman's Rho",
          show_rownames = FALSE,
          show_colnames = FALSE,
          scale="none",
          annotation_col = annotation_data)
dev.off()
#############################################################################################################################################
### DE
#############################################################################################################################################

### +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
### All samples in experiment inflammation
### +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

### Empty list to store the results
my_full_mixlm_results <- list()

### subset
ix <- which(colnames(colitis_metadata) %in% c("QC","Plate","Timepoint","Study","Type","Fluid","Colitis","Patient","Colitis_Timepoint","Study_Type"))
annotation_data <- colitis_metadata[,ix]
### subset for inflammation panel
ix <- which(protein_metadata$Experiment_Short %in% "I")
colitis_data_matrix_I <- colitis_data_matrix[ix,]
### subset
ix <- which(! annotation_data$Study %in% "control")
colitis_data_matrix_I <- colitis_data_matrix_I[,ix]
annotation_data <- annotation_data[ix,]
### subset (remove warnings from experiment)
ix <- which(! annotation_data$QC %in% "Warning")
colitis_data_matrix_I <- colitis_data_matrix_I[,ix]
annotation_data <- annotation_data[ix,]

### full model (Very complex, considers study, timepoint and plate as covariates)
design <- model.matrix( ~ 0 + Colitis + Study_Type + Timepoint + Plate, data = annotation_data )
colnames(design)[1:2] <- c("No","Yes")
colnames(design) <- make.names(colnames(design))
# comparisons
contr.matrix <- makeContrasts( 
  "Colitis" = No - Yes,
  "C_vs_A" = TimepointC,
  levels = colnames(design) )
### STATs
block2 <- as.numeric(as.factor( annotation_data$Patient ))
dupcor2 <- duplicateCorrelation(colitis_data_matrix_I, design, block=block2)
vfit <- lmFit(colitis_data_matrix_I, design, block=block2, correlation=dupcor2$consensus)
vfit <- contrasts.fit(vfit, contrasts=contr.matrix)
efit <- eBayes(vfit, robust = TRUE) ; summary(decideTests(efit))

my_full_mixlm_results[[1]] <- topTable(efit, coef="Colitis", n=Inf, adjust.method="BH")
my_full_mixlm_results[[1]]$Contrast <- "No - Yes"
my_full_mixlm_results[[1]]$Data <- "All Studies"
my_full_mixlm_results[[1]]$Model  <- "~ 0 + Colitis + Study_Type + Timepoint + Plate"
my_full_mixlm_results[[1]]$Panel  <- "Inflammation"

my_full_mixlm_results[[2]] <- topTable(efit, coef="C_vs_A", n=Inf, adjust.method="BH")
my_full_mixlm_results[[2]]$Contrast <- "C - A"
my_full_mixlm_results[[2]]$Data <- "All Studies"
my_full_mixlm_results[[2]]$Model  <- "~ 0 + Colitis + Study_Type + Timepoint + Plate"
my_full_mixlm_results[[2]]$Panel  <- "Inflammation"

#summary(unlist(colitis_data_matrix_I[rownames(colitis_data_matrix_I) %in% "CXCL9...112",annotation_data$Timepoint %in% "A"]))
#summary(unlist(colitis_data_matrix_I[rownames(colitis_data_matrix_I) %in% "CXCL9...112",annotation_data$Timepoint %in% "C"]))#

### time point specific model ###
### Vert simple, compares colitis between both timepoints
design <- model.matrix( ~ 0 + Colitis_Timepoint , data = annotation_data )
colnames(design) <- gsub("Colitis_Timepoint","",colnames(design))
colnames(design) <- make.names(colnames(design))
### comparisons
contr.matrix <- makeContrasts( 
  "Colitis_A" = Yes___A - No___A,
  "Colitis_C" = Yes___C - No___C,
  levels = colnames(design) )
### STATs
block2 <- as.numeric(as.factor( annotation_data$Patient ))
dupcor2 <- duplicateCorrelation(colitis_data_matrix_I, design, block=block2)
vfit <- lmFit(colitis_data_matrix_I, design, block=block2, correlation=dupcor2$consensus)
vfit <- contrasts.fit(vfit, contrasts=contr.matrix)
efit <- eBayes(vfit, robust = TRUE) ; summary(decideTests(efit))

my_full_mixlm_results[[3]] <- topTable(efit, coef="Colitis_A", n=Inf, adjust.method="BH")
my_full_mixlm_results[[3]]$Contrast <- "Yes___A - No___A"
my_full_mixlm_results[[3]]$Data <- "All Studies"
my_full_mixlm_results[[3]]$Model  <- "~ 0 + Colitis_Timepoint"
my_full_mixlm_results[[3]]$Panel  <- "Inflammation"

my_full_mixlm_results[[4]] <- topTable(efit, coef="Colitis_C", n=Inf, adjust.method="BH")
my_full_mixlm_results[[4]]$Contrast <- "Yes___C - No___C"
my_full_mixlm_results[[4]]$Data <- "All Studies"
my_full_mixlm_results[[4]]$Model  <- "~ 0 + Colitis_Timepoint"
my_full_mixlm_results[[4]]$Panel  <- "Inflammation"

### Because there was no immmediate signal between colitis, we explore variance profiles in a complex model
### VP is used for the main covariates.
### form & run
form <- ~ (1|Study_Type) + (1|Patient) + (1|Plate) + (1|Timepoint) + (1|Colitis) 
variance_exprs_matrix <- fitExtractVarPartModel(colitis_data_matrix_I, form, annotation_data)
### Figure
pdf(file="u01_olink_colitis_figures/VP_olink_inflammation_experiment_all_covariates_model.pdf",width = 4, height = 4)
plotVarPart( sortCols( variance_exprs_matrix[,] , decreasing=FALSE) ) + 
  ggtitle('Olink Target 96 Inflammation(v.3022)') + 
  ylim(0,75) + 
  theme(axis.text.x = element_text(angle = 45,hjust = 1)) + 
  coord_flip() + 
  theme(plot.title = element_text(size=12)) + 
  ylab("Variance explained (%)") +
  theme(text = element_text(size=rel(4.8)),
        legend.text=element_text(size=rel(3.0)),
        axis.title.x = element_text(size=rel(3.0)),
        plot.title=element_text(size=rel(2)) )
dev.off()
###

my_full_mixlm_results[[5]] <- variance_exprs_matrix
my_full_mixlm_results[[5]]$Data <- "All Studies"
my_full_mixlm_results[[5]]$Model <- "~ (1|Study_Type) + (1|Patient) + (1|Plate) + (1|Timepoint) + (1|Colitis)"
my_full_mixlm_results[[5]]$Panel  <- "Inflammation"

### +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
### +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
### study correction using non-parametric bayes adjustment
### +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
library(sva)
Batch = as.numeric(as.factor(as.character(annotation_data$Study_Type)))
modcombat = model.matrix(~1, data=annotation_data)
colitis_data_matrix_I_study_corrected = ComBat(dat=colitis_data_matrix_I, batch=Batch, mod=modcombat, par.prior=TRUE, prior.plots=FALSE)
### +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
### +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

### time point specific model
design <- model.matrix( ~ 0 + Colitis_Timepoint , data = annotation_data )
colnames(design) <- gsub("Colitis_Timepoint","",colnames(design))
colnames(design) <- make.names(colnames(design))
### comparisons
contr.matrix <- makeContrasts( 
  "Colitis_A" = Yes___A - No___A,
  "Colitis_C" = Yes___C - No___C,
  levels = colnames(design) )
### STATs
block2 <- as.numeric(as.factor( annotation_data$Patient ))
dupcor2 <- duplicateCorrelation(colitis_data_matrix_I_study_corrected, design, block=block2)
vfit <- lmFit(colitis_data_matrix_I_study_corrected, design, block=block2, correlation=dupcor2$consensus)
vfit <- contrasts.fit(vfit, contrasts=contr.matrix)
efit <- eBayes(vfit, robust = TRUE) ; summary(decideTests(efit))

my_full_mixlm_results[[6]] <- topTable(efit, coef="Colitis_A", n=Inf, adjust.method="BH")
my_full_mixlm_results[[6]]$Contrast <- "Yes___A - No___A"
my_full_mixlm_results[[6]]$Data <- "Batch corrected for Studies"
my_full_mixlm_results[[6]]$Model  <- "~ 0 + Colitis_Timepoint"
my_full_mixlm_results[[6]]$Panel  <- "Inflammation"

i<-7
my_full_mixlm_results[[i]] <- topTable(efit, coef="Colitis_C", n=Inf, adjust.method="BH")
my_full_mixlm_results[[i]]$Contrast <- "Yes___C - No___C"
my_full_mixlm_results[[i]]$Data <- "Batch corrected for Studies"
my_full_mixlm_results[[i]]$Model  <- "~ 0 + Colitis_Timepoint"
my_full_mixlm_results[[i]]$Panel  <- "Inflammation"

### study correction
form <- ~ (1|Study_Type) + (1|Patient) + (1|Plate) + (1|Timepoint) + (1|Colitis) 
variance_exprs_matrix <- fitExtractVarPartModel(colitis_data_matrix_I_study_corrected, form, annotation_data)
###
pdf(file="u01_olink_colitis_figures/VP_olink_inflammation_experiment_all_covariates_model_CORRECTED_4_study.pdf",width = 4, height = 4)
plotVarPart( sortCols( variance_exprs_matrix[,] , decreasing=FALSE) ) + 
  ggtitle('Olink Target 96 Inflammation(v.3022)') + 
  ylim(0,75) + theme(axis.text.x = element_text(angle = 45,hjust = 1)) + 
  coord_flip() + 
  ylab("Variance explained (%)") +
  theme(text = element_text(size=rel(4.8)),
        legend.text=element_text(size=rel(3.0)),
        axis.title.x = element_text(size=rel(3.0)),
        plot.title=element_text(size=rel(2)) )
dev.off()
###
i<-8
my_full_mixlm_results[[i]] <- variance_exprs_matrix
my_full_mixlm_results[[i]]$Data <- "All Studies"
my_full_mixlm_results[[i]]$Model <- "~ (1|Study_Type) + (1|Patient) + (1|Plate) + (1|Timepoint) + (1|Colitis)"
my_full_mixlm_results[[i]]$Panel  <- "Inflammation"

### +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
### All samples in experiment inflammation/Oncology
### +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

### subset
ix <- which(colnames(colitis_metadata) %in% c("QC","Plate","Timepoint","Study","Type","Fluid","Colitis","Patient","Colitis_Timepoint","Study_Type"))
annotation_data <- colitis_metadata[,ix]
### subset for inflammation panel
ix <- which(protein_metadata$Experiment_Short %in% "IO")
colitis_data_matrix_I <- colitis_data_matrix[ix,]
### subset
ix <- which(! annotation_data$Study %in% "control")
colitis_data_matrix_I <- colitis_data_matrix_I[,ix]
annotation_data <- annotation_data[ix,]
### subset (remove warnings from experiment)
ix <- which(! annotation_data$QC %in% "Warning")
colitis_data_matrix_I <- colitis_data_matrix_I[,ix]
annotation_data <- annotation_data[ix,]

### full model
design <- model.matrix( ~ 0 + Colitis + Study_Type + Timepoint + Plate, data = annotation_data )
colnames(design)[1:2] <- c("No","Yes")
colnames(design) <- make.names(colnames(design))
# comparisons
contr.matrix <- makeContrasts( 
  "Colitis" = No - Yes,
  "C_vs_A" = TimepointC,
  levels = colnames(design) )
### STATs
block2 <- as.numeric(as.factor( annotation_data$Patient ))
dupcor2 <- duplicateCorrelation(colitis_data_matrix_I, design, block=block2)
vfit <- lmFit(colitis_data_matrix_I, design, block=block2, correlation=dupcor2$consensus)
vfit <- contrasts.fit(vfit, contrasts=contr.matrix)
efit <- eBayes(vfit, robust = TRUE) ; summary(decideTests(efit))

i<-9
my_full_mixlm_results[[i]] <- topTable(efit, coef="Colitis", n=Inf, adjust.method="BH")
my_full_mixlm_results[[i]]$Contrast <- "No - Yes"
my_full_mixlm_results[[i]]$Data <- "All Studies"
my_full_mixlm_results[[i]]$Model  <- "~ 0 + Colitis + Study_Type + Timepoint + Plate"
my_full_mixlm_results[[i]]$Panel  <- "IO"
i<-10
my_full_mixlm_results[[i]] <- topTable(efit, coef="C_vs_A", n=Inf, adjust.method="BH")
my_full_mixlm_results[[i]]$Contrast <- "C - A"
my_full_mixlm_results[[i]]$Data <- "All Studies"
my_full_mixlm_results[[i]]$Model  <- "~ 0 + Colitis + Study_Type + Timepoint + Plate"
my_full_mixlm_results[[i]]$Panel  <- "IO"
#summary(unlist(colitis_data_matrix_I[rownames(colitis_data_matrix_I) %in% "CXCL9...112",annotation_data$Timepoint %in% "A"]))
#summary(unlist(colitis_data_matrix_I[rownames(colitis_data_matrix_I) %in% "CXCL9...112",annotation_data$Timepoint %in% "C"]))#

### time point specific model
design <- model.matrix( ~ 0 + Colitis_Timepoint , data = annotation_data )
colnames(design) <- gsub("Colitis_Timepoint","",colnames(design))
colnames(design) <- make.names(colnames(design))
### comparisons
contr.matrix <- makeContrasts( 
  "Colitis_A" = Yes___A - No___A,
  "Colitis_C" = Yes___C - No___C,
  levels = colnames(design) )
### STATs
block2 <- as.numeric(as.factor( annotation_data$Patient ))
dupcor2 <- duplicateCorrelation(colitis_data_matrix_I, design, block=block2)
vfit <- lmFit(colitis_data_matrix_I, design, block=block2, correlation=dupcor2$consensus)
vfit <- contrasts.fit(vfit, contrasts=contr.matrix)
efit <- eBayes(vfit, robust = TRUE) ; summary(decideTests(efit))

i<-11
my_full_mixlm_results[[i]] <- topTable(efit, coef="Colitis_A", n=Inf, adjust.method="BH")
my_full_mixlm_results[[i]]$Contrast <- "Yes___A - No___A"
my_full_mixlm_results[[i]]$Data <- "All Studies"
my_full_mixlm_results[[i]]$Model  <- "~ 0 + Colitis_Timepoint"
my_full_mixlm_results[[i]]$Panel  <- "IO"

i<-12
my_full_mixlm_results[[i]] <- topTable(efit, coef="Colitis_C", n=Inf, adjust.method="BH")
my_full_mixlm_results[[i]]$Contrast <- "Yes___C - No___C"
my_full_mixlm_results[[i]]$Data <- "All Studies"
my_full_mixlm_results[[i]]$Model  <- "~ 0 + Colitis_Timepoint"
my_full_mixlm_results[[i]]$Panel  <- "IO"

###
form <- ~ (1|Study_Type) + (1|Patient) + (1|Plate) + (1|Timepoint) + (1|Colitis) 
variance_exprs_matrix <- fitExtractVarPartModel(colitis_data_matrix_I, form, annotation_data)
###
pdf(file="u01_olink_colitis_figures/VP_olink_IO_experiment_all_covariates_model.pdf",width = 4, height = 4)
plotVarPart( sortCols( variance_exprs_matrix[,] , decreasing=FALSE) ) + 
  ggtitle('Olink Target 96 Immuno-Oncology(v.3112)') + ylim(0,75) + 
  theme(axis.text.x = element_text(angle = 45,hjust = 1)) + 
  coord_flip() + 
  ylab("Variance explained (%)") +
  theme(text = element_text(size=rel(4.8)),
        legend.text=element_text(size=rel(3.0)),
        axis.title.x = element_text(size=rel(3.0)),
        plot.title=element_text(size=rel(2)) )
dev.off()
###

i<-13
my_full_mixlm_results[[i]] <- variance_exprs_matrix
my_full_mixlm_results[[i]]$Data <- "All Studies"
my_full_mixlm_results[[i]]$Model <- "~ (1|Study_Type) + (1|Patient) + (1|Plate) + (1|Timepoint) + (1|Colitis)"
my_full_mixlm_results[[i]]$Panel  <- "IO"

### study correction
library(sva)
Batch = as.numeric(as.factor(as.character(annotation_data$Study_Type)))
modcombat = model.matrix(~1, data=annotation_data)
colitis_data_matrix_I_study_corrected = ComBat(dat=colitis_data_matrix_I, batch=Batch, mod=modcombat, par.prior=TRUE, prior.plots=FALSE)
###

### time point specific model
design <- model.matrix( ~ 0 + Colitis_Timepoint , data = annotation_data )
colnames(design) <- gsub("Colitis_Timepoint","",colnames(design))
colnames(design) <- make.names(colnames(design))
### comparisons
contr.matrix <- makeContrasts( 
  "Colitis_A" = Yes___A - No___A,
  "Colitis_C" = Yes___C - No___C,
  levels = colnames(design) )
### STATs
block2 <- as.numeric(as.factor( annotation_data$Patient ))
dupcor2 <- duplicateCorrelation(colitis_data_matrix_I_study_corrected, design, block=block2)
vfit <- lmFit(colitis_data_matrix_I_study_corrected, design, block=block2, correlation=dupcor2$consensus)
vfit <- contrasts.fit(vfit, contrasts=contr.matrix)
efit <- eBayes(vfit, robust = TRUE) ; summary(decideTests(efit))

i<-14
my_full_mixlm_results[[i]] <- topTable(efit, coef="Colitis_A", n=Inf, adjust.method="BH")
my_full_mixlm_results[[i]]$Contrast <- "Yes___A - No___A"
my_full_mixlm_results[[i]]$Data <- "Batch corrected for Studies"
my_full_mixlm_results[[i]]$Model  <- "~ 0 + Colitis_Timepoint"
my_full_mixlm_results[[i]]$Panel  <- "IO"

i<-15
my_full_mixlm_results[[i]] <- topTable(efit, coef="Colitis_C", n=Inf, adjust.method="BH")
my_full_mixlm_results[[i]]$Contrast <- "Yes___C - No___C"
my_full_mixlm_results[[i]]$Data <- "Batch corrected for Studies"
my_full_mixlm_results[[i]]$Model  <- "~ 0 + Colitis_Timepoint"
my_full_mixlm_results[[i]]$Panel  <- "IO"

### study correction
form <- ~ (1|Study_Type) + (1|Patient) + (1|Plate) + (1|Timepoint) + (1|Colitis) 
variance_exprs_matrix <- fitExtractVarPartModel(colitis_data_matrix_I_study_corrected, form, annotation_data)
###
pdf(file="u01_olink_colitis_figures/VP_olink_IO_experiment_all_covariates_model_CORRECTED_4_study.pdf",width = 4, height = 4)
plotVarPart( sortCols( variance_exprs_matrix[,] , decreasing=FALSE) ) + 
  ggtitle('Olink Target 96 Immuno-Oncology(v.3112)') + ylim(0,75) + 
  theme(axis.text.x = element_text(angle = 45,hjust = 1)) + 
  coord_flip() + 
  ylab("Variance explained (%)") +
  theme(text = element_text(size=rel(4.8)),
        legend.text=element_text(size=rel(3.0)),
        axis.title.x = element_text(size=rel(3.0)),
        plot.title=element_text(size=rel(2)) )
dev.off()
###
i<-16
my_full_mixlm_results[[i]] <- variance_exprs_matrix
my_full_mixlm_results[[i]]$Data <- "All Studies"
my_full_mixlm_results[[i]]$Model <- "~ (1|Study_Type) + (1|Patient) + (1|Plate) + (1|Timepoint) + (1|Colitis)"
my_full_mixlm_results[[i]]$Panel  <- "IO"

#############################################################################################################################################
### +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
### +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
### REVIEWED STATUS 15 March 2021
### All samples in experiment 1 ###  per study 
### To evalute each study independently we write a for loop with a series of comparisons to be done PER STUDY
###
### +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
### +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#############################################################################################################################################

identical(colnames(colitis_data_matrix),colitis_metadata$sample_id)
table(colitis_metadata$Study_Type,colitis_metadata$Colitis)

my_studies <- c("09-155","13-075, 17-162, 13-169","CTLA-4")
my_experiment <- levels(as.factor(as.character(protein_metadata$Experiment_Short)))
my_full_mixlm_results_SL <- list()
count=1

for (ixx in 1:length(my_experiment)) {
###
ix <- which(protein_metadata$Experiment_Short %in% my_experiment[ixx] )
colitis_data_matrix_I <- colitis_data_matrix[ix,]
###
ix <- which(colnames(colitis_metadata) %in% c("QC","Plate","Timepoint","Study","Study_Type","Type","Fluid","Colitis","sample_id","Patient"))
annotation_data <- colitis_metadata[,ix]
###
for (iyy in 1:length(my_studies)) {
###
ix <- which(annotation_data$Study %in% my_studies[iyy]) 
annotation_data_v2 <- annotation_data[ix,]
colitis_data_matrix_I_v2 <- colitis_data_matrix_I[,ix]
###
print(identical(colnames(colitis_data_matrix_I_v2),annotation_data_v2$sample_id))
print("DE run 1")
### Model
design <- model.matrix( ~ 0 + Colitis + Timepoint, data = annotation_data_v2 )
colnames(design)[1:2] <- c("No","Yes")
colnames(design) <- make.names(colnames(design))
### Comp
contr.matrix <- makeContrasts( 
  "Colitis" = No - Yes,
  "B_vs_A" = TimepointB,
  "C_vs_A" = TimepointC,
  "D_vs_A" = TimepointD,
  levels = colnames(design) )
### STATs
block2 <- as.numeric(as.factor( annotation_data_v2$Patient ))
dupcor2 <- duplicateCorrelation(colitis_data_matrix_I_v2, design, block=block2)
vfit <- lmFit(colitis_data_matrix_I_v2, design, block=block2, correlation=dupcor2$consensus)
vfit <- contrasts.fit(vfit, contrasts=contr.matrix)
efit <- eBayes(vfit, robust = TRUE) ; summary(decideTests(efit))
###
my_comparisons <- colnames(summary(decideTests(efit)))
print(count)
###
for ( iii in 1:length(my_comparisons)){
my_full_mixlm_results_SL[[count]] <- topTable(efit, coef=my_comparisons[iii], n=Inf, adjust.method="BH")
my_full_mixlm_results_SL[[count]]$Contrast <- my_comparisons[iii]
my_full_mixlm_results_SL[[count]]$Data <- my_studies[iyy]
my_full_mixlm_results_SL[[count]]$Model  <- "~ 0 + Colitis + Timepoint"
my_full_mixlm_results_SL[[count]]$Panel  <- my_experiment[ixx]
count <- count+1 }
print("DE run 2")
###
annotation_data_v2$Colitis_Timepoint <- paste(annotation_data_v2$Colitis,annotation_data_v2$Timepoint,sep="___")
design <- model.matrix( ~ 0 + Colitis_Timepoint, data = annotation_data_v2 )
colnames(design) <- gsub("Colitis_Timepoint","",colnames(design))
colnames(design) <- make.names(colnames(design))
###
contr.matrix <- makeContrasts( 
  "Colitis_A" = Yes___A - No___A,
  "Colitis_C" = Yes___C - No___C,
  levels = colnames(design) )
###
block2 <- as.numeric(as.factor( annotation_data_v2$Patient ))
dupcor2 <- duplicateCorrelation(colitis_data_matrix_I_v2, design, block=block2)
vfit <- lmFit(colitis_data_matrix_I_v2, design, block=block2, correlation=dupcor2$consensus)
vfit <- contrasts.fit(vfit, contrasts=contr.matrix)
efit <- eBayes(vfit, robust = TRUE) ; summary(decideTests(efit))
###
my_comparisons <- colnames(summary(decideTests(efit)))
###
print(count)
for ( iii in 1:length(my_comparisons)){
  my_full_mixlm_results_SL[[count]] <- topTable(efit, coef=my_comparisons[iii], n=Inf, adjust.method="BH")
  my_full_mixlm_results_SL[[count]]$Contrast <- my_comparisons[iii]
  my_full_mixlm_results_SL[[count]]$Data <- my_studies[iyy]
  my_full_mixlm_results_SL[[count]]$Model  <- "~ 0 + Colitis_Timepoint"
  my_full_mixlm_results_SL[[count]]$Panel  <- my_experiment[ixx]
  count <- count+1 }
###
print("parvar")
print(count)
### full model
form <- ~ (1|Patient) + (1|Plate) + (1|Timepoint) + (1|Colitis) 
variance_exprs_matrix <- fitExtractVarPartModel(colitis_data_matrix_I_v2, form, annotation_data_v2)
### plot
my_full_mixlm_results_SL[[count]] <- variance_exprs_matrix
my_full_mixlm_results_SL[[count]]$Data <- my_studies[iyy]
my_full_mixlm_results_SL[[count]]$Model <- paste(form)
my_full_mixlm_results_SL[[count]]$Panel  <- my_experiment[ixx]
count=count+1
my_full_mixlm_results_SL[[count]] <- plotVarPart( sortCols( variance_exprs_matrix[,] , decreasing=FALSE) ) + ggtitle(paste(my_experiment[ixx],my_studies[iyy],sep="___")) + ylim(0,75) + theme(axis.text.x = element_text(angle = 45,hjust = 1)) + coord_flip() + theme(plot.title = element_text(size=12)) + ylab("Variance explained (%)")
count=count+1
#print(count)
}
}
### lapply(my_full_mixlm_results_SL, function(x) { dim(x) })

### label lists
names(my_full_mixlm_results_SL) <- rep(c(rep("DE1",4),rep("DE2",2),"varpar","vplot"),6)
names(my_full_mixlm_results) <- rep(c(rep("DE1",2),rep("DE2",2),"varpar",rep("DE2_batch",2),"varpar_batch"),2)

### clean up 

#############################################################################################################################################
### RESULTS AND FIGURES
#############################################################################################################################################

###
length(my_full_mixlm_results) ### Previous models 16
length(my_full_mixlm_results_SL) ### New Models 48
###

### add scores and labels to these objects
my_numbers <- grep("DE",names(my_full_mixlm_results)) 
for (i in 1:length(my_numbers) ) {
  my_full_mixlm_results[[my_numbers[i]]]$nLog10FDR <- -log10(my_full_mixlm_results[[my_numbers[i]]]$adj.P.Val)
  my_full_mixlm_results[[my_numbers[i]]]$nLog10pvalue <- -log10(my_full_mixlm_results[[my_numbers[i]]]$P.Value)
  my_full_mixlm_results[[my_numbers[i]]]$Protein <- rownames(my_full_mixlm_results[[my_numbers[i]]])
  my_full_mixlm_results[[my_numbers[i]]]$Label <-  names(my_full_mixlm_results)[my_numbers[i]]
}

### add scores
my_numbers <- grep("DE",names(my_full_mixlm_results_SL))
for (i in 1:length(my_numbers) ) {
  my_full_mixlm_results_SL[[my_numbers[i]]]$nLog10FDR <- -log10(my_full_mixlm_results_SL[[my_numbers[i]]]$adj.P.Val)
  my_full_mixlm_results_SL[[my_numbers[i]]]$nLog10pvalue <- -log10(my_full_mixlm_results_SL[[my_numbers[i]]]$P.Value)
  my_full_mixlm_results_SL[[my_numbers[i]]]$Protein <- rownames(my_full_mixlm_results_SL[[my_numbers[i]]])
  my_full_mixlm_results_SL[[my_numbers[i]]]$Label <-  names(my_full_mixlm_results_SL)[my_numbers[i]]
}

### +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
### All data, model all co-variates
### +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

### Select only DE analysis
my_numbers <- grep("DE",names(my_full_mixlm_results)) 
### merge
data_for_figure <- do.call(rbind,my_full_mixlm_results[my_numbers])
data_for_figure$logFC[data_for_figure$Contrast %in% "No - Yes"] <- -data_for_figure$logFC[data_for_figure$Contrast %in% "No - Yes"]
data_for_figure$Contrast[data_for_figure$Contrast %in% "No - Yes"] <- "Colitis: Yes vs No"
data_for_figure$Contrast[data_for_figure$Contrast %in% "C - A"] <- "Time: C vs A"
data_for_figure$Contrast[data_for_figure$Contrast %in% "Yes___A - No___A"] <- "Colitis: Yes vs No (A)"
data_for_figure$Contrast[data_for_figure$Contrast %in% "Yes___C - No___C"] <- "Colitis: Yes vs No (C)"
data_for_figure$Significance <- "ND"
data_for_figure$Significance [data_for_figure$P.Value<0.05] <-"Nominal"
data_for_figure$Significance [data_for_figure$adj.P.Val<0.05] <-"FDR"
data_for_figure$Label[grep("batch",data_for_figure$Label)] <- "Batch"
data_for_figure$Label[grep("DE",data_for_figure$Label)] <- ""
data_for_figure$Panel[grep("Inflammation",data_for_figure$Panel)] <- "I"

data_for_figure$Comparisons <- paste(data_for_figure$Contrast,data_for_figure$Label,data_for_figure$Panel)
data_for_figure$Comparisons_Simple <- paste(data_for_figure$Contrast,data_for_figure$Label)

pdf(file="u01_olink_colitis_figures/all_data_models_DE_results_corrected_combined_models_PVALo5.pdf",width = 8,height = 12)
ggplot(data_for_figure[data_for_figure$P.Value<0.05,], aes(Protein, Comparisons_Simple, color= logFC, size=nLog10FDR, shape=Significance )) + 
  geom_point()+ #geom_point(shape=16)+ 
  geom_text(aes(label=round(logFC,1)),size=2,color="black") +
  scale_colour_gradient2(low='firebrick', mid="white",high='steelblue', name='LogFC') + #,breaks = c(-2.8, 0, 1, 5)) +
  theme_bw() +rotate_x_text(angle=25) + coord_flip()+ 
  labs(x ='', y='', title='') 
dev.off()

pdf(file="u01_olink_colitis_figures/all_data_models_DE_results_corrected_combined_models_FDRo5.pdf",width = 8,height = 12)
ggplot(data_for_figure[data_for_figure$adj.P.Val<0.05,], aes(Protein, Comparisons_Simple, color= logFC, size=nLog10FDR, shape=Significance )) + 
  geom_point()+ #geom_point(shape=16)+ 
  geom_text(aes(label=round(logFC,1)),size=2,color="black") +
  scale_colour_gradient2(low='firebrick', mid="white",high='steelblue', name='LogFC') + #,breaks = c(-2.8, 0, 1, 5)) +
  theme_bw() +rotate_x_text(angle=25) + coord_flip()+ 
  labs(x ='', y='', title='') 
dev.off()

#I
head(my_full_mixlm_results[[5]])
head(my_full_mixlm_results[[8]])
#IO
head(my_full_mixlm_results[[13]])
head(my_full_mixlm_results[[16]])

### top results for colitis
### top10

my_top_proteins <- list(
vp1=rownames(my_full_mixlm_results[[5]])[order(my_full_mixlm_results[[5]]$Colitis,decreasing = TRUE)] [1:5],
vp2=rownames(my_full_mixlm_results[[8]])[order(my_full_mixlm_results[[8]]$Colitis,decreasing = TRUE)] [1:5],
vp3=rownames(my_full_mixlm_results[[13]])[order(my_full_mixlm_results[[13]]$Colitis,decreasing = TRUE)] [1:5],
vp4=rownames(my_full_mixlm_results[[16]])[order(my_full_mixlm_results[[16]]$Colitis,decreasing = TRUE)] [1:5],
de=data_for_figure$Protein[data_for_figure$Contrast %in% c("Colitis: Yes vs No","Colitis: Yes vs No (A)","Colitis: Yes vs No (C)") & data_for_figure$P.Value<0.05] )

#data_for_figure[data_for_figure$Contrast %in% c("Colitis: Yes vs No","Colitis: Yes vs No (A)","Colitis: Yes vs No (C)") & data_for_figure$adj.P.Val<0.3,] # 1 hit

my_top_proteins$vp1 <- tidyr::separate(data.frame(my_top_proteins$vp1), 1, sep="\\.\\.\\.", c("a","b","c"))$a
my_top_proteins$vp2 <- tidyr::separate(data.frame(my_top_proteins$vp2), 1, sep="\\.\\.\\.", c("a","b","c"))$a
my_top_proteins$vp3 <- tidyr::separate(data.frame(my_top_proteins$vp3), 1, sep="\\.\\.\\.", c("a","b","c"))$a
my_top_proteins$vp4 <- tidyr::separate(data.frame(my_top_proteins$vp4), 1, sep="\\.\\.\\.", c("a","b","c"))$a
my_top_proteins$de <- tidyr::separate(data.frame(my_top_proteins$de), 1, sep="\\.\\.\\.", c("a","b","c"))$a


### Reduce(intersect,my_top_proteins[1:4]) ### 3 proteins ###  "TNF" "CXCL9" "CXCL10"
### All proteins significant for colitis at some point in time
my_top_proteins_global <- unique(sort(unlist(my_top_proteins)))
 
### +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
### study specific results
### +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

### selec the DE experiments
my_numbers <- grep("DE",names(my_full_mixlm_results_SL))

data_for_figure <- do.call(rbind,my_full_mixlm_results_SL[my_numbers])
data_for_figure$logFC[data_for_figure$Contrast %in% "Colitis"] <- -data_for_figure$logFC[data_for_figure$Contrast %in% "Colitis"]
data_for_figure$Contrast[data_for_figure$Contrast %in% "Colitis"] <- "Colitis: Yes - No (all)"
data_for_figure$Contrast[data_for_figure$Contrast %in% "Colitis_A"] <- "Colitis: Yes - No (A)"
data_for_figure$Contrast[data_for_figure$Contrast %in% "Colitis_C"] <- "Colitis: Yes - No (C)"
data_for_figure$Contrast[data_for_figure$Contrast %in% "B_vs_A"] <- "Time: (B) - (A) "
data_for_figure$Contrast[data_for_figure$Contrast %in% "C_vs_A"] <- "Time: (C) - (A)"
data_for_figure$Contrast[data_for_figure$Contrast %in% "D_vs_A"] <- "Time: (D) - (A)"
data_for_figure$Significance <- "ND"
data_for_figure$Significance [data_for_figure$P.Value<0.05] <-"Nominal"
data_for_figure$Significance [data_for_figure$adj.P.Val<0.05] <-"FDR"

data_for_figure$Comparisons <- paste(data_for_figure$Contrast,data_for_figure$Data,data_for_figure$Panel)
data_for_figure$Comparisons_Simple <- paste(data_for_figure$Contrast,data_for_figure$Data)
data_for_figure <- data_for_figure[data_for_figure$P.Value<0.05,]

### head(data_for_figure)
pdf(file="u01_olink_colitis_figures/individual_data_models_DE_results_PVALo5.pdf",width = 14,height = 15)
ggplot(data_for_figure[data_for_figure$P.Value<0.05,], aes(Protein, Comparisons_Simple, color= logFC, size=nLog10FDR, shape=Significance )) + 
  geom_point()+ #geom_point(shape=16)+ 
  geom_text(aes(label=round(logFC,1)),size=2,color="black") +
  scale_colour_gradient2(low='firebrick', mid="white",high='steelblue', name='LogFC') + #,breaks = c(-2.8, 0, 1, 5)) +
  theme_bw() +rotate_x_text(angle=25) + coord_flip()+ 
  labs(x ='', y='', title='') 
dev.off()

pdf(file="u01_olink_colitis_figures/individual_data_models_DE_results_FDRo5.pdf",width = 14,height = 15)
ggplot(data_for_figure[data_for_figure$adj.P.Val<0.05,], aes(Protein, Comparisons_Simple, color= logFC, size=nLog10FDR, shape=Significance )) + 
  geom_point()+ #geom_point(shape=16)+ 
  geom_text(aes(label=round(logFC,1)),size=2,color="black") +
  scale_colour_gradient2(low='firebrick', mid="white",high='steelblue', name='LogFC') + #,breaks = c(-2.8, 0, 1, 5)) +
  theme_bw() +rotate_x_text(angle=25) + coord_flip()+ 
  labs(x ='', y='', title='') 
dev.off()

### +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
### olink heat maps all tests
### +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

my_top_proteins_local <- data_for_figure$Protein[data_for_figure$Contrast %in% c("Colitis: Yes - No (A)","Colitis: Yes - No (all)","Colitis: Yes - No (C)") & data_for_figure$P.Value<0.03]
my_top_proteins_local <- tidyr::separate(data.frame(my_top_proteins_local), 1, sep="\\.\\.\\.", c("a","b","c"))$a
my_top_proteins_local <- unique(sort(my_top_proteins_local))

list(a=my_top_proteins_global,b=my_top_proteins_local)
Reduce(intersect,list(a=my_top_proteins_global,b=my_top_proteins_local))

my_top_proteins <- unique(sort(unlist(list(a=my_top_proteins_global,b=my_top_proteins_local))))

###

ix <- which(colnames(colitis_metadata) %in% c("QC","Timepoint","Study","Colitis"))
annotation_data <- colitis_metadata[,ix]
rownames(annotation_data) <- colitis_metadata$sample_id

ix <- which( ! annotation_data$Study %in% "control" )
### iy <- which(protein_metadata$Experiment_Short %in% "I")

colitis_data_matrix_I <- as.matrix(colitis_data_matrix[,ix])
annotation_data <- annotation_data[ix,]

### grep( "IL-15RA", rownames(colitis_data_matrix_I))
my_position<-list()
for( i in 1:length(my_top_proteins)){ my_position[[i]] <- grep( my_top_proteins[i], rownames(colitis_data_matrix_I))   }
my_position <- sort(unlist(my_position))

library(sva)
colitis_data_matrix_I <- as.matrix(colitis_data_matrix_I)
Batch = as.numeric(as.factor(as.character(annotation_data$Study)))
modcombat = model.matrix(~1, data=annotation_data)
colitis_data_matrix_I_study_corrected = ComBat(dat=colitis_data_matrix_I, batch=Batch, mod=modcombat, par.prior=TRUE, prior.plots=FALSE)

###
colitis_data_matrix_I_study_corrected <- colitis_data_matrix_I_study_corrected[my_position,] 
# comparisons
design <- model.matrix( ~ 0 + Colitis, data = annotation_data ) ; colnames(design) <- c("No","Yes")
contr.matrix <- makeContrasts(  "Colitis" = Yes - No, levels = colnames(design) )
vfit <- lmFit(colitis_data_matrix_I_study_corrected, design) ; vfit <- contrasts.fit(vfit, contrasts=contr.matrix)
efit <- eBayes(vfit, robust = TRUE) ; summary(decideTests(efit))
### final result
final_result <- topTable(efit, coef="Colitis", n=Inf, adjust.method="BH")
###
head(final_result)

###
final_result$Comparisons <- " "
final_result$Proteins <- final_result$ID
final_result <- final_result[ order(final_result$Proteins), ]
final_result$nLog10FDR <- -log10(final_result$adj.P.Val)
rownames(final_result) <- 1:nrow(final_result) 
final_result <- final_result[-c(9:10),] 

pdf(file="u01_olink_colitis_figures/norm_batch_corrected_data_models_DE_lms_results_FDRo5.pdf",width = 8,height = 4)
ggplot(final_result[which(final_result$adj.P.Val<0.05),], aes(Proteins, Comparisons, color= logFC, size=nLog10FDR )) + 
  geom_point()+ #geom_point(shape=16)+ 
  geom_text(aes(label=round(logFC,1)),size=2,color="black") +
  scale_colour_gradient2(low='firebrick', mid="white",high='steelblue', name='LogFC') + #,breaks = c(-2.8, 0, 1, 5)) +
  theme_bw() +rotate_x_text(angle=45) + 
  labs(x ='', y='', title="Colitis: Yes vs No") 
dev.off()

###
ix <- which(rownames(colitis_data_matrix_I_study_corrected) %in% final_result$ID[final_result$adj.P.Val<0.05 & final_result$logFC > 0.5])

out_matrix <- pheatmap( colitis_data_matrix_I_study_corrected[ix,],
          main="Study Corrected NPX values",
          show_rownames = TRUE,
          show_colnames = FALSE,
          scale="row",
          annotation_col = annotation_data)

annotation_data$HCluster <- as.character(cutree(out_matrix$tree_col,k=4))
annotation_data$HCluster <- paste("HCluster",annotation_data$HCluster,sep="_")

annotation_data$Study[annotation_data$Study %in% "09-155"] <- "Study1"
annotation_data$Study[annotation_data$Study %in% "13-075, 17-162, 13-169"] <- "Study2"
annotation_data$Study[annotation_data$Study %in% "CTLA-4"] <- "Study3"

# annotation_data <- annotation_data[order(annotation_data$HCluster),]
iy <- c( grep( "HCluster_2",annotation_data$HCluster),grep( "HCluster_1",annotation_data$HCluster),grep( "HCluster_3",annotation_data$HCluster),grep( "HCluster_4",annotation_data$HCluster) )
annotation_data <- annotation_data[iy,]
iy <- match( rownames(annotation_data),colnames(colitis_data_matrix_I_study_corrected ) )

annotation_data_col_save <- annotation_data

annotation_row <- protein_metadata
annotation_row <- annotation_row [ which(annotation_row$Marker %in% rownames(colitis_data_matrix_I_study_corrected[ix,]) ) , ]
annotation_row <- annotation_row[,c("QC1","Experiment")]  
annotation_row$QC1 <- round(as.numeric(as.character(annotation_row$QC1)),2)
annotation_row$Experiment[annotation_row$Experiment %in% "Olink Target 96 Inflammation(v.3022)"] <- "I"
annotation_row$Experiment[annotation_row$Experiment %in% "Olink Target 96 Immuno-Oncology(v.3112)"] <- "IO"
  
annotation_colors <- list ( HCluster = c( HCluster_1="#781C6DFF",HCluster_2="#000004FF",HCluster_3="#ED6925FF",HCluster_4="#FCFFA4FF"),
                              Colitis = c( No="steelblue",Yes="firebrick"),
                              Timepoint = c( A="grey20",B="grey40",C="grey60",D="grey80"),
                              Study = c( Study1="green",Study2="orange",Study3="purple"),
                              QC = c( Pass="black",Warning="red"),
                              Experiment = c(I="grey10",IO="grey50") )

colitis_data_matrix_I_study_corrected <- colitis_data_matrix_I_study_corrected[ix,iy]
colitis_data_matrix_I_study_corrected <- colitis_data_matrix_I_study_corrected[unique(rownames(colitis_data_matrix_I_study_corrected)),]

png(file="u01_olink_colitis_figures/pheatmap_clustered_significant_markers_colitis.png", width = 10, height = 7.5,units = "in", res = 300)
pheatmap( colitis_data_matrix_I_study_corrected,
          main="Study Corrected NPX values",
          cluster_cols = FALSE,
          show_rownames = TRUE,
          show_colnames = FALSE,
          scale="row",
          gaps_col = c(as.numeric(cumsum(table(annotation_data$HCluster)[c("HCluster_2","HCluster_1","HCluster_3","HCluster_4")]))) ,
          border_color=FALSE,
          fontsize=10,
          fontsize_row = 12,
          angle_col = 90,
          annotation_col = annotation_data,
          annotation_row = annotation_row,
          annotation_colors = annotation_colors )
dev.off()
dev.off()

### prepare for violin plots
my_melted_df <- melt(colitis_data_matrix_I_study_corrected) ### only the top proteins
colnames(my_melted_df) <- c("Protein","sample_id","value")
my_melted_df <- merge(my_melted_df,colitis_metadata,by="sample_id")
ix <- which(my_melted_df$Protein %in% "TNF...87")

## ggplot(data=my_melted_df) + aes(x=Colitis,y=value) + geom_boxplot() + facet_wrap(~Protein,scales = "free") + theme_bw()

my_melted_df$Protein <- as.character(my_melted_df$Protein)
my_melted_df$Protein <- tidyr::separate(data.frame(my_melted_df$Protein), 1, sep="\\.\\.\\.", c("a","b","c"))$a
my_melted_df$Protein <- as.factor(my_melted_df$Protein)
table(my_melted_df$Protein)
my_melted_df$Protein <- factor(my_melted_df$Protein, levels = c("CXCL10","CXCL9","IFN-gamma","TNF","MMP-10","IL-17C","4E-BP1"))  

library(ggbeeswarm)
my_comparisons <- list( c("Yes","No") )

#pdf(file="u01_olink_colitis_figures/boxplots_clustered_significant_markers_colitis.pdf ",width = 8,height = 7)
ggviolin(data=my_melted_df, x="Colitis", y="value", color="Colitis",fill="Colitis",alpha=0.5) + theme_bw() +
  facet_wrap(~Protein,scales = "free",nrow=2) +
  geom_boxplot(width=0.2) + geom_quasirandom(varwidth = TRUE,alpha=0.3,size=0.5) +
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test", paired=FALSE, label="p.signif", hide.ns = TRUE) +
  theme(axis.text.x = element_text(angle = 60,hjust = 1)) + theme(legend.position="bottom") + theme(axis.title.x=element_blank(),
                                                                                                    axis.text.x=element_blank(),
                                                                                                    axis.ticks.x=element_blank(),
                                                                                                    strip.text=element_text(size=rel(1.5),face="bold"),
                                                                                                    legend.text = element_text(size=rel(1.5),face="bold"),
                                                                                                    legend.title = element_text(size=rel(1.5),face="bold"),
                                                                                                    axis.text.y = element_text(size=rel(1.5),face="bold"),
                                                                                                    axis.title.y = element_text(size=rel(1.5),face="bold")) +
  scale_fill_manual(values = c("firebrick","steelblue")) +
  scale_color_manual(values = c("firebrick","steelblue")) +
  ggtitle("Study Corrected NPX values") + labs(x = "Colitis", y = "NPX") 
#dev.off()

my_melted_df$Status <- paste(my_melted_df$Timepoint,my_melted_df$Colitis,sep="___")
table(my_melted_df$Status)

my_melted_df$Status <- gsub("A_","Baseline_",my_melted_df$Status)
my_melted_df$Status <- gsub("D_","Post_",my_melted_df$Status)
my_melted_df$Status <- gsub("B_","Colitis_",my_melted_df$Status)
my_melted_df$Status <- gsub("C_","Week6_",my_melted_df$Status)

my_melted_df$Status <- factor(my_melted_df$Status, levels=c("Baseline___No","Baseline___Yes","Week6___No","Week6___Yes","Colitis___Yes","Post___Yes"))

my_comparisons <- list( c("Baseline___No","Baseline___Yes"),
                        c("Baseline___No","Week6___No"),
                        c("Baseline___Yes","Week6___Yes"),
                        c("Week6___No","Week6___Yes"),
                        c("Week6___No","Colitis___Yes"))

### Figure for markers from heatmap
### Model using adjusted data for study
### Top Markers after various models for DE
### +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
pdf(file="u01_olink_colitis_figures/boxplots_clustered_significant_markers_colitis_adjusted.model.pdf ",width = 18,height = 7)
ggviolin(data=my_melted_df, x="Status", y="value", color="Colitis",fill="Colitis",alpha=0.5) + theme_bw() +
  facet_wrap(~Protein,scales = "free",ncol=7) +
  geom_boxplot(width=0.2) + geom_quasirandom(varwidth = TRUE,alpha=0.3,size=0.5) +
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test", paired=FALSE, label="p.signif", hide.ns = FALSE) +
  theme(axis.text.x = element_text(angle = 45,hjust = 1)) + theme(legend.position="bottom") + theme(#axis.title.x=element_blank(),
                                                                                                    #axis.text.x=element_blank(),
                                                                                                    #axis.ticks.x=element_blank(),
                                                                                                    strip.text=element_text(size=rel(1.5),face="bold"),
                                                                                                    legend.text = element_text(size=rel(1.5),face="bold"),
                                                                                                    legend.title = element_text(size=rel(1.5),face="bold"),
                                                                                                    axis.text.y = element_text(size=rel(1.5),face="bold"),
                                                                                                    axis.title.y = element_text(size=rel(1.5),face="bold")) +
  scale_fill_manual(values = c("firebrick","steelblue")) +
  scale_color_manual(values = c("firebrick","steelblue")) +
  ggtitle("") + labs(x = "", y = "Study Adj. NPX values") +
  theme(text = element_text(size=rel(4.8)),
        legend.text = element_text(size=rel(4.8)),
        legend.title = element_text(size=rel(4.8)),
        plot.subtitle = element_text(size=rel(4.8)),
        axis.title.x = element_text(size=rel(4.8)),
        axis.text.y = element_text(size=rel(4.8)),
        strip.text = element_text(size=rel(4.8)),
        axis.title.y = element_text(size=rel(4.8)),
        plot.title=element_text(size=rel(4.8)) ) +
  theme(plot.margin=unit(c(t=5,r=5,b=5,l=100),"pt"))
dev.off()
### +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

### +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
### olink percent vars varpar
### +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

ix <- c(order(my_full_mixlm_results$varpar$Residuals)[1:5],
  order(my_full_mixlm_results$varpar$Colitis,decreasing = TRUE)[1:5],
  order(my_full_mixlm_results$varpar$Study_Type,decreasing = TRUE)[1:5],
  order(my_full_mixlm_results$varpar$Timepoint,decreasing = TRUE)[1:5] )
  
#pdf(file="u01_olink_colitis_figures/plotPercentBars_global_colitis.pdf ",width = 5,height = 8)
plotPercentBars( my_full_mixlm_results$varpar[ix,c(1:6)] )
#dev.off()

ix <- c(order(my_full_mixlm_results$varpar_batch$Residuals)[1:5],
        order(my_full_mixlm_results$varpar_batch$Colitis,decreasing = TRUE)[1:5],
        order(my_full_mixlm_results$varpar_batch$Study_Type,decreasing = TRUE)[1:5],
        order(my_full_mixlm_results$varpar_batch$Timepoint,decreasing = TRUE)[1:5] )

#pdf(file="u01_olink_colitis_figures/plotPercentBars_global_colitis_batch.pdf ",width = 5,height = 8)
plotPercentBars( my_full_mixlm_results$varpar_batch[ix,c(1:6)] )
#dev.off()

#############################################################################################################################################
### 
#############################################################################################################################################

#colitis_data_matrix
#colitis_metadata
#protein_metadata

### another checkpoint after analysis
### save.image("u01_RData_files/temporal_analysis.RData")
# load("u01_RData_files/temporal_analysis.RData")

#############################################################################################################################################
### DE between clusters defined by cytokines
#############################################################################################################################################

### corrected olink data
dim(colitis_data_matrix_I_study_corrected)

### order 2,1,3,4 for increase

ix <- which( colnames(colitis_data_matrix) %in% rownames(annotation_data) )
colitis_data_matrix_I <- as.matrix(colitis_data_matrix[,ix])

annotation_data <- annotation_data[ match( colnames(colitis_data_matrix_I),rownames(annotation_data) ) , ]
identical(colnames(colitis_data_matrix_I),rownames(annotation_data))

design <- model.matrix( ~ 0 + HCluster, data = annotation_data ) ; colnames(design) <- c("H1","H2","H3","H4")
contr.matrix <- makeContrasts(  "H_2v1" = H2 - H1,
                                "H_3v1" = H3 - H1,
                                "H_4v1" = H4 - H1,
                                "H_3v2" = H3 - H2,
                                "H_4v2" = H4 - H2,
                                "H_4v3" = H4 - H3,
                                levels = colnames(design) )
vfit <- lmFit(colitis_data_matrix_I_study_corrected, design) ; vfit <- contrasts.fit(vfit, contrasts=contr.matrix)
efit <- eBayes(vfit, robust = TRUE) ; summary(decideTests(efit))

### final result
de_hclusters <- topTable(efit, coef="Colitis", n=Inf, adjust.method="BH")

de_hclusters <- list()
for ( i in 1:ncol(summary(decideTests(efit))) ) { de_hclusters[[i]] <- topTable(efit, coef=i,n=Inf,adjust.method="BH") 
de_hclusters[[i]]$Comparison <- colnames(summary(decideTests(efit)))[i]
de_hclusters[[i]]$Protein <- rownames( topTable(efit, coef=i,n=Inf,adjust.method="BH") ) }
###
de_hclusters <- do.call(rbind,de_hclusters)
rownames(de_hclusters) <- NULL

### NOTE: this verifies that 4 of these proteins are indeed deferentially expressed between these clusters.
### The table is simple, thus can be stored for further analysis if PIs decide it's interesting.
###
de_hclusters_filtered <- de_hclusters[which(de_hclusters$adj.P.Val<0.05),]
###

#############################################################################################################################################
### Results breakdown FOR improved and complementary FIGURES
#############################################################################################################################################

### interesting proteins between local and global models
list(a=my_top_proteins_global,b=my_top_proteins_local)
Reduce(intersect,list(a=my_top_proteins_global,b=my_top_proteins_local))
my_top_proteins <- unique(sort(unlist(list(a=my_top_proteins_global,b=my_top_proteins_local))))

###
ix <- which(colnames(colitis_metadata) %in% c("QC","Timepoint","Study","Colitis"))
annotation_data <- colitis_metadata[,ix]
rownames(annotation_data) <- colitis_metadata$sample_id
ix <- which( ! annotation_data$Study %in% "control" )
### iy <- which(protein_metadata$Experiment_Short %in% "I")
colitis_data_matrix_I <- as.matrix(colitis_data_matrix[,ix])
annotation_data <- annotation_data[ix,]

my_position <- list()
for( i in 1:length(my_top_proteins)){ my_position[[i]] <- grep( my_top_proteins[i], rownames(colitis_data_matrix_I))   }
my_position <- unique(sort(unlist(my_position)))

library(sva)
colitis_data_matrix_II <- as.matrix(colitis_data_matrix_I)
Batch = as.numeric(as.factor(as.character(annotation_data$Study)))
modcombat = model.matrix(~1, data=annotation_data)
colitis_data_matrix_I_study_corrected = ComBat(dat=colitis_data_matrix_II, batch=Batch, mod=modcombat, par.prior=TRUE, prior.plots=FALSE)

out_matrix <- pheatmap( colitis_data_matrix_I_study_corrected[my_position,],
                        main="Study Corrected NPX values",
                        show_rownames = TRUE,
                        show_colnames = FALSE,
                        scale="row",
                        annotation_col = annotation_data)

annotation_data$HCluster <- as.character(cutree(out_matrix$tree_col,k=4))
annotation_data$HCluster <- paste("HCluster",annotation_data$HCluster,sep="_")

annotation_data$Study[annotation_data$Study %in% "09-155"] <- "Study1"
annotation_data$Study[annotation_data$Study %in% "13-075, 17-162, 13-169"] <- "Study2"
annotation_data$Study[annotation_data$Study %in% "CTLA-4"] <- "Study3"

# annotation_data <- annotation_data[order(annotation_data$HCluster),]
iy <- c( grep( "HCluster_2",annotation_data$HCluster),grep( "HCluster_1",annotation_data$HCluster),grep( "HCluster_3",annotation_data$HCluster),grep( "HCluster_4",annotation_data$HCluster) )
annotation_data <- annotation_data[iy,]
iy <- match( rownames(annotation_data),colnames(colitis_data_matrix_I_study_corrected ) )

annotation_row <- protein_metadata
annotation_row <- annotation_row [ which(annotation_row$Marker %in% rownames(colitis_data_matrix_I_study_corrected[,]) ) , ]
annotation_row <- annotation_row[,c("QC1","Experiment")]  
annotation_row$QC1 <- round(as.numeric(as.character(annotation_row$QC1)),2)
annotation_row$Experiment[annotation_row$Experiment %in% "Olink Target 96 Inflammation(v.3022)"] <- "I"
annotation_row$Experiment[annotation_row$Experiment %in% "Olink Target 96 Immuno-Oncology(v.3112)"] <- "IO"

annotation_colors <- list ( HCluster = c( HCluster_1="#781C6DFF",HCluster_2="#000004FF",HCluster_3="#ED6925FF",HCluster_4="#FCFFA4FF"),
                            Colitis = c( No="steelblue",Yes="firebrick"),
                            Timepoint = c( A="grey20",B="grey40",C="grey60",D="grey80"),
                            Study = c( Study1="green",Study2="orange",Study3="purple"),
                            QC = c( Pass="black",Warning="red"),
                            Experiment = c(I="grey10",IO="grey50") )

colitis_data_matrix_I_study_corrected <- colitis_data_matrix_I_study_corrected[,iy]

#png(file="u01_olink_colitis_figures/pheatmap_clustered_significant_markers_colitis_all_models_significant.png", width = 10, height = 4,units = "in", res = 300)
pheatmap( colitis_data_matrix_I_study_corrected[my_position,],
          main="Study Corrected NPX values",
          cluster_cols = FALSE,
          show_rownames = TRUE,
          show_colnames = FALSE,
          scale="row",
          gaps_col = c(as.numeric(cumsum(table(annotation_data$HCluster)[c("HCluster_2","HCluster_1","HCluster_3","HCluster_4")]))) ,
          border_color=FALSE,
          fontsize=6,
          angle_col = 90,
          annotation_col = annotation_data,
          annotation_row = annotation_row,
          annotation_colors = annotation_colors )
dev.off()
dev.off()

### prepare for violin plots
my_melted_df <- melt(colitis_data_matrix_I_study_corrected[my_position,])
colnames(my_melted_df) <- c("Protein","sample_id","value")
my_melted_df <- merge(my_melted_df,colitis_metadata,by="sample_id")
#ix <- which(my_melted_df$Protein %in% "TNF...87")
## ggplot(data=my_melted_df) + aes(x=Colitis,y=value) + geom_boxplot() + facet_wrap(~Protein,scales = "free") + theme_bw()

my_melted_df$Protein <- as.character(my_melted_df$Protein)
my_melted_df$Protein <- tidyr::separate(data.frame(my_melted_df$Protein), 1, sep="\\.\\.\\.", c("a","b","c"))$a
my_melted_df$Protein <- as.factor(my_melted_df$Protein)
table(my_melted_df$Protein)
### my_melted_df$Protein <- factor(my_melted_df$Protein, levels = c("CXCL10","CXCL9","IFN-gamma","TNF","MMP-10","IL-17C","CXCL1"))  

library(ggbeeswarm)
my_comparisons <- list( c("Yes","No") )

#pdf(file="u01_olink_colitis_figures/boxplots_clustered_significant_markers_colitis_all_models_significant.pdf ",width = 16,height = 14)
ggviolin(data=my_melted_df, x="Colitis", y="value", color="Colitis",fill="Colitis",alpha=0.5) + theme_bw() +
  facet_wrap(~Protein,scales = "free",ncol=7,nrow=4) +
  geom_boxplot(width=0.2) + geom_quasirandom(varwidth = TRUE,alpha=0.3,size=0.5) +
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test", paired=FALSE, label="p.signif", hide.ns = FALSE) +
  theme(axis.text.x = element_text(angle = 60,hjust = 1)) + theme(legend.position="bottom") + theme(axis.title.x=element_blank(),
                                                                                                    axis.text.x=element_blank(),
                                                                                                    axis.ticks.x=element_blank(),
                                                                                                    strip.text=element_text(size=rel(1.5),face="bold"),
                                                                                                    legend.text = element_text(size=rel(1.5),face="bold"),
                                                                                                    legend.title = element_text(size=rel(1.5),face="bold"),
                                                                                                    axis.text.y = element_text(size=rel(1.5),face="bold"),
                                                                                                    axis.title.y = element_text(size=rel(1.5),face="bold")) +
  scale_fill_manual(values = c("firebrick","steelblue")) +
  scale_color_manual(values = c("firebrick","steelblue")) +
  ggtitle("") + labs(x = "Colitis", y = "NPX") 
#dev.off()


my_melted_df$Status <- paste(my_melted_df$Timepoint,my_melted_df$Colitis,sep="___")
table(my_melted_df$Status)

my_melted_df$Status <- factor(my_melted_df$Status, levels=c("A___No","A___Yes","B___Yes","C___No","C___Yes","D___Yes"))

my_comparisons <- list( c("A___No","A___Yes"),
                        c("A___No","C___No"),
                        c("A___Yes","C___Yes"),
                        c("C___No","C___Yes"),
                        c("C___No","B___Yes"))

#pdf(file="u01_olink_colitis_figures/boxplots_clustered_significant_markers_colitis_time_all_models_significant.pdf ",width = 24,height = 14)
ggviolin(data=my_melted_df, x="Status", y="value", color="Colitis",fill="Colitis",alpha=0.5) + theme_bw() +
  facet_wrap(~Protein,scales = "free",ncol=7,nrow=4) +
  geom_boxplot(width=0.2) + geom_quasirandom(varwidth = TRUE,alpha=0.3,size=0.5) +
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test", paired=FALSE, label="p.signif", hide.ns = FALSE) +
  theme(axis.text.x = element_text(angle = 60,hjust = 1)) + theme(legend.position="bottom") + theme(#axis.title.x=element_blank(),
    #axis.text.x=element_blank(),
    #axis.ticks.x=element_blank(),
    strip.text=element_text(size=rel(1.5),face="bold"),
    legend.text = element_text(size=rel(1.5),face="bold"),
    legend.title = element_text(size=rel(1.5),face="bold"),
    axis.text.y = element_text(size=rel(1.5),face="bold"),
    axis.title.y = element_text(size=rel(1.5),face="bold")) +
  scale_fill_manual(values = c("firebrick","steelblue")) +
  scale_color_manual(values = c("firebrick","steelblue")) +
  ggtitle("") + labs(x = "Colitis", y = "Adj.NPX") 
#dev.off()

#############################################################################################################################################
### PD1 review
#############################################################################################################################################

### selection
ix <- which(colnames(colitis_metadata) %in% c("QC","Timepoint","Study","Colitis"))
annotation_data <- colitis_metadata[,ix]
rownames(annotation_data) <- colitis_metadata$sample_id
ix <- which( ! annotation_data$Study %in% "control" )
### iy <- which(protein_metadata$Experiment_Short %in% "I")
colitis_data_matrix_I <- as.matrix(colitis_data_matrix[,ix])
annotation_data <- annotation_data[ix,]

### prepare for violin plots
my_melted_df <- melt(colitis_data_matrix_I[,])
colnames(my_melted_df) <- c("Protein","sample_id","value")
my_melted_df <- merge(my_melted_df,colitis_metadata,by="sample_id")
#ix <- which(my_melted_df$Protein %in% "TNF...87")
## ggplot(data=my_melted_df) + aes(x=Colitis,y=value) + geom_boxplot() + facet_wrap(~Protein,scales = "free") + theme_bw()

my_melted_df$Protein <- as.character(my_melted_df$Protein)
my_melted_df$Protein <- tidyr::separate(data.frame(my_melted_df$Protein), 1, sep="\\.\\.\\.", c("a","b","c"))$a
my_melted_df$Protein <- as.factor(my_melted_df$Protein)
table(my_melted_df$Protein)

# grep("PD", my_melted_df$Protein,value=TRUE)
# data_for_figure ### Results from various mixed linear models.
my_melted_df$Status <- paste(my_melted_df$Timepoint,my_melted_df$Colitis,sep="___")
my_melted_df$Status <- factor(my_melted_df$Status, levels=c("A___No","A___Yes","B___Yes","C___No","C___Yes","D___Yes"))

my_comparisons <- list( c("A___No","A___Yes"),
                        c("A___No","C___No"),
                        c("A___Yes","C___Yes"),
                        c("C___No","C___Yes"),
                        c("C___No","B___Yes"))

ix <- which( my_melted_df$Protein %in% c("PDCD1","PD-L1","PD-L2") )

#pdf(file="u01_olink_colitis_figures/boxplots_clustered_significant_markers_colitis_time_PD1.pdf ",width = 12,height = 14)
ggviolin(data=my_melted_df[ix,], x="Status", y="value", color="Colitis",fill="Colitis",alpha=0.5) + theme_bw() +
  facet_wrap(~Protein+Study,scales = "free") +
  geom_boxplot(width=0.2) + geom_quasirandom(varwidth = TRUE,alpha=0.3,size=0.5) +
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test", paired=FALSE, label="p.signif", hide.ns = FALSE) +
  theme(axis.text.x = element_text(angle = 60,hjust = 1)) + theme(legend.position="bottom") + theme(#axis.title.x=element_blank(),
    #axis.text.x=element_blank(),
    #axis.ticks.x=element_blank(),
    strip.text=element_text(size=rel(1.5),face="bold"),
    legend.text = element_text(size=rel(1.5),face="bold"),
    legend.title = element_text(size=rel(1.5),face="bold"),
    axis.text.y = element_text(size=rel(1.5),face="bold"),
    axis.title.y = element_text(size=rel(1.5),face="bold")) +
  scale_fill_manual(values = c("firebrick","steelblue")) +
  scale_color_manual(values = c("firebrick","steelblue")) +
  ggtitle("") + labs(x = "Colitis", y = "RAW NPX Values") 
#dev.off()

#############################################################################################################################################
### Significant proteins, figures withtout correction
#############################################################################################################################################

### selection
ix <- which(colnames(colitis_metadata) %in% c("QC","Timepoint","Study","Colitis"))
annotation_data <- colitis_metadata[,ix]
rownames(annotation_data) <- colitis_metadata$sample_id
ix <- which( ! annotation_data$Study %in% "control" )
### iy <- which(protein_metadata$Experiment_Short %in% "I")
colitis_data_matrix_I <- as.matrix(colitis_data_matrix[,ix])
annotation_data <- annotation_data[ix,]

final <- final_result$ID[final_result$adj.P.Val<0.05 & final_result$logFC > 0.5]

my_position <- list()
for( i in 1:length(final)){ my_position[[i]] <- grep( final[i], rownames(colitis_data_matrix_I))   }
my_position <- unique(sort(unlist(my_position)))

### prepare for violin plots
my_melted_df <- melt(colitis_data_matrix_I[my_position,])
colnames(my_melted_df) <- c("Protein","sample_id","value")
my_melted_df <- merge(my_melted_df,colitis_metadata,by="sample_id")
#ix <- which(my_melted_df$Protein %in% "TNF...87")
## ggplot(data=my_melted_df) + aes(x=Colitis,y=value) + geom_boxplot() + facet_wrap(~Protein,scales = "free") + theme_bw()

my_melted_df$Protein <- as.character(my_melted_df$Protein)
my_melted_df$Protein <- tidyr::separate(data.frame(my_melted_df$Protein), 1, sep="\\.\\.\\.", c("a","b","c"))$a
my_melted_df$Protein <- as.factor(my_melted_df$Protein)
table(my_melted_df$Protein)

# grep("PD", my_melted_df$Protein,value=TRUE)
# data_for_figure ### Results from various mixed linear models.
my_melted_df$Status <- paste(my_melted_df$Timepoint,my_melted_df$Colitis,sep="___")
my_melted_df$Status <- factor(my_melted_df$Status, levels=c("A___No","A___Yes","B___Yes","C___No","C___Yes","D___Yes"))

my_comparisons <- list( c("A___No","A___Yes"),
                        c("A___No","C___No"),
                        c("A___Yes","C___Yes"),
                        c("C___No","C___Yes"),
                        c("C___No","B___Yes"))

ix <- which(my_melted_df$Study %in% names(table(my_melted_df$Study))[1])
#pdf(file="u01_olink_colitis_figures/boxplots_clustered_significant_markers_colitis_time_without_any_mod_study1.pdf ",width = 16,height = 7)
ggviolin(data=my_melted_df[ix,], x="Status", y="value", color="Colitis",fill="Colitis",alpha=0.5) + theme_bw() +
  facet_wrap(~Protein,scales = "free",ncol=7) +
  geom_boxplot(width=0.2) + geom_quasirandom(varwidth = TRUE,alpha=0.3,size=0.5) +
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test", paired=FALSE, label="p.signif", hide.ns = FALSE) +
  theme(axis.text.x = element_text(angle = 60,hjust = 1)) + theme(legend.position="bottom") + theme(#axis.title.x=element_blank(),
    #axis.text.x=element_blank(),
    #axis.ticks.x=element_blank(),
    strip.text=element_text(size=rel(1.5),face="bold"),
    legend.text = element_text(size=rel(1.5),face="bold"),
    legend.title = element_text(size=rel(1.5),face="bold"),
    axis.text.y = element_text(size=rel(1.5),face="bold"),
    axis.title.y = element_text(size=rel(1.5),face="bold")) +
  scale_fill_manual(values = c("firebrick","steelblue")) +
  scale_color_manual(values = c("firebrick","steelblue")) +
  ggtitle(paste("Study:",names(table(my_melted_df$Study))[1])) + labs(x = "Colitis", y = "RAW NPX Values") 
#dev.off()

ix <- which(my_melted_df$Study %in% names(table(my_melted_df$Study))[2])
#pdf(file="u01_olink_colitis_figures/boxplots_clustered_significant_markers_colitis_time_without_any_mod_study2.pdf ",width = 16,height = 7)
ggviolin(data=my_melted_df[ix,], x="Status", y="value", color="Colitis",fill="Colitis",alpha=0.5) + theme_bw() +
  facet_wrap(~Protein,scales = "free",ncol=7) +
  geom_boxplot(width=0.2) + geom_quasirandom(varwidth = TRUE,alpha=0.3,size=0.5) +
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test", paired=FALSE, label="p.signif", hide.ns = FALSE) +
  theme(axis.text.x = element_text(angle = 60,hjust = 1)) + theme(legend.position="bottom") + theme(#axis.title.x=element_blank(),
    #axis.text.x=element_blank(),
    #axis.ticks.x=element_blank(),
    strip.text=element_text(size=rel(1.5),face="bold"),
    legend.text = element_text(size=rel(1.5),face="bold"),
    legend.title = element_text(size=rel(1.5),face="bold"),
    axis.text.y = element_text(size=rel(1.5),face="bold"),
    axis.title.y = element_text(size=rel(1.5),face="bold")) +
  scale_fill_manual(values = c("firebrick","steelblue")) +
  scale_color_manual(values = c("firebrick","steelblue")) +
  ggtitle(paste("Study:",names(table(my_melted_df$Study))[2])) + labs(x = "Colitis", y = "RAW NPX Values") 
#dev.off()

ix <- which(my_melted_df$Study %in% names(table(my_melted_df$Study))[3])
#pdf(file="u01_olink_colitis_figures/boxplots_clustered_significant_markers_colitis_time_without_any_mod_study3.pdf ",width = 16,height = 7)
ggviolin(data=my_melted_df[ix,], x="Status", y="value", color="Colitis",fill="Colitis",alpha=0.5) + theme_bw() +
  facet_wrap(~Protein,scales = "free",ncol=7) +
  geom_boxplot(width=0.2) + geom_quasirandom(varwidth = TRUE,alpha=0.3,size=0.5) +
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test", paired=FALSE, label="p.signif", hide.ns = FALSE) +
  theme(axis.text.x = element_text(angle = 60,hjust = 1)) + theme(legend.position="bottom") + theme(#axis.title.x=element_blank(),
    #axis.text.x=element_blank(),
    #axis.ticks.x=element_blank(),
    strip.text=element_text(size=rel(1.5),face="bold"),
    legend.text = element_text(size=rel(1.5),face="bold"),
    legend.title = element_text(size=rel(1.5),face="bold"),
    axis.text.y = element_text(size=rel(1.5),face="bold"),
    axis.title.y = element_text(size=rel(1.5),face="bold")) +
  scale_fill_manual(values = c("firebrick","steelblue")) +
  scale_color_manual(values = c("firebrick","steelblue")) +
  ggtitle(paste("Study:",names(table(my_melted_df$Study))[3])) + labs(x = "Colitis", y = "RAW NPX Values") 
#dev.off()

#############################################################################################################################################
#############################################################################################################################################
### Repeat models per study
### Here, we use a different approach. We repeat the analysis above to verify results, using a similar yet different modeling strategy
### same matrix, separated by study.
#############################################################################################################################################
#############################################################################################################################################

table(colitis_metadata$Colitis,colitis_metadata$Timepoint)
dim(data_for_figure)
table(data_for_figure$Contrast)
table(data_for_figure$Data[data_for_figure$Contrast %in% "Colitis: Yes - No (C)"][data_for_figure$adj.P.Val[data_for_figure$Contrast %in% "Colitis: Yes - No (C)"]<0.05])
table(data_for_figure$Contrast %in% "Colitis: Yes - No (C)" & data_for_figure$Data %in% "09-155" & data_for_figure$adj.P.Val<0.05)
table(data_for_figure$Contrast %in% "Colitis: Yes - No (C)" & data_for_figure$Data %in% "09-155" & data_for_figure$P.Value<0.05) ### 7
table(data_for_figure$Contrast %in% "Colitis: Yes - No (A)" & data_for_figure$Data %in% "13-075, 17-162, 13-169" & data_for_figure$adj.P.Val<0.05)
table(data_for_figure$Contrast %in% "Colitis: Yes - No (C)" & data_for_figure$Data %in% "CTLA-4" & data_for_figure$adj.P.Val<0.05)

dim(colitis_data_matrix)
dim(colitis_metadata)
dim(protein_metadata)

### split

ix <- which(colnames(colitis_metadata) %in% c("QC","Timepoint","Study","Colitis"))
annotation_data <- colitis_metadata[,ix]
rownames(annotation_data) <- colitis_metadata$sample_id
ix <- which( ! annotation_data$Study %in% "control" )
### iy <- which(protein_metadata$Experiment_Short %in% "I")
colitis_data_matrix_I <- as.matrix(colitis_data_matrix[,ix])
annotation_data <- annotation_data[ix,]

annotation_data$Study_simple <- as.character(annotation_data$Study)
annotation_data$Study_simple[annotation_data$Study_simple %in% "09-155"] <- "A"
annotation_data$Study_simple[annotation_data$Study_simple %in% "13-075, 17-162, 13-169"] <- "B"
annotation_data$Study_simple[annotation_data$Study_simple %in% "CTLA-4"] <- "C"

annotation_data$Condition <- paste( annotation_data$Study_simple,annotation_data$Colitis,annotation_data$Timepoint,sep="___")

table(annotation_data$Colitis[annotation_data$Study_simple %in% "A"],annotation_data$Timepoint[annotation_data$Study_simple %in% "A"])
table(annotation_data$Colitis[annotation_data$Study_simple %in% "B"],annotation_data$Timepoint[annotation_data$Study_simple %in% "B"])
table(annotation_data$Colitis[annotation_data$Study_simple %in% "C"],annotation_data$Timepoint[annotation_data$Study_simple %in% "C"])

### time point specific model
design <- model.matrix( ~ 0 + Condition , data = annotation_data )
colnames(design) <- gsub("Condition","",colnames(design))
colnames(design) <- make.names(colnames(design))
### comparisons
contr.matrix <- makeContrasts( 
  "Study_A_Baseline_Yes_vs_No" = A___Yes___A - A___No___A,
  "Study_A_Week6_Yes_vs_No" = A___Yes___C - A___No___C,
  
  "Study_B_Baseline_Yes_vs_No" = B___Yes___A - B___No___A,
  "Study_B_Week6_Yes_vs_No" = B___Yes___C - B___No___C,
  
  "Study_C_Baseline_Yes_vs_No" = C___Yes___A - C___No___A,
  "Study_C_Week6_Yes_vs_No" = C___Yes___C - C___No___C,

  levels = colnames(design) )
### STATs
#block2 <- as.numeric(as.factor( annotation_data$Patient ))
#dupcor2 <- duplicateCorrelation(colitis_data_matrix_I, design, block=block2)
#vfit <- lmFit(colitis_data_matrix_I, design, block=block2, correlation=dupcor2$consensus)
vfit <- lmFit(colitis_data_matrix_I, design)
vfit <- contrasts.fit(vfit, contrasts=contr.matrix)
efit <- eBayes(vfit, robust = TRUE) ; summary(decideTests(efit))

lmfreq_results <- list()
for ( i in 1:ncol(summary(decideTests(efit))) ) { lmfreq_results[[i]] <- topTable(efit, coef=i,n=Inf,adjust.method="BH") 
lmfreq_results[[i]]$Comparison <- colnames(summary(decideTests(efit)))[i]
lmfreq_results[[i]]$Protein <- rownames( topTable(efit, coef=i,n=Inf,adjust.method="BH") ) }
lmfreq_results <- do.call(rbind,lmfreq_results)
rownames(lmfreq_results) <- NULL
###
###### from single study comparisons
write.csv(file="statistics_single_comparisons.csv",lmfreq_results)
###### significant results from study corrected models
write.csv(file="statistics_adjusted_models_comparisons.csv",final_result)

###
#############################################################################################################################################
#############################################################################################################################################
### Based on sugestions from Sacha, We'll merge both panels, treating them as independed measurements
### IO + I panel merge
#############################################################################################################################################
#############################################################################################################################################
###

ix <- which(colnames(colitis_metadata) %in% c("QC","Timepoint","Study","Colitis"))
annotation_data <- colitis_metadata[,ix]
rownames(annotation_data) <- colitis_metadata$sample_id
ix <- which( ! annotation_data$Study %in% "control" )
### iy <- which(protein_metadata$Experiment_Short %in% "I")
colitis_data_matrix_I <- as.matrix(colitis_data_matrix[,ix])
annotation_data <- annotation_data[ix,]

colitis_data_matrix_I <- colitis_data_matrix_I[order(rownames(colitis_data_matrix_I)),]
### identify duplicated proteins
my_proteins <- rownames(colitis_data_matrix_I)
my_proteins <- tidyr::separate(data.frame(my_proteins), 1, sep="\\.\\.\\.", c("a","b","c"))$a
my_proteins <- my_proteins[which(duplicated(my_proteins))]

corr.res.imol<-list()
for ( i in 1:length(my_proteins)){
#print(grep(my_proteins[i],rownames(colitis_data_matrix_I)))
#print(cor.test(colitis_data_matrix_I[a,],colitis_data_matrix_I[b,],method="spearman")$estimate)
a<-grep(my_proteins[i],rownames(colitis_data_matrix_I))[1]
b<-grep(my_proteins[i],rownames(colitis_data_matrix_I))[2]
corr.res.imol[[i]] <- data.frame( rho = cor.test(colitis_data_matrix_I[a,],colitis_data_matrix_I[b,],method="spearman")$estimate, 
                                  molecule = grep(my_proteins[i],rownames(colitis_data_matrix_I),value=TRUE)[1],
                                  posa=a,
                                  posb=b,
                                  average_a=median(colitis_data_matrix_I[a,]),
                                  average_b=median(colitis_data_matrix_I[b,])) }

corr.res.imol <- do.call(rbind,corr.res.imol)
corr.res.imol$symbol <- my_proteins
corr.res.imol <- corr.res.imol[order(corr.res.imol$rho,decreasing=TRUE),]
### to average these:
my_proteins <- corr.res.imol$symbol[corr.res.imol$rho>0.88]

average_between_dup_markers <- list()
for ( i in 1:length(my_proteins)){
  a<-grep(my_proteins[i],rownames(colitis_data_matrix_I))[1]
  b<-grep(my_proteins[i],rownames(colitis_data_matrix_I))[2]
  average_between_dup_markers[[i]] <- colMedians(colitis_data_matrix_I[a:b,])  }

average_between_dup_markers <- do.call(cbind,average_between_dup_markers)
colnames(average_between_dup_markers) <- my_proteins

average_between_dup_markers <- as.data.frame(t(average_between_dup_markers))
colnames(average_between_dup_markers) <- colnames(colitis_data_matrix_I)

dim(average_between_dup_markers)
#corr.res.imol
#rownames(average_between_dup_markers)
summary(colitis_data_matrix_I[116,])
dim(corr.res.imol)
corr.res.imol$average_b[43] <- 1.1360
corr.res.imol$average_b[42] <- 0.8159
corr.res.imol$average_b[41] <- 0.11921
corr.res.imol$average_b[39] <- 1.1984

my_proteins <- c(corr.res.imol$posa[corr.res.imol$rho>0.88],corr.res.imol$posb[corr.res.imol$rho>0.88])
colitis_data_matrix_I <- colitis_data_matrix_I[-my_proteins,]
test <- rbind(colitis_data_matrix_I,average_between_dup_markers)
colitis_data_matrix_I <- test ; rm(test,average_between_dup_markers,my_proteins)

### rename
colitis_data_matrix_iio <- colitis_data_matrix_I
corr.res.imol$difference <- abs(corr.res.imol$average_a-corr.res.imol$average_b)

rownames(corr.res.imol) <- NULL
write.csv(file="correlation_duplicated_variables.csv",corr.res.imol)

### CORRELATION BETWEEN I AND IO PANELS
#############################################################################################################################################
pdf(file="u01_olink_colitis_figures/correlation_between_duplicated_panel_proteins_summary.pdf",width = 8,height = 6)
ggplot(data=corr.res.imol,aes(x=average_a,y=average_b) ) + theme_bw() +
geom_point(shape=16, colour = "black", aes(size = 1)) +
  geom_point(shape=16, colour = "white", aes(size = 0.8*max(rho))) +
  geom_point(shape=16, aes(size = 0.81*rho,color=rho)) + 
  scale_color_viridis()+ theme_classic2()+ 
  labs(x ='NPX IO', y='NPX I', title='Proteins common to IO/I panels',colour="Rho",size="Rho") + 
  geom_label_repel(aes(label=symbol), force = 10, max.iter = 5000,size=3) +
  theme(text = element_text(size=rel(4.8)),
      legend.text=element_text(size=rel(4.8)),
      axis.title.x = element_text(size=rel(4.8)),
      plot.title=element_text(size=rel(4.8)) )
dev.off()
#############################################################################################################################################

#############################################################################################################################################
### REVIEW clustering strategy and markers 
#############################################################################################################################################

### Review proteins in common between corrected and not corrected data
### 

list(a=my_top_proteins_global,b=my_top_proteins_local)
Reduce(intersect,list(a=my_top_proteins_global,b=my_top_proteins_local))
my_top_proteins <- unique(sort(unlist(list(a=my_top_proteins_global,b=my_top_proteins_local))))

### Subset of annotation data for figure

ix <- which(colnames(colitis_metadata) %in% c("QC","Timepoint","Study","Colitis"))
annotation_data <- colitis_metadata[,ix]
rownames(annotation_data) <- colitis_metadata$sample_id
ix <- which( ! annotation_data$Study %in% "control" )
annotation_data <- annotation_data[ix,]

### Select first the top 7 proteins showed before
### These are comming from study adjusted data.

final_result$Proteins <- tidyr::separate(data.frame(final_result$Proteins), 1, sep="\\.\\.\\.", c("a","b","c"))$a
my_top_proteins <- final_result$Proteins[final_result$adj.P.Val<0.05 & final_result$logFC>0.5]
my_top_proteins <- unique(sort(my_top_proteins))

### NOTE:
### use annotation data from original clustering ###

library(sva)
colitis_data_matrix_II <- as.matrix(colitis_data_matrix_iio)
Batch = as.numeric(as.factor(as.character(annotation_data$Study)))
modcombat = model.matrix(~1, data=annotation_data)
colitis_data_matrix_I_study_corrected = ComBat(dat=colitis_data_matrix_II, batch=Batch, mod=modcombat, par.prior=TRUE, prior.plots=FALSE)

ix <- which(rownames(colitis_data_matrix_I_study_corrected) %in% final_result$Proteins[final_result$adj.P.Val<0.05 & final_result$logFC > 0.5])

out_matrix <- pheatmap( colitis_data_matrix_I_study_corrected[ix,],
                        main="Study Corrected NPX values",
                        show_rownames = TRUE,
                        show_colnames = FALSE,
                        scale="row",
                        annotation_col = annotation_data)

annotation_data$HCluster <- as.character(cutree(out_matrix$tree_col,k=4))
annotation_data$HCluster <- paste("HCluster",annotation_data$HCluster,sep="_")

iy <- c( grep( "HCluster_2",annotation_data$HCluster),grep( "HCluster_1",annotation_data$HCluster),grep( "HCluster_3",annotation_data$HCluster),grep( "HCluster_4",annotation_data$HCluster) )
annotation_data <- annotation_data[iy,]

iy <- match(rownames(annotation_data), colnames(colitis_data_matrix_I_study_corrected ) )

identical(rownames(annotation_data),colnames(colitis_data_matrix_I_study_corrected )[iy])
identical(colnames(colitis_data_matrix_iio),colnames(colitis_data_matrix_I_study_corrected ))

annotation_row <- protein_metadata
annotation_row <- annotation_row [ which(annotation_row$Marker %in% rownames(colitis_data_matrix_I_study_corrected[,]) ) , ]
annotation_row <- annotation_row[,c("QC1","Experiment")]  
annotation_row$QC1 <- round(as.numeric(as.character(annotation_row$QC1)),2)
annotation_row$Experiment[annotation_row$Experiment %in% "Olink Target 96 Inflammation(v.3022)"] <- "I"
annotation_row$Experiment[annotation_row$Experiment %in% "Olink Target 96 Immuno-Oncology(v.3112)"] <- "IO"

annotation_colors <- list ( HCluster = c( HCluster_1="#781C6DFF",HCluster_2="#000004FF",HCluster_3="#ED6925FF",HCluster_4="#FCFFA4FF"),
                            Colitis = c( No="steelblue",Yes="firebrick"),
                            Timepoint = c( A="grey20",B="grey40",C="grey60",D="grey80"),
                            Study = c( `CTLA-4`="green",`13-075, 17-162, 13-169`="orange",`09-155`="purple"),
                            QC = c( Pass="black",Warning="red"),
                            Experiment = c(I="grey10",IO="grey50") )

### To compare with raw pnx values we want to know the position of these 7 markers in the original matrix
### To do this, we use a regex expression using a trick with paste and grep functions

### prepare exact match regex expression
my_top_proteins <- paste( "^",my_top_proteins,"$" , sep="")
### fish proteins position in matrix with the regex expression
my_position <- list()
for( i in 1:length(my_top_proteins)){ my_position[[i]] <- grep( my_top_proteins[i], rownames(colitis_data_matrix_iio))   }
my_position <- unique(sort(unlist(my_position)))

### sort matrices based on clusters
colitis_data_matrix_I_study_corrected <- colitis_data_matrix_I_study_corrected[,iy]
colitis_data_matrix_iio <- colitis_data_matrix_iio[,iy]

png(file="u01_olink_colitis_figures/pheatmap_clustered_significant_markers_colitis_all_models_significant_study_corrected.png", width = 10, height = 3.5,units = "in", res = 300)
pheatmap( colitis_data_matrix_I_study_corrected[my_position,],
          color = colorRampPalette(c("steelblue", "white","firebrick"))(255),
          main="Study Corrected NPX values",
          cluster_cols = FALSE,
          show_rownames = TRUE,
          show_colnames = FALSE,
          scale="row",
          gaps_col = c(as.numeric(cumsum(table(annotation_data$HCluster)[c("HCluster_2","HCluster_1","HCluster_3","HCluster_4")]))) ,
          border_color=FALSE,
          fontsize=5,
          angle_col = 90,
          annotation_col = annotation_data,
          annotation_row = annotation_row,
          annotation_colors = annotation_colors )
dev.off()
dev.off()

png(file="u01_olink_colitis_figures/pheatmap_clustered_significant_markers_colitis_all_models_significant_raw_npx.png", width = 10, height = 3.5,units = "in", res = 300)
pheatmap( colitis_data_matrix_iio[my_position,],
          color = colorRampPalette(c("steelblue", "white","firebrick"))(255),
          main="Original NPX values",
          cluster_cols = FALSE,
          show_rownames = TRUE,
          show_colnames = FALSE,
          scale="row",
          gaps_col = c(as.numeric(cumsum(table(annotation_data$HCluster)[c("HCluster_2","HCluster_1","HCluster_3","HCluster_4")]))) ,
          border_color=FALSE,
          fontsize=5,
          angle_col = 90,
          annotation_col = annotation_data,
          annotation_row = annotation_row,
          annotation_colors = annotation_colors )
dev.off()
dev.off()

### prepare for violin plots
my_melted_df <- melt(colitis_data_matrix_I_study_corrected[my_position,])
colnames(my_melted_df) <- c("Protein","sample_id","value")
my_melted_df <- merge(my_melted_df,colitis_metadata,by="sample_id")
#ix <- which(my_melted_df$Protein %in% "TNF...87")
## ggplot(data=my_melted_df) + aes(x=Colitis,y=value) + geom_boxplot() + facet_wrap(~Protein,scales = "free") + theme_bw()
my_melted_df$Protein <- as.character(my_melted_df$Protein)
my_melted_df$Protein <- tidyr::separate(data.frame(my_melted_df$Protein), 1, sep="\\.\\.\\.", c("a","b","c"))$a
my_melted_df$Protein <- as.factor(my_melted_df$Protein)
table(my_melted_df$Protein)
### my_melted_df$Protein <- factor(my_melted_df$Protein, levels = c("CXCL10","CXCL9","IFN-gamma","TNF","MMP-10","IL-17C","CXCL1"))  
library(ggbeeswarm)
my_melted_df$Status <- paste(my_melted_df$Timepoint,my_melted_df$Colitis,sep="___")
my_melted_df$Status <- gsub("A_","Baseline_",my_melted_df$Status)
my_melted_df$Status <- gsub("D_","Post_",my_melted_df$Status)
my_melted_df$Status <- gsub("B_","Colitis_",my_melted_df$Status)
my_melted_df$Status <- gsub("C_","Week6_",my_melted_df$Status)

my_melted_df$Status <- factor(my_melted_df$Status, levels=c("Baseline___No","Baseline___Yes","Week6___No","Week6___Yes","Colitis___Yes","Post___Yes"))

my_comparisons <- list( c("Baseline___No","Baseline___Yes"),
                        c("Baseline___No","Week6___No"),
                        c("Baseline___Yes","Week6___Yes"),
                        c("Week6___No","Week6___Yes"),
                        c("Week6___No","Colitis___Yes"))


pdf(file="u01_olink_colitis_figures/boxplots_clustered_significant_markers_colitis_time_combination_models_significant_study_corrected.pdf ",width = 18,height = 7)
ggviolin(data=my_melted_df, x="Status", y="value", color="Colitis",fill="Colitis",alpha=0.5) + theme_bw() +
  facet_wrap(~Protein,scales = "free",ncol=7) +
  geom_boxplot(width=0.2) + geom_quasirandom(varwidth = TRUE,alpha=0.3,size=0.5) +
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test", paired=FALSE, label="p.signif", hide.ns = FALSE) +
  theme(axis.text.x = element_text(angle = 45,hjust = 1)) + theme(legend.position="bottom") + theme(#axis.title.x=element_blank(),
    #axis.text.x=element_blank(),
    #axis.ticks.x=element_blank(),
    strip.text=element_text(size=rel(1.5),face="bold"),
    legend.text = element_text(size=rel(1.5),face="bold"),
    legend.title = element_text(size=rel(1.5),face="bold"),
    axis.text.y = element_text(size=rel(1.5),face="bold"),
    axis.title.y = element_text(size=rel(1.5),face="bold")) +
  scale_fill_manual(values = c("firebrick","steelblue")) +
  scale_color_manual(values = c("firebrick","steelblue")) +
  ggtitle("") + labs(x = "", y = "Study Adj. NPX values") +
  theme(text = element_text(size=rel(4.8)),
        legend.text = element_text(size=rel(4.8)),
        legend.title = element_text(size=rel(4.8)),
        plot.subtitle = element_text(size=rel(4.8)),
        axis.title.x = element_text(size=rel(4.8)),
        axis.text.y = element_text(size=rel(4.8)),
        strip.text = element_text(size=rel(4.8)),
        axis.title.y = element_text(size=rel(4.8)),
        plot.title=element_text(size=rel(4.8)) ) +
  theme(plot.margin=unit(c(t=5,r=5,b=5,l=100),"pt"))
dev.off()

### This needs to be TRUE
identical(rownames(colitis_data_matrix_iio), rownames(colitis_data_matrix_I_study_corrected) )

### prepare for violin plots
my_melted_df <- melt(as.matrix(colitis_data_matrix_iio[my_position,]))
colnames(my_melted_df) <- c("Protein","sample_id","value")
my_melted_df <- merge(my_melted_df,colitis_metadata,by="sample_id")
#ix <- which(my_melted_df$Protein %in% "TNF...87")
## ggplot(data=my_melted_df) + aes(x=Colitis,y=value) + geom_boxplot() + facet_wrap(~Protein,scales = "free") + theme_bw()
my_melted_df$Protein <- as.character(my_melted_df$Protein)
my_melted_df$Protein <- tidyr::separate(data.frame(my_melted_df$Protein), 1, sep="\\.\\.\\.", c("a","b","c"))$a
my_melted_df$Protein <- as.factor(my_melted_df$Protein)
table(my_melted_df$Protein)
### my_melted_df$Protein <- factor(my_melted_df$Protein, levels = c("CXCL10","CXCL9","IFN-gamma","TNF","MMP-10","IL-17C","CXCL1"))  
library(ggbeeswarm)
my_melted_df$Status <- paste(my_melted_df$Timepoint,my_melted_df$Colitis,sep="___")
my_melted_df$Status <- gsub("A_","Baseline_",my_melted_df$Status)
my_melted_df$Status <- gsub("D_","Post_",my_melted_df$Status)
my_melted_df$Status <- gsub("B_","Colitis_",my_melted_df$Status)
my_melted_df$Status <- gsub("C_","Week6_",my_melted_df$Status)

my_melted_df$Status <- factor(my_melted_df$Status, levels=c("Baseline___No","Baseline___Yes","Week6___No","Week6___Yes","Colitis___Yes","Post___Yes"))

my_comparisons <- list( c("Baseline___No","Baseline___Yes"),
                        c("Baseline___No","Week6___No"),
                        c("Baseline___Yes","Week6___Yes"),
                        c("Week6___No","Week6___Yes"),
                        c("Week6___No","Colitis___Yes"))

pdf(file="u01_olink_colitis_figures/boxplots_clustered_significant_markers_colitis_time_all_models_significant_npx_raw.pdf ",width = 18,height = 7)
ggviolin(data=my_melted_df, x="Status", y="value", color="Colitis",fill="Colitis",alpha=0.5) + theme_bw() +
  facet_wrap(~Protein,scales = "free",ncol=7) +
  geom_boxplot(width=0.2) + geom_quasirandom(varwidth = TRUE,alpha=0.3,size=0.5) +
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test", paired=FALSE, label="p.signif", hide.ns = FALSE) +
  theme(axis.text.x = element_text(angle = 45,hjust = 1)) + theme(legend.position="bottom") + theme(#axis.title.x=element_blank(),
    #axis.text.x=element_blank(),
    #axis.ticks.x=element_blank(),
    strip.text=element_text(size=rel(1.5),face="bold"),
    legend.text = element_text(size=rel(1.5),face="bold"),
    legend.title = element_text(size=rel(1.5),face="bold"),
    axis.text.y = element_text(size=rel(1.5),face="bold"),
    axis.title.y = element_text(size=rel(1.5),face="bold")) +
  scale_fill_manual(values = c("firebrick","steelblue")) +
  scale_color_manual(values = c("firebrick","steelblue")) +
  ggtitle("") + labs(x = "", y = "Raw NPX values") +
  theme(text = element_text(size=rel(4.8)),
        legend.text = element_text(size=rel(4.8)),
        legend.title = element_text(size=rel(4.8)),
        plot.subtitle = element_text(size=rel(4.8)),
        axis.title.x = element_text(size=rel(4.8)),
        axis.text.y = element_text(size=rel(4.8)),
        strip.text = element_text(size=rel(4.8)),
        axis.title.y = element_text(size=rel(4.8)),
        plot.title=element_text(size=rel(4.8)) ) +
  theme(plot.margin=unit(c(t=5,r=5,b=5,l=100),"pt"))
dev.off()

#############################################################################################################################################
### Volcano
#############################################################################################################################################

ix <- which(colnames(colitis_metadata) %in% c("QC","Timepoint","Study","Colitis"))
annotation_data <- colitis_metadata[,ix]
rownames(annotation_data) <- colitis_metadata$sample_id
ix <- which( ! annotation_data$Study %in% "control" )

colitis_data_matrix_I <- as.matrix(colitis_data_matrix[,ix])
annotation_data <- annotation_data[ix,]

library(sva)
colitis_data_matrix_I <- as.matrix(colitis_data_matrix_I)
Batch = as.numeric(as.factor(as.character(annotation_data$Study)))
modcombat = model.matrix(~1, data=annotation_data)
colitis_data_matrix_I_study_corrected = ComBat(dat=colitis_data_matrix_I, batch=Batch, mod=modcombat, par.prior=TRUE, prior.plots=FALSE)
colitis_data_matrix_I_study_corrected <- colitis_data_matrix_I_study_corrected

# comparisons
design <- model.matrix( ~ 0 + Colitis, data = annotation_data ) ; colnames(design) <- c("No","Yes")
contr.matrix <- makeContrasts(  "Colitis" = Yes - No, levels = colnames(design) )
vfit <- lmFit(colitis_data_matrix_I_study_corrected, design) ; vfit <- contrasts.fit(vfit, contrasts=contr.matrix)
efit <- eBayes(vfit, robust = TRUE) ; summary(decideTests(efit))
### final result
final_result <- topTable(efit, coef="Colitis", n=Inf, adjust.method="BH")
###
head(final_result)
###
final_result$ID <- rownames(final_result) # mod
final_result$Comparisons <- " "
final_result$Proteins <- tidyr::separate(data.frame(final_result$ID), 1, sep="\\.\\.\\.", c("a","b","c"))$a
final_result <- final_result[ order(final_result$Proteins), ]
final_result$nLog10FDR <- -log10(final_result$adj.P.Val)

ix<-which(!duplicated(final_result$Proteins))
final_result2 <- final_result[ix,]

#pdf(file="u01_olink_colitis_figures/volcano_clustered_significant_markers_colitis.pdf ",width = 6,height = 8)
ggplot(final_result2)+aes(logFC,nLog10FDR)+
  geom_point() + 
  geom_point(data=final_result2[which(final_result2$adj.P.Val>0.05),],color="grey50") + 
  theme_bw() + 
  geom_vline(xintercept = 0, linetype = 2, color = 'black') +
  geom_vline(xintercept = 0.5, linetype = 2, color = 'red') +
  geom_hline(yintercept = 1.3, linetype = 2, color = 'red') +
  geom_label_repel(data=final_result2[which(final_result2$adj.P.Val<0.1 & final_result2$logFC>0.3),],aes(label=Proteins))
#dev.off()

annotation_row <- protein_metadata
annotation_row <- annotation_row [ which(annotation_row$Marker %in% rownames(colitis_data_matrix_I_study_corrected[ix,]) ) , ]
annotation_row <- annotation_row[,c("QC1","Experiment")]  
annotation_row$QC1 <- round(as.numeric(as.character(annotation_row$QC1)),2)
annotation_row$Experiment[annotation_row$Experiment %in% "Olink Target 96 Inflammation(v.3022)"] <- "I"
annotation_row$Experiment[annotation_row$Experiment %in% "Olink Target 96 Immuno-Oncology(v.3112)"] <- "IO"

annotation_colors <- list ( HCluster = c( HCluster_1="#781C6DFF",HCluster_2="#000004FF",HCluster_3="#ED6925FF",HCluster_4="#FCFFA4FF"),
                            Colitis = c( No="steelblue",Yes="firebrick"),
                            Timepoint = c( A="grey20",B="grey40",C="grey60",D="grey80"),
                            Study = c( Study1="green",Study2="orange",Study3="purple"),
                            QC = c( Pass="black",Warning="red"),
                            Experiment = c(I="grey10",IO="grey50") )

rownames(colitis_data_matrix_I_study_corrected) <- tidyr::separate(data.frame(rownames(colitis_data_matrix_I_study_corrected)), 1, sep="\\.\\.\\.", c("a","b","c"))$a

ix <- which(rownames(colitis_data_matrix_I_study_corrected) %in% final_result2$Proteins[which(final_result2$adj.P.Val<0.1 & final_result2$logFC>0.5)])

colitis_data_matrix_I_study_corrected <- colitis_data_matrix_I_study_corrected[ , match(rownames(annotation_data_col_save), colnames(colitis_data_matrix_I_study_corrected)) ]
colitis_data_matrix_I_study_corrected <- colitis_data_matrix_I_study_corrected[ix,]
colitis_data_matrix_I_study_corrected <- colitis_data_matrix_I_study_corrected[unique(rownames(colitis_data_matrix_I_study_corrected)),]
colitis_data_matrix_I_study_corrected <- colitis_data_matrix_I_study_corrected[-which(duplicated(rownames(colitis_data_matrix_I_study_corrected))),]

#png(file="u01_olink_colitis_figures/pheatmap_clustered_significant_markers_colitis.png", width = 10, height = 4,units = "in", res = 300)
pheatmap( colitis_data_matrix_I_study_corrected,
          color = colorRampPalette(c("steelblue", "white","firebrick"))(255),
          main="Study Corrected NPX values",
          cluster_cols = FALSE,
          show_rownames = TRUE,
          show_colnames = FALSE,
          scale="row",
          gaps_col = c(as.numeric(cumsum(table(annotation_data_col_save$HCluster)[c("HCluster_2","HCluster_1","HCluster_3","HCluster_4")]))) ,
          border_color=FALSE,
          fontsize=6,
          angle_col = 90,
          annotation_col = annotation_data_col_save,
          #annotation_row = annotation_row,
          annotation_colors = annotation_colors )
#dev.off()
dev.off()

### prepare for violin plots
my_melted_df <- melt(colitis_data_matrix_I_study_corrected)
colnames(my_melted_df) <- c("Protein","sample_id","value")
my_melted_df <- merge(my_melted_df,colitis_metadata,by="sample_id")

my_melted_df$Protein <- as.character(my_melted_df$Protein)
my_melted_df$Protein <- as.factor(my_melted_df$Protein)
my_melted_df$Protein <- factor(my_melted_df$Protein, levels = c("CXCL10","CXCL9","IFN-gamma","TNF","MMP-10","IL-17C","IL-17A","CXCL1"))  

library(ggbeeswarm)
my_comparisons <- list( c("Yes","No") )

#pdf(file="u01_olink_colitis_figures/boxplots_clustered_significant_markers_colitis.pdf ",width = 8,height = 7)
ggviolin(data=my_melted_df, x="Colitis", y="value", color="Colitis",fill="Colitis",alpha=0.5) + theme_bw() +
  facet_wrap(~Protein,scales = "free") +
  geom_boxplot(width=0.2) + geom_quasirandom(varwidth = TRUE,alpha=0.3,size=0.5) +
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test", paired=FALSE, label="p.signif", hide.ns = TRUE) +
  theme(axis.text.x = element_text(angle = 60,hjust = 1)) + theme(legend.position="bottom") + theme(axis.title.x=element_blank(),
                                                                                                    axis.text.x=element_blank(),
                                                                                                    axis.ticks.x=element_blank(),
                                                                                                    strip.text=element_text(size=rel(1.5),face="bold"),
                                                                                                    legend.text = element_text(size=rel(1.5),face="bold"),
                                                                                                    legend.title = element_text(size=rel(1.5),face="bold"),
                                                                                                    axis.text.y = element_text(size=rel(1.5),face="bold"),
                                                                                                    axis.title.y = element_text(size=rel(1.5),face="bold")) +
  scale_fill_manual(values = c("firebrick","steelblue")) +
  scale_color_manual(values = c("firebrick","steelblue")) +
  ggtitle("") + labs(x = "Colitis", y = "NPX") 
#dev.off()

my_melted_df$Status <- paste(my_melted_df$Timepoint,my_melted_df$Colitis,sep="___")
my_melted_df$Status <- factor(my_melted_df$Status, levels=c("A___No","A___Yes","B___Yes","C___No","C___Yes","D___Yes"))

my_comparisons <- list( c("A___No","A___Yes"),
                        c("A___No","C___No"),
                        c("A___Yes","C___Yes"),
                        c("C___No","C___Yes"),
                        c("C___No","B___Yes"))

#pdf(file="u01_olink_colitis_figures/boxplots_clustered_significant_markers_colitis_time.pdf ",width = 18,height = 7)
ggviolin(data=my_melted_df, x="Status", y="value", color="Colitis",fill="Colitis",alpha=0.5) + theme_bw() +
  facet_wrap(~Protein,scales = "free",ncol=8) +
  geom_boxplot(width=0.2) + geom_quasirandom(varwidth = TRUE,alpha=0.3,size=0.5) +
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test", paired=FALSE, label="p.signif", hide.ns = FALSE) +
  theme(axis.text.x = element_text(angle = 60,hjust = 1)) + theme(legend.position="bottom") + theme(#axis.title.x=element_blank(),
    #axis.text.x=element_blank(),
    #axis.ticks.x=element_blank(),
    strip.text=element_text(size=rel(1.5),face="bold"),
    legend.text = element_text(size=rel(1.5),face="bold"),
    legend.title = element_text(size=rel(1.5),face="bold"),
    axis.text.y = element_text(size=rel(1.5),face="bold"),
    axis.title.y = element_text(size=rel(1.5),face="bold")) +
  scale_fill_manual(values = c("firebrick","steelblue")) +
  scale_color_manual(values = c("firebrick","steelblue")) +
  ggtitle("") + labs(x = "Colitis", y = "Adj.NPX") 
#dev.off()


#############################################################################################################################################
### full model to show
#############################################################################################################################################
###
ix <- which(colnames(colitis_metadata) %in% c("QC","Timepoint","Study","Colitis"))
annotation_data <- colitis_metadata[,ix]
rownames(annotation_data) <- colitis_metadata$sample_id

ix <- which( ! annotation_data$Study %in% "control" )
colitis_data_matrix_I <- as.matrix(colitis_data_matrix[,ix])
annotation_data <- annotation_data[ix,]

library(sva)
colitis_data_matrix_I <- as.matrix(colitis_data_matrix_I)
Batch = as.numeric(as.factor(as.character(annotation_data$Study)))
modcombat = model.matrix(~1, data=annotation_data)
colitis_data_matrix_I_study_corrected = ComBat(dat=colitis_data_matrix_I, batch=Batch, mod=modcombat, par.prior=TRUE, prior.plots=FALSE)

annotation_data$CT <- paste(annotation_data$Colitis,annotation_data$Timepoint,sep="___")

# comparisons
design <- model.matrix( ~ 0 + CT, data = annotation_data ) ; colnames(design) <- gsub("CT","",colnames(design))
contr.matrix <- makeContrasts(  "Colitis_A" = Yes___A - No___A, 
                                "Colitis_C" = Yes___C - No___C, 
                                "Not_Colitis_C_A" = No___C - No___A, 
                                "Colitis_B_A" = Yes___B - Yes___A,
                                "Colitis_C_A" = Yes___C - Yes___A,
                                "Colitis_D_A" = Yes___D - Yes___A,
                                "Colitis_C_B" = Yes___C - Yes___B,
                                "Colitis_D_B" = Yes___D - Yes___B,
                                "Colitis_D_C" = Yes___D - Yes___C,
                                levels = colnames(design) )

vfit <- lmFit(colitis_data_matrix_I_study_corrected, design) ; vfit <- contrasts.fit(vfit, contrasts=contr.matrix)
efit <- eBayes(vfit, robust = TRUE) ; summary(decideTests(efit))
### final result
final_result_corrected_model<-list()
for ( iii in 1:length(colnames(summary(decideTests(efit))))){
  final_result_corrected_model[[iii]] <- topTable(efit, coef=colnames(summary(decideTests(efit)))[iii], n=Inf, adjust.method="BH")
  final_result_corrected_model[[iii]]$Contrast <- colnames(summary(decideTests(efit)))[iii]
  final_result_corrected_model[[iii]]$Proteins <- rownames(final_result_corrected_model[[iii]]) }

final_result_corrected_model <- do.call(rbind,final_result_corrected_model)

final_result_corrected_model$Contrast_long <- final_result_corrected_model$Contrast

final_result_corrected_model$Contrast[final_result_corrected_model$Contrast %in% "Colitis_A"] <- "Baseline: (Colitis vs No Colitis)"
final_result_corrected_model$Contrast[final_result_corrected_model$Contrast %in% "Colitis_C"] <- "Week 6: (Colitis vs No Colitis)"
final_result_corrected_model$Contrast[final_result_corrected_model$Contrast %in% "Not_Colitis_C_A"] <- "No colitis: (Week 6 vs Baseline)"
final_result_corrected_model$Contrast[final_result_corrected_model$Contrast %in% "Colitis_B_A"] <- "Colitis: (Colitis vs Baseline)"
final_result_corrected_model$Contrast[final_result_corrected_model$Contrast %in% "Colitis_C_A"] <- "Colitis: (Week6 vs Baseline)"
final_result_corrected_model$Contrast[final_result_corrected_model$Contrast %in% "Colitis_D_A"] <- "Colitis: (Post vs Baseline)"
final_result_corrected_model$Contrast[final_result_corrected_model$Contrast %in% "Colitis_D_B"] <- "Colitis: (Post vs Colitis)"
final_result_corrected_model$Contrast[final_result_corrected_model$Contrast %in% "Colitis_C_B"] <- "Colitis: (Week 6 vs Colitis)"
final_result_corrected_model$Contrast[final_result_corrected_model$Contrast %in% "Colitis_D_C"] <- "Colitis: (Post vs Week 6)"

final_result_corrected_model$nLog10FDR <- -log10(final_result_corrected_model$adj.P.Val)
#final_result_corrected_model$Significance 

table(final_result_corrected_model$Contrast)

final_result_corrected_model$Contrast <- factor(final_result_corrected_model$Contrast, levels = c("Baseline: (Colitis vs No Colitis)","Week 6: (Colitis vs No Colitis)", "Colitis: (Colitis vs Baseline)","Colitis: (Week 6 vs Colitis)",
                                                                                                  "Colitis: (Week6 vs Baseline)","No colitis: (Week 6 vs Baseline)", "Colitis: (Post vs Baseline)","Colitis: (Post vs Colitis)","Colitis: (Post vs Week 6)") )


out_mat <- pivot_wider(data=final_result_corrected_model[,c("Proteins","Contrast","logFC")], names_from=Proteins, values_from=logFC) # %>% unnest() 
out_mat <- as.data.frame(out_mat)
rownames(out_mat) <- out_mat$Contrast
colnames(out_mat)
###
out <- pheatmap(as.data.frame(out_mat[,-1]))
dev.off()
dev.off()
###
testXYY <- final_result_corrected_model
testXYY$Contrast <- factor(testXYY$Contrast, levels=out$tree_row$labels[out$tree_row$order])
testXYY$Proteins <- factor(testXYY$Proteins, levels=out$tree_col$labels[out$tree_col$order])
###
pdf(file="u01_olink_colitis_figures/mixed_linear_model_reviewed_results_study_correction.pdf",width = 4.5,height = 7.5)
ggplot(testXYY[testXYY$adj.P.Val<0.05,], aes(Proteins, Contrast, color= logFC, size=nLog10FDR )) + 
  geom_point()+ #geom_point(shape=16)+ 
  geom_text(aes(label=round(logFC,1)),size=2,color="black") +
  scale_colour_gradient2(low='steelblue', mid="white",high='firebrick', name='LogFC') + #,breaks = c(-2.8, 0, 1, 5)) +
  theme_bw() +rotate_x_text(angle=25) + coord_flip()+ 
  labs(x ='', y='', title='Study Adjusted Model') 
dev.off()
###

rm(testXYY,out,out_mat)

#pdf(file="u01_olink_colitis_figures/mixed_linear_model_reviewed_results_study_correction_nominal_significance.pdf",width = 8,height = 12)
#ggplot(final_result_corrected_model[final_result_corrected_model$P.Value<0.05,], aes(Proteins, Contrast, color= logFC, size=nLog10FDR )) + 
#  geom_point()+ #geom_point(shape=16)+ 
#  geom_text(aes(label=round(logFC,1)),size=2,color="black") +
#  scale_colour_gradient2(low='steelblue', mid="white",high='firebrick', name='LogFC') + #,breaks = c(-2.8, 0, 1, 5)) +
#  theme_bw() +rotate_x_text(angle=25) + coord_flip()+ 
#  labs(x ='', y='', title='') 
#dev.off()

#############################################################################################################################################
### SAMPLES TABLE
#############################################################################################################################################

dim(colitis_data_matrix_I)
dim(annotation_data)
###
annotation_data <- annotation_data[order(annotation_data$ID),]
colitis_data_matrix_I <- colitis_data_matrix_I[,order(colnames(colitis_data_matrix_I))]
###
identical(rownames(annotation_data),colnames(colitis_data_matrix_I))
###
### table(annotation_data$Study,annotation_data$CT)
my_temporal_matrix <- as.data.frame(table(annotation_data$Study,annotation_data$CT))
my_temporal_matrix$Var2 <- as.character(my_temporal_matrix$Var2)
my_temporal_matrix$Var2[which(my_temporal_matrix$Var2 %in% "No___A")] <- "No Col. Baseline"
my_temporal_matrix$Var2[which(my_temporal_matrix$Var2 %in% "Yes___A")] <- "Col. Baseline"
my_temporal_matrix$Var2[which(my_temporal_matrix$Var2 %in% "Yes___B")] <- "Colitis"
my_temporal_matrix$Var2[which(my_temporal_matrix$Var2 %in% "Yes___C")] <- "Col. Week 6"
my_temporal_matrix$Var2[which(my_temporal_matrix$Var2 %in% "No___C")] <- "No Col. Week 6"
my_temporal_matrix$Var2[which(my_temporal_matrix$Var2 %in% "Yes___D")] <- "Col. Post"

my_temporal_matrix$Var2 <- factor(my_temporal_matrix$Var2, levels=c("No Col. Baseline","Col. Baseline","No Col. Week 6","Col. Week 6","Colitis","Col. Post") )

pdf(file="u01_olink_colitis_figures/mixed_linear_model_reviewed_samples_table_final.pdf",width = 4,height = 3.5)
ggplot(data=my_temporal_matrix,aes(x=Var1,y=Var2,fill=Freq)) + geom_tile() + 
  scale_fill_viridis() + 
  geom_text(aes(label=Freq)) +
  theme_bw() +rotate_x_text(angle=25) +
  labs(x ='', y='', title='Samples Table',fill="#n") 
dev.off()

#############################################################################################################################################
### full model to show (per study)
#############################################################################################################################################

identical(rownames(annotation_data), colnames(colitis_data_matrix_I))
annotation_data <- annotation_data[ match(colnames(colitis_data_matrix_I), rownames(annotation_data)), ]
identical(rownames(annotation_data), colnames(colitis_data_matrix_I))

study_levels <- levels(as.factor(annotation_data$Study))
final_result_per_study_model <- list()

for ( i in 1:length(study_levels)) {
  
  annotation_data_X2  <- annotation_data[which(annotation_data$Study %in% study_levels[i]),]
  colitis_data_matrix_IX2  <- colitis_data_matrix_I[,which(annotation_data$Study %in% study_levels[i])]
  
# comparisons
design <- model.matrix( ~ 0 + CT, data = annotation_data_X2 ) ; colnames(design) <- gsub("CT","",colnames(design))
contr.matrix <- makeContrasts(  "Colitis_A" = Yes___A - No___A, 
                                "Colitis_C" = Yes___C - No___C, 
                                "Not_Colitis_C_A" = No___C - No___A, 
                                "Colitis_B_A" = Yes___B - Yes___A,
                                "Colitis_C_A" = Yes___C - Yes___A,
                                "Colitis_D_A" = Yes___D - Yes___A,
                                "Colitis_C_B" = Yes___C - Yes___B,
                                "Colitis_D_B" = Yes___D - Yes___B,
                                "Colitis_D_C" = Yes___D - Yes___C,
                                levels = colnames(design) )

vfit <- lmFit(colitis_data_matrix_IX2, design) ; vfit <- contrasts.fit(vfit, contrasts=contr.matrix)
efit <- eBayes(vfit, robust = TRUE) ; summary(decideTests(efit))
### final result
temporal_results_file<-list()
for ( iii in 1:length(colnames(summary(decideTests(efit))))){
  temporal_results_file[[iii]] <- topTable(efit, coef=colnames(summary(decideTests(efit)))[iii], n=Inf, adjust.method="BH")
  temporal_results_file[[iii]]$Contrast <- colnames(summary(decideTests(efit)))[iii]
  temporal_results_file[[iii]]$Proteins <- rownames(temporal_results_file[[iii]]) }

temporal_results_file <- do.call(rbind,temporal_results_file)
temporal_results_file$Contrast_long <- temporal_results_file$Contrast
temporal_results_file$Contrast[temporal_results_file$Contrast %in% "Colitis_A"] <- "Baseline: (Colitis vs No Colitis)"
temporal_results_file$Contrast[temporal_results_file$Contrast %in% "Colitis_C"] <- "Week 6: (Colitis vs No Colitis)"
temporal_results_file$Contrast[temporal_results_file$Contrast %in% "Not_Colitis_C_A"] <- "No colitis: (Week 6 vs Baseline)"
temporal_results_file$Contrast[temporal_results_file$Contrast %in% "Colitis_B_A"] <- "Colitis: (Colitis vs Baseline)"
temporal_results_file$Contrast[temporal_results_file$Contrast %in% "Colitis_C_A"] <- "Colitis: (Week6 vs Baseline)"
temporal_results_file$Contrast[temporal_results_file$Contrast %in% "Colitis_D_A"] <- "Colitis: (Post vs Baseline)"
temporal_results_file$Contrast[temporal_results_file$Contrast %in% "Colitis_D_B"] <- "Colitis: (Post vs Colitis)"
temporal_results_file$Contrast[temporal_results_file$Contrast %in% "Colitis_C_B"] <- "Colitis: (Week 6 vs Colitis)"
temporal_results_file$Contrast[temporal_results_file$Contrast %in% "Colitis_D_C"] <- "Colitis: (Post vs Week 6)"
temporal_results_file$study <- study_levels[i]
final_result_per_study_model[[i]]<- temporal_results_file
}

final_result_per_study_model <- do.call(rbind,final_result_per_study_model)
final_result_per_study_model$nLog10FDR <- -log10(final_result_per_study_model$adj.P.Val)


pdf(file="u01_olink_colitis_figures/mixed_linear_model_reviewed_results_per_study_part1.pdf",width = 4.55,height = 6)
ix <- which(final_result_per_study_model$adj.P.Val<0.05 & final_result_per_study_model$study %in% study_levels[1])
ggplot(final_result_per_study_model[ix,], aes(Proteins, Contrast, color= logFC, size=nLog10FDR )) + 
  geom_point()+ #geom_point(shape=16)+ 
  geom_text(aes(label=round(logFC,1)),size=2,color="black") +
  scale_colour_gradient2(low='steelblue', mid="white",high='firebrick', name='LogFC') + #,breaks = c(-2.8, 0, 1, 5)) +
  theme_bw() +rotate_x_text(angle=25) + coord_flip()+ 
  labs(x ='', y='', title=paste("Study:",study_levels[1])) 
dev.off()


pdf(file="u01_olink_colitis_figures/mixed_linear_model_reviewed_results_per_study_part2.pdf",width = 8,height = 6)
ix <- which(final_result_per_study_model$adj.P.Val<0.05 & final_result_per_study_model$study %in% study_levels[2])
ggplot(final_result_per_study_model[ix,], aes(Proteins, Contrast, color= logFC, size=nLog10FDR )) + 
  geom_point()+ #geom_point(shape=16)+ 
  geom_text(aes(label=round(logFC,1)),size=2,color="black") +
  scale_colour_gradient2(low='steelblue', mid="white",high='firebrick', name='LogFC') + #,breaks = c(-2.8, 0, 1, 5)) +
  theme_bw() +rotate_x_text(angle=25) + coord_flip()+ 
  labs(x ='', y='', title=paste("Study:",study_levels[2])) +
  theme(text=element_text(size=21)) +
  theme(plot.margin=unit(c(t=5,r=5,b=5,l=60),"pt"))
dev.off()

pdf(file="u01_olink_colitis_figures/mixed_linear_model_reviewed_results_per_study_part3.pdf",width = 4,height = 4)
ix <- which(final_result_per_study_model$adj.P.Val<0.05 & final_result_per_study_model$study %in% study_levels[3])
ggplot(final_result_per_study_model[ix,], aes(Proteins, Contrast, color= logFC, size=nLog10FDR )) + 
  geom_point()+ #geom_point(shape=16)+ 
  geom_text(aes(label=round(logFC,1)),size=2,color="black") +
  scale_colour_gradient2(low='steelblue', mid="white",high='firebrick', name='LogFC') + #,breaks = c(-2.8, 0, 1, 5)) +
  theme_bw() +rotate_x_text(angle=25) + coord_flip()+ 
  labs(x ='', y='', title=paste("Study:",study_levels[3])) 
dev.off()

final_result_corrected_model$direction <- final_result_corrected_model$logFC > 0
final_result_corrected_model$direction[final_result_corrected_model$direction %in% "TRUE"] <- "up"
final_result_corrected_model$direction[final_result_corrected_model$direction %in% "FALSE"] <- "down"

final_result_corrected_model$significance <- final_result_corrected_model$adj.P.Val<0.05
final_result_corrected_model$significance[final_result_corrected_model$significance %in% "TRUE"] <- "sig"
final_result_corrected_model$significance[final_result_corrected_model$significance %in% "FALSE"] <- "not.sig"

ix <-which(final_result_corrected_model$Contrast %in% "Colitis: (Colitis vs Baseline)")
#pdf(file="u01_olink_colitis_figures/mixed_linear_model_reviewed_results_study_correction_volcano1.pdf",width = 8,height = 12)
ggplot(final_result_corrected_model[ix,], aes(logFC,nLog10FDR,color=direction)) + 
  geom_point()+ 
  scale_color_manual(values = c("firebrick","steelblue")) +
  geom_vline(xintercept = 0, linetype = 2, color = 'black') +
  geom_hline(yintercept = 1.3, linetype = 2, color = 'red') +
  geom_point(data=final_result_corrected_model[final_result_corrected_model$significance %in% "not.sig",],color="grey60")+ 
  theme_bw() + labs(x ='', y='', title="Colitis: (Colitis vs Baseline)")+ 
  geom_label_repel(data=final_result_corrected_model[ix[1:10],],aes(label=Proteins),nudge_x = -2.5,force = 10)+
  theme(legend.position = "none")
#dev.off()

ix <-which(final_result_corrected_model$Contrast %in% "Colitis: (Week6 vs Baseline)")
#pdf(file="u01_olink_colitis_figures/mixed_linear_model_reviewed_results_study_correction_volcano2.pdf",width = 6,height = 6)
ggplot(final_result_corrected_model[ix,], aes(logFC,nLog10FDR,color=direction)) + 
  geom_point()+ 
  scale_color_manual(values = c("firebrick","steelblue")) +
  geom_vline(xintercept = 0, linetype = 2, color = 'black') +
  geom_hline(yintercept = 1.3, linetype = 2, color = 'red') +
  geom_point(data=final_result_corrected_model[final_result_corrected_model$significance %in% "not.sig",],color="grey60")+ 
  theme_bw() + labs(x ='', y='', title="Colitis: (Week6 vs Baseline)")+ 
  geom_label_repel(data=final_result_corrected_model[ix[1:10],],aes(label=Proteins),nudge_x = -2.5,force = 10)+
  theme(legend.position = "none")
#dev.off()

ix <-which(final_result_corrected_model$Contrast %in% "No colitis: (Week 6 vs Baseline)")
#pdf(file="u01_olink_colitis_figures/mixed_linear_model_reviewed_results_study_correction_volcano3.pdf",width = 6,height = 6)
ggplot(final_result_corrected_model[ix,], aes(logFC,nLog10FDR,color=direction)) + 
  geom_point()+ 
  scale_color_manual(values = c("firebrick","steelblue")) +
  geom_vline(xintercept = 0, linetype = 2, color = 'black') +
  geom_hline(yintercept = 1.3, linetype = 2, color = 'red') +
  geom_point(data=final_result_corrected_model[final_result_corrected_model$significance %in% "not.sig",],color="grey60")+ 
  theme_bw() + labs(x ='', y='', title="No colitis: (Week 6 vs Baseline)")+ 
  geom_label_repel(data=final_result_corrected_model[ix[1:10],],aes(label=Proteins),nudge_x = -2.5,force = 10)+
  theme(legend.position = "none")
#dev.off()

ix <-which(final_result_corrected_model$Contrast %in% "Colitis: (Post vs Baseline)")
#pdf(file="u01_olink_colitis_figures/mixed_linear_model_reviewed_results_study_correction_volcano4.pdf",width = 6,height = 6)
ggplot(final_result_corrected_model[ix,], aes(logFC,nLog10FDR,color=direction)) + 
  geom_point()+ 
  scale_color_manual(values = c("firebrick","steelblue")) +
  geom_vline(xintercept = 0, linetype = 2, color = 'black') +
  geom_hline(yintercept = 1.3, linetype = 2, color = 'red') +
  geom_point(data=final_result_corrected_model[final_result_corrected_model$significance %in% "not.sig",],color="grey60")+ 
  theme_bw() + labs(x ='', y='', title="Colitis: (Post vs Baseline)")+ 
  geom_label_repel(data=final_result_corrected_model[ix[1:10],],aes(label=Proteins),nudge_x = -2.5,force = 10)+
  theme(legend.position = "none")
#dev.off()

#############################################################################################################################################
### stack barplot
#############################################################################################################################################

stacked_bar_plot <- final_result_corrected_model

stacked_bar_plot$logFC_dir <- stacked_bar_plot$logFC > 0 
stacked_bar_plot$logFC_dir[stacked_bar_plot$logFC > 0] <- "Up"
stacked_bar_plot$logFC_dir[stacked_bar_plot$logFC < 0] <- "Down"
stacked_bar_plot$sig <- stacked_bar_plot$adj.P.Val<0.05

stacked_bar_plot <- data.frame(table(stacked_bar_plot$Contrast, stacked_bar_plot$logFC_dir, stacked_bar_plot$sig))
colnames(stacked_bar_plot) <- c('Contrast', 'logFC_dir','sig', 'count')

regional_dmr_summary <- stacked_bar_plot %>%  filter(sig=='TRUE') %>% dplyr::select(-sig)
regional_dmr_summary$count[regional_dmr_summary$logFC_dir=='Down'] <- 0 - regional_dmr_summary$count[regional_dmr_summary$logFC_dir=='Down']

regional_dmr_summary$Contrast <- as.character(regional_dmr_summary$Contrast) 
regional_dmr_summary$Contrast <- factor(regional_dmr_summary$Contrast, levels = c("Baseline: (Colitis vs No Colitis)","Week 6: (Colitis vs No Colitis)", "Colitis: (Colitis vs Baseline)","Colitis: (Week 6 vs Colitis)",
                                                                                                  "Colitis: (Week6 vs Baseline)","No colitis: (Week 6 vs Baseline)", "Colitis: (Post vs Baseline)","Colitis: (Post vs Colitis)","Colitis: (Post vs Week 6)") )

#pdf(file="u01_olink_colitis_figures/mixed_linear_model_reviewed_results_study_correction_summary.pdf",width = 5,height = 3)
ggplot(regional_dmr_summary) + 
  aes(x=count, y=Contrast, fill=logFC_dir, label=abs(count)) + geom_bar(stat='identity', width=0.5) +
  theme_classic() + theme(axis.text.y = element_text(face='bold'),
                          axis.text.x = element_text(face='bold'),
                          axis.title = element_text(face = "bold"),
                          legend.title= element_text(face='bold'),
                          strip.text=element_text(face='bold'),
                          plot.title=element_text(face='bold')) + 
  scale_fill_manual(values=c('steelblue', 'firebrick'), name='LogFC Direction') +
  labs(x='no.DEGs(FDR<0.05)', y='Contrast') + theme(legend.position="bottom")
#dev.off()

rm(regional_dmr_summary,stacked_bar_plot)
#dev.off()

#############################################################################################################################################
### Upset plot for model using study adjusted data
#############################################################################################################################################

my_list <- list()
test <- final_result_corrected_model[which(final_result_corrected_model$adj.P.Val<0.05),]
test$Contrast <- as.character(test$Contrast)
for (i in 1:length(levels(as.factor(test$Contrast)))){
  ix <- which(test$Contrast %in% levels(as.factor(test$Contrast))[i])
  my_list[[i]] <- test$Proteins[ix]
  names(my_list)[i] <- levels(as.factor(test$Contrast))[i] }

rm(test)

pdf(file="u01_olink_colitis_figures/mixed_linear_model_reviewed_results_study_correction_upset.pdf",width = 8,height = 4)
upset(fromList(my_list), order.by = c("freq"), #"degree",
      number.angles = 0, nsets=10, 
      mainbar.y.label = "Intersection",
      empty.intersections = "on", 
      sets.bar.color = "black", #
      cutoff = 5, #
      nintersects = 30,
      matrix.color = "darkred", keep.order = TRUE,
      main.bar.color = "darkblue", #mb.ratio = c(0.7, 0.55),
      text.scale = c(3, 2, 1.5, 1.5, 2, 2),
      point.size = 4, line.size = 2)
dev.off()

#############################################################################################################################################
### Upset plot for models  PER study  
#############################################################################################################################################

my_list <- list()
test <- final_result_per_study_model[which(final_result_per_study_model$adj.P.Val<0.05),]
test$Contrast <- as.character(test$Contrast)
test$study <- as.character(test$study)
test$condition <- paste( test$study, test$Contrast, sep=" " ) 

for (i in 1:length(levels(as.factor(test$condition)))){
  ix <- which(test$condition %in% levels(as.factor(test$condition))[i])
  my_list[[i]] <- test$Proteins[ix]
  names(my_list)[i] <- levels(as.factor(test$condition))[i]
  } ; rm(test)
  
names(my_list)

pdf(file="u01_olink_colitis_figures/mixed_linear_model_reviewed_results_PER_study_upset.pdf",width = 12,height = 8)
upset(fromList(my_list), order.by = c("freq"), #"degree",
      number.angles = 0, nsets=10, 
      mainbar.y.label = "Intersection",
      empty.intersections = "on", 
      sets.bar.color = "black", #
      cutoff = 5, #
      nintersects = 30,
      matrix.color = "darkred", keep.order = TRUE,
      main.bar.color = "darkblue", #mb.ratio = c(0.7, 0.55),
      text.scale = c(3, 2, 1.5, 1.5, 2, 2),
      point.size = 4, line.size = 2)
dev.off()

#############################################################################################################################################
### Violin plot selection for ADJUSTED DATA
#############################################################################################################################################

my_list <- list()
test <- final_result_corrected_model[which(final_result_corrected_model$adj.P.Val<0.05),]
test$Contrast <- as.character(test$Contrast)
for (i in 1:length(levels(as.factor(test$Contrast)))){
  ix <- which(test$Contrast %in% levels(as.factor(test$Contrast))[i])
  my_list[[i]] <- test$Proteins[ix]
  names(my_list)[i] <- levels(as.factor(test$Contrast))[i] }

rm(test)

library(UpSetR)
library(ComplexHeatmap)
m = make_comb_mat(my_list)
set_size(m)

to_test <- list( C0 = extract_comb(m,'1000'), #1
                 C1 = extract_comb(m,'0100') ,#2
                 C2 = extract_comb(m,'0010') ,#3
                 C3 = extract_comb(m,'0001') ,#4
                 C5 = extract_comb(m,'1110') ) #5
### names(to_test) <- names(my_list)

my_proteins_to_test <- unique(sort(unlist(to_test[5])))


ix <- which(rownames(colitis_data_matrix_I_study_corrected) %in% my_proteins_to_test)
for_plot_heatmap <- colitis_data_matrix_I_study_corrected[ix,]
annotation_data$ID <- rownames(annotation_data)

#identical(colnames(for_plot_heatmap),annotation_data$ID)
#dim(for_plot_heatmap)
#dim(annotation_data)

melted_for_plot_heatmap <- melt(for_plot_heatmap)
colnames(melted_for_plot_heatmap) <- c("Protein","ID","value")

melted_for_plot_heatmap <- merge(melted_for_plot_heatmap,annotation_data,by="ID")

table(melted_for_plot_heatmap$CT)
table(melted_for_plot_heatmap$Protein)
### melted_for_plot_heatmap$status <- as.character(melted_for_plot_heatmap$CT)

melted_for_plot_heatmap$Status <- paste(melted_for_plot_heatmap$Timepoint,melted_for_plot_heatmap$Colitis,sep="___")
table(melted_for_plot_heatmap$Status)

melted_for_plot_heatmap$Status <- gsub("A_","Baseline_",melted_for_plot_heatmap$Status)
melted_for_plot_heatmap$Status <- gsub("D_","Post_",melted_for_plot_heatmap$Status)
melted_for_plot_heatmap$Status <- gsub("B_","Colitis_",melted_for_plot_heatmap$Status)
melted_for_plot_heatmap$Status <- gsub("C_","Week6_",melted_for_plot_heatmap$Status)

melted_for_plot_heatmap$Status <- factor(melted_for_plot_heatmap$Status, levels=c("Baseline___No","Baseline___Yes","Week6___No","Week6___Yes","Colitis___Yes","Post___Yes"))

my_comparisons <- list( c("Baseline___No","Baseline___Yes"),
                        c("Baseline___No","Week6___No"),
                        c("Baseline___Yes","Week6___Yes"),
                        c("Week6___No","Week6___Yes"),
                        c("Week6___No","Colitis___Yes"))

melted_for_plot_heatmap$Protein <- as.character(melted_for_plot_heatmap$Protein)
melted_for_plot_heatmap$Protein <- factor(melted_for_plot_heatmap$Protein, levels=my_proteins_to_test)

pdf(file="u01_olink_colitis_figures/mixed_linear_model_reviewed_results_study_correction_violin_plots_intersection_colitis.pdf", width = 12, height = 7)
ggviolin(data=melted_for_plot_heatmap, x="Status", y="value", color="Colitis",fill="Colitis",alpha=0.5) + theme_bw() +
  facet_wrap(~Protein,scales = "free",ncol=4) +
  geom_boxplot(width=0.2) + geom_quasirandom(varwidth = TRUE,alpha=0.3,size=0.5) +
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test", paired=FALSE, label="p.signif", hide.ns = FALSE) +
  theme(axis.text.x = element_text(angle = 45,hjust = 1)) + theme(legend.position="bottom") + theme(#axis.title.x=element_blank(),
    #axis.text.x=element_blank(),
    #axis.ticks.x=element_blank(),
    strip.text=element_text(size=rel(1.5),face="bold"),
    legend.text = element_text(size=rel(1.5),face="bold"),
    legend.title = element_text(size=rel(1.5),face="bold"),
    axis.text.y = element_text(size=rel(1.5),face="bold"),
    axis.title.y = element_text(size=rel(1.5),face="bold")) +
  scale_fill_manual(values = c("firebrick","steelblue")) +
  scale_color_manual(values = c("firebrick","steelblue")) +
  ggtitle("") + labs(x = "", y = "Study Adj. NPX values") +
  theme(text = element_text(size=rel(4.8)),
        legend.text = element_text(size=rel(4.8)),
        legend.title = element_text(size=rel(4.8)),
        plot.subtitle = element_text(size=rel(4.8)),
        axis.title.x = element_text(size=rel(4.8)),
        axis.text.y = element_text(size=rel(4.8)),
        strip.text = element_text(size=rel(4.8)),
        axis.title.y = element_text(size=rel(4.8)),
        plot.title=element_text(size=rel(4.8)) ) +
  theme(plot.margin=unit(c(t=5,r=5,b=5,l=100),"pt"))
dev.off()

#############################################################################################################################################
### Wilcox test to verify linear models
#############################################################################################################################################
###
identical(colitis_metadata$sample_id,colnames(colitis_data_matrix))
### copy and subset
annotation_data <- colitis_metadata
###
ix <- which(! annotation_data$Study %in% "control")
colitis_data_matrix_I <- colitis_data_matrix[,ix]
annotation_data <- annotation_data[ix,]
### subset (remove warnings from experiment)
ix <- which(! annotation_data$QC %in% "Warning")
colitis_data_matrix_I <- colitis_data_matrix_I[,ix]
annotation_data <- annotation_data[ix,]
###
identical(annotation_data$sample_id,colnames(colitis_data_matrix_I))
###

dim(annotation_data)
dim(colitis_data_matrix_I)

### loop for wilcox/man whitney
my_results_wilcox <- data.frame()
for (x in 1:nrow(colitis_data_matrix_I)) {
 for ( i in 1:length(levels(as.factor(annotation_data$Study)))) {
  a <- which(annotation_data$Colitis_Timepoint %in% "Yes___C" & annotation_data$Study %in% levels(as.factor(annotation_data$Study))[i])
  b <- which(annotation_data$Colitis_Timepoint %in% "No___C" & annotation_data$Study %in% levels(as.factor(annotation_data$Study))[i])
  test <- data.frame( p.value = wilcox.test( as.numeric(colitis_data_matrix_I[x,a]),as.numeric(colitis_data_matrix_I[x,b]))$p.value,
                      Protein = rownames(colitis_data_matrix_I)[x],
                      TimePoint = "C",
                      Study = levels(as.factor(annotation_data$Study))[i] )
  my_results_wilcox <- rbind(my_results_wilcox,test)
  
  a <- which(annotation_data$Colitis_Timepoint %in% "Yes___A" & annotation_data$Study %in% levels(as.factor(annotation_data$Study))[i])
  b <- which(annotation_data$Colitis_Timepoint %in% "No___A" & annotation_data$Study %in% levels(as.factor(annotation_data$Study))[i])
  test <- data.frame( p.value = wilcox.test( as.numeric(colitis_data_matrix_I[x,a]),as.numeric(colitis_data_matrix_I[x,b]))$p.value,
                      Protein = rownames(colitis_data_matrix_I)[x],
                      TimePoint = "A",
                      Study = levels(as.factor(annotation_data$Study))[i] )
  my_results_wilcox <- rbind(my_results_wilcox,test)
  
}}

table(my_results_wilcox$p.value<0.05)
table(my_results_wilcox$p.value<0.01)

table(my_results_wilcox$Study[which(my_results_wilcox$p.value<0.05)])
my_results_wilcox[which(my_results_wilcox$p.value<0.01),]

test <- list()
for (i in 1:nrow(my_results_wilcox)){test[[i]] <- p.adjust(my_results_wilcox$p.value[i], n = 92, method = 'fdr') }
my_results_wilcox$adj.p.value <- unlist(test)
table(my_results_wilcox$adj.p.value<0.1)

my_results_wilcox <- my_results_wilcox[order(my_results_wilcox$p.value),]

write.csv(file="wilcox.test.results.csv",my_results_wilcox)

my_wilcox_top_proteins <- unique(sort(my_results_wilcox$Protein[which(my_results_wilcox$p.value < 0.05)]))

#############################################################################################################################################
### Wilcox test results
#############################################################################################################################################

test<- my_results_wilcox[which(my_results_wilcox$p.value<0.05),]

test$Condition <- paste(test$TimePoint,test$Study,sep="___")

out_mat <- pivot_wider(data=test[,c("Protein","Condition","p.value")], names_from=Protein, values_from=p.value) # %>% unnest() 
out_mat <- as.data.frame(out_mat)
rownames(out_mat) <- out_mat$Condition
colnames(out_mat)
out_mat[is.na(out_mat)]<-0
###
out <- pheatmap(as.data.frame(out_mat[,-1]))
dev.off()
dev.off()
###
testXYY <- test
out$tree_row$labels
testXYY$Condition
###
out$tree_row$labels<-gsub("A_",'Baseline',out$tree_row$labels)
testXYY$Condition<-gsub("A_",'Baseline',testXYY$Condition)

out$tree_row$labels<-gsub("C_",'Week6',out$tree_row$labels)
testXYY$Condition<-gsub("C_",'Week6',testXYY$Condition)

out$tree_row$labels<-gsub("__",'\n',out$tree_row$labels)
testXYY$Condition<-gsub("__",'\n',testXYY$Condition)

testXYY$Condition <- factor(testXYY$Condition, levels=out$tree_row$labels[out$tree_row$order])
testXYY$Protein <- factor(testXYY$Protein, levels=out$tree_col$labels[out$tree_col$order])

pdf(file="u01_olink_colitis_figures/p.value.significant.markers.wilcox.pdf", width = 16, height = 5)
ggplot(testXYY, aes(Protein, Condition, color= p.value, size=p.value)) + 
  geom_point()+ #geom_point(shape=16)+ 
  #geom_text(aes(label=round(logFC,1)),size=2,color="black") +
  scale_colour_viridis() +
  theme_bw() +rotate_x_text(angle=45) + 
  labs(x ='', y='', title='Colitis vs Non Colitis') + theme(text = element_text(size=rel(4.8)), 
                                                            legend.text = element_text(size=rel(4.8)),
                                                            plot.title=element_text(size=rel(4.8)),
                                                            axis.text.x = element_text(size=rel(2.8)))
dev.off()


#############################################################################################################################################
### Violin plot figures
#############################################################################################################################################

#head(my_results_wilcox)
#table(my_wilcox_top_proteins %in%  my_top_proteins_global )
#table(my_wilcox_top_proteins %in%  my_top_proteins_local )
#table(my_top_proteins_global %in% my_top_proteins_local)
#my_results_wilcox$Protein[1:5]

my_position <- c(which(my_results_wilcox$TimePoint %in% "A" & my_results_wilcox$Study %in% my_studies[1] & my_results_wilcox$p.value<0.05)[1:2],
which(my_results_wilcox$TimePoint %in% "A" & my_results_wilcox$Study %in% my_studies[2] & my_results_wilcox$p.value<0.05 )[1:2],
which(my_results_wilcox$TimePoint %in% "A" & my_results_wilcox$Study %in% my_studies[3] & my_results_wilcox$p.value<0.05)[1:2],
which(my_results_wilcox$TimePoint %in% "C" & my_results_wilcox$Study %in% my_studies[1] & my_results_wilcox$p.value<0.05)[1:2],
which(my_results_wilcox$TimePoint %in% "C" & my_results_wilcox$Study %in% my_studies[2] & my_results_wilcox$p.value<0.05)[1:2],
which(my_results_wilcox$TimePoint %in% "C" & my_results_wilcox$Study %in% my_studies[3] & my_results_wilcox$p.value<0.05)[1:2])


my_results_wilcox$Protein[my_position]

ix <- which(rownames(colitis_data_matrix_iio) %in% c("CD28","IL2...27","FGF-21","TNF","IL10","MMP-10","CCL11"))
rownames(colitis_data_matrix_iio)[ix]

### prepare for violin plots
my_melted_df <- melt(as.matrix(colitis_data_matrix_iio[ix,]))
colnames(my_melted_df) <- c("Protein","sample_id","value")
my_melted_df <- merge(my_melted_df,colitis_metadata,by="sample_id")
my_melted_df$Protein <- as.character(my_melted_df$Protein)
my_melted_df$Protein <- tidyr::separate(data.frame(my_melted_df$Protein), 1, sep="\\.\\.\\.", c("a","b","c"))$a
my_melted_df$Protein <- as.factor(my_melted_df$Protein)
table(my_melted_df$Protein)

#ix <- which(my_melted_df$Study %in% "13-075, 17-162, 13-169")
#my_melted_df <-my_melted_df[ix,]

my_melted_df$Status <- paste(my_melted_df$Timepoint,my_melted_df$Colitis,sep="___")
my_melted_df$Status <- gsub("A_","Baseline_",my_melted_df$Status)
my_melted_df$Status <- gsub("D_","Post_",my_melted_df$Status)
my_melted_df$Status <- gsub("B_","Colitis_",my_melted_df$Status)
my_melted_df$Status <- gsub("C_","Week6_",my_melted_df$Status)

my_melted_df$Status <- factor(my_melted_df$Status, levels=c("Baseline___No","Baseline___Yes","Week6___No","Week6___Yes","Colitis___Yes","Post___Yes"))

my_comparisons <- list( c("Baseline___No","Baseline___Yes"),
                        #c("Baseline___No","Week6___No"),
                        #c("Baseline___Yes","Week6___Yes"),
                        c("Week6___No","Week6___Yes")
                        #c("Week6___No","Colitis___Yes")
                        )

pdf(file="u01_olink_colitis_figures/boxplots_clustered_significant_markers_colitis_wilcox_npx_raw.pdf ",width = 18,height = 16)
ggviolin(data=my_melted_df, x="Status", y="value", color="Colitis",fill="Colitis",alpha=0.5) + theme_bw() +
  facet_wrap(~Study+Protein,scales = "free",ncol=7) +
  geom_boxplot(width=0.2) + geom_quasirandom(varwidth = TRUE,alpha=0.3,size=0.5) +
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test", paired=FALSE, label="p.signif", hide.ns = FALSE) +
  theme(axis.text.x = element_text(angle = 45,hjust = 1)) + theme(legend.position="bottom") + theme(
    #axis.title.x=element_blank(),
    #axis.text.x=element_blank(),
    #axis.ticks.x=element_blank(),
    strip.text=element_text(size=rel(1.5),face="bold"),
    legend.text = element_text(size=rel(1.5),face="bold"),
    legend.title = element_text(size=rel(1.5),face="bold"),
    axis.text.y = element_text(size=rel(1.5),face="bold"),
    axis.title.y = element_text(size=rel(1.5),face="bold")) +
  scale_fill_manual(values = c("firebrick","steelblue")) +
  scale_color_manual(values = c("firebrick","steelblue")) +
  ggtitle("") + labs(x = "", y = "Raw NPX values") +
  theme(text = element_text(size=rel(4.8)),
        legend.text = element_text(size=rel(4.8)),
        legend.title = element_text(size=rel(4.8)),
        plot.subtitle = element_text(size=rel(4.8)),
        axis.title.x = element_text(size=rel(4.8)),
        axis.text.y = element_text(size=rel(4.8)),
        strip.text = element_text(size=rel(4.8)),
        axis.title.y = element_text(size=rel(4.8)),
        plot.title=element_text(size=rel(4.8)) ) +
  theme(plot.margin=unit(c(t=5,r=5,b=5,l=100),"pt"))
dev.off()

#############################################################################################################################################
### 
#############################################################################################################################################

#write.csv(file="stats_final_result_reviewed_models.csv",final_result_corrected_model)
#table( !(final_result_corrected_model$adj.P.Val<0.05 & final_result_corrected_model$Contrast %in% "No colitis: (Week 6 vs Baseline)") & (final_result_corrected_model$adj.P.Val<0.05 & final_result_corrected_model$Contrast %in% "Colitis: (Week6 vs Baseline)"))
#ix<-which( !(final_result_corrected_model$adj.P.Val<0.05 & final_result_corrected_model$Contrast %in% "No colitis: (Week 6 vs Baseline)") & (final_result_corrected_model$adj.P.Val<0.05 & final_result_corrected_model$Contrast %in% "Colitis: (Week6 vs Baseline)"))
#final_result_corrected_model[ix,]

#############################################################################################################################################
### 
#############################################################################################################################################

#############################################################################################################################################
### 
#############################################################################################################################################

#############################################################################################################################################
### STORAGE
#############################################################################################################################################
### save.image("u01_RData_files/temporal_analysis.RData")
### load("u01_RData_files/temporal_analysis.RData")
#############################################################################################################################################
### The end. (for now)
#############################################################################################################################################