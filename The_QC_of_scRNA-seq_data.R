#QC and selecting cells for further analysis
library(Seurat)
library(cowplot)
library(patchwork)
library(dplyr)
library(ggplot2)
library(stringr)
library(data.table)

help1 = c("E18_5WT","E18_5HM")
file1 =str_c(file2_dir,help1,"/outs/filtered_feature_bc_matrix.h5")
print(file1[1])

sample1 <- Read10X_h5(file1[1],use.names = TRUE, unique.features = TRUE)
sample1 <- CreateSeuratObject(counts = sample1, project = "sample1")
print(sample1)
sample1[["percent.mt"]] <- PercentageFeatureSet(sample1, pattern = "^mt-")
VlnPlot(sample1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

sample1 <- subset(sample1, subset = nFeature_RNA > 200 
		  & nFeature_RNA < 6000 & percent.mt < 10)
sample1@meta.data$batch <- help1[1]
saveRDS(sample1,str_c(save_dat,help1[1],"_qc.rds"))


sample1 <- Read10X_h5(file1[2],use.names = TRUE, unique.features = TRUE)
sample1 <- CreateSeuratObject(counts = sample1, project = "sample1")
print(sample1)
sample1[["percent.mt"]] <- PercentageFeatureSet(sample1, pattern = "^mt-")
VlnPlot(sample1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

sample1 <- subset(sample1, subset = nFeature_RNA > 200 
		  & nFeature_RNA < 6000 & percent.mt < 10)
sample1@meta.data$batch <- help1[2]
saveRDS(sample1,str_c(save_dat,help1[2],"_qc.rds"))



