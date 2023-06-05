library(Seurat)
library(monocle3)
library(monocle)
library(cowplot)
library(patchwork)
library(dplyr)
library(ggplot2)
library(stringr)
library(data.table)


# pseudotime

E18_Mesen <- readRDS(file = 'E18_Mesen.rds')

Mesen<- subset(E18_Mesen, idents = c('ASM', 'MyoFB'))
Mesen$newtype <- Idents(Mesen)
expdata = Mesen@assays$RNA@counts
cell_metadata=Mesen@meta.data
gene_annotation =as.data.frame(row.names(expdata))
rownames(gene_annotation) =row.names(expdata)
colnames(gene_annotation) ='id'
gene_annotation$gene_short_name = gene_annotation$id

cds <-new_cell_data_set(expdata,cell_metadata = cell_metadata,gene_metadata = gene_annotation)

reducedDims(cds)$PCA <-as.data.frame(Mesen@reductions$Embeded_z0.6@cell.embeddings)

reducedDims(cds)$UMAP <-as.data.frame(Mesen@reductions$umap0.6@cell.embeddings)

cds <-cluster_cells(cds,reduction_method ='UMAP')

cds <-learn_graph(cds,use_partition =FALSE)
cds <- order_cells(cds)
cds <-order_cells(cds,root_cells =rownames(subset(cds@colData,desc_0.6=='9')))

pdf(file = 'pseudotime.pdf', height = 6,width = 7)
plot_cells(cds,
           color_cells_by ="pseudotime",
           label_cell_groups= T,
           label_leaves=T,
           label_branch_points=T,
           graph_label_size=1.5,
           label_groups_by_cluster = T)
dev.off()




