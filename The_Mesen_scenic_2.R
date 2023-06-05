#running pyscenic and got the E18_Mesen.loom file
library(SCENIC)
library(SCopeLoomR)
library(AUCell)
library(Seurat)
library(pheatmap)
library(SummarizedExperiment)
setwd("./Mesenchymal/result/")
E18_Mesen <- readRDS(file = 'E18_Mesen_pure.rds')

pyScenicLoomFile <- file.path( "E18_Mesen.loom")
loom <- open_loom(pyScenicLoomFile, mode="r")

# Read information from loom file:
regulons_incidMat <- get_regulons(loom, column.attr.name='Regulons')
regulons <- regulonsToGeneLists(regulons_incidMat)
regulons$`Sox5 (16g)`
regulonsAUC <- get_regulons_AUC(loom, column.attr.name='RegulonsAUC')
regulonsAucThresholds <- get_regulon_thresholds(loom)
tail(regulonsAucThresholds[order(as.numeric(names(regulonsAucThresholds)))])
embeddings <- get_embeddings(loom)
close_loom(loom)

sub_regulonAUC <- regulonsAUC[,match(colnames(E18_Mesen),colnames(regulonsAUC))]
dim(sub_regulonAUC)

#make sure whether identical
identical(colnames(sub_regulonAUC), colnames(E18_Mesen))

cellClusters <- data.frame(row.names = colnames(E18_Mesen), 
                           seurat_clusters = as.character(E18_Mesen$newtype))
cellTypes <- data.frame(row.names = colnames(E18_Mesen), 
                        celltype = E18_Mesen$newtype)
head(cellTypes)
head(cellClusters)
sub_regulonAUC[1:4,1:4] 
save(sub_regulonAUC,cellTypes,cellClusters,E18_Mesen,
     file = 'E18_Mesen_for_rss_and_visual.Rdata')
load('E18_Mesen_for_rss_and_visual.Rdata')
regulonsToPlot = c('Etv1(+)')
regulonsToPlot
E18_Mesen@meta.data = cbind(E18_Mesen@meta.data ,t(assay(sub_regulonAUC[regulonsToPlot,])))
Idents(E18_Mesen) <- E18_Mesen$newtype
table(Idents(E18_Mesen) )
VlnPlot(E18_Mesen, pt.size = 0.0,split.by = 'batch', features =  regulonsToPlot, ncol = 3 )+ NoLegend()


head(cellTypes)
selectedResolution <- "celltype"
cellsPerGroup <- split(rownames(cellTypes), 
                       cellTypes[,selectedResolution]) 
sub_regulonAUC <- sub_regulonsAUC[onlyNonDuplicatedExtended(rownames(sub_regulonsAUC)),] # 去除extened regulons
dim(sub_regulonAUC)
# Calculate average expression:
regulonActivity_byGroup <- sapply(cellsPerGroup,
                                  function(cells) 
                                    rowMeans(getAUC(sub_regulonAUC)[,cells]))

# Scale expression:
regulonActivity_byGroup_Scaled <- t(scale(t(regulonActivity_byGroup),
                                          center = T, scale=T)) 

# Scale processing of the same regulon in different clusters
dim(regulonActivity_byGroup_Scaled)
regulonActivity_byGroup_Scaled=regulonActivity_byGroup_Scaled[]
regulonActivity_byGroup_Scaled=na.omit(regulonActivity_byGroup_Scaled)
pheatmap:pheatmap(regulonActivity_byGroup_Scaled)# Calculate average expression:
rss <- calcRSS(AUC=getAUC(sub_regulonAUC), 
               cellAnnotation=cellTypes[colnames(sub_regulonAUC), 
                                        selectedResolution]) 
rss=na.omit(rss) 
rssPlot <- plotRSS(rss)
pdf(file = 'Mesen_scenic.pdf', height = 10, width = 4)
rssPlot$plot
dev.off()


