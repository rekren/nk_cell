set_libs <- function() {
  library(Seurat);library(SeuratDisk);library(dplyr);library(tidyverse); library(ggplot2); library(reshape2);library(monocle3);library(clustree);library(SeuratWrappers);library(dittoSeq);
  
}
set_libs()

pbmc_combined <- LoadH5Seurat(file = "mapped_integrated_pbmc_cohort.h5seurat")

DefaultAssay(pbmc_combined)='prediction.score.celltype.l2'
# What are the cell names of all NK cells?
WhichCells(pbmc_combined, idents = c("NK", "NK_CD56bright","NK Proliferating"))

NK_only <- subset(pbmc_combined, idents = c("NK", "NK_CD56bright","NK Proliferating"))
SaveH5Seurat(NK_only,verbose = T,filename = "NK_40.h5seurat")

rm(list=setdiff(ls(), "NK_only"))
gc()

######################
NK_only <- LoadH5Seurat("NK_40.h5seurat")
DefaultAssay(NK_only) <- "integrated"
NK_only <- FindVariableFeatures(NK_only,assay = "integrated")
NK_only <- ScaleData(NK_only)
NK_only <- RunPCA(NK_only)
NK_only <- FindNeighbors(NK_only,force.recalc = T)
NK_only <- FindClusters(NK_only,algorithm = 3)
NK_only <- RunUMAP(NK_only,dims = 1:10)

## We compute the Progeny (https://doi.org/doi:10.18129/B9.bioc.progeny) activity 
## scores and add them to our Seurat object as a new assay called Progeny. 
NK_only <- progeny(NK_only, scale=FALSE, organism="Human", top=500, perm=1, 
                return_assay = TRUE)

## We can now directly apply Seurat functions in our Progeny scores. 
## For instance, we scale the pathway activity scores. 
NK_only <- Seurat::ScaleData(NK_only, assay = "progeny") 

## We transform Progeny scores into a data frame to better handling the results
progeny_scores_df <- 
  as.data.frame(t(GetAssayData(NK_only, slot = "scale.data", 
                               assay = "progeny"))) %>%
  rownames_to_column("Cell") %>%
  gather(Pathway, Activity, -Cell) 

## We create a data frame with the specification of the cells that belong to 
## each cluster to match with the Progeny scores.
CellsClusters <- data.frame(Cell = names(Idents(NK_only)), 
                            CellType = as.character(Idents(NK_only)),
                            stringsAsFactors = FALSE)

## We match Progeny scores with the cell clusters.
progeny_scores_df <- inner_join(progeny_scores_df, CellsClusters)

## We summarize the Progeny scores by cell population
summarized_progeny_scores <- progeny_scores_df %>% 
  group_by(Pathway, CellType) %>%
  summarise(avg = mean(Activity), std = sd(Activity))

# Here we plot the different pathway activities for the different cell populations
## We prepare the data for the plot
summarized_progeny_scores_df <- summarized_progeny_scores %>%
  dplyr::select(-std) %>%   
  spread(Pathway, avg) %>%
  data.frame(row.names = 1, check.names = FALSE, stringsAsFactors = FALSE) 

paletteLength = 100
myColor = colorRampPalette(c("Darkblue", "white","red"))(paletteLength)

progenyBreaks = c(seq(min(summarized_progeny_scores_df), 0, 
                      length.out=ceiling(paletteLength/2) + 1),
                  seq(max(summarized_progeny_scores_df)/paletteLength, 
                      max(summarized_progeny_scores_df), 
                      length.out=floor(paletteLength/2)))
progeny_hmap = pheatmap(t(summarized_progeny_scores_df[,-1]),fontsize=14, 
                        fontsize_row = 10, 
                        color=myColor, breaks = progenyBreaks, 
                        # main = "PROGENy (500)",
                        angle_col = "45",
                        treeheight_col = 25,  border_color = NA,cluster_cols = T)







DefaultAssay(NK_only) <- "integrated"
Idents(NK_only) <- NK_only$condition
# MM_features <-(FindMarkers(NK_only, ident.1 = "MM", ident.2 = "HD",test.use = "roc",only.pos = T ))
# HD_features <-(FindMarkers(NK_only, ident.1 = "HD", ident.2 = "MM",test.use = "roc",only.pos = T))

markers <- FindAllMarkers(NK_only, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10
DoHeatmap(NK_only, features = top10$gene,) + NoLegend()
markers %>%
  group_by(cluster) %>%
  top_n(n = 20, wt = avg_log2FC) -> top20
top20hm_HD_MM <- DoHeatmap(NK_only,features = top20$gene)+ NoLegend()
Idents(NK_only) <- NK_only$seurat_clusters
top20hm <-DoHeatmap(NK_only,features = top20$gene)+ NoLegend()
top20hm_HD_MM + top20hm


markers <- FindAllMarkers(NK_only, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

p0 <- dittoDimPlot(NK_only,var = "condition",do.raster = T)
p1 <- dittoDimPlot(NK_only,var = "seurat_clusters", split.by = "condition",do.contour = T,do.label = T,do.raster = T)

p2 <- dittoBarPlot(NK_only,var = "condition",group.by = "seurat_clusters",x.labels.rotate = 45,
                   x.reorder = c(1,2,11,12,13,14,15,16,17,18,3,4,5,6,7,8,9,10))
p0 | (p1/p2)
p1|top20hm_HD_MM|p2


dittoFreqPlot(NK_only,var= "condition",sample.by = "seurat_clusters",group.by = "seurat_clusters")
# heatmap_genes <- top20$gene
# dittoHeatmap(NK_only, heatmap_genes,
#              annot.by = c("condition", "seurat_clusters"))


top20hm_HD_MM | feature_p

feature_p <- FeaturePlot(NK_only, features = c("NKG7","CD3D","IL7R","CD8A","FCGR3A","GZMK","GZMB","CX3CR1", "CD7","adt_CD16","adt_CD335-NKp46","adt_CD38"))




p2 <- dittoBarPlot(object = NK_only,group.by = "seurat_clusters",var = "condition")
# (p1|p2)/(feature_p|mon_plot)
top20hm_HD_MM|top20hm




Idents(NK_only) <- NK_only$seurat_clusters
# gs <-row.names(top_features)
# p1 <-DoHeatmap(NK_only, features = gs)
# 
# rm(list=setdiff(ls(), "NK_only"))
# gc()
# DotPlot(NK_only,features = top20$gene)
# DoHeatmap(NK_only,features = gs,angle = 45)
# dittoDotPlot(NK_only,vars = top20$gene,group.by = "condition")
# p2 <-dittoDotPlot(NK_only,vars = gs,group.by = "condition")
# dittoDotPlot(NK_only,vars = top20$gene,group.by = "seurat_clusters",split.by = "condition",assay = "RNA")
# # dittoHeatmap(NK_only,vars = top20$gene,group.by = "seurat_clusters",split.by = "condition",assay = "RNA")
# p4 <-dittoDotPlot(NK_only,vars = gs,group.by = "seurat_clusters")

markers <- FindAllMarkers(NK_only, only.pos = TRUE, test.use = "MAST", min.pct = 0.25, logfc.threshold = 0.25)
(p1+p2)/(p3+p4)
top10 <- markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
top50 <- features%>% top_n(n = 50, wt = power)
DoHeatmap(NK_only, features = gs, group.by = c("condition"), label = T)


plot(table(NK_only$seurat_clusters))
plot(prop.table(table(NK_only$seurat_clusters)))


Idents(NK_only) <- NK_only$condition

DotPlot(NK_only,cluster.idents = T, features = sort(unique(top10$gene))) + RotatedAxis()









