###################################################
# library(Seurat);library(ggplot2);library(patchwork);
# # Load in the RNA UMI matrix
# 
# # Note that this dataset also contains ~5% of mouse cells, which we can use as negative controls
# # for the protein measurements. For this reason, the gene expression matrix has HUMAN_ or MOUSE_
# # appended to the beginning of each gene.
# cbmc.rna <- as.sparse(read.csv(file = "/home/ruchan/Desktop/work/demo/GSE100866_CBMC_8K_13AB_10X-RNA_umi.csv.gz", sep = ",", 
#                                header = TRUE, row.names = 1))
# # To make life a bit easier going forward, we're going to discard all but the top 100 most
# # highly expressed mouse genes, and remove the 'HUMAN_' from the CITE-seq prefix
# cbmc.rna <- CollapseSpeciesExpressionMatrix(cbmc.rna)
# # Load in the ADT UMI matrix
# cbmc.adt <- as.sparse(read.csv(file = "/home/ruchan/Desktop/work/demo/GSE100866_CBMC_8K_13AB_10X-ADT_umi.csv.gz", sep = ",", 
#                                header = TRUE, row.names = 1))
# # When adding multimodal data to Seurat, it's okay to have duplicate feature names. Each set of
# # modal data (eg. RNA, ADT, etc.) is stored in its own Assay object.  One of these Assay objects
# # is called the 'default assay', meaning it's used for all analyses and visualization.  To pull
# # data from an assay that isn't the default, you can specify a key that's linked to an assay for
# # feature pulling.  To see all keys for all objects, use the Key function.  Lastly, we observed
# # poor enrichments for CCR5, CCR7, and CD10 - and therefore remove them from the matrix
# # (optional)
# cbmc.adt <- cbmc.adt[setdiff(rownames(x = cbmc.adt), c("CCR5", "CCR7", "CD10")), ]
# 
# cbmc <- CreateSeuratObject(counts = cbmc.rna)
# # standard log-normalization
# cbmc <- NormalizeData(cbmc)
# # choose ~1k variable features
# cbmc <- FindVariableFeatures(cbmc)
# # standard scaling (no regression)
# cbmc <- ScaleData(cbmc)
# # Run PCA, select 13 PCs for tSNE visualization and graph-based clustering
# cbmc <- RunPCA(cbmc, verbose = FALSE)
# ElbowPlot(cbmc, ndims = 50)
# cbmc <- FindNeighbors(cbmc, dims = 1:25)
# cbmc <- FindClusters(cbmc, resolution = 0.8)
# cbmc <- RunTSNE(cbmc, dims = 1:25, method = "FIt-SNE")
# 
# # Find the markers that define each cluster, and use these to annotate the clusters, we use
# # max.cells.per.ident to speed up the process
# cbmc.rna.markers <- FindAllMarkers(cbmc, max.cells.per.ident = 100, min.diff.pct = 0.3, only.pos = TRUE)
# # Note, for simplicity we are merging two CD14+ Monocyte clusters (that differ in expression of
# # HLA-DR genes) and NK clusters (that differ in cell cycle stage)
# DimPlot(cbmc, label = TRUE) + NoLegend()
# new.cluster.ids <- c("Memory CD4 T", "CD14+ Mono", "Naive CD4 T", "NK", "CD14+ Mono", "Mouse", "B", 
#                      "CD8 T", "CD16+ Mono", "T/Mono doublets", "NK", "CD34+", "Multiplets", "Mouse", "Eryth", "Mk", 
#                      "Mouse", "DC", "pDCs")
# names(new.cluster.ids) <- levels(cbmc)
# cbmc <- RenameIdents(cbmc, new.cluster.ids)
# DimPlot(cbmc, label = TRUE) + NoLegend()
# # We will define an ADT assay, and store raw counts for it
# 
# # If you are interested in how these data are internally stored, you can check out the Assay
# # class, which is defined in objects.R; note that all single-cell expression data, including RNA
# # data, are still stored in Assay objects, and can also be accessed using GetAssayData
# cbmc[["ADT"]] <- CreateAssayObject(counts = cbmc.adt)
# 
# # Now we can repeat the preprocessing (normalization and scaling) steps that we typically run
# # with RNA, but modifying the 'assay' argument.  For CITE-seq data, we do not recommend typical
# # LogNormalization. Instead, we use a centered log-ratio (CLR) normalization, computed
# # independently for each feature.  This is a slightly improved procedure from the original
# # publication, and we will release more advanced versions of CITE-seq normalizations soon.
# cbmc <- NormalizeData(cbmc, assay = "ADT", normalization.method = "CLR")
# cbmc <- ScaleData(cbmc, assay = "ADT")
# # in this plot, protein (ADT) levels are on top, and RNA levels are on the bottom
# FeaturePlot(cbmc, features = c("adt_CD3", "adt_CD11c", "adt_CD8", "adt_CD16", "CD3E", "ITGAX", "CD8A", 
#                                "FCGR3A"), min.cutoff = "q05", max.cutoff = "q95", ncol = 4)
# RidgePlot(cbmc, features = c("adt_CD3", "adt_CD11c", "adt_CD8", "adt_CD16"), ncol = 2)
# # Draw ADT scatter plots (like biaxial plots for FACS). Note that you can even 'gate' cells if
# # desired by using HoverLocator and FeatureLocator
# FeatureScatter(cbmc, feature1 = "adt_CD19", feature2 = "adt_CD3")
# # view relationship between protein and RNA
# FeatureScatter(cbmc, feature1 = "adt_CD3", feature2 = "CD3E")
# # Let's plot CD4 vs CD8 levels in T cells
# tcells <- subset(cbmc, idents = c("Naive CD4 T", "Memory CD4 T", "CD8 T"))
# FeatureScatter(tcells, feature1 = "adt_CD4", feature2 = "adt_CD8")
# # # Let's look at the raw (non-normalized) ADT counts. You can see the values are quite high,
# # particularly in comparison to RNA values. This is due to the significantly higher protein copy
# # number in cells, which significantly reduces 'drop-out' in ADT data
# FeatureScatter(tcells, feature1 = "adt_CD4", feature2 = "adt_CD8", slot = "counts")
# 
# saveRDS(cbmc, file = "../output/cbmc.rds")
# 
# 
# # Downsample the clusters to a maximum of 300 cells each (makes the heatmap easier to see for
# # small clusters)
# cbmc.small <- subset(cbmc, downsample = 300)
# 
# # Find protein markers for all clusters, and draw a heatmap
# adt.markers <- FindAllMarkers(cbmc.small, assay = "ADT", only.pos = TRUE)
# DoHeatmap(cbmc.small, features = unique(adt.markers$gene), assay = "ADT", angle = 90) + NoLegend()
# # You can see that our unknown cells co-express both myeloid and lymphoid markers (true at the
# # RNA level as well). They are likely cell clumps (multiplets) that should be discarded. We'll
# # remove the mouse cells now as well
# cbmc <- subset(cbmc, idents = c("Multiplets", "Mouse"), invert = TRUE)
# 
# ###
#  ###Cluster directly on protein levels
# ###
# # Because we're going to be working with the ADT data extensively, we're going to switch the
# # default assay to the 'ADT' assay.  This will cause all functions to use ADT data by default,
# # rather than requiring us to specify it each time
# DefaultAssay(cbmc) <- "ADT"
# cbmc <- RunPCA(cbmc, features = rownames(cbmc), reduction.name = "pca_adt", reduction.key = "pca_adt_", 
#                verbose = FALSE)
# DimPlot(cbmc, reduction = "pca_adt")
# 
# 
# # Before we recluster the data on ADT levels, we'll stash the RNA cluster IDs for later
# cbmc[["rnaClusterID"]] <- Idents(cbmc)
# 
# # Now, we rerun tSNE using the PCA only on ADT (protein) levels.
# cbmc <- RunTSNE(cbmc, dims = 1:10, reduction = "pca_adt", reduction.key = "adtTSNE_", reduction.name = "tsne_adt")
# cbmc <- FindNeighbors(cbmc, features = rownames(cbmc), dims = NULL)
# cbmc <- FindClusters(cbmc, resolution = 0.2, graph.name = "ADT_snn")
# 
# # We can compare the RNA and protein clustering, and use this to annotate the protein clustering
# # (we could also of course use FindMarkers)
# clustering.table <- table(Idents(cbmc), cbmc$rnaClusterID)
# clustering.table
# 
# 
# 
# new.cluster.ids <- c("CD4 T", "CD14+ Mono", "NK", "B", "CD8 T", "NK", "CD34+", "T/Mono doublets", 
#                      "CD16+ Mono", "pDCs", "B")
# names(new.cluster.ids) <- levels(cbmc)
# cbmc <- RenameIdents(cbmc, new.cluster.ids)
# 
# tsne_rnaClusters <- DimPlot(cbmc, reduction = "tsne_adt", group.by = "rnaClusterID") + NoLegend()
# tsne_rnaClusters <- tsne_rnaClusters + ggtitle("Clustering based on scRNA-seq") + theme(plot.title = element_text(hjust = 0.5))
# tsne_rnaClusters <- LabelClusters(plot = tsne_rnaClusters, id = "rnaClusterID", size = 4)
# 
# tsne_adtClusters <- DimPlot(cbmc, reduction = "tsne_adt", pt.size = 0.5) + NoLegend()
# tsne_adtClusters <- tsne_adtClusters + ggtitle("Clustering based on ADT signal") + theme(plot.title = element_text(hjust = 0.5))
# tsne_adtClusters <- LabelClusters(plot = tsne_adtClusters, id = "ident", size = 4)
# 
# # Note: for this comparison, both the RNA and protein clustering are visualized on a tSNE
# # generated using the ADT distance matrix.
# wrap_plots(list(tsne_rnaClusters, tsne_adtClusters), ncol = 2)
# 
# ###The ADT-based clustering yields similar results, but with a few differences
# # Clustering is improved for CD4/CD8 T cell populations, based on the robust ADT data for CD4, CD8, CD14, and CD45RA
# # However, some clusters for which the ADT data does not contain good distinguishing protein markers (i.e. Mk/Ery/DC) lose separation
# # You can verify this using FindMarkers at the RNA level, as well
# # 
# # 
# 
# tcells <- subset(cbmc, idents = c("CD4 T", "CD8 T"))
# FeatureScatter(tcells, feature1 = "CD4", feature2 = "CD8")
# 
# RidgePlot(cbmc, features = c("adt_CD11c", "adt_CD8", "adt_CD16", "adt_CD4", "adt_CD19", "adt_CD14"), 
#           ncol = 2)
# 
# 
# saveRDS(cbmc, file = "../output/cbmc_multimodal.rds") 
# 
# 
# 
# library(Seurat);library(ggplot2);library(patchwork);library(tidyverse);
# nk <- readRDS("/home/ruchan/Desktop/work/RobjectsAndCo/NKonly.RDS")
# nk <- readRDS("NKonly.RDS")
# # nk <- readRDS("/home/ruchan/Desktop/work/RobjectsAndCo/NKonly.RDS")
# nk 
# # nk <-FindVariableFeatures(nk)
# nk
# ElbowPlot(nk,ndims = 50,reduction = "pca")
# DimPlot(nk,reduction = "umap", label = TRUE,pt.size = 2) + NoLegend()
# new.cluster.ids <- c("4","4","4","4","1","2","4","1","3","4","3","4","2","4","4","_15_","1")
# names(new.cluster.ids) <- levels(nk)
# nk <- RenameIdents(nk, new.cluster.ids)
# nk.rna.markers <- FindAllMarkers(nk, only.pos = TRUE,test.use = "MAST",min.pct = 0.25, logfc.threshold = 0.25)
# nk.small <- subset(nk, downsample = 300)
# 
# DoHeatmap(nk.small,features = unique(nk.rna.markers$gene),assay = 'SCT')
# top10 <- nk.rna.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
# DoHeatmap(nk, features = unique(top10$gene),assay = 'SCT')
# 
# 
# nk <- NormalizeData(nk, assay = "ADT", normalization.method = "CLR")
# nk <- ScaleData(nk, assay = "ADT")

# Weighted nearest neighbor analysis DEMO#
# devtools::install_github('satijalab/seurat-data')
library(Seurat)
library(SeuratData)
library(cowplot)
library(dplyr)
InstallData("bmcite")
bm <- LoadData(ds = "bmcite")
# We first perform pre-processing and dimensional reduction on both assays independently. 
# We use standard normalization, but you can also use SCTransform or any alternative method.
DefaultAssay(bm) <- 'RNA'
bm <- NormalizeData(bm) %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA()

DefaultAssay(bm) <- 'ADT'
# we will use all ADT features for dimensional reduction
# we set a dimensional reduction name to avoid overwriting the 
VariableFeatures(bm) <- rownames(bm[["ADT"]])
bm <- NormalizeData(bm, normalization.method = 'CLR', margin = 2) %>% 
  ScaleData() %>% RunPCA(reduction.name = 'apca')
# For each cell, we calculate its closest neighbors in the dataset based on a weighted 
# combination of RNA and protein similarities. The cell-specific modality weights and 
# multimodal neighbors are calculated in a single function, which takes ~2 minutes to 
# run on this dataset. We specify the dimensionality of each modality (similar to specifying 
# the number of PCs to include in scRNA-seq clustering), but you can vary these settings to see 
# that small changes have minimal effect on the overall results.
# 
# Identify multimodal neighbors. These will be stored in the neighbors slot,
# and can be accessed using bm[['weighted.nn']]
# The WNN graph can be accessed at bm[["wknn"]],
# and the SNN graph used for clustering at bm[["wsnn"]]
# Cell-specific modality weights can be accessed at bm$RNA.weight
bm <- FindMultiModalNeighbors(
  bm, reduction.list = list("pca", "apca"), 
  dims.list = list(1:30, 1:18), modality.weight.name = "RNA.weight"
)

# We can now use these results for downstream analysis, such as visualization 
# and clustering. For example, we can create a UMAP visualization of the data based 
# on a weighted combination of RNA and protein data We can also perform graph-based 
# clustering and visualize these results on the UMAP, alongside a set of cell annotations.

bm <- RunUMAP(bm, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
bm <- FindClusters(bm, graph.name = "wsnn", algorithm = 3, resolution = 2, verbose = FALSE)


p1 <- DimPlot(bm, reduction = 'wnn.umap', label = TRUE, repel = TRUE, label.size = 2.5) + NoLegend()
p2 <- DimPlot(bm, reduction = 'wnn.umap', group.by = 'celltype.l2', label = TRUE, repel = TRUE, label.size = 2.5) + NoLegend()
p1 + p2

# We can also compute UMAP visualization based on only the RNA and protein data and compare. 
# We find that the RNA analysis is more informative than the ADT analysis in identifying 
# progenitor states (the ADT panel contains markers for differentiated cells), 
# while the converse is true of T cell states (where the ADT analysis outperforms RNA).

bm <- RunUMAP(bm, reduction = 'pca', dims = 1:30, assay = 'RNA', 
              reduction.name = 'rna.umap', reduction.key = 'rnaUMAP_')
bm <- RunUMAP(bm, reduction = 'apca', dims = 1:18, assay = 'ADT', 
              reduction.name = 'adt.umap', reduction.key = 'adtUMAP_')
p3 <- DimPlot(bm, reduction = 'rna.umap', group.by = 'celltype.l2', label = TRUE, 
              repel = TRUE, label.size = 2.5) + NoLegend()
p4 <- DimPlot(bm, reduction = 'adt.umap', group.by = 'celltype.l2', label = TRUE, 
              repel = TRUE, label.size = 2.5) + NoLegend()
p3 + p4

# We can visualize the expression of canonical marker genes and proteins on the multimodal UMAP, 
# which can assist in verifying the provided annotations:

p5 <- FeaturePlot(bm, features = c("adt_CD45RA","adt_CD16","adt_CD161"),
                    reduction = 'wnn.umap', max.cutoff = 2, 
                    cols = c("lightgrey","darkgreen"), ncol = 3)
p6 <- FeaturePlot(bm, features = c("rna_TRDC","rna_MPO","rna_AVP"), 
                  reduction = 'wnn.umap', max.cutoff = 3, ncol = 3)
p5 / p6

# Finally, we can visualize the modality weights that were learned for each cell. 
# Each of the populations with the highest RNA weights represent progenitor cells,
# while the populations with the highest protein weights represent T cells. This is in 
# line with our biological expectations, as the antibody panel does not contain markers 
# that can distinguish between different progenitor populations.

VlnPlot(bm, features = "RNA.weight", group.by = 'celltype.l2', sort = TRUE, pt.size = 0.1) +
  NoLegend()
