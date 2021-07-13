library(Seurat);library(SeuratDisk);library(dplyr);library(tidyverse); library(ggplot2); library(reshape2);library(monocle3);library(clustree);library(SeuratWrappers);library(CIPR);library(escape);library(dittoSeq);

NK_only <- LoadH5Seurat("NK_40.h5seurat")
#####################
DefaultAssay(NK_only) <- "integrated"
NK_only <- FindVariableFeatures(NK_only,assay = "integrated")
NK_only <- ScaleData(NK_only)
NK_only <- RunPCA(NK_only)
NK_only <- FindNeighbors(NK_only,force.recalc = T)
NK_only <- JackStraw(NK_only, assay = "integrated",num.replicate = 100)
NK_only <- ScoreJackStraw(NK_only,dims = 1:20,reduction = "pca")
JackStrawPlot(NK_only)

NK_only <- FindClusters(NK_only,algorithm = 3)
NK_only <- RunUMAP(NK_only,dims = 1:10)
Idents(NK_only) <- NK_only$seurat_clusters
# clustree(NK_only,assay = "integrated")

clustree_overlay(NK_only,red_dim = "umap", x_value = "umap1", y_value = "umap2",
                 use_colour = "points", alt_colour = "blue")

clustree_overlay(NK_only,red_dim = "tsne", x_value = "tsne1", y_value = "tsne2",
                 use_colour = "points", alt_colour = "blue")
#### Create a Cell DataSet (CDS) object from Seurat object
# for Monocle3 running on the already calculated UMAPs
#### Create a Monocle CDS Object
# Project PC dimensions to whole data set
Reductions(NK_only)
my.so <- ProjectDim(NK_only,assay = "integrated", reduction = "umap")
# Create an expression matrix
expression_matrix <- my.so@assays$integrated@data
# Get cell metadata
cell_metadata <- my.so@meta.data
if (all.equal(colnames(expression_matrix), rownames(cell_metadata))) {
  print(sprintf("Cell identifiers match"))
} else {
  print(sprintf("Cell identifier mismatch - %i cells in expression matrix, %i cells in metadata",
                ncol(expression_matrix), nrow(cell_metadata)))
  print("If the counts are equal, sort differences will throw this error")
}
# get gene annotations
gene_annotation <- data.frame(gene_short_name = rownames(my.so@assays$integrated), row.names = rownames(my.so@assays$integrated))
if (all.equal(rownames(expression_matrix), rownames(gene_annotation))) {
  print(sprintf("Gene identifiers all match"))
} else {
  print(sprintf("Gene identifier mismatch - %i genes in expression matrix, %i gene in gene annotation",
                nrow(expression_matrix), nrow(gene_annotation)))
  print("If the counts are equal, sort differences will throw this error")
}
# Seurat-derived CDS
my.cds <- new_cell_data_set(expression_matrix,
                            cell_metadata = cell_metadata,
                            gene_metadata = gene_annotation)
# Transfer Seurat embeddings
# Note that these may be calculated on the Integrated object, not the counts
#   and thus will involve fewer genes
reducedDim(my.cds, type = "PCA") <- my.so@reductions$pca@cell.embeddings 
my.cds@preprocess_aux$prop_var_expl <- my.so@reductions$pca@stdev
plot_pc_variance_explained(my.cds)
# Transfer Seurat UMAP embeddings
my.cds@int_colData@listData$reducedDims$UMAP <- my.so@reductions$umap@cell.embeddings
plot_cells(my.cds,label_groups_by_cluster = T,group_label_size = 5)
# Copy cluster info from Seurat
my.cds@clusters$UMAP_so$clusters <- my.so@meta.data$seurat_clusters
my.cds <- cluster_cells(my.cds, reduction_method = "UMAP")
# Fix from https://gitmemory.com/cole-trapnell-lab
rownames(my.cds@principal_graph_aux$UMAP$dp_mst) <- NULL
colnames(my.cds@int_colData@listData$reducedDims$UMAP) <- NULL
my.cds <- learn_graph(my.cds, use_partition = TRUE,learn_graph_control=list("geodesic_distance_ratio"=1/2))
my.cds <- order_cells(my.cds)
mon_plot <- plot_cells(
  cds = my.cds,
  color_cells_by = "pseudotime",
  show_trajectory_graph = T
)

FeaturePlot(NK_only,split.by = "condition",label = T,features = c("NCAM1","PRF1","CD7","GNLY","CD3D","SDC1"))
p1 <- dittoDimPlot(object = NK_only,var = "seurat_clusters")
p2 <- dittoDimPlot(object = NK_only,var = "condition")
p3 <- dittoDimHex(object = NK_only)
p4 <- dittoBarPlot(object = NK_only,group.by = "seurat_clusters",var = "condition")
p5 <- dittoBarPlot(object = NK_only,group.by = "seurat_clusters",var = "source")
# dittoHeatmap(NK_only ,genes = top10$gene,order.by = "seurat_clusters",
#              annot.by = c("condition","seurat_clusters"), 
#              fontsize = 8, 
#              cluster_cols = F,)


dittoFreqPlot(object = NK_only,group.by = "condition",var = "seurat_clusters")


Idents(NK_only) <- NK_only$seurat_clusters
DimPlot(NK_only)




##### Looking for Cell-cycle markers and their possible effects on PCs
# ##### Not included into the analysis, just personally wanted to see outcome
# s.genes <- (cc.genes.updated.2019$s.genes)
# g2m.genes <- (cc.genes.updated.2019$g2m.genes)
# NK_only <- CellCycleScoring(object = NK_only, s.features = s.genes,g2m.features = g2m.genes,set.ident = T)
# NK_only<- RunPCA(NK_only, features = c(s.genes, g2m.genes))
# DimPlot(NK_only,reduction = "pca")
# 
# NK_only <- ScaleData(NK_only, vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(NK_only))
# # Now, a PCA on the variable genes no longer returns components associated with cell cycle
# NK_only <- RunPCA(NK_only, features = VariableFeatures(NK_only), nfeatures.print = 10)
# # When running a PCA on only cell cycle genes, cells no longer separate by cell-cycle phase
# NK_only <- RunPCA(NK_only, features = c(s.genes, g2m.genes))
# DimPlot(NK_only,reduction = "pca")
# 
# ### Alternate Cell-cycle regressing
# NK_only$CC.Difference <- NK_only$S.Score - NK_only$G2M.Score
# NK_only <- ScaleData(NK_only, vars.to.regress = "CC.Difference", features = rownames(NK_only))
# 
# # cell cycle effects strongly mitigated in PCA
# NK_only <- RunPCA(NK_only, features = VariableFeatures(NK_only), nfeatures.print = 10)
# 
# # when running a PCA on cell cycle genes, actively proliferating cells remain distinct from G1
# # cells however, within actively proliferating cells, G2M and S phase cells group together
# NK_only <- RunPCA(NK_only, features = c(s.genes, g2m.genes))
# DimPlot(NK_only)
# DimPlot(NK_only,reduction = "pca")

#Splitting the data into HD and MM samples 
Idents(NK_only) <- NK_only$condition
NK_MM <-subset(x = NK_only, idents = "MM")
NK_HD <-subset(x = NK_only, idents = "HD")

############################################################################
Reductions(NK_HD)
my.so <- ProjectDim(NK_HD,assay = "integrated", reduction = "umap")
# Create an expression matrix
expression_matrix <- my.so@assays$integrated@data
# Get cell metadata
cell_metadata <- my.so@meta.data
if (all.equal(colnames(expression_matrix), rownames(cell_metadata))) {
  print(sprintf("Cell identifiers match"))
} else {
  print(sprintf("Cell identifier mismatch - %i cells in expression matrix, %i cells in metadata",
                ncol(expression_matrix), nrow(cell_metadata)))
  print("If the counts are equal, sort differences will throw this error")
}
# get gene annotations
gene_annotation <- data.frame(gene_short_name = rownames(my.so@assays$integrated), row.names = rownames(my.so@assays$integrated))
if (all.equal(rownames(expression_matrix), rownames(gene_annotation))) {
  print(sprintf("Gene identifiers all match"))
} else {
  print(sprintf("Gene identifier mismatch - %i genes in expression matrix, %i gene in gene annotation",
                nrow(expression_matrix), nrow(gene_annotation)))
  print("If the counts are equal, sort differences will throw this error")
}
# Seurat-derived CDS
my.cds <- new_cell_data_set(expression_matrix,
                            cell_metadata = cell_metadata,
                            gene_metadata = gene_annotation)
# Transfer Seurat embeddings
# Note that these may be calculated on the Integrated object, not the counts
#   and thus will involve fewer genes
reducedDim(my.cds, type = "PCA") <- my.so@reductions$pca@cell.embeddings 
my.cds@preprocess_aux$prop_var_expl <- my.so@reductions$pca@stdev
plot_pc_variance_explained(my.cds)
# Transfer Seurat UMAP embeddings
my.cds@int_colData@listData$reducedDims$UMAP <- my.so@reductions$umap@cell.embeddings
plot_cells(my.cds,label_groups_by_cluster = T,group_label_size = 5)
# Copy cluster info from Seurat
my.cds@clusters$UMAP_so$clusters <- my.so@meta.data$seurat_clusters
my.cds <- cluster_cells(my.cds, reduction_method = "UMAP")
# Fix from https://gitmemory.com/cole-trapnell-lab
rownames(my.cds@principal_graph_aux$UMAP$dp_mst) <- NULL
colnames(my.cds@int_colData@listData$reducedDims$UMAP) <- NULL
my.cds <- learn_graph(my.cds, use_partition = TRUE,learn_graph_control=list("geodesic_distance_ratio"=1/2))
my.cds <- order_cells(my.cds)
mon_plot <- plot_cells(
  cds = my.cds,
  color_cells_by = "pseudotime",
  show_trajectory_graph = T
)





Reductions(NK_MM)
my.so <- ProjectDim(NK_MM,assay = "integrated", reduction = "umap")
# Create an expression matrix
expression_matrix <- my.so@assays$integrated@data
# Get cell metadata
cell_metadata <- my.so@meta.data
if (all.equal(colnames(expression_matrix), rownames(cell_metadata))) {
  print(sprintf("Cell identifiers match"))
} else {
  print(sprintf("Cell identifier mismatch - %i cells in expression matrix, %i cells in metadata",
                ncol(expression_matrix), nrow(cell_metadata)))
  print("If the counts are equal, sort differences will throw this error")
}
# get gene annotations
gene_annotation <- data.frame(gene_short_name = rownames(my.so@assays$integrated), row.names = rownames(my.so@assays$integrated))
if (all.equal(rownames(expression_matrix), rownames(gene_annotation))) {
  print(sprintf("Gene identifiers all match"))
} else {
  print(sprintf("Gene identifier mismatch - %i genes in expression matrix, %i gene in gene annotation",
                nrow(expression_matrix), nrow(gene_annotation)))
  print("If the counts are equal, sort differences will throw this error")
}
# Seurat-derived CDS
my.cds <- new_cell_data_set(expression_matrix,
                            cell_metadata = cell_metadata,
                            gene_metadata = gene_annotation)
# Transfer Seurat embeddings
# Note that these may be calculated on the Integrated object, not the counts
#   and thus will involve fewer genes
reducedDim(my.cds, type = "PCA") <- my.so@reductions$pca@cell.embeddings 
my.cds@preprocess_aux$prop_var_expl <- my.so@reductions$pca@stdev
plot_pc_variance_explained(my.cds)
# Transfer Seurat UMAP embeddings
my.cds@int_colData@listData$reducedDims$UMAP <- my.so@reductions$umap@cell.embeddings
plot_cells(my.cds,label_groups_by_cluster = T,group_label_size = 5)
# Copy cluster info from Seurat
my.cds@clusters$UMAP_so$clusters <- my.so@meta.data$seurat_clusters
my.cds <- cluster_cells(my.cds, reduction_method = "UMAP")
# Fix from https://gitmemory.com/cole-trapnell-lab
rownames(my.cds@principal_graph_aux$UMAP$dp_mst) <- NULL
colnames(my.cds@int_colData@listData$reducedDims$UMAP) <- NULL
my.cds <- learn_graph(my.cds, use_partition = TRUE,learn_graph_control=list("geodesic_distance_ratio"=1/2))
my.cds <- order_cells(my.cds)
mon_plot_2 <- plot_cells(
  cds = my.cds,
  color_cells_by = "pseudotime",
  show_trajectory_graph = T
)
