setwd("/home/rekren/work/")

library(Seurat)
library(SeuratDisk)
library(SeuratWrappers)
library(dplyr)
#
######## Loading untouched cohort data ######## 
pbmc.combined <-LoadH5Seurat(file = "/home/rekren/save/resources/raw_pbmc_combined.h5seurat")
######## Doing preps for RPCA data integration for batch-effect prevention, CCA integration is not favorable for big dataset ######## 
pbmc.list <- SplitObject(pbmc.combined,split.by="batch")
# normalize and identify variable features for each dataset independently
pbmc.list <- lapply(X = pbmc.list, FUN = function(x) {
   x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})
# select features that are repeatedly variable across datasets for integration run PCA on each
# dataset using these features
features <- SelectIntegrationFeatures(object.list = pbmc.list)
pbmc.list <- lapply(X = pbmc.list, FUN = function(x) {
   x <- ScaleData(x, features = features, verbose = T)
   x <- RunPCA(x, features = features, verbose = T)
})
######## Integration of different batch of data ########     
immune.anchors <- FindIntegrationAnchors(object.list = pbmc.list, anchor.features = features, reduction = "rpca")
# this command creates an 'integrated' data assay
pbmc.combined <- IntegrateData(anchorset = immune.anchors)
# Now we can run a single integrated analysis on all cells!
# Run the standard workflow for visualization and clustering
#
# specify that we will perform downstream analysis on the corrected data note that the original
# unmodified data still resides in the 'RNA' assay
DefaultAssay(pbmc.combined) <- "integrated"
pbmc.combined <- ScaleData(pbmc.combined, verbose = FALSE)
pbmc.combined <- RunPCA(pbmc.combined, verbose = FALSE)
pbmc.combined <- RunUMAP(pbmc.combined, reduction = "pca", dims = 1:30)
pbmc.combined <- FindNeighbors(pbmc.combined, reduction = "pca",dims = 1:30)
pbmc.combined <- FindClusters(pbmc.combined, resolution = 0.5)
SaveH5Seurat(pbmc.combined,verbose = T,filename = 'integrated_pbmc_cohort.h5seurat')
# Save an object to a file
saveRDS(pbmc.combined, file = "integrated_pbmc_cohort.rds")
# Restore the object
###### Preps for WNN ######
DefaultAssay(pbmc.combined) <- 'RNA'
pbmc.combined <- NormalizeData(pbmc.combined) %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA()
#
DefaultAssay(pbmc.combined) <- 'ADT'
# we will use all ADT features for dimensional reduction
# we set a dimensional reduction name to avoid overwriting 
VariableFeatures(pbmc.combined) <- rownames(pbmc.combined[["ADT"]])
pbmc.combined <- NormalizeData(pbmc.combined, normalization.method = 'CLR', margin = 2) %>%
  ScaleData() %>% RunPCA(reduction.name = 'apca')
#
# Identify multimodal neighbors. These will be stored in the neighbors slot, 
# and can be accessed using pbmc.combined[['weighted.nn']]
# The WNN graph can be accessed at pbmc.combined[["wknn"]], 
# and the SNN graph used for clustering at pbmc.combined[["wsnn"]]
# Cell-specific modality weights can be accessed at pbmc.combined$RNA.weight
pbmc.combined <- FindMultiModalNeighbors(
   pbmc.combined, reduction.list = list("pca", "apca"),
   dims.list = list(1:30, 1:18), modality.weight.name = "RNA.weight"
 )
pbmc.combined <- RunUMAP(pbmc.combined, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
pbmc.combined <- FindClusters(pbmc.combined, graph.name = "wsnn", algorithm = 3, resolution = 2, verbose = FALSE)
SaveH5Seurat(pbmc.combined,verbose = T,filename = 'integrated_wnn_cohort.h5seurat')
saveRDS(pbmc.combined, file = "wnn_integrated_pbmc_cohort.rds")

######## Reference mapping ########
rm(list=setdiff(ls(), "pbmc.combined"))
gc()

pbmc.combined <- LoadH5Seurat(file = "/home/rekren/work/integrated_wnn_cohort.h5seurat")
reference <- LoadH5Seurat(file = "/home/rekren/save/resources/reference_pbmc_multimodal.h5seurat")
if (!all(c("nCount_RNA", "nFeature_RNA") %in% c(colnames(x = pbmc.combined[[]])))) {
  calcn <- as.data.frame(x = Seurat:::CalcN(object = pbmc.combined))
  colnames(x = calcn) <- paste(
    colnames(x = calcn),
    "RNA",
    sep = '_'
  )
  pbmc.combined <- AddMetaData(
    object = pbmc.combined,
    metadata = calcn
  )
  rm(calcn)
}

# Calculate percent mitochondrial genes if the query contains genes
# matching the regular expression "^MT-"
if (any(grepl(pattern = '^MT-', x = rownames(x = pbmc.combined)))) {
  pbmc.combined <- PercentageFeatureSet(
    object = pbmc.combined,
    pattern = '^MT-',
    col.name = 'percent.mt',
    assay = "RNA"
  )
}

# Filter cells based on the thresholds for nCount_RNA and nFeature_RNA
# you set in the app
cells.use <- pbmc.combined[["nCount_RNA", drop = TRUE]] <= 7500 &
  pbmc.combined[["nCount_RNA", drop = TRUE]] >= 500 &
  pbmc.combined[["nFeature_RNA", drop = TRUE]] <= 2000 &
  pbmc.combined[["nFeature_RNA", drop = TRUE]] >= 50

# If the query contains mitochondrial genes, filter cells based on the
# thresholds for percent.mt you set in the app
if ("percent.mt" %in% c(colnames(x = pbmc.combined[[]]))) {
  cells.use <- pbmc.combined[["percent.mt", drop = TRUE]] <= 25 &
    pbmc.combined[["percent.mt", drop = TRUE]] >= 1
}

# Remove filtered cells from the query
pbmc.combined <- pbmc.combined[, cells.use]

# Preprocess with SCTransform
pbmc.combined <- SCTransform(
  object = pbmc.combined,
  assay = "RNA",
  method = 'glmGamPoi')

# Find anchors between query and reference
anchors <- FindTransferAnchors(
  reference = reference,
  query = pbmc.combined,
  normalization.method = "SCT",
  reference.reduction = "spca",
  dims = 1:50,verbose = T, recompute.residuals = F
)


pbmc.combined <- MapQuery(
  anchorset = anchors,
  query = pbmc.combined,
  reference = reference,
  refdata = list(
    celltype.l1 = "celltype.l1",
    celltype.l2 = "celltype.l2",
    predicted_ADT = "ADT"
  ),
  reference.reduction = "spca", 
  reduction.model = "wnn.umap",verbose = T
)

SaveH5Seurat(pbmc.combined,verbose = T,filename = 'mapped_integrated_pbmc_cohort.h5seurat')
saveRDS(pbmc.combined, file = "annotated_wnn_integrated_pbmc_cohort.rds")
#save.image(file = "full_envr.RData")
#quit()
