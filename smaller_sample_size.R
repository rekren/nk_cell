set_libs <- function() {
  library(Seurat);library(ggplot2);library(patchwork);library(tidyverse); library(monocle3);library(SeuratDisk);library(SeuratWrappers)
}
set_libs()
samples<- c(
  '/mnt/SERVER-CRCT-STORAGE/CRCT13/99PoleTechno/BioInformatique/CITE_seq/2021_CellRanger_Citeseq_NC01/H21b/outs/filtered_feature_bc_matrix',
  '/mnt/SERVER-CRCT-STORAGE/CRCT13/99PoleTechno/BioInformatique/CITE_seq/2021_CellRanger_Citeseq_NC01/H21m/outs/filtered_feature_bc_matrix',
  '/mnt/SERVER-CRCT-STORAGE/CRCT13/99PoleTechno/BioInformatique/CITE_seq/2021_CellRanger_Citeseq_NC01/H22b/outs/filtered_feature_bc_matrix',
  '/mnt/SERVER-CRCT-STORAGE/CRCT13/99PoleTechno/BioInformatique/CITE_seq/2021_CellRanger_Citeseq_NC01/H22m/outs/filtered_feature_bc_matrix',
  '/mnt/SERVER-CRCT-STORAGE/CRCT13/99PoleTechno/BioInformatique/CITE_seq/2021_CellRanger_Citeseq_NC01/H23b/outs/filtered_feature_bc_matrix',
  '/mnt/SERVER-CRCT-STORAGE/CRCT13/99PoleTechno/BioInformatique/CITE_seq/2021_CellRanger_Citeseq_NC01/H23m/outs/filtered_feature_bc_matrix',
  '/mnt/SERVER-CRCT-STORAGE/CRCT13/99PoleTechno/BioInformatique/CITE_seq/2021_CellRanger_Citeseq_NC01/H26b/outs/filtered_feature_bc_matrix',
  '/mnt/SERVER-CRCT-STORAGE/CRCT13/99PoleTechno/BioInformatique/CITE_seq/2021_CellRanger_Citeseq_NC01/H26m/outs/filtered_feature_bc_matrix',
  '/mnt/SERVER-CRCT-STORAGE/CRCT13/99PoleTechno/BioInformatique/CITE_seq/2021_CellRanger_Citeseq_NC01/T22497b/outs/filtered_feature_bc_matrix',
  '/mnt/SERVER-CRCT-STORAGE/CRCT13/99PoleTechno/BioInformatique/CITE_seq/2021_CellRanger_Citeseq_NC01/T22497m/outs/filtered_feature_bc_matrix',
  '/mnt/SERVER-CRCT-STORAGE/CRCT13/99PoleTechno/BioInformatique/CITE_seq/2021_CellRanger_Citeseq_NC01/T23491b/outs/filtered_feature_bc_matrix',
  '/mnt/SERVER-CRCT-STORAGE/CRCT13/99PoleTechno/BioInformatique/CITE_seq/2021_CellRanger_Citeseq_NC01/T23491m/outs/filtered_feature_bc_matrix',
  '/mnt/SERVER-CRCT-STORAGE/CRCT13/99PoleTechno/BioInformatique/CITE_seq/2021_CellRanger_Citeseq_NC01/T23819b/outs/filtered_feature_bc_matrix',
  '/mnt/SERVER-CRCT-STORAGE/CRCT13/99PoleTechno/BioInformatique/CITE_seq/2021_CellRanger_Citeseq_NC01/T23819m/outs/filtered_feature_bc_matrix',
  '/mnt/SERVER-CRCT-STORAGE/CRCT13/99PoleTechno/BioInformatique/CITE_seq/2021_CellRanger_Citeseq_NC01/T23833b/outs/filtered_feature_bc_matrix',
  '/mnt/SERVER-CRCT-STORAGE/CRCT13/99PoleTechno/BioInformatique/CITE_seq/2021_CellRanger_Citeseq_NC01/T23833m/outs/filtered_feature_bc_matrix')

name_of_samples <- basename(dirname(dirname(samples)))

i <-1
#### For the first sample, run those manually. Then, there will come `repeater`
tmp <- Seurat::Read10X(data.dir =samples[i],unique.features = T,cell.column = 1,gene.column = 2);

tmp_seurat <- CreateSeuratObject(counts = tmp[["Gene Expression"]], min.cells = 3, min.features = 200, project = name_of_samples[i])
tmp_seurat <- NormalizeData(tmp_seurat)
###### Cleaning TotalSeqB tag from the prospective ADT data
rownames(tmp[["Antibody Capture"]])<- stringr::str_replace_all(rownames(tmp[["Antibody Capture"]]),
                                                               setNames("", "_*TotalSeqB"));
###### Adding ADT modality to Seurat object 
tmp_seurat[["ADT"]] <- CreateAssayObject(tmp[["Antibody Capture"]][, colnames(x = tmp_seurat)])
tmp_seurat <- NormalizeData(tmp_seurat, assay = "ADT", normalization.method = "CLR",
                            project = name_of_samples[i])
###### Embedding additional infos
tmp_seurat$condition = ifelse(grepl(pattern="H",x=name_of_samples[i]),"HD","MM")
tmp_seurat$source = ifelse(grepl(pattern="m",x=name_of_samples[i]),"BoneMarrow","Blood")
rm(tmp)
####   Only required for for 1st iteration, manually run for the first one up to here    
pbmc.combined <- tmp_seurat
i <- i+1

repeat {
  
  # Run the inner loop for once manually then  repeat it with repeater
  ###### Reading one by one, each element from the list of samples into Seurat
  tmp <- Seurat::Read10X(data.dir =samples[i],unique.features = T,cell.column = 1,gene.column = 2);
  
  tmp_seurat <- CreateSeuratObject(counts = tmp[["Gene Expression"]], min.cells = 3, min.features = 200, project = name_of_samples[i])
  tmp_seurat <- NormalizeData(tmp_seurat)
  ###### Cleaning TotalSeqB tag from the prospective ADT data
  rownames(tmp[["Antibody Capture"]])<- stringr::str_replace_all(rownames(tmp[["Antibody Capture"]]),
                                                                 setNames("", "_*TotalSeqB"));
  ###### Adding ADT modality to Seurat object 
  tmp_seurat[["ADT"]] <- CreateAssayObject(tmp[["Antibody Capture"]][, colnames(x = tmp_seurat)])
  tmp_seurat <- NormalizeData(tmp_seurat, assay = "ADT", normalization.method = "CLR",
                              project = name_of_samples[i])
  ###### Embedding additional infos
  tmp_seurat$condition = ifelse(grepl(pattern="H",x=name_of_samples[i]),"HD","MM")
  tmp_seurat$source = ifelse(grepl(pattern="m",x=name_of_samples[i]),"BoneMarrow","Blood")
  rm(tmp)
  ####   Only required for for 1st iteration, manually run for the first one up to here    
  # pbmc.combined <- tmp_seurat
  
  pbmc.combined <- merge(pbmc.combined, y = tmp_seurat)
  
  #Manually increase for the first iteration
  i <- i+1
  if(i > length(samples)) {
    break
  } 
} 

table(pbmc.combined$orig.ident)        ## Each sample's cell numbers
sum(table(pbmc.combined$orig.ident))   ## 244511 Cell
table(pbmc.combined$source)            ## Blood  114379  BoneMarrow 130132 
table(pbmc.combined$condition)         ## 105850 138661 

pbmc.combined[["percent.mt"]] <- Seurat::PercentageFeatureSet(pbmc.combined, pattern = "^MT-")

pbmc.combined <- subset(pbmc.combined, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 25)

saveRDS(object = pbmc.combined,file = 'untouched_pbmc_combined.RDS')
SaveH5Seurat(pbmc.combined,verbose = T,filename = 'untouched_pbmc_combined.h5seurat')

# reference <- LoadH5Seurat('pbmc_multimodal.h5seurat')

# >>>>> burada kaldın dün >>>>

table(pbmc.combined$orig.ident)       
sum(table(pbmc.combined$orig.ident))  
table(pbmc.combined$source)           
table(pbmc.combined$condition)        


pbmc.combined <- Seurat::NormalizeData(pbmc.combined)

# Normally for smaller cell numbers do Canonical-Correlation Analysis (CCA) 
# su.list <- Seurat::SplitObject(pbmc.combined,split.by="Origin")
# su.anchors <- Seurat::FindIntegrationAnchors(object.list=su.list,dims=1:15)
# rm(list=setdiff(ls(), "su.anchors"))

saveRDS()
gc()

# For big datasets (i.e. several hundred thousand single-cells)
# follow Reciprocal PCA (RPCA) data integration
rm(list=setdiff(ls(), "pbmc.combined"))

pbmc.list <- SplitObject(pbmc.combined,split.by="orig.ident")
pbmc.list <- lapply(X = pbmc.list, FUN = function(x) {
  x <- NormalizeData(x, verbose = T)
  x <- FindVariableFeatures(x, verbose = T)
})

features <- SelectIntegrationFeatures(object.list = pbmc.list)
pbmc.list <- lapply(X = pbmc.list, FUN = function(x) {
  x <- ScaleData(x, features = features, verbose = T)
  x <- RunPCA(x, features = features, verbose = T)
})
anchors <- FindIntegrationAnchors(object.list = pbmc.list, reduction = "rpca", 
                                  dims = 1:50)
pbmc.combined <- Seurat::IntegrateData(anchorset = anchors, dims = 1:15)

Seurat::DefaultAssay(pbmc.combined) <- "integrated"


# Run the standard workflow for visualization and clustering
pbmc.combined <- Seurat::ScaleData(pbmc.combined, verbose = T)
pbmc.combined <- Seurat::RunPCA(pbmc.combined,assay ="integrated",  npcs=15)


pbmc.combined <- Seurat::FindNeighbors(pbmc.combined,assay = "integrated", dims = 1:10)
pbmc.combined <- Seurat::FindClusters(pbmc.combined, resolution = 0.5)


pbmc.combined <- Seurat::RunUMAP(pbmc.combined, dims = 1:10)
# saveRDS(pbmc.combined,"untouched_cohort.RDS")


