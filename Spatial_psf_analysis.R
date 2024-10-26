library(Seurat)
library(SeuratData)
library(ggplot2)
library(patchwork)
library(dplyr)
library(biomaRt)
library(psf)

spatial_psf_analysis <- function(spatial_obj, pathway_collection, gene_symbol_to_entrez, nthreads = 1, return_only_shiny_vars = TRUE) {
  
  ### Exp data Normalization
  cat("Exp SCT transform (normalization)\n")
  spatial_obj <- SCTransform(spatial_obj, assay = "Spatial", verbose = FALSE)
  
  ### Exp clustering
  cat("Exp clustering\n")
  spatial_obj <- RunPCA(spatial_obj, assay = "SCT", verbose = FALSE)
  spatial_obj <- FindNeighbors(spatial_obj, reduction = "pca", dims = 1:30)
  spatial_obj <- FindClusters(spatial_obj, verbose = FALSE, resolution = 0.8)
  spatial_obj <- RunUMAP(spatial_obj, reduction = "pca", dims = 1:30)
  
  
  ### Exp fold change calculation
  spatial_sct_matrix <- as.matrix(Seurat::Assays(spatial_obj, slot = "SCT")@counts)
  spatial_sct_matrix_fc <- (spatial_sct_matrix + 1)/rowMeans(spatial_sct_matrix + 1)
  
  rownames(spatial_sct_matrix_fc) <- gene_symbol_to_entrez[rownames(spatial_sct_matrix_fc)]
  spatial_sct_matrix_fc <- spatial_sct_matrix_fc[which(!is.na(rownames(spatial_sct_matrix_fc))),]
  
  ### Pathway signal flow analysis
  cat("Pathway activity calculation\n")
  spatial_psf <- run_psf(entrez.fc = spatial_sct_matrix_fc, kegg.collection = pathway_collection, 
                         calculate.significance = F, ncores = nthreads)
  
  
  spatial_psf_mat <- Reduce(rbind, 
                            lapply(names(spatial_psf), function(x) {
                              print(x)
                              pathway <- spatial_psf[[x]]
                              
                              activity_mat <- pathway$psf_activities[pathway$sink.nodes,,drop = F]
                              
                              rownames(activity_mat) <- paste0(rownames(activity_mat), "; ", x)
                              
                              activity_mat
                            })   
  )
  
  
  spatial_psf_obj <- CreateSeuratObject(counts = spatial_psf_mat, meta.data = spatial_obj@meta.data, assay = "Spatial")
  
  cat("Patwahy activity clustering\n")
  ### Patwhay activity data normalization
  # spatial_psf_obj <- NormalizeData(spatial_psf_obj)
  
  ### Pathway activity clustering
  spatial_psf_obj <- FindVariableFeatures(spatial_psf_obj)
  spatial_psf_obj <- ScaleData(spatial_psf_obj)
  spatial_psf_obj <- RunPCA(spatial_psf_obj, verbose = FALSE, approx=TRUE)
  spatial_psf_obj <- FindNeighbors(spatial_psf_obj, dims = 1:30)
  spatial_psf_obj <- FindClusters(spatial_psf_obj, resolution = 0.8, verbose = FALSE)
  spatial_psf_obj <- RunUMAP(spatial_psf_obj, dims = 1:30)
  
  spatial_psf_obj@images <- spatial_obj@images
  
  idents <- as.character(sort(unique(as.numeric(as.character(Idents(spatial_psf_obj))))))
  
  ### Finding Cluster specific markers
  cat("Identification of cluster specific features\n")
  ident_markers <- lapply(idents, function(x) {
    markers <- FindMarkers(spatial_psf_obj, ident.1 = x, logfc.threshold = 0.25, test.use = "wilcox")
    
    if("image" %in% names(pathway_collection[[1]]$attrs)) {
      sink_annot <- kegg_sink_to_process[gsub("-", "_", rownames(markers)), c("Pathway_name", "Sink", "Process", "Pathway")]
      colnames(sink_annot)[4] <- "Dowstream_pathway"
      cbind(markers, sink_annot)
    } else {
      markers$Pathway_name <- sapply(rownames(markers), function(y) {gsub("-", "_", unlist(strsplit(y, split = "; "))[2])})
      markers$Sink <- sapply(rownames(markers), function(y) {
        graph::nodeData(pathway_collection[[gsub("-", "_", unlist(strsplit(y, split = "; "))[2])]]$graph, n = unlist(strsplit(y, split = "; "))[1], attr = "label")
      })
    }
    
  })
  
  names(ident_markers) <- idents
  
  if(return_only_shiny_vars) {
    return(list(spatial_psf_collection = spatial_psf, psf_ident_markers = ident_markers, psf_mat = spatial_psf_mat, 
         spatial_image = spatial_psf_obj@images$slice1@image, coords = GetTissueCoordinates(object = spatial_psf_obj, scale = "lowres"), 
         meta.data = spatial_psf_obj@meta.data))
  } else {
    return(list(spatial_obj = spatial_obj, spatial_psf_obj = spatial_psf_obj, spatial_psf_collection = spatial_psf,
                psf_ident_markers = ident_markers))
  }
  
}