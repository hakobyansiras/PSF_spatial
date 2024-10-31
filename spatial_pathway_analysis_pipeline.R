library(Seurat)
library(SeuratData)
library(ggplot2)
library(patchwork)
library(dplyr)
library(biomaRt)
library(psf)
library(networkD3)
library(data.table)
library(DT)
library(miniUI)
library(shiny)

### Spatial clustering and pathway activity analysis function
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
                psf_ident_markers = ident_markers, psf_mat = spatial_psf_mat))
  }
  
}


## ggplot colors function to color sankey diagram nodes
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

## Silhouette plot function
custom_fviz_silhouette <- function (sil.obj, label = FALSE, print.summary = TRUE, clust_order = NULL, clust_colors, ...) {
  df <- as.data.frame(sil.obj, stringsAsFactors = TRUE)
  df <- df[order(df$cluster_num, -df$sil_width), ]
  # if (!is.null(rownames(df))) {
  #   df$name <- factor(rownames(df), levels = rownames(df))
  # } else {
  #   df$name <- as.factor(1:nrow(df)) 
  # }
  df$name <- as.factor(1:nrow(df)) 
  if(!is.null(clust_order)) {
    df$cluster <- factor(df$cluster, levels = clust_order) # as.factor(df$cluster)
  } else {
    df$cluster <- as.factor(df$cluster)
  }
  mapping <- aes_string(x = "name", y = "sil_width", color = "cluster", 
                        fill = "cluster")
  p <- ggplot(df, mapping) + geom_bar(stat = "identity") + 
    labs(y = "Silhouette width Si", x = "", 
         title = paste0("Clusters silhouette plot ", "\n Average silhouette width: ", round(mean(df$sil_width), 2))) + 
    ggplot2::ylim(c(NA, 1)) + geom_hline(yintercept = mean(df$sil_width), linetype = "dashed", color = "red") + 
    scale_fill_manual(values=clust_colors) + scale_colour_manual(values=clust_colors) 
  p <- ggpubr::ggpar(p, ...)
  if (!label) 
    p <- p + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
  else if (label) 
    p <- p + theme(axis.text.x = element_text(angle = 45))
  ave <- tapply(df$sil_width, df$cluster, mean)
  n <- tapply(df$cluster, df$cluster, length)
  sil.sum <- data.frame(cluster = names(ave), size = n, ave.sil.width = round(ave, 2), stringsAsFactors = TRUE)
  if (print.summary) 
    print(sil.sum)
  p
}

## Silhouette score calcualtion function
silhouette_calc <- function(s_obj, clust_name_convert, clust_colors, clust_order) {
  
  # silhouette metric
  library(factoextra)
  library(cluster, quietly = TRUE)
  reduction <- "umap"
  dims <- 1:2
  dist.matrix <- dist(x = Embeddings(object = s_obj[[reduction]])[, dims])
  clusters <- s_obj$seurat_clusters
  sil <- silhouette(x = as.numeric(x = as.factor(x = clusters)), dist = dist.matrix)
  sil[,"cluster"] <- sil[,"cluster"] - 1
  sil[,"neighbor"] <- sil[,"neighbor"] - 1
  
  sil <- data.frame(sil)
  sil[,"cluster"] <- clust_name_convert[as.character(sil[,"cluster"])]
  sil$cluster_num <- as.integer(gsub("[A-Z]", "", sil$cluster))
  sil[,"neighbor"] <- clust_name_convert[as.character(sil[,"neighbor"])]
  
  custom_fviz_silhouette(sil, clust_order =  clust_order, clust_colors = clust_colors)
}


## sankey plot building function
sankey_builder <- function(links, node_colors = NULL) {
  
  links <- links[which(links$value > 0),]
  
  nodes <- data.frame(
    name=c(as.character(links$source), 
           as.character(links$target)) %>% unique()
  )
  
  if(is.null(node_colors)) {
    color_scale_js <- NULL
  } else {
    color_scale_js <- paste0(
      'd3.scaleOrdinal()',
      '.domain([', paste0('"', nodes$name, '"', collapse = ", "), '])',
      '.range([', paste0('"', node_colors[nodes$name], '"', collapse = ", "), '])'
    )
  }
  
  
  
  links$IDsource <- match(links$source, nodes$name)-1 
  links$IDtarget <- match(links$target, nodes$name)-1
  
  # Make the Network
  if(is.null(color_scale_js)) {
    sn <- sankeyNetwork(Links = links, Nodes = nodes,
                        Source = "IDsource", Target = "IDtarget",
                        Value = "value", NodeID = "name",
                        sinksRight=FALSE, fontSize = 20)
  } else {
    sn <- sankeyNetwork(Links = links, Nodes = nodes,
                        Source = "IDsource", Target = "IDtarget",
                        Value = "value", NodeID = "name", colourScale = color_scale_js,
                        sinksRight=FALSE, fontSize = 20)
  }
  
  
  sn$x$nodes <-
    sn$x$nodes %>% 
    mutate(is_source_node = name %in% links$source)
  
  htmlwidgets::onRender(
    sn,
    '
  function(el,x) {
  d3.select(el)
    .selectAll(".node text")
    .filter(function(d) { return d.is_source_node; })
    .attr("x", x.options.nodeWidth - 16)
    .attr("text-anchor", "end");
  
  d3.select(el)
    .selectAll(".node text")
    .filter(function(d) { return !d.is_source_node; })
    .attr("x", x.options.nodeWidth)
    .attr("text-anchor", "start");
  }
  '
  )
  
}

### loading saptial dataset
spatial_melanoma <- Load10X_Spatial("/data/spacial_melanoma_psf/spatial_data", 
                                    filename = "CytAssist_FFPE_Human_Skin_Melanoma_filtered_feature_bc_matrix.h5")


### Downloading gene symbol to entrez id conversion data
# ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl")
# gene_symbol_to_entrez <- getBM(attributes=c('entrezgene_id', 'hgnc_symbol'), mart = ensembl)
# gene_symbol_to_entrez <- gene_symbol_to_entrez[which(!is.na(gene_symbol_to_entrez$entrezgene_id)),]
# gene_symbol_to_entrez <- gene_symbol_to_entrez[which(!duplicated(gene_symbol_to_entrez$entrezgene_id)),]
# gene_symbol_to_entrez <- setNames(object = gene_symbol_to_entrez$entrezgene_id, nm = gene_symbol_to_entrez$hgnc_symbol)

## load the same conversion vector which was used in the paper to have an identical results
load("gene_symbol_to_entrez.RData")

## loading curated 40 KEGG pathway from psf package 
load(system.file("extdata", "kegg_curated_40_signalings.RData", package="psf"))
load(system.file("extdata", "kegg_sink_to_process.RData", package="psf"))


### Running spatial data clustering and pathway analysis
psf_spatial <- spatial_psf_analysis(spatial_obj = spatial_melanoma, pathway_collection = kegg_curated_40_signalings, 
                                    gene_symbol_to_entrez = gene_symbol_to_entrez, nthreads = 30)

## Setting manual labels
exp_clust_name_convert <- setNames(object = c("E14", "E1", "E6", "E2", "E3", "E15", "E4", "E12", "E8", "E9", "E5", "E11", "E10", "E7", "E13"), 
                                   nm = c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14"))

psf_clust_name_convert <- setNames(object = c("P1", "P4", "P6", "P3", "P2", "P7", "P5", "P8", "P9"),
                                   nm = c("0", "1", "2", "3", "4", "5", "6", "7", "8"))


psf_spatial$spatial_obj@meta.data$Cluster <- exp_clust_name_convert[as.character(psf_spatial$spatial_obj@meta.data$seurat_clusters)]
psf_spatial$spatial_obj@meta.data$Cluster <- factor(psf_spatial$spatial_obj@meta.data$Cluster, levels = paste0("E", 1:15))

psf_spatial$spatial_obj@meta.data$spot_colors <- setNames(object = gg_color_hue(15), nm = as.character(0:14))[as.character(psf_spatial$spatial_obj@meta.data$seurat_clusters)]


psf_spatial$spatial_psf_obj@meta.data$Cluster <- psf_clust_name_convert[as.character(psf_spatial$spatial_psf_obj@meta.data$seurat_clusters)]
psf_spatial$spatial_psf_obj@meta.data$Cluster <- factor(psf_spatial$spatial_psf_obj@meta.data$Cluster, levels = paste0("P", 1:9))

psf_spatial$spatial_psf_obj@meta.data$spot_colors <- setNames(object = gg_color_hue(9), nm = as.character(0:8))[as.character(psf_spatial$spatial_psf_obj@meta.data$seurat_clusters)]

## generating silhouette plots for expression and PSF clusters
silhouette_calc(s_obj = psf_spatial$spatial_obj, clust_name_convert = exp_clust_name_convert, 
                clust_colors = setNames(nm = unique(psf_spatial$spatial_obj$Cluster), 
                                        object = unique(psf_spatial$spatial_obj$spot_colors)),
                clust_order = paste0("E", 1:15)
)

silhouette_calc(s_obj = psf_spatial$spatial_psf_obj, clust_name_convert = psf_clust_name_convert, 
                clust_colors = setNames(nm = psf_spatial$spatial_psf_obj$Cluster,
                                        object = psf_spatial$spatial_psf_obj$spot_colors),
                clust_order = paste0("P", 1:9)
)

## Loading spot annotations from Schmid et al paper (https://doi.org/10.3390/cimb46050284)
load("Schmidt_annotation.RData")

psf_spatial$spatial_obj@meta.data$cancer_curated <- Schmidt_annotation[rownames(psf_spatial$spatial_obj@meta.data),"cancer_curated"]
psf_spatial$spatial_obj@meta.data$general_som_clusts <- Schmidt_annotation[rownames(psf_spatial$spatial_obj@meta.data),"general_som_clusts"]


## generating UMAP plot for expression clustering
exp_umap <- DimPlot(psf_spatial$spatial_obj, label = TRUE, pt.size = 0.3, 
                    label.size = 5, reduction = "umap", group.by = "Cluster",
                    cols = setNames(object = psf_spatial$spatial_obj@meta.data$spot_colors,
                                    nm = psf_spatial$spatial_obj@meta.data$Cluster)) + 
  ggtitle("Exp UMAP") +
  labs(x = "UMAP 1", 
       y = "UMAP 2") + 
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=16,face="bold"),
        # legend.text = element_blank()
  )

## generating UMAP plot for PSF clustering
psf_umap <- DimPlot(psf_spatial$spatial_psf_obj, label = TRUE, pt.size = 0.3, 
                    label.size = 5, reduction = "umap", group.by = "Cluster",
                    cols = setNames(object = psf_spatial$spatial_psf_obj@meta.data$spot_colors,
                                    nm = psf_spatial$spatial_psf_obj@meta.data$Cluster)) + 
  ggtitle("PSF UMAP") +
  labs(x = "UMAP 1", 
       y = "UMAP 2") + 
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=16,face="bold"),
        # legend.text = element_blank()
  )


## generating spatial plot with expression cluster annotation
exp_groups <- SpatialDimPlot(psf_spatial$spatial_obj, label = TRUE, label.size = 6, group.by = "Cluster",
                             cols = setNames(object = psf_spatial$spatial_obj@meta.data$spot_colors,
                                             nm = psf_spatial$spatial_obj@meta.data$Cluster), repel = T
) + theme(legend.position="none")

## generating spatial plot with PSF cluster annotation
psf_groups <- SpatialDimPlot(psf_spatial$spatial_psf_obj, label = TRUE, label.size = 6, group.by = "Cluster",
                             cols = setNames(object = psf_spatial$spatial_psf_obj@meta.data$spot_colors,
                                             nm = psf_spatial$spatial_psf_obj@meta.data$Cluster)
) + theme(legend.position="none")


## generating spatial plot for manually cuatrated general melanoma subtypes and cancer cell type annotations from Shcmidt. et al.
som_groups <- SpatialDimPlot(psf_spatial$spatial_obj, label = TRUE, label.size = 5, group.by = "general_som_clusts") + theme(legend.position="none")
cancer_groups <- SpatialDimPlot(psf_spatial$spatial_obj, label = TRUE, label.size = 5, group.by = "cancer_curated", repel = T) + theme(legend.position="none")

## rendering expression umap and spatial lot
exp_umap + exp_groups
## rendering PSF umap and spatial lot
psf_umap + psf_groups

## rendering Manually cuatrated general melanoma subtypes annotations from Shcmidt. et al., epxression and PSF clusters
som_groups + exp_groups + psf_groups

## rendering Manually cuatrated general melanoma subtypes and cancer cell type annotations from Shcmidt. et al.
som_groups + cancer_groups



## Combining PSF and exp clustering metadata information for sankey diagrams
spatial_metadata <- psf_spatial$spatial_obj@meta.data
spatial_metadata$psf_clusts <- as.character(psf_spatial$spatial_psf_obj@meta.data[rownames(spatial_metadata), "Cluster"])


exp_to_psf <- Reduce(rbind,
                     lapply(as.character(sort(unique(spatial_metadata$psf_clusts))), function(x) {
                       data.frame(source = as.character(names(table(spatial_metadata$Cluster[which(spatial_metadata$psf_clusts == x)]))),
                                  target = x,
                                  value = as.numeric(unname(table(spatial_metadata$Cluster[which(spatial_metadata$psf_clusts == x)]))) 
                       )
                     })
)


som_general_to_exp <- Reduce(rbind,
                             lapply(sort(unique(psf_spatial$spatial_obj@meta.data$general_som_clusts)), function(x) {
                               data.frame(source = x, 
                                          target = as.character(names(table(psf_spatial$spatial_obj$Cluster[which(psf_spatial$spatial_obj$general_som_clusts == x)]))),
                                          value = as.numeric(unname(table(psf_spatial$spatial_obj$Cluster[which(psf_spatial$spatial_obj$general_som_clusts == x)]))) 
                               )
                             })
)


### setting manual node colors to match with spatial clustering colors
node_colors <- c(setNames(nm = unique(psf_spatial$spatial_psf_obj$Cluster), object = unique(psf_spatial$spatial_psf_obj$spot_colors)),
                 setNames(nm = unique(psf_spatial$spatial_obj$Cluster), object = unique(psf_spatial$spatial_obj$spot_colors)),
                 setNames(nm = sort(unique(psf_spatial$spatial_obj$general_som_clusts)), object = gg_color_hue(length(sort(unique(psf_spatial$spatial_obj$general_som_clusts)))))
)


### rendering sankey diagram
sankey_builder(links = exp_to_psf, node_colors = node_colors)
sankey_builder(links = rbind(som_general_to_exp, exp_to_psf), node_colors = node_colors)


### plotting average idend pathway activity fold change values across all pathways
spatial_psf_matrix <- as.matrix(psf_spatial$spatial_psf_obj@assays$Spatial@data)

psf_mean_FC_by_ident <- Reduce(rbind,
                               lapply(levels(psf_spatial$spatial_psf_obj$Cluster), function(x) {
                                 
                                 clust_cells <- rownames(psf_spatial$spatial_psf_obj@meta.data)[which(as.character(psf_spatial$spatial_psf_obj@meta.data$Cluster) == x)]
                                 
                                 mean_fc = log2(rowMeans(spatial_psf_matrix[,clust_cells])/rowMeans(spatial_psf_matrix[,setdiff(colnames(spatial_psf_matrix), clust_cells)]))
                                 
                                 data.frame(Cluster = x, mean_fc_values = unname(mean_fc),
                                            pathway = names(mean_fc))
                                 
                               })
)

psf_mean_FC_by_ident$pathway <- sapply(psf_mean_FC_by_ident$pathway, function(x) {unlist(strsplit(x, split = "; "))[2]})

psf_mean_FC_by_ident$pathway <- gsub("_pathway", "", psf_mean_FC_by_ident$pathway)


ggplot(psf_mean_FC_by_ident, aes(x=Cluster, y=mean_fc_values, fill=Cluster)) +
  geom_boxplot() +
  geom_hline(yintercept = 0, colour = "red", linewidth = 0.2) +
  facet_wrap(vars(pathway), nrow = 6, ncol = 8) +
  scale_fill_manual(values = setNames(object = psf_spatial$spatial_psf_obj@meta.data$spot_colors,
                                        nm = psf_spatial$spatial_psf_obj@meta.data$Cluster)) +
  # scale_fill_viridis(discrete = TRUE, alpha=0.6, option="A") +
  # theme_ipsum() +
  labs(title = "PSF FC between clusters", x = "Cluster", y = "Mean log2 PSF FC", fill = "Ident\n") +
  theme(
    legend.position="none",
    plot.title = element_text(size=18, face = 'bold'),
    axis.title=element_text(size=14,face="bold"),
    strip.text = element_text(size = 11, face = "bold")
  )


### calculating average PSF activities for each PSF cluster
ident_coldata <- data.frame(cell_id = rownames(psf_spatial$spatial_psf_obj@meta.data), 
                            group = as.character(psf_spatial$spatial_psf_obj@meta.data$Cluster))

spatial_melanoma_mean_psf <- lapply(psf_spatial$spatial_psf_collection, function(x) {
  print(x$attrs$title)
  mean_PSF_mat <- sapply(sort(unique(ident_coldata$group)), function(y) {
    clust_cells <- ident_coldata$cell_id[which(ident_coldata$group == y)]
    rowMeans(x$psf_activities[,clust_cells])
  })
  
  mean_exp_mat <- sapply(sort(unique(ident_coldata$group)), function(y) {
    clust_cells <- ident_coldata$cell_id[which(ident_coldata$group == y)]
    rowMeans(x$exp_fc[,clust_cells])
  })
  
  x$psf_activities <- mean_PSF_mat
  x$exp_fc <- mean_exp_mat
  x
})


## Plotting pathways
plot_pathway(spatial_melanoma_mean_psf$Hippo_signaling_pathway, plot_type = "kegg", use_old_images = T, color_nodes = "psf_activities", sample_id = "P3")

plot_pathway(spatial_melanoma_mean_psf$NF_kappa_B_signaling_pathway, plot_type = "kegg", use_old_images = T, color_nodes = "exp_fc", sample_id = c("P2", "P5"), multi_color_nodes = T)


### manual spot selection and comparison of selected spot groups

### Customized linked dimplot shiny app for manual spot selection
modified_linked_dimplot <-function (object, dims = 1:2, reduction = NULL, image = NULL, 
                                    group.by = NULL, alpha = c(0.1, 1), combine = TRUE) 
{
  ui <- miniPage(gadgetTitleBar(title = "LinkedDimPlot", left = miniTitleBarButton(inputId = "reset", label = "Reset")), 
                 miniContentPanel(
                   actionButton("record_selection", label = "Record"),
                   fillRow(
                     plotOutput(outputId = "spatialplot", height = "100%", click = clickOpts(id = "spclick", clip = TRUE),hover = hoverOpts(id = "sphover", delay = 10, nullOutside = TRUE)), 
                     plotOutput(outputId = "dimplot", height = "100%", 
                                brush = brushOpts(id = "brush", delay = 10, clip = TRUE, resetOnNew = FALSE), 
                                click = clickOpts(id = "dimclick",clip = TRUE), 
                                hover = hoverOpts(id = "dimhover", delay = 10, nullOutside = TRUE)), height = "97%"
                   ), 
                   verbatimTextOutput(outputId = "info")
                   
                 )
  )
  image <- image %||% Seurat:::DefaultImage(object = object)
  cells.use <- Cells(x = object[[image]])
  reduction <- reduction %||% DefaultDimReduc(object = object)
  dims <- dims[1:2]
  dims <- paste0(Key(object = object[[reduction]]), dims)
  group.by <- group.by %||% "ident"
  group.data <- FetchData(object = object, vars = group.by, 
                          cells = cells.use)
  coords <- GetTissueCoordinates(object = object[[image]])
  embeddings <- Embeddings(object = object[[reduction]])[cells.use, 
                                                         dims]
  plot.data <- cbind(coords, group.data, embeddings)
  plot.data$selected_ <- FALSE
  Idents(object = object) <- group.by
  server <- function(input, output, session) {
    click <- reactiveValues(pt = NULL, invert = FALSE)
    plot.env <- reactiveValues(data = plot.data, alpha.by = NULL, selected_cells = c(), clust_num = 0)
    observeEvent(eventExpr = input$done, handlerExpr = {
      plots <- list(plot.env$spatialplot, plot.env$dimplot)
      if (combine) {
        plots <- wrap_plots(plots, ncol = 2)
      }
      stopApp(returnValue = plots)
    })
    observeEvent(eventExpr = input$reset, handlerExpr = {
      click$pt <- NULL
      click$invert <- FALSE
      session$resetBrush(brushId = "brush")
    })
    observeEvent(eventExpr = input$brush, handlerExpr = click$pt <- NULL)
    observeEvent(eventExpr = input$spclick, handlerExpr = {
      click$pt <- input$spclick
      click$invert <- TRUE
      
      clicked <- nearPoints(df = plot.data, coordinfo = if (click$invert) {
        Seurat:::InvertCoordinate(x = click$pt)
      }
      else {
        click$pt
      }, threshold = 10, maxpoints = 1)
      
      
      plot.env$selected_cells <- c(plot.env$selected_cells, rownames(clicked))
      
    })
    
    observeEvent(input$record_selection, {
      
      write(plot.env$selected_cells, file = paste0("cluster", plot.env$clust_num))
      
      plot.env$clust_num = plot.env$clust_num + 1
      plot.env$selected_cells <- c()
      
    })
    
    observeEvent(eventExpr = input$dimclick, handlerExpr = {
      click$pt <- input$dimclick
      click$invert <- FALSE
    })
    observeEvent(eventExpr = c(input$brush, # input$spclick, 
                               input$dimclick), handlerExpr = {
                                 plot.env$data <- if (is.null(x = input$brush)) {
                                   clicked <- nearPoints(df = plot.data, coordinfo = if (click$invert) {
                                     Seurat:::InvertCoordinate(x = click$pt)
                                   }
                                   else {
                                     click$pt
                                   }, threshold = 10, maxpoints = 1)
                                   if (nrow(x = clicked) == 1) {
                                     cell.clicked <- rownames(x = clicked)
                                     group.clicked <- plot.data[cell.clicked, group.by, 
                                                                drop = TRUE]
                                     idx.group <- which(x = plot.data[[group.by]] == 
                                                          group.clicked)
                                     plot.data[idx.group, "selected_"] <- TRUE
                                     plot.data
                                   }
                                   else {
                                     plot.data
                                   }
                                 }
                                 else if (input$brush$outputId == "dimplot") {
                                   brushedPoints(df = plot.data, brush = input$brush, 
                                                 allRows = TRUE)
                                 }
                                 else if (input$brush$outputId == "spatialplot") {
                                   brushedPoints(df = plot.data, brush = Seurat:::InvertCoordinate(x = input$brush), 
                                                 allRows = TRUE)
                                 }
                                 plot.env$alpha.by <- if (any(plot.env$data$selected_)) {
                                   "selected_"
                                 }
                                 else {
                                   NULL
                                 }
                               })
    output$spatialplot <- renderPlot(expr = {
      plot.env$spatialplot <- SingleSpatialPlot(data = plot.env$data, 
                                                image = object[[image]], col.by = group.by, pt.size.factor = 1.6, 
                                                crop = TRUE, alpha.by = plot.env$alpha.by) + 
        scale_alpha_ordinal(range = alpha) + NoLegend()
      plot.env$spatialplot
    })
    output$dimplot <- renderPlot(expr = {
      plot.env$dimplot <- SingleDimPlot(data = plot.env$data, 
                                        dims = dims, col.by = group.by, alpha.by = plot.env$alpha.by) + 
        scale_alpha_ordinal(range = alpha) + guides(alpha = "none")
      plot.env$dimplot
    })
    output$info <- renderPrint(expr = {
      cell.hover <- rownames(x = nearPoints(df = plot.data, 
                                            coordinfo = if (is.null(x = input[["sphover"]])) {
                                              input$dimhover
                                            }
                                            else {
                                              Seurat:::InvertCoordinate(x = input$sphover)
                                            }, threshold = 10, maxpoints = 1))
      if (length(x = cell.hover) == 1) {
        paste(cell.hover, paste("Group:", plot.data[cell.hover, 
                                                    group.by, drop = TRUE]), collapse = "<br />")
      }
      else {
        NULL
      }
    })
    
  }
  runGadget(app = ui, server = server)
}

# modified_linked_dimplot(psf_spatial$spatial_psf_obj)


### load previosly manually slected spot groups
load("psf_manual_clusters.RData")

manual_cluster_labels <- setNames(object = rep("non_selected", length(rownames(psf_spatial$spatial_psf_obj@meta.data))), nm = rownames(psf_spatial$spatial_psf_obj@meta.data))
manual_cluster_labels[names(psf_manual_clusters)] <- unname(psf_manual_clusters)


psf_spatial$spatial_psf_obj@meta.data$manual_selection <- manual_cluster_labels[rownames(psf_spatial$spatial_psf_obj@meta.data)]


cluster_selection_list <- lapply(unique(manual_cluster_labels), function(x) {
  names(manual_cluster_labels)[which(manual_cluster_labels == x)]
})

names(cluster_selection_list) <- unique(manual_cluster_labels)
names(cluster_selection_list) <- gsub("cluster", "c", names(cluster_selection_list))


### selecteing spot groups of interest
manual_selection_subset <- cluster_selection_list[c("c4", "c2", "c18", "c12", "c10", "c13", "c16", "c20")]
manual_selection_subset_vec <- unlist(lapply(names(manual_selection_subset), function(x) {
  setNames(object = rep(x, length(manual_selection_subset[[x]])), nm = manual_selection_subset[[x]])
}))


### visualizing selected spot groups
psf_spatial$spatial_psf_obj@meta.data$custom_selection <- as.character(psf_spatial$spatial_psf_obj@meta.data$seurat_clusters)
psf_spatial$spatial_psf_obj@meta.data[names(manual_selection_subset_vec),"custom_selection"] <- unname(manual_selection_subset_vec)

cluster_selection_list_subset <- sapply(unique(psf_spatial$spatial_psf_obj@meta.data$custom_selection), function(x) {
  rownames(psf_spatial$spatial_psf_obj@meta.data)[which(psf_spatial$spatial_psf_obj@meta.data$custom_selection == x)]
})


### Spatial plot of 
SpatialDimPlot(psf_spatial$spatial_psf_obj, label = F, label.size = 6,
               cols.highlight =  c(gg_color_hue(9), c("#fcc244", "#fcc244", "#e33a14", "#e33a14", "#006400", "#006400", "#fcc244", "#006400")),
               cells.highlight = cluster_selection_list_subset[order(names(cluster_selection_list_subset))])



### performingDetecting significantly differentially activated pathways branches between border and core spot groups for each cluster
psf_clust_comparisons <- list(c_12_10 = c("c12", "c10"),
                              c_18_4 = c("c18", "c4"),
                              c_20_10 = c("c20", "c10"),
                              c_13_16 = c("c13", "c16")
)

psf_spatial$spatial_psf_obj@meta.data$manual_selection <- gsub("cluster", "c", psf_spatial$spatial_psf_obj@meta.data$manual_selection)

psf_adjesent_clust_comparison <- lapply(psf_clust_comparisons, function(y) {
  
  clust_0 <- psf_spatial$psf_mat[,rownames(psf_spatial$spatial_psf_obj@meta.data)[which(psf_spatial$spatial_psf_obj@meta.data$manual_selection == y[1])]]
  clust_1 <- psf_spatial$psf_mat[,rownames(psf_spatial$spatial_psf_obj@meta.data)[which(psf_spatial$spatial_psf_obj@meta.data$manual_selection == y[2])]]
  
  clust_0_and_1_psf <- cbind(clust_0, clust_1)
  
  gr <- factor(c(rep("clus0", times=ncol(clust_0)), rep("clus1", times=ncol(clust_1))), levels=c("clus0", "clus1"))
  
  p.val <- apply(clust_0_and_1_psf, 1, function(x) {
    fit <- lm(x~gr)
    summary(fit)$coefficients[2,4]
  })
  
  stat_df <- data.frame(sink = rownames(clust_0_and_1_psf), log2_fc = log2(rowMeans(clust_0)/rowMeans(clust_1)),  p_val = p.val, fdr = p.adjust(p.val, method = "fdr"))
  
  stat_df <- stat_df[order(stat_df$fdr),]
  
  stat_df[which(stat_df$fdr < 0.05),]
  
})

### annotating significant branches
psf_adjesent_clust_comparison <- lapply(psf_adjesent_clust_comparison, function(x) {
  if(nrow(x) > 0) {
    cbind(x, sink_id = kegg_sink_to_process[x$sink, "sink_id"], 
          sink_name = kegg_sink_to_process[x$sink, "Sink"],
          pathway = kegg_sink_to_process[x$sink, "Pathway_name"],
          process = kegg_sink_to_process[x$sink, "Process"], 
          downstream_pathway = kegg_sink_to_process[x$sink, "Pathway"])
  } else {
    x
  }
})


### Detecting receptor ligand connections between KEGG signaling pathways
determine.input.nodes <- function(pathway) {
  dfs <- graphnel_to_df(pathway)
  return(setdiff(dfs$edge_table$from, dfs$edge_table$to))
}

kegg_pathway_connection_list <- lapply(names(kegg_curated_40_signalings), function(x) {
  
  sink_genes <- graph::nodeData(kegg_curated_40_signalings[[x]]$graph, 
                                n = kegg_curated_40_signalings[[x]]$sink.nodes, attr = "genes")
  
  connections <- lapply(names(kegg_curated_40_signalings), function(y) {
    
    input_genes <- graph::nodeData(kegg_curated_40_signalings[[y]]$graph, 
                                   n = determine.input.nodes(kegg_curated_40_signalings[[y]]), attr = "genes")
    
    Reduce(rbind, lapply(names(sink_genes), function(z) {
      Reduce(rbind, lapply(names(input_genes), function(h) {
        if(length(intersect(sink_genes[[z]], input_genes[[h]])) > 0) {
          data.frame(sink = paste0(z, "; ", x), input_node = h, input_pathway = y)
          # paste0(z, " > ", h)
        }
      }))
      
    }))   
    
  })
  
  connections[!sapply(connections, is.null)]
  
  Reduce(rbind, connections)
  
})

kegg_pathway_connection_df <- Reduce(rbind, kegg_pathway_connection_list)


### Checking if there are interacting significant pathways in adjacent border spot groups of tow clusters
adjacent_border_pairs <- list(
  c("c_18_4", "c_12_10"),
  c("c_13_16", "c_20_10")
)


lapply(adjacent_border_pairs, function(x) {
  filtered_connection_df <- kegg_pathway_connection_df[which(kegg_pathway_connection_df$sink %in% rownames(psf_adjesent_clust_comparison[[x[1]]])),]
  filtered_connection_df[which(filtered_connection_df$input_pathway %in% psf_adjesent_clust_comparison[[x[2]]]$pathway),]
})


lapply(adjacent_border_pairs, function(x) {
  filtered_connection_df <- kegg_pathway_connection_df[which(kegg_pathway_connection_df$sink %in% rownames(psf_adjesent_clust_comparison[[x[2]]])),]
  filtered_connection_df[which(filtered_connection_df$input_pathway %in% psf_adjesent_clust_comparison[[x[1]]]$pathway),]
})
