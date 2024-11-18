# PSF Spatial Browser

This repository contains the code for running the PSF spatial browser locally. The PSF spatial browser is an interactive R Shiny application designed to explore the results of topology-aware pathway analysis of spatial transcriptomic data using the Pathway Signal Flow (PSF) algorithm.

## Installation

To install the PSF spatial browser, follow these steps:

1.  Clone the GitHub repository:

    ``` bash
    git clone https://github.com/hakobyansiras/PSF_spatial.git
    ```

2.  Install the necessary R libraries:

    ``` r
    remotes::install_github("hakobyansiras/psf")
    if (!require("BiocManager", quietly = TRUE)) {install.packages("BiocManager")}
    BiocManager::install("biomaRt")
    install.packages(c("ggplot2", "patchwork", "dplyr", "data.table", "DT", "miniUI", "shiny", "Seurat", "psf", "magick", "shinyjs", "visNetwork", "biomaRt", "SeuratData"))
    ```

3.  Load the required R functions:

    ``` r
    source("Spatial_psf_analysis.R")
    source("psf_spatial_browser.R")
    ```

## Getting Started

### Loading Spatial Dataset

Load your spatial dataset using the `Load10X_Spatial` function:

``` r
spatial_melanoma <- Load10X_Spatial("/data/spacial_melanoma_psf/spatial_data",
                                    filename = "CytAssist_FFPE_Human_Skin_Melanoma_filtered_feature_bc_matrix.h5")
```

### Downloading Gene Symbol to Entrez ID Conversion Data

You can download gene symbol to Entrez ID conversion data or use the preloaded version:

``` r
# ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl")
# gene_symbol_to_entrez <- getBM(attributes=c('entrezgene_id', 'hgnc_symbol'), mart = ensembl)
# gene_symbol_to_entrez <- gene_symbol_to_entrez[which(!is.na(gene_symbol_to_entrez$entrezgene_id)),]
# gene_symbol_to_entrez <- gene_symbol_to_entrez[which(!duplicated(gene_symbol_to_entrez$entrezgene_id)),]
# gene_symbol_to_entrez <- setNames(object = gene_symbol_to_entrez$entrezgene_id, nm = gene_symbol_to_entrez$hgnc_symbol)
load("gene_symbol_to_entrez.RData")
```

### Loading KEGG Signaling Pathways

Load the KEGG signaling pathways:

``` r
load(system.file("extdata", "kegg_curated_40_signalings.RData", package="psf"))
load(system.file("extdata", "kegg_sink_to_process.RData", package="psf"))
```

### Running Pathway Analysis

Run the pathway analysis using the `spatial_psf_analysis` function:

``` r
psf_spatial <- spatial_psf_analysis(spatial_obj = spatial_melanoma, pathway_collection = kegg_curated_40_signalings,
                                    gene_symbol_to_entrez = gene_symbol_to_entrez, nthreads = 30)
```

### Launching the App

Launch the PSF spatial browser application:

``` r
run_psf_spatial_browser(psf_saptial_results = psf_spatial)
```

## App Usage

To see a demonstration of how to use the PSF spatial browser, you can watch the short demo video below: [PSF Spatial Browser Demo](https://youtu.be/lHTgYBA374o)
