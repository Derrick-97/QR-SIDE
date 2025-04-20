# QR-SIDE
Robust spatial cell-type deconvolution with qualitative reference for spatial transcriptomics

Package: QRSIDE

Version: 1.0.0

Title: Robust Spatial Cell-Type Deconvolution with Qualitative Reference for Spatial Transcriptomics
Description: QR-SIDE is a package for cell type deconvolution on spatial transcriptomics datasets developed by Jin Liu's lab. It has the capability to effectively leverage a large number of non-marker genes as well as "qualitative" information about marker genes without using a reference dataset. QR-SIDE can quantify the spatial heterogeneity of individual marker genes using a mixture Poisson regression model, and simultaneously perform spatial clustering by specifying the latent spot-separable topics. Using a Potts model, QR-SIDE integrates spatial information to promote spatial continuity for the identified spatial domains.

# Dependency
License: GPL-3 

Encoding: UTF-8

Depends: 
    R (>= 4.2.3)
    
Imports:
    BiocSingular (>= 1.14.0),
    MASS (>= 7.3-60),
    Rcpp (>= 1.0.12),
    STdeconvolve (>= 1.3.1),
    SingleCellExperiment (>= 1.20.1),
    clue (>= 0.3-65),
    combinat,
    mclust (>= 6.1),
    scater (>= 1.26.1),
    scran (>= 1.26.2),
    SpatialDecon
    
LinkingTo: 
    Rcpp,
    BH,
    RcppArmadillo
    
SystemRequirements: C++11

Suggests: 
    testthat (>= 3.0.0)
    
RoxygenNote: 7.3.1

# Installation Guide
1. Download the following files to your local machine:
    SpatialDecon_1.0.tar.gz
    QRSIDE_1.0.0.tar.gz
   
2. Install SpatialDecon
    Open ​​R​​ or ​​RStudio​​ and run:
```r
#Replace '/path/to/' with the actual directory containing SpatialDecon_1.0.tar.gz
install.packages(
  pkgs = "/path/to/SpatialDecon_1.0.tar.gz",
  repos = NULL,
  type = "source"
)
```
3. Install QRSIDE
In the same R session, run:
```r
#Replace '/path/to/' with the actual directory containing QRSIDE_1.0.0.tar.gz
install.packages(
  pkgs = "/path/to/QRSIDE_1.0.0.tar.gz",
  repos = NULL,
  type = "source"
)
```

# Troubleshooting
Install required dependencies manually. For example:
```r
install.packages("dplyr")  # CRAN package
if (!require("BiocManager")) install.packages("BiocManager")
BiocManager::install("SingleCellExperiment")  # Bioconductor package
```

# Usage
QR-SIDE​​ takes the following inputs to perform spatial domain deconvolution:

Key Inputs:

sp_expr​​: A ​​spot-by-gene matrix​​ (spatial transcriptomics data in matrix/dataframe format).
​​sp_pos​​: A ​​2D spatial coordinate matrix​​ (spot locations in X-Y coordinates).
​​top_DEGs​​: A ​​list of differentially expressed genes​​ (cell-type marker genes).
​​Num_topic​​ (int): Number of spatial domains to infer.
​​Num_HVG​​ (int): Number of highly variable genes (HVGs) to include in training.
​​dim_embed​​ (int): Latent dimension for hierarchical factor modeling.
​​top_marker_num​​ (int): Only use the ​​top n marker genes​​ per cell type from top_DEGs.
​​fixed_marker_list​​ (logical):

FALSE → Use top top_marker_num genes per cell type.
TRUE → Use all genes in top_DEGs.

For a quick start example, see the 'tutorial/MOB.ipynb'

