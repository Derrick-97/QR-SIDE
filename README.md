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
Additional_repositories: 
  https://bioconductor.org/packages/release/bioc
RoxygenNote: 7.3.1

# Installation Guide
1. Download the following files to your local machine:
    SpatialDecon_1.0.tar.gz
    QRSIDE_1.0.0.tar.gz
   
2. Install SpatialDecon
    Open ​​R​​ or ​​RStudio​​ and run:

    # Replace '/path/to/' with the actual directory containing SpatialDecon_1.0.tar.gz
    install.packages(
      pkgs = "/path/to/SpatialDecon_1.0.tar.gz",
      repos = NULL,
      type = "source"
    )
3. Install QRSIDE
In the same R session, run:

# Replace '/path/to/' with the actual directory containing QRSIDE_1.0.0.tar.gz
install.packages(
  pkgs = "/path/to/QRSIDE_1.0.0.tar.gz",
  repos = NULL,
  type = "source"
)

# Troubleshooting
Install required dependencies manually. For example:

  install.packages("dplyr")  # CRAN package
  if (!require("BiocManager")) install.packages("BiocManager")
  BiocManager::install("SingleCellExperiment")  # Bioconductor package


# Usage
The inputs of QR-SIDE include a Spot by gene SRT matrix `sp_expr`, 2D spatial coordination `sp_pos`, and differentially expressed gene list `top_DEGs`. `Num_topic` is the number of spatial domains. `Num_HVG` is the number of highly variable genes involved in training. `dim_embed` is the factor dimension in the hierachical factor models. If you only want top n genes of each cell type in `top_DEGs` to be involved, you can specific the value of `top_marker_num = n` and turn `fixed_marker_list` to **FALSE** to focus on the top n marker genes only.

