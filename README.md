# QR-SIDE
Robust spatial cell-type deconvolution with qualitative reference for spatial transcriptomics

QR-SIDE is a package for cell type deconvolution on spatial transcriptomics datasets developed by Jin Liu's lab. It has the capability to effectively leverage a large number of non-marker genes as well as “qualitative” information about marker genes without using a reference dataset. QR-SIDE can quantify the spatial heterogeneity of individual marker genes using a mixture Poisson regression mode, and simultaneously perform spatial clustering by specifying the latent spot-separable topics. Using a Potts model, QR-SIDE integrates spatial information to promote spatial continuity for the identified spatial domains. 

# Dependency
The required dependencies for QR-SIDE is listed as following:

`clue : 0.3-65`

`mclust : 6.1`

`STdeconvolve : 1.3.1` 

`scater : 1.26.1` 

`scran : 1.26.2` 

`SingleCellExperiment : 1.20.1` 

`BiocSingular : 1.14.0` 

`MASS : 7.3-60.0.1`

`Rcpp : 1.0.12`

Additionally, the package `SpatialDecon` is also required. The installation package `SpatialDecon.tar.gz` is included in the R folder.

# Usage
The inputs of QR-SIDE include a Spot by gene SRT matrix `sp_expr`, 2D spatial coordination `sp_pos`, and differentially expressed gene list `top_DEGs`. `Num_topic` is the number of spatial domains. `Num_HVG` is the number of highly variable genes involved in training. `dim_embed` is the factor dimension in the hierachical factor models. If you only want top n genes of each cell type in `top_DEGs` to be involved, you can specific the value of `top_marker_num = n` and turn `fixed_marker_list` to **FALSE** to focus on the top n marker genes only.

`out=QR_SIDE(sp_expr, sp_pos, top_DEGs, Num_topic, Num_HVG, dim_embed, top_marker_num, fixed_marker_list)`
