# QR-SIDE
Robust spatial cell-type deconvolution with qualitative reference for spatial transcriptomics

QR-SIDE is a package for cell type deconvolution on spatial transcriptomics datasets developed by Jin Liu's lab. It has the capability to effectively leverage a large number of non-marker genes as well as “qualitative” information about marker genes without using a reference dataset. QR-SIDE can quantify the spatial heterogeneity of individual marker genes using a mixture Poisson regression mode, and simultaneously perform spatial clustering by specifying the latent spot-separable topics. Using a Potts model, QR-SIDE integrates spatial information to promote spatial continuity for the identified spatial domains. 

# Dependency
The required dependencies for QR-SIDE is listed as following:
dplyr : 1.1.4 
clue : 0.3-65 
mclust : 6.1 
STdeconvolve : 1.3.1 
scater : 1.26.1 
scran : 1.26.2 
scuttle : 1.8.4 
SingleCellExperiment : 1.20.1 
SummarizedExperiment : 1.28.0 
Biobase : 2.58.0 
GenomicRanges : 1.50.2 
GenomeInfoDb : 1.34.9 
IRanges : 2.32.0 
S4Vectors : 0.36.2 
BiocGenerics : 0.44.0 
MatrixGenerics : 1.10.0 
matrixStats : 1.1.0 
BiocSingular : 1.14.0 
MASS : 7.3-60.0.1 
Rcpp : 1.0.12 
Seurat : 5.0.1 
SeuratObject : 5.0.1 
sp : 2.1-3 
data.table : 1.15.0 

Additionally, the package `SpatialDecon` is also required. The installation package `SpatialDecon.tar.gz` is included in the R folder.

# Dependency
The inputs of QR-SIDE include sp_expr, sp_pos

sp_expr=filtered_matrix # Spot by gene SRT matrix

sp_pos=aligned_pos #Spot by 2, 2D spatial coordination

top_DEGs=top_DEGs #DEGs dataframe from Seurat sorted by 'avg_log2FC'

Num_topic=7 #Num of topic domains

Num_HVG=1000 #num of highly variable genes used in QR_SIDE

dim_embed=20 #dimension of embeddings in hierachical factor models

top_marker_num=4 #Num of markers of each cell type from "top_DEGs" used in QR_SIDE 

fixed_marker_list=FALSE #whether QR-SIDE uses the all genes in "top_DEGs" as markers. If FALSE, "top_marker_num" must be provided

out=QR_SIDE(sp_expr, sp_pos, top_DEGs, Num_topic, Num_HVG, dim_embed, top_marker_num, fixed_marker_list)
