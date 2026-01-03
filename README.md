# zebrafish SpatioTemperal Expression Profiles (zSTEP)

### Yang Dong<sup>†</sup>,  Tao Cheng<sup>†</sup>, Xiang Liu<sup>†</sup>,  ..., Jie Liao*, Xiaohui Fan*, Pengfei Xu*
[![DOI](https://img.shields.io/badge/DOI-10.1186/s13059--025--03917--8-blue)](https://link.springer.com/article/10.1186/s13059-025-03917-8) [![DOI](https://img.shields.io/badge/DOI-10.1101/2024.07.01.601472-yellowgreen)](https://www.biorxiv.org/content/10.1101/2024.07.01.601472v1)
<p align="justify">zSTEP is a zebrafish spatial gene expression atlas that provides a powerful platform for visualizing gene expression patterns in real zebrafish embryonic cartography and investigating spatial cell-cell interactions within selected regions. As a result, zSTEP holds significant potential for investigating a wide range of developmental events during development. </p>
<p align="justify">The generation of zSTEP involves performing bulk RNA-seq on serial cryosections taken along the left-right axis, implementation of Palette to infer spatial gene expression using Stereo-seq data<sup>1</sup> as the reference and projection of the ST spots onto the zebrafish live images<sup>2</sup>.</p><br>

<p align="center">
  <img width="800"  src="https://github.com/ldo2zju/zSTEP/blob/main/images/DreAM.png">
</p>

# Palette
<p align="justify">Palette is a bioinformatics pipeline designed to infer spatial gene expression from bulk RNA-seq data. Only spatial transcriptomics data is required as the reference, and bulk RNA-seq data is used as the input. The pipeline incorporates imputation, smoothing and expression amplification, resulting in accurate predictions of spatial expression patterns.</p><br>
<p align="center">
  <img width="1200"  src="https://github.com/ldo2zju/zSTEP/blob/main/images/Palette pipeline_Final.png">
</p>

## Requirements and Installation
Palette is a R-programming pipeline, and employs BayesSpace<sup>3</sup> for spatial clustering and MuSiC<sup>4</sup> for deconvolution. The following codes show the installation of these two R packages.

```{r, eval = FALSE}
# install BiocManager if necessary
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager", repos = "https://cran.R-project.org")

# install BayesSpace
BiocManager::install("BayesSpace")

# install devtools if necessary
install.packages('devtools')

# install MuSiC
devtools::install_github('xuranw/MuSiC')
```

## Tutorial
Here we provide several step-by-step tutorials using _Drosophila_ Stereo-seq data<sup>5</sup>, MERFISH ST reference<sup>6</sup> and Visium reference<sup>7</sup> as examples. The data is available in the _test_ folder.
* [Inferring spatial gene expression patterns from _Drosophila_ pseudo bulk data](https://rawcdn.githack.com/ldo2zju/zSTEP/21de0ea373cbdd03ff04fa7688c00898e577cb0f/analysis/Palette_Drosophila_tutorial_update.html)
* [Inferring spatial gene expression patterns from mouse hypothalamus bulk (MERFISH as ST reference)](https://rawcdn.githack.com/ldo2zju/zSTEP/37773124553cd7262b3be06ddcdc6d333576b6bd/analysis/Palette_mouse_hypothalamus_tutorial(MERFISH).html)
* [Inferring spatial gene expression patterns from melanoma pseudo bulk data (Visium as ST reference)](https://rawcdn.githack.com/ldo2zju/zSTEP/37773124553cd7262b3be06ddcdc6d333576b6bd/analysis/Palette_melanoma_tutorial(Visium).html)


## About
If you have any questions, please feel free to contact the following authors responsible for the Palette design, Dr. Yang Dong (yang.dongau@gmail.com), Dr. Tao Cheng (chengtao2ldo@zju.edu.cn) or Dr. Jie Liao (liaojie@zju.edu.cn).

## References
1. Liu, C. et al. Spatiotemporal mapping of gene expression landscapes and developmental trajectories during zebrafish embryogenesis. _Dev Cell_ 57, 1284-1298 e1285 (2022).
2. Shah, G. et al. Multi-scale imaging and analysis identify pan-embryo cell dynamics of germlayer formation in zebrafish. _Nat Commun_ 10, 5753 (2019).
3. Zhao, E. et al. Spatial transcriptomics at subspot resolution with BayesSpace. _Nature Biotechnology_ 39, 1375-+ (2021).
4. Wang, X., Park, J., Susztak, K., Zhang, N.R. & Li, M. Bulk tissue cell type deconvolution with multi-subject single-cell expression reference. _Nat Commun_ 10, 380 (2019).
5. Wang, M. et al. High-resolution 3D spatiotemporal transcriptomic maps of developing Drosophila embryos and larvae. _Dev Cell_ 57, 1271-1283 e1274 (2022).
6. Moffitt, J. et al. Molecular, spatial, and functional single-cell profiling of the hypothalamic preoptic region. _Science_ 362, 6416 (2018).
7. Thrane, K. et al. Spatially resolved transcriptomics enables dissection of genetic heterogeneity in stage III cutaneous malignant melanoma. _Cancer Res_ 78, 5970–5979 (2018).

## Citation
<p align="justify">Our work has been published on <em>Genome Biology</em> and can be cited as follows.</p>

Dong, Y. et al. Unravelling the progression of the zebrafish primary body axis with reconstructed spatiotemporal transcriptomics. _Genome Biol_ (2026) https://doi.org/10.1186/s13059-025-03917-8



