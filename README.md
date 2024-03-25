# _Danio rerio_ Asymmetrical Maps (<em>Dre</em>AM)

### Yang Dong<sup>†</sup>,  Tao Cheng<sup>†</sup>, Xiang Liu<sup>†</sup>,  ..., Jie Liao*, Xiaohui Fan*, Pengfei Xu*
<p align="justify"><em>Dre</em>AM is a zebrafish spatial gene expression atlas that provides a powerful platform for visualizing gene expression patterns in real zebrafish embryonic cartography and investigating spatial cell-cell interactions within selected regions. As a result, <em>Dre</em>AM holds significant potential for investigating a wide range of developmental events during development. </p>
<p align="justify">The generation of <em>Dre</em>AM involves performing bulk RNA-seq on serial cryosections taken along the left-right axis, implementation of Palette to infer spatial gene expression using Stereo-seq data<sup>1</sup> as the reference and projection of the ST spots onto the zebrafish live images<sup>2</sup>.</p><br>

<p align="center">
  <img width="700"  src="https://github.com/ldo2zju/DreAM/blob/main/images/DreAM.png">
</p>

# Palette
<p align="justify">Palette is a bioinformatics pipeline designed to infer spatial gene expression from bulk RNA-seq data. Only spatial transcriptomics data is required as the reference, and bulk RNA-seq data is used as the input. The pipeline incorporates imputation, smoothing and expression amplification, resulting in accurate predictions of spatial expression patterns.</p><br>
<p align="center">
  <img width="1200"  src="https://github.com/ldo2zju/DreAM/blob/main/images/Palette%20workflow.png">
</p>

## Requirements and Installation


## Tutorial
Here we provide a step-by-step tutorial



## About
<p align="justify">If you have any questions, please feel free to contact the following authors responsible for the Palette design, Dr. Yang Dong (yang.dongau@gmail.com), Dr. Tao Cheng (chengtao2ldo@zju.edu.cn) or Dr. Jie Liao (liaojie@zju.edu.cn).</p>

## Reference
1. Liu, C. et al. Spatiotemporal mapping of gene expression landscapes and developmental trajectories during zebrafish embryogenesis. _Dev Cell_ 57, 1284-1298 e1285 (2022).
2. Shah, G. et al. Multi-scale imaging and analysis identify pan-embryo cell dynamics of germlayer formation in zebrafish. _Nat Commun_ 10, 5753 (2019).



