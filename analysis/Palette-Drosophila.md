Inferring spatial gene expression patterns from *Drosophila* pseudo bulk
data
================
Yang Dong
2024-03-26

<style type="text/css">
  body{
  font-size: 12pt;
}
</style>

### Setup

``` r
library(limma)
library(Seurat)
library(readr)
library(tidyverse)
library(magrittr)
library(pander)
library(SingleCellExperiment)
library(reshape2)
library(BayesSpace)
library(MuSiC)
library(ggplot2)
library(ggpubr)
library(Matrix.utils)

theme_set(theme_bw())
panderOptions("big.mark", ",")
panderOptions("table.split.table", Inf)
panderOptions("table.style", "rmarkdown")
if (interactive()) setwd(here::here("analysis"))
```

### Summary

Here we take 2 neighboring ST slices (slice 4 and slice 5) from
*Drosophila* ST data. Slice 4 will be used as the ST reference, while
slice 5 will be converted a pseudo bulk data. Palette will then be
implemented to infer the expression patterns from the pseudo bulk data
using the slice 4 as the ST reference. The inferred expression patterns
will then be compared with the real expression of slice 5 ST data.

### Load data

Load the neighboring *Drosophila* ST slice.

``` r
Drosophila_slice <- readRDS(here::here("test","Drosophila_slice.rds"))
```

### *Drosophila* ST slices

Here we visualize the ST slice and each ST spot is colored with spot
annotation. The two ST slices show similar cell type distribution.

``` r
# Assign color to each cell type
color_df <- data.frame(unique(Drosophila_slice$annotation))
colnames(color_df)[1] <- "annotation"
color_df$color <- c("#FFF799","#2ED9FF","#683b79","#009200","#6e9BC5","#F091A0","#FEAF16","#ff66ff","#AEF359","#C0C0C0")
x <- Drosophila_slice@meta.data %>%
  left_join(color_df)
Drosophila_slice@meta.data$color <- x$color

slice_4 <- subset(Drosophila_slice, subset = slice_ID == "E14-16h_a_S04")
slice_5 <- subset(Drosophila_slice, subset = slice_ID == "E14-16h_a_S05")

# Visualize ST slice
g1 <- ggplot(slice_4@meta.data %>% as.data.frame(),aes(x=new_x, y=new_y, color=color)) +
  geom_point(size=0.5) +
  theme_classic() +
  xlim(-25, 25) +
  ylim(-25, 25) +
    coord_fixed(ratio = 1, xlim = NULL, ylim = NULL, expand = TRUE, clip = "on") +
  scale_color_manual(values = levels(slice_4@meta.data$color %>% as.factor())) +
  theme_classic() +
  theme(axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        legend.position="none")
g2 <- ggplot(slice_5@meta.data %>% as.data.frame(),aes(x=new_x, y=new_y, color=color)) +
  geom_point(size=0.5) +
  theme_classic() +
  xlim(-25, 25) +
  ylim(-25, 25) +
    coord_fixed(ratio = 1, xlim = NULL, ylim = NULL, expand = TRUE, clip = "on") +
  scale_color_manual(values = levels(slice_5@meta.data$color %>% as.factor())) +
  theme_classic() +
  theme(axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        legend.position="none")

ggarrange(g1, g2, ncol = 2)
```

<div class="figure" style="text-align: center">

<img src="1.-Palette-pipeline--Drosophila-slices-as-example-_files/figure-gfm/Drosophila ST slices-1.png" alt="The adjacent _Drosophila_ ST slices." width="60%" />
<p class="caption">
The adjacent *Drosophila* ST slices.
</p>

</div>

<p>
 
</p>

### Construct pseudo bulk data and select highly expressed genes

Aggregate expression of ST data to generate pseudo bulk. Slice 5 is
convert to pseudo bulk, while the aggregated expression of ST data
(slice 4) is used for selecting highly expressed genes.

``` r
# Convert both to pseudo bulk data
slice_sce <- SingleCellExperiment(assays = list(counts = Drosophila_slice@assays$RNA@counts), 
                           colData = Drosophila_slice@meta.data)
groups <- colData(slice_sce)[, c("slice_ID")]
slice_bulk <- aggregate.Matrix(t(counts(slice_sce)), 
                       groupings = groups, fun = "sum") %>% 
  as.data.frame() %>% t()

# Slice 5 as pseudo bulk
slice_5_bulk <- slice_bulk[,2, drop = FALSE] %>% set_colnames("slice_5")
# Slice 4 as aggregated ST data for selecting highly expressed genes
slice_4_ST_sum <- slice_bulk[,1, drop = FALSE] %>% set_colnames("slice_4_ST")

# Select highly expressed genes in ST data
slice_4_ST_sum_high_exp <- slice_4_ST_sum %>%
  as.data.frame() %>%
  dplyr::filter(rownames(.) %in% rownames(slice_5_bulk)) %>%
  dplyr::filter(.$slice_4_ST > 30)

# Reorder genes of pseudo bulk data (slice 5) based on the ST data (slice 4)
slice_5_bulk_reorder <- slice_5_bulk[match(rownames(slice_4_ST_sum_high_exp), rownames(slice_5_bulk)),,drop = FALSE] %>%
  as.data.frame()

# Define the ratio of the pseudo bulk (slice 5) and the aggregate ST (slice 4) as gene filter
gene_filter <- (slice_5_bulk_reorder/slice_4_ST_sum_high_exp) %>%
  set_colnames("ratio") %>%
  .[is.finite(rowSums(.)),, drop = FALSE]
```

We use density plot to visualize the expression ratios of the pseudo
bulk (slice 5) and the aggregated ST (slice 4). We take the top 2000
genes closest to the mean and median ratios for the following analyses.
These genes are highly expressed in both pseudo bulk and ST data.

``` r
ggplot(gene_filter, aes(x=ratio)) + 
  geom_density() + theme_classic() +
  theme(text = element_text(size = 20)) 
```

<div class="figure" style="text-align: center">

<img src="1.-Palette-pipeline--Drosophila-slices-as-example-_files/figure-gfm/unnamed-chunk-2-1.png" alt="Density plot showing the expression ratio of pseudo bulk and aggregated ST data." width="40%" />
<p class="caption">
Density plot showing the expression ratio of pseudo bulk and aggregated
ST data.
</p>

</div>

``` r
# The difference of gene ratio to the mean and median ratio
gene_filter$differ <- abs(gene_filter$ratio - (mean(gene_filter$ratio)+median(gene_filter$ratio))/2)

# Select the top 2000 genes closer to the mean and median for spatial clustering and deconvolution
high_exp_genes <- gene_filter %>%
  .[order(.$differ),] %>%
   top_n(.,-2000)
```

<p>
 
</p>

### Spatial clustering of ST data

BayesSpace is employed for spatial clustering on ST data (slice 4). Only
the selected 2000 highly expressed genes are used for spatial
clustering. The cluster identities are then assigned to each ST spot.

``` r
# Convert coordinates into pixel coordinates for BayesSpace
metadata <-  slice_4@meta.data
metadata$col <-  metadata$new_x %>% as.numeric() %>% as.integer()
metadata$row <-  metadata$new_y %>% as.numeric() %>% as.integer()

# Construct sce object
sce <- SingleCellExperiment(assays = list(counts = slice_4@assays$RNA@counts, logcounts = slice_4@assays$RNA@data),
                            colData=metadata)

# Only take highly expressed genes for spatial clustering
sce_filtered <- sce[rownames(high_exp_genes),]
sce_filtered <- spatialPreprocess(sce_filtered, platform="ST",
                              n.PCs=4,n.HVGs = 2000, log.normalize=TRUE)
# Select the number of clusters
p <- qTune(sce_filtered, qs=seq(2, 10), platform="ST")

# Clustering with BayesSpace
set.seed(149)
spatial_cluster <- spatialCluster(p, q=5, platform="ST", d=7,
                           init.method="kmeans", model="t", gamma=2,
                           nrep=10000, burn.in=100,
                           save.chain=TRUE)

# Visualizing spatial clusters
clusterPlot(spatial_cluster, color="black") +
  theme_bw() +
  xlab("") +
  ylab("") + 
  theme(text = element_text(size = 20)) 
```

<div class="figure" style="text-align: center">

<img src="1.-Palette-pipeline--Drosophila-slices-as-example-_files/figure-gfm/unnamed-chunk-3-1.png" alt="Spatial Clustering of ST data." width="40%" />
<p class="caption">
Spatial Clustering of ST data.
</p>

</div>

``` r
slice_4$spatial_cluster <- spatial_cluster$spatial.cluster
Idents(slice_4) <- slice_4$spatial_cluster

# Check the proportions of each cluster
prop_table <- prop.table(table(slice_4$spatial_cluster)) %>% t()
colnames(prop_table) <- c("cluster 1","cluster 2","cluster 3","cluster 4","cluster 5")
rownames(prop_table) <- "Proportion"
knitr::kable(prop_table,
             caption = "The proportions of each cluster.")
```

|            | cluster 1 | cluster 2 | cluster 3 | cluster 4 | cluster 5 |
|:-----------|----------:|----------:|----------:|----------:|----------:|
| Proportion | 0.3997195 | 0.1640954 | 0.1823282 | 0.1107994 | 0.1430575 |

The proportions of each cluster.

<p>
 
</p>

### Pseudo bulk deconvolution

MuSiC is employed for deconvolution. Also, only the selected 2000 highly
expressed genes are used for deconvolution. We can then achieve the
abundances of each cluster in pseudo bulk data.

``` r
# Pre-process data for MuSiC
bulk.est <- Biobase::ExpressionSet(assayData = slice_bulk %>% as.matrix())
sc.eset <- ExpressionSet(assayData = as.matrix(GetAssayData(slice_4[rownames(high_exp_genes),])), phenoData = new("AnnotatedDataFrame", slice_4@meta.data))

slice_4$sample_id <- colnames(slice_4)
sce <- SingleCellExperiment(assays = list(counts = slice_4[rownames(high_exp_genes),]@assays$RNA@counts, logcounts = slice_4[rownames(high_exp_genes),]@assays$RNA@data),
                            colData=slice_4@meta.data)

# MuSiC deconvolution
Est.pro = music_prop(bulk.mtx = exprs(bulk.est), sc.sce = sce,  clusters = "spatial_cluster", samples = "sample_id")
proportion <- Est.pro$Est.prop.weighted
proportion_reorder <- proportion[, c("1","2","3","4","5")]

# The estimated proportion of each cluster in pseudo bulk data
prop_table_2 <- proportion_reorder[2,] %>% t()
colnames(prop_table_2) <- c("cluster 1","cluster 2","cluster 3","cluster 4","cluster 5")
rownames(prop_table_2) <- "Proportion"
knitr::kable(prop_table_2,
             caption = "The estimated proportions of each cluster in pseudo bulk data (slice 5).")
```

|            | cluster 1 | cluster 2 | cluster 3 | cluster 4 | cluster 5 |
|:-----------|----------:|----------:|----------:|----------:|----------:|
| Proportion |  0.396305 | 0.2290089 | 0.2436111 | 0.0656715 | 0.0654036 |

The estimated proportions of each cluster in pseudo bulk data (slice 5).

<p>
 
</p>

### Variable factor calculation

Here we define a variable factor to adjust the expression difference
between pseudo bulk and ST data. The variable factor is defined as the
ratio of bulk expression matrix to the pseudo bulk matrix. Each gene has
its own variable factor to adjust its expression.

``` r
# The average gene expression of each cluster
cluster.averages <- AverageExpression(slice_4) %>%
  as.data.frame() %>%
  as.matrix()

cluster.averages_filtered <- cluster.averages %>%
  as.data.frame() %>%
  rownames_to_column() %>%
  dplyr::filter(.$rowname %in% rownames(slice_bulk)) %>%
  column_to_rownames(var= "rowname")  %>%
  as.matrix() %>%
  .[,c("RNA.1","RNA.2","RNA.3","RNA.4","RNA.5")]

# The abundance of each cluster in pseudo bulk 
abundance <- proportion_reorder %>%
  .[2,,drop=FALSE] %>%
  as.data.frame() %>%
  as.matrix()


slice_5_bulk_reorder <- slice_5_bulk[match(rownames(cluster.averages_filtered), rownames(slice_5_bulk)),,drop=FALSE] %>%
  as.data.frame()

# Define variable factor
variable_factor <- (slice_5_bulk_reorder/(cluster.averages_filtered %*% (abundance %>% t()))) %>%
  .[is.finite(rowSums(.)),, drop = FALSE] %>%
  as.data.frame() %>%
  set_colnames("var")

variable_factor_matrix <- rep(variable_factor,5) %>%
  as.data.frame() %>%
  set_colnames(c("1","2","3","4","5"))
rownames(variable_factor_matrix) <- rownames(variable_factor)
```

    ## [1] "The top 5 rows of variable factor"

    ##                          var
    ## 128up               497.0545
    ## 14-3-3epsilon       364.0360
    ## 14-3-3zeta          355.6674
    ## 140up               386.9003
    ## 18SrRNA-Psi:CR41602 344.8970

The adjusted matrix is the overall gene expression of each cluster,
which is obtained by taking the dot product of the cluster expression
matrix of ST slice and the variable factor.

``` r
adjusted_matrix <- ((cluster.averages_filtered[rownames(variable_factor_matrix),] * variable_factor_matrix)) 
colnames(adjusted_matrix) <- c("cluster 1","cluster 2","cluster 3","cluster 4","cluster 5")
```

    ## [1] "The top 5 rows of adjusted matrix"

    ##                     cluster 1 cluster 2   cluster 3   cluster 4 cluster 5
    ## 128up                30.25231  45.50767   32.064688   64.999345  24.00725
    ## 14-3-3epsilon       482.38149 428.18683 1022.029042 1264.198563 451.50266
    ## 14-3-3zeta          724.20988 644.44031  840.021008  854.037742 752.64577
    ## 140up                10.88988  11.16867    7.669921    4.940966  16.02148
    ## 18SrRNA-Psi:CR41602  51.66063 147.46913   52.623960  149.876118 129.43712

<p>
 
</p>

### Estimate expression of each spot

The expression of each spot is estimated through a loop algorithm. In
each interaction of the loop, the procedure begins by selecting one
random spot and its neighboring spots. The number of neighboring spots
can be changed depending on the data. The expression of spots belonging
to the same cluster is aggregated to form a pseudo bulk data called
regional cluster ST (RST). Here we assume that the ratio of RST to the
entire ST data is equal to the ratio of the adjusted RST to the adjusted
matrix. We define the ratio of RST to the entire ST as regional factor
K. The adjust regional ST equals to the dot product of adjusted matrix
and regional factor K. The expression of adjusted RST is then evenly
allocated into the selected spots of this cluster. After thousands of
interactions, the average expression of each spot is almost stable and
considered as the output estimated expression.

``` r
ST <- slice_4_ST_sum %>%
    as.data.frame() %>%
    rownames_to_column() %>%
    dplyr::select(Gene=rowname, Counts = 'slice_4_ST')


#### Nearest 8 spots
output_final <- data.frame(matrix(ncol = 0, nrow = length(rownames(variable_factor_matrix))))
rownames(output_final) <- rownames(variable_factor_matrix)
for (i in 1:5000) {
  # Select one random spot
  set.seed(i)
  slice_4.random.obj <- slice_4[, sample(colnames(slice_4), size = 1, replace=F)]
  spot <- slice_4.random.obj@meta.data[,c("new_x","new_y")] %>% as.data.frame()
  df <- slice_4@meta.data[,c("new_x","new_y")] %>% as.data.frame()
  # Calculate the distance from other spots to the selected spot
  df$distance = sqrt((df$new_x - spot$new_x)^2 + (df$new_y - spot$new_y)^2)
  # Here we select the nearest 8 spots. The number of the nearest spots can be changed here. 
  df <- df[order(df$distance),][1:9,]
  slice_4.random.neighbor <-  slice_4[,rownames(df)]
  output_region <- data.frame(matrix(ncol = 0, nrow = length(rownames(variable_factor_matrix))))
  rownames(output_region) <- rownames(variable_factor_matrix)
  # For each cluster in the selecting spots
  for (k in as.numeric(levels(as.factor(slice_4.random.neighbor$spatial_cluster)))){
    # Select the spots belonging to the same cluster
    cluster <- subset(slice_4.random.neighbor, subset = spatial_cluster == k)
    # Sum the expression of these spots from the same cluster
    RST <- colSums(t((cluster@assays$RNA@counts))) %>%
      as.data.frame() %>%
      rownames_to_column() %>%
      dplyr::select(Gene=rowname, Counts = '.')
    # Define a regional factor K as the ratio of RST to ST
    K <- RST$Counts/ST$Counts %>%
      as.data.frame()
    rownames(K) <- ST$Gene
    K_filtered <- K %>%
      rownames_to_column() %>%
      dplyr::select(Gene=rowname, factor = '.') %>%
      dplyr::filter(.$Gene %in% rownames(variable_factor_matrix)) %>%
      column_to_rownames(var= "Gene")  %>%
      as.matrix()
    # Adjusted regional bulk equals to the dot product of adjusted matrix and regional factor K
    region_bulk <- (adjusted_matrix[,k] %>% as.data.frame()) * K_filtered
    # Assign expression to each spot
    each_spot <- data.frame(region_bulk/length(cluster$orig.ident))
    output <- cbind(each_spot, rep(each_spot[1],  length(cluster$orig.ident)))
    output <- output[1:(length(output)-1)]
    colnames(output) <- rownames(cluster@meta.data)
    output_region <- output_region %>%
      cbind(output)
    }
  output_final <- output_final %>%
    cbind(output_region)
}

# Calculate the average expression of each spot
cell_list <- as.list(rownames(slice_4@meta.data))

spot_average_counts_final <- data.frame(matrix(ncol = 0, nrow = length(rownames(variable_factor_matrix))))
rownames(spot_average_counts_final) <- rownames(variable_factor_matrix)
spot_rep_times <- data.frame()

for ( i in 1:length(cell_list)) {
  spot <- colnames(output_final) == cell_list[i]
  spot_counts <- output_final[, spot, drop = FALSE]
  spot_rep_times <- spot_rep_times %>%
    rbind(length(colnames(spot_counts)))
  spot_average_counts <-  data.frame(rowSums(spot_counts))/length(spot_counts)
  colnames(spot_average_counts) <- cell_list[i]
  spot_average_counts_final <- spot_average_counts_final %>%
    cbind(spot_average_counts)
}

# Covert Palette-inferred data into seurat object
slice_4_adjusted <- CreateSeuratObject(spot_average_counts_final,assay = "RNA",meta.data = metadata)
slice_4_adjusted@images$image <- slice_4@images$image
saveRDS(slice_4_adjusted, here::here("test","slice_4_adjusted.rds"))
```

<p>
 
</p>

### Palette performance assessment

We assess the Palette performance through comparing the expression
patterns of Palette implemented data to original ST data and *in situ*
hybridization images. Here we show the expression patterns of marker
genes CG5171 (amnioserosa) and Idgf6 (hemolymph). Compared the original
ST data (Left), the Palette-implemented data (Right) significantly
increases the signal-to-noise ratio and exhibits more specific
expression patterns, which are more similar to the bona fide patterns
observed by *in situ* hybridization.

``` r
# Construct sce object for expression visualization 
metadata <-  slice_4_adjusted@meta.data
metadata$col <-  slice_4_adjusted@meta.data$col
metadata$row <-  slice_4_adjusted@meta.data$row


sce_S4 <- SingleCellExperiment(assays = list(counts = slice_4_adjusted@assays$RNA@counts, logcounts = slice_4_adjusted@assays$RNA@data),
                            colData= metadata)

sce_S4 <- spatialPreprocess(sce_S4, platform="ST",
                              n.PCs=4,n.HVGs = 2000, log.normalize=FALSE,skip.PCA=TRUE)


p1 <- featurePlot(
  sce_S4,
  "CG5171",
  assay.type = "logcounts",
  diverging = FALSE,
  low = NULL,
  high = "darkblue",
  mid = NULL,
  color = NULL,
  platform = NULL,
  is.enhanced = NULL,
)

p2 <- featurePlot(
  sce_S4,
  "Idgf6",
  assay.type = "logcounts",
  diverging = FALSE,
  low = NULL,
  high = "darkblue",
  mid = NULL,
  color = NULL,
  platform = NULL,
  is.enhanced = NULL,
)

slice_5_metadata <-  slice_5@meta.data
slice_5_metadata$col <-  slice_5_metadata$new_x %>% as.numeric() %>% as.integer()
slice_5_metadata$row <-  slice_5_metadata$new_y %>% as.numeric() %>% as.integer()


sce_S5 <- SingleCellExperiment(assays = list(counts = slice_5@assays$RNA@counts, logcounts = slice_5@assays$RNA@data),
                            colData=slice_5_metadata)

sce_S5 <- spatialPreprocess(sce_S5, platform="ST",
                              n.PCs=4,n.HVGs = 2000, log.normalize=FALSE,skip.PCA=TRUE)
markers <- FindAllMarkers(slice_5)
p3 <- featurePlot(
  sce_S5,
  "CG5171",
  assay.type = "logcounts",
  diverging = FALSE,
  low = NULL,
  high = "darkblue",
  mid = NULL,
  color = NULL,
  platform = NULL,
  is.enhanced = NULL,
)

p4 <- featurePlot(
  sce_S5,
  "Idgf6",
  assay.type = "logcounts",
  diverging = FALSE,
  low = NULL,
  high = "darkblue",
  mid = NULL,
  color = NULL,
  platform = NULL,
  is.enhanced = NULL,
)

ggarrange(p3, p1, p4, p2, ncol = 2, nrow = 2)
```

<div class="figure" style="text-align: center">

<img src="1.-Palette-pipeline--Drosophila-slices-as-example-_files/figure-gfm/unnamed-chunk-4-1.png" alt="Comparsion of expression patterns." width="70%" />
<p class="caption">
Comparsion of expression patterns.
</p>

</div>

<p>
 
</p>

### Session information

``` r
devtools::session_info()
```

    ## ─ Session info ───────────────────────────────────────────────────────────────
    ##  setting  value
    ##  version  R version 4.3.1 (2023-06-16)
    ##  os       macOS Sonoma 14.2
    ##  system   aarch64, darwin20
    ##  ui       X11
    ##  language (EN)
    ##  collate  en_US.UTF-8
    ##  ctype    en_US.UTF-8
    ##  tz       Asia/Shanghai
    ##  date     2024-03-27
    ##  pandoc   3.1.1 @ /Applications/RStudio.app/Contents/Resources/app/quarto/bin/tools/ (via rmarkdown)
    ## 
    ## ─ Packages ───────────────────────────────────────────────────────────────────
    ##  package              * version   date (UTC) lib source
    ##  abind                  1.4-5     2016-07-21 [1] CRAN (R 4.3.0)
    ##  assertthat             0.2.1     2019-03-21 [1] CRAN (R 4.3.0)
    ##  backports              1.4.1     2021-12-13 [1] CRAN (R 4.3.0)
    ##  BayesSpace           * 1.10.1    2023-06-11 [1] Bioconductor
    ##  beachmat               2.16.0    2023-05-08 [1] Bioconductor
    ##  beeswarm               0.4.0     2021-06-01 [1] CRAN (R 4.3.0)
    ##  Biobase              * 2.60.0    2023-05-08 [1] Bioconductor
    ##  BiocFileCache          2.8.0     2023-05-08 [1] Bioconductor
    ##  BiocGenerics         * 0.46.0    2023-06-04 [1] Bioconductor
    ##  BiocNeighbors          1.18.0    2023-05-08 [1] Bioconductor
    ##  BiocParallel           1.34.2    2023-05-28 [1] Bioconductor
    ##  BiocSingular           1.16.0    2023-05-08 [1] Bioconductor
    ##  bit                    4.0.5     2022-11-15 [1] CRAN (R 4.3.0)
    ##  bit64                  4.0.5     2020-08-30 [1] CRAN (R 4.3.0)
    ##  bitops                 1.0-7     2021-04-24 [1] CRAN (R 4.3.0)
    ##  blob                   1.2.4     2023-03-17 [1] CRAN (R 4.3.0)
    ##  bluster                1.10.0    2023-05-08 [1] Bioconductor
    ##  broom                  1.0.5     2023-06-09 [1] CRAN (R 4.3.0)
    ##  cachem                 1.0.8     2023-05-01 [1] CRAN (R 4.3.0)
    ##  callr                  3.7.3     2022-11-02 [1] CRAN (R 4.3.0)
    ##  car                    3.1-2     2023-03-30 [1] CRAN (R 4.3.0)
    ##  carData                3.0-5     2022-01-06 [1] CRAN (R 4.3.0)
    ##  class                  7.3-22    2023-05-03 [1] CRAN (R 4.3.1)
    ##  cli                    3.6.1     2023-03-23 [1] CRAN (R 4.3.0)
    ##  cluster                2.1.4     2022-08-22 [1] CRAN (R 4.3.1)
    ##  coda                   0.19-4    2020-09-30 [1] CRAN (R 4.3.0)
    ##  codetools              0.2-19    2023-02-01 [1] CRAN (R 4.3.1)
    ##  colorspace             2.1-0     2023-01-23 [1] CRAN (R 4.3.0)
    ##  corpcor                1.6.10    2021-09-16 [1] CRAN (R 4.3.0)
    ##  cowplot                1.1.1     2020-12-30 [1] CRAN (R 4.3.0)
    ##  crayon                 1.5.2     2022-09-29 [1] CRAN (R 4.3.0)
    ##  curl                   5.1.0     2023-10-02 [1] CRAN (R 4.3.1)
    ##  data.table             1.14.8    2023-02-17 [1] CRAN (R 4.3.0)
    ##  DBI                    1.1.3     2022-06-18 [1] CRAN (R 4.3.0)
    ##  dbplyr                 2.3.4     2023-09-26 [1] CRAN (R 4.3.1)
    ##  DelayedArray           0.26.7    2023-07-30 [1] Bioconductor
    ##  DelayedMatrixStats     1.22.6    2023-09-03 [1] Bioconductor
    ##  deldir                 1.0-9     2023-05-17 [1] CRAN (R 4.3.0)
    ##  devtools               2.4.5     2022-10-11 [1] CRAN (R 4.3.0)
    ##  digest                 0.6.33    2023-07-07 [1] CRAN (R 4.3.0)
    ##  DirichletReg           0.7-1     2021-05-18 [1] CRAN (R 4.3.0)
    ##  doParallel             1.0.17    2022-02-07 [1] CRAN (R 4.3.0)
    ##  dotCall64              1.1-0     2023-10-17 [1] CRAN (R 4.3.1)
    ##  dplyr                * 1.1.4     2023-11-17 [1] CRAN (R 4.3.1)
    ##  dqrng                  0.3.1     2023-08-30 [1] CRAN (R 4.3.0)
    ##  e1071                  1.7-13    2023-02-01 [1] CRAN (R 4.3.0)
    ##  edgeR                  3.42.4    2023-06-04 [1] Bioconductor
    ##  ellipsis               0.3.2     2021-04-29 [1] CRAN (R 4.3.0)
    ##  EpiDISH              * 2.16.0    2023-05-08 [1] Bioconductor
    ##  evaluate               0.22      2023-09-29 [1] CRAN (R 4.3.1)
    ##  fansi                  1.0.5     2023-10-08 [1] CRAN (R 4.3.1)
    ##  farver                 2.1.1     2022-07-06 [1] CRAN (R 4.3.0)
    ##  fastmap                1.1.1     2023-02-24 [1] CRAN (R 4.3.0)
    ##  filelock               1.0.2     2018-10-05 [1] CRAN (R 4.3.0)
    ##  fitdistrplus           1.1-11    2023-04-25 [1] CRAN (R 4.3.0)
    ##  forcats              * 1.0.0     2023-01-29 [1] CRAN (R 4.3.0)
    ##  foreach                1.5.2     2022-02-02 [1] CRAN (R 4.3.0)
    ##  Formula                1.2-5     2023-02-24 [1] CRAN (R 4.3.0)
    ##  fs                     1.6.3     2023-07-20 [1] CRAN (R 4.3.0)
    ##  future                 1.33.0    2023-07-01 [1] CRAN (R 4.3.0)
    ##  future.apply           1.11.0    2023-05-21 [1] CRAN (R 4.3.0)
    ##  generics               0.1.3     2022-07-05 [1] CRAN (R 4.3.0)
    ##  GenomeInfoDb         * 1.36.4    2023-10-08 [1] Bioconductor
    ##  GenomeInfoDbData       1.2.10    2023-06-27 [1] Bioconductor
    ##  GenomicRanges        * 1.52.1    2023-10-08 [1] Bioconductor
    ##  GGally                 2.1.2     2021-06-21 [1] CRAN (R 4.3.0)
    ##  ggbeeswarm             0.7.2     2023-04-29 [1] CRAN (R 4.3.0)
    ##  ggplot2              * 3.4.4     2023-10-12 [1] CRAN (R 4.3.1)
    ##  ggpubr               * 0.6.0     2023-02-10 [1] CRAN (R 4.3.0)
    ##  ggrepel                0.9.4     2023-10-13 [1] CRAN (R 4.3.1)
    ##  ggridges               0.5.4     2022-09-26 [1] CRAN (R 4.3.0)
    ##  ggsignif               0.6.4     2022-10-13 [1] CRAN (R 4.3.0)
    ##  globals                0.16.2    2022-11-21 [1] CRAN (R 4.3.0)
    ##  glue                   1.6.2     2022-02-24 [1] CRAN (R 4.3.0)
    ##  goftest                1.2-3     2021-10-07 [1] CRAN (R 4.3.0)
    ##  gridExtra              2.3       2017-09-09 [1] CRAN (R 4.3.0)
    ##  grr                    0.9.5     2016-08-26 [1] CRAN (R 4.3.0)
    ##  gtable                 0.3.4     2023-08-21 [1] CRAN (R 4.3.0)
    ##  here                   1.0.1     2020-12-13 [1] CRAN (R 4.3.0)
    ##  hms                    1.1.3     2023-03-21 [1] CRAN (R 4.3.0)
    ##  htmltools              0.5.6.1   2023-10-06 [1] CRAN (R 4.3.1)
    ##  htmlwidgets            1.6.2     2023-03-17 [1] CRAN (R 4.3.0)
    ##  httpuv                 1.6.11    2023-05-11 [1] CRAN (R 4.3.0)
    ##  httr                   1.4.7     2023-08-15 [1] CRAN (R 4.3.0)
    ##  ica                    1.0-3     2022-07-08 [1] CRAN (R 4.3.0)
    ##  igraph                 1.5.1     2023-08-10 [1] CRAN (R 4.3.0)
    ##  IRanges              * 2.34.1    2023-06-22 [1] Bioconductor
    ##  irlba                  2.3.5.1   2022-10-03 [1] CRAN (R 4.3.1)
    ##  iterators              1.0.14    2022-02-05 [1] CRAN (R 4.3.0)
    ##  jsonlite               1.8.7     2023-06-29 [1] CRAN (R 4.3.0)
    ##  KernSmooth             2.23-22   2023-07-10 [1] CRAN (R 4.3.0)
    ##  knitr                  1.44      2023-09-11 [1] CRAN (R 4.3.0)
    ##  labeling               0.4.3     2023-08-29 [1] CRAN (R 4.3.0)
    ##  later                  1.3.1     2023-05-02 [1] CRAN (R 4.3.0)
    ##  lattice                0.21-9    2023-10-01 [1] CRAN (R 4.3.1)
    ##  lazyeval               0.2.2     2019-03-15 [1] CRAN (R 4.3.0)
    ##  leiden                 0.4.3     2022-09-10 [1] CRAN (R 4.3.0)
    ##  lifecycle              1.0.3     2022-10-07 [1] CRAN (R 4.3.0)
    ##  limma                * 3.56.2    2023-06-04 [1] Bioconductor
    ##  listenv                0.9.0     2022-12-16 [1] CRAN (R 4.3.0)
    ##  lmtest                 0.9-40    2022-03-21 [1] CRAN (R 4.3.0)
    ##  locfdr                 1.1-8     2015-07-15 [1] CRAN (R 4.3.0)
    ##  locfit                 1.5-9.8   2023-06-11 [1] CRAN (R 4.3.0)
    ##  lubridate            * 1.9.3     2023-09-27 [1] CRAN (R 4.3.1)
    ##  magrittr             * 2.0.3     2022-03-30 [1] CRAN (R 4.3.0)
    ##  MASS                   7.3-60    2023-05-04 [1] CRAN (R 4.3.1)
    ##  Matrix               * 1.7-0     2024-01-17 [1] R-Forge (R 4.3.1)
    ##  Matrix.utils         * 0.9.8     2020-02-26 [1] CRAN (R 4.3.1)
    ##  MatrixGenerics       * 1.12.3    2023-07-30 [1] Bioconductor
    ##  MatrixModels           0.5-2     2023-07-10 [1] CRAN (R 4.3.0)
    ##  matrixStats          * 1.0.0     2023-06-02 [1] CRAN (R 4.3.0)
    ##  maxLik                 1.5-2     2021-07-26 [1] CRAN (R 4.3.0)
    ##  mclust                 6.0.0     2022-10-31 [1] CRAN (R 4.3.0)
    ##  mcmc                   0.9-7     2020-03-21 [1] CRAN (R 4.3.0)
    ##  MCMCpack               1.6-3     2022-04-13 [1] CRAN (R 4.3.0)
    ##  memoise                2.0.1     2021-11-26 [1] CRAN (R 4.3.0)
    ##  metapod                1.8.0     2023-04-25 [1] Bioconductor
    ##  mime                   0.12      2021-09-28 [1] CRAN (R 4.3.0)
    ##  miniUI                 0.1.1.1   2018-05-18 [1] CRAN (R 4.3.0)
    ##  miscTools              0.6-28    2023-05-03 [1] CRAN (R 4.3.0)
    ##  munsell                0.5.0     2018-06-12 [1] CRAN (R 4.3.0)
    ##  MuSiC                * 1.0.0     2023-06-28 [1] Github (xuranw/MuSiC@0a3e3af)
    ##  nlme                   3.1-163   2023-08-09 [1] CRAN (R 4.3.0)
    ##  nnls                 * 1.5       2023-09-11 [1] CRAN (R 4.3.0)
    ##  pander               * 0.6.5     2022-03-18 [1] CRAN (R 4.3.0)
    ##  parallelly             1.36.0    2023-05-26 [1] CRAN (R 4.3.0)
    ##  patchwork              1.1.3     2023-08-14 [1] CRAN (R 4.3.0)
    ##  pbapply                1.7-2     2023-06-27 [1] CRAN (R 4.3.0)
    ##  pillar                 1.9.0     2023-03-22 [1] CRAN (R 4.3.0)
    ##  pkgbuild               1.4.2     2023-06-26 [1] CRAN (R 4.3.0)
    ##  pkgconfig              2.0.3     2019-09-22 [1] CRAN (R 4.3.0)
    ##  pkgload                1.3.3     2023-09-22 [1] CRAN (R 4.3.1)
    ##  plotly                 4.10.3    2023-10-21 [1] CRAN (R 4.3.1)
    ##  plyr                   1.8.9     2023-10-02 [1] CRAN (R 4.3.1)
    ##  png                    0.1-8     2022-11-29 [1] CRAN (R 4.3.0)
    ##  polyclip               1.10-6    2023-09-27 [1] CRAN (R 4.3.1)
    ##  prettyunits            1.2.0     2023-09-24 [1] CRAN (R 4.3.1)
    ##  processx               3.8.2     2023-06-30 [1] CRAN (R 4.3.0)
    ##  profvis                0.3.8     2023-05-02 [1] CRAN (R 4.3.0)
    ##  progressr              0.14.0    2023-08-10 [1] CRAN (R 4.3.0)
    ##  promises               1.2.1     2023-08-10 [1] CRAN (R 4.3.0)
    ##  proxy                  0.4-27    2022-06-09 [1] CRAN (R 4.3.0)
    ##  ps                     1.7.5     2023-04-18 [1] CRAN (R 4.3.0)
    ##  purrr                * 1.0.2     2023-08-10 [1] CRAN (R 4.3.0)
    ##  quadprog             * 1.5-8     2019-11-20 [1] CRAN (R 4.3.0)
    ##  quantreg               5.97      2023-08-19 [1] CRAN (R 4.3.0)
    ##  R6                     2.5.1     2021-08-19 [1] CRAN (R 4.3.0)
    ##  RANN                   2.6.1     2019-01-08 [1] CRAN (R 4.3.0)
    ##  RColorBrewer           1.1-3     2022-04-03 [1] CRAN (R 4.3.0)
    ##  Rcpp                   1.0.11    2023-07-06 [1] CRAN (R 4.3.0)
    ##  RcppAnnoy              0.0.21    2023-07-02 [1] CRAN (R 4.3.0)
    ##  RCurl                  1.98-1.12 2023-03-27 [1] CRAN (R 4.3.0)
    ##  readr                * 2.1.4     2023-02-10 [1] CRAN (R 4.3.0)
    ##  remotes                2.4.2.1   2023-07-18 [1] CRAN (R 4.3.0)
    ##  reshape                0.8.9     2022-04-12 [1] CRAN (R 4.3.0)
    ##  reshape2             * 1.4.4     2020-04-09 [1] CRAN (R 4.3.0)
    ##  reticulate             1.34.0    2023-10-12 [1] CRAN (R 4.3.1)
    ##  rhdf5                  2.44.0    2023-05-08 [1] Bioconductor
    ##  rhdf5filters           1.12.1    2023-05-08 [1] Bioconductor
    ##  Rhdf5lib               1.22.1    2023-09-10 [1] Bioconductor
    ##  rlang                  1.1.1     2023-04-28 [1] CRAN (R 4.3.0)
    ##  rmarkdown              2.25      2023-09-18 [1] CRAN (R 4.3.1)
    ##  ROCR                   1.0-11    2020-05-02 [1] CRAN (R 4.3.0)
    ##  rprojroot              2.0.3     2022-04-02 [1] CRAN (R 4.3.0)
    ##  RSQLite                2.3.1     2023-04-03 [1] CRAN (R 4.3.0)
    ##  rstatix                0.7.2     2023-02-01 [1] CRAN (R 4.3.0)
    ##  rstudioapi             0.15.0    2023-07-07 [1] CRAN (R 4.3.0)
    ##  rsvd                   1.0.5     2021-04-16 [1] CRAN (R 4.3.0)
    ##  Rtsne                  0.16      2022-04-17 [1] CRAN (R 4.3.0)
    ##  S4Arrays               1.2.0     2023-10-26 [1] Bioconductor
    ##  S4Vectors            * 0.38.2    2023-09-24 [1] Bioconductor
    ##  sandwich               3.0-2     2022-06-15 [1] CRAN (R 4.3.0)
    ##  ScaledMatrix           1.8.1     2023-05-08 [1] Bioconductor
    ##  scales                 1.2.1     2022-08-20 [1] CRAN (R 4.3.0)
    ##  scater                 1.28.0    2023-04-25 [1] Bioconductor
    ##  scattermore            1.2       2023-06-12 [1] CRAN (R 4.3.0)
    ##  scran                  1.28.2    2023-07-23 [1] Bioconductor
    ##  sctransform            0.4.1     2023-10-19 [1] CRAN (R 4.3.1)
    ##  scuttle                1.10.3    2023-10-15 [1] Bioconductor
    ##  sessioninfo            1.2.2     2021-12-06 [1] CRAN (R 4.3.0)
    ##  Seurat               * 4.4.0     2023-09-28 [1] CRAN (R 4.3.1)
    ##  SeuratObject         * 5.0.1     2023-11-17 [1] CRAN (R 4.3.1)
    ##  shiny                  1.7.5.1   2023-10-14 [1] CRAN (R 4.3.1)
    ##  SingleCellExperiment * 1.22.0    2023-05-08 [1] Bioconductor
    ##  sp                     2.1-1     2023-10-16 [1] CRAN (R 4.3.1)
    ##  spam                   2.9-1     2022-08-07 [1] CRAN (R 4.3.0)
    ##  SparseM                1.81      2021-02-18 [1] CRAN (R 4.3.0)
    ##  sparseMatrixStats      1.12.2    2023-07-02 [1] Bioconductor
    ##  spatstat.data          3.0-1     2023-03-12 [1] CRAN (R 4.3.0)
    ##  spatstat.explore       3.2-3     2023-09-07 [1] CRAN (R 4.3.0)
    ##  spatstat.geom          3.2-7     2023-10-20 [1] CRAN (R 4.3.1)
    ##  spatstat.random        3.2-1     2023-10-21 [1] CRAN (R 4.3.1)
    ##  spatstat.sparse        3.0-2     2023-06-25 [1] CRAN (R 4.3.0)
    ##  spatstat.utils         3.0-3     2023-05-09 [1] CRAN (R 4.3.0)
    ##  statmod                1.5.0     2023-01-06 [1] CRAN (R 4.3.0)
    ##  stringi                1.7.12    2023-01-11 [1] CRAN (R 4.3.0)
    ##  stringr              * 1.5.0     2022-12-02 [1] CRAN (R 4.3.0)
    ##  SummarizedExperiment * 1.30.2    2023-06-11 [1] Bioconductor
    ##  survival               3.5-7     2023-08-14 [1] CRAN (R 4.3.0)
    ##  tensor                 1.5       2012-05-05 [1] CRAN (R 4.3.0)
    ##  tibble               * 3.2.1     2023-03-20 [1] CRAN (R 4.3.0)
    ##  tidyr                * 1.3.0     2023-01-24 [1] CRAN (R 4.3.0)
    ##  tidyselect             1.2.0     2022-10-10 [1] CRAN (R 4.3.0)
    ##  tidyverse            * 2.0.0     2023-02-22 [1] CRAN (R 4.3.0)
    ##  timechange             0.2.0     2023-01-11 [1] CRAN (R 4.3.0)
    ##  TOAST                * 1.14.0    2023-05-08 [1] Bioconductor
    ##  tzdb                   0.4.0     2023-05-12 [1] CRAN (R 4.3.0)
    ##  urlchecker             1.0.1     2021-11-30 [1] CRAN (R 4.3.0)
    ##  usethis                2.2.2     2023-07-06 [1] CRAN (R 4.3.0)
    ##  utf8                   1.2.3     2023-01-31 [1] CRAN (R 4.3.0)
    ##  uwot                   0.1.16    2023-06-29 [1] CRAN (R 4.3.0)
    ##  vctrs                  0.6.4     2023-10-12 [1] CRAN (R 4.3.1)
    ##  vipor                  0.4.5     2017-03-22 [1] CRAN (R 4.3.0)
    ##  viridis                0.6.4     2023-07-22 [1] CRAN (R 4.3.0)
    ##  viridisLite            0.4.2     2023-05-02 [1] CRAN (R 4.3.0)
    ##  withr                  2.5.1     2023-09-26 [1] CRAN (R 4.3.1)
    ##  xfun                   0.40      2023-08-09 [1] CRAN (R 4.3.0)
    ##  xgboost                1.7.5.1   2023-03-30 [1] CRAN (R 4.3.0)
    ##  xtable                 1.8-4     2019-04-21 [1] CRAN (R 4.3.0)
    ##  XVector                0.40.0    2023-05-08 [1] Bioconductor
    ##  yaml                   2.3.7     2023-01-23 [1] CRAN (R 4.3.0)
    ##  zlibbioc               1.46.0    2023-05-08 [1] Bioconductor
    ##  zoo                    1.8-12    2023-04-13 [1] CRAN (R 4.3.0)
    ## 
    ##  [1] /Library/Frameworks/R.framework/Versions/4.3-arm64/Resources/library
    ## 
    ## ──────────────────────────────────────────────────────────────────────────────
