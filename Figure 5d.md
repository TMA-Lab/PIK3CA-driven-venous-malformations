Fig.5d: GO enrichment analysis of upregulated DEGs unique for Pik3ca-1
or Pik3ca-2, and those shared with FOXO1 signature
================
Marle Kraft
2023-11-20

- [Packages](#packages)
- [Cluster Pik3ca-1 data frame
  generation](#cluster-pik3ca-1-data-frame-generation)
- [Cluster Pik3ca-1 matrix and heatmap
  generation](#cluster-pik3ca-1-matrix-and-heatmap-generation)
- [Cluster Pik3ca-2 data frame
  generation](#cluster-pik3ca-2-data-frame-generation)
- [Cluster Pik3ca-2 matrix and heatmap
  generation](#cluster-pik3ca-2-matrix-and-heatmap-generation)

## Packages

``` r
library(pheatmap)
```

## Cluster Pik3ca-1 data frame generation

``` r
# Create a matrix with P-values for cluster Pik3ca-1
# note 1: As input the log10(P-value) was used at data input and the scale labels were changed retrospectively
# note 2: some GO terms were abbriviated to meet figure size restrictions
# note 3: data can be found in Fig.5 source data files


Pik3ca_1 <- data.frame(
  GO_Term = c(
    "lysosome organization",
    "vacuole organization",
    "ER organization",
    "cellular response to oxygen levels",
    "glycosaminoglycan metabolic process",
    "reactive nitrogen species metabolic process",
    "negative regulation of protein catabolic process",
    "cellular lipid catabolic process",
    "fatty acid catabolic process",
    "monocarboxylic acid catabolic process",
    "carbohydrate derivative catabolic process",
    "ECM organization",
    "carbohydrate catabolic process",
    "carbohydrate catabolic process",
    "monosaccharide metabolic process",
    "ATP metabolic process",
    "glycolytic process",
    "ATP generation from ADP",
    "purine NTP metabolic process",
    "NDP phosphorylation",
    "nucleotide phosphorylation",
    "pyruvate metabolic process",
    "ribonucleotide biosynthetic process",
    "glucose metabolic process"
  ),
  
  Pik3ca_1_Unique = c(5.314145792,
                      4.863088441,
                      4.264238154,
                      4.324614466,
                      3.78372493,
                      5.030135119,
                      4.411711966,
                      5.824711287,
                      5.363960776,
                      4.398365687,
                      4.814680359,
                      5.211589017,
                      6.994311111,
                      6.994311111,
                      4.279479097,
                      2.237505632,
                      2.903091336,
                      2.872655906,
                      2.175055083,
                      2.438573169,
                      2.343522733,
                      0,
                      0,
                      0),
  Pik3ca_1_Shared_with_FOXO1 = c(0,
                                 0,
                                 0,
                                 0,
                                 0,
                                 0,
                                 0,
                                 0,
                                 0,
                                 2.236823601,
                                 3.294026221,
                                 9.979707641,
                                 9.935716,
                                 9.935716,
                                 7.68624374,
                                 6.918570151,
                                 8.011508588,
                                 7.965477109,
                                 7.929240042,
                                 8.409739816,
                                 8.242841508,
                                 9.03440919,
                                 6.313978907,
                                 6.206231307)
)
```

## Cluster Pik3ca-1 matrix and heatmap generation

``` r
matrix_data <- as.matrix(Pik3ca_1[, c("Pik3ca_1_Unique", "Pik3ca_1_Shared_with_FOXO1")])

# Create a row names vector for the heatmap
rownames(matrix_data) <- Pik3ca_1$GO_Term

# Create the heatmap
p <- pheatmap(
  matrix_data,
  color = colorRampPalette(c("white", "red"))(50),
  main = "GO terms of PI3K signature",
  fontsize = 6, 
  cellwidth = 15, 
  cellheight = 7,  
  display_numbers = FALSE, 
  show_colnames = TRUE,  
  show_rownames = TRUE,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  labels_col = c("Unique", "Shared with FOXO1")
)
```

![](Fig.5d_files/figure-gfm/unnamed-chunk-2-1.png)<!-- --> \# Save
cluster Pik3ca-1 heatmap as a PDF

``` r
pdf("GO terms of PI3K_Pik3ca-1.pdf")
print(p)  
dev.off()  
```

    ## quartz_off_screen 
    ##                 2

## Cluster Pik3ca-2 data frame generation

``` r
# Create a matrix with P-values for cluster Pik3ca-2
# note 1: As input the log10(P-value) was used at data input and the scale labels were changed retrospectively
# note 2: some GO terms were abbriviated to meet figure size restrictions
# note 3: data can be found in Fig.5 source data files

Pik3ca_2 <- data.frame(
  GO_Term = c("pos. regulation of GTPase activity",
              "regulation of exocytosis",
              "actin filament bundle assembly",
              "actin filament bundle organization",
              "reg. of cell junction assembly",
              "EC migration",
              "cell-substrate adhesion",
              "tube morphogenesis",
              "blood vessel development",
              "angiogenesis",
              "blood vessel morphogenesis",
              "ECM organization",
              "extracellular structure organization",
              "external encapsulating structure organization",
              "regulation of cell adhesion",
              "regulation of cell migration",
              "regulation of angiogenesis",
              "carbohydrate metabolic process",
              "mononuclear cell differentiation",
              "nucleotide metabolic process",
              "lymphocyte differentiation",
              "purine nucleoside triphosphate metabolic process",
              "leukocyte differentiation",
              "leukocyte cell-cell adhesion",
              "regulation of leukocyte differentiation",
              "regulation of T cell activation",
              "purine ribonucleotide metabolic process",
              "cellular carbohydrate metabolic process",
              "T cell activation"
              ),
  
  Pik3ca_2_Unique = c(4.54846623,4.023897929, 4.201628236, 4.138362232,5.296642616,5.768456284,5.910855063,7.736673199,5.836149506,6.2502825,5.921519499,6.805741493,6.78575173,6.765841003,5.324664118,7.385690803,2.242725014,2.868092638,0,0,0,0,0,0,0,0,0,0,0),
  Pik3ca_2_Shared_with_FOXO1 = c(0,0,0,0,2.087076424,4.279117064,3.610777383,15.57724815,16.51542274,14.18611281,15.21265324,11.43829521,11.41687143,11.39552214,8.264132614,8.8392506,5.023171268,5.93129934,7.436019425,7.763027114,6.365607815,5.71725927,6.675887904,3.965747559,4.557138561,4.265697718,6.904425353,5.097551873,5.217015829)
)
```

## Cluster Pik3ca-2 matrix and heatmap generation

``` r
matrix_data <- as.matrix(Pik3ca_2[, c("Pik3ca_2_Unique", "Pik3ca_2_Shared_with_FOXO1")])

# Create a row names vector for the heatmap
rownames(matrix_data) <- Pik3ca_2$GO_Term

# Create the heatmap
p <- pheatmap(
  matrix_data,
  color = colorRampPalette(c("white", "red"))(50),
  main = "GO terms of PI3K signature",
  fontsize = 6, 
  cellwidth = 15, 
  cellheight = 7,  
  display_numbers = FALSE, 
  show_colnames = TRUE,  
  show_rownames = TRUE,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  labels_col = c("Unique", "Shared with FOXO1")
)
```

![](Fig.5d_files/figure-gfm/unnamed-chunk-5-1.png)<!-- --> \# Save
cluster Pik3ca-2 heatmap as a PDF

``` r
pdf("GO terms of PI3K_Pik3ca-2.pdf")
print(p)  
dev.off()  
```

    ## quartz_off_screen 
    ##                 2
