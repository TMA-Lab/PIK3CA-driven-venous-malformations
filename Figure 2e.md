Fig.2e: GO enrichment analysis using top Vein and Artery cluster markers
================
Marle Kraft
2023-10-02

- [Packages](#packages)
- [Data frame generation](#data-frame-generation)
- [Cluster Vein/Artery matrix and heatmap
  generation](#cluster-veinartery-matrix-and-heatmap-generation)

## Packages

``` r
library(pheatmap)
```

## Data frame generation

``` r
#This code creates a heatmap of selected GO terms, where each cell represents the P-value for a specific GO term in the clusters "Vein" and "arterial"Artery"

# Create a matrix with P-values 
# note 1: As input the log10(P-value) was used at data input and the scale labels were changed retrospectively
# note 2: some GO terms were abbriviated to meet figure size restrictions
# note 3: Grouping of GO terms in the figure was performed manually
# note 4: data can be found in Fig.2 source data files


Artery_and_Vein <- data.frame(
  GO_Term = c("reg. of vasoconstriction",
              "vasoconstriction",
              "pos. reg. of vasoconstriction",
              "reg. of tube size",
              "reg. of tube diameter",
              "blood vessel diameter maintenance",
              "actin filament bundle assembly",
              "endothelium development",
              "blood vessel morphogenesis",
              "blood vessel development",
              "vasculature development",
              "blood vessel remodeling",
              "artery development",
              "artery morphogenesis",
              "regulation of muscle cell differentiation",
              "regulation of cell-substrate adhesion",
              "regulation of cell-matrix adhesion",
              "pos. reg. of cell-substrate adhesion",
              "regulation of blood vessel EC migration",
              "integrin-mediated signaling pathway",
              "pos. reg. of EC migration",
              "reg. of EC migration",
              "EC migration",
              "blood vessel EC migration",
              "pos. Reg. of leukocyte migration",
              "reg. of leukocyte migration",
              "reg. of leukocyte migration",
              "heterophilic cell-cell adhesion",
              "leukocyte adhesion to vascular EC",
              "cellular extravasation",
              "response to toxic substance",
              "response to hydrogen peroxide",
              "response to ROS",
              "neg. reg. of extrinsic apoptotic signaling pathway",
              "reg. of extrinsic apoptotic signaling pathway"
  ),
  Vein = c(3.40811240508341,3.02024705904035,0,3.17554475442736,3.18684499824345,3.18684499824345,2.31069951007139,7.07646145051875,14.7818663064868,17.6335783867311,18.1814050174235,5.56971030691044,7.66582212995808,7.45025964099538,5.48922034600584,4.71323526586547,4.8990702946901,6.5166948934222,4.40302392441232,4.21472116555849,4.09879360665064,7.42641631517567,7.45646486057162,6.22364302963367,0,0,0,0,0,0,0,0,0,0,0),
  Artery = c(5.059816769,4.559416001,4.586248537,4.669351539,4.683475819,4.683475819,3.634684329,2.98289072904831,4.26423815371653,3.60371944451099,3.76869803166972,0,0,0,0,0,0,0,0,0,0,0,0,0,7.300827615,6.039349703,4.832713101,4.687364667,4.586248537,3.861583822,3.182540748,4.324643439,4.454754013,3.526106488, 3.71889433243673)
)
```

## Cluster Vein/Artery matrix and heatmap generation

``` r
matrix_data <- as.matrix(Artery_and_Vein[, c("Artery", "Vein")])

# Create a row names vector for the heatmap
rownames(matrix_data) <- Artery_and_Vein$GO_Term

# Create the heatmap
p <- pheatmap(
  matrix_data,
  color = colorRampPalette(c("white", "red"))(50),
  main = "GO terms of Vein/Artery signature",
  fontsize = 6, 
  cellwidth = 15, 
  cellheight = 7,  
  display_numbers = FALSE, 
  show_colnames = TRUE,  
  show_rownames = TRUE,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  labels_col = c("Vein", "Artery")
)
```

![](Figure-2e_files/figure-gfm/unnamed-chunk-2-1.png)

# Save cluster Artey/Vein heatmap as a PDF

``` r
pdf("GO terms of Artery and Vein.pdf")
print(p)  
dev.off()  
```
