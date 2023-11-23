Fig.3d and 4b: Heatmap and Expression profile of zonation marker
expression
================
Marle Kraft
2023-10-02

- [Packages](#packages)
- [Load dataset](#load-dataset)
- [WT](#wt)
  - [Subset dataset for WT
    trajectory](#subset-dataset-for-wt-trajectory)
  - [Calculate and draw trajectory for WT
    dataset](#calculate-and-draw-trajectory-for-wt-dataset)
  - [Calculate gene importance (gimp) construct trajectory heatmap for
    WT dataset (Fig.
    3d)](#calculate-gene-importance-gimp-construct-trajectory-heatmap-for-wt-dataset-fig-3d)
  - [Heatmap of Sox17 and Lrg1 for WT dataset
    (Fig.4b)](#heatmap-of-sox17-and-lrg1-for-wt-dataset-fig4b)
  - [Save heatmap source data WT dataset
    (Fig.4b)](#save-heatmap-source-data-wt-dataset-fig4b)
  - [Expression profile of Sox17 and Lrg1 for WT dataset
    (Fig.4b)](#expression-profile-of-sox17-and-lrg1-for-wt-dataset-fig4b)
- [Mutant](#mutant)
  - [Subset dataset for mutant
    trajectory](#subset-dataset-for-mutant-trajectory)
  - [Calculate and draw trajectory for mutant
    dataset](#calculate-and-draw-trajectory-for-mutant-dataset)
  - [Calculate gene importance (gimp) construct trajectory heatmap for
    mutant dataset (Fig.
    3d)](#calculate-gene-importance-gimp-construct-trajectory-heatmap-for-mutant-dataset-fig-3d)
  - [Heatmap of Sox17 and Lrg1 for mutant dataset
    (Fig.4b)](#heatmap-of-sox17-and-lrg1-for-mutant-dataset-fig4b)
  - [Save heatmap source data mutant dataset
    (Fig.4b)](#save-heatmap-source-data-mutant-dataset-fig4b)
  - [Expression profile of Sox17 and Lrg1 for WT dataset
    (Fig.4b)](#expression-profile-of-sox17-and-lrg1-for-wt-dataset-fig4b-1)

## Packages

``` r
library(Seurat)
library(ggsci)
library(SCORPIUS)
library(ggplot2)
library(xlsx)
library(dplyr)
library(tibble)
library(reshape2)
```

## Load dataset

``` r
set.seed(100)

# load data
#note: this file should include metadata information about cluster "$Annotation" and genotype "$group" (see documentation for Figure 3a)

#BEC.integrated<-readRDS("Your/File/Directory/"Integrated_BEC.rds")

```

# WT

## Subset dataset for WT trajectory

``` r
Idents(BEC.integrated)<-BEC.integrated$group
BEC_wt<-subset(BEC.integrated,  idents = "WT")
Idents(BEC_wt)<- BEC_wt$Annotation

#exclude Tip-cell cluster from analysis
BEC_wt_subset<-subset(BEC_wt, idents=c("Artery","aCap","vCap","Vein"))
```

## Calculate and draw trajectory for WT dataset

``` r
# set cluster colors
bec_cols_wt <- c("#E377C2FF", "#7F7F7FFF","#D62728FF","#1F77B4FF")
bec_cols_wt <- setNames(bec_cols_wt,as.character(sort(unique(BEC_wt_subset$Annotation))))

expression <- t(as.matrix(BEC_wt_subset@assays$RNA@data))
groups <- BEC_wt_subset$Annotation
space <- as.matrix(Embeddings(BEC_wt_subset,"umap"))
traj <- infer_trajectory(space,k=4)

BEC_wt_trajectory <- UMAPPlot(BEC_wt_subset,cols=bec_cols_wt, pt.size=3)+
  geom_path(data=data.frame(traj$path), aes(x=Comp1, y=Comp2), size=2, alpha=0.7)

print(BEC_wt_trajectory)
```

![](Figure-3d-and-4b_files/figure-gfm/unnamed-chunk-3-1.png)<!-- -->

## Calculate gene importance (gimp) construct trajectory heatmap for WT dataset (Fig. 3d)

``` r
# calculate gene importance
gimp_wt <- gene_importances(expression, traj$time, num_permutations = 0, ntree = 10)
gene_sel_wt <- gimp_wt[1:100,]
expr_sel_wt <- expression[,gene_sel_wt$gene]

Heatmap_wt<-draw_trajectory_heatmap(
  expr_sel_wt, traj$time, progression_group=groups,show_labels_row=T,scale_features = TRUE,
  progression_group_palette = setNames(bec_cols_wt,as.character(sort(unique(BEC_wt_subset$Annotation))))
)
```

![](Figure-3d-and-4b_files/figure-gfm/unnamed-chunk-4-1.png)<!-- -->

``` r
# Save as Pdf
ggsave("Fig3d_wt.pdf", Heatmap_wt, width = 8, height = 12)
```

## Heatmap of Sox17 and Lrg1 for WT dataset (Fig.4b)

``` r
expr_sel_wt_markers <- FetchData(BEC_wt_subset, vars = c("Sox17","Lrg1"))

Heatmap_wt_markers<-draw_trajectory_heatmap(
  expr_sel_wt_markers, traj$time, progression_group=groups,show_labels_row=T,scale_features = TRUE,
  progression_group_palette = setNames(bec_cols_wt,as.character(sort(unique(BEC_wt_subset$Annotation))))
)
```

![](Figure-3d-and-4b_files/figure-gfm/unnamed-chunk-5-1.png)<!-- -->

``` r
ggsave("Fig4b_wt_markers.pdf", Heatmap_wt_markers, width = 6, height = 3)
```

## Save heatmap source data WT dataset (Fig.4b)

``` r
#create data frame for marker expression and pseudotime information for each cluster
common_cells <- intersect(colnames(BEC_wt_subset), rownames(expr_sel_wt_markers))
cluster_info <- Idents(BEC_wt_subset)[common_cells]
source_data <- data.frame(expr_sel_wt_markers, cluster = cluster_info)
source_data <- data.frame(expr_sel_wt_markers, pseudotime = traj$time)

#save data frame
write.xlsx(source_data,file="Fig4b_wt_sourcedata.xlsx")
```

## Expression profile of Sox17 and Lrg1 for WT dataset (Fig.4b)

``` r
# create data frame for relative gene expression 
expr_wt <-AverageExpression(object = BEC_wt_subset,group.by = "ident",slot = "scale.data")
expr_wt <- as.data.frame(expr_wt)
expr_wt <- expr_wt %>%
  rownames_to_column(var = "gene")

# select genes
marker_gene_names <- c("Lrg1","Sox17")

# create data frame for selected genes
expr_wt <- expr_wt[expr_wt$gene %in% marker_gene_names, ]
order_index <- match(expr_wt$gene, marker_gene_names)
expr_wt <- expr_wt[order(order_index), ]
expr_wt <- melt(expr_wt, id.vars = "gene")

# save data frame
write.xlsx(expr_wt,file="Fig4b_Expressionprofile_wt.xlsx")

# plot expression profile
Expressionprofile_wt<-ggplot(expr_wt, aes(x = variable, y = value, color = gene, group = gene)) +
  geom_line() +
  labs(title = "Gene Expression Plot", x = "Clusters", y = "Expression Value") +
  theme_minimal() + theme(axis.text.x = element_text(angle = 90))

# save as Pdf
ggsave("Fig4b_Expressionprofile_wt.pdf", Expressionprofile_wt, width = 4, height = 3)
```

# Mutant

## Subset dataset for mutant trajectory

``` r
Idents(BEC.integrated)<-BEC.integrated$group
BEC_mutant<-subset(BEC.integrated, idents="mutation")
Idents(BEC_mutant)<-BEC_mutant$Annotation

#exclude Tip-cell and proliferating Pik3ca-3 cluster from analysis
BEC_mutant_subset<-subset(BEC_mutant,idents=c("Artery","aCap","vCap","Pik3ca-1","Pik3ca-2"))
```

## Calculate and draw trajectory for mutant dataset

``` r
bec_cols_mutant <- c("#E377C2FF", "#7F7F7FFF","#D62728FF","#BCBD22FF","#2CA02CFF")

bec_cols_mutant <- setNames(bec_cols_mutant,as.character(sort(unique(BEC_mutant_subset$Annotation))))

expression <- t(as.matrix(BEC_mutant_subset@assays$RNA@data))
groups <- BEC_mutant_subset$Annotation
space <- as.matrix(Embeddings(BEC_mutant_subset,"umap"))
traj <- infer_trajectory(space,k=4)
BEC_mutant_trajectory <- UMAPPlot(BEC_mutant_subset,cols=bec_cols_mutant, pt.size=3)+
  geom_path(data=data.frame(traj$path), aes(x=Comp1, y=Comp2), size=2, alpha=0.7)

print(BEC_mutant_trajectory)
```

![](Figure-3d-and-4b_files/figure-gfm/unnamed-chunk-9-1.png)<!-- -->

## Calculate gene importance (gimp) construct trajectory heatmap for mutant dataset (Fig. 3d)

``` r
# calculate gene importance
gimp_mutant <- gene_importances(expression, traj$time, num_permutations = 0, ntree = 10)
gene_sel_mutant <- gimp_mutant[1:100,]
expr_sel_mutant <- expression[,gene_sel_mutant$gene]

Heatmap_mutant<-draw_trajectory_heatmap(
  expr_sel_mutant, traj$time, progression_group=groups,show_labels_row=T,scale_features = TRUE,
  progression_group_palette = setNames(bec_cols_mutant,as.character(sort(unique(BEC_mutant_subset$Annotation))))
)
```

![](Figure-3d-and-4b_files/figure-gfm/unnamed-chunk-10-1.png)<!-- -->

``` r
# Save as Pdf
ggsave("Fig3d_mutant.pdf", Heatmap_mutant, width = 8, height = 12)
```

## Heatmap of Sox17 and Lrg1 for mutant dataset (Fig.4b)

``` r
expr_sel_mutant_markers <- FetchData(BEC_mutant_subset, vars = c("Sox17","Lrg1"))

Heatmap_mutant_markers<-draw_trajectory_heatmap(
  expr_sel_mutant_markers, traj$time, progression_group=groups,show_labels_row=T,scale_features = TRUE,
  progression_group_palette = setNames(bec_cols_mutant,as.character(sort(unique(BEC_mutant_subset$Annotation))))
)
```

![](Figure-3d-and-4b_files/figure-gfm/unnamed-chunk-11-1.png)<!-- -->

``` r
ggsave("Fig4b_mutant_markers.pdf", Heatmap_mutant_markers, width = 6, height = 3)
```

## Save heatmap source data mutant dataset (Fig.4b)

``` r
#create data frame for marker expression and pseudotime information for each cluster
common_cells <- intersect(colnames(BEC_mutant_subset), rownames(expr_sel_mutant_markers))
cluster_info <- Idents(BEC_mutant_subset)[common_cells]
source_data_mutant <- data.frame(expr_sel_mutant_markers, cluster = cluster_info)
source_data_mutant <- data.frame(expr_sel_mutant_markers, pseudotime = traj$time)

#save data frame
write.xlsx(source_data_mutant,file="Fig4b_mutant_sourcedata.xlsx")
```

## Expression profile of Sox17 and Lrg1 for WT dataset (Fig.4b)

``` r
# create data frame for relative gene expression 
expr_mutant <-AverageExpression(object = BEC_mutant_subset,group.by = "ident",slot = "scale.data")
expr_mutant <- as.data.frame(expr_mutant)
expr_mutant <- expr_mutant %>%
  rownames_to_column(var = "gene")

# select genes
marker_gene_names <- c("Lrg1","Sox17")

# create data frame for "selected genes
expr_mutant <- expr_mutant[expr_mutant$gene %in% marker_gene_names, ]
order_index <- match(expr_mutant$gene, marker_gene_names)
expr_mutant <- expr_mutant[order(order_index), ]
expr_mutant <- melt(expr_mutant, id.vars = "gene")

# save data frame
write.xlsx(expr_mutant,file="Fig4b_Expressionprofile_mutant.xlsx")

# plot expression profile
Expressionprofile_mutant<-ggplot(expr_mutant, aes(x = variable, y = value, color = gene, group = gene)) +
  geom_line() +
  labs(title = "Gene Expression Plot", x = "Clusters", y = "Expression Value") +
  theme_minimal() + theme(axis.text.x = element_text(angle = 90))

# save as Pdf
ggsave("Fig4b_Expressionprofile_mutant.pdf", Expressionprofile_mutant, width = 4, height = 3)
```
