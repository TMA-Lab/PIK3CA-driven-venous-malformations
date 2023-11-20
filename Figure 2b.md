Fig.2b: UMAP representation of 321 dermal BECs from control mice 
================
Marle Kraft
2023-10-02

- <a href="#packages" id="toc-packages">Packages</a>
- <a href="#data-processing---control-becs"
  id="toc-data-processing---control-becs">Data processing - control
  BECs</a>
- <a href="#umap-clustering-representation---control-becs"
  id="toc-umap-clustering-representation---control-becs">UMAP clustering
  representation - control BECs</a>
- <a href="#deg-analysis---control-becs"
  id="toc-deg-analysis---control-becs">DEG analysis - control BECs</a>
- <a href="#removal-of-contamination-and-re-analysis---control-becs"
  id="toc-removal-of-contamination-and-re-analysis---control-becs">Removal
  of contamination and re-analysis - control BECs</a>
- <a href="#deg-analysis-of-cleaned-dataset---control-becs"
  id="toc-deg-analysis-of-cleaned-dataset---control-becs">DEG analysis of
  cleaned dataset - control BECs</a>
- <a href="#fig-2b-plot---control-becs"
  id="toc-fig-2b-plot---control-becs">Fig. 2b plot - control BECs</a>
- <a href="#save-final-dataset---control-becs"
  id="toc-save-final-dataset---control-becs">Save final dataset - control
  BECs</a>

## Packages

``` r
library(Seurat)
library(xlsx)
library(ggplot2)
```

## Data processing - control BECs

``` r
# load data
# bec_complete<-readRDS("/Your/File/Directory/bec_complete.rds")

# subset control BECs only
Idents(bec_complete)<-bec_complete$plate
control<-subset(bec_complete, idents=c("895s","901s","903s"))

# batch effect correction
control.list <- SplitObject(object =control, split.by = "plate")

# prior to finding anchors, we perform standard preprocessing, and identify variable features individually for each plate
for (i in 1:length(control.list)) {
  control.list[[i]] <- NormalizeData(control.list[[i]], verbose = FALSE,scale.factor = 1000000)
  control.list[[i]] <- FindVariableFeatures(control.list[[i]], selection.method = "vst",
                                        nfeatures = 3000, verbose = FALSE)
}

# finding integration anchors
control.anchors <- FindIntegrationAnchors(object.list =control.list, anchor.features = 3000)

# integrate data
control.batchCorr<- IntegrateData(anchorset = control.anchors)

# run the standard workflow for visualization and clustering
DefaultAssay(control.batchCorr) <- "integrated"
control.batchCorr<- ScaleData(object = control.batchCorr, verbose = FALSE)

# dimensionality reduction by pca
control.batchCorr<- RunPCA(object = control.batchCorr, npcs = 30, verbose = FALSE)

# run umap visualization
control.batchCorr<- RunUMAP(object = control.batchCorr, reduction = "pca", dims = 1:10)
control.batchCorr <- FindNeighbors(control.batchCorr, reduction = "pca", dims = 1:10)
control.batchCorr <- FindClusters(control.batchCorr, resolution = 0.8)
```

## UMAP clustering representation - control BECs

``` r
# plot UMAP representation
DimPlot(control.batchCorr, reduction = "umap")
```

![](Fig.2b_control_Analysis_files/figure-gfm/UMAP-1.png)<!-- -->

## DEG analysis - control BECs

``` r
# run DEG analysis to find cluster markers   
control.markers <- FindAllMarkers(control.batchCorr, only.pos = TRUE, logfc.threshold = 0.25, assay="RNA")

# save cluster marker list
write.table(control.markers, "control.markers.txt", sep="\t")
```

## Removal of contamination and re-analysis - control BECs

cluster 4 likely to be contamination and was removed from the further
analysis

``` r
# remove cluster 4 and re-analyse
clean_BEC_control<-subset(control.batchCorr, idents=c("0","1","2","3"))

#re-analysis
DefaultAssay(clean_BEC_control) <- "integrated"
clean_BEC_control<- ScaleData(object = clean_BEC_control, verbose = FALSE)
clean_BEC_control<- RunPCA(object = clean_BEC_control, npcs = 30, verbose = FALSE)
clean_BEC_control<- RunUMAP(object = clean_BEC_control, reduction = "pca", dims = 1:10)
clean_BEC_control <- FindNeighbors(clean_BEC_control, reduction = "pca", dims = 1:10)
clean_BEC_control<- FindClusters(clean_BEC_control, resolution = 0.8)
```

## DEG analysis of cleaned dataset - control BECs

``` r
# DEG analysis
clean_BEC_control.markers <- FindAllMarkers(clean_BEC_control, only.pos = TRUE, logfc.threshold = 0.25, assay="RNA")

# save cluster markers
write.table(clean_BEC_control.markers, "clean_BEC_control.markers.txt", sep="\t")
```

## Fig. 2b plot - control BECs

``` r
#plot UMAP representation
 mycolor=c("#E377C2FF","#7F7F7FFF","#D62728FF","#1F77B4FF")

new.cluster.ids <- c("aCap","vCap","Vein","Artery")
names(new.cluster.ids) <- levels(clean_BEC_control)
clean_BEC_control<- RenameIdents(clean_BEC_control, new.cluster.ids)
levels(clean_BEC_control) <- c("Artery","aCap","vCap","Vein")

DimPlot(clean_BEC_control, reduction = "umap",label=TRUE,cols=mycolor,pt.size = 5, label.size = FALSE)
```

![](Fig.2b_control_Analysis_files/figure-gfm/Fig.2b-1.png)<!-- -->

``` r
  theme(
    legend.text=element_text(size=40),
    axis.text.x = element_text(size=40),
    axis.text.y = element_text(size=40),
    text = element_text(size=40)
  )
```

## Save final dataset - control BECs

``` r
# Save as an RDS file
saveRDS(clean_BEC_control, "clean_BEC_control.rds")
```
