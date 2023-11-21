Fig.3b: Pik3ca transcript status visualized in combined UMAP clustering.
================
Marle Kraft
2023-10-02

- [Packages](#packages)
- [Load data](#load-data)
- [Extraction of Pik3caH1047R transcript
  information](#extraction-of-pik3cah1047r-transcript-information)
- [Visualization (WT and PIK3CA
  status)](#visualization-wt-and-pik3ca-status)
- [Visualization (PIK3CA status
  alone)](#visualization-pik3ca-status-alone)

## Packages

``` r
library(Seurat)
```

    ## Loading required package: SeuratObject

    ## Loading required package: sp

    ## 'SeuratObject' was built under R 4.3.0 but the current version is
    ## 4.3.2; it is recomended that you reinstall 'SeuratObject' as the ABI
    ## for R may have changed

    ## 'SeuratObject' was built with package 'Matrix' 1.6.1.1 but the current
    ## version is 1.6.2; it is recomended that you reinstall 'SeuratObject' as
    ## the ABI for 'Matrix' may have changed

    ## 
    ## Attaching package: 'SeuratObject'

    ## The following object is masked from 'package:base':
    ## 
    ##     intersect

``` r
library(RColorBrewer)
library(ggplot2)
```

## Load data

``` r
#BEC.integrated<-readRDS("Your/File/Directory/"Integrated_BEC.rds")
BEC.integrated<-readRDS("/Users/markr691/Desktop/19042023_manuscript/Finalize manuscript/BEC.dat.rds")
```

## Extraction of Pik3caH1047R transcript information

``` r
#The single-cell sequence data were aligned to the Gallus gallus PIK3CA sequence (National Center for Biotechnology Information [NCBI] nucleotide sequence ID: NM_001004410) and its transcript was named "Chicken_PIK3CA" in our dataset
#This chunk of code identifies the Pik3caH1047R transcript in the transgenic mice, encoded by G. gallus Pik3ca

#Annotate Chicken_PIK3CA positive cells and extract their umap embedding information
Chicken.PIK3CA.pos<-row.names(BEC.integrated@meta.data[BEC.integrated$Chicken_PIK3CA>0,])
umap_embedding<-Embeddings(BEC.integrated,reduction = "umap")
umap.PIK3CA.pos<-data.frame(umap_embedding[Chicken.PIK3CA.pos,],PIK3CA_status="PIK3CA_positive")

#Annotate Chicken_PIK3CA negative mutant and WT cells and extract their umap embedding information
Chicken.PIK3CA.neg<-row.names(BEC.integrated@meta.data[(BEC.integrated$Chicken_PIK3CA==0) &(BEC.integrated$group!="WT"),])
umap.PIK3CA.neg<-data.frame(umap_embedding[Chicken.PIK3CA.neg,],PIK3CA_status="PIK3CA_negative")

WT.cells<-row.names(BEC.integrated@meta.data[(BEC.integrated$group=="WT"),])
umap.WT<-data.frame(umap_embedding[WT.cells,],PIK3CA_status="WT")

#combine PIK3CA pos/neg and WT cells
umap.PIK3CA.pos.neg.WT<-rbind(umap.WT,umap.PIK3CA.neg,umap.PIK3CA.pos)
umap.PIK3CA.pos.neg.WT$PIK3CA_status<-factor(umap.PIK3CA.pos.neg.WT$PIK3CA_status,levels=c("WT","PIK3CA_negative","PIK3CA_positive"))
```

## Visualization (WT and PIK3CA status)

``` r
# plot
ggplot(umap.PIK3CA.pos.neg.WT,aes(x=UMAP_1, y=UMAP_2,col=PIK3CA_status))+
  geom_point(size=2)+
  scale_color_manual(values = c("grey","#377EB8","#E41A1C"))+
  theme_bw()+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        legend.title = element_blank(),
        text = element_text(size=20)
  )
```

![](Figure-3b_files/figure-gfm/unnamed-chunk-3-1.png)<!-- -->

## Visualization (PIK3CA status alone)

``` r
#In the figure we do not display WT cells:

umap.PIK3CA.pos.neg_mutant<-rbind(umap.PIK3CA.neg,umap.PIK3CA.pos)
umap.PIK3CA.pos.neg_mutant$PIK3CA_status_mutant<-factor(umap.PIK3CA.pos.neg_mutant$PIK3CA_status,levels=c("PIK3CA_negative","PIK3CA_positive"))

ggplot(umap.PIK3CA.pos.neg_mutant,aes(x=UMAP_1, y=UMAP_2,col=PIK3CA_status_mutant))+
  geom_point(size=2)+
  scale_color_manual(values = c("#377EB8","#E41A1C"))+
  theme_bw()+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        legend.title = element_blank(),
        text = element_text(size=20)
  )
```

![](Figure-3b_files/figure-gfm/unnamed-chunk-4-1.png)<!-- -->
