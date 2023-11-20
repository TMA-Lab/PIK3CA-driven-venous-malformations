Fig.2c_control_markers
================
Marle Kraft
2023-10-02

- [Packages](#packages)

## Packages

``` r
library(Seurat)
library(xlsx)
library(ggplot2)
```

``` r
# load data
clean_BEC_control<-readRDS("/Users/markr691/Desktop/VM_manuscript/VMsourcecode/clean_BEC_control.rds")
```

``` r
# note: The label text was changed to "Relative Expression" and "% Expression" retrospectively in Adobe Illustrator.

DefaultAssay(clean_BEC_control)<-"RNA"

levels(clean_BEC_control) <- c("Vein","vCap","aCap","Artery")

Bubbleplot<-DotPlot(clean_BEC_control, features = c("Vegfc","Bmx","Notch1","Ly6a","Hey1","Cldn5","Sox17","Dll4","Itga1","Cd36","Nrp1","Aqp7","Cxcl12","Tgfbr3","Pltp","Scarb1","Myc","Angpt2","Col15a1","Nrp2","Emcn","Plvap","Vwf","Apoe","Selp","Sele","Bmp4","Nr2f2","Ackr1","Icam1","Klf2","Vcam1"), scale = TRUE, cols=c("RdYlBu"), dot.scale=5) + RotatedAxis() + theme(axis.text.x = element_text(angle = 90, vjust =-0.05, hjust = -0.05)) + scale_x_discrete(position="top")

Bubbleplot<-Bubbleplot + theme(axis.title.x=element_blank()) + theme(axis.title.y=element_blank()) + theme(legend.position="bottom") + theme(legend.text = element_text(size=7)) + theme(legend.title = element_text(size=10))

ggsave("Figure2c.pdf", Bubbleplot, width = 8, height = 4)
```
