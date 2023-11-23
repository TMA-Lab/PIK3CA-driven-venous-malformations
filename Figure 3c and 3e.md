Fig.3c and 3e: Visualizing gene expression of selected venous, aterial
and proliferation marker genes
================
Marle Kraft
2023-10-02

- [Packages](#packages)
- [Load dataset](#load-dataset)
- [Stacked Violine Plot of known venous identity and proliferation
  markers (Fig.
  3c)](#stacked-violine-plot-of-known-venous-identity-and-proliferation-markers-fig-3c)
- [Heatmap: Vein vs Artery (Fig.3e)](#heatmap-vein-vs-artery-fig3e)

## Packages

``` r
library(Seurat)
library(patchwork)
library(ggplot2)
library(dplyr)
library(tibble)
library(RColorBrewer)
```

## Load dataset

``` r
# load data
#note: this file should include metadata information about cluster "$Annotation" and genotype "$group" (see documentation for Figure 3a)

#BEC.integrated<-readRDS("Your/File/Directory/"Integrated_BEC.rds")
BEC.integrated<-readRDS("/Users/markr691/BEC.integrated.rds")

#set cluster colors
mycolor<-c("#E377C2FF", "#7F7F7FFF","#D62728FF","#1F77B4FF","#BCBD22FF","#2CA02CFF","#8C564BFF")

#exclude Tip cell-cluster from analysis below
BEC_subset<-subset(BEC.integrated, idents=c("Artery","aCap","vCap","Vein","Pik3ca-1","Pik3ca-2","Pik3ca-3"))
```

## Stacked Violine Plot of known venous identity and proliferation markers (Fig. 3c)

``` r
#Prepare Stacked Violine Plot

#note1: boarder lines were removed manually in Adobe Illustrator 
#note2: pt.size can be adjusted to show data points

## remove the x-axis text and tick
## plot.margin to adjust the white space between each plot.
## ... pass any arguments to VlnPlot in Seurat
modify_vlnplot<- function(obj,
                          feature,
                          pt.size = 0,
                          plot.margin = unit(c(-0.75, 0, -0.75, 0), "cm"),
                          ...) {
  p<- VlnPlot(obj, features = feature, pt.size = pt.size, cols = mycolor,...)  +
    xlab("") + ylab(feature) + ggtitle("") +
    theme(legend.position = "none",
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.title.y = element_text(size = rel(1), angle = 0),
          axis.text.y = element_text(size = rel(1)),
          plot.margin = plot.margin )
  return(p)
}

## extract the max value of the y axis
extract_max<- function(p){
  ymax<- max(ggplot_build(p)$layout$panel_scales_y[[1]]$range$range)
  return(ceiling(ymax))
}

## main function
StackedVlnPlot<- function(obj, features,
                          pt.size = 0,
                          plot.margin = unit(c(-0.75, 0, -0.75, 0), "cm"),
                          ...) {

  plot_list<- purrr::map(features, function(x) modify_vlnplot(obj = obj,feature = x, ...))

  # Add back x-axis title to bottom plot. patchwork is going to support this?
  plot_list[[length(plot_list)]]<- plot_list[[length(plot_list)]] +
    theme(axis.text.x=element_text(), axis.ticks.x = element_line()) + theme(axis.text.x = element_text(angle = 45))

  # change the y-axis tick to only max value
  ymaxs<- purrr::map_dbl(plot_list, extract_max)
  plot_list<- purrr::map2(plot_list, ymaxs, function(x,y) x +
                            scale_y_continuous(breaks = c(y)) +
                            expand_limits(y = y))

  p<- patchwork::wrap_plots(plotlist = plot_list, ncol = 1)
  return(p)
}

# Plot Stacked Violine Plot

features=c("Nrp2","Vwf","Nr2f2","Ackr1","Cdk1","Mki67")

Violine_plot_integrated<-StackedVlnPlot(BEC_subset, features)

ggsave("Violine_plot_integrated.pdf", Violine_plot_integrated, width = 4, height = 8)
```

## Heatmap: Vein vs Artery (Fig.3e)

``` r
VeinvsArtery.markers <- FindMarkers(BEC.integrated, ident.1="Artery",ident.2="Vein",logfc.threshold = 0.25, assay="RNA")

Artery.select<-VeinvsArtery.markers[VeinvsArtery.markers$avg_log2FC>0,]
Vein.select<-VeinvsArtery.markers[VeinvsArtery.markers$avg_log2FC<0,]

top50_Artery <- Artery.select %>% top_n(50, -p_val)
top50_Vein <- Vein.select %>% top_n(50, -p_val)

top100 <- rbind(top50_Artery, top50_Vein)

top100 <- top100 %>%
  rownames_to_column(var = "gene")

# Subset Artery, Vein, Pik3ca-1 and Pik3ca-2 

BEC_subset_heatmap<-subset(BEC.integrated, idents=c("Artery","Vein","Pik3ca-1","Pik3ca-2"))

mycolors2<-c("#E377C2FF","#1F77B4FF","#BCBD22FF","#2CA02CFF")

ArteryvsVein_heatmap<-DoHeatmap(object = BEC_subset_heatmap,features = top100$gene,group.colors = mycolors2) +
scale_fill_gradientn(colors = rev(RColorBrewer::brewer.pal(n = 10, name = "RdBu")))

ggsave("Fig3e_ArteryvsVein_heatmap.pdf", ArteryvsVein_heatmap, width = 8, height = 12)
```
