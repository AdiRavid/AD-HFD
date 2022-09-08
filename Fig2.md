---
title: "Figure 2"
output:
  html_document:
    keep_md: yes
  pdf_document: default
---


Setup:


Panel a

```r
obj <- load_as_seurat('ALL')
DimPlot(obj, group.by='CellType', label=TRUE) + NoLegend()
```

![](/Volumes/habib-lab/adi.ravid/AD-HFD/AD-HFD/Fig2_files/figure-html/unnamed-chunk-2-1.png)<!-- -->
Panels b, c

```r
plot_cell_type('DGs')
```

![](/Volumes/habib-lab/adi.ravid/AD-HFD/AD-HFD/Fig2_files/figure-html/unnamed-chunk-3-1.png)<!-- -->
Panel d

```r
#OR
```
Panels e,f

```r
plot_cell_type('Microglia')
```

![](/Volumes/habib-lab/adi.ravid/AD-HFD/AD-HFD/Fig2_files/figure-html/unnamed-chunk-5-1.png)<!-- -->
Panels e,f

```r
plot_cell_type('Astrocytes')
```

![](/Volumes/habib-lab/adi.ravid/AD-HFD/AD-HFD/Fig2_files/figure-html/unnamed-chunk-6-1.png)<!-- -->
Panels e,f

```r
plot_cell_type('Oligodendrocytes')
```

![](/Volumes/habib-lab/adi.ravid/AD-HFD/AD-HFD/Fig2_files/figure-html/unnamed-chunk-7-1.png)<!-- -->
Panels k,l,m

```r
volcs <- lapply(names(data.names.map[3:length(data.names.map)]), function(cell.type){
    obj <- load_as_seurat(data.names.map[[cell.type]])
    degs <- cd_hfd_de(obj)
    return(VolcanoPlot(degs, tname=paste(cell.type, '5xFAD','HFD vs CD', sep=' ,')))
})
cowplot::plot_grid(plotlist=volcs, nrow=1)
```

![](/Volumes/habib-lab/adi.ravid/AD-HFD/AD-HFD/Fig2_files/figure-html/unnamed-chunk-8-1.png)<!-- -->

