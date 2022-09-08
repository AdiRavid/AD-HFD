Figure 2
================

Setup:

Panel a

``` r
obj <- load_as_seurat('ALL')
DimPlot(obj, group.by='CellType', label=TRUE) + NoLegend()
```

![](/Volumes/habib-lab/adi.ravid/AD-HFD/AD-HFD/Fig2_files/figure-gfm/a-1.png)<!-- -->
Panels b, c

``` r
plot_cell_type('DGs')
```

![](/Volumes/habib-lab/adi.ravid/AD-HFD/AD-HFD/Fig2_files/figure-gfm/b-c-1.png)<!-- -->
Panel d

``` r
#OR
```

Panels e,f

``` r
plot_cell_type('Microglia')
```

![](/Volumes/habib-lab/adi.ravid/AD-HFD/AD-HFD/Fig2_files/figure-gfm/e-f-1.png)<!-- -->
Panels g,h

``` r
plot_cell_type('Astrocytes')
```

![](/Volumes/habib-lab/adi.ravid/AD-HFD/AD-HFD/Fig2_files/figure-gfm/g-h-1.png)<!-- -->
Panels i,j

``` r
plot_cell_type('Oligodendrocytes')
```

![](/Volumes/habib-lab/adi.ravid/AD-HFD/AD-HFD/Fig2_files/figure-gfm/i-j-1.png)<!-- -->
Panels k,l,m

``` r
volcs <- lapply(names(data.names.map[3:length(data.names.map)]), function(cell.type){
    obj <- load_as_seurat(data.names.map[[cell.type]])
    degs <- cd_hfd_de(obj)
    return(VolcanoPlot(degs, tname=paste(cell.type, '5xFAD','HFD vs CD', sep=' ,')))
})
cowplot::plot_grid(plotlist=volcs, nrow=1)
```

![](/Volumes/habib-lab/adi.ravid/AD-HFD/AD-HFD/Fig2_files/figure-gfm/k-l-m-1.png)<!-- -->
