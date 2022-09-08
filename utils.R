data.names.map <- c(all = 'ALL', DGs = 'DG',  Microglia = 'MG', Astrocytes = 'AST', Oligodendrocytes = 'OLG')
params.palette <- c(`WT-CD` = '#aedaea', `WT-HFD` = '#89aeb7', `5xFAD-CD` = '#f7ccc9', `5xFAD-HFD` = '#cd8d82')

load_as_seurat <- function(name){
    require(dplyr)
    exp <- Read10X(file.path(here::here(), '..', 'data', name))
    obj <- CreateSeuratObject(counts = exp)
    meta.data <- read.csv(file.path(here::here(), '..', 'data', name, 'metadata.csv'))
    rownames(meta.data) <- meta.data$X
    obj <- AddMetaData(obj, dplyr::select(meta.data, -X))
    umap <- read.csv(file.path(here::here(), '..', 'data', name, 'umap.csv'))
    rownames(umap) <- umap$X
    obj[['umap']] <- CreateDimReducObject(embeddings = as.matrix(dplyr::select(umap, -X)), key = "UMAP_")
    obj <- NormalizeData(obj) %>% FindVariableFeatures(., nfeatures=2500) %>% ScaleData()
    return(obj)
}

cd_hfd_de <- function(obj){
    obj <- subset(obj, subset=Genotype == '5xFAD')
    obj <- subset(obj, subset=orig.ident == 640, invert=TRUE)
    Idents(obj) <- obj$Diet

    future::plan('multicore', workers = 4)
    m <- FindMarkers(obj, assay='RNA', test.use='MAST', latent.vars=c('Batch'),
                         ident.1='HFD', ident.2='CD', logfc.threshold=0)
    future::plan('sequential')

    m$gene <- rownames(m)
    m <- m[-grep('^Gm|Rik$', m$gene), ]
    m$qvalue <- p.adjust(m$p_val, 'fdr')
    return(m)
}

plot_cell_type <- function(cell.type){
    obj <- load_as_seurat(data.names.map[[cell.type]])
    p1 <- DimPlot(obj, group.by='State', label=TRUE) + NoLegend()
    p2 <- freq_box_plot(subset(obj, subset=orig.ident %in% c(636, 640), invert=TRUE), group.by='State')
    cowplot::plot_grid(plotlist=list(p1, p2), nrow=1, rel_widths=c(0.25, 0.75))
}

create_freq_df <- function(obj, group.by, param, batch = 'Batch') {
    require(dplyr); require(tidyr)

    if (group.by == 'active.ident') {
        group.by.col <- obj@active.ident
    } else {
        group.by.col <- obj@meta.data[, group.by]
    }

    df <- data.frame(orig.ident = obj$orig.ident, group.by = group.by.col) %>%
        count(orig.ident, group.by, .drop = FALSE, name = 'n') %>%
        complete(orig.ident, group.by, fill = list(n = 0)) %>%
        group_by(orig.ident, .drop = FALSE) %>%
        mutate(freq = n / sum(n) * 100)


    param.df <- data.frame(orig.ident = obj$orig.ident, param = obj@meta.data[, param], batch = obj@meta.data[, batch]) %>%
        distinct()

    df <- df %>% left_join(param.df, by = 'orig.ident')
    return(df)
}

freq_box_plot <- function(obj, group.by = "active.ident", param = "Model", x.lab = "Clusters",
                          y.lab = "Frequency", scales = "fixed", legend.position = "right", ...) {
    require(dplyr); require(rstatix); require(ggpubr)

    if (length(unique(obj@meta.data[, param])) == 1){
        return(NULL)
    }

    df <- create_freq_df(obj, group.by = group.by, param = param)

    pwc <- tryCatch({
        pwc <- df %>%
            group_by(group.by) %>%
            pairwise_wilcox_test(freq ~ param, pool.sd = FALSE, p.adjust.method = "bonferroni")
        pwc <- pwc %>%
            add_xy_position(step.increase = 0)
        pwc$y.position <- pwc$y.position + 0.01
        pwc
    }, error = function(cond){
        return(NULL)
    }
    )

    p <- ggplot(df, aes(x = param, y = freq)) +
        geom_boxplot(aes(fill = param, color = param), size = 0.5, outlier.shape = NA) +
        geom_dotplot(aes(color = param), binaxis = "y", stackdir = "center", dotsize = 0.25, fill = "white",
                     stroke = 0.5, position = position_dodge(width = 0.8), stackratio = 1.1) +
        facet_wrap(~group.by, scales = scales, nrow = 1, labeller = labeller(group.by = label_wrap_gen(10))) +
        labs(x = x.lab, y = y.lab, fill = param, color = param) +
        theme(panel.background = element_rect(fill = "white", color = "gray"),
              strip.background = element_rect(fill = "white", color = "gray"),
              strip.text = element_text(size = 10, family = "sans"), text = element_text(size = 10, family = "sans"),
              axis.ticks.x = element_blank(), axis.text.x = element_blank(),
              axis.ticks.y = element_line(size = 0.5), legend.text = element_text(size = 10, family = "sans"),
              legend.key = element_rect(fill = "white"), legend.position = legend.position) +
        scale_fill_manual(values = params.palette[unique(df$param)]) +
        scale_color_manual(values = params.palette[unique(df$param)])

    if (! is.null(pwc)){
        p <- p +
            stat_pvalue_manual(pwc, label = "p.adj", hide.ns = T, step.increase = 0.02,
                               bracket.size = 0.1, tip.length = 0.01, label.size = 2.5) +
            scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0.05, 0.15)))
    }

    return(p)
}

VolcanoPlot <- function(m, FC = 0.25, pval = 0.01, tname = '') {
    require(ggrepel)
    eps <- 1e-500

    volc <- ggplot(m, aes(x = avg_log2FC, y = -log(qvalue + eps))) +
        geom_point(size = 0.5) +
        ggtitle(label = '', subtitle = tname) +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              panel.background = element_blank(), axis.line = element_line(colour = "black"))

    volc <- volc +
        geom_text_repel(data = m[m$avg_log2FC > FC & m$qvalue < pval,], aes(label = gene), size = 3, color = "steelblue4") +
        geom_text_repel(data = m[m$avg_log2FC < -FC & m$qvalue < pval,], aes(label = gene), size = 3, color = "plum4")
    return(volc)
}


