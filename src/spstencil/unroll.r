library(magick)
library(imager)
# library(STutility)
library(zeallot)
library(magrittr)
library(dplyr)
library(dbscan)
library(ggplot2)

setwd("/share//nas1/gaoy/Project/BMK240416-BZ180-ZX01-05/Dig/")
source("CreateBmkObject_v2.R")
# matrix_path = "../BMK_2_BSTMatrix_analysis/BMK_2_BSTViewer/S115.BSTViewer_project/subdata/L7_heAuto/",  #矩阵文件目录
# 重新选择roi后，使用更新的subdata结果
single.seurat <- CreateBmkObject(
    matrix_path = "../ROI/S115.BSTViewer_project/roi_0/subdata/L7_heAuto/", # 矩阵文件目录
    png_path = "../BMK_2_BSTMatrix_analysis/BMK_2_BSTViewer/S115.BSTViewer_project/he_roi_small.png", # png格式图片路径
    min.cells = 5, # 一个基因至少在n个细胞中表达才被保留，可自行调整，默认值5
    min.features = 100, # 一个细胞至少有n个基因才被保留，可自行调整，默认值100 #作图点的半径，细胞分割数据必须指定
    type = "S3000" # 数据类型，必须指定，可选参数有'S1000', 'S2000', 'S3000', 'cell_split'
)

spots <- single.seurat@images[["sample1"]]@coordinates[, c("imagecol", "imagerow")]
zoom <- single.seurat@images$sample1@scale.factors$spot
zoom
spots$imagecol <- spots$imagecol * zoom
spots$imagerow <- spots$imagerow * zoom

x <- image_read("./S115.BSTViewer_project/S115_he_roi_small.png") %>% as.raster()
xy.all <- as.data.frame(which(x == "#000000ff", arr.ind = T))
xy.all$layer <- 1
plot(xy.all[, 2:1], xlim = c(0, 1024), ylim = c(1026, 0), cex = 0.1, col = xy.all$layer)
points(spots, cex = 0.3, pch = 19, col = "gray")

p <- ggplot(spots, aes(imagecol, 2000 - imagerow)) +
    geom_point(size = 1, color = "orange") +
    theme_void()

p <- ggplot(xy.all, aes(col, 2000 - row)) +
    geom_point(size = 0.1, color = "black") +
    theme_void()

set1 <- spots[, c("imagecol", "imagerow")]
set2 <- xy.all[, 2:1]
mindists <- apply(set1, 1, function(x) {
    which.min(sqrt(colSums((t(set2) - x)^2)))
})

gg <- cbind(xy.all[mindists, ], spots)
ggplot() +
    geom_segment(data = gg, aes(x = col, xend = imagecol, y = row, yend = imagerow)) +
    scale_y_reverse()

library(dbscan)

# Add a unique id for each spot
xy.all$id <- paste0("id", 1:nrow(xy.all))
rownames(xy.all) <- xy.all$id

# Find k nearest neighbors
y <- dbscan::kNN(x = xy.all[, 1:2], k = 5, sort = FALSE)

# Create a data.frame in long format specifying the connections with their associated distances
adj <- rbind(
    data.frame(from = xy.all$id, to = xy.all$id[y$id[, 1]], d = y$dist[, 1]),
    data.frame(from = xy.all$id, to = xy.all$id[y$id[, 2]], d = y$dist[, 2]),
    data.frame(from = xy.all$id, to = xy.all$id[y$id[, 3]], d = y$dist[, 3]),
    data.frame(from = xy.all$id, to = xy.all$id[y$id[, 4]], d = y$dist[, 4]),
    data.frame(from = xy.all$id, to = xy.all$id[y$id[, 5]], d = y$dist[, 5])
)

# Convert data.frame into an adjecency matrix
adj.mat <- reshape2::acast(data = adj, formula = from ~ to, value.var = "d", fill = 0)

# Convert adjacency matrix into igraph object
ig <- igraph::graph_from_adjacency_matrix(adj.mat)

end.points <- names(sort(rowMeans(y$dist), decreasing = TRUE))[1:2]

rownames(xy.all) <- xy.all$id
gg <- data.frame(xy.all)
gg$hl <- ifelse(gg$id %in% end.points, gg$id, "")
p <- ggplot() +
    geom_point(data = gg, aes(col, row)) +
    geom_text(data = gg, aes(col, row, label = hl)) +
    theme_void()

sp.full <- igraph::shortest_paths(
    graph = ig,
    from = end.points[1],
    to = end.points[2]
)

xy.all <- xy.all[names(sp.full$vpath[[1]]), ]
xy.all$ord <- 1:nrow(xy.all)

p <- ggplot() +
    geom_point(data = xy.all, aes(col, 2000 - row)) +
    geom_label(data = xy.all[seq(1, nrow(xy.all), length.out = 40), ], aes(col, 2000 - row, label = ord)) +
    theme_void()

library(LearnGeom)

adjust.coords <- function(x, y, P) {
    avg <- c(mean(c(x[1], y[1])), mean(c(x[2], y[2])))
    df <- (P - avg)
    x_new <- x - df
    y_new <- y - df
    return(c(x_new, y_new))
}
sider <- function(x, y, P) {
    (P[1] - x[1]) * (y[2] - x[2]) - (P[2] - x[2]) * (y[1] - x[1])
}
angle <- function(x, y, P) {
    Angle(A = P, B = c(mean(c(x[1], y[1])), mean(c(x[2], y[2]))), C = y)
}
sign_angle <- function(x, y, P) {
    angle(x, y, P) * sign(sider(x, y, P))
}
expand.range <- function(x, exp.factor = 5, maxval = 1e4) {
    y <- c(max(1, x - exp.factor), min(x + exp.factor, maxval))
    return(y)
}

library(pbmcapply)
set1 <- spots[, c("imagecol", "imagerow")]
set2 <- xy.all[, 2:1]
dists <- as.matrix(apply(set1, 1, function(x) {
    sqrt(colSums((t(set2) - x)^2))
}))

spots.list <- do.call(rbind, pbmcapply::pbmclapply(colnames(dists), function(s) {
    P <- spots[s, ]
    dists.subset <- dists[, s]
    dists.subset <- dists.subset[dists.subset < 250]
    checks <- sapply(names(dists.subset), function(b) {
        xs <- expand.range(xy.all[b, "ord"] %>% as.numeric(),
            exp.factor = 5, maxval = nrow(xy.all)
        )
        xs <- xy.all[xs, ]
        ls <- adjust.coords(x = as.numeric(xs[1, 2:1]), y = as.numeric(xs[2, 2:1]), P = as.numeric(P))
        x <- ls[1:2]
        y <- ls[3:4]
        y <- sign_angle(x = x, y = y, P = as.numeric(P))
        135 > y && y > 45
    })
    if (sum(checks) > 0) {
        bl.subset <- dists.subset[checks]
        return(data.frame(dist = bl.subset, id = names(bl.subset), spot = s))
    } else {
        return(NULL)
    }
}))

spots.summarized <- spots.list %>%
    group_by(spot) %>%
    top_n(n = 1, wt = -dist)

g <- merge(xy.all, spots.summarized, by = "id", all = T)
g <- cbind(g, spots[g$spot, ])
g <- na.omit(g)

im <- image_read("../BMK_2_BSTMatrix_analysis/BMK_2_BSTViewer/S115.BSTViewer_project/he_roi_small.png") %>% as.raster()
gr <- grid::rasterGrob(im, width = unit(1, "npc"), height = unit(1, "npc"), interpolate = TRUE)

p <- ggplot() +
    annotation_custom(gr, -Inf, Inf, -Inf, Inf) +
    geom_point(data = g, aes(col, dim(im)[1] - row), color = "red", size = 0.1) +
    geom_point(data = g, aes(imagecol, dim(im)[1] - imagerow), color = "blue", shape = 21, size = 0.1) +
    geom_segment(data = g, aes(x = col, xend = imagecol, y = dim(im)[1] - row, yend = dim(im)[1] - imagerow), color = "black", size = 0.3) +
    theme_void() +
    scale_x_continuous(limits = c(0, dim(im)[2]), expand = c(0, 0)) +
    scale_y_continuous(limits = c(0, dim(im)[1]), expand = c(0, 0))

single.seurat <- NormalizeData(single.seurat, normalization.method = "LogNormalize", scale.factor = 10000, assay = "Spatial")
single.seurat <- ScaleData(single.seurat, model.use = "negbinom", do.scale = T, do.center = T, assay = "Spatial")
single.seurat <- SCTransform(single.seurat, assay = "Spatial", verbose = FALSE, return.only.var.genes = FALSE)

proximal <- c("Car1", "Mettl7b", "Emp1", "Fabp2")
mid <- c("Retnlb", "Sprr2a2")
distal <- c("Prdx6", "Tgm3", "Ly6g", "Eno3")
ftrs <- c(proximal, mid, distal)

ggs <- cbind(g, t(as.matrix(single.seurat@assays$SCT@data[ftrs, g$spot])))
plotlist <- lapply(ftrs, function(ftr) {
    ggs <- ggs[order(ggs[, ftr], decreasing = F), ]
    p <- ggplot(ggs, aes_string("ord", "dist",
        z = paste0("`", ftr, "`"),
        color = paste0("`", ftr, "`")
    )) +
        stat_summary_hex(binwidth = c(20, 3)) +
        scale_fill_gradientn(colours = c("lightgray", "mistyrose", "red", "darkred")) +
        theme_void() +
        theme(
            plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), units = "cm")
        ) +
        labs(x = "position along colon", y = "thickness", title = ftr)
})
p <- cowplot::plot_grid(plotlist = plotlist, ncol = 2)

ftrs <- c("Hmgcs2", "Ang4", "B4galt1")

ggs <- cbind(g, t(as.matrix(single.seurat@assays$SCT@data[ftrs, g$spot])))
plotlist <- lapply(ftrs, function(ftr) {
    ggs <- ggs[order(ggs[, ftr], decreasing = F), ]
    p <- ggplot(ggs, aes_string("ord", "dist",
        z = paste0("`", ftr, "`"),
        color = paste0("`", ftr, "`")
    )) +
        stat_summary_hex(binwidth = c(20, 3)) +
        scale_fill_gradientn(colours = c("lightgray", "mistyrose", "red", "darkred")) +
        theme_void() +
        theme(
            plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), units = "cm")
        ) +
        labs(x = "position along colon", y = "thickness", title = ftr)
})
p <- cowplot::plot_grid(plotlist = plotlist, ncol = 1)

proximal <- c("Car1", "Mettl7b", "Emp1", "Fabp2")
mid <- c("Retnlb", "Sprr2a2")
distal <- c("Prdx6", "Tgm3", "Ly6g", "Eno3")
ftrs <- c(proximal, mid, distal)

ggs <- reshape2::melt(cbind(g, t(single.seurat@assays$SCT@scale.data[ftrs, g$spot])), measure.vars = ftrs)
rs <- ggs %>%
    mutate(bin = ntile(ord, 50)) %>%
    group_by(bin, variable) %>%
    summarize(avg = mean(value)) %>%
    group_by(variable) %>%
    mutate(avg.norm = scales::rescale(avg, c(-1, 1)))

p <- ggplot() +
    geom_tile(data = rs, aes(bin, variable, fill = avg.norm)) +
    scale_fill_gradientn(colours = RColorBrewer::brewer.pal(n = 11, name = "RdBu") %>% rev()) +
    theme_minimal() +
    theme(panel.grid = element_blank(), axis.text.x = element_blank()) +
    labs(x = "", fill = "expression")