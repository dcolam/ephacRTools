library(ephacRTools)
library(ggplot2)
library(hexSticker)
library(ggsci)

mypal <- pal_aaas("default", alpha = 0.9)(9)
mypal
#> [1] "#E64B35B2" "#4DBBD5B2" "#00A087B2" "#3C5488B2" "#F39B7FB2" "#8491B4B2"
#> [7] "#91D1C2B2" "#DC0000B2" "#7E6148B2"

scales::show_col(mypal)

data("se_pn")

se_pn <- reducedDim.Cellwise(se_pn, assayList = c("Minima", "Maxima"))

melted_pn <- sechm::meltSE(se_pn, features=row.names(se_pn), assayName = c("Minima", "Maxima"))
melted_pn <- reshape2::melt(melted_pn, measure.vars = c("Minima", "Maxima"))

melted_pn$value <- melted_pn$value*10^10
p <- ggplot(subset(melted_pn), aes(x = V_Clamp, y = value, color=variable, group=interaction(variable))) +
  stat_summary(geom='errorbar',fun.data=mean_se, linewidth=0.5, alpha=0.6) +
  stat_summary(geom='line', fun = "mean", linewidth=0.5, alpha=1) + scale_color_aaas(alpha=0.7)
p1 <- p + theme_void() + theme_transparent() + guides(color="none", ) + facet_grid(~cluster.umap)+ theme(strip.background = element_blank(),
                                                                                                         strip.text.x = element_blank())

p1


sticker(p1, package="ephacRTools", p_size=6, s_x=1, s_y=.75, s_width=1.3, s_height=1,
        filename="inst/figures/ggplot2.svg", h_fill = "#008280E5", h_color = "#38595c", white_around_sticker = T)

redDim.melt <- reducedDims(se_pn)$UMAP
redDim.melt$Well <- row.names(redDim.melt)
melted_pn <- merge(melted_pn, redDim.melt, by="Well")


p2 <- ggplot(data=melted_pn, aes(x=UMAP1, y=UMAP2, color=cluster.tsne)) + geom_point(size=0.01, alpha=0.05)+ theme_void() + theme_transparent() + guides(color="none")+ theme(strip.background = element_blank(),
                                                                                                                                                                            strip.text.x = element_blank())
p2

p <- ggpubr::ggarrange(p1, p2, nrow=2)

sticker(p1, package="ephacRTools", p_size=20, s_x=1, s_y=.75, s_width=1.3, s_height=1,
        filename="inst/figures/ggplot2.png", h_fill = "#3a838a", h_color = "#38595c")
