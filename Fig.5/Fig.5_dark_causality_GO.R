
#[1] Plot Dark causality strength
library(ggplot2)
library(patchwork)
load("CRC_GSE146771_AGR2_genes_dark_causality.RData")
# remove causality == 0
causality_dark[causality_dark == 0] <- NA

#AGR2_to_KRT8
t=causality_dark[,c(1,2)]
t <- subset(t, !is.na(t[,2]))
p1 <- ggplot(t, aes(x = Time, y = AGR2_to_KRT8)) +
  geom_line(color = "skyblue", size = 1) +            # Add line (light blue)
  geom_point(color = "darkorange", size = 3) +          # Add points (light red)
  labs(
    x = "Time",
    y = "Dark causality strength"
  ) +
  theme(
    panel.background = element_blank(),     # Remove gray background
    plot.background = element_blank(),      # Optional: remove full plot background
    panel.grid.major = element_blank(),     # Remove major grid lines
    panel.grid.minor = element_blank(),     # Remove minor grid lines
    axis.line = element_line(color = "black"),  # Keep black axis lines
    axis.text = element_text(color = "black"),  # Axis text color
    axis.title = element_text(color = "black")  # Axis title color
  )

#AGR2_to_EPCAM
t=causality_dark[,c(1,3)]
t <- subset(t, !is.na(t[,2]))
p2 <- ggplot(t, aes(x = Time, y = AGR2_to_EPCAM)) +
  geom_line(color = "skyblue", size = 1) +            # Add line (light blue)
  geom_point(color = "darkorange", size = 3) +          # Add points (light red)
  labs(
    x = "Time",
    y = "Dark causality strength"
  ) +
  theme(
    panel.background = element_blank(),     # Remove gray background
    plot.background = element_blank(),      # Optional: remove full plot background
    panel.grid.major = element_blank(),     # Remove major grid lines
    panel.grid.minor = element_blank(),     # Remove minor grid lines
    axis.line = element_line(color = "black"),  # Keep black axis lines
    axis.text = element_text(color = "black"),  # Axis text color
    axis.title = element_text(color = "black")  # Axis title color
  )

#AGR2_to_CLDN4
t=causality_dark[,c(1,4)]
t <- subset(t, !is.na(t[,2]))
p3 <- ggplot(t, aes(x = Time, y = AGR2_to_CLDN4)) +
  geom_line(color = "skyblue", size = 1) +            # Add line (light blue)
  geom_point(color = "darkorange", size = 3) +          # Add points (light red)
  labs(
    x = "Time",
    y = "Dark causality strength"
  ) +
  theme(
    panel.background = element_blank(),     # Remove gray background
    plot.background = element_blank(),      # Optional: remove full plot background
    panel.grid.major = element_blank(),     # Remove major grid lines
    panel.grid.minor = element_blank(),     # Remove minor grid lines
    axis.line = element_line(color = "black"),  # Keep black axis lines
    axis.text = element_text(color = "black"),  # Axis text color
    axis.title = element_text(color = "black")  # Axis title color
  )

#KRT8_to_AGR2
t=causality_dark[,c(1,5)]
t <- subset(t, !is.na(t[,2]))
p4 <- ggplot(t, aes(x = Time, y = KRT8_to_AGR2)) +
  geom_line(color = "skyblue", size = 1) +            # Add line (light blue)
  geom_point(color = "darkorange", size = 3) +          # Add points (light red)
  labs(
    x = "Time",
    y = "Dark causality strength"
  ) +
  theme(
    panel.background = element_blank(),     # Remove gray background
    plot.background = element_blank(),      # Optional: remove full plot background
    panel.grid.major = element_blank(),     # Remove major grid lines
    panel.grid.minor = element_blank(),     # Remove minor grid lines
    axis.line = element_line(color = "black"),  # Keep black axis lines
    axis.text = element_text(color = "black"),  # Axis text color
    axis.title = element_text(color = "black")  # Axis title color
  )

#EPCAM_to_AGR2
t=causality_dark[,c(1,6)]
t <- subset(t, !is.na(t[,2]))
p5 <- ggplot(t, aes(x = Time, y = EPCAM_to_AGR2)) +
  geom_line(color = "skyblue", size = 1) +            # Add line (light blue)
  geom_point(color = "darkorange", size = 3) +          # Add points (light red)
  labs(
    x = "Time",
    y = "Dark causality strength"
  ) +
  theme(
    panel.background = element_blank(),     # Remove gray background
    plot.background = element_blank(),      # Optional: remove full plot background
    panel.grid.major = element_blank(),     # Remove major grid lines
    panel.grid.minor = element_blank(),     # Remove minor grid lines
    axis.line = element_line(color = "black"),  # Keep black axis lines
    axis.text = element_text(color = "black"),  # Axis text color
    axis.title = element_text(color = "black")  # Axis title color
  )

#CLDN4_to_AGR2
t=causality_dark[,c(1,7)]
t <- subset(t, !is.na(t[,2]))
p6 <- ggplot(t, aes(x = Time, y = CLDN4_to_AGR2)) +
  geom_line(color = "skyblue", size = 1) +            # Add line (light blue)
  geom_point(color = "darkorange", size = 3) +          # Add points (light red)
  labs(
    x = "Time",
    y = "Dark causality strength"
  ) +
  theme(
    panel.background = element_blank(),     # Remove gray background
    plot.background = element_blank(),      # Optional: remove full plot background
    panel.grid.major = element_blank(),     # Remove major grid lines
    panel.grid.minor = element_blank(),     # Remove minor grid lines
    axis.line = element_line(color = "black"),  # Keep black axis lines
    axis.text = element_text(color = "black"),  # Axis text color
    axis.title = element_text(color = "black")  # Axis title color
  )

#[2] Plot GO enrichment analysis
library("clusterProfiler")
library("org.Hs.eg.db")
library("DOSE")
library("ggplot2")
library("enrichR")
library("GOplot")
library("enrichplot")
library(ReactomePA)

top_20_genes <- readRDS("top_20_genes.rds")
#output results of top20
x13=top_20_genes[["C5_GSE166555_xtoy"]][,1]
x14=top_20_genes[["C5_GSE166555_ytox"]][,1]
GSE166555_subtype1_C5 = union(x13,x14)
x15=top_20_genes[["C11_GSE166555_xtoy"]][,1]
x16=top_20_genes[["C11_GSE166555_ytox"]][,1]
GSE166555_subtype2_C11 = union(x15,x16)
x17=top_20_genes[["C13_GSE166555_xtoy"]][,1]
x18=top_20_genes[["C13_GSE166555_ytox"]][,1]
GSE166555_subtype3_C13 = union(x17,x18)
x19=top_20_genes[["C14_GSE166555_xtoy"]][,1]
x20=top_20_genes[["C14_GSE166555_ytox"]][,1]
GSE166555_subtype4_C14 = union(x19,x20)
x21=top_20_genes[["C18_GSE166555_xtoy"]][,1]
x22=top_20_genes[["C18_GSE166555_ytox"]][,1]
GSE166555_subtype5_C18 = union(x21,x22)
x23=top_20_genes[["C29_GSE166555_xtoy"]][,1]
x24=top_20_genes[["C29_GSE166555_ytox"]][,1]
GSE166555_subtype6_C29 = union(x23,x24)

gene_symbol_data=c( "GSE166555_subtype1_C5","GSE166555_subtype2_C11","GSE166555_subtype3_C13",
                   "GSE166555_subtype4_C14","GSE166555_subtype5_C18","GSE166555_subtype6_C29")

kk_6subtypes=list()
plot_6subtypes=list()
for(i in gene_symbol_data) {
  ###genes
  genes_ = get(i)
  entrezIDs <- mget(genes_, org.Hs.egSYMBOL2EG, ifnotfound=NA)# mapIds(org.Hs.eg.db,keys = genes,keytype = "SYMBOL",column = "ENTREZID")
  genes <- as.character(entrezIDs)

  ###GO分析
  kk <- enrichGO(gene = genes,OrgDb = org.Hs.eg.db, pvalueCutoff =0.05, qvalueCutoff = 0.2, ont="all",readable =T)
  plot <- dotplot(kk, showCategory = 5,split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale='free')

  kk_6subtypes <- c(kk_6subtypes,list(kk))
  plot_6subtypes <- c(plot_6subtypes,list(plot))
  print(i)
  #if(i=="GSE146771_subtype2_C17"){break;}
}
#save(kk_6subtypes,plot_6subtypes,file="GO_plot_6_GSE146771.RData")
#load("GO_plot_6_GSE146771..RData")

p7 <- plot_6subtypes[[1]]
p8 <- plot_6subtypes[[2]]
p9 <- plot_6subtypes[[3]]
p10 <- plot_6subtypes[[4]]
p11 <- plot_6subtypes[[5]]
p12 <- plot_6subtypes[[6]]

#[3] merge plots
library(patchwork)
#1-6
combined_p1 <- (p1 | p2 | p3  ) / ( p4 | p5 | p6 ) + plot_annotation(tag_levels = 'a')
print(combined_p1)
ggsave(
  filename = "Combined_plots_dark.png",  # File name of the output image
  plot = combined_p1,                  # The ggplot/patchwork object to save
  width = 22,                          # Width of the image in inches
  height = 14,                         # Height of the image in inches
  dpi = 300                            # Resolution in DPI (dots per inch)
)
#7-12
combined_p2 <- (p7 | p8 | p9  ) /(p10 | p11 | p12  ) + plot_annotation(tag_levels =list( c("g", "h", "i", "j", "k", "l")))
print(combined_p2)
ggsave(
  filename = "Combined_plots_GO.png",  # File name of the output image
  plot = combined_p2,                  # The ggplot/patchwork object to save
  width = 22,                          # Width of the image in inches
  height = 14,                         # Height of the image in inches
  dpi = 300                            # Resolution in DPI (dots per inch)
)
# merge
library(png)
library(grid)
# read PNG
img1_raster <- png::readPNG("Combined_plots_dark.png")
img2_raster <- png::readPNG("Combined_plots_GO.png")
# to grid grob
grob1 <- grid::rasterGrob(img1_raster, interpolate = TRUE)
grob2 <- grid::rasterGrob(img2_raster, interpolate = TRUE)
library(cowplot)
#  grob to plot
plot1 <- ggdraw() + draw_grob(grob1)
plot2 <- ggdraw() + draw_grob(grob2)
combined_p3 <- plot1 / plot2
print(combined_p3)
ggsave(
  filename = "Combined_plots_dark_GO.pdf",  # File name of the output image
  plot = combined_p3,                  # The ggplot/patchwork object to save
  width = 22,                          # Width of the image in inches
  height = 28,                         # Height of the image in inches
  dpi = 300                            # Resolution in DPI (dots per inch)
)



