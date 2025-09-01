
library(ggplot2)
library(ggrepel)
library(patchwork)

#[1] DE for CRC_GSE146771
data_ <- read.csv("CRC_GSE146771_Smartseq2_DE_cluster_14_17_9.csv")
# read data
head(data_)
c_14 = data_[1:1888,]
c_17 = data_[1889:4441,]
c_9 = data_[4442:6978,]
#c_14
c_dt = c_14
de_genes <- data.frame(
  gene = c_dt[,2],
  logFC = c_dt[,3],
  adjusted.pvalue = c_dt[,5]
)
# add: -log10(p-value)
de_genes$logadj.p <- -log10(de_genes$adjusted.pvalue)
data_c14 = de_genes
#c_17
c_dt = c_17
de_genes <- data.frame(
  gene = c_dt[,2],
  logFC = c_dt[,3],
  adjusted.pvalue = c_dt[,5]
)
# add: -log10(p-value)
de_genes$logadj.p <- -log10(de_genes$adjusted.pvalue)
data_c17 = de_genes
#c_9
c_dt = c_9
de_genes <- data.frame(
  gene = c_dt[,2],
  logFC = c_dt[,3],
  adjusted.pvalue = c_dt[,5]
)
# add: -log10(p-value)
de_genes$logadj.p <- -log10(de_genes$adjusted.pvalue)
data_c9 = de_genes
####logFC is 2

plot_de_fig <- function(dat) {
  #import gene id
  m=ggplot(data=dat, aes(x=logFC, y =-log10(adjusted.pvalue))) +
  # divide data into 4 parts,-log10(0.05)=1.30103,min(abs(data$logFC)=0.25
  geom_point(data=subset(data, data$logadj.p < 1.3), color="gray",size =1) +
  geom_point(data=subset(data, data$logadj.p > 1.3 & abs(data$logFC) < 0.25), color="gray",size =1) +
  geom_point(data=subset(data, data$logadj.p > 1.3 & logFC > 0.25), color="#F4A261",size =1) +
  geom_point(data=subset(data, data$logadj.p > 1.3 & logFC < -0.25), color="#2a9d8f",size =1) +
  ## add line
  geom_hline(yintercept = 1.3, linetype = "dashed", color = "#C6B1B2")+
  geom_vline(xintercept = c(-0.25, 0.25), linetype = "dashed", color = "#C6B1B2")+
  #set background white
  theme_bw()+
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        axis.text.x = element_text(size = 16),  # x-axis scale font size of the number
        axis.text.y = element_text(size = 16)   # y-axis scale font size of the number
        )+ 
  ## modify axis
  labs(x="Log2(Fold Change)",y="-Log10(adjusted p-value)")+
  ## remove caption 
  theme(legend.position='none')
  ## add text information gene of logFC>5
  #geom_text_repel(data=subset(data, abs(logFC) > 5), aes(label=id),color="black",alpha = 0.8)

  return (m)
}

# [2] DE for EMTAB8107
data_ <- read.csv("CRC_EMTAB8107_DE_cluster_6_8_9.csv")
# read data
head(data_)
c_6 = data_[1:1671,]
c_8 = data_[1672:2986,]
c_9 = data_[2987:4486,]
#c_6
c_dt = c_6
de_genes <- data.frame(
  gene = c_dt[,2],
  logFC = c_dt[,3],
  adjusted.pvalue = c_dt[,5]
)
# add: -log10(p-value)
de_genes$logadj.p <- -log10(de_genes$adjusted.pvalue)
data_c6 = de_genes
#c_8
c_dt = c_8
de_genes <- data.frame(
  gene = c_dt[,2],
  logFC = c_dt[,3],
  adjusted.pvalue = c_dt[,5]
)
# add: -log10(p-value)
de_genes$logadj.p <- -log10(de_genes$adjusted.pvalue)
data_c8 = de_genes
#c_9
c_dt = c_9
de_genes <- data.frame(
  gene = c_dt[,2],
  logFC = c_dt[,3],
  adjusted.pvalue = c_dt[,5]
)
# add: -log10(p-value)
de_genes$logadj.p <- -log10(de_genes$adjusted.pvalue)
data_c9 = de_genes

m1=plot_de_fig(dat=data_c14)
m2=plot_de_fig(dat=data_c17)
m3=plot_de_fig(dat=data_c9)
m4=plot_de_fig(dat=data_c6)
m5=plot_de_fig(dat=data_c8)
m6=plot_de_fig(dat=data_c9)

#[3] single cell pseudotime analysis
load("CRC_GSE146771_data_DE_C14_17_9.RData")
load("CRC_EMTAB8107_data_DE_C6_8_9.RData")
#filter cells QC: [cell count < 1, gene number < 10]
filted_cells_genes <- function(mat)
{
  # filter：[count < 1]
  row_zeros <- apply(mat, 1, function(x) sum(x == 0))
  a=length(colnames(mat))#
  #a
  #summary(row_zeros)
  #row_zeros
  row_zeros_sorted = sort(row_zeros,decreasing = TRUE)
  #row_zeros_sorted
  threshold = a-1
  filted_genes=names(row_zeros_sorted[row_zeros_sorted > threshold])
  #filted_genes
  #aa=mat[filted_genes,]
  ##filter cells QC: [gene number < 10]
  col_zeros <- apply(mat, 2, function(x) sum(x == 0))
  b=length(rownames(mat))#
  #b
  #summary(col_zeros)
  #col_zeros
  col_zeros_sorted = sort(col_zeros,decreasing = TRUE)
  #col_zeros_sorted
  threshold = b-10
  filted_cells=names(col_zeros_sorted[col_zeros_sorted > threshold])
  #filted_cells
  #bb=mat[filted_cells,]
  #save
  t=setdiff(rownames(mat),filted_genes)
  tt=mat[t,]
  return(tt)
}
CRC_GSE146771_data_DE_C14_filtered = filted_cells_genes(CRC_GSE146771_data_DE_C14)
CRC_GSE146771_data_DE_C17_filtered = filted_cells_genes(CRC_GSE146771_data_DE_C17)
CRC_GSE146771_data_DE_C9_filtered = filted_cells_genes(CRC_GSE146771_data_DE_C9)
rm(list = c("CRC_GSE146771_data_DE_C14","CRC_GSE146771_data_DE_C17","CRC_GSE146771_data_DE_C9"))
CRC_EMTAB8107_data_DE_C6_filtered = filted_cells_genes(CRC_EMTAB8107_data_DE_C6)
CRC_EMTAB8107_data_DE_C8_filtered = filted_cells_genes(CRC_EMTAB8107_data_DE_C8)
CRC_EMTAB8107_data_DE_C9_filtered = filted_cells_genes(CRC_EMTAB8107_data_DE_C9)
rm(list = c("CRC_EMTAB8107_data_DE_C6","CRC_EMTAB8107_data_DE_C8","CRC_EMTAB8107_data_DE_C9"))

#if (!requireNamespace("BiocManager", quietly = TRUE))
# install.packages("BiocManager")
#BiocManager::install("TSCAN")
#browseVignettes("TSCAN")
library(TSCAN)
pseudotime_cells <- function(mat)
{
  procdata <-preprocess(mat,cvcutoff = 0.1)#cvcutoff = 0.5,1.0
  lpsmclust <- exprmclust(procdata)
  lpsorder <- TSCANorder(lpsmclust)
  plot_ <- plotmclust(lpsmclust,show_cell_names = F)
  pseudotime_cell_order  <- lpsorder
  # sorted cells using the pseudotime results
  data_ <- mat[,pseudotime_cell_order]
  result = list(data_,plot_)
  return(result)
}
## ----CRC_GSE146771----------
#CRC_GSE146771_data_DE_C14_17_9_filtered
a=rownames(CRC_GSE146771_data_DE_C14_filtered)
b=rownames(CRC_GSE146771_data_DE_C17_filtered)
c=rownames(CRC_GSE146771_data_DE_C9_filtered)
intersection <- intersect(intersect(a, b), c)
matrix1=CRC_GSE146771_data_DE_C14_filtered[intersection,]
matrix2=CRC_GSE146771_data_DE_C17_filtered[intersection,]
matrix3=CRC_GSE146771_data_DE_C9_filtered[intersection,]
merged_matrix <- do.call(cbind, list(matrix1, matrix2, matrix3))
CRC_GSE146771_data_DE_C14_17_9_filtered = merged_matrix
#C14 17 9
result1 = pseudotime_cells(CRC_GSE146771_data_DE_C14_filtered)
CRC_GSE146771_data_DE_C14_filtered_pseudotime = result1[[1]]
print(result1[[2]])
result2 = pseudotime_cells(CRC_GSE146771_data_DE_C17_filtered)
CRC_GSE146771_data_DE_C17_filtered_pseudotime = result2[[1]]
print(result2[[2]])
result3 = pseudotime_cells(CRC_GSE146771_data_DE_C9_filtered)
CRC_GSE146771_data_DE_C9_filtered_pseudotime = result3[[1]]
print(result3[[2]])
#
m7=result1[[2]] # C14
m8=result2[[2]] # C17
m9=result3[[2]] # C9

## ----CRC_EMTAB8107----------
#CRC_EMTAB8107_data_DE_C6_8_9_filtered
a=rownames(CRC_EMTAB8107_data_DE_C6_filtered)
b=rownames(CRC_EMTAB8107_data_DE_C8_filtered)
c=rownames(CRC_EMTAB8107_data_DE_C9_filtered)
intersection <- intersect(intersect(a, b), c)
matrix1=CRC_EMTAB8107_data_DE_C6_filtered[intersection,]
matrix2=CRC_EMTAB8107_data_DE_C8_filtered[intersection,]
matrix3=CRC_EMTAB8107_data_DE_C9_filtered[intersection,]
merged_matrix <- do.call(cbind, list(matrix1, matrix2, matrix3))
CRC_EMTAB8107_data_DE_C6_8_9_filtered = merged_matrix
#C6 8 9
result1 = pseudotime_cells(CRC_EMTAB8107_data_DE_C6_filtered)
CRC_EMTAB8107_data_DE_C6_filtered_pseudotime = result1[[1]]
print(result1[[2]])
result2 = pseudotime_cells(CRC_EMTAB8107_data_DE_C8_filtered)
CRC_EMTAB8107_data_DE_C8_filtered_pseudotime = result2[[1]]
print(result2[[2]])
result3 = pseudotime_cells(CRC_EMTAB8107_data_DE_C9_filtered)
CRC_EMTAB8107_data_DE_C9_filtered_pseudotime = result3[[1]]
print(result3[[2]])
#
m10=result1[[2]] # C6
m11=result2[[2]] # C8
m12=result3[[2]] # C9

combined_plot <- (m1 + m2 + m3) / (m4 + m5 + m6) /  (m7 + m8 + m9) / (m10 + m11 + m12) + plot_annotation(tag_levels = 'a')
combined_plot
#save
ggsave(combined_plot,file="CRC_GSE146771_EMTAB8107_cluster_14_17_9_6_8_9_DE_filtered_pseudotime.pdf",width = 14,height = 17)#save pdf



