
#######Inference of CRC driver genes based on dark causality#######
#####[1]load: Different CRC cancer subtypes (malignant cell subpopulations)
load("CRC_GSE146771_data_DE_C14_17_9.RData")
load("CRC_EMTAB8107_data_DE_C6_8_9.RData")
load("CRC_GSE166555_data_DE_C5_11_13_14_18_29.RData")

#Quality Control (QC) of Malignant Cell Populations: Assessing the number of 
#expressed cells and genes [Filtering criteria: Cells with <1 gene, Genes expressed in <10 cells]
filted_cells_genes <- function(mat)
{
#Check the expression of genes across different cells [Filtering criteria: Retain genes expressed in ≥1 cell]
row_zeros <- apply(mat, 1, function(x) sum(x == 0))
a=length(colnames(mat))# all cell numbers
#a
#summary(row_zeros)
#row_zeros
row_zeros_sorted = sort(row_zeros,decreasing = TRUE)
#row_zeros_sorted
threshold = a-1
filted_genes=names(row_zeros_sorted[row_zeros_sorted > threshold])
#filted_genes
#aa=mat[filted_genes,]
#Check the expression levels of different genes in cells [Filtering criteria: Retain genes detected in ≥10 cells]
col_zeros <- apply(mat, 2, function(x) sum(x == 0))
b=length(rownames(mat))# all gene numbers
#b
#summary(col_zeros)
#col_zeros
col_zeros_sorted = sort(col_zeros,decreasing = TRUE)
#col_zeros_sorted
threshold = b-10
filted_cells=names(col_zeros_sorted[col_zeros_sorted > threshold])
#filted_cells
#bb=mat[filted_cells,]
#save filtered gene
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
CRC_GSE166555_data_DE_C5_filtered = filted_cells_genes(CRC_GSE166555_data_DE_C5)
CRC_GSE166555_data_DE_C11_filtered = filted_cells_genes(CRC_GSE166555_data_DE_C11)
CRC_GSE166555_data_DE_C13_filtered = filted_cells_genes(CRC_GSE166555_data_DE_C13)
CRC_GSE166555_data_DE_C14_filtered = filted_cells_genes(CRC_GSE166555_data_DE_C14)
CRC_GSE166555_data_DE_C18_filtered = filted_cells_genes(CRC_GSE166555_data_DE_C18)
CRC_GSE166555_data_DE_C29_filtered = filted_cells_genes(CRC_GSE166555_data_DE_C29)
rm(list = c("CRC_GSE166555_data_DE_C5","CRC_GSE166555_data_DE_C11","CRC_GSE166555_data_DE_C13","CRC_GSE166555_data_DE_C14","CRC_GSE166555_data_DE_C18","CRC_GSE166555_data_DE_C29"))


#####[2]Cell Ordering: Single-cell Pseudotime Analysis
#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("TSCAN")
#browseVignettes("TSCAN")
#library(TSCAN)
#data(lpsdata)
#procdata <- preprocess(lpsdata)
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
result = pseudotime_cells(CRC_GSE146771_data_DE_C14_17_9_filtered)
print(result[[2]])
#save the overall figure
library(ggrepel)
library(ggplot2)
ggsave(result[[2]],file="1.pdf",width = 8,height = 6)#save pdf
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
#save 3 single figure
library(patchwork)
m1=result1[[2]] # C14
m2=result2[[2]] # C17
m3=result3[[2]] # C9
# Combine and display the three plots
combined_plot <- m1 + m2 + m3 + plot_annotation(tag_levels = 'A')
# display a composite figure
combined_plot
ggsave(combined_plot,file="2.pdf",width = 24,height = 6)# save pdf
#save data
save(CRC_GSE146771_data_DE_C14_filtered_pseudotime,CRC_GSE146771_data_DE_C17_filtered_pseudotime,CRC_GSE146771_data_DE_C9_filtered_pseudotime, file = "CRC_GSE146771_data_DE_C14_17_9_filtered_pseudotime.RData")
load("CRC_GSE146771_data_DE_C14_17_9_filtered_pseudotime.RData")

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
result = pseudotime_cells(CRC_EMTAB8107_data_DE_C6_8_9_filtered)
print(result[[2]])
#
library(ggrepel)
library(ggplot2)
ggsave(result[[2]],file="1.pdf",width = 8,height = 6)#pdf
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
library(patchwork)
m1=result1[[2]] # C14
m2=result2[[2]] # C17
m3=result3[[2]] # C9
# 
combined_plot <- m1 + m2 + m3 + plot_annotation(tag_levels = 'A')
# 
combined_plot
ggsave(combined_plot,file="2.pdf",width = 24,height = 6)#pdf
#
save(CRC_EMTAB8107_data_DE_C6_filtered_pseudotime,CRC_EMTAB8107_data_DE_C8_filtered_pseudotime,CRC_EMTAB8107_data_DE_C9_filtered_pseudotime, file = "CRC_EMTAB8107_data_DE_C6_8_9_filtered_pseudotime.RData")
load("CRC_EMTAB8107_data_DE_C6_8_9_filtered_pseudotime.RData")

## ----CRC_GSE166555----------
#CRC_GSE166555_data_DE_C5_11_13_14_18_29_filtered
a=rownames(CRC_GSE166555_data_DE_C5_filtered)
b=rownames(CRC_GSE166555_data_DE_C11_filtered)
c=rownames(CRC_GSE166555_data_DE_C13_filtered)
d=rownames(CRC_GSE166555_data_DE_C14_filtered)
e=rownames(CRC_GSE166555_data_DE_C18_filtered)
f=rownames(CRC_GSE166555_data_DE_C29_filtered)
intersection <- Reduce(intersect, list(a, b, c, d, e, f))
matrix1=CRC_GSE166555_data_DE_C5_filtered[intersection,]
matrix2=CRC_GSE166555_data_DE_C11_filtered[intersection,]
matrix3=CRC_GSE166555_data_DE_C13_filtered[intersection,]
matrix4=CRC_GSE166555_data_DE_C14_filtered[intersection,]
matrix5=CRC_GSE166555_data_DE_C18_filtered[intersection,]
matrix6=CRC_GSE166555_data_DE_C29_filtered[intersection,]
merged_matrix <- do.call(cbind, list(matrix1, matrix2, matrix3,matrix4, matrix5, matrix6))
CRC_GSE166555_data_DE_C5_11_13_14_18_29_filtered = merged_matrix
result = pseudotime_cells(CRC_GSE166555_data_DE_C5_11_13_14_18_29_filtered)
print(result[[2]])
#
library(ggrepel)
library(ggplot2)
ggsave(result[[2]],file="1.pdf",width = 8,height = 6)#
#C5 11 13 14 18 29
result1 = pseudotime_cells(CRC_GSE166555_data_DE_C5_filtered)
CRC_GSE166555_data_DE_C5_filtered_pseudotime = result1[[1]]
print(result1[[2]])
result2 = pseudotime_cells(CRC_GSE166555_data_DE_C11_filtered)
CRC_GSE166555_data_DE_C11_filtered_pseudotime = result2[[1]]
print(result2[[2]])
result3 = pseudotime_cells(CRC_GSE166555_data_DE_C13_filtered)
CRC_GSE166555_data_DE_C13_filtered_pseudotime = result3[[1]]
print(result3[[2]])
result4 = pseudotime_cells(CRC_GSE166555_data_DE_C14_filtered)
CRC_GSE166555_data_DE_C14_filtered_pseudotime = result4[[1]]
print(result4[[2]])
result5 = pseudotime_cells(CRC_GSE166555_data_DE_C18_filtered)
CRC_GSE166555_data_DE_C18_filtered_pseudotime = result5[[1]]
print(result5[[2]])
result6 = pseudotime_cells(CRC_GSE166555_data_DE_C29_filtered)
CRC_GSE166555_data_DE_C29_filtered_pseudotime = result6[[1]]
print(result6[[2]])
#
library(patchwork)
m1=result1[[2]] # C5
m2=result2[[2]] # C11
m3=result3[[2]] # C13
m4=result4[[2]] # C14
m5=result5[[2]] # C18
m6=result6[[2]] # C29
# 
combined_plot <- (m1 + m2 + m3)/(m4 + m5 + m6) + plot_annotation(tag_levels = 'A')
# 
combined_plot
ggsave(combined_plot,file="2.pdf",width = 24,height = 12)#pdf
#
save(CRC_GSE166555_data_DE_C5_filtered_pseudotime,CRC_GSE166555_data_DE_C11_filtered_pseudotime,
     CRC_GSE166555_data_DE_C13_filtered_pseudotime,CRC_GSE166555_data_DE_C14_filtered_pseudotime,
     CRC_GSE166555_data_DE_C18_filtered_pseudotime,CRC_GSE166555_data_DE_C29_filtered_pseudotime,
     file = "CRC_GSE166555_data_DE_C5_11_13_14_18_29_filtered_pseudotime.RData")
load("CRC_GSE166555_data_DE_C5_11_13_14_18_29_filtered_pseudotime.RData")


#####[3]Prior Information Fusion to identify target marker genes (tmg) for malignant cell subpopulations
#Differential Expression; Gene Mutation Information (Mutation Matrix, CNV, SNP); Cancer Cell Markers
load("CRC_GSE146771_C14_17_9_EMTAB8107_C6_8_9_GSE166555_C5_11_13_14_18_29_tmg.RData")


#####[4] Dark causal inference
## Function Implementation：import function from dark_causal_infer_fun.R
source("dark_causal_infer_fun.R")

#(1)CRC_GSE146771 ...
# test
load("CRC_GSE146771_data_DE_C14_17_9_filtered_pseudotime.RData")
load("CRC_GSE146771_C14_17_9_EMTAB8107_C6_8_9_GSE166555_C5_11_13_14_18_29_tmg.RData")
x_gene_ = CRC_GSE146771_data_DE_C9_filtered_pseudotime#CRC_GSE146771_data_DE_C17_filtered_pseudotime #CRC_GSE146771_data_DE_C14_filtered_pseudotime
y_gene_names_ = CRC_GSE146771_C9_target_marker_gene#CRC_GSE146771_C17_target_marker_gene #CRC_GSE146771_C14_target_marker_gene

#(2)CRC_EMTAB8107 ...
# test
load("CRC_EMTAB8107_data_DE_C6_8_9_filtered_pseudotime.RData")
load("CRC_GSE146771_C14_17_9_EMTAB8107_C6_8_9_GSE166555_C5_11_13_14_18_29_tmg.RData")
x_gene_ = CRC_EMTAB8107_data_DE_C6_filtered_pseudotime#
y_gene_names_ = CRC_EMTAB8107_C6_target_marker_gene

#(3)CRC_GSE166555 ...
# test
load("CRC_GSE166555_data_DE_C5_11_13_14_18_29_filtered_pseudotime.RData")
load("CRC_GSE146771_C14_17_9_EMTAB8107_C6_8_9_GSE166555_C5_11_13_14_18_29_tmg.RData")
x_gene_ = CRC_GSE166555_data_DE_C5_filtered_pseudotime#
y_gene_names_ = CRC_GSE166555_C5_target_marker_gene

#Optimal Parameters
#Randomly select 100 genes for testing to determine the optimal parameters E and tau; initial values set as E=3 and tau=2
set.seed(123)  # Set a random seed to ensure reproducibility of results

# Application:
Dark_causal_inference <- function(t=0) {
  # load data
  #CRC_GSE146771_data_DE_C9_filtered_pseudotime #CRC_GSE146771_data_DE_C17_filtered_pseudotime #CRC_GSE146771_data_DE_C14_filtered_pseudotime
  #CRC_EMTAB8107_data_DE_C9_filtered_pseudotime#CRC_EMTAB8107_data_DE_C8_filtered_pseudotime#CRC_EMTAB8107_data_DE_C6_filtered_pseudotime#
  # CRC_GSE166555_data_DE_C13_filtered_pseudotime#CRC_GSE166555_data_DE_C11_filtered_pseudotime#CRC_GSE166555_data_DE_C5_filtered_pseudotime#
  x_gene = CRC_GSE166555_data_DE_C29_filtered_pseudotime#CRC_GSE166555_data_DE_C18_filtered_pseudotime#CRC_GSE166555_data_DE_C14_filtered_pseudotime
  #CRC_GSE146771_C9_target_marker_gene #CRC_GSE146771_C17_target_marker_gene #CRC_GSE146771_C14_target_marker_gene
  #CRC_EMTAB8107_C9_target_marker_gene#CRC_EMTAB8107_C8_target_marker_gene#CRC_EMTAB8107_C6_target_marker_gene#
  #CRC_GSE166555_C13_target_marker_gene#CRC_GSE166555_C11_target_marker_gene#CRC_GSE166555_C5_target_marker_gene
  y_gene_names = CRC_GSE166555_C29_target_marker_gene#CRC_GSE166555_C18_target_marker_gene#CRC_GSE166555_C14_target_marker_gene

   # regulation strength with dark causality
  pnd_reg_str <- matrix(nrow=length(rownames(x_gene))*length(y_gene_names),ncol=5)
  colnames(pnd_reg_str) <- c("Reg_pair","Dark")
  Xi <-sequence(length(rownames(x_gene)))
  Yj <-sequence(length(y_gene_names))
  cnt <- 0
  for (j in Yj) {
    y_gene_name <- y_gene_names[j]
    y_gene_exp <- x_gene[y_gene_name,]
    for (i in Xi) {
      x_gene_name <- rownames(x_gene)[i]
      x_gene_exp <- x_gene[x_gene_name,]
      cnt <- cnt + 1

      # X->Y or Y->X regulation relationships
      reg_pair <- paste0(x_gene_name,"_reg_",y_gene_name)
      #reg_pair <- paste0(y_gene_name,"_reg_",x_gene_name)
      
      # 
      result <- dark_causality(input_X=x_gene_exp,input_Y=y_gene_exp, embedding_E = 3, delay_tau = 2, dist_metric = "euclidean", horizon_h = 1, is_weighted = TRUE, show_progress = FALSE)
      #result <- dark_causality(input_X=y_gene_exp,input_Y=x_gene_exp, embedding_E = 3, delay_tau = 2, dist_metric = "euclidean", horizon_h = 1, is_weighted = TRUE, show_progress = FALSE)
      
      d_reg_str[cnt,1] <- reg_pair
      d_reg_str[cnt,5] <- result[["dark"]]
      #if(cnt == 10){break;}
    }
    #if(cnt == 10){break;}
  }
  return(d_reg_str)
}
#Parallel Computing Execution
library(parallel)
detectCores(logical = F)  # 4
mc <- getOption("mc.cores", 4)
system.time({ D_causality_null_1 <- mclapply(list(1), Dark_causal_inference, mc.cores = mc);});
stopCluster(mc)


#####[5] Dark causal inference CRC Driver Gene Results and Comparative Analysis
# load saved varibles
#(1)CRC_GSE146771 ; (2)CRC_EMTAB8107; (3)CRC_GSE166555
#save
saveRDS(top_20_genes, file = "top_20_genes.rds")


#####[6] Inferenced CRC Driver Gene analysis
top_20_genes <- readRDS("top_20_genes.rds")
tmp_t20= c(top_20_genes[["C14_GSE146771_xtoy"]][1:5,1],top_20_genes[["C14_GSE146771_ytox"]][1:5,1],
  top_20_genes[["C17_GSE146771_xtoy"]][1:5,1],top_20_genes[["C17_GSE146771_ytox"]][1:5,1],
  top_20_genes[["C9_GSE146771_xtoy"]][1:5,1],top_20_genes[["C9_GSE146771_ytox"]][1:5,1],
  top_20_genes[["C6_EMTAB8107_xtoy"]][1:5,1],top_20_genes[["C6_EMTAB8107_ytox"]][1:5,1],
  top_20_genes[["C8_EMTAB8107_xtoy"]][1:5,1],top_20_genes[["C8_EMTAB8107_ytox"]][1:5,1],
  top_20_genes[["C9_EMTAB8107_xtoy"]][1:5,1],top_20_genes[["C9_EMTAB8107_ytox"]][1:5,1],
  top_20_genes[["C5_GSE166555_xtoy"]][1:5,1],top_20_genes[["C5_GSE166555_ytox"]][1:5,1],
  top_20_genes[["C11_GSE166555_xtoy"]][1:5,1],top_20_genes[["C11_GSE166555_ytox"]][1:5,1],
  top_20_genes[["C13_GSE166555_xtoy"]][1:5,1],top_20_genes[["C13_GSE166555_ytox"]][1:5,1],
  top_20_genes[["C14_GSE166555_xtoy"]][1:5,1],top_20_genes[["C14_GSE166555_ytox"]][1:5,1],
  top_20_genes[["C18_GSE166555_xtoy"]][1:5,1],top_20_genes[["C18_GSE166555_ytox"]][1:5,1],
  top_20_genes[["C29_GSE166555_xtoy"]][1:5,1],top_20_genes[["C29_GSE166555_ytox"]][1:5,1])

sort(table(tmp_t20))

#output top20 genes
x1=top_20_genes[["C14_GSE146771_xtoy"]][,1]
write.csv(x1, file = "C14_GSE146771_xtoy.csv")
x2=top_20_genes[["C14_GSE146771_ytox"]][,1]
write.csv(x2, file = "C14_GSE146771_ytox.csv")
x3=top_20_genes[["C17_GSE146771_xtoy"]][,1]
write.csv(x3, file = "C17_GSE146771_xtoy.csv")
x4=top_20_genes[["C17_GSE146771_ytox"]][,1]
write.csv(x4, file = "C17_GSE146771_ytox.csv")
x5=top_20_genes[["C9_GSE146771_xtoy"]][,1]
write.csv(x5, file = "C9_GSE146771_xtoy.csv")
x6=top_20_genes[["C9_GSE146771_ytox"]][,1]
write.csv(x6, file = "C9_GSE146771_ytox.csv")

x7=top_20_genes[["C6_EMTAB8107_xtoy"]][,1]
write.csv(x7, file = "C6_EMTAB8107_xtoy.csv")
x8=top_20_genes[["C6_EMTAB8107_ytox"]][,1]
write.csv(x8, file = "C6_EMTAB8107_ytox.csv")
x9=top_20_genes[["C8_EMTAB8107_xtoy"]][,1]
write.csv(x9, file = "C8_EMTAB8107_xtoy.csv")
x10=top_20_genes[["C8_EMTAB8107_ytox"]][,1]
write.csv(x10, file = "C8_EMTAB8107_ytox.csv")
x11=top_20_genes[["C9_EMTAB8107_xtoy"]][,1]
write.csv(x11, file = "C9_EMTAB8107_xtoy.csv")
x12=top_20_genes[["C9_EMTAB8107_ytox"]][,1]
write.csv(x12, file = "C9_EMTAB8107_ytox.csv")

x13=top_20_genes[["C5_GSE166555_xtoy"]][,1]
write.csv(x13, file = "C5_GSE166555_xtoy.csv")
x14=top_20_genes[["C5_GSE166555_ytox"]][,1]
write.csv(x14, file = "C5_GSE166555_ytox.csv")
x15=top_20_genes[["C11_GSE166555_xtoy"]][,1]
write.csv(x15, file = "C11_GSE166555_xtoy.csv")
x16=top_20_genes[["C11_GSE166555_ytox"]][,1]
write.csv(x16, file = "C11_GSE166555_ytox.csv")
x17=top_20_genes[["C13_GSE166555_xtoy"]][,1]
write.csv(x17, file = "C13_GSE166555_xtoy.csv")
x18=top_20_genes[["C13_GSE166555_ytox"]][,1]
write.csv(x18, file = "C13_GSE166555_ytox.csv")

x19=top_20_genes[["C14_GSE166555_xtoy"]][,1]
write.csv(x19, file = "C14_GSE166555_xtoy.csv")
x20=top_20_genes[["C14_GSE166555_ytox"]][,1]
write.csv(x20, file = "C14_GSE166555_ytox.csv")
x21=top_20_genes[["C18_GSE166555_xtoy"]][,1]
write.csv(x21, file = "C18_GSE166555_xtoy.csv")
x22=top_20_genes[["C18_GSE166555_ytox"]][,1]
write.csv(x22, file = "C18_GSE166555_ytox.csv")
x23=top_20_genes[["C29_GSE166555_xtoy"]][,1]
write.csv(x23, file = "C29_GSE166555_xtoy.csv")
x24=top_20_genes[["C29_GSE166555_ytox"]][,1]
write.csv(x24, file = "C29_GSE166555_ytox.csv")








