
## ---- enrichment analysis---
library("clusterProfiler")
library("org.Hs.eg.db")
library("DOSE")
library(ReactomePA)
library("ggplot2")
library("enrichR")
library("GOplot")
library("enrichplot")
#library(miRspongeR)

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("miRspongeR")
install.packages("miRspongeR")

#####(1) Data preparation
top_20_genes <- readRDS("top_20_genes.rds")

#output top20
x1=top_20_genes[["C14_GSE146771_xtoy"]][,1]
x2=top_20_genes[["C14_GSE146771_ytox"]][,1]
GSE146771_subtype1_C14 = union(x1,x2)
x3=top_20_genes[["C17_GSE146771_xtoy"]][,1]
x4=top_20_genes[["C17_GSE146771_ytox"]][,1]
GSE146771_subtype2_C17 = union(x3,x4)
x5=top_20_genes[["C9_GSE146771_xtoy"]][,1]
x6=top_20_genes[["C9_GSE146771_ytox"]][,1]
GSE146771_subtype3_C9 = union(x5,x6)

x7=top_20_genes[["C6_EMTAB8107_xtoy"]][,1]
x8=top_20_genes[["C6_EMTAB8107_ytox"]][,1]
EMTAB8107_subtype1_C6 = union(x7,x8)
x9=top_20_genes[["C8_EMTAB8107_xtoy"]][,1]
x10=top_20_genes[["C8_EMTAB8107_ytox"]][,1]
EMTAB8107_subtype2_C8 = union(x9,x10)
x11=top_20_genes[["C9_EMTAB8107_xtoy"]][,1]
x12=top_20_genes[["C9_EMTAB8107_ytox"]][,1]
EMTAB8107_subtype3_C9 = union(x11,x12)

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

gene_symbol_data=c("GSE146771_subtype1_C14","GSE146771_subtype2_C17","GSE146771_subtype3_C9",
                   "EMTAB8107_subtype1_C6","EMTAB8107_subtype2_C8","EMTAB8107_subtype3_C9",
                   "GSE166555_subtype1_C5","GSE166555_subtype2_C11","GSE166555_subtype3_C13",
                   "GSE166555_subtype4_C14","GSE166555_subtype5_C18","GSE166555_subtype6_C29")

kk_12subtypes=list()
plot_12subtypes=list()
for(i in gene_symbol_data) {
  ###genes
  genes_ = get(i)
  entrezIDs <- mget(genes_, org.Hs.egSYMBOL2EG, ifnotfound=NA)# mapIds(org.Hs.eg.db,keys = genes,keytype = "SYMBOL",column = "ENTREZID")
  genes <- as.character(entrezIDs)
  
  ###GO分析
  #kk <- enrichGO(gene = genes,OrgDb = org.Hs.eg.db, pvalueCutoff =0.05, qvalueCutoff = 0.2, ont="all",readable =T)
  #plot <- dotplot(kk, showCategory = 5,split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale='free')
  
  ###KEGG分析
  #kk <- enrichKEGG(gene = genes, organism = "hsa",keyType = "kegg", pvalueCutoff =0.05, qvalueCutoff =0.2)
  #plot <- barplot(kk, drop = TRUE, showCategory = 10)
  
  ###Reactome分析
  #kk <- enrichPathway(gene = genes,organism = "human",pvalueCutoff = 0.05,qvalueCutoff = 0.2)
  #plot <- barplot(kk, drop = TRUE, showCategory = 10)
  
  ###DO分析
  kk = enrichDO(gene = genes, pvalueCutoff=0.05)
  plot <- barplot(kk, drop = TRUE, showCategory = 10)
  
  #
  kk_12subtypes <- c(kk_12subtypes,list(kk))
  plot_12subtypes <- c(plot_12subtypes,list(plot))
  print(i)
  #if(i=="GSE146771_subtype2_C17"){break;}
}
#save(kk_12subtypes,plot_12subtypes,file="kk_GO_plot.RData")
#save(kk_12subtypes,plot_12subtypes,file="kk_KEGG_plot.RData")
#save(kk_12subtypes,plot_12subtypes,file="kk_Reactome_plot.RData")
save(kk_12subtypes,plot_12subtypes,file="kk_DO_plot.RData")

#(2) analysis plots for GO, KEGG, Reactome, DO
#GO analysis
load("kk_GO_plot.RData")
library(patchwork)
combined_cnet <- (plot_12subtypes[[1]] | plot_12subtypes[[2]] | plot_12subtypes[[3]]) /
  (plot_12subtypes[[4]] | plot_12subtypes[[5]] | plot_12subtypes[[6]]) + plot_annotation(tag_levels = 'a')
print(combined_cnet)
ggsave(filename = "Combined_plots_GO_1.pdf",  
       plot = combined_cnet,            
       width = 25,                       
       height =23,                      
       dpi = 300                     
       
)
combined_cnet <- (plot_12subtypes[[7]] | plot_12subtypes[[8]] | plot_12subtypes[[9]]) /
  (plot_12subtypes[[10]] | plot_12subtypes[[11]] | plot_12subtypes[[12]]) + plot_annotation(tag_levels = 'a')
print(combined_cnet)
ggsave(filename = "Combined_plots_GO_2.pdf", 
       plot = combined_cnet,            
       width = 18,  #18           
       height = 15, #11                   
       dpi = 300                  
)

#KEGG analysis
load("kk_KEGG_plot.RData")
# 1 3 6 7 8 9 10 11 
library(patchwork)
combined_cnet <- (plot_12subtypes[[1]] | plot_12subtypes[[3]] | plot_12subtypes[[6]]) /
  (plot_12subtypes[[7]] | plot_12subtypes[[8]] | plot_12subtypes[[9]]) / 
  (plot_12subtypes[[10]] | plot_12subtypes[[11]] ) + plot_annotation(tag_levels = 'a')
print(combined_cnet)
ggsave(filename = "Combined_plots_KEGG_1.pdf",
       plot = combined_cnet,       
       width = 25,                     
       height =23,                
       dpi = 300             
       
)

#Reactome analysis
load("kk_Reactome_plot.RData")
# 1 3 7 8 9 12 
combined_cnet <- (plot_12subtypes[[1]] | plot_12subtypes[[3]] | plot_12subtypes[[7]]) /
  (plot_12subtypes[[8]] | plot_12subtypes[[9]] | plot_12subtypes[[12]]) + plot_annotation(tag_levels = 'a')
#print(combined_cnet)
ggsave(filename = "combined_plots_Reactome_1.pdf",
       plot = combined_cnet,             
       width = 18,  #18                    
       height = 15, #11       
       dpi = 300                
)

#DO analysis
load("kk_DO_plot.RData")
# 1 2 3 5 6 7 9 10 11 
combined_cnet <- (plot_12subtypes[[1]] | plot_12subtypes[[2]] | plot_12subtypes[[3]]) /
  (plot_12subtypes[[5]] | plot_12subtypes[[6]] | plot_12subtypes[[7]]) /
  (plot_12subtypes[[9]] | plot_12subtypes[[10]] | plot_12subtypes[[11]]) + plot_annotation(tag_levels = 'a')
#print(combined_cnet)
ggsave(filename = "combined_plots_DO_1.pdf",  
       plot = combined_cnet,             
       width = 25,  #18                    
       height = 23, #11                     
       dpi = 300              
)


