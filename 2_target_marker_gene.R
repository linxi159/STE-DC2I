
#####[1]load
#CRC_GSE146771. CRC_EMTAB8107. CRC_GSE166555
CRC_GSE146771_DE <- read.csv("CRC_GSE146771_Smartseq2_DE_cluster_14_17_9.csv", header = TRUE, sep = ",")
CRC_EMTAB8107_DE <- read.csv("CRC_EMTAB8107_DE_cluster_6_8_9.csv", header = TRUE, sep = ",")
CRC_GSE166555_DE <- read.csv("CRC_GSE166555_DE_cluster_5_11_13_14_18_29.csv", header = TRUE, sep = ",")

CRC_GSE146771_DE_C14 <- CRC_GSE146771_DE[1:1887,c(2,3)]#malignant cell cluster C14: Differentially expressed genes
CRC_GSE146771_DE_C14_ <- CRC_GSE146771_DE_C14[order(abs(CRC_GSE146771_DE_C14$log2FC),decreasing = TRUE), ]
CRC_GSE146771_DE_C17 <- CRC_GSE146771_DE[1888:4439,c(2,3)]#malignant cell cluster C17: Differentially expressed genes
CRC_GSE146771_DE_C17_ <- CRC_GSE146771_DE_C17[order(abs(CRC_GSE146771_DE_C17$log2FC),decreasing = TRUE), ]
CRC_GSE146771_DE_C9 <- CRC_GSE146771_DE[4440:6974,c(2,3)]#malignant cell cluster C9: Differentially expressed genes
CRC_GSE146771_DE_C9_ <- CRC_GSE146771_DE_C9[order(abs(CRC_GSE146771_DE_C9$log2FC),decreasing = TRUE), ]

CRC_EMTAB8107_DE_C6 <- CRC_EMTAB8107_DE[1:1670,c(2,3)]#malignant cell cluster C6: Differentially expressed genes
CRC_EMTAB8107_DE_C6_ <- CRC_EMTAB8107_DE_C6[order(abs(CRC_EMTAB8107_DE_C6$log2FC),decreasing = TRUE), ]
CRC_EMTAB8107_DE_C8 <- CRC_EMTAB8107_DE[1671:2984,c(2,3)]#malignant cell cluster C8: Differentially expressed genes
CRC_EMTAB8107_DE_C8_ <- CRC_EMTAB8107_DE_C8[order(abs(CRC_EMTAB8107_DE_C8$log2FC),decreasing = TRUE), ]
CRC_EMTAB8107_DE_C9 <- CRC_EMTAB8107_DE[2985:4483,c(2,3)]#malignant cell cluster C9: Differentially expressed genes
CRC_EMTAB8107_DE_C9_ <- CRC_EMTAB8107_DE_C9[order(abs(CRC_EMTAB8107_DE_C9$log2FC),decreasing = TRUE), ]

CRC_GSE166555_DE_C5 <- CRC_GSE166555_DE[7486:8832,c(2,3)]#malignant cell cluster C5: Differentially expressed genes
CRC_GSE166555_DE_C5_ <- CRC_GSE166555_DE_C5[order(abs(CRC_GSE166555_DE_C5$log2FC),decreasing = TRUE), ]
CRC_GSE166555_DE_C11 <- CRC_GSE166555_DE[1:1345,c(2,3)]#malignant cell cluster C11: Differentially expressed genes
CRC_GSE166555_DE_C11_ <- CRC_GSE166555_DE_C11[order(abs(CRC_GSE166555_DE_C11$log2FC),decreasing = TRUE), ]
CRC_GSE166555_DE_C13 <- CRC_GSE166555_DE[1346:2677,c(2,3)]#malignant cell cluster C13: Differentially expressed genes
CRC_GSE166555_DE_C13_ <- CRC_GSE166555_DE_C13[order(abs(CRC_GSE166555_DE_C13$log2FC),decreasing = TRUE), ]
CRC_GSE166555_DE_C14 <- CRC_GSE166555_DE[2678:4881,c(2,3)]#malignant cell cluster C14: Differentially expressed genes
CRC_GSE166555_DE_C14_ <- CRC_GSE166555_DE_C14[order(abs(CRC_GSE166555_DE_C14$log2FC),decreasing = TRUE), ]
CRC_GSE166555_DE_C18 <- CRC_GSE166555_DE[4882:6516,c(2,3)]#malignant cell cluster C18: Differentially expressed genes
CRC_GSE166555_DE_C18_ <- CRC_GSE166555_DE_C18[order(abs(CRC_GSE166555_DE_C18$log2FC),decreasing = TRUE), ]
CRC_GSE166555_DE_C29 <- CRC_GSE166555_DE[6517:7485,c(2,3)]#malignant cell cluster C29: Differentially expressed genes
CRC_GSE166555_DE_C29_ <- CRC_GSE166555_DE_C29[order(abs(CRC_GSE166555_DE_C29$log2FC),decreasing = TRUE), ]

#####[2]top-k gene
k=10
CRC_GSE146771_C14_topk  = CRC_GSE146771_DE_C14_[1:k,1]
CRC_GSE146771_C17_topk  = CRC_GSE146771_DE_C17_[1:k,1]
CRC_GSE146771_C9_topk  = CRC_GSE146771_DE_C9_[1:k,1]
CRC_EMTAB8107_C6_topk = CRC_EMTAB8107_DE_C6_[1:k,1]
CRC_EMTAB8107_C8_topk = CRC_EMTAB8107_DE_C8_[1:k,1]
CRC_EMTAB8107_C9_topk = CRC_EMTAB8107_DE_C9_[1:k,1]
CRC_GSE166555_C5_topk = CRC_GSE166555_DE_C5_[1:k,1]
CRC_GSE166555_C11_topk = CRC_GSE166555_DE_C11_[1:k,1]
CRC_GSE166555_C13_topk = CRC_GSE166555_DE_C13_[1:k,1]
CRC_GSE166555_C14_topk = CRC_GSE166555_DE_C14_[1:k,1]
CRC_GSE166555_C18_topk = CRC_GSE166555_DE_C18_[1:k,1]
CRC_GSE166555_C29_topk = CRC_GSE166555_DE_C29_[1:k,1]

#####[3]marker gene
CRC_Cell_marker_gene <- read.csv("Cell_marker_Human_Seq_CRC_gene.csv", header = TRUE, sep = ",")
CRC_mutation_top20gene_TCGA <- read.csv("CRC_mutation_top20gene_TCGA.csv", header = TRUE, sep = ",")
#并集
CRC_marker_gene = union(CRC_Cell_marker_gene[,1],CRC_mutation_top20gene_TCGA[,1])

#####[4]target marker gene
#CRC_GSE146771: target marker gene
tt = CRC_GSE146771_C14_topk
intersect(tt,CRC_Cell_marker_gene[,1]) #交集
intersect(tt,CRC_mutation_top20gene_TCGA[,1])
CRC_GSE146771_C14_target_marker_gene = intersect(tt,CRC_marker_gene)
tt = CRC_GSE146771_C17_topk
intersect(tt,CRC_Cell_marker_gene[,1])
intersect(tt,CRC_mutation_top20gene_TCGA[,1])
CRC_GSE146771_C17_target_marker_gene = intersect(tt,CRC_marker_gene)
tt = CRC_GSE146771_C9_topk
intersect(tt,CRC_Cell_marker_gene[,1])
intersect(tt,CRC_mutation_top20gene_TCGA[,1])
CRC_GSE146771_C9_target_marker_gene = intersect(tt,CRC_marker_gene)

#CRC_EMTAB8107: target marker gene
tt = CRC_EMTAB8107_C6_topk
intersect(tt,CRC_Cell_marker_gene[,1])
intersect(tt,CRC_mutation_top20gene_TCGA[,1])
CRC_EMTAB8107_C6_target_marker_gene = intersect(tt,CRC_marker_gene)
tt = CRC_EMTAB8107_C8_topk
intersect(tt,CRC_Cell_marker_gene[,1])
intersect(tt,CRC_mutation_top20gene_TCGA[,1])
CRC_EMTAB8107_C8_target_marker_gene = intersect(tt,CRC_marker_gene)
tt = CRC_EMTAB8107_C9_topk
intersect(tt,CRC_Cell_marker_gene[,1])
intersect(tt,CRC_mutation_top20gene_TCGA[,1])
CRC_EMTAB8107_C9_target_marker_gene = intersect(tt,CRC_marker_gene)

#CRC_GSE166555: target marker gene
tt = CRC_GSE166555_C5_topk
intersect(tt,CRC_Cell_marker_gene[,1])
intersect(tt,CRC_mutation_top20gene_TCGA[,1])
CRC_GSE166555_C5_target_marker_gene = intersect(tt,CRC_marker_gene)
tt = CRC_GSE166555_C11_topk
intersect(tt,CRC_Cell_marker_gene[,1])
intersect(tt,CRC_mutation_top20gene_TCGA[,1])
CRC_GSE166555_C11_target_marker_gene = intersect(tt,CRC_marker_gene)
tt = CRC_GSE166555_C13_topk
intersect(tt,CRC_Cell_marker_gene[,1])
intersect(tt,CRC_mutation_top20gene_TCGA[,1])
CRC_GSE166555_C13_target_marker_gene = intersect(tt,CRC_marker_gene)

tt = CRC_GSE166555_C14_topk
intersect(tt,CRC_Cell_marker_gene[,1])
intersect(tt,CRC_mutation_top20gene_TCGA[,1])
CRC_GSE166555_C14_target_marker_gene = intersect(tt,CRC_marker_gene)
tt = CRC_GSE166555_C18_topk
intersect(tt,CRC_Cell_marker_gene[,1])
intersect(tt,CRC_mutation_top20gene_TCGA[,1])
CRC_GSE166555_C18_target_marker_gene = intersect(tt,CRC_marker_gene)
tt = CRC_GSE166555_C29_topk
intersect(tt,CRC_Cell_marker_gene[,1])
intersect(tt,CRC_mutation_top20gene_TCGA[,1])
CRC_GSE166555_C29_target_marker_gene = intersect(tt,CRC_marker_gene)

#save data
save(CRC_GSE146771_C14_target_marker_gene,CRC_GSE146771_C17_target_marker_gene,CRC_GSE146771_C9_target_marker_gene,
     CRC_EMTAB8107_C6_target_marker_gene,CRC_EMTAB8107_C8_target_marker_gene,CRC_EMTAB8107_C9_target_marker_gene,
     CRC_GSE166555_C5_target_marker_gene,CRC_GSE166555_C11_target_marker_gene,CRC_GSE166555_C13_target_marker_gene,
     CRC_GSE166555_C14_target_marker_gene,CRC_GSE166555_C18_target_marker_gene,CRC_GSE166555_C29_target_marker_gene,
     file = "CRC_GSE146771_C14_17_9_EMTAB8107_C6_8_9_GSE166555_C5_11_13_14_18_29_tmg.RData")
load("CRC_GSE146771_C14_17_9_EMTAB8107_C6_8_9_GSE166555_C5_11_13_14_18_29_tmg.RData")



