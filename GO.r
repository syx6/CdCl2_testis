library(ReactomePA)
library(tidyverse)
library(data.table)
library(org.Mm.eg.db)
library(clusterProfiler)
library(biomaRt)
library(enrichplot)


Undiff_CDH <- SPG_diff_CDH %>% dplyr::filter(celltype == 'Undifferentiated SSCs')

Undiff_CDH_gene_sort <- Undiff_CDH %>% 
  arrange(desc(avg_log2FC))

Undiff_CDH_gene_list_sort <- Undiff_CDH_gene_sort$avg_log2FC 
names(Undiff_CDH_gene_list_sort) <- Undiff_CDH_gene_sort$gene1
head(Undiff_CDH_gene_list_sort)


diff_CDH <- SPG_diff_CDH %>% dplyr::filter(celltype == 'Differentiated SSCs')

diff_CDH_gene_sort <- diff_CDH %>% 
  arrange(desc(avg_log2FC))

diff_CDH_gene_list_sort <- diff_CDH_gene_sort$avg_log2FC 
names(diff_CDH_gene_list_sort) <- diff_CDH_gene_sort$gene1
head(diff_CDH_gene_list_sort)

diff_CDH_gsea_go <- gseGO(
  geneList     = diff_CDH_gene_list_sort,
  OrgDb        = org.Mm.eg.db,
  keyType      = "SYMBOL",
  ont          = "BP",
  #nPerm        = 1000, #置换检验的次数，默认为1000
  minGSSize    = 10,
  maxGSSize    = 500,
  pvalueCutoff = 1,
  verbose      = FALSE,
  seed = FALSE, 
  by = "fgsea")


gseaplot2(Undiff_CDH_gsea_go, 
          title = Undiff_CDH_gsea_go$Description[944], geneSetID = 944)

gseaplot2(diff_CDH_gsea_go, 
          title = diff_CDH_gsea_go$Description[3682], geneSetID = 3682)

gseaplot2(diff_CDH_gsea_go, 
          title = diff_CDH_gsea_go$Description[173], geneSetID = 173)

Undiff_CDH_GO <- Undiff_CDH_gene_sort %>%
  dplyr::filter(sig != 'non') %>%
  dplyr::filter(celltype == 'Undifferentiated SSCs') %>%
  dplyr::pull(gene, name=avg_log2FC) %>% 
  bitr(fromType = "SYMBOL",
       toType = c("ENSEMBL", "ENTREZID"),
       OrgDb = org.Mm.eg.db) %>% 
  .$ENTREZID %>%
  enrichGO(gene = .,
           OrgDb = org.Mm.eg.db,
           ont = "BP",
           readable = TRUE,
           pvalueCutoff = 0.05,
           qvalueCutoff = 0.05)


Undiff_CDH_GO <- Undiff_CDH_gene_sort %>%
  dplyr::filter(sig != 'non') %>%
  dplyr::filter(celltype == 'Undifferentiated SSCs') %>%
  dplyr::pull(gene, name=avg_log2FC) %>% 
  bitr(fromType = "SYMBOL",
       toType = c("ENSEMBL", "ENTREZID"),
       OrgDb = org.Mm.eg.db) %>% 
  .$ENTREZID %>%
  enrichGO(gene = .,
           OrgDb = org.Mm.eg.db,
           ont = "BP",
           readable = TRUE,
           pvalueCutoff = 0.05,
           qvalueCutoff = 0.05)


