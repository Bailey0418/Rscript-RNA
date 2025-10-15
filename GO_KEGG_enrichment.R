setwd("")
library(clusterProfiler)
library(org.Rn.eg.db)
gene <- read.csv("de_deseq2_6vs6.csv")
gene <- gene[abs(gene$log2FoldChange)>0.5&gene$pvalue<0.05,]
#GO enrich
GO = enrichGO(gene$X, OrgDb = org.Rn.eg.db, keyType = 'SYMBOL', ont = "ALL", pvalueCutoff = 0.05, qvalueCutoff = 1)
write.csv(GO@result,"GO_enrichment_results.csv")
barplot(GO)
dotplot(GO)

#KEGG enrich
gene_KEGG <- bitr(gene$X,fromType = 'SYMBOL',toType = 'ENTREZID', OrgDb = org.Rn.eg.db, drop = T)
KEGG_database <- "rno"
KEGG = enrichKEGG(gene_KEGG$ENTREZID,organism = 'rno', pvalueCutoff = 0.05, qvalueCutoff = 1, use_internal_data = TRUE)
write.csv(KEGG@result,"KEGG_enrichment_results.csv")
barplot(KEGG)
dotplot(KEGG)

#可视化
library(enrichplot)
KEGG2 = pairwise_termsim(KEGG)
treeplot(GO2)
emapplot(GO2,showCategory = 15, layout = "kk")
cnetplot(GO2)
heatplot(KEGG2)
upsetplot(GO2)

#可视化
KEGG <- read.csv("KEGG_enrichment_results.csv")
library(ggplot2)
library(stringr)
library(forcats)
KEGG <- KEGG[KEGG$pvalue<0.05,]
KEGG <- KEGG[order(KEGG$FoldEnrichment, KEGG$Count, KEGG$pvalue), ]
#KEGG <- KEGG[c(1:20),]
p2 <- ggplot(KEGG, aes(x = FoldEnrichment, y = fct_reorder(subcategory, Count))) +
  geom_point(aes(size = Count, fill = pvalue), shape = 21, color = 'black') +  # 点的大小和颜色
  scale_fill_gradient(low = '#E27371', high = '#5D82A7') +  # 颜色渐变
  labs(title = 'KEGG Pathway Enrichment', y = 'KEGG Pathway', x = 'GeneRatio') +  # 标签
  guides(fill = guide_colorbar(reverse = TRUE, order = 1)) +  # 图例
  theme_bw()  # 主题
ggsave('KEGG_bubble.pdf', p2,width=12,height=9)

GO <- read.csv("GO_enrichment_results.csv")
top_GO <- GO[which(GO$pvalue<0.05),]
top_GO <- top_GO %>%
  group_by(ONTOLOGY) %>%
  slice_max(order_by = GeneRatio, n = 5, with_ties = FALSE) %>%
  ungroup() %>%
  mutate(Description = fct_reorder(Description, GeneRatio))
p2 <- ggplot(top_GO, aes(x = GeneRatio, y = fct_reorder(Description, Count)))+ #将term顺序按照GeneRatio进行排序
  geom_point(aes(size = Count, fill = pvalue), shape = 21, color = 'black')+ #分面
  facet_grid(ONTOLOGY~., scale = 'free_y', space = 'free_y')+
  scale_fill_gradient(low = '#E27371', high = '#5D82A7')+ #修改气泡图颜色
  labs(title = 'GO Enrichment', y = 'GO term', x = 'GeneRatio')+ #标题修改
  guides(fill = guide_colorbar(reverse = T,order = 1))+
  theme_bw()#主题
ggsave("GO_bubble.pdf",p2,width = 12,height = 9)
