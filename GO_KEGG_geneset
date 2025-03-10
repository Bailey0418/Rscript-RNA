#1.批量获得GO术语对应的基因集合
library(clusterProfiler)
library(org.Hs.eg.db) 
go_terms <- c("GO:0001525", "GO:0060055", "GO:0061042", "GO:0101023", "GO:2000772","GO:0090398","GO:2000774","GO:2000773") 
gene_list <- list()
for (go in go_terms) {
  go_genes <- AnnotationDbi::select(
    org.Hs.eg.db,
    keys = go,
    keytype = "GO",
    columns = "SYMBOL"
  )$SYMBOL
  
  gene_list[[go]] <- go_genes
}

unique_genes <- unique(unlist(gene_list))
cat("去重后的基因数量:", length(unique_genes), "\n")
head(unique_genes)
write.csv(unique_genes, "GO_unique_genes.csv", row.names = FALSE)

#2.批量获得KEGG通路对应的基因集合
library(clusterProfiler)
library(org.Hs.eg.db)
library(KEGGREST)
kegg_terms <- c("hsa04218","hsa04115","hsa04110")
gene_list_kegg <- list()

for (kegg in kegg_terms) {
  # 使用 keggREST 包从 KEGG 数据库获取基因
  kegg_genes <- KEGGREST::keggGet(kegg)[[1]]$GENE
  # 提取 KEGG 基因 ID
  kegg_gene_ids <- gsub(".*:(.*)", "\\1", kegg_genes)
  # 将 KEGG 基因 ID 转换为基因符号
  gene_symbols <- mapIds(org.Hs.eg.db, keys = kegg_gene_ids, keytype = "ENTREZID", column = "SYMBOL")
  # 移除 NA 值（那些没有找到对应符号的基因）
  gene_symbols <- na.omit(gene_symbols)
  # 将基因符号保存到列表中
  gene_list_kegg[[kegg]] <- gene_symbols
}

unique_kegg_genes <- unique(unlist(gene_list_kegg))
# 输出基因数和前几个基因
cat("去重后的基因数量:", length(unique_kegg_genes), "\n")
head(unique_kegg_genes)
write.csv(unique_kegg_genes, "KEGG_unique_genes.csv", row.names = FALSE)
