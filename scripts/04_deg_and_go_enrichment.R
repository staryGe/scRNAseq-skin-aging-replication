# ==============================================================================
# Script 04: DEG Analysis & GO Pathway Enrichment for Mac4 Subcluster
# Project: GSE336400 Replication (Mouse Skin Aging scRNA-seq)
# ==============================================================================

# 1. 加载依赖包 ----------------------------------------------------------------
library(Seurat)
library(tidyverse)
library(clusterProfiler)  # 用于 GO 富集分析
library(org.Mm.eg.db)     # 小鼠基因注释数据库
library(ggrepel)          # 用于火山图标签防重叠

# 确保目录存在
if (!dir.exists("plots")) dir.create("plots")
if (!dir.exists("processed_data")) dir.create("processed_data")

# 2. 读取 03 脚本保存的巨噬细胞数据 -------------------------------------------
message("--> 正在加载巨噬细胞数据集 (03_seu_macrophage.rds)...")
seu_mac <- readRDS("processed_data/03_seu_macrophage.rds")

# 设主 Idents 为 mac_subcluster
Idents(seu_mac) <- "mac_subcluster"

# 3. 寻找 Mac4 亚群的差异表达基因 (DEGs vs 其他所有 Mac 亚群) --------------------
message("--> 正在计算 Mac4 亚群的差异表达基因 (FindMarkers)...")
mac4_degs <- FindMarkers(
  seu_mac, 
  ident.1 = "Mac4", 
  min.pct = 0.25, 
  logfc.threshold = 0.25
)

# 整理 DEG 表格
mac4_degs <- mac4_degs %>% 
  rownames_to_column(var = "gene") %>% 
  mutate(
    change = case_when(
      p_val_adj < 0.05 & avg_log2FC > 0.5 ~ "Up",
      p_val_adj < 0.05 & avg_log2FC < -0.5 ~ "Down",
      TRUE ~ "NotSig"
    )
  )

# 保存 DEG 表格
write.csv(mac4_degs, file = "processed_data/04_Mac4_DEGs.csv", row.names = FALSE)
message("--> Mac4 差异基因表格已保存至: processed_data/04_Mac4_DEGs.csv")

# 4. 绘制 Mac4 差异基因火山图 --------------------------------------------------
top_genes <- mac4_degs %>% 
  filter(change == "Up") %>% 
  arrange(p_val_adj, desc(avg_log2FC)) %>% 
  head(10) %>% 
  pull(gene)

p_volcano <- ggplot(mac4_degs, aes(x = avg_log2FC, y = -log10(p_val_adj), color = change)) +
  geom_point(alpha = 0.6, size = 1.5) +
  scale_color_manual(values = c("Up" = "#d73027", "Down" = "#4575b4", "NotSig" = "grey70")) +
  geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed", color = "gray50") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "gray50") +
  geom_text_repel(
    data = filter(mac4_degs, gene %in% top_genes),
    aes(label = gene),
    size = 3.5,
    max.overlaps = 20
  ) +
  theme_classic() +
  labs(title = "Mac4 Differentially Expressed Genes (Volcano Plot)",
       x = "Average log2 Fold Change", y = "-log10 (Adjusted P-value)")

ggsave("plots/04_Mac4_DEG_Volcano.png", plot = p_volcano, width = 7, height = 6)

# 5. GO (Biological Process) 功能富集分析 --------------------------------------
message("--> 正在提取 Mac4 显著上调基因进行 GO 富集分析...")
up_genes <- mac4_degs %>% 
  filter(change == "Up") %>% 
  pull(gene)

# 将 Symbol 转换成 Entrez ID
gene_ids <- bitr(up_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db)

ego <- enrichGO(
  gene          = gene_ids$ENTREZID,
  OrgDb         = org.Mm.eg.db,
  keyType       = "ENTREZID",
  ont           = "BP",              # Biological Process
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.2
)

# 保存 GO 富集对象
saveRDS(ego, file = "processed_data/04_Mac4_GO_results.rds")

# 6. 绘制 GO 富集分析条形图 ---------------------------------------------------
if (!is.null(ego) && nrow(as.data.frame(ego)) > 0) {
  p_go <- dotplot(ego, showCategory = 12) + 
    ggtitle("GO Biological Process Enrichment for Mac4") +
    theme_bw()
  
  ggsave("plots/04_Mac4_GO_Enrichment.png", plot = p_go, width = 8, height = 6)
  message("--> GO 富集分析气泡图已保存至: plots/04_Mac4_GO_Enrichment.png")
} else {
  message("--> 未检测到显著富集的 GO 通路，请检查阈值设置。")
}

message("\n🎉🎉🎉 恭喜！整个复现流程（Script 01 ~ 04）已全部运行完毕！")