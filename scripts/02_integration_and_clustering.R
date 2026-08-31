# ==============================================================================
# Script 02: Data Integration, Clustering, and Major Cell Type Annotation
# Project: GSE336400 Replication (Mouse Skin Aging scRNA-seq)
# ==============================================================================

# 1. 加载依赖包 ----------------------------------------------------------------
library(Seurat)
library(tidyverse)
library(harmony) # 用于多样本去批次整合

# 确保输出目录存在
if (!dir.exists("plots")) dir.create("plots")
if (!dir.exists("processed_data")) dir.create("processed_data")

# 2. 读取步骤 1 保存的质控后数据 ----------------------------------------------
message("--> 正在加载 QC 后的数据...")
seu_qc <- readRDS("processed_data/01_seu_qc.rds")

# 3. 标准化与高变基因选择 ------------------------------------------------------
message("--> 正在进行数据标准化与高变基因筛选...")
seu_qc <- NormalizeData(seu_qc, normalization.method = "LogNormalize", scale.factor = 10000)
seu_qc <- FindVariableFeatures(seu_qc, selection.method = "vst", nfeatures = 2000)

# 4. 数据 Scaling 与 PCA 降维 --------------------------------------------------
message("--> 正在执行 ScaleData 与 PCA 降维...")
seu_qc <- ScaleData(seu_qc, features = VariableFeatures(seu_qc))
seu_qc <- RunPCA(seu_qc, npcs = 30, verbose = FALSE)

# 5. 使用 Harmony 进行多样本去批次整合 ---------------------------------------
message("--> 正在运行 Harmony 去批次整合 (按样本 GSM_ID 整合)...")
seu_integrated <- RunHarmony(seu_qc, group.by.vars = "GSM_ID", dims.use = 1:20)

# 6. UMAP 降维与全局聚类 ------------------------------------------------------
message("--> 正在运行 UMAP 降维与图形聚类...")
seu_integrated <- RunUMAP(seu_integrated, reduction = "harmony", dims = 1:20)
seu_integrated <- FindNeighbors(seu_integrated, reduction = "harmony", dims = 1:20)

# resolution = 0.5 适合用于划分主要的大细胞类群
seu_integrated <- FindClusters(seu_integrated, resolution = 0.5)

# 7. 可视化与保存图件至 plots/ 目录 -------------------------------------------
# 图 A: 全局 Cluster 聚类图
p_umap_cluster <- DimPlot(seu_integrated, reduction = "umap", label = TRUE, pt.size = 0.3) + 
  ggtitle("Global Clusters (Harmony Integrated)")

# 图 B: 按年龄拆分的 UMAP 图（观察不同年龄段细胞的分布状态）
p_umap_age <- DimPlot(seu_integrated, reduction = "umap", group.by = "Age", split.by = "Age", pt.size = 0.2) + 
  ggtitle("Cell Distribution across Ages (8wk vs 26wk vs 52wk)")

# 图 C: 检查皮肤组织主要细胞类型的经典 Marker
# Ptprc: 免疫细胞 | Adgre1, Cd68: 巨噬细胞 | Col1a1: 成纤维细胞 | Krt5, Krt14: 角质形成细胞 | Pecam1: 内皮细胞
marker_genes <- c("Ptprc", "Adgre1", "Cd68", "Col1a1", "Dcn", "Krt5", "Krt14", "Pecam1")
p_markers <- DotPlot(seu_integrated, features = marker_genes) + 
  RotatedAxis() + 
  ggtitle("Major Cell Type Marker Gene Expression")

# 保存所有图片到 plots 文件夹
ggsave("plots/02_UMAP_global_clusters.png", plot = p_umap_cluster, width = 8, height = 6)
ggsave("plots/02_UMAP_split_by_age.png", plot = p_umap_age, width = 14, height = 5)
ggsave("plots/02_Marker_genes_dotplot.png", plot = p_markers, width = 10, height = 5)

message("--> 所有 UMAP 及 Marker 表达图已成功保存至 plots/ 目录！")

# 8. 保存整合与聚类后的 R 对象 -------------------------------------------------
saveRDS(seu_integrated, file = "processed_data/02_seu_integrated.rds")
message("\n--> 步骤 2 完成！RDS 文件已成功保存至: processed_data/02_seu_integrated.rds")

####################################手动印上标签###################################

# 1. 将画图的主 Active Identity 设回纯数字的 seurat_clusters
Idents(seu_integrated) <- "seurat_clusters"

# 2. 构造图例专用的标签映射（数字_细胞类型）
legend_labels <- c(
  "0"  = "0_Keratinocytes_1",
  "1"  = "1_Keratinocytes_2",
  "2"  = "2_Fibroblasts_1",
  "3"  = "3_Keratinocytes_3",
  "4"  = "4_Keratinocytes_4",
  "5"  = "5_Keratinocytes_5",
  "6"  = "6_Fibroblasts_2",
  "7"  = "7_Fibroblasts_3",
  "8"  = "8_Macrophage",      # 👈 核心巨噬细胞
  "9"  = "9_Immune_Other",
  "10" = "10_Endothelial",
  "11" = "11_Keratinocytes_6",
  "12" = "12_Keratinocytes_7",
  "13" = "13_Pericytes/SMC",
  "14" = "14_Myeloid_Other",
  "15" = "15_Keratinocytes_8",
  "16" = "16_Fibroblasts_4",
  "17" = "17_Keratinocytes_9"
)

# 3. 绘图：图中纯数字，右侧图例显示 编号_细胞类型
p_umap_perfect <- DimPlot(
  seu_integrated, 
  reduction = "umap", 
  label = TRUE,           # 图中央只印纯数字
  repel = TRUE, 
  pt.size = 0.3
) + 
  scale_color_discrete(labels = legend_labels) + 
  ggtitle("Global Cell Types (Harmony Integrated)")

# 4. 保存图片
ggsave("plots/02_UMAP_annotated.png", plot = p_umap_perfect, width = 10, height = 6)

# 5. 【修复处】使用 unname() 消除向量名称冲突，成功写入 meta.data
seu_integrated$cell_type <- unname(legend_labels[as.character(seu_integrated$seurat_clusters)])

# 6. 保存 RDS 文件
saveRDS(seu_integrated, file = "processed_data/02_seu_integrated.rds")

message("--> 完美解决！图片已更新至 plots/02_UMAP_annotated.png，RDS 文件已成功保存！")