# ==============================================================================
# Script 03: Macrophage Subclustering & Mac4 Identification
# Project: GSE336400 Replication (Mouse Skin Aging scRNA-seq)
# ==============================================================================

# 1. 加载依赖包 ----------------------------------------------------------------
library(Seurat)
library(tidyverse)
library(harmony)

# 确保目录存在
if (!dir.exists("plots")) dir.create("plots")
if (!dir.exists("processed_data")) dir.create("processed_data")

# 2. 读取 02 脚本保存的大对象并提取巨噬细胞 ------------------------------------
message("--> 正在加载全局整合数据 (02_seu_integrated.rds)...")
seu_integrated <- readRDS("processed_data/02_seu_integrated.rds")

message("--> 正在提取 Cluster 8 (Macrophage) 细胞...")
# 提取 8 号 Cluster (巨噬细胞)
seu_mac <- subset(seu_integrated, idents = "8")

# 3. 对巨噬细胞子集重新跑 Standard Workflow -----------------------------------
# 注意：提取子集后，必须重新寻找高变基因、重新 Scale 和跑 PCA，以展现亚群间的微小差异
message("--> 正在对巨噬细胞子集进行重新标准化与 PCA 降维...")
seu_mac <- FindVariableFeatures(seu_mac, selection.method = "vst", nfeatures = 2000)
seu_mac <- ScaleData(seu_mac, features = VariableFeatures(seu_mac))
seu_mac <- RunPCA(seu_mac, npcs = 20, verbose = FALSE)

# 4. 重新运行 Harmony 去批次与二次 UMAP 降维 ----------------------------------
message("--> 正在对巨噬细胞子集运行 Harmony 与 UMAP 聚类...")
seu_mac <- RunHarmony(seu_mac, group.by.vars = "GSM_ID", dims.use = 1:15)
seu_mac <- RunUMAP(seu_mac, reduction = "harmony", dims = 1:15)
seu_mac <- FindNeighbors(seu_mac, reduction = "harmony", dims = 1:15)

# 使用较低的 resolution (如 0.3) 来切分出 3~5 个精细的 Macrophage 亚群
seu_mac <- FindClusters(seu_mac, resolution = 0.3)

# 5. 重命名 Mac 亚群为 Mac1, Mac2, Mac3, Mac4... -------------------------------
# 将纯数字 cluster 映射为标准的亚群命名
mac_cluster_names <- paste0("Mac", as.numeric(levels(seu_mac$seurat_clusters)) + 1)
names(mac_cluster_names) <- levels(seu_mac$seurat_clusters)
seu_mac <- RenameIdents(seu_mac, mac_cluster_names)
seu_mac$mac_subcluster <- Idents(seu_mac)

# 6. 出图 1: 巨噬细胞亚群 UMAP 图 ---------------------------------------------
p_mac_umap <- DimPlot(seu_mac, reduction = "umap", label = TRUE, pt.size = 1.0) + 
  ggtitle("Macrophage Subclusters (UMAP)")

ggsave("plots/03_Macrophage_UMAP_subclusters.png", plot = p_mac_umap, width = 7, height = 5)

# 7. 出图 2: Mac4 关键特异 Marker 表达气泡图 -----------------------------------
# 衰老/促炎巨噬细胞（Mac4）经典 Marker：S100a8, S100a9, Cd14, Il1b, Tnf
mac4_markers <- c("S100a8", "S100a9", "Cd14", "Il1b", "Tnf", "Adgre1", "Cd68")
p_mac4_bubble <- DotPlot(seu_mac, features = mac4_markers) + 
  RotatedAxis() + 
  ggtitle("Macrophage Subcluster Marker Expression (Focus on Mac4)")

ggsave("plots/03_Mac4_Marker_Bubble.png", plot = p_mac4_bubble, width = 8, height = 5)

# 8. 出图 3: 不同年龄段中巨噬细胞亚群的细胞比例堆叠图 ------------------------
# 验证 Mac4 是否在 52wk (老年组) 中显著富集
prop_data <- seu_mac@meta.data %>%
  group_by(Age, mac_subcluster) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(Age) %>%
  mutate(percentage = count / sum(count) * 100)

p_prop <- ggplot(prop_data, aes(x = Age, y = percentage, fill = mac_subcluster)) +
  geom_bar(stat = "identity", position = "stack", width = 0.6) +
  theme_classic() +
  labs(title = "Macrophage Subcluster Proportions Across Ages",
       x = "Age Group", y = "Percentage (%)", fill = "Subcluster") +
  scale_fill_brewer(palette = "Set2")

ggsave("plots/03_Mac4_Proportion_by_Age.png", plot = p_prop, width = 6, height = 5)

message("--> 3 张巨噬细胞亚群分析图件已成功保存至 plots/ 目录！")

# 9. 保存巨噬细胞专用 RDS 数据对象 --------------------------------------------
saveRDS(seu_mac, file = "processed_data/03_seu_macrophage.rds")
message("\n--> 步骤 3 完成！RDS 已成功保存至: processed_data/03_seu_macrophage.rds")

##################依据年龄分组，重新绘图######################################

# 1. 将 Age 字段设为有序因子 (Factor)
seu_mac$Age <- factor(seu_mac$Age, levels = c("8wk", "26wk", "52wk"))

# 2. 重新计算比例
prop_data <- seu_mac@meta.data %>%
  group_by(Age, mac_subcluster) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(Age) %>%
  mutate(percentage = count / sum(count) * 100)

# 3. 重新绘图
p_prop_ordered <- ggplot(prop_data, aes(x = Age, y = percentage, fill = mac_subcluster)) +
  geom_bar(stat = "identity", position = "stack", width = 0.5) +
  theme_classic() +
  labs(title = "Macrophage Subcluster Proportions Across Ages",
       x = "Age Group", y = "Percentage (%)", fill = "Subcluster") +
  scale_fill_brewer(palette = "Set2")

# 4. 覆盖保存
ggsave("plots/03_Mac4_Proportion_by_Age.png", plot = p_prop_ordered, width = 6, height = 5)
message("--> 横坐标已按 8wk -> 26wk -> 52wk 完美重新排序并保存！")