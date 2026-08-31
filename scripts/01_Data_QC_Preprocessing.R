# ==============================================================================
# Script 01: Data QC and Preprocessing
# Project: GSE336400 Replication (Mouse Skin Aging scRNA-seq)
# ==============================================================================

# 1. 加载依赖包 ----------------------------------------------------------------
library(Seurat)
library(tidyverse)

# 创建用于保存结果和中间数据的文件夹（如果不存在的话）
if (!dir.exists("processed_data")) dir.create("processed_data")
if (!dir.exists("results")) dir.create("results")

# 2. 批量读取数据并注入 Metadata -----------------------------------------------
sample_dirs <- list.dirs("data", recursive = FALSE, full.names = TRUE)
sample_dirs <- sample_dirs[grep("GSM", sample_dirs)]

seurat_list <- list()

for (dir in sample_dirs) {
  folder_name <- basename(dir)
  message("--> 正在读取样本: ", folder_name)
  
  # 读取 10x 格式三件套
  cts <- Read10X(data.dir = dir)
  
  # 创建 Seurat 对象 (初步过滤极低表达的基因和细胞)
  obj <- CreateSeuratObject(
    counts = cts,
    project = folder_name,
    min.cells = 3,
    min.features = 200
  )
  
  # 自动解析文件夹名称提取元数据信息
  # 文件夹格式样例: GSM9834280_8wk_ctrl
  name_parts <- unlist(strsplit(folder_name, "_"))
  obj$GSM_ID    <- name_parts[1]
  obj$Age       <- name_parts[2]
  obj$Treatment <- name_parts[3]
  obj$Group     <- paste(name_parts[2], name_parts[3], sep = "_")
  
  # 计算小鼠线粒体基因比例 (小鼠基因前缀为 ^mt-)
  obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^mt-")
  
  seurat_list[[folder_name]] <- obj
}

# 3. 合并 6 个样本为单一 Seurat 大对象 ------------------------------------------
seu_obj <- merge(
  x = seurat_list[[1]],
  y = seurat_list[2:length(seurat_list)],
  add.cell.ids = names(seurat_list)
)

message("\n=== 合并完成！原始数据概况 ===")
print(seu_obj)

# 4. 绘制 QC 质控图并保存 ------------------------------------------------------
qc_plot <- VlnPlot(
  seu_obj, 
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
  group.by = "Group", 
  ncol = 3, 
  pt.size = 0.1
)

# 将 QC 图保存至 results 目录
ggsave("plots/01_QC_vlnplot_before_filtering.png", plot = qc_plot, width = 12, height = 5)
message("--> 质控图已保存至: results/01_QC_vlnplot_before_filtering.png")

# 5. 执行细胞质量过滤 (QC Filtering) -------------------------------------------
# 标准单细胞小鼠质控阈值：
# - nFeature_RNA > 300 (过滤空液滴)
# - nFeature_RNA < 6000 (过滤潜在的双细胞 Doublets)
# - percent.mt < 15% (过滤破损/濒死细胞)

seu_qc <- subset(
  seu_obj,
  subset = nFeature_RNA > 300 & nFeature_RNA < 6000 & percent.mt < 15
)

# 输出过滤前后的细胞总数对比
n_before <- ncol(seu_obj)
n_after  <- ncol(seu_qc)
message(sprintf("\n=== 质量过滤结果 SUMMARY ==="))
message(sprintf("过滤前细胞总数: %d", n_before))
message(sprintf("过滤后细胞总数: %d", n_after))
message(sprintf("剔除细胞比例: %.2f%%", (n_before - n_after) / n_before * 100))

# 5.1 绘制并保存过滤后的 QC 质控图（用于对比） ---------------------------------
qc_plot_after <- VlnPlot(
  seu_qc, 
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
  group.by = "Group", 
  ncol = 3, 
  pt.size = 0.1
)

# 保存过滤后的质控图
ggsave("plots/01_QC_vlnplot_after_filtering.png", plot = qc_plot_after, width = 12, height = 5)
message("--> 过滤后的质控图已保存至: results/01_QC_vlnplot_after_filtering.png")


# 6. 保存过滤后的 R 数据对象 (.rds) --------------------------------------------
saveRDS(seu_qc, file = "processed_data/01_seu_qc.rds")
message("\n--> 步骤 1 完成！RDS 文件已成功保存至: processed_data/01_seu_qc.rds")