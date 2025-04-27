# 加载包
library(edgeR)
library(pheatmap)

# 获取文件列表
file_list <- list.files(
  path = "tumor-transcriptome-demo",
  pattern = "\\.txt$",
  recursive = TRUE,
  full.names = TRUE
)

# 通过第一个文件获取基因信息
fc_colnames <- c("Geneid", "Chr", "Start", "End", "Strand", "Length", "Counts")


first_file <- read.delim(
  file_list[1],
  header = FALSE,
  skip = 2,
  col.names = fc_colnames
)
gene_ids <- first_file$Geneid
num_genes <- length(gene_ids)

# 初始化矩阵
count_matrix <- matrix(
  nrow = num_genes, 
  ncol = length(file_list),
  dimnames = list(gene_ids, NULL)
)

# 读取样本
for (i in seq_along(file_list)) {
  data <- read.delim(
    file_list[i],
    header = FALSE,
    skip = 2,
    col.names = fc_colnames
  )
  
 
  count_matrix[, i] <- data$Counts
}

# 设置列名（样本ID）
sample_ids <- sapply(file_list, function(f) {
  parts <- unlist(strsplit(f, .Platform$file.sep))
  cancer_type <- parts[length(parts)-1]
  sample_name <- tools::file_path_sans_ext(parts[length(parts)])
  paste(cancer_type, sample_name, sep = "_")
})
colnames(count_matrix) <- sample_ids

# 计算logCPM
dge <- DGEList(counts = count_matrix)
cpm <- edgeR::cpm(dge,log=F)
logcpm <- log10(cpm+1)

# 计算Z-score
z_scores <- (logcpm - rowMeans(logcpm))/apply(logcpm,1,sd)
z_scores[is.na(z_scores)] <- 0
z_scores[is.infinite(z_scores)] <- 0

# 按照癌症类型排序样本
cancer_types <- factor(sapply(colnames(z_scores), function(x) {
  strsplit(x, "_")[[1]][1]
}))
annotation_col <- data.frame(CancerType = cancer_types)
rownames(annotation_col) <- colnames(z_scores)
cancer_order <- c("COAD", "READ", "ESCA")
annotation_col$CancerType <- factor(annotation_col$CancerType, levels = cancer_order)
type_sorted <- order(annotation_col$CancerType)
z_sorted <- z_scores[, type_sorted]

# 截断
z_clipped <- z_scores
z_clipped[z_clipped > 2] <- 2
z_clipped[z_clipped < -2] <- -2

# 绘图
pheatmap(
  z_clipped,
  annotation_col = annotation_col,
  cluster_cols = FALSE,
  show_rownames = FALSE,
  show_colnames = FALSE,
  treeheight_row = 0,
  clustering_distance_rows = "euclidean",
  clustering_method = "complete",
  main = "Z-score of logCPM Expression",
  color = colorRampPalette(c("blue", "white", "red"))(100)
)


