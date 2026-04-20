# 安装对应包
# install.packages(c("survival", "survminer"))

# 加载必要的包
library(survival)
library(survminer)
library(ggplot2)
library(dplyr)
library(tidyr)

setwd('D:/毕设/生存分析')
# 定义数据集名称和对应的文件路径
datasets <- c("TCGA_LIHC", "CHCC", "ICGC_LIRI")
file_paths <- c(
  "D:/毕设/mime/data/Mime_TCGA_LIHC.txt",
  "D:/毕设/mime/data/Mime_CHCC.txt",
  "D:/毕设/mime/data/Mime_ICGC_LIRI.txt"
)
names(file_paths) <- datasets

# 定义目标基因
gene_of_interest <- "COBLL1"

# 用于存储每个数据集的生存拟合结果
fit_list <- list()
data_list <- list()

for (ds in datasets) {
  cat("\nProcessing", ds, "...\n")
  
  # 读取数据（第一列为ID，作为行名；列名包含OS.time、OS及基因）
  data_raw <- read.table(file_paths[ds], header = TRUE, row.names = 1, check.names = FALSE)
    
  # 提取临床信息（生存时间与状态）
  clin <- data_raw[, c("OS.time", "OS")]
  
  # 提取基因表达矩阵：除 OS.time 和 OS 外的所有列
  expr <- data_raw[, !(colnames(data_raw) %in% c("OS.time", "OS"))]
  
    
  # 提取该基因的表达值（向量）
  gene_exp <- expr[[gene_of_interest]]
  
  # 构建生存数据框
  surv_data <- data.frame(
    OS.time = clin$OS.time,
    OS = clin$OS,
    exp = gene_exp,
    stringsAsFactors = FALSE
  )
  
  # 剔除无效生存数据
  surv_data <- surv_data[!is.na(surv_data$OS.time) & surv_data$OS.time > 0, ]
  surv_data <- surv_data[!is.na(surv_data$OS), ]
  
  if (nrow(surv_data) == 0) {
    warning(ds, " 无有效生存数据，跳过")
    next
  }
  
  cut <- surv_cutpoint(surv_data, time = "OS.time", event = "OS", variables = "exp")
  # 查看最佳切点
  summary(cut)
  # 生成包含分组的新数据框
  cut_data <- surv_categorize(cut)
  surv_data$group <- factor(cut_data$exp, levels = c("low", "high"))
  
  # 拟合KM曲线
  fit <- survfit(Surv(OS.time, OS) ~ group, data = surv_data)
  
  fit_list[[ds]] <- fit
  data_list[[ds]] <- surv_data
}

# 检查是否有成功处理的数据集
if (length(fit_list) == 0) {
  stop("没有成功处理任何数据集，请检查文件路径和基因名")
}


# 分别绘制每个数据集的生存曲线并保存为PDF
for (ds in names(fit_list)) {
  surv_plot <- ggsurvplot(
    fit_list[[ds]],
    data = data_list[[ds]],
    pval = TRUE,
    pval.coord = c(800, 0.15),        # 手动指定p值位置
    conf.int = TRUE,                   # 添加置信区间
    risk.table = TRUE,
    risk.table.height = 0.25,
    fontsize = 4,                      # risk table字体
    title = paste("Survival analysis of", gene_of_interest, "in", ds),
    xlab = "Survival Time (Days)",
    ylab = "Overall Survival Probability",
    legend.title = "Expression Groups",
    legend.labs = c("Low expression", "High expression"),
    palette = c("#2E9FDF", "#E7B800"),
    ggtheme = theme_bw(base_size = 14)
  )
  # 保存为PDF
  combined <- arrange_ggsurvplots(list(surv_plot), print = FALSE, ncol = 1, nrow = 1)
  pdf(paste0("survival_", ds, ".pdf"), width =8, height = 10)
  print(combined)
  dev.off()
}
