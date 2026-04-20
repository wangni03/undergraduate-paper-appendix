# 准备工作
rm(list = ls())
setwd('D:/毕设/mime')

library(Mime1)
library(survminer)
library(Cairo)
library(snowfall)
library(randomForestSRC)
library(survivalROC)

sfInit(parallel = TRUE, cpus = 8)
#genelist <- read.table("data/targets_common_new.txt", header = FALSE, sep = " ")
genelist = read.table("data/senescence_genelist.txt")
genelist <- unlist(genelist)


# 数据预处理
Dataset1 <- read.table("data/Mime_TCGA_LIHC.txt", header = TRUE)
Dataset2 <- read.table("data/Mime_CHCC.txt", header = TRUE)
Dataset3 <- read.table("data/Mime_ICGC_LIRI.txt", header = TRUE)

Dataset1 <- Dataset1[, c(1, 3, 2, 4:ncol(Dataset1))]
Dataset1 <- Dataset1[Dataset1[, 2] != 0, ]
Dataset2 <- Dataset2[Dataset2[, 2] != 0, ]
Dataset3 <- Dataset3[Dataset3[, 2] != 0, ]

list_train_vali_Data <- list(TCGA_LIHC = Dataset1,
                             CHCC = Dataset2,
                             ICGC_LIRI = Dataset3)

# 构建模型
res <- ML.Dev.Prog.Sig(train_data = list_train_vali_Data$TCGA_LIHC,
                       list_train_vali_Data = list_train_vali_Data,
                       unicox.filter.for.candi = T,
                       unicox_p_cutoff = 0.01,
                       candidate_genes = genelist,
                       mode = 'all',
                       nodesize =10,seed = 5201314)                       


# 绘制所有模型的C-index热图
CairoPDF('Heatmap_C_index.pdf', width = 7, height = 12)
cindex_dis_all(
  res,                                  # 模型结果对象
  validate_set = names(list_train_vali_Data)[-1],  # 验证集名称（除第一个之外）
  order = names(list_train_vali_Data),               # 排序顺序
  width = 0.35                         # 热图单元格宽度
)
dev.off()


# 绘制生存曲线
survplot <- vector("list", 3)
for (i in 1:3) {
  # 为每个数据集生成生存曲线
  print(
    survplot[[i]] <- rs_sur(
      res,
      model_name = "StepCox[forward] + Enet[α=0.3]",  # 使用的模型
      dataset = names(list_train_vali_Data)[i],   # 数据集名称
      median.line = "hv",                         # 中位生存线类型
      cutoff = 0.5,                              # 风险分组的截断值
      conf.int = TRUE,                           # 显示置信区间
      xlab = "Day",                             # X轴标签
      pval.coord = c(1000, 0.9)                  # p值位置
    )
  )
}
# 合并生存曲线图
CairoPDF('survival_curve.pdf', width = 18, height = 6)
aplot::plot_list(gglist = survplot, ncol = 3)
dev.off()

# 计算各模型在不同时间点的AUC
all.auc.1y <- cal_AUC_ml_res(res.by.ML.Dev.Prog.Sig = res,train_data = list_train_vali_Data[["TCGA_LIHC"]],
                             inputmatrix.list = list_train_vali_Data,mode = 'all',AUC_time = 1,
                             auc_cal_method="KM")
all.auc.2y <- cal_AUC_ml_res(res.by.ML.Dev.Prog.Sig = res,train_data = list_train_vali_Data[["TCGA_LIHC"]],
                             inputmatrix.list = list_train_vali_Data,mode = 'all',AUC_time = 2,
                             auc_cal_method="KM")
all.auc.3y <- cal_AUC_ml_res(res.by.ML.Dev.Prog.Sig = res,train_data = list_train_vali_Data[["TCGA_LIHC"]],
                             inputmatrix.list = list_train_vali_Data,mode = 'all',AUC_time = 3,
                             auc_cal_method="KM")

# 1年AUC绘图
CairoPDF('1year_AUC.pdf', width = 10, height = 12)
auc_dis_all(all.auc.1y, dataset = names(list_train_vali_Data), validate_set = names(list_train_vali_Data)[-1],
            order = names(list_train_vali_Data), width = 0.35, year = 1)
dev.off()

# 绘制最佳模型的ROC曲线（1年）
CairoPDF('1year_AUC_by_best_model.pdf', width = 10, height = 10)
roc_vis(all.auc.1y, model_name = "StepCox[forward] + Enet[α=0.3]", dataset = names(list_train_vali_Data),
        order = names(list_train_vali_Data), anno_position = c(0.65, 0.55), year = 1)
dev.off()

# 2年AUC绘图
CairoPDF('2year_AUC.pdf', width = 10, height = 12)
auc_dis_all(all.auc.2y, dataset = names(list_train_vali_Data), validate_set = names(list_train_vali_Data)[-1],
            order = names(list_train_vali_Data), width = 0.35, year = 2)
dev.off()

# 绘制最佳模型的ROC曲线（2年）
CairoPDF('2year_AUC_by_best_model.pdf', width = 10, height = 10)
roc_vis(all.auc.2y, model_name = "StepCox[forward] + Enet[α=0.3]", dataset = names(list_train_vali_Data),
        order = names(list_train_vali_Data), anno_position = c(0.65, 0.55), year = 2)
dev.off()

# 3年AUC绘图
CairoPDF('3year_AUC.pdf', width = 10, height = 12)
auc_dis_all(all.auc.3y, dataset = names(list_train_vali_Data), validate_set = names(list_train_vali_Data)[-1],
            order = names(list_train_vali_Data), width = 0.35, year = 3)
dev.off()

# 3年最佳模型ROC曲线
CairoPDF('3year_AUC_by_best_model.pdf', width = 10, height = 10)
roc_vis(all.auc.3y, model_name = "StepCox[forward] + Enet[α=0.3]", dataset = names(list_train_vali_Data),
        order = names(list_train_vali_Data), anno_position = c(0.65, 0.55), year = 3)
dev.off()



rs.glioma.lgg.gbm <- cal_RS_pre.prog.sig(use_your_own_collected_sig = F,type.sig = c('LGG','GBM','Glioma'),
                                        list_input_data = list_train_vali_Data)
pdf('HR_com.pdf',width=30,height=10)
HR_com(rs.glioma.lgg.gbm,
       res,
       model_name="StepCox[forward] + Enet[α=0.3]",
       dataset=names(list_train_vali_Data),
       type = "categorical")
dev.off()

# 与已发表模型比较
cc.glioma.lgg.gbm <- cal_cindex_pre.prog.sig(use_your_own_collected_sig = F,type.sig = c('Glioma','LGG','GBM'),
                                            list_input_data = list_train_vali_Data)

pdf('Cindex_comparison.pdf', width = 11, height = 12)
cindex_comp(cc.glioma.lgg.gbm,
            res,
            model_name="StepCox[forward] + Enet[α=0.3]",
            dataset=names(list_train_vali_Data))
dev.off()

#AUC比较
auc.glioma.lgg.gbm.1 <- cal_auc_pre.prog.sig(use_your_own_collected_sig = F,
                                            type.sig = c('Glioma','LGG','GBM'),
                                            list_input_data = list_train_vali_Data,AUC_time = 1,
                                            auc_cal_method = 'KM')

pdf('1year_AUC_comparison.pdf', width = 11, height = 12)
auc_comp(auc.glioma.lgg.gbm.1,
         all.auc.1y,
         model_name="StepCox[forward] + Enet[α=0.3]",
         dataset=names(list_train_vali_Data))
dev.off()

auc.glioma.lgg.gbm.2 <- cal_auc_pre.prog.sig(use_your_own_collected_sig = F,
                                            type.sig = c('Glioma','LGG','GBM'),
                                            list_input_data = list_train_vali_Data,AUC_time = 2,
                                            auc_cal_method = 'KM')
pdf('2year_AUC_comparison.pdf', width = 11, height = 12)
auc_comp(auc.glioma.lgg.gbm.2,
         all.auc.2y,
         model_name="StepCox[forward] + Enet[α=0.3]",
         dataset=names(list_train_vali_Data))
dev.off()

auc.glioma.lgg.gbm.3 <- cal_auc_pre.prog.sig(use_your_own_collected_sig = F,
                                            type.sig = c('Glioma','LGG','GBM'),
                                            list_input_data = list_train_vali_Data,AUC_time = 3,
                                            auc_cal_method = 'KM')
pdf('3year_AUC_comparison.pdf', width = 11, height = 12)
auc_comp(auc.glioma.lgg.gbm.3,
         all.auc.3y,
         model_name="StepCox[forward] + Enet[α=0.3]",
         dataset=names(list_train_vali_Data))
dev.off()


# res.feature.all <- ML.Corefeature.Prog.Screen(InputMatrix = list_train_vali_Data$TCGA_LIHC,
#                                             candidate_genes = genelist,
#                                             mode = "all_without_SVM",nodesize =15,seed = 1234)


# #通过不同方法过滤的基因的翻转图：
# pdf(file="core_feature_select.pdf",width=20,height=15)
# core_feature_select(res.feature.all)
# dev.off()

# # 绘制通过不同方法过滤的基因的排名
# pdf(file="core_feature_rank.pdf",width=10,height=15)
# core_feature_rank(res.feature.all, top=20)
# dev.off()

# #我们随机选择前两个基因来分析它们的相关性
# dataset_col<-c("#3182BDFF","#E6550DFF",'green')
# corplot <- list()
# for (i in c(1:3)) {
#   print(corplot[[i]]<-cor_plot(list_train_vali_Data[[i]],
#                                dataset=names(list_train_vali_Data)[i],
#                                color = dataset_col[i],
#                                feature1="PON1",
#                                feature2="IVD",
#                                method="pearson"))
# }
# pdf(file="core_feature_corplot.pdf",width=30,height=10)
# aplot::plot_list(gglist=corplot,ncol=3)
# dev.off()

# # 根据不同数据集中特定基因的中位表达水平绘制患者的生存曲线
# survplot <- vector("list",3) 
# for (i in c(1:3)) {
#   print(survplot[[i]]<-core_feature_sur("PON1", 
#                                         InputMatrix=list_train_vali_Data[[i]],
#                                         dataset = names(list_train_vali_Data)[i],
#                                         color=c("blue","red"),
#                                         median.line = "hv",
#                                         cutoff = 0.5,
#                                         conf.int = T,
#                                         xlab="Day",pval.coord=c(1000,0.9)))
# }
# pdf(file="survival_gene_plot.pdf",width=30,height=10)
# aplot::plot_list(gglist=survplot,ncol=3)
# dev.off()