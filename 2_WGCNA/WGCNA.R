library(dynamicTreeCut)
library(fastcluster)
library(WGCNA)
library(stringr)
library(ggplot2)
options(stringsAsFactors = FALSE) 
# 打开多线程
enableWGCNAThreads()
#加载数据
setwd('D:/毕设/wgcna')


exprMat <- "data/2_TCGA-LIHC_logCPM_matrix.txt"
type = "unsigned"
corType = "pearson"
corFnc = ifelse(corType=="pearson", cor, bicor)
maxPOutliers = ifelse(corType=="pearson",1,0.05)
# 关联样品性状的二元变量时，设置
robustY = ifelse(corType=="pearson",T,F)
##导入数据##
dataExpr <- read.table(exprMat, sep='\t', row.names=1, header=T, 
                     quote="", comment="", check.names=F)

dim(dataExpr)
head(dataExpr)[,1:8]

#数据筛选
## 筛选中位绝对偏差前75%的基因，至少MAD大于0.01
## 筛选后会降低运算量，也会失去部分信息
m.mad <- apply(dataExpr,1,mad)
dataExprVar <- dataExpr[which(m.mad > 
                 max(quantile(m.mad,na.rm=TRUE, probs=seq(0, 1, 0.25))[2],0.01)),]

#转换为样品在行，基因在列的矩阵
dataExpr <- as.data.frame(t(dataExprVar))

#检测缺失值
gsg = goodSamplesGenes(dataExpr, verbose = 3)

if (!gsg$allOK){
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes)>0) 
    printFlush(paste("Removing genes:", 
                     paste(names(dataExpr)[!gsg$goodGenes], collapse = ",")));
  if (sum(!gsg$goodSamples)>0) 
    printFlush(paste("Removing samples:", 
                     paste(rownames(dataExpr)[!gsg$goodSamples], collapse = ",")));
  # Remove the offending genes and samples from the data:
  dataExpr = dataExpr[gsg$goodSamples, gsg$goodGenes]
}

nGenes = ncol(dataExpr)
nSamples = nrow(dataExpr)

dim(dataExpr)
head(dataExpr)[,1:8]
#                                A1BG A1BG-AS1  A1CF    A2M A2M-AS1  A2ML1 A4GALT
# TCGA-DD-AACI-01A-11R-A41C-07 -2.116   -0.505 6.332 11.610   1.938 -4.602  3.577
# TCGA-UB-AA0V-01A-11R-A38B-07  1.431   -0.774 8.030 12.591   1.623 -3.249  1.832
# TCGA-G3-AAV3-01A-11R-A37K-07  0.010   -0.746 7.403  7.338   0.699 -4.602  2.903
# TCGA-EP-A3RK-01A-11R-A22L-07 -2.559   -1.501 6.300 11.650   0.485 -4.602  3.966
# TCGA-2Y-A9H9-01A-21R-A39D-07  0.400    0.335 7.380 11.643   0.867 -4.602  1.955
# TCGA-MI-A75H-01A-11R-A32O-07  0.244   -0.644 7.905  8.874   0.401 -4.602  2.605
#                               A4GNT
# TCGA-DD-AACI-01A-11R-A41C-07 -4.602
# TCGA-UB-AA0V-01A-11R-A38B-07 -1.916
# TCGA-G3-AAV3-01A-11R-A37K-07 -4.018
# TCGA-EP-A3RK-01A-11R-A22L-07 -4.602
# TCGA-2Y-A9H9-01A-21R-A39D-07 -3.645
# TCGA-MI-A75H-01A-11R-A32O-07 -4.602

#软阈值筛选
## 查看是否有离群样品
sampleTree = hclust(dist(dataExpr), method = "average")
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="")

#构建无标度网络
powers = c(c(1:10), seq(from = 12, to=30, by=2))
sft = pickSoftThreshold(dataExpr, powerVector=powers, 
                        networkType=type, verbose=5)

# pickSoftThreshold: will use block size 2225.
#  pickSoftThreshold: calculating connectivity for given powers...
#    ..working on genes 1 through 2225 of 20101
#    ..working on genes 2226 through 4450 of 20101
#    ..working on genes 4451 through 6675 of 20101
#    ..working on genes 6676 through 8900 of 20101
#    ..working on genes 8901 through 11125 of 20101
#    ..working on genes 11126 through 13350 of 20101
#    ..working on genes 13351 through 15575 of 20101
#    ..working on genes 15576 through 17800 of 20101
#    ..working on genes 17801 through 20025 of 20101
#    ..working on genes 20026 through 20101 of 20101
#    Power SFT.R.sq  slope truncated.R.sq  mean.k. median.k.  max.k.
# 1      1  0.00156  0.108          0.944 2.63e+03  2.59e+03 4720.00
# 2      2  0.49000 -1.830          0.935 5.76e+02  5.19e+02 1840.00
# 3      3  0.72100 -2.370          0.974 1.70e+02  1.31e+02  905.00
# 4      4  0.77600 -2.520          0.981 6.20e+01  3.83e+01  510.00
# 5      5  0.78900 -2.500          0.986 2.65e+01  1.25e+01  312.00
# 6      6  0.80100 -2.360          0.983 1.29e+01  4.51e+00  203.00
# 7      7  0.80500 -2.180          0.977 6.92e+00  1.75e+00  138.00
# 8      8  0.82200 -1.970          0.966 4.03e+00  7.29e-01   96.50
# 9      9  0.91900 -1.700          0.995 2.50e+00  3.27e-01   69.60
# 10    10  0.92000 -1.680          0.994 1.64e+00  1.53e-01   57.10
# 11    12  0.92100 -1.610          0.994 7.97e-01  3.69e-02   40.20
# 12    14  0.92100 -1.560          0.995 4.37e-01  1.00e-02   29.40
# 13    16  0.92700 -1.520          0.996 2.61e-01  2.94e-03   22.20
# 14    18  0.92700 -1.480          0.992 1.66e-01  9.11e-04   17.00
# 15    20  0.92300 -1.450          0.990 1.11e-01  2.88e-04   13.30
# 16    22  0.90500 -1.450          0.964 7.68e-02  9.60e-05   10.60
# 17    24  0.91500 -1.420          0.963 5.47e-02  3.23e-05    8.50
# 18    26  0.40400 -1.780          0.333 3.99e-02  1.11e-05    6.90
# 19    28  0.40900 -1.760          0.335 2.98e-02  3.91e-06    5.72
# 20    30  0.41100 -1.720          0.324 2.26e-02  1.39e-06    4.82

power = sft$powerEstimate
power
# [1] 9

# 无满足条件的power时选用
# 无向网络在power小于15或有向网络power小于30内，没有一个power值可以使
# 无标度网络图谱结构R^2达到0.8，平均连接度较高如在100以上，可能是由于
# 部分样品与其他样品差别太大。这可能由批次效应、样品异质性或实验条件对
# 表达影响太大等造成。可以通过绘制样品聚类查看分组信息和有无异常样品。
# 如果这确实是由有意义的生物变化引起的，也可以使用下面的经验power值。
if (is.na(power)){
  power = ifelse(nSamples<20, ifelse(type == "unsigned", 9, 18),
          ifelse(nSamples<30, ifelse(type == "unsigned", 8, 16),
          ifelse(nSamples<40, ifelse(type == "unsigned", 7, 14),
          ifelse(type == "unsigned", 6, 12))       
          )
          )
}



##一步法网络构建：One-step network construction and module detection##
# power: 上一步计算的软阈值
# maxBlockSize: 计算机能处理的最大模块的基因数量 (默认5000)；
# corType: pearson or bicor
# numericLabels: 返回数字是颜色作为模块的名字，后面可以再转换为数字
# saveTOMs：最耗费时间的计算，存储起来，供后续使用
# mergeCutHeight: 合并模块的阈值，越大模块越少
net = blockwiseModules(dataExpr, power = power, maxBlockSize = nGenes,
                       TOMType = type, minModuleSize = 30,
                       reassignThreshold = 0, mergeCutHeight = 0.25,
                       numericLabels = FALSE, pamRespectsDendro = FALSE,
                       saveTOMs=TRUE, corType = corType, 
                       maxPOutliers=maxPOutliers, loadTOMs=TRUE,
                       saveTOMFileBase = paste0(exprMat, ".tom"),
                       verbose = 3)
table(net$colors)

## 灰色的为**未分类**到模块的基因。
# Convert labels to colors for plotting
moduleLabels = net$colors
moduleColors = labels2colors(moduleLabels)
# Plot the dendrogram and the module colors underneath
pdf(file="merged_dynamic.pdf",width=10,height=10) 
plotDendroAndColors(net$dendrograms[[1]], moduleColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()

#关联表型数据
trait <- "data/TCGA-LIHC_clinical.txt"
# 读入表型数据，不是必须的
if(trait != "") {
  traitData <- read.table(file=trait, sep='\t', header=T, row.names=1,
                          check.names=FALSE, comment='',quote="")
  traitData$ajcc_pathologic_stage <- factor(traitData$ajcc_pathologic_stage, 
                   levels = c("Stage I", "Stage II",c("Stage IIIA", "Stage IIIB", "Stage III", "Stage IIIC"),c("Stage IVA", "Stage IVB", "Stage IV")),
                   labels = c(1, 2, 3, 3, 3, 3, 4, 4, 4))
  traitData$ajcc_pathologic_t <- factor(traitData$ajcc_pathologic_t,
                   levels = c("T1", "T2a","T2b","T2","T3a","T3b","T3", "T4"),
                   labels = c(1, 2, 2, 2, 3, 3, 3, 4))
  # traitData$tumor_grade <- factor(traitData$atumor_grade,
  #                  levels = c("G1", "G2","G3","G4"),
  #                  labels = c(1, 2, 3, 4))
  traitData$gender <- factor(traitData$gender,
                   levels = c("female", "male"),
                   labels = c(0,1))
  sampleName = rownames(dataExpr)
  traitData = traitData[match(str_sub(sampleName,1,12), rownames(traitData)), ]
}


### 模块与表型数据关联
MEs = net$MEs
MEs_col = MEs
colnames(MEs_col) = paste0("ME", labels2colors(
  as.numeric(str_replace_all(colnames(MEs),"ME",""))))
MEs_col = orderMEs(MEs_col)

if (corType=="pearson") {
  modTraitCor = cor(MEs_col, traitData, use = 'p')
  modTraitP = corPvalueStudent(modTraitCor, nSamples)
} else {
  modTraitCorP = bicorAndPvalue(MEs_col, traitData, robustY=robustY)
  modTraitCor = modTraitCorP$bicor
  modTraitP   = modTraitCorP$p
}

# signif表示保留几位小数
pdf(file="Module-trait_relationships.pdf",width=10,height=10) 
textMatrix = paste(signif(modTraitCor, 2), "\n(", signif(modTraitP, 1), ")", sep = "")
dim(textMatrix) = dim(modTraitCor)
par(mar = c(6, 8.5, 3, 3))
labeledHeatmap(Matrix = modTraitCor, xLabels = colnames(traitData), 
               yLabels = colnames(MEs_col), 
               cex.lab = 0.5, 
               ySymbols = colnames(MEs_col), colorLabels = FALSE, 
               colors = blueWhiteRed(50), 
               textMatrix = textMatrix, setStdMargins = FALSE, 
               cex.text = 0.5, zlim = c(-1,1),
               main = paste("Module-trait relationships"))
dev.off()


## 模块内基因与表型数据关联

# 性状跟模块虽然求出了相关性，可以挑选最相关的那些模块来分析，
# 但是模块本身仍然包含非常多的基因，还需进一步的寻找最重要的基因。
# 所有的模块都可以跟基因算出相关系数，所有的连续型性状也可以跟基因的表达值算出相关系数。

### 计算模块与基因的相关性矩阵

if (corType=="pearson") {
  geneModuleMembership = as.data.frame(cor(dataExpr, MEs_col, use = "p"))
  MMPvalue = as.data.frame(corPvalueStudent(
             as.matrix(geneModuleMembership), nSamples))
} else {
  geneModuleMembershipA = bicorAndPvalue(dataExpr, MEs_col, robustY=robustY)
  geneModuleMembership = geneModuleMembershipA$bicor
  MMPvalue   = geneModuleMembershipA$p
}


# 计算性状与基因的相关性矩阵

## 只有连续型性状才能进行计算，如果是离散变量，在构建样品表时就转为0-1矩阵。

if (corType=="pearson") {
  geneTraitCor = as.data.frame(cor(dataExpr, traitData, use = "p"))
  geneTraitP = as.data.frame(corPvalueStudent(
             as.matrix(geneTraitCor), nSamples))
} else {
  geneTraitCorA = bicorAndPvalue(dataExpr, traitData, robustY=robustY)
  geneTraitCor = as.data.frame(geneTraitCorA$bicor)
  geneTraitP   = as.data.frame(geneTraitCorA$p)
}

# 最后把两个相关性矩阵联合起来,指定感兴趣模块进行分析
nGenes = ncol(dataExpr)
nSamples = nrow(dataExpr)



#  ajcc_pathologic_stage
ajcc_pathologic_stage = as.data.frame(traitData$ajcc_pathologic_stage)
names(ajcc_pathologic_stage) = "ajcc_pathologic_stage"
modNames = substring(names(MEs), 3)
geneModuleMembership = as.data.frame(cor(dataExpr, MEs, use = "p"))
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership 
 ), nSamples))
names(geneModuleMembership) = paste("MM", modNames, sep="")
names(MMPvalue) = paste("p.MM", modNames, sep="")
# 计算性状和基因表达量之间的相关性（GS）
geneTraitSignificance = as.data.frame(cor(dataExpr, ajcc_pathologic_stage, use = "p"))
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), 
  nSamples))
names(geneTraitSignificance) = paste("GS.", names(ajcc_pathologic_stage), sep="")
names(GSPvalue) = paste("p.GS.", names(ajcc_pathologic_stage), sep="")


#yellow 模块
module = "4" #4:yellow;6:red
column = match(module, modNames)
module='yellow'
moduleGenes = moduleColors==module
table(moduleGenes)
yellow_module<-as.data.frame(dimnames(data.frame(dataExpr))[[2]][moduleGenes]) 
names(yellow_module)="genename"

#moduleGenes
#FALSE  TRUE 
#19086  1015 


#筛选hub基因
MM<-abs(geneModuleMembership[moduleGenes,column])
GS<-abs(geneTraitSignificance[moduleGenes, 1])
c_1<-as.data.frame(cbind(MM,GS))
rownames(c_1)=yellow_module$genename
head(c_1)

#   MM         GS
# ABCA9  0.5207030 0.23055613
# ABHD3  0.3930192 0.04355719
# ABLIM3 0.3976121 0.02815002
# ACAD8  0.4048236 0.13664815
# ACBD7  0.4103006 0.07397600
# ACER1  0.3884292 0.08782678


c_1$group <- abs(c_1$MM)>0.6&abs(c_1$GS)>0.15
write.csv(c_1, "hubgene_MMGS_yellow_ajcc_pathologic_stage.csv")


pdf("MM vs GS_yellow_ajcc_pathologic_stage.pdf",width = 7,height = 7)
ggplot(data=c_1, aes(x=MM, y=GS, color=group))+geom_point(size=1.5)+scale_colour_manual(values=c("grey60", "red"))+ theme_bw()+  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+  labs(x="Module Membership in yellow module", y="Gene significance for ajcc_pathologic_stage",title = "Module membership vs. gene significance ")+theme(axis.title.x =element_text(size=14), axis.title.y=element_text(size=14),axis.text = element_text(size = 12),axis.text.x = element_text(colour = "black"),axis.text.y = element_text(colour = "black"),plot.title = element_text(hjust = 0.5,size = 16,face = "bold"),plot.margin = unit(rep(2,4),'lines')) +theme(legend.position = 'none')+geom_hline(aes(yintercept=0.15),colour="#5B9BD5",lwd=1,linetype=1)+geom_vline(aes(xintercept=0.6),colour="#5B9BD5",lwd=1,linetype=1)
dev.off()


#red 模块
module = "6" #4:yellow;6:red
column = match(module, modNames)
module='red'
moduleGenes = moduleColors==module
table(moduleGenes)

# moduleGenes
# FALSE  TRUE 
# 19292   809 

red_module<-as.data.frame(dimnames(data.frame(dataExpr))[[2]][moduleGenes]) 
names(red_module)="genename"

#筛选hub基因
MM<-abs(geneModuleMembership[moduleGenes,column])
GS<-abs(geneTraitSignificance[moduleGenes, 1])
c_2<-as.data.frame(cbind(MM,GS))
rownames(c_2)=red_module$genename
head(c_2)

#              MM        GS
# A1BG  0.5773017 0.2164258
# A1CF  0.7396354 0.1606311
# AADAC 0.7820721 0.2225708
# AASS  0.6317908 0.1607554
# ABAT  0.8471267 0.2259109
# ABCA6 0.7938410 0.2123513


c_2$group <- abs(c_2$MM)>0.6&abs(c_2$GS)>0.15
write.csv(c_2, "hubgene_MMGS_red_ajcc_pathologic_stage.csv")


pdf("MM vs GS_red_ajcc_pathologic_stage.pdf",width = 7,height = 7)
ggplot(data=c_2, aes(x=MM, y=GS, color=group))+geom_point(size=1.5)+scale_colour_manual(values=c("grey60", "red"))+ theme_bw()+  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+  labs(x="Module Membership in red module", y="Gene significance for ajcc_pathologic_stage",title = "Module membership vs. gene significance ")+theme(axis.title.x =element_text(size=14), axis.title.y=element_text(size=14),axis.text = element_text(size = 12),axis.text.x = element_text(colour = "black"),axis.text.y = element_text(colour = "black"),plot.title = element_text(hjust = 0.5,size = 16,face = "bold"),plot.margin = unit(rep(2,4),'lines')) +theme(legend.position = 'none')+geom_hline(aes(yintercept=0.15),colour="#5B9BD5",lwd=1,linetype=1)+geom_vline(aes(xintercept=0.6),colour="#5B9BD5",lwd=1,linetype=1)
dev.off()





#  ajcc_pathologic_t
ajcc_pathologic_t = as.data.frame(traitData$ajcc_pathologic_t)
names(ajcc_pathologic_t) = "ajcc_pathologic_t"
modNames = substring(names(MEs), 3)
geneModuleMembership = as.data.frame(cor(dataExpr, MEs, use = "p"))
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership 
 ), nSamples))
names(geneModuleMembership) = paste("MM", modNames, sep="")
names(MMPvalue) = paste("p.MM", modNames, sep="")
# 计算性状和基因表达量之间的相关性（GS）
geneTraitSignificance = as.data.frame(cor(dataExpr, ajcc_pathologic_t, use = "p"))
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), 
  nSamples))
names(geneTraitSignificance) = paste("GS.", names(ajcc_pathologic_t), sep="")
names(GSPvalue) = paste("p.GS.", names(ajcc_pathologic_t), sep="")


#yellow 模块
module = "4" #4:yellow;6:red
column = match(module, modNames)
module='yellow'
moduleGenes = moduleColors==module
table(moduleGenes)
moduleGenes

# FALSE  TRUE 
# 19086  1015 

yellow_module<-as.data.frame(dimnames(data.frame(dataExpr))[[2]][moduleGenes]) 
names(yellow_module)="genename"

#筛选hub基因
MM<-abs(geneModuleMembership[moduleGenes,column])
GS<-abs(geneTraitSignificance[moduleGenes, 1])
c_3<-as.data.frame(cbind(MM,GS))
rownames(c_3)=yellow_module$genename
head(c_3)

#               MM         GS
# ABCA9  0.5207030 0.22219210
# ABHD3  0.3930192 0.03391066
# ABLIM3 0.3976121 0.03957260
# ACAD8  0.4048236 0.14484596
# ACBD7  0.4103006 0.10203411
# ACER1  0.3884292 0.09926795

c_3$group <- abs(c_3$MM)>0.6&abs(c_3$GS)>0.15
write.csv(c_3, "hubgene_MMGS_yellow_ajcc_pathologic_t.csv")


pdf("MM vs GS_yellow_ajcc_pathologic_t.pdf",width = 7,height = 7)
ggplot(data=c_3, aes(x=MM, y=GS, color=group))+geom_point(size=1.5)+scale_colour_manual(values=c("grey60", "red"))+ theme_bw()+  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+  labs(x="Module Membership in yellow module", y="Gene significance for ajcc_pathologic_t",title = "Module membership vs. gene significance ")+theme(axis.title.x =element_text(size=14), axis.title.y=element_text(size=14),axis.text = element_text(size = 12),axis.text.x = element_text(colour = "black"),axis.text.y = element_text(colour = "black"),plot.title = element_text(hjust = 0.5,size = 16,face = "bold"),plot.margin = unit(rep(2,4),'lines')) +theme(legend.position = 'none')+geom_hline(aes(yintercept=0.15),colour="#5B9BD5",lwd=1,linetype=1)+geom_vline(aes(xintercept=0.6),colour="#5B9BD5",lwd=1,linetype=1)
dev.off()


#red 模块
module = "6" #4:yellow;6:red
column = match(module, modNames)
module='red'
moduleGenes = moduleColors==module
table(moduleGenes)

# moduleGenes
# FALSE  TRUE 
# 19292   809 

red_module<-as.data.frame(dimnames(data.frame(dataExpr))[[2]][moduleGenes]) 
names(red_module)="genename"

#筛选hub基因
MM<-abs(geneModuleMembership[moduleGenes,column])
GS<-abs(geneTraitSignificance[moduleGenes, 1])
c_4<-as.data.frame(cbind(MM,GS))
rownames(c_4)=red_module$genename
head(c_4)

#   MM        GS
# A1BG  0.5773017 0.2085144
# A1CF  0.7396354 0.1461486
# AADAC 0.7820721 0.1972142
# AASS  0.6317908 0.1744571
# ABAT  0.8471267 0.2097060
# ABCA6 0.7938410 0.1650807

c_4$group <- abs(c_4$MM)>0.6&abs(c_4$GS)>0.15
write.csv(c_4, "hubgene_MMGS_red_ajcc_pathologic_t.csv")


pdf("MM vs GS_red_ajcc_pathologic_t.pdf",width = 7,height = 7)
ggplot(data=c_4, aes(x=MM, y=GS, color=group))+geom_point(size=1.5)+scale_colour_manual(values=c("grey60", "red"))+ theme_bw()+  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+  labs(x="Module Membership in red module", y="Gene significance for ajcc_pathologic_t",title = "Module membership vs. gene significance ")+theme(axis.title.x =element_text(size=14), axis.title.y=element_text(size=14),axis.text = element_text(size = 12),axis.text.x = element_text(colour = "black"),axis.text.y = element_text(colour = "black"),plot.title = element_text(hjust = 0.5,size = 16,face = "bold"),plot.margin = unit(rep(2,4),'lines')) +theme(legend.position = 'none')+geom_hline(aes(yintercept=0.15),colour="#5B9BD5",lwd=1,linetype=1)+geom_vline(aes(xintercept=0.6),colour="#5B9BD5",lwd=1,linetype=1)
dev.off()




#tumor_grade
tumor_grade = as.data.frame(traitData$tumor_grade)
names(tumor_grade) = "tumor_grade"
modNames = substring(names(MEs), 3)
geneModuleMembership = as.data.frame(cor(dataExpr, MEs, use = "p"))
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership 
 ), nSamples))
names(geneModuleMembership) = paste("MM", modNames, sep="")
names(MMPvalue) = paste("p.MM", modNames, sep="")
# 计算性状和基因表达量之间的相关性（GS）
geneTraitSignificance = as.data.frame(cor(dataExpr, tumor_grade, use = "p"))
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), 
  nSamples))
names(geneTraitSignificance) = paste("GS.", names(tumor_grade), sep="")
names(GSPvalue) = paste("p.GS.", names(tumor_grade), sep="")


#yellow 模块
module = "4" #4:yellow;6:red
column = match(module, modNames)
module='yellow'
moduleGenes = moduleColors==module
table(moduleGenes)

# moduleGenes
# FALSE  TRUE 
# 19086  1015 

yellow_module<-as.data.frame(dimnames(data.frame(dataExpr))[[2]][moduleGenes]) 
names(yellow_module)="genename"

#筛选hub基因
MM<-abs(geneModuleMembership[moduleGenes,column])
GS<-abs(geneTraitSignificance[moduleGenes, 1])
c_5<-as.data.frame(cbind(MM,GS))
rownames(c_5)=yellow_module$genename
head(c_5)

#               MM         GS
# ABCA9  0.5207030 0.20226059
# ABHD3  0.3930192 0.12285449
# ABLIM3 0.3976121 0.08853754
# ACAD8  0.4048236 0.19589894
# ACBD7  0.4103006 0.09273654
# ACER1  0.3884292 0.15592343

c_5$group <- abs(c_5$MM)>0.6&abs(c_5$GS)>0.15
write.csv(c_5, "hubgene_MMGS_yellow_tumor_grade.csv")


pdf("MM vs GS_yellow_tumor_grade.pdf",width = 7,height = 7)
ggplot(data=c_5, aes(x=MM, y=GS, color=group))+geom_point(size=1.5)+scale_colour_manual(values=c("grey60", "red"))+ theme_bw()+  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+  labs(x="Module Membership in yellow module", y="Gene significance for tumor_grade",title = "Module membership vs. gene significance ")+theme(axis.title.x =element_text(size=14), axis.title.y=element_text(size=14),axis.text = element_text(size = 12),axis.text.x = element_text(colour = "black"),axis.text.y = element_text(colour = "black"),plot.title = element_text(hjust = 0.5,size = 16,face = "bold"),plot.margin = unit(rep(2,4),'lines')) +theme(legend.position = 'none')+geom_hline(aes(yintercept=0.15),colour="#5B9BD5",lwd=1,linetype=1)+geom_vline(aes(xintercept=0.6),colour="#5B9BD5",lwd=1,linetype=1)
dev.off()


#red 模块
module = "6" #4:yellow;6:red
column = match(module, modNames)
module='red'
moduleGenes = moduleColors==module
table(moduleGenes)

# moduleGenes
# FALSE  TRUE 
# 19292   809 

red_module<-as.data.frame(dimnames(data.frame(dataExpr))[[2]][moduleGenes]) 
names(red_module)="genename"

#筛选hub基因
MM<-abs(geneModuleMembership[moduleGenes,column])
GS<-abs(geneTraitSignificance[moduleGenes, 1])
c_6<-as.data.frame(cbind(MM,GS))
rownames(c_6)=red_module$genename
head(c_6)

#              MM         GS
# A1BG  0.5773017 0.09996296
# A1CF  0.7396354 0.11447216
# AADAC 0.7820721 0.17240161
# AASS  0.6317908 0.22986579
# ABAT  0.8471267 0.25313439
# ABCA6 0.7938410 0.20154306

c_6$group <- abs(c_6$MM)>0.6&abs(c_6$GS)>0.15
write.csv(c_6, "hubgene_MMGS_red_tumor_grade.csv")


pdf("MM vs GS_red_tumor_grade.pdf",width = 7,height = 7)
ggplot(data=c_6, aes(x=MM, y=GS, color=group))+geom_point(size=1.5)+scale_colour_manual(values=c("grey60", "red"))+ theme_bw()+  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+  labs(x="Module Membership in red module", y="Gene significance for tumor_grade",title = "Module membership vs. gene significance ")+theme(axis.title.x =element_text(size=14), axis.title.y=element_text(size=14),axis.text = element_text(size = 12),axis.text.x = element_text(colour = "black"),axis.text.y = element_text(colour = "black"),plot.title = element_text(hjust = 0.5,size = 16,face = "bold"),plot.margin = unit(rep(2,4),'lines')) +theme(legend.position = 'none')+geom_hline(aes(yintercept=0.15),colour="#5B9BD5",lwd=1,linetype=1)+geom_vline(aes(xintercept=0.6),colour="#5B9BD5",lwd=1,linetype=1)
dev.off()
