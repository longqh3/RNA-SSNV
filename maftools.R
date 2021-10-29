# 修改R包路径
.libPaths(c("/home/lqh/R/x86_64-conda_cos6-linux-gnu-library/3.6"))
# 设置待分析癌症类型
CANCER_TYPE = "LUSC"

# 载入所需的R包
library("survival") # 生存分析相关包
library("survminer") # 生存曲线绘制相关包
library(maftools) # maftools主要分析包
library(mclust) # 癌症驱动基因分析包
library(readxl) # 读取excel数据文件
library(data.table)
library(tidyr)

# 设置工作文件夹存放结果
setwd("/home/lqh/test/clinical/LUSC")

### 批量读取所有maf文件+临床资料文件 ###
# Clinical data
# 读取临床信息
# Read in clinical file

# 引入第三方标准数据CDR数据集信息
# 添加PFI等信息进入临床信息内
# 添加DSS、DFI、PFI信息
CDR_gdc_clin = read_excel("/home/lqh/Codes/Data/TCGA_file_list/TCGA-CDR-SupplementalTableS1.xlsx",
                          sheet = "TCGA-CDR") %>% as.data.frame()
names(CDR_gdc_clin)[names(CDR_gdc_clin) == 'bcr_patient_barcode'] <- 'Tumor_Sample_Barcode' # 修改"bcr_patient_barcode"列名为规范形式
names(CDR_gdc_clin)[names(CDR_gdc_clin) == 'age_at_initial_pathologic_diagnosis'] <- 'age' # 修改"age_at_initial_pathologic_diagnosis"列名为规范形式
# 选择对应癌症项目数据(暂且以添加为主，未来再以其为主进行分析)
CDR_gdc_clin = subset(CDR_gdc_clin, type==CANCER_TYPE,
                      select = c("Tumor_Sample_Barcode", "gender", "race", "ajcc_pathologic_tumor_stage", "histological_type", "age", "treatment_outcome_first_course",
                                 "OS", "OS.time", "DSS", "DSS.time", "DFI", "DFI.time", "PFI", "PFI.time"))

# 去除结局变量相关空值
CDR_gdc_clin<-subset(CDR_gdc_clin,(OS.time!="NA")|(OS.time!="")) # OS删除空值
# 去除年龄相关空值
CDR_gdc_clin<-subset(CDR_gdc_clin,(age!="NA"))


# Convert value
# 批量修改性别、年龄信息
# male -> 1, female -> 0
CDR_gdc_clin$gender <- factor(CDR_gdc_clin$gender,    #转为factor形式，替换标签
                          levels=c("FEMALE", "MALE"),
                          labels=c(0, 1))
# age<65 -> 0, age≥65 -> 1
CDR_gdc_clin$age_index = ifelse(CDR_gdc_clin$age<65, 0, 1)

# 数据类型转换
CDR_gdc_clin$OS.time = as.numeric(as.character(CDR_gdc_clin$OS.time))
CDR_gdc_clin$OS = as.numeric(CDR_gdc_clin$OS)
CDR_gdc_clin$PFI.time = as.numeric(as.character(CDR_gdc_clin$PFI.time))
CDR_gdc_clin$PFI = as.numeric(CDR_gdc_clin$PFI)
CDR_gdc_clin$age = as.numeric(as.character(CDR_gdc_clin$age))

# 读取对应maf文件信息
RNA_maf_info = maftools::read.maf("/home/lqh/Codes/Python/Integrative_Analysis_Bioinformatics_Pipeline/results/LUSC/RNA/RNA_somatic_mutation/VcfAssembly_new/DNA_only_negative_maftools_updated.txt",
                                  clinicalData = CDR_gdc_clin,
                                  isTCGA = TRUE)
WXS_maf_info = maftools::read.maf("/home/lqh/Codes/Data/TCGA_maf_files/TCGA-LUSC",
                                  clinicalData = CDR_gdc_clin,
                                  isTCGA = TRUE)

### 突变的可视化 ###
# MAF文件汇总统计图
plotmafSummary(maf=WXS_maf_info, rmOutlier=TRUE, addStat="median", dashboard=TRUE, titvRaw = FALSE)
plotmafSummary(maf=RNA_maf_info, rmOutlier=TRUE, addStat="median", dashboard=TRUE, titvRaw = FALSE)
# Oncoplot(瀑布图)
# 通过top=20设定要展示的频率前20的突变
oncoplot(maf=RNA_maf_info, top=50, borderCol=NULL) # 有待进一步调整，详见函数帮助
oncoplot(maf=WXS_maf_info, top=20, borderCol=NULL) # 有待进一步调整，详见函数帮助

### 突变的数据分析 ###
# 生存分析
# 使用函数mafSurvival可以利用突变和临床数据绘制Kaplan-Meier曲线（KM曲线）进行生存分析。
# 参数genes用于选定基因，time设定临床数据中保存随访时间的字段，Status设定临床数据中存放生存状态的字段。

# Predict gene-sets associated with survival
# 应用生存分析探索top 20或30的突变基因对应突变状态与OS、PFI长度的关系
prog_geneset = survGroup(maf = RNA_maf_info, top = 50,
                         geneSetSize = 1,
                         time = "OS.time",
                         Status = "OS",
                         verbose = FALSE)
print(prog_geneset)
prog_geneset = survGroup(maf = RNA_maf_info, top = 30,
                         geneSetSize = 1,
                         time = "PFI.time",
                         Status = "PFI",
                         verbose = FALSE)
print(prog_geneset)
prog_geneset = survGroup(maf = WXS_maf_info, top = 50,
                         geneSetSize = 1,
                         time = "OS.time",
                         Status = "OS",
                         verbose = FALSE)
print(prog_geneset)
prog_geneset = survGroup(maf = WXS_maf_info, top = 50,
                         geneSetSize = 1,
                         time = "PFI.time",
                         Status = "PFI",
                         verbose = FALSE)
print(prog_geneset)
# 针对某给定基因进行单因素生存分析结果图绘制，并给出log-rank检验P值
# RNA as compare
mafSurvGroup(maf = RNA_maf_info, geneSet = c("TP53"),
             time = "OS.time",
             Status = "OS")

mafSurvGroup(maf = RNA_maf_info, geneSet = c("COPA"),
             time = "PFI.time",
             Status = "PFI")
fit <- survfit(Surv(PFI.time, PFI) ~ PRKDC, data=RNA_maf_info)
p1 <- ggsurvplot(fit)
p1

# DNA as base line
mafSurvGroup(maf = WXS_maf_info, geneSet = c("BIRC6"),
             time = "OS.time",
             Status = "OS")
# Test single gene
mafSurvival(maf=WXS_maf_info, genes="WHSC1L1",
            time = "OS.time",
            Status = "OS", isTCGA=TRUE)

mafSurvival(maf=WXS_maf_info, genes="PIK3CA",
            time = "OS.time",
            Status = "OS", isTCGA=TRUE)