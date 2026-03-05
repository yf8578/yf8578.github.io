---
title: "【分析】850k甲基化数据分析-meffil"
date: 2026-03-05
permalink: /posts/2026/03/分析850k甲基化数据分析-meffil/
excerpt: 本文详细介绍了使用 meffil R包进行 850K 甲基化芯片（EWAS）数据的上游处理流程，包括原始 IDAT 文件的读取、细胞比例评估、质量控制（QC）与不合格样本过滤、探针标准化以及最终的表观基因组关联分析模型构建。
tags:
  - Bioinformatics
  - Methylation
  - EWAS
  - R
---

# https://mp.weixin.qq.com/s/Vw0o0Im8Z9RE0bXsrnByFw

![cover_image](/images/blog/b736879e3c846fdc78e31b5ff68a171a.jpg)

# 【分析】850k甲基化数据分析-meffil

 

哈喽👋，大家好，我是摸鱼凡。

上期简单介绍了一下EWAS的相关知识，那么我们应该如何分析呢？有没有一个标准的分析流程？这么多工具又要如何选择呢？

这些曾是我纠结的问题，因为能用的工具实在是太多了（ChAMP、minfi..）也没个标准，我试试这个试试那个，总是觉得做的不行。

后来看了一篇文章做的MWAS，里面用的是meffil做的，能从原始数据开始处理，生成质检质控报告等，对于计算的MWAS结果再使用bacon进行矫正，后面又使用methylGSA去做GO和KEGG富集，有基因组数据还做了meQTL......

自己也在公众号搜索过，大多数都是ChAMP还有minfi，基本上没有meffil的分析流程。其实这个meffil用起来也挺方便的，内置数据集能够估算免疫细胞比例，检测、矫正混杂因素，做EWAS分析。

**作者也贴心写了很详细的教程，https://github.com/perishky/meffil/wiki**

# meffil安装

```bash
install.packages("BiocManager")  
BiocManager::install("SmartSVA")  
BiocManager::install("illuminaio")  
BiocManager::install("limma")  
BiocManager::install("DNAcopy")  
BiocManager::install("preprocessCore")  
  
install.packages("markdown")  
install.packages("knitr")  
install.packages("remotes")  
library(remotes)  
install_github("perishky/meffil")
```

在后面分析的时候会出现类似下面的报错（参考：https://github.com/bioconda/bioconda-recipes/issues/31420）

`"ERROR; return code from pthread_create() is 22" for normalize.quantiles in bioconductor-preprocesscore`

对于这个问题下面两种方法重新安装一下 `preprocessCore` 就好了

```bash
BiocManager::install("preprocessCore", configure.args=c(preprocessCore = "--disable-threading"))  
# or  
R CMD INSTALL --configure-args="--disable-threading" preprocessCore_1.39.1.tar.gz
```

# meffil分析流程

这个图也是对meffil功能的简单介绍了，还是比较全面的。

## 原始数据读入

```bash
# ==============================================================================  
# 0. 环境设置与路径  
# ==============================================================================  
library(meffil)  
library(dplyr)  
library(data.table)  
  
# 设置多核并行 (建议设为 10-20)  
# options(mc.cores = 10)  
  
# 路径定义  
idat_dir   <- "/data/work/data"  
pheno_file <- "/data/work/00_data/pheno/138_meth_sample_group_lookup_table__with_immune_new_0731_v1.csv"  
output_dir <- "/data/work/01_analysis_meffil/01_QC"  
  
if(!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)  
# ......  
# ==============================================================================  
# 构建 IDAT 索引 (Samplesheet)  
# ==============================================================================  
message("🔍 正在扫描 IDAT 文件...")  
  
# meffil 会自动去子文件夹找文件，并生成包含 Sample_Name, Sentrix_ID, Basename 的表格  
samplesheet<- meffil.create.samplesheet(path = idat_dir, recursive = TRUE)  
  
message(paste0("✅ 扫描到 IDAT 文件组数: ", nrow(samplesheet_meffil)))  
# samplesheet_meffil 的 Sample_Name 列通常是 "2xxxxx_Rxxxxx" 这种格式
```

到这里是创建一个样本表为了后面QC做准备。

```bash
# 细胞类型参考集  
meffil.list.cell.type.references()  
# [1] "andrews and bakulski cord blood" "blood gse167998"                  
# [3] "blood gse35069"                  "blood gse35069 chen"              
# [5] "blood gse35069 complete"         "blood idoloptimized"              
# [7] "blood idoloptimized epic"        "combined cord blood"              
# [9] "cord blood gse68456"             "gervin and lyle cord blood"       
# [11] "guintivano dlpfc"                "saliva gse147318"                 
# [13] "saliva gse48472"         
# 根据自身需求进行选择  
  
# 读取原始数据开始质控  
## 这一步时间会比较长，耐心等等就好  
qc.objects <- meffil.qc(samplesheet, cell.type.reference="blood gse35069 complete", verbose=TRUE)  
# 保存qc对象  
save(qc.objects,file="qc.objects.Robj")
```

## QC设置

处理完原始数据之后就需要做一些质控了

```bash
# 定义过滤参数 (这些是常规严格阈值，如果剔除太多可适当放宽)  
qc.parameters <- meffil.qc.parameters(  
    beadnum.samples.threshold             = 0.1,  # 珠子数不够的探针比例 > 10% 则剔除样本  
    detectionp.samples.threshold          = 0.1,  # P值不达标探针比例 > 10% 则剔除样本  
    detectionp.cpgs.threshold             = 0.1,  # (探针过滤) P值不达标样本 > 10% 则剔除探针  
    beadnum.cpgs.threshold                = 0.1,  # (探针过滤) 珠子数不够样本 > 10% 则剔除探针  
    sex.outlier.sd                        = 5,    # 性别预测偏离度  
    snp.concordance.threshold             = 0.95, # SNP 一致性  
    sample.genotype.concordance.threshold = 0.9  
)  
  
# 开始质控  
qc.summary <- meffil.qc.summary(  
        qc.objects,  
        parameters = qc.parameters,  
        genotypes=genotypes  
)  
  
save(qc.summary, file="qcsummary.Robj")  
# 生成QC报告  
meffil.qc.report(qc.summary, output.file="qc-report.html")
```

## 移除不合格样本

这也是我比较喜欢的一个功能，根据刚才我们设置的qc参数，会对样本情况进行统计，看哪些样本不合格，并且将其剔除。

```bash
outlier <- qc.summary$bad.samples  
table(outlier$issue)  
index <- outlier$issue %in% c("Control probe (dye.bias)",   
                              "Methylated vs Unmethylated",  
                              "X-Y ratio outlier",  
                              "Low bead numbers",  
                              "Detection p-value",  
                              "Sex mismatch",  
                              "Genotype mismatch",  
                              "Control probe (bisulfite1)",  
                              "Control probe (bisulfite2)")  
  
outlier <- outlier[index,]  
  
length(qc.objects)  
# 剔除不合格样本  
qc.objects <- meffil.remove.samples(qc.objects, outlier$sample.name)  
length(qc.objects)  
save(qc.objects,file="qc.objects.clean.Robj")  
# 剔除样本之后还可以再生成新的质控报告检查一下
```

在我的数据中显示有几个样本预测性别和收集的表型中性别不符，可以输出不和个样本名进行检查。

下图为qc报告中检测性别的图，可以看到有三个红点，就是我那三个样本性别有误的样本。

## 标准化

接下来，我们需要估计使用多少个主成分（PCs）来调整甲基化水平以消除技术影响。我们可以采用10折交叉验证来估计在拟合n个主成分后的残差方差。随着主成分数量的增加，残差应该会持续减小。生成的图表类似于一条肘形曲线。其核心思路是选择残差突然下降时对应的主成分数量。有时可能会出现多个“肘点”，这时选择主成分数量最多的那个。

```bash
y <- meffil.plot.pc.fit(qc.objects)  
ggsave(y$plot,filename="pc.fit.pdf",height=6,width=6)  
  
# 根据绘制的图来选择PC  
pc <- 10
```

现在进行功能标准化。我们在qc.summary中识别出的由于检测分数不佳而导致的劣质CpG也会被移除。

```bash
# 标准化  
norm.objects <- meffil.normalize.quantiles(qc.objects, number.pcs=pc)  
save(norm.objects,file=paste("norm.obj.pc",pc,".Robj",sep=""))  
  
norm.beta <- meffil.normalize.samples(norm.objects, cpglist.remove=qc.summary$bad.cpgs$name)  
save(norm.beta,file=paste("norm.beta.pc",pc,".Robj",sep=""))
```

对于标准化结果，同样也能生成一个报告

在标准化报告中，将对控制矩阵和变异性最大的探针进行主成分分析。此外，还将计算前10个主成分（PC）与样本表中指定的批次变量之间的关联。这些关联检验可用于识别可能的离群值。例如，如果“载玻片（Slide）”是你的批次变量之一，它会为每个载玻片提供一个p值，而不是一个整体的p值。质量不佳的载玻片可以在标准化后被识别并移除。在norm.parameters变量中，你可以设置批次变量（如载玻片、培养板、组织等）、从控制矩阵中提取的主成分数量、从标准化β值中提取的主成分数量、可变探针的数量以及关联检验中使用的p值阈值。重要的是，要将批次变量编码为因子，以便研究主成分与批次变量（如载玻片、sentrix\_row、sentrix\_col、性别及其他批次）之间的关联。你可以通过以下方式进行检查：

```bash
str(norm.objects[[1]]$samplesheet)  
  
#You change it by running a loop  
  
for (i in 1:length(norm.objects)){  
norm.objects[[i]]$samplesheet$Slide<-as.factor(norm.objects[[i]]$samplesheet$Slide)  
norm.objects[[i]]$samplesheet$Sex<-as.factor(norm.objects[[i]]$samplesheet$Sex)  
norm.objects[[i]]$samplesheet$sentrix_row<-as.factor(norm.objects[[i]]$samplesheet$sentrix_row)  
norm.objects[[i]]$samplesheet$sentrix_col<-as.factor(norm.objects[[i]]$samplesheet$sentrix_col)  
}  
  
batch_var<-c("Slide", "plate","Sex") #PLEASE EDIT THIS LINE  
norm.parameters <- meffil.normalization.parameters(  
        norm.objects,  
        variables=batch_var,  
        control.pcs=1:10,  
        batch.pcs=1:10,  
        batch.threshold=0.01  
)  
  
pcs <- meffil.methylation.pcs(norm.beta,probe.range=20000)  
save(pcs,file="pcs.norm.beta.Robj")  
norm.summary <- meffil.normalization.summary(norm.objects, pcs=pcs,parameters=norm.parameters)  
meffil.normalization.report(norm.summary, output.file="normalization-report.html")
```

提取细胞比例

```bash
cc<-t(sapply(qc.objects, function(obj) obj$cell.counts$counts))  
cc<-data.frame(IID=row.names(cc),cc)  
write.table(cc,"cellcounts.txt",sep="\t",row.names=F,col.names=T,quote=F)
```

## 过滤XY染色体探针

```bash
> meffil.list.featuresets()  
[1] "450k"            "common"          "epic"            "epic2"            
[5] "mouse"           "450k:epic"       "450k:epic2"      "epic:epic2"       
[9] "450k:epic:epic2"  
# 选择对应平台  
featureset<-"epic"  
autosomal.sites <- meffil.get.autosomal.sites(featureset)  
autosomal.sites <- intersect(autosomal.sites, rownames(norm.beta))  
norm.beta <- norm.beta[autosomal.sites,]
```

## EWAS

Meffil 实现了四种不同的表观基因组关联研究（EWAS）模型：

1. 无=无协变量

2. 全部=用户协变量

3. sva=除用户协变量外的替代变量

4. isva=除用户协变量外的独立替代变量

```bash
load(beta_file) # methylation betas: CpGs in rows, samples in cols  
                # assume that it loads a matrix called 'norm.beta'  
covs <- read.table(covs_file, he=T, stringsAsFactors=FALSE) #covariate file: samples in rows, covs in cols  
                # assume that the rownames of 'covs' correspond to the column names of 'norm.beta'  
phen <- read.table(phen_file, he=T, stringsAsFactors=FALSE) #phenotype file: samples in rows, phenotypes in cols   
                # assume that the rownames of 'phen' correspond to the column names of 'norm.beta'  
                # assume that the first column of 'phen' is the phenotype of interest for the EWAS  
  
phen<-na.omit(phen)  
  
norm.beta <- norm.beta[, colnames(norm.beta) %in% rownames(phen)]  
phen <- phen[match(colnames(norm.beta), rownames(phen)), , drop=FALSE]  
stopifnot(identical(rownames(phen), colnames(norm.beta)))  
      
covs <- covs[match(colnames(norm.beta), rownames(covs)), , drop=FALSE]  
stopifnot(identical(rownames(covs), colnames(norm.beta)))  
  
set.seed(123456) #to make siva and sva reproducible you need to set a seed.      
ewas.ret <- meffil.ewas(norm.beta, variable=phen[,1], covariates=covs, isva=F)  
save(ewas.ret, file="myewas.Robj")  
  
# This function sets your default model and your significance threshold. It will show statistics for all four models for all CpGs that pass significance threshold with the default model.  
  
ewas.parameters <- meffil.ewas.parameters(sig.threshold=1e-7,  ## EWAS p-value threshold  
                                          max.plots=100, ## plot at most 100 CpG sites  
                                          qq.inflation.method="median",  ## measure inflation using median  
                                          model="sva") ## select default EWAS model;   
  
  
ewas.summary<-meffil.ewas.summary(ewas.ret,norm.beta,parameters=ewas.parameters)                                
  
meffil.ewas.report(ewas.summary, output.file=paste(report_file,".ewas.report.html",sep=""))
```

 

