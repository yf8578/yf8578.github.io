---
title: '【分析】850k甲基化数据分析-meffil'
date: 2026-03-05
permalink: /posts/2026/03/850k-methylation-meffil/
tags:
  - Bioinformatics
  - Methylation
  - EWAS
  - R
---

哈喽👋，大家好，我是摸鱼凡。

上期简单介绍了一下EWAS的相关知识，那么我们应该如何分析呢？有没有一个标准的分析流程？这么多工具又要如何选择呢？

这些曾是我纠结的问题，因为能用的工具实在是太多了（ChAMP、minfi..）也没个标准，我试试这个试试那个，总是觉得做的不行。

后来看了一篇文章做的MWAS，里面用的是meffil做的，能从原始数据开始处理，生成质检质控报告等，对于计算的MWAS结果再使用bacon进行矫正，后面又使用methylGSA去做GO和KEGG富集，有基因组数据还做了meQTL......

自己也在公众号搜索过，大多数都是ChAMP还有minfi，基本上没有meffil的分析流程。其实这个meffil用起来也挺方便的，内置数据集能够估算免疫细胞比例，检测、矫正混杂因素，做EWAS分析。

作者也贴心写了很详细的教程，[https://github.com/perishky/meffil/wiki](https://github.com/perishky/meffil/wiki)

## meffil安装

```R
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

`"ERROR; return code from pthread_create() is 22" for normalize.quantiles in bioconductor-preprocesscore`

对于这个问题下面两种方法重新安装一下 `preprocessCore` 就好了

```R
BiocManager::install("preprocessCore", configure.args=c(preprocessCore = "--disable-threading"))  
# or  
R CMD INSTALL --configure-args="--disable-threading" preprocessCore_1.39.1.tar.gz
```

## meffil分析流程

### 原始数据读入

```R
# ==============================================================================  
# 0. 环境设置与路径  
# ==============================================================================  
library(meffil)  
library(dplyr)  
library(data.table)  
  
# 设置多核并行 (建议设为 10-20)  
# options(mc.cores = 10)  
  
# 路径定义  
idat_dir   <- "/data/work/data"  
pheno_file <- "/data/work/00_data/pheno/138_meth_sample_group_lookup_table__with_immune_new_0731_v1.csv"  
output_dir <- "/data/work/01_analysis_meffil/01_QC"  
  
if(!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)  
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

```R
# 细胞类型参考集  
meffil.list.cell.type.references()  
# [1] "andrews and bakulski cord blood" "blood gse167998"                 
# [3] "blood gse35069"                  "blood gse35069 chen"             
# [5] "blood gse35069 complete"         "blood idoloptimized"             
# [7] "blood idoloptimized epic"        "combined cord blood"             
# [9] "cord blood gse68456"             "gervin and lyle cord blood"      
# [11] "guintivano dlpfc"                "saliva gse147318"                
# [13] "saliva gse48472"                 
# 根据自身需求进行选择  
  
# 读取原始数据开始质控  
## 这一步时间会比较长，耐心等等就好  
qc.objects <- meffil.qc(samplesheet, cell.type.reference="blood gse35069 complete", verbose=TRUE)  
# 保存qc对象  
save(qc.objects,file="qc.objects.Robj")
```

### QC设置

处理完原始数据之后就需要做一些质控了

```R
# 定义过滤参数 (这些是常规严格阈值，如果剔除太多可适当放宽)  
qc.parameters <- meffil.qc.parameters(  
    beadnum.samples.threshold             = 0.1,  # 珠子数不够的探针比例 > 10% 则剔除样本  
    detectionp.samples.threshold          = 0.1,  # P值不达标探针比例 > 10% 则剔除样本  
    detectionp.cpgs.threshold             = 0.1,  # (探针过滤) P值不达标样本 > 10% 则剔除探针  
    beadnum.cpgs.threshold                = 0.1,  # (探针过滤) 珠子数不够样本 > 10% 则剔除探针  
    sex.outlier.sd                        = 5,    # 性别预测偏离度  
    snp.concordance.threshold             = 0.95, # SNP 一致性  
    sample.genotype.concordance.threshold = 0.9  
)  
  
# 开始质控  
qc.summary <- meffil.qc.summary(  
        qc.objects,  
        parameters = qc.parameters,  
        genotypes=genotypes  
)  
  
save(qc.summary, file="qcsummary.Robj")  
# 生成QC报告  
meffil.qc.report(qc.summary, output.file="qc-report.html")
```

### 移除不合格样本

根据刚才我们设置的qc参数，会对样本情况进行统计，看哪些样本不合格，并且将其剔除。

```R
outlier <- qc.summary$bad.samples  
table(outlier$issue)  
index <- outlier$issue %in% c("Control probe (dye.bias)",   
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

### 标准化

接下来，我们需要估计使用多少个主成分（PCs）来调整甲基化水平以消除技术影响。我们可以采用10折交叉验证来估计在拟合n个主成分后的残差方差。

```R
y <- meffil.plot.pc.fit(qc.objects)  
ggsave(y$plot,filename="pc.fit.pdf",height=6,width=6)  
  
# 根据绘制的图来选择PC  
pc <- 10
```

现在进行功能标准化。我们在qc.summary中识别出的由于检测分数不佳而导致的劣质CpG也会被移除。

```R
# 标准化  
norm.objects <- meffil.normalize.quantiles(qc.objects, number.pcs=pc)  
save(norm.objects,file=paste("norm.obj.pc",pc,".Robj",sep=""))  
  
norm.beta <- meffil.normalize.samples(norm.objects, cpglist.remove=qc.summary$bad.cpgs$name)  
save(norm.beta,file=paste("norm.beta.pc",pc,".Robj",sep=""))
```

提取细胞比例：

```R
cc<-t(sapply(qc.objects, function(obj) obj$cell.counts$counts))  
cc<-data.frame(IID=row.names(cc),cc)  
write.table(cc,"cellcounts.txt",sep="\t",row.names=F,col.names=T,quote=F)
```

### 过滤XY染色体探针

```R
> meffil.list.featuresets()  
[1] "450k"            "common"          "epic"            "epic2"           
[5] "mouse"           "450k:epic"       "450k:epic2"      "epic:epic2"      
[9] "450k:epic:epic2"  
# 选择对应平台  
featureset<-"epic"  
autosomal.sites <- meffil.get.autosomal.sites(featureset)  
autosomal.sites <- intersect(autosomal.sites, rownames(norm.beta))  
norm.beta <- norm.beta[autosomal.sites,]
```

### EWAS

Meffil 实现了四种不同的表观基因组关联研究（EWAS）模型：
1. 无=无协变量
2. 全部=用户协变量
3. sva=除用户协变量外的替代变量
4. isva=除用户协变量外的独立替代变量

```R
load(beta_file) # methylation betas: CpGs in rows, samples in cols  
                # assume that it loads a matrix called 'norm.beta'  
covs <- read.table(covs_file, he=T, stringsAsFactors=FALSE) #covariate file: samples in rows, covs in cols  
                # assume that the rownames of 'covs' correspond to the column names of 'norm.beta'  
phen <- read.table(phen_file, he=T, stringsAsFactors=FALSE) #phenotype file: samples in rows, phenotypes in cols   
                # assume that the rownames of 'phen' correspond to the column names of 'norm.beta'  
                # assume that the first column of 'phen' is the phenotype of interest for the EWAS  
  
phen<-na.omit(phen)  
  
norm.beta <- norm.beta[, colnames(norm.beta) %in% rownames(phen)]  
phen <- phen[match(colnames(norm.beta), rownames(phen)), , drop=FALSE]  
stopifnot(identical(rownames(phen), colnames(norm.beta)))  
    
covs <- covs[match(colnames(norm.beta), rownames(covs)), , drop=FALSE]  
stopifnot(identical(rownames(covs), colnames(norm.beta)))  
  
set.seed(123456) #to make siva and sva reproducible you need to set a seed.     
ewas.ret <- meffil.ewas(norm.beta, variable=phen[,1], covariates=covs, isva=F)  
save(ewas.ret, file="myewas.Robj")  
  
# This function sets your default model and your significance threshold. It will show statistics for all four models for all CpGs that pass significance threshold with the default model.  
  
ewas.parameters <- meffil.ewas.parameters(sig.threshold=1e-7,  ## EWAS p-value threshold  
                                          max.plots=100, ## plot at most 100 CpG sites  
                                          qq.inflation.method="median",  ## measure inflation using median  
                                          model="sva") ## select default EWAS model;   
  
  
ewas.summary<-meffil.ewas.summary(ewas.ret,norm.beta,parameters=ewas.parameters)                               
  
meffil.ewas.report(ewas.summary, output.file=paste(report_file,".ewas.report.html",sep=""))
```
