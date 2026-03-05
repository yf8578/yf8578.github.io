---
title: "RNA-seq上游：从FASTQ到表达矩阵"
date: 2026-03-05
permalink: /posts/2026/03/RNA-seq-上游---从FASTQ到表达矩阵/
excerpt: RNA-seq 基础数据处理全流程指南。从测序下机获取的 FASTQ 文件开始，详细演示了使用 fastp 进行质控与接头切除、使用 HISAT2 过滤 rRNA 并进行基因组比对、最终通过 featureCounts 生成基因表达定量矩阵（TPM）的标准实战代码。
tags:
  - Bioinformatics
  - RNA-seq
---

# https://mp.weixin.qq.com/s/3S14HI845LhahnNGO3p8sA

![cover_image](/images/blog/70c223f2f02c42b680c0e4f1599b02f6.jpg)

# RNA-seq 上游 | 从FASTQ到表达矩阵

 

我是摸鱼凡，摸鱼时间写了RNA-seq上游教程...

本科毕设就是一个RNA-seq上有的流程，当时也测了不少软件，跑了好多pipeline，下面先介绍一下基本概念，之后是自己做过的一个流程。

还想要深入了解的话，记得之前看过一篇综述，很是详细（见下图）。

![https://doi.org/10.1038/s41576-019-0150-2]( "null")

再说回来，RNA-seq的上游主要目的是从测序得到的FASTQ经过一系列处理后得到基因的定量文件（看研究目的哈，也可以输出可变剪接，转录本等定量结果，不用怕，处理流程大差不差的，基本上都是要经过质控，比对，定量这几个步骤。）。

# RNA-seq介绍

RNA-seq（RNA sequencing）是一种基于二代测序技术（Next-Generation Sequencing, NGS）研究转录组学的方法，能够快速获取给定时刻的某一基因组中RNA的种类和数量，从而实现定性与定量分析。

随着二代测序技术的发展，RNA-seq的应用场景也不断拓展，在生物学和临床医学等领域都有广泛应用。如今，RNA-seq在诊断、预后和治疗各种疾病（包括癌症和传染病）的应用中具有广阔前景。

RNA-seq数据分析的基础是下机数据的质量检查与控制（简称“质检、质控”），即对测序质量总体上有所了解，并进行质量控制，保留满足后续分析的高质量数据。不同的生物信息学工具(如FastQC)能够对测序数据整体及评分，如Phred质量评分、读长长度分布、GC含量、接头含量、重复读长等。在需要去除接头时，可以使用Cutadapt，Fastp和Trimmomatic等工具。

下一步，使用splice-aware算法将原始数据reads映射到人类参考基因组，例如STAR、TopHat2或HISAT2。此时，必须根据研究类型和表型调整一些重要的变量，如基因组的参考版本、参考基因组的注释文件等。为了丰富参考基因组注释，可以通过融合多个数据库得到较为全面的背景数据集。

在将读数映射到基因组后，有一些技术和生物学偏好会影响灵敏度阈值。比对结果中3’末端基因的偏好性可能表明RNA的降解，或者表明数据来自3’测定(例如，寡聚T起始的3’RNA-seq)，可以使用RSeQC等工具来评估BAM文件中基因的覆盖率。

RNA-seq最广泛的应用就是用来评估基因和转录本的表达，这一应用主要是基于比对到转录组区间内的读长的数量。例如HTSeq，FeatureCounts和Salmon，这些工具可以量化特定基因特征中的映射读数。其中的一些偏差，如基因长度或GC含量可能会影响量化的过程，并对差异表达分析(DEA)产生负面影响。目前已经有几种方法，能够减少这些偏差。包括基于基因长度和文库大小(每次复制的读取总数)将读长计数归一化。常用的方法包括对单端读长以RPKM(每千个碱基的转录每百万映射读取的读长)为单位，双端读长以FPKM（每千个碱基的转录每百万映射读取的片段）作为单位，和以TPM（每百万映射读取的读长）作为单位。

基因表达水平定量之后的分析可使用相应的软件包来实现。差异表达分析可选择NOIseq识别RNA-seq数据中的误差来源，并使用相应的方法来标准化这些误差，同时也可以使用COMBAT和ARSyN等及进行批次矫正。此外还可使用edgeR、DESeq、baySeq和EBSeq等进行分析。可变剪切分析可以使用BASIS、CuffDiff2、rSeqDiff等进行分析…

许多独立的研究都已经证实，选择不同的方法会对结果有一定影响，而且没有哪一种方法能够适用于所有的数据，所以，在具体分析的时候需要使用多个软件进行相互验证。

# RNA-seq pipeline

关于上游的流程已经有很多工具可以选择，实际使用中可以根据自己的需求选择合适的工具组合。

![https://doi.org/10.1038/s41467-017-00050-4]( "null")

自己本科时期用的一套流程（其实现在也在用）如下：

![]( "null")

## 流程主体

### Raw Data QC

#### step1 Cut adapter

这一步主要是对FASTQ数据进行质控去除接头的，这里使用的是非常非常好用的`fastp`，编写代码时如果知道接头可以指定接头序列，不知道的话可以设置让其自动检测。

`fastp`
```bash
fastp \\  
-i $fq1 -I $fq2  \\  
-o $outdir_step1_fastp/N1.r1.fq -O $outdir_step1_fastp/N1.r2.fq  \\  
-h $outdir_step1_fastp/cutadapter_fastp.html -j $outdir_step1_fastp/cutadapter_fastp.json \\  
--adapter_sequence="接头序列"\\  
--adapter_sequence_r2="接头序列"\\  
--disable_trim_poly_g \\  
--disable_quality_filtering  \\  
--disable_length_filtering \\  
--dont_eval_duplication \\  
--thread=12  -c \\  
--failed_out $outdir_step1_fastp/${id}_cutadapter.failed_out.fq.gz &&
```

#### step2 Quality control

这里的质控真的是字面意义上的质量控制，是对FASTQ数据中测得的比较低质量的base/reads进行处理。

```bash
fastp \\  
-i $outdir_step1_fastp/N1.r1.fq -I $outdir_step1_fastp/N1.r2.fq  \\  
-o $outdir_step1_fastp/N3.r1.fq -O $outdir_step1_fastp/N3.r2.fq  \\  
-h $outdir_step1_fastp/quality_filter.fastp.html -j $outdir_step1_fastp/quality_filter_.fastp.json \\  
--disable_adapter_trimming \\  
--thread=12 \\  
--disable_trim_poly_g \\  
--qualified_quality_phred=30 \\ ##设置过滤质量值  
--unqualified_percent_limit=50 \\  
--disable_length_filtering \\  
--n_base_limit=10 \\  
-p  -c \\  
--failed_out $outdir_step1_fastp/${id}_quality_failed_out.fq.gz &&
```

#### step3 Length control

这里是对长度进行控制，因为当时做的是游离核酸，本身长度就比较短，所以这里需要再特地做一下长度控制。

```bash
/fastp \\  
-i $outdir_step1_fastp/N3.r1.fq -I $outdir_step1_fastp/N3.r2.fq  \\  
-o $outdir_step1_fastp/${id}.r1.fq.gz -O $outdir_step1_fastp/${id}.r2.fq.gz  \\  
-h $outdir_step1_fastp/${id}_lengthfilter.fastp.html -j $outdir_step1_fastp/${id}_lengthfilter.fastp.json \\  
--disable_adapter_trimming \\  
--thread=12 \\  
--disable_trim_poly_g \\  
--length_required=指定 \\  #设置过滤长度  
--disable_quality_filtering  \\  
--dont_eval_duplication \\  
-p  -c \\  
--failed_out $outdir_step1_fastp/${id}_length_failed_out.fq.gz &&  
  
  
rm $outdir_step1_fastp/N3.r1.fq  
rm $outdir_step1_fastp/N3.r2.fq
```

**上述几段代码可以合在一起做质控，之前拆开写是需要检查每个部分的结果，流程跑通就可以再优化啦。**

### rRNA deletion

看情况吧，可能大多数是不需要在这里去rRNA的。之前处理的数据会有很多rRNA，进而导致测序数据比对率非常低，因此这里是在生信层面进行去除的一步操作。如果需要的话还要自己去制作一个rRNA的索引哈。

```bash
hisat2 -k 10 -p 12 --no-unal \\  
-x  rRNA_indedx -S $outdir_step2_rRNA/${id}_accepted_hits.sam \\  
--fr --dta --avoid-pseudogene \\  
--un-conc-gz $outdir_step2_rRNA/${id}_non_rRNA --summary-file $outdir_step2_rRNA/${id}_rRNA_result.txt \\  
-1  $outdir_step1_fastp/$id.r1.fq.gz -2 $outdir_step1_fastp/$id.r2.fq.gz &&  
  
rm $outdir_step2_rRNA/${id}_accepted_hits.sam &&  
mv $outdir_step2_rRNA/${id}_non_rRNA.1 $outdir_step2_rRNA/${id}_non_rRNA.1.fq.gz &&  
mv $outdir_step2_rRNA/${id}_non_rRNA.2 $outdir_step2_rRNA/${id}_non_rRNA.2.fq.gz &&
```

### mapping

太牛了，已经做到比对了（手动点赞）！！

这里比对工具有很多可以选择`STAR`、`HISAT2`、`bowtie`....这些工具都是比较经典的，根据自己需求选择就好。

`STAR`
`HISAT2`
`bowtie`
```bash
hisat2 -k 10 -p 12 --no-unal \\  
-x  index_path  -S $outdir_step3_hisat2/${id}_accepted_hits.sam \\  
--novel-splicesite-outfile $outdir_step3_hisat2/${id}_junctions.bed  \\  
--fr --dta --avoid-pseudogene \\  
--un-conc-gz $outdir_step3_hisat2/${id}_unmapped --summary-file $outdir_step3_hisat2/${id}_hisat2_result.txt \\  
-1  $outdir_step2_rRNA/${id}_non_rRNA.1.fq.gz -2 $outdir_step2_rRNA/${id}_non_rRNA.2.fq.gz &&  
samtools view -@ 4 -Su $outdir_step3_hisat2/${id}_accepted_hits.sam | samtools sort -@ 4 -o $outdir_step3_hisat2/${id}_accepted_hits.sorted.sam -O SAM &&   
samtools view -@ 4 -bS $outdir_step3_hisat2/${id}_accepted_hits.sorted.sam > $outdir_step3_hisat2/${id}_accepted_hits.sorted.bam && samtools index $outdir_step3_hisat2/${id}_accepted_hits.sorted.bam &&   
  
rm $outdir_step3_hisat2/${id}_accepted_hits.sam $outdir_step3_hisat2/${id}_accepted_hits.sorted.sam
```

### Quantification

比对之后我们得到的是二进制的`bam`文件，这一步看自己需要的是什么矩阵，再选择相应的工具。

`bam`

此处使用featureCounts对上一步的比对结果进行定量，速度还是非常快的。

```bash
featureCounts \  
-T 12 -t gene  -g gene_id -B -C  -p   -O --minOverlap 10 \  
-a  GCF_000001405.40_GRCh38.p14_genomic.gtf \  
-o ${outdir_step4_featureCounts_1}/$id.featurecounts.txt ${dir}/result/$id/step3_hisat2/${id}_accepted_hits.sorted.bam
```

![]( "null")

之后便可得到各基因定量的结果，**但要注意的是，这里输出的计数是reads数目（整数），需要根据自己的需要将各基因的表达量进行标准化**。

# 结束

一个RNA-seq上游的分析就结束啦，每个样本都会有一个定量文件，通过代码将其整合起来就好了（见下图），最终样式是第一列为基因ID，之后为每个样本及其对应的表达量。

![]( "null")

你可能已经发现了，怎么这里面基因的表达量跟刚刚说的整数不一样呢？？！！

**回答：这里的表达量是经过标准化之后的TPM，可以不是整数喔。下游分析的时候会有软件要求输入必须是整数，不过这里的分析软件也是可以根据自己需要选择的，没那么古板。**  
上述代码仅供参考，具体参数需要根据实际进行调整~~~

有任何疑问和建议，欢迎交流探讨~~~

我是摸鱼凡，摸鱼时间写生信教程，希望能为大家求学路上撑一把伞

# 参考

1. Sahraeian S M E, Mohiyuddin M, Sebra R, et al. Gaining comprehensive biological insight into the transcriptome by performing a broad-spectrum RNA-seq analysis[J]. Nature communications, 2017, 8(1): 1-15.

2. Conesa A, Madrigal P, Tarazona S, et al. A survey of best practices for RNA-seq data analysis[J]. Genome biology, 2016, 17(1): 1-19.

3. Marco-Puche G, Lois S, Benítez J, et al. RNA-Seq perspectives to improve clinical diagnosis[J]. Frontiers in genetics, 2019, 10: 1152.

 

