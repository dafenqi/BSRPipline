# BSRPipline
 Analysis Pipline Of Bulked Segregant RNA-Seq
BSR模块使用说明

一 关联分析
1 SNP-index是近年来发表的一种通过混池间的基因型频率差异进行标记关联分析的方法，主要是寻找混池之间基因型频率的显著差异，用Δ(SNP-index)统计。标记SNP与性状关联度越强，Δ(SNP-index)越接近于1。
2 G value方法，Magwene等人提出了另一种确定来自NGS-BSA的QTL统计学显着性的方法。基于观察到的和预期的等位基因深度计算每个SNP的修正的G统计量并使用三倍体平滑核平滑该值。使用平滑的G统计量或G'。该方法允许降低噪声，同时还解决SNP之间的连锁不平衡。此外，由于G'接近于正常分布的对数，因此可以使用G'的零分布的非参数估计来估计每个SNP的p值。这提供了清晰且易于理解的结果以及多个测试校正的选项。

二 流程使用方法

2.1 功能描述
检测SNP，筛选候选区间
2.2 流程准备工作
1 配置 run_BSR.sh 执行脚本
-R 参考基因组文件
-t bam/sam 顺序按突变型，野生型
-i bam文件夹 顺序按突变型，野生型
-o 输出文件夹SNP
-b bam文件名 去掉.bam
-n 样品名
-gff gtf文件路径
--scriptsdir 对应脚本路径
--toolsdir 对应工具路径
