# PlantNGSTools - Plant RNA-seq analysis tool


**Latest version:** [![Latest Version](https://img.shields.io/github/release/biomarble/PlantNGSTools.svg?style=flat?maxAge=86400)](https://github.com/biomarble/PlantNGSTools/releases)


## Tutorial

[https://blog.ugeneyun.cn/software/PlantNGSTools.manual.html](https://blog.ugeneyun.cn/software/PlantNGSTools.manual.html)



## Installation



```
install.packages('devtools')
devtools::install_github('biomarble/PlantNGSTools')
```



## Module


|用途|模块名| 备注                                                 |
|-|-|-|
|DEG analysis ***with*** replicates<br>有生物学重复的差异表达基因分析|DEGAnalysis_DESeq2|[DESeq2](https://doi.org/10.1186/s13059-014-0550-8)|
|DEG analysis ***without*** replicates<br>没有生物学重复的差异表达基因分析|DEGAnalysis_EBSeq|[EBSeq](https://doi.org/10.1093/bioinformatics/btt087)|
|Gene Ontology (GO) enrichment analysis<br>GO 功能富集|GOEnrich|[topGO](https://doi.org/10.1093/bioinformatics/btl140)|
|Gene Ontology (GO) enrichment analysis with eggNOG-Mapper annotation<br>基于[eggNOG-Mapper](http://eggnog-mapper.embl.de/)注释的GO富集|GOEnrich_eggnog| [topGO](https://doi.org/10.1093/bioinformatics/btl140)|
|Gene Ontology (GO) enrichment analysis with PANNZER2 annotation<br>基于[PANNZER2](http://ekhidna2.biocenter.helsinki.fi/sanspanz/)注释的GO富集|GOEnrich_pannzer2| [topGO](https://doi.org/10.1093/bioinformatics/btl140)|
|Gene Ontology (GO) enrichment analysis with custom annotation<br/>自主注释的GO富集|GOEnrich_customMapping<br>GOEnrich_customTable|[topGO](https://doi.org/10.1093/bioinformatics/btl140) |
|KEGG Pathway enrichment analysis <br>KEGG Pathway富集|KEGGenrich|-|
|KEGG Pathway enrichment analysis with Blastkoala annotation <br>基于[BlastKOALA](https://www.kegg.jp/blastkoala/)的KEGG Pathway富集|KEGGenrich_blastkoala| - |
|KEGG Pathway enrichment analysis with custom pathway annotation <br>基于自主注释的KEGG Pathway富集|KEGGenrich_customTable| - |
|DEG Volcano Plot<br>差异表达分析火山图|VolcanoPlot|[ggplot2](https://ggplot2.tidyverse.org/)|
|DEG MA Plot<br>差异表达分析MA图|MAPlot|[ggplot2](https://ggplot2.tidyverse.org/)|
|Bubble plot for KEGG enrichment<br>KEGG富集气泡图|KEGGbubble|[ggplot2](https://ggplot2.tidyverse.org/)|
|Bubble plot for GO enrichment<br>GO富集气泡图|GObubble|[ggplot2](https://ggplot2.tidyverse.org/)<br>|
|Bar plot for secondary GO Terms<br>GO二级节点[^1]条形图|GOBar|[ggplot2](https://ggplot2.tidyverse.org/)<br>|
|Trait statistics<br>表型常用统计量[^2]|multiTraitStat|-|
|Trait plot<br>表型统计图[^3]| multiTraitPlot|-|
|Expression matrix combining<br>表达量矩阵按组合并|matrixGroup|-|

[^1]: GO的二级节点是MF、BP、CC三个主节点的直接子节点
[^2]: 最大/小值，均值，中位数，偏度，峰度，shannon多态性指数，变异系数
[^3]: 直方图，相关性散点图，相关性显著性,相关系数

## Tips

- 任何疑问、提交Bug，请在[Issues](https://github.com/biomarble/PlantNGSTools/issues)反馈

- 扫码关注公众号，不定期发布培训课程：<br>

![qrcode.png](./qrcode.png)
