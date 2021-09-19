# PlantNGSTools - Plant RNA-seq analysis tool



## Tutorial

[https://blog.ugeneyun.cn/software/PlantNGSTools.manual.html](https://blog.ugeneyun.cn/software/PlantNGSTools.manual.html)



## Installation



```
install.packages('devtools')
devtools::install_github('biomarble/PlantNGSTools')
```



## Module


|Usage|Module| Method                                                 |
|-|-|-|
|DEG analysis ***with*** replicates<br>有生物学重复的差异表达基因分析|DEGAnalysis_DESeq2|[DESeq2](https://doi.org/10.1186/s13059-014-0550-8)|
|DEG analysis ***without*** replicates<br>没有生物学重复的差异表达基因分析|DEGAnalysis_EBSeq|[EBSeq](https://doi.org/10.1093/bioinformatics/btt087)|
|Gene Ontology (GO) enrichment analysis<br>GO 功能富集|GOEnrich|[topGO](https://doi.org/10.1093/bioinformatics/btl140)|
|Gene Ontology (GO) enrichment analysis with eggNOG-Mapper annotation<br>基于[eggNOG-Mapper](http://eggnog-mapper.embl.de/)注释的GO富集|GOEnrich_eggnog| [topGO](https://doi.org/10.1093/bioinformatics/btl140)|
|Gene Ontology (GO) enrichment analysis with PANNZER2 annotation<br>基于[PANNZER2](http://ekhidna2.biocenter.helsinki.fi/sanspanz/)注释的GO富集|GOEnrich_pannzer2| [topGO](https://doi.org/10.1093/bioinformatics/btl140)|
|Gene Ontology (GO) enrichment analysis with custom annotation<br/>自主注释的GO富集|GOEnrich_customMapping<br>GOEnrich_customTable|[topGO](https://doi.org/10.1093/bioinformatics/btl140) |
|KEGG Pathway enrichment analysis <br>KEGG Pathway富集|KEGGenrich|-|
|KEGG Pathway enrichment analysis with Blastkoala annotation <br>基于[BlastKOALA](https://www.kegg.jp/blastkoala/)的KEGG Pathway富集|KEGGenrich_blastkoala| - |
|差异表达分析火山图<br>DEG Volcano Plot|VolcanoPlot|[ggplot2](https://ggplot2.tidyverse.org/)|
|差异表达分析MA图<br>DEG MA Plot|MAPlot|[ggplot2](https://ggplot2.tidyverse.org/)|
|KEGG富集气泡图<br>Bubble plot for KEGG enrichment|KEGGbubble|[ggplot2](https://ggplot2.tidyverse.org/)|
|GO富集气泡图<br>Bubble plot for GO enrichment|GObubble|[ggplot2](https://ggplot2.tidyverse.org/)<br>|




## Tips

- 任何疑问、提交Bug，请在[Issues](https://github.com/biomarble/PlantNGSTools/issues)反馈

- 扫码关注公众号，不定期发布培训课程：<br>

![qrcode.png](./qrcode.png)
