colorScheme = function(colorsheme = "ryb", n = 100) {
    if (colorsheme == "ryb")
        cc = colorRampPalette(c('blue', 'yellow', 'red'))
    if (colorsheme == "rwb")
        cc = colorRampPalette(c('blue', 'white', 'red'))
    if (colorsheme == "ryg")
        cc = colorRampPalette(c('green', 'yellow', 'red'))
    if (colorsheme == "rwg")
        cc = colorRampPalette(c('green', 'white', 'red'))
    if (colorsheme == "rw")
        cc = colorRampPalette(c('white', 'red'))
    if (colorsheme == "bw")
        cc = colorRampPalette(c('white', 'blue'))
    if (colorsheme == "gw")
        cc = colorRampPalette(c('white', 'green'))
    return(cc(n))
}
#' @param dataset a dataframe of KEGG Enrichment output, with these columns:'ID', 'KEGGPathway', 'Count', 'all', 'P'
#'
#' @param MainTitle main title
#' @param top  top n kegg pathways used to plot
#' @param ColorScheme  "ryb", 'rwb', 'ryg', 'rwg', 'rw', 'bw', 'gw'
#' @param ColorReverse  whether to reverse color order
#' @param xangle   angle of x axis text
#'
#' @title KEGG enrichment bubble plot.
#' @description  do KEGG pathway enrichment bubble plot using KEGGenrich result.
#' @import ggplot2
#' @importFrom magrittr %>%
#' @importFrom scales label_wrap
#' @export
KEGGbubble = function(dataset,
                      MainTitle = "",
                      TrimName = F,
                      NameWidth = 50,
                      top = 10,
                      ColorScheme = "rw",
                      ColorReverse = F,xangle=45) {
    checkParams(ColorScheme,
                c("ryb", 'rwb', 'ryg', 'rwg', 'rw', 'bw', 'gw'),
                'ColorScheme')
    col = colorScheme(colorsheme = ColorScheme, n = 100)
    if (ColorReverse) {
        col = rev(col)
    }
    col = col[15:85]

    colnames(dataset) = c('ID', 'KEGGPathway', 'Count', 'all', 'P')

    if(TrimName){
        dataset$KEGGPathway=strtrim(dataset$KEGGPathway,NameWidth)
        if(length(unique(dataset$KEGGPathway))!=length(dataset$KEGGPathway)){
            stop(call. = F,"Pathway名称重复\n请设置参数TrimName=F，或调大参数NameWidth后重试！",'\n')
        }
    }

    k = sort(dataset$P,
             index.return = T,
             decreasing = F)$ix
    dataset = dataset[k, ]
    if (nrow(dataset) > top) {
        dataset = dataset[1:top, ]
    }
    dataset$P = -log10(dataset$P)
    dataset$GeneRatio = dataset$Count / dataset$all * 100
    k = sort(dataset$Count,
             index.return = T,
             decreasing = F)$ix
    dataset$KEGGPathway=unlist(lapply(dataset$KEGGPathway, function(x) {n=strsplit(x, split=";");n[[1]][1]}))
    dataset$KEGGPathway = factor(dataset$KEGGPathway, levels = dataset$KEGGPathway[k])
    if (max(dataset$Count) < 10) {
        Xbreaks=seq(0, ceiling(max(dataset$Count)), by = 1)
    } else if (max(dataset$Count) < 30) {
        Xbreaks=seq(0, ceiling(max(dataset$Count)), by = 5)
    } else{
        Xbreaks=seq(0, ceiling(max(dataset$Count)), by = 10)
    }
    if(length(Xbreaks)>10){
        Xbreaks=pretty(c(0,ceiling(max(dataset$Count))),n=10)
    }
    drawbreaks = pretty(c(floor(min(
        dataset$GeneRatio
    )), ceiling(max(
        dataset$GeneRatio
    ))), min.n = 3)
    g = ggplot(dataset, aes(Count, KEGGPathway)) +
        geom_point(
            mapping = aes(size = GeneRatio, color = P),
            show.legend = T,
            stroke = 1
        ) +
        theme_bw() +
        scale_size(breaks = drawbreaks ,
                   labels = drawbreaks,
                   range = c(5, 10)) +
        theme(
            plot.title = element_text(hjust = 0.5, colour = "black"),
            axis.text.y = element_text(size = 13, colour = "black"),
            axis.title = element_text(size = 14, colour = "black"),
            axis.text.x = element_text(size=13,angle = xangle,colour='black',hjust = 1),
            legend.position = "right"
        ) +
        scale_color_gradientn(colours = col) +
        labs(
            title = MainTitle,
            size = "Gene Ratio (%)",
            y = "",
            x = 'Count',
            color = expression(-log[10] * italic(P))
        ) +
        scale_y_discrete(position = "left",labels=scales::label_wrap(NameWidth)) +
        scale_x_continuous(limits = c(0, max(dataset$Count)),
                           breaks =Xbreaks)
    return(g)
}

#' GO enrichment bubble plot
#'
#' @param dataset a dataframe of GOEnrich output, with these columns: GO.ID,Term,Annotated,Significant,AllDEG,Pvalue,Class,gene
#' @param MainTitle Plot Title
#' @param onlySig whether plot only significant Terms
#' @param useFDR  whether use FDR to define significant
#' @param cut     p value significant threshold, if useFDR is set, FDR pvalue is used
#' @param top     top n terms used to plot in each BP/MF/CC
#' @param TrimTerm  whether to trim long term name to specific length
#' @param TermWidth Term longest width, if TrimTerm is set, long Terms will be trimmed, otherwise wrapped into specific length
#' @param ColorScheme  Color scheme used to plot pvalue, 'ryb', 'rwb', 'ryg', 'rwg', 'rw'(default), 'bw', 'gw'
#' @param ColorReverse  whether reverse the color order
#'
#' @import ggplot2
#' @importFrom  dplyr group_by
#' @importFrom  dplyr slice_max
#' @importFrom  dplyr filter
#' @importFrom  scales label_wrap
#' @importFrom magrittr %>%
#' @export
GObubble = function(dataset,
                    MainTitle = "",
                    onlySig = T,
                    useFDR=F,
                    cut = 0.05,
                    top = 10,
                    TrimTerm=T,
                    TermWidth=60,
                    ColorScheme = "rw",
                    ColorReverse = F) {
    checkParams(ColorScheme,
                c("ryb", 'rwb', 'ryg', 'rwg', 'rw', 'bw', 'gw'),
                'ColorScheme')

    if(useFDR){
        dataset$Pvalue=p.adjust(dataset$Pvalue,method="BH")
    }
    if (onlySig) {
        dataset = dataset %>% filter(Pvalue < cut)
    }
    dataset$Pvalue = -log10(dataset$Pvalue)
    dataset = dataset %>% group_by(Class) %>%
        slice_max(order_by = Pvalue, n = top) %>%as.data.frame


    if(TrimTerm){
        #   dataset$Term2=sapply(strsplit(dataset$Term,split = SplitSep),function(x) x[1])
        dataset$Term=strtrim(dataset$Term,TermWidth)
        if(length(unique(dataset$Term))!=length(dataset$Term)){
            stop(call. = F,"Term名称重复\n请设置参数TrimTerm=F，或调大参数TermWidth后重试！",'\n')
        }
    }

    if(nrow(dataset%>%filter(Class=='BP'))>top){
        warning(call. = F,paste("BP结果存在相等的p值，参与绘图的BP Term可能大于",top))
    }
    if(nrow(dataset%>%filter(Class=='MF'))>top){
        warning(call. = F,paste("MF结果存在相等的p值，参与绘图的MF Term可能大于",top))
    }
    if(nrow(dataset%>%filter(Class=='CC'))>top){
        warning(call. = F,paste("CC结果存在相等的p值，参与绘图的CC Term可能大于",top))
    }
    dataset$GeneRatio = dataset$Significant / dataset$AllDEG * 100
    order <- order(dataset$GeneRatio)
    dataset$Term <-factor(dataset$Term, levels = dataset$Term[order])
    xLabel = 'Count'
    yLabel = ''
    sizeLabel = 'Gene Ratio (%)'
    colorLabel = expression('-log'[10] * 'P-value')
    DotRange = c(2, 10)
    col = colorScheme(colorsheme = ColorScheme, n = 100)
    if (ColorReverse) {
        col = rev(col)
    }
    col = col[15:85]
    colorbreaks = pretty(dataset$Pvalue, 5)
    sizebreaks = pretty(dataset$GeneRatio, 5)
    g <- ggplot(data = dataset, mapping = aes(Significant, Term)) +
        theme_bw() +
        theme(aspect.ratio = 1.1)+
        geom_point(aes(size = GeneRatio, color = Pvalue),
                   show.legend = T,
                   stroke = 1) +
        facet_wrap(~Class, scales = 'free_y',ncol=1,strip.position = 'right') +
        scale_size(
            name = sizeLabel,
            breaks = sizebreaks,
            labels = sizebreaks,
            range = DotRange,
            guide = "legend"
        ) +
        scale_color_gradientn(colours = col, breaks = colorbreaks) +
        scale_y_discrete(position = "left",labels=scales::label_wrap(TermWidth)) +
        labs(
            title = MainTitle,
            x = xLabel,
            y = yLabel,
            size = sizeLabel,
            color = colorLabel
        ) +
        theme(
            plot.title = element_text(hjust = 0.5),panel.spacing = unit(0, "cm"),

            axis.text = element_text(colour = "black", size = 11),
            axis.text.y = element_text(lineheight = 0.65),
            legend.text = element_text(size = 10),
            legend.title = element_text(size = 10),
            legend.position = "right",
            strip.background = element_blank(),
            strip.text = element_text(size = 15, face = 'bold')
        )
    return(g)
}

#' @param dat  deg analysis result dataframe with ID,FCcut,log2FoldChange,Pvalue or FDR
#'
#' @param useFDR whether to use FDR to define significance
#' @param MainTitle main title of plot
#' @param cut   P value significant cut , if useFDR is set, FDR is used
#' @param FCcut  FoldChange cut significant cut
#' @param xlim   x axis limits
#' @param showlabel  whether to show top significant gene names
#' @param showlabel.num  top number
#' @param point.size  point size
#'
#' @title Volcano plot for DEG analysis results.
#' @description  do volcano plot.
#' @import ggplot2
#' @importFrom  dplyr filter
#' @importFrom  dplyr select
#' @importFrom  dplyr mutate
#' @importFrom magrittr %>%
#' @importFrom ggrepel geom_text_repel
#' @export
VolcanoPlot = function(dat,
                       useFDR = T,
                       MainTitle = "",
                       cut = 0.05,
                       FCcut = 2,
                       xlim = 10,
                       showlabel = "no",
                       showlabel.num = 10,
                       point.size=2) {
    checkParams(showlabel,
                c("allSig", 'topPvalue', 'topFC', "no"),
                'showlabel')
    dat=na.omit(dat)
    if(!('Pvalue' %in% colnames(dat)) ){
        if('FDR' %in% colnames(dat)  ){
            dat$Pvalue=dat$FDR
        }else{
            stop(call. = F,"错误：未找到FDR或Pvalue列",'\n')
        }
    }
    if (useFDR) {
        input <- dat %>%
            mutate(sig = ifelse(
                dat$FDR < cut & abs(dat$log2FoldChange) >= log2(FCcut),
                ifelse(dat$log2FoldChange >= 0,
                       "Up",
                       "Down"),
                "Non-Significant"
            ))
    } else{
        input <- dat %>%
            mutate(sig = ifelse(
                dat$Pvalue < cut & abs(dat$log2FoldChange) >= log2(FCcut),
                ifelse(dat$log2FoldChange >= 0,
                       "Up",
                       "Down"),
                "Non-Significant"
            ))
    }
    input = input %>%
        mutate(lim = ifelse(
            dat$log2FoldChange > xlim,
            'up-censored',
            ifelse(dat$log2FoldChange < (-xlim),
                   'down-censored',
                   'origin')
        ))

    input[input$log2FoldChange > xlim, 'log2FoldChange'] = xlim
    input[input$log2FoldChange < (-xlim), 'log2FoldChange'] = (-xlim)

    shapevector = c(
        "down-censored" = "\u25C4",
        "origin" = "\u25CF",
        "up-censored" = "\u25BA"
    )
    colorVector = c("Down" = "blue",
                    "Non-Significant" = "black",
                    "Up" = "red")
    g <-
        ggplot(input, aes(log2FoldChange, -log10(Pvalue))) + geom_point(aes(col = sig, shape = lim), size = point.size)
    g <- g + theme_bw()
    g <-
        g + scale_shape_manual(values = shapevector, guide = 'none')
    g <- g + scale_color_manual(values = colorVector)
    g <- g + geom_vline(xintercept = 0,
                        color = 'grey',
                        size = 1)
    g <-
        g + theme(axis.text = element_text(size = 10, colour = "black"))
    g <-
        g + theme(legend.position = "right", legend.title = element_blank())
    g <-
        g + labs(title = MainTitle, x = "log2(Fold Change)", y = "-log10(P-value)")
    g <- g + theme(plot.title = element_text(hjust = 0.5))

    if (showlabel != "no") {
        if (showlabel == "allSig") {
            labeldat = input[input$sig != "Non-Significant", ]
        } else if (showlabel == "topPvalue") {
            labeldat = head(input[order(input$Pvalue), ], showlabel.num)
        } else if (showlabel == "topFC") {
            labeldat = head(input[order(abs(input$log2FoldChange), decreasing = T), ], showlabel.num)
        } else{
            stop('ERROR')
        }
        g <-
            g + geom_text_repel(
                data = labeldat,
                aes(label = rownames(labeldat), color = sig),
                size = 3,
                max.overlaps = 100,
                show.legend = F
            )
    }
    return(g)
}

#' @param dat  deg analysis result dataframe with ID,FCcut,log2FoldChange,Pvalue or FDR
#'
#' @param useFDR whether to use FDR to define significance
#' @param MainTitle main title of plot
#' @param cut   P value significant cut , if useFDR is set, FDR is used
#' @param FCcut  FoldChange cut significant cut
#' @param ylim    y axis limits
#' @param showlabel  whether to show top significant gene names
#' @param showlabel.num  top number
#' @param point.size  point size
#'
#' @title MA plot for DEG analysis results.
#' @description  do MA plot.
#' @import ggplot2
#' @importFrom  dplyr filter
#' @importFrom  dplyr select
#' @importFrom  dplyr mutate
#' @importFrom magrittr %>%
#' @importFrom ggrepel geom_text_repel
#' @export
MAPlot = function(dat,
                  useFDR = T,
                  MainTitle = "",
                  cut = 0.05,
                  FCcut = 2,
                  ylim = 10,
                  showlabel = "no",
                  showlabel.num = 10,
                  point.size=2) {
    dat=na.omit(dat)
    if(!('Pvalue' %in% colnames(dat)) ){
        if('FDR' %in% colnames(dat)  ){
            dat$Pvalue=dat$FDR
        }else{
            stop(call. = F,"错误：未找到FDR或Pvalue列",'\n')
        }
    }
    if (useFDR) {
        input <- dat %>%
            mutate(sig = ifelse(
                dat$FDR < cut & abs(dat$log2FoldChange) >= log2(FCcut),
                ifelse(dat$log2FoldChange >= 0,
                       "Up",
                       "Down"),
                "Non-Significant"
            ))
    } else{
        input <- dat %>%
            mutate(sig = ifelse(
                dat$Pvalue < cut & abs(dat$log2FoldChange) >= log2(FCcut),
                ifelse(dat$log2FoldChange >= 0,
                       "Up",
                       "Down"),
                "Non-Significant"
            ))
    }
    input = input %>%
        mutate(lim = ifelse(
            dat$log2FoldChange > ylim,
            'up-censored',
            ifelse(dat$log2FoldChange < (-ylim),
                   'down-censored',
                   'origin')
        ))

    shapevector = c(
        "down-censored" = "\u25BC",
        "origin" = "\u25CF",
        "up-censored" = "\u25B2"
    )
    colorvector = c("Down" = "blue",
                    "Non-Significant" = "black",
                    "Up" = "red")
    input[input$log2FoldChange > ylim, 'log2FoldChange'] = ylim
    input[input$log2FoldChange < (-ylim), 'log2FoldChange'] = -ylim
    g <-
        ggplot(input, aes(log10(MeanExpression), log2FoldChange)) + geom_point(aes(color =
                                                                                       sig, shape = lim), size = point.size) +
        theme_bw() +
        scale_shape_manual(values = shapevector, guide = 'none') +
        scale_color_manual(values = colorvector) +
        geom_hline(yintercept = 0,
                   color = 'grey',
                   size = 1) +
        theme(axis.text = element_text(size = 10, colour = "black")) +
        theme(legend.position = "right", legend.title = element_blank()) +
        labs(title = MainTitle, x = "Normalized Expression", y = "log2(Fold Change)") +
        theme(plot.title = element_text(hjust = 0.5))

    if (showlabel != "no") {
        if (showlabel == "allSig") {
            labeldat = input[input$sig != "Non-Significant", ]
        } else if (showlabel == "topPvalue") {
            labeldat = head(input[order(input$Pvalue), ], showlabel.num)
        } else if (showlabel == "topFC") {
            labeldat = head(input[order(abs(input$log2FoldChange), decreasing = T), ], showlabel.num)
        } else{
            stop('ERROR')
        }
        g <-
            g + geom_text_repel(
                data = labeldat,
                aes(label = rownames(labeldat), color = sig),
                size = 3,
                max.overlaps = 100,
                show.legend = F
            )
    }
    return(g)
}



#' @param data  input trait dataframe, with traits as columns
#'
#' @title Frequently used stats for trait analysis.
#' @description  Frequently used stats for trait analysis.
#' @importFrom  fBasics basicStats
#' @importFrom  vegan diversity
#' @export
multiTraitStat<-function(data){

    stats=as.data.frame(t(basicStats(data)[c('nobs','NA','Minimum','Maximum','Mean','Median','Stdev','Skewness','Kurtosis'),]))

    DiversityIndex=sapply(apply(data,2,na.omit),diversity)
    CV=stats$Stdev/stats$Mean
    out=data.frame(stats,CV,DiversityIndex)
    return(out)
}

#' @param data input trait dataframe, with traits as columns
#'
#' @title  trait Plot.
#' @description  Trait plot with histograms/correlation/density.
#' @export
multiTraitPlot<-function(data){
    custom_cor <- function(x, y,cutoff=0.05) {
        usr <- par("usr"); on.exit(par(usr))
        par(usr = c(0, 1, 0, 1))
        ct <- cor.test(x,y,method="pearson")
        sig <- symnum(ct$p.value,
                      corr = FALSE,
                      na = FALSE,
                      cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1),
                      symbols = c("***", "**", "*", ".", "NS"))
        r <- ct$estimate
        rt <- format(r, digits=2)[1]
        cex=1.5
        if(ct$p.value<cutoff){
            color='red'
        }else{color='grey50'}
        text(.5,.4,rt,cex=cex,col =color)
        text(.5,.6,sig,cex=cex,col =color)
    }

    custom_smooth <- function(x,y) {
        points(x,y,col='black',pch=20)
        abline(lm(y~x))
    }
    custom_hist <- function(x, ...)
    {
        usr <- par("usr"); on.exit(par(usr))
        par(usr = c(usr[1:2], 0, 1.5) )
        h <- hist(x, plot = FALSE)
        breaks <- h$breaks; nB <- length(breaks)
        y <- h$counts; y <- y/max(y)
        rect(breaks[-nB], 0, breaks[-1], y, col = "green", ...)
        lines(density(x,na.rm=T),lwd=1,col="chocolate3")
    }
    pairs(data,upper.panel=custom_cor,diag.panel=custom_hist,lower.panel=custom_smooth)
}

#' @title get All Children  .
#' @description get all children  .
#' @importFrom GO.db GOBPCHILDREN
#' @importFrom GO.db GOMFCHILDREN
#' @importFrom GO.db GOCCCHILDREN
#' @importFrom AnnotationDbi mget
#' @import topGO
getAllChildren <- function(goids,class='BP')
{
  if(class == "BP"){
    ans <- unique(unlist(AnnotationDbi::mget(goids, GOBPCHILDREN), use.names=FALSE))
  }else if(class == "CC"){
    ans <- unique(unlist(AnnotationDbi::mget(goids, GOCCCHILDREN), use.names=FALSE))
  }else if(class == "MF"){
    ans <- unique(unlist(AnnotationDbi::mget(goids, GOMFCHILDREN), use.names=FALSE))
  }else{
    stop('error\n')
  }
  ans <- ans[!is.na(ans)]
}

#' @param deglist DEG name vector
#'
#' @param taxonid taxon id , use GOdbInfo() to check available ids
#' @param customMapping    path of the custom GO annotation file, which annotates each gene in a row, separate GOs with comma(,)
#' @param eggnog   path of the custom GO annotation using eggNOG database
#' @param pannzer2  path of the custom GO annotation using pannzer2 database
#' @param customTable path of the custom GO annotation file, which annotates two-column pairwise gene-GO per row
#'
#' @title GO barplot .
#' @description Barplot of secondary GO terms.
#' @importFrom dplyr select
#' @importFrom reshape2 melt
#' @importFrom magrittr %>%
#' @import ggplot2
#' @import topGO
#' @export
GOBar=function(deglist,taxonid=NULL,customMapping=NULL,eggnog=NULL,pannzer2=NULL,customTable=NULL){
  if(!is.null(taxonid)){
    checkParams(taxonid, names(godb), 'taxon')
    geneID2GO = godb[[taxonid]][['db']]
  }else  if(!is.null(customMapping)){
    geneID2GO<-readMappings(file=customMapping)
  }else  if(!is.null(eggnog)){
    data = read_delim(
        eggnog,
        comment = '##',
        delim = "\t",
        na = '-',
        col_names = T,
        col_types =cols(.default = col_character()))%>%
       dplyr::select(c('#query', 'GOs'))%>% na.omit()
    geneID2GO <- strsplit(data$GOs, ",")
    names(geneID2GO) = data$`#query`
  }else  if(!is.null(pannzer2)){
    geneID2GO = PANNZERres2GOdb(pannzer2)
  }else  if(!is.null(customTable)){
    geneID2GOtable=read.delim(customTable,header=T,sep="\t")
    checkParams(colnames(geneID2GOtable),c('GeneID','GOID'),string = "geneID2GOtable列名错误")
    geneID2GO=lapply(unique(geneID2GOtable$GeneID), function (x) geneID2GOtable$GOID[geneID2GOtable$GeneID==x])
    names(geneID2GO)=unique(geneID2GOtable$GeneID)
  }else{
    stop("error no GO database selected!\n")
  }
  GO2geneID <- topGO::inverseList(geneID2GO)
  allgenes=names(geneID2GO)
  universe = factor(as.integer(allgenes %in% deglist))
  names(universe) = allgenes
  bp_terms <- getAllChildren("GO:0008150",'BP')
  cc_terms <- getAllChildren("GO:0005575",'CC')
  mf_terms <- getAllChildren("GO:0003674",'MF')

  bpObj <-new("topGOdata",nodeSize = 1,ontology = 'BP',allGenes = universe,annot = annFUN.gene2GO,gene2GO = geneID2GO)
  mfObj <-new("topGOdata",nodeSize = 1,ontology = 'MF',allGenes = universe,annot = annFUN.gene2GO,gene2GO = geneID2GO)
  ccObj <-new("topGOdata",nodeSize = 1,ontology = 'CC',allGenes = universe,annot = annFUN.gene2GO,gene2GO = geneID2GO)

  allbpGO = topGO::genesInTerm(bpObj,bp_terms)
  allccGO = topGO::genesInTerm(ccObj,cc_terms)
  allmfGO = topGO::genesInTerm(mfObj,mf_terms)

  stat1=data.frame(ID=AnnotationDbi::Term(names(allbpGO)),
                   allGene=sapply(names(allbpGO),
                                  function(x){length(allbpGO[[x]])}),
                   DEG=sapply(names(allbpGO),
                              function(x){length(intersect(deglist,allbpGO[[x]]))}),
                   Class='BP'
  )
  stat2=data.frame(ID=AnnotationDbi::Term(names(allmfGO)),
                   allGene=sapply(names(allmfGO),
                                  function(x){length(allmfGO[[x]])}),
                   DEG=sapply(names(allmfGO),
                              function(x){length(intersect(deglist,allmfGO[[x]]))}),
                   Class='MF'
  )
  stat3=data.frame(ID=AnnotationDbi::Term(names(allccGO)),
                   allGene=sapply(names(allccGO),
                                  function(x){length(allccGO[[x]])}),
                   DEG=sapply(names(allccGO),
                              function(x){length(intersect(deglist,allccGO[[x]]))}),
                   Class='CC'
  )
  plotdata=rbind(stat1,stat2,stat3)%>%as.data.frame()%>%reshape2::melt(id.vars=c('ID','Class'))

  ord=order(plotdata[plotdata$variable=='allGene','value'],decreasing = T)
  plotdata$ID=factor(plotdata$ID,levels=plotdata[plotdata$variable=='allGene','ID'][ord])
  plotdata[plotdata$value==0,'value']=1
  usecolors=c("#A6CEE3", "#1F78B4" ,"#B2DF8A" ,"#33A02C", "#FB9A99", "#E31A1C")
  names(usecolors)=c('BP DEG','BP allGene','CC DEG','CC allGene','MF DEG','MF allGene')
  g=ggplot(plotdata,aes(ID,value))+
    geom_bar(aes(fill=paste(Class,variable)),stat='identity',position = 'dodge')+
    labs(x="",y="",fill="")+
    facet_grid(.~Class,scales = 'free',space='free')+
    scale_y_continuous(trans='log10',expand=c(0,0))+
    scale_fill_manual(values=usecolors)+
    theme_bw()+
    theme(
      strip.background = element_blank(),
      strip.placement = 'outside',
      legend.position="right",
      axis.text.x = element_text(angle=45,hjust = 1)
    )+
    guides(
      fill = guide_legend(ncol = 2,byrow = T,label.position = 'left')
    )

  return(g)
}
