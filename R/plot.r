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
#' @title KEGG enrichment bubble plot.
#' @description  do KEGG pathway enrichment bubble plot using KEGGenrich result.
#' @import ggplot2
#' @importFrom magrittr %>%
#' @export
KEGGbubble = function(dataset,
                      MainTitle = "",
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
        scale_y_discrete(position = "left") +
        scale_x_continuous(limits = c(0, max(dataset$Count)),
                           breaks =Xbreaks)
    return(g)
}

#' @title GO enrichment bubble plot.
#' @description  do GO  enrichment bubble plot using GOEnrich result.
#' @import ggplot2
#' @importFrom  dplyr group_by
#' @importFrom  dplyr slice_max
#' @importFrom  dplyr filter
#' @importFrom magrittr %>%
#' @export
GObubble = function(dataset,
                    MainTitle = "",
                    onlySig = T,
                    useFDR=F,
                    cut = 0.05,
                    top = 10,
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
    #   dataset$Term=strtrim(dataset$Term,maxTermLen)
    dataset$Pvalue = -log10(dataset$Pvalue)
    dataset = dataset %>% group_by(Class) %>% slice_max(order_by = Pvalue, n = top) %>%
        as.data.frame

    dataset$GeneRatio = dataset$Significant / dataset$AllDEG * 100
    order <- order(dataset$GeneRatio)
    dataset$Term <-
        factor(dataset$Term, levels = dataset$Term[order])
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
        geom_point(aes(size = GeneRatio, color = Pvalue),
                   show.legend = T,
                   stroke = 1) +
        facet_grid(Class ~ ., scales = "free", space = "free") +
        scale_size(
            name = sizeLabel,
            breaks = sizebreaks,
            labels = sizebreaks,
            range = DotRange,
            guide = "legend"
        ) +
        scale_color_gradientn(colours = col, breaks = colorbreaks) +
        scale_y_discrete(position = "left") +
        labs(
            title = MainTitle,
            x = xLabel,
            y = yLabel,
            size = sizeLabel,
            color = colorLabel
        ) +
        theme_bw() +
        theme(
            plot.title = element_text(hjust = 0.5),
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

#' @title Frequently used stats for trait.
#' @description  Frequently used stats for trait.
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
