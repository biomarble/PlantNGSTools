#' @title GOdbInfo.
#' @description  Show GO database information.
#' @importFrom magrittr %>%
#' @export
GOdbInfo = function() {
    info = sapply(names(godb), function(x)
        godb[[x]][['version']]) %>% as.data.frame()
    info2 = sapply(names(godb), function(x)
        godb[[x]][['update']]) %>% as.data.frame()
    info = cbind(rownames(info), info, info2) %>% as.data.frame()
    colnames(info) = c('taxon', 'ReferenceGenomeVersion', 'update')
    print(info, row.names = F)
}
#' @title KEGGdbInfo.
#' @description  Show KEGG database information.
#' @importFrom magrittr %>%
#' @export
KEGGdbInfo = function() {
    info = sapply(names(keggdb), function(x)
        unique(keggdb[[x]][['Version']])) %>% as.data.frame()
    info2 = sapply(names(keggdb), function(x)
        unique(keggdb[[x]][['update']])) %>% as.data.frame()
    info = cbind(rownames(info), info, info2) %>% as.data.frame()
    colnames(info) = c('taxon', 'ReferenceGenomeVersion', 'update')
    print(info, row.names = F)
}
fetchGOs = function(ensembldataset, useArchive = F) {
    if (useArchive) {
        ensembl = useMart(host = 'nov2020-plants.ensembl.org',
                          biomart = "plants_mart",
                          dataset = ensembldataset)
    } else{
        ensembl <-
            useEnsemblGenomes(biomart = "plants_mart", dataset = ensembldataset)
    }
    res = getBM(attributes = c('ensembl_gene_id', 'go_id'),
                mart = ensembl) %>% filter(go_id != "")
    geneID2GO = lapply(unique(res$ensembl_gene_id), function (x)
        paste(res$go_id[res$ensembl_gene_id == x]))
    names(geneID2GO) = unique(res$ensembl_gene_id)
    return(geneID2GO)
}

fetchPathways = function(taxon, ensembldataset, useArchive = F) {
    parseName = function(input,
                         sep = ":",
                         retind = 2) {
        res = strsplit(input, split = sep) %>% unlist() %>% matrix(., ncol = 2, byrow =
                                                                       T)
        return(res[, retind])
    }
    taxonPathway = function(taxon) {
        usepathway = keggList('pathway', organism = taxon)
        pathwayID = names(usepathway) %>% parseName()
        usepathway = cbind(
            pathwayID,
            unname(usepathway) %>% stringi::stri_replace_last_regex(pattern = " - ", replacement = '::') %>%
                parseName(sep = "::", retind = 1)
        ) %>% as.data.frame()
        colnames(usepathway) = c('PathwayID', 'Pathway')
        l = keggLink("pathway", taxon)
        id = names(l)
        path = l %>% unname() %>% parseName()
        res = left_join(as.data.frame(cbind(
            ID = id, PathwayID = path
        )), usepathway, by = "PathwayID")
        return(res)
    }
    id2entrez <- function(query, ensembldataset) {
        if (useArchive) {
            ensembl = useMart(host = 'nov2020-plants.ensembl.org',
                              biomart = "plants_mart",
                              dataset = ensembldataset)
        } else{
            ensembl <-
                useEnsemblGenomes(biomart = "plants_mart", dataset = ensembldataset)
        }
        id2entrez = getBM(
            attributes = c('ensembl_gene_id', 'entrezgene_id'),
            mart = ensembl
        ) %>% filter(!is.na(entrezgene_id))
        id = id2entrez[match(query %>% parseName(), id2entrez$entrezgene_id), 1]
        res = cbind(entrezID = query, GeneID = id) %>% as.data.frame()
        return(res)
    }

    pathways = taxonPathway(taxon)
    res = pathways
    print(head(res))
    entrezIDorNot = "NULL"
    while (!(entrezIDorNot %in% c('y', 'n'))) {
        entrezIDorNot <- readline(prompt = "is it Entrez ID ? (y/n)")
    }
    if (entrezIDorNot == "y") {
        ID2Entrez = id2entrez(pathways$ID, ensembldataset)
        res = left_join(pathways, ID2Entrez, by = c("ID" = "entrezID")) %>%
            filter(!is.na(GeneID)) %>% select(c('GeneID', 'ID' , 'PathwayID', 'Pathway')) %>%
            unique()
        head(res)
    } else{
        res$GeneID = res$ID %>% parseName()
        res = res %>% as.data.frame() %>% select(c('GeneID', 'ID' , 'PathwayID', 'Pathway')) %>%
            unique()
    }
    print(head(res))
    return(res)
}

#' @title KEGG pathway enrichment tools.
#' @description  do KEGG pathway enrichment by a list of degs.
#' @import R2HTML
#' @importFrom magrittr %>%
#' @export
KEGGenrich <- function(deglist,
                       taxon,
                       outdir = NULL,
                       outprefix,
                       fdr = FALSE) {
    checkParams(taxon, names(keggdb), 'taxon')

    if (is.null(outdir)) {
        outdir = getwd()
    }
    if (!dir.exists(outdir)) {
        dir.create(outdir)
    }
    pathinfo = keggdb[[taxon]]
    pathway = unique(pathinfo[, c('PathwayID', 'Pathway')])
    deggene = unique(deglist)
    allgene = unique(pathinfo$GeneID)
    deginfo = pathinfo[pathinfo[, 1] %in% deggene, ]
    allinfo = pathinfo[pathinfo[, 1] %in% allgene, ]
    k = length(unique(deginfo[, 1]))
    N = length(unique(allinfo[, 1]))
    p <- data.frame()
    urlColor <- NULL
    geneLink <- NULL
    for (i in seq(1, nrow(pathway))) {
        pid = pathway[i, 'PathwayID']
        allInPath = allinfo[allinfo[, 'PathwayID'] %in% pid, ]
        m = length(unique(allInPath[, 1]))
        degInPath = deginfo[deginfo[, 'PathwayID'] %in% pid, ]
        x = length(unique(degInPath[, 1]))
        p[i, 1] = x
        p[i, 2] = m
        p[i, 3] = phyper(x - 1, m, N - m, k, lower.tail = FALSE)
        urlColor[i] = apply(as.matrix(paste(
            "/", t(degInPath[, 'ID']), '%09red', sep = ""
        )), 2, paste, collapse = "")
        Link = paste('<a href=',
                     shQuote(
                         paste('https://www.genome.jp/entry/', degInPath[, 'ID'], sep = "")
                     ),
                     '/',
                     '>',
                     degInPath[, 'GeneID'],
                     '</a>',
                     sep = "")
        if (x == 0) {
            geneLink[i] = "--"
        } else{
            geneLink[i] = apply(as.matrix(Link), 2, paste, collapse = ", ")
        }
    }
    output = cbind(pathway, p)
    colnames(output) = c('ID',
                         'Pathway',
                         'DEGsInPathway',
                         'GenesInPathway',
                         'Pvalue')
    if (fdr) {
        fdr = p.adjust(p[, 3], method = "BH")
        output = cbind(output, FDR = fdr)
    }
    ind = order(output[, "Pvalue"])
    output = output[ind, ]
    urlColor = urlColor[ind]
    DEGs = geneLink[ind]
    output2 = cbind(output, DEGs)
    ind = output$DEGsInPathway > 0
    output = output[ind, ]
    urlColor = urlColor[ind]
    DEGs = geneLink[ind]
    output2 = output2[ind, ]
    write.table(
        output,
        file = paste(outdir, "/", outprefix, ".tsv", sep = ""),
        quote = F,
        row.names = F,
        col.names = T,
        sep = "\t"
    )
    urlTitles <- output[, 'ID']
    url = paste(
        'https://www.genome.jp/kegg-bin/show_pathway?',
        gsub("ko", "map", urlTitles),
        urlColor,
        sep = ""
    )
    htmlOUT <-
        transform(output2,
                  ID = paste('<a href = ', shQuote(url), '/', '>', urlTitles, '</a>'))
    target <-
        HTMLInitFile(
            Title = "KEGG Pathway Ernichment Results",
            outdir = outdir,
            filename = paste0(outprefix),
            BackGroundColor = "#f8f8f8",
            useGrid = T,
            useLaTeX = F,
            HTMLframe = F
        )
    write(
        "<style>table{border-collapse:collapse;width:1280px;}a:link,a:visited {color: #3366FF;text-decoration: none; }a:hover,a:active {color: #ff3366;text-decoration: none; };</style>",
        target,
        append = T
    )
    HTML(
        htmlOUT,
        file = target,
        innerBorder = 1,
        row.names = F,
        digits = 3
    )
    HTMLEndFile()
    cat('Enrichment Done! Result files saved:\n',
        outdir,
        '/',
        outprefix,
        '.tsv\n',
        sep = "")
    cat(outdir, '/', outprefix, '.html\n', sep = "")
    return(output)
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
                      ColorReverse = F) {
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
        xstep = 1
    } else if (max(dataset$Count) < 30) {
        xstep = 5
    } else{
        xstep = 10
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
            axis.text = element_text(size = 13, colour = "black"),
            axis.title = element_text(size = 14, colour = "black"),
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
                           breaks = seq(0, ceiling(max(dataset$Count)), by = xstep))
    return(g)
}

#' @title GO enrichment .
#' @description  do GO enrichment by a list of DEGs.
#' @import ggplot2
#' @import GO.db
#' @importFrom  dplyr select
#' @import topGO
#' @export
GOEnrich = function(deglist,
                    taxon,
                    fdr = FALSE,
                    outdir = NULL,
                    outprefix) {
    checkParams(taxon, names(godb), 'taxon')
    if (is.null(outdir)) {
        outdir = getwd()
    }
    if (!dir.exists(outdir)) {
        dir.create(outdir)
    }
    geneID2GO = godb[[taxon]][['db']]
    version = godb[[taxon]][['version']]
    GOenrichsub <- function(deglist, geneID2GO, class, output) {
        cat('Processing', class, 'class...')
        geneList <- rep(1, length(geneID2GO))
        names(geneList) <- names(geneID2GO)
        geneList[match(deglist, names(geneList))] = 0
        sampleGOdata <-
            suppressMessages(
                new(
                    "topGOdata",
                    nodeSize = 1,
                    ontology = class,
                    allGenes = geneList,
                    annot = annFUN.gene2GO,
                    gene2GO = geneID2GO,
                    geneSel = function(allScore) {
                        return(allScore < 0.01)
                    }
                )
            )
        cat('...')
        result <-
            suppressMessages(runTest(
                sampleGOdata,
                algorithm = "elim",
                statistic = "fisher"
            ))
        allRes <-
            GenTable(
                sampleGOdata,
                Pvalue = result,
                orderBy = "Pvalue",
                topNodes = attributes(result)$geneData[4],
                numChar = 1000
            )
        cat('...')

        gene = NULL
        GOs = sampleGOdata@graph@nodeData@data
        for (i in allRes$GO.ID) {
            e = paste0('g=names(GOs$`', i, '`$genes)')
            eval(parse(text = e))

            d = intersect(g, deglist)
            e = paste0('gene[\'', i, '\']=toString(d)')
            eval(parse(text = e))
        }
        allRes$gene = gene
        allRes = cbind(allRes,
                       AllDEG = length(geneList[geneList == 0]),
                       Class = class)
        cat('done!\n')
        return(allRes)
    }
    cat('genome version:', version, '\n')
    BP.res = GOenrichsub(deglist,
                         geneID2GO,
                         'BP',
                         paste0(outdir, '/', outprefix))
    MF.res = GOenrichsub(deglist,
                         geneID2GO,
                         'MF',
                         paste0(outdir, '/', outprefix))
    CC.res = GOenrichsub(deglist,
                         geneID2GO,
                         'CC',
                         paste0(outdir, '/', outprefix))
    all = rbind(BP.res, MF.res, CC.res)
    all$Pvalue = as.numeric(all$Pvalue)
    all = all %>% select(
        c(
            'GO.ID',
            'Term',
            'Annotated',
            'Significant',
            'AllDEG',
            'Pvalue',
            'Class',
            'gene'
        )
    )
    write.table(
        all[all$Pvalue < 0.05,],
        file = paste0(outdir, '/', outprefix, ".GO.significant.tsv"),
        sep = "\t",
        quote = FALSE,
        col.names = TRUE,
        row.names = FALSE
    )
    write.table(
        all,
        file = paste0(outdir, '/', outprefix, ".GO.all.tsv"),
        sep = "\t",
        quote = FALSE,
        col.names = TRUE,
        row.names = FALSE
    )
    cat(
        '\nGO enrichment Done!\nResult files saved at:\n',
        outdir,
        '/',
        outprefix,
        '.GO.significant.tsv\n',
        sep = ""
    )
    cat(outdir, '/', outprefix, '.GO.all.tsv\n', sep = "")
    return(all)
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
                    pCut = 0.05,
                    top = 10,
                    ColorScheme = "rw",
                    ColorReverse = F) {
    checkParams(ColorScheme,
                c("ryb", 'rwb', 'ryg', 'rwg', 'rw', 'bw', 'gw'),
                'ColorScheme')

    if (onlySig) {
        dataset = dataset %>% filter(Pvalue < pCut)
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

#' @title DEG analysis using samples with replicates.
#' @description  do DEG analysis using samples with replicates by DESeq2.
#' @import DESeq2
#' @importFrom  dplyr filter
#' @importFrom  dplyr select
#' @importFrom magrittr %>%
#' @importFrom tibble rownames_to_column
#' @export
DEGAnalysis_DESeq2 <-
    function (countMatrix,
              group,
              useFDR = F,
              cut = 0.05,
              FCcut = 2,
              control = NULL,
              treat = NULL,
              outdir = NULL) {
        checkGroup(countMatrix, group)
        if (is.null(outdir)) {
            outdir = getwd()
        }
        if (!dir.exists(outdir)) {
            dir.create(outdir)
        }
        colnames(group) = c('sample', 'group')
        allgroups = unique(group$group)

        while (!(control %in% allgroups)) {
            print(group)
            cat('可用分组名称：', paste(allgroups), '\n')
            control <- readline(prompt = "输入对照组名称: ")
        }
        allgroups = allgroups[!(allgroups %in% control)]
        if (length(allgroups) == 1) {
            treat = allgroups[1]
            cat('实验组名称（自动选择）：', treat, '\n')
        } else{
            while (!(treat %in% allgroups)) {
                cat('可用分组名称：', paste(allgroups), '\n')
                treat <- readline(prompt = "输入实验组名称: ")
            }
        }
        meta.data = group %>% filter(group %in% c(control, treat))
        meta.data$group = as.factor(meta.data$group)
        #    cat(meta.data$sample)
        countMatrix = countMatrix %>% select(meta.data$sample)
        dds <-
            DESeqDataSetFromMatrix(countMatrix, colData = meta.data, ~ group)
        dds <- dds[rowSums(counts(dds)) > 1, ]
        dds <- suppressMessages(DESeq(dds))
        res <- results(dds, contrast = c("group", treat, control))
        out = select(as.data.frame(res),
                     'baseMean',
                     'log2FoldChange',
                     'pvalue',
                     'padj')
        out = rownames_to_column(out, 'ID')
        out = filter(out, !is.na(padj))
        colnames(out) = c('ID', 'MeanExpression', 'log2FoldChange', 'Pvalue', 'FDR')
        if (useFDR) {
            diff = filter(out, FDR < cut, abs(log2FoldChange) >= abs(log2(FCcut)))
        } else{
            diff = filter(out, Pvalue < cut, abs(log2FoldChange) >= abs(log2(FCcut)))
        }
        write.table(
            out,
            paste0(outdir, '/', control, '.vs.', treat, '.AllGene.tsv'),
            row.names = F,
            quote = FALSE,
            sep = '\t'
        )
        write.table(
            diff,
            paste0(outdir, '/', control, '.vs.', treat, '.DEG.tsv'),
            row.names = F,
            quote = FALSE,
            sep = '\t'
        )
        write.table(
            diff$ID,
            paste0(outdir, '/', control, '.vs.', treat, '.DEGlist.tsv'),
            row.names = F,
            quote = FALSE,
            sep = '\t'
        )
        cat('差异基因分析结果见以下目录:\n')
        cat(paste0(outdir, '/', control, '.vs.', treat, '.AllGene.tsv'),
            '\n')
        cat(paste0(outdir, '/', control, '.vs.', treat, '.DEG.tsv'),
            '\n')
        cat(paste0(outdir, '/', control, '.vs.', treat, '.DEGlist.tsv'),
            '\n')
        filepath = NULL
        filepath$allgene = paste0(outdir, '/', control, '.vs.', treat, '.AllGene.tsv')
        filepath$deg = paste0(outdir, '/', control, '.vs.', treat, '.DEG.tsv')
        return(filepath)
    }
#' @title DEG analysis using samples without replicates.
#' @description  do DEG analysis using samples without replicates by EBSeq.
#' @import EBSeq
#' @importFrom  dplyr filter
#' @importFrom  dplyr select
#' @importFrom magrittr %>%
#' @importFrom tibble rownames_to_column
#' @export
DEGAnalysis_EBSeq <-
    function (countMatrix,
              cut = 0.05,
              FCcut = 2,
              outdir = NULL,
              control = NULL,
              treat = NULL
    ) {
        if (is.null(outdir)) {
            outdir = getwd()
        }
        if (!dir.exists(outdir)) {
            dir.create(outdir)
        }
        allgroups = unique(colnames(countMatrix))
        group = cbind(sample = allgroups, group = allgroups) %>% as.data.frame()

        while (!(control %in% allgroups)) {
            cat('可用样本名：', paste(allgroups), '\n')
            control <- readline(prompt = "输入对照样本名称: ")
        }
        allgroups = allgroups[!(allgroups %in% control)]
        if (length(allgroups) == 1) {
            treat = allgroups[1]
            cat('实验样本名称（自动选择）：', treat, '\n')
        } else{
            while (!(treat %in% allgroups)) {
                cat('可用样本名：', paste(allgroups), '\n')
                treat <- readline(prompt = "输入实验样本名称: ")
            }
        }
        meta.data = group %>% filter(group %in% c(control, treat))
        meta.data$group = as.factor(meta.data$group)
        countMatrix = countMatrix %>% select(meta.data$sample)
        countMatrix = as.matrix(countMatrix)
        Sizes = MedianNorm(countMatrix)
        EBOut = EBTest(
            countMatrix,
            Conditions = meta.data$group,
            sizeFactors = Sizes,
            maxround = 5
        )
        FCCut = 1 / FCcut
        EBDERes = GetDEResults(EBOut, FDR = cut, Threshold_FC = FCCut)
        FDR = 1 - EBOut$PPMat[, 'PPDE']
        GeneFC <- PostFC(EBOut)
        outFC = GeneFC$PostFC
        if (GeneFC$Direction == paste0(control, ' Over ', treat)) {
            outFC = 1 / outFC
        }
        out.remain = cbind(log2(outFC), FDR)

        nFilt = length(EBOut$AllZeroIndex)
        out.filt = cbind(rep(NA, nFilt), rep(NA, nFilt))
        row.names(out.filt) = names(EBOut$AllZeroIndex)
        out = rbind(out.remain, out.filt)
        colnames(out) = c('log2(FoldChange)', 'FDR')
        out = rownames_to_column(as.data.frame(out), 'ID')
        diff = out[match(EBDERes$DEfound, out$ID), ]

        write.table(
            out,
            paste0(outdir, '/', control, '.vs.', treat, '.AllGene.tsv'),
            row.names = F,
            quote = FALSE,
            sep = '\t'
        )
        write.table(
            diff,
            paste0(outdir, '/', control, '.vs.', treat, '.DEG.tsv'),
            row.names = F,
            quote = FALSE,
            sep = '\t'
        )
        write.table(
            diff$ID,
            paste0(outdir, '/', control, '.vs.', treat, '.DEGlist.tsv'),
            row.names = F,
            quote = FALSE,
            sep = '\t'
        )
        cat('差异基因分析结果见以下目录:\n')
        cat(paste0(outdir, '/', control, '.vs.', treat, '.AllGene.tsv'),
            '\n')
        cat(paste0(outdir, '/', control, '.vs.', treat, '.DEG.tsv'),
            '\n')
        cat(paste0(outdir, '/', control, '.vs.', treat, '.DEGlist.tsv'),
            '\n')
        return(diff)
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
                       showlabel.num = 10) {
    checkParams(showlabel,
                c("allSig", 'topPvalue', 'topFC', "no"),
                'showlabel')
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
        ggplot(input, aes(log2FoldChange, -log10(Pvalue))) + geom_point(aes(col =
                                                                                sig, shape = lim), size = 2)
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


checkParams = function(value, candidate, string) {
    if (!(value %in% candidate)) {
        stop(call. = F,
            '\n',
            string,
            '参数错误！\n可输入的值为:',
            toString(candidate),
            '\n而您的输入是:',
            value,
            '\n'
        )
    }
}

checkGroup = function(matrix, group) {
    samples = colnames(matrix)
    samples2 = group[, 1]
    if (all(samples != samples2)) {
        stop(
            '\nCount矩阵中定义的样本名与Group文件中定义的样本名不对应！\n矩阵样本名为:',
            toString(samples),
            '\n而Group文件的样本名为:',
            toString(samples2),
            '\n'
        )
    }
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
                  showlabel.num = 10) {
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
                                                                                       sig, shape = lim), size = 2) +
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
#' @title merge column values according group information.
#' @description  merge column values according group information.
#' @export
matrixGroup = function(matrix, groupInfo, method = "mean") {
    checkParams(method,
                c("mean", 'median', 'sum'),
                'method')
    checkGroup(matrix, groupInfo)
    matrix = t(as.data.frame(apply(count, 1, function(x)
        tapply(as.double(x), groupInfo[, 2], method))))
    return(matrix)
}
