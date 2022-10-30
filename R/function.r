#' @title GOdbInfo.
#' @description  Show GO database information.
#' @importFrom magrittr %>%
#' @export
GOdbInfo = function() {
    info = sapply(names(godb), function(x)
        godb[[x]][['version']]) %>% as.data.frame()
    info2 = sapply(names(godb), function(x)
        godb[[x]][['update']]) %>% as.data.frame()
    demoid=sapply(names(keggdb), function(x) head(keggdb[[x]][['GeneID']],n=1))%>% as.data.frame()

    info = cbind(rownames(info), info, info2,demoid) %>% as.data.frame()
    colnames(info) = c('taxon', 'ReferenceGenomeVersion', 'update','IDformat')
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
    demoid=sapply(names(keggdb), function(x) head(keggdb[[x]][['GeneID']],n=1))%>% as.data.frame()

    info = cbind(rownames(info), info, info2,demoid) %>% as.data.frame()
    colnames(info) = c('taxon', 'ReferenceGenomeVersion', 'update','IDformat')
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

#' @param countMatrix  input count matrix, rows are genes and columns are samples
#'
#' @param group   group dataframe, first column: sample ID, second column: group name
#' @param useFDR  whether to use FDR to define significant DEG
#' @param cut  significant P value threshold, is useFDR is set , FDR is used
#' @param FCcut Fold Change threshold
#' @param outdir  output file directory
#' @param control control group name, used as denominator in Fold Change calculation
#' @param treat  treat group name, used as numerator  in Fold Change calculation
#'
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
              control = "Manually_select",
              treat = "Manually_select",
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
        exps=counts(dds, normalized=TRUE)
        out = filter(out, !is.na(padj))

        x=exps[match(rownames(out),rownames(exps)),match(meta.data$sample,colnames(exps))]
        colnames(out) = c('MeanExpression', 'log2FoldChange', 'Pvalue', 'FDR')
        out=cbind(x,out)
        out = rownames_to_column(out, 'ID')
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
#' @param countMatrix  input count matrix, rows are genes and columns are samples
#'
#' @param cut  significant threshold of FDR
#' @param FCcut Fold Change threshold
#' @param outdir  output file directory
#' @param control control sample name, used as denominator in Fold Change calculation
#' @param treat  treat sample name, used as numerator  in Fold Change calculation
#'
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
              control = "Manually_select",
              treat = "Manually_select"
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
        colnames(out) = c('log2FoldChange', 'FDR')
        normexp=EBOut$DataNorm[match(rownames(out),rownames(EBOut$DataNorm)),]
        MeanExpression=rowMeans(normexp)
        out=cbind(normexp,MeanExpression,out)
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
        filepath = NULL
        filepath$allgene = paste0(outdir, '/', control, '.vs.', treat, '.AllGene.tsv')
        filepath$deg = paste0(outdir, '/', control, '.vs.', treat, '.DEG.tsv')
        return(filepath)
    }

checkParams = function(value, candidate, string) {
    if (!(all(value %in% candidate))) {
        stop(call. = F,
            '\n',
            string,
            '参数错误！\n可输入的值为:',
            toString(candidate),
            '\n而您的输入是:',
            paste(value,collapse = ','),
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

#' @param matrix input dataframe
#'
#' @param groupInfo group dataframe, with two column : ID,Group
#' @param method  grouping method , "mean", 'median', 'sum'
#'
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

