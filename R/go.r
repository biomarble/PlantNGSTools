#' @title GO enrichment .
#' @description  do GO enrichment by a list of DEGs.
#' @export
GOEnrich = function(deglist,
                    taxon,
                    fdr = FALSE,
                    outdir = NULL,
                    outprefix) {
    checkParams(taxon, names(godb), 'taxon')

    geneID2GO = godb[[taxon]][['db']]
    version = godb[[taxon]][['version']]
    cat('genome version:', version, '\n')
    res=GOenrich_common(deglist,geneID2GO,fdr=fdr,outdir=outdir,outprefix=outprefix)
    return(res)
}

#' @title GO enrichment using pannzer2 result.
#' @description  do GO enrichment by a list of DEGs using pannzer2 result.
#' @export
GOEnrich_pannzer2 = function(deglist,
                    pannzerfile,
                    fdr = FALSE,
                    outdir = NULL,
                    outprefix=NULL) {
    geneID2GO = PANNZERres2GOdb(pannzerfile)
    res=GOenrich_common(deglist,geneID2GO,fdr=fdr,outdir=outdir,outprefix=outprefix)
    return(res)
}

#' @title GO enrichment using custom table  annotation.
#' @description  do GO enrichment by a list of DEGs using custom table annotation.
#' @export
GOEnrich_customTable = function(deglist,
                           tableFile,
                           fdr = FALSE,
                           outdir = NULL,
                           outprefix=NULL) {
    geneID2GOtable=read.table(tableFile,header=T,sep="\t")
    checkParams(colnames(geneID2GOtable),c('GeneID','GOID'),string = "geneID2GOtable列名错误")

    geneID2GO=lapply(unique(geneID2GOtable$GeneID), function (x) geneID2GOtable$GOID[geneID2GOtable$GeneID==x])
    names(geneID2GO)=unique(geneID2GOtable$GeneID)

    res=GOenrich_common(deglist,geneID2GO,fdr=fdr,outdir=outdir,outprefix=outprefix)
    return(res)
}

#' @title GO enrichment using custom mapping file.
#' @description  do GO enrichment by a list of DEGs using custom mapping  file.
#' @importFrom topGO readMappings
#' @export
GOEnrich_customMapping = function(deglist,
                           mappingfile,
                           fdr = FALSE,
                           outdir = NULL,
                           outprefix=NULL) {
    geneID2GO<-readMappings(file=mappingfile)
    res=GOenrich_common(deglist,geneID2GO,fdr=fdr,outdir=outdir,outprefix=outprefix)
    return(res)
}

#' @title GO enrichment  sub.
#' @description  do GO enrichment by a list of DEGs.
#' @import ggplot2
#' @import GO.db
#' @importFrom  dplyr select
#' @import topGO
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

#' @title convert PANNZER GO annotation to godb addon.
#' @description  convert PANNZER GO annotation(GO.out file) to geneID2GO format.
#' @importFrom magrittr %>%
PANNZERres2GOdb=function(PANNZERfile){
    dat=read.table(PANNZERfile,header=T,sep="\t",)
    GOcontent=cbind(GID=dat$qpid,GOID=paste0('GO:', sprintf("%07d", dat$goid)))%>%as.data.frame()
    geneID2GO=lapply(unique(GOcontent$GID), function (x) GOcontent$GOID[GOcontent$GID==x])
    names(geneID2GO)=unique(GOcontent$GID)
    return(geneID2GO)
}

#' @title GO enrichment using  annotation.
#' @description  do GO enrichment by a list of DEGs using  annotation.
#' @import ggplot2
#' @import GO.db
#' @importFrom  dplyr select
#' @import topGO
GOenrich_common = function(deglist,
                           geneID2GO,
                           fdr = FALSE,
                           outdir = NULL,
                           outprefix=NULL) {
    if (is.null(outdir)) {
        outdir = getwd()
    }
    if (!dir.exists(outdir)) {
        dir.create(outdir)
    }
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
    all=all[!is.na(all$Pvalue),]
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
