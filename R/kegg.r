#' @title KEGG pathway enrichment tools.
#' @description  do KEGG pathway enrichment by a list of degs.
#' @export
KEGGenrich <- function(deglist,
                       taxon,
                       outdir = NULL,
                       outprefix,
                       useFDR = FALSE) {
    checkParams(taxon, names(keggdb), 'taxon')
    pathinfo = keggdb[[taxon]]
    res=KEGGenrich_common(deglist,pathinfo,outdir,outprefix,fdr=useFDR)
    return(res)
}

#' @title KEGG pathway enrichment tools for blastkoala ko annotation.
#' @description  do KEGG pathway enrichment by a list of degs. Using blastkoala method to convert from GeneID to KEGG ID.
#' @export
KEGGenrich_blastkoala <- function(deglist,
                       blastkoalafile,
                       taxonid='ko',
                       outdir = NULL,
                       outprefix,
                       useFDR = FALSE) {
    pathinfo = KAAS2Keggdb(blastkoalafile,taxonid)
    res=KEGGenrich_common(deglist,pathinfo,outdir,outprefix,fdr=useFDR)
    return(res)
}

#' @title KEGG pathway enrichment tools commmon.
#' @description  do KEGG pathway enrichment by a list of degs common.
#' @import R2HTML
#' @importFrom magrittr %>%
KEGGenrich_common=function(deglist,
                  pathinfo,
                  outdir = NULL,
                  outprefix,
                  fdr = FALSE){
    if (is.null(outdir)) {
        outdir = getwd()
    }
    if (!dir.exists(outdir)) {
        dir.create(outdir)
    }
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

#' @title convert KEGG kaas annotation to keggdb addon.
#' @description  convert KEGG kaas annotation to keggdb addon.
#' @importFrom magrittr %>%
#' @importFrom dplyr left_join
#' @importFrom dplyr select
#' @importFrom dplyr filter
#' @import KEGGREST
KAAS2Keggdb=function(KAASfile,taxonid='ko'){

    koquery=read.delim(KAASfile,sep="\t",header=F,na.strings = "",col.names = c('GeneID','ID'))%>%
        filter(!is.na(ID))

    koquery$ID=paste0("ko:",koquery$ID)
    pathwayName = keggList('pathway', organism = taxonid)
    pathwayName = cbind(
        PathwayID=names(pathwayName) %>% gsub('path:','',x = .),
        Pathway=unname(pathwayName)
    ) %>% as.data.frame()

    kopathway=keggLink("pathway",taxonid)
    kopathway=data.frame(ID=names(kopathway),PathwayID=unname(kopathway))
    kopathway=kopathway[grep(pattern = paste0('path:',taxonid),x = kopathway$PathwayID ),]
    kopathway$PathwayID=gsub('path:','',kopathway$PathwayID)

    if(taxonid!='ko'){
        gid2ko=keggLink('ko',taxonid)
        kopathway$ID=gid2ko[kopathway$ID]
    }


    res=left_join(kopathway,koquery,by="ID")%>%left_join(pathwayName,by="PathwayID")%>%select('GeneID','ID','PathwayID','Pathway')%>%filter(!is.na(GeneID))
    return(res)
}
