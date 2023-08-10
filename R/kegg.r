#' @param deglist  DEG vector
#'
#' @param taxon  taxon id , use KEGGdbInfo() to check available ids
#' @param outdir output directory
#' @param outprefix output file prefix
#' @param useFDR  whether to add FDR column
#'
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

#' @param deglist  DEG vector
#'
#' @param taxonid  KEGG pathway taxon id prefix , for example, osa is for Oryza Sativa.  available ids : https://www.genome.jp/kegg-bin/find_org_www?mode=abbr&obj=mode.map
#' @param outdir output directory
#' @param outprefix output file prefix
#' @param useFDR  whether to add FDR column
#' @param blastkoalafile blastkoala output file path,  Gene ID to KEGG ID pair by row
#'
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

#'
#' @title KEGG pathway enrichment tools commmon.
#' @description  do KEGG pathway enrichment by a list of degs common.
#' @import R2HTML
#' @importFrom magrittr %>%
#' @importFrom dplyr select
#' @importFrom dplyr mutate
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
    allurlColor <- NULL
    allgeneLink <- NULL
    options(scipen = 9)
    for (i in seq(1, nrow(pathway))) {
        pid = pathway[i, 'PathwayID']
        allInPath = unique(allinfo[allinfo[, 'PathwayID'] %in% pid, ])
        m = length(unique(allInPath[, 1]))
        degInPath = unique(deginfo[deginfo[, 'PathwayID'] %in% pid, ])
        x = length(unique(degInPath[, 1]))
        p[i, 1] = x
        p[i, 2] = m
        p[i, 3] = phyper(x - 1, m, N - m, k, lower.tail = FALSE)
        urlColor[i] = apply(as.matrix(paste(
            "/", t(unique(degInPath[, 'ID'])), "%09red", sep = ""
        )), 2, paste, collapse = "")

        Link = paste("<a target=\"_blank\" href=\"",
                     paste("https://www.genome.jp/entry/", degInPath[, "ID"], sep = ""),
                     "/",
                     ">",
                     degInPath[, "GeneID"],
                     "\"</a>",
                     sep = "")

        if (x == 0) {
            geneLink[i] = "--"
        } else{
            geneLink[i] = apply(as.matrix(Link), 2, paste, collapse = ", ")
        }

        geneLink[i]=paste('<button class="btn btn-link "  data-toggle="popover" data-html="true" data-trigger="click" title="DEG List" data-content=\'',geneLink[i],'\'> ',x,'</button>')

        aLink = paste("<a target=\"_blank\" href=\"",
                     paste("https://www.genome.jp/entry/", allInPath[, "ID"], sep = ""),
                     "/",
                     ">",
                     allInPath[, "GeneID"],
                     "\"</a>",
                     sep = "")
        if (m == 0) {
            allgeneLink[i] = "--"
        } else{
            allgeneLink[i] = apply(as.matrix(aLink), 2, paste, collapse = ", ")
        }
        allgeneLink[i]=paste('<button class="btn btn-link"  data-toggle="popover" data-trigger="click" data-html="true" title="All Gene List" data-content=\'',allgeneLink[i],'\'> ',m,'</button>')
    }
    output = cbind(pathway, p)
    colnames(output) = c('ID',
                         'Pathway',
                         'DEGsInPathway',
                         'GenesInPathway',
                         'Pvalue')
    if (fdr) {
        output = cbind(output, FDR = p.adjust(p[, 3], method = "BH"))
    }
    ind = order(output[, "Pvalue"])
    output = output[ind, ]
    urlColor = urlColor[ind]
    DEGs = geneLink[ind]
    allGenes = allgeneLink[ind]
    output2 = cbind(output, DEGs,allGenes)

    ind = output$DEGsInPathway > 0
    output = output[ind, ]
    urlColor = urlColor[ind]
    output2 = output2[ind, ]

   if (fdr) {
        output2 = output2%>%select(ID,Pathway,DEGsInPathway=DEGs,GenesInPathway=allGenes,Pvalue,FDR)
    }else{
        output2 = output2%>%select(ID,Pathway,DEGsInPathway=DEGs,GenesInPathway=allGenes,Pvalue)
    }

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
                  ID = paste('<a  target="_blank" href = \"', url, '\" />', urlTitles, '  </a>'))
    target <-
        HTMLInitFile(
            Title = "KEGG Pathway Ernichment Results",
            outdir = outdir,
            filename = paste0(outprefix),
            BackGroundColor = "#f8f8f8",
            useGrid = F,
            useLaTeX = F,
            HTMLframe = F
        )
 #   write( #insert global css
 #       "<style>table{border-collapse:collapse;width:1280px;}</style>",
 #       target,
 #       append = T
 #   )
    write(
        '<link rel="stylesheet" href="https://cdn.staticfile.org/twitter-bootstrap/4.3.1/css/bootstrap.min.css">
        <script src="https://cdn.staticfile.org/jquery/3.2.1/jquery.min.js"></script>
        <script src="https://cdn.staticfile.org/popper.js/1.15.0/umd/popper.min.js"></script>
        <script src="https://cdn.staticfile.org/twitter-bootstrap/4.3.1/js/bootstrap.min.js"></script>
        <script>$(document).ready(function(){$(\'[data-toggle="popover"]\').popover();});</script>
        <script>$(\'body\').on(\'click\',function (e) {$(\'[data-toggle="popover"]\').each(function () {if (!$(this).is(e.target) && $(this).has(e.target).length === 0 && $(\'.popover\').has(e.target).length === 0) {$(this).popover(\'hide\');}});});</script>
        ' ,
        target,
        append = T
    )
    HTML(
        htmlOUT,
        classtable='"table table-hover table-striped table-bordered"',
        classfirstline = '"thead-dark"',
        file = target,
        Border = 0,
        innerBorder = 0,
        digits=3,nsmall=3,
        row.names = FALSE,
        sortableDF =FALSE,
        decimal.mark="."
    )
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
#' @importFrom dplyr slice
#' @importFrom dplyr group_by
#' @importFrom dplyr ungroup
#' @import KEGGREST
KAAS2Keggdb=function(KAASfile,taxonid='ko'){

    koquery=read.delim(KAASfile,sep="\t",header=F,na.strings = "",col.names = c('GeneID','ID'))%>%
    filter(!is.na(ID))%>%mutate(KOID=paste0('ko:',ID),.keep='unused')

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
       kopathway$KOID=gid2ko[kopathway$ID]
    }else{
       kopathway$KOID=kopathway$ID
    }
    res=left_join(kopathway,koquery,by="KOID")%>%left_join(pathwayName,by="PathwayID")%>%select('GeneID','ID','KOID','PathwayID','Pathway')%>%
        filter(!is.na(GeneID))%>%group_by(GeneID,KOID,PathwayID)%>%slice(1)%>%ungroup()%>%as.data.frame()

    return(res)
}
