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
                       alllist=NULL,
                       taxon,
                       outdir = NULL,
                       outprefix,cutP=0.05,
                       useFDR = FALSE) {
    checkParams(taxon, names(keggdb), 'taxon')
    pathinfo = keggdb[[taxon]]
    if(!is.null(alllist)) pathinfo=pathinfo%>%filter(GeneID%in%alllist)
    res=KEGGenrich_common(deglist,pathinfo,outdir,outprefix,fdr=useFDR,cutP = cutP)
    return(res)
}

#' @param deglist  DEG vector
#'
#' @param tableFile  custom kegg to map file
#' @param outdir output directory
#' @param outprefix output file prefix
#' @param useFDR  whether to add FDR column
#'
#' @title KEGG pathway enrichment tools using custom annotation file.
#' @description  do KEGG pathway enrichment by a list of degs using custom annotation file.
#' @export
KEGGenrich_customTable <- function(deglist,alllist=NULL,
                       tableFile,
                       outdir = NULL,
                       outprefix,cutP=0.05,
                       useFDR = FALSE) {
    gid2keggpathway=read.delim(tableFile,header=T,sep="\t")
    colNum=ncol(gid2keggpathway)
    if(colNum==3){
        checkParams(colnames(gid2keggpathway),c('GeneID','KOID','PathwayID'),string = "tableFile列名错误")
        gid2keggpathway=unique(gid2keggpathway)
        dbpathwayInfo=keggList('pathway')
        gid2keggpathway$Pathway=dbpathwayInfo[gid2keggpathway$PathwayID]
    }else{
        checkParams(colnames(gid2keggpathway),c('GeneID','KOID','PathwayID','Pathway'),string = "tableFile列名错误")
    }
    if(!is.null(alllist)) gid2keggpathway=gid2keggpathway%>%filter(GeneID%in%alllist)
    res=KEGGenrich_common(deglist,gid2keggpathway,outdir,outprefix,fdr=useFDR,cutP)
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
KEGGenrich_blastkoala <- function(deglist,alllist=NULL,
                       blastkoalafile,
                       taxonid='ko',
                       outdir = NULL,
                       outprefix,
                       useFDR = FALSE) {
    pathinfo = KAAS2Keggdb(blastkoalafile,taxonid)
    if(!is.null(alllist)) pathinfo=pathinfo%>%filter(GeneID%in%alllist)
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
#' @importFrom dplyr group_by
#' @importFrom dplyr summarise
#' @importFrom dplyr left_join
#' @importFrom dplyr inner_join
#' @importFrom dplyr slice_head
#' @importFrom tidyr separate_rows
KEGGenrich_common=function(deglist,
                  pathinfo,
                  outdir = NULL,
                  outprefix,
                  cutP,
                  fdr = FALSE){
    if (is.null(outdir)) {
        outdir = getwd()
    }
    if (!dir.exists(outdir)) {
        dir.create(outdir)
    }
    pathway = unique(pathinfo[, c('PathwayID', 'Pathway')])

    keypref <- gsub("[0-9].*", "", pathway$PathwayID[1])
    if(keypref != 'ko' && keypref !='map'){
        temp_dir_path <- paste0(dirname(tempdir()),'/PlantNGSTools')
        if(!dir.exists(temp_dir_path)) dir.create(temp_dir_path)
        ko2idfile=paste0(temp_dir_path,"/",keypref,'.ko2id')
        if(file.exists(ko2idfile)){
            ko2id=read.csv(ko2idfile,header=T,col.names = c('ko','gid','path'))
        }else{
            ko2id=KEGGREST::keggLink(keypref,'ko')
            ko2id=data.frame(ko=names(ko2id),gid=ko2id)
            ko2id$ko=gsub('ko:','',ko2id$ko)
            gid2path=KEGGREST::keggLink('pathway',keypref)
            gid2path=data.frame(path=gsub("path:","",gid2path),gid=names(gid2path))
            ko2id=inner_join(ko2id,gid2path,by=c('gid'='gid'))
            write.csv(file = ko2idfile,ko2id,quote=F,row.names = F)
        }
    }else{
        ko2id=pathinfo%>%select(ko=KOID,gid=KOID,path=PathwayID)
    }
    ko2id=ko2id%>%group_by(ko,path)%>%slice_head(n=1)%>%ungroup()

    deggene = unique(deglist)
    allgene = unique(pathinfo$GeneID)
    deginfo = pathinfo[pathinfo$GeneID %in% deggene, ]%>%
        group_by(KOID,PathwayID,Pathway) %>%
        summarise(GeneID = paste(GeneID, collapse = ","), .groups = "drop")
    allinfo = pathinfo[pathinfo$GeneID %in% allgene, ]%>%
        group_by(KOID,PathwayID,Pathway) %>%
        summarise(GeneID = paste(GeneID, collapse = ","), .groups = "drop")
    k = length(unique(deginfo$KOID))
    N = length(unique(allinfo$KOID))
    p <- data.frame()
    urlColor <- NULL
    geneLink <- NULL
    allurlColor <- NULL
    allgeneLink <- NULL
    DEGdetailDat<-data.frame()
    ALLdetailDat<-data.frame()
    options(scipen = 9)
    get_first_two_elements <- function(string) {
        elements <- unlist(strsplit(string, ","))
        if (length(elements) > 2) {
            paste(elements[1:2], collapse = ",")
        } else {
            paste(elements, collapse = ",")
        }
    }

    showString=function(koid,x){
        first_elements <- sapply(x, function(string) {
            elements <- unlist(strsplit(string, ","))
            if (length(elements) > 1) {
                paste0(paste(get_first_two_elements(string), sep = ","),"...")
            } else {
                get_first_two_elements(string)
            }
        })
        res=paste0(koid,' (',first_elements,')')
        return(res)
    }


    for (i in seq(1, nrow(pathway))) {
        pid = pathway[i, 'PathwayID']
        allInPath = allinfo%>%filter(PathwayID==pid)
        m = length(unique(allInPath$KOID))
        degInPath=deginfo%>%filter(PathwayID==pid)
        DEGdetailDat=rbind(DEGdetailDat,degInPath%>%separate_rows(GeneID,sep=","))
        ALLdetailDat=rbind(ALLdetailDat,allInPath%>%separate_rows(GeneID,sep=","))
        x = length(unique(degInPath$KOID))

        degShowString=showString(degInPath$KOID,degInPath$GeneID)
        allShowString=showString(allInPath$KOID,allInPath$GeneID)

        degInPath=degInPath%>%left_join(ko2id,by=c('KOID'='ko','PathwayID'='path'))
        p[i, 1] = x
        p[i, 2] = m
        p[i, 3] = phyper(x - 1, m, N - m, k, lower.tail = FALSE)
        urlColor[i] = apply(as.matrix(paste(
            "/", t(unique(degInPath$gid)),  sep = ""
        )), 2, paste, collapse = "")
        Link = paste("<a target=\"_blank\" href=\"https://www.genome.jp/entry/",
                      degInPath$KOID,
                      "\"/>",
                     degShowString,
                     "</a>",sep=""
                     )

        if (x == 0) {
            geneLink[i] = "--"
        } else{
            geneLink[i] = apply(as.matrix(Link), 2, paste, collapse = "<br>")
        }

        geneLink[i]=paste('<button class="btn btn-link "  data-toggle="popover" data-html="true" data-trigger="click" title="DEG List" data-content=\'',geneLink[i],'\'> ',x,'</button>')

        aLink = paste("<a target=\"_blank\" href=\"https://www.genome.jp/entry/",
                     allInPath$KOID,
                     "\"/>",
                     allShowString,
                     "</a>",sep=""
        )
        if (m == 0) {
            allgeneLink[i] = "--"
        } else{
            allgeneLink[i] = apply(as.matrix(aLink), 2, paste, collapse = "<br>")
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
    detailDat=left_join(DEGdetailDat,output%>%dplyr::select(-Pathway),by=c('PathwayID'='ID'))%>%
        dplyr::select(PathwayID,Pathway,Pvalue,KOID,GeneID)%>%unique()

    write.table(
        detailDat,
        file = paste(outdir, "/", outprefix, ".DetailDEGList.tsv", sep = ""),
        quote = F,
        row.names = F,
        col.names = T,
        sep = "\t"
    )
    write.table(
        output,
        file = paste(outdir, "/", outprefix, ".all.tsv", sep = ""),
        quote = F,
        row.names = F,
        col.names = T,
        sep = "\t"
    )
    if(fdr){
        sigout=output%>%filter(FDR<cutP)
    }else{
        sigout=output%>%filter(Pvalue<cutP)
    }
    write.table(
        sigout,
        file = paste(outdir, "/", outprefix, ".significant.tsv", sep = ""),
        quote = F,
        row.names = F,
        col.names = T,
        sep = "\t"
    )
    if(keypref == 'map' || keypref== 'ko'){
        urlTitles <- gsub("^[a-z]+", "ko", output$ID)
    }else{
        urlTitles <-output$ID
    }

    url = paste(
        'https://www.genome.jp/kegg-bin/show_pathway?',
        urlTitles,
        urlColor,
        sep = ""
    )
    htmlOUT <-
        transform(output2,
                  ID = paste('<a  target="_blank" href = \"', url, '\" />', output$ID, '  </a>'))
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
        '<link rel="stylesheet" href="https://cdn.bootcss.com/twitter-bootstrap/4.3.1/css/bootstrap.min.css">
        <script src="https://cdn.bootcss.com/jquery/3.2.1/jquery.min.js"></script>
        <script src="https://cdn.bootcss.com/popper.js/1.15.0/umd/popper.min.js"></script>
        <script src="https://cdn.bootcss.com/twitter-bootstrap/4.3.1/js/bootstrap.min.js"></script>
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
