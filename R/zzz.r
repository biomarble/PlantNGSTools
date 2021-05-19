##' @importFrom utils packageDescription
##' @importFrom topGO groupGOTerms
.onAttach <- function(libname, pkgname="PlantNGSTools") {
    pkgVersion <- packageDescription(pkgname, fields="Version")
    author<-packageDescription(pkgname,fields="Maintainer")
    msg <- paste0(pkgname, " v", pkgVersion, " load successfully!\n",
		  "Author:",author,"\n"
#               ,   "More Packages: https://biomarble.github.io\n"
    )
    packageStartupMessage(msg)
    where <- match(paste("package:", pkgname, sep=""), search())
    suppressMessages(groupGOTerms(where))
}
