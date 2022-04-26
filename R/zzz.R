#Package version checking
.onAttach <- function(libname = find.package("dispRity"), pkgname = "dispRity") {
    packageStartupMessage(paste0("       --- dispRity package ---\nThis is the CRAN release version (1.6.9) of the package.\nFor news, vignettes and future releases,\nvisit https://github.com/TGuillerme/dispRity\n"))
}