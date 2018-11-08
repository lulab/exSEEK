if (tolower(Sys.getenv('CI')) != 'true') {
    if (file.exists('~/.Rprofile')) source('~/.Rprofile')

    options(knitr.package.root.dir = normalizePath('./'))

    options('defaultPackages' = unique(c(getOption('defaultPackages'), 'magrittr')));


    # setdiff2 <- function(x,y){setdiff(y,x)}
}
