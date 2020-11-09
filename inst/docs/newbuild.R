# Modeled after roxygenizator.R
# Run this from the command line with Rscript with
# one trailing argument that will be the name of your
# new BTYD package, such as:
# Rscript newbuild.R BTYD

args <- commandArgs(trailingOnly = TRUE)
library("purrr")

# Whatever you named your version of the BTYD package
# (e.g. BTYD2) assign it to the `where` variable here:
where <- args[1]

# And `what` may be either `where` or `tarball`
# depending on whether you want to R CMD check the
# installed build or the tarball produced by
# R CMD build
what <- args[2]

# make sure the `vignettes` subdirectory is there
# and that it has a `figure` subdirectory of its own:
foo <- setwd(file.path(here::here(), where))
dir.create(file.path('vignettes', 'figure'))
devtools::document()
setwd(foo)
rm(foo)
foo <- setwd(file.path(here::here(), where, 'man'))

# a helper for getting the package version
getPackageVersion <- function(where) {
  x <- readLines(file.path(here::here(), 
                           where, 
                           'DESCRIPTION'))
  v <- 'Version: '
  gsub(pattern = v, 
       replacement = '', 
       x[grepl(v, x)])
}
# Obliterate the current source tarball
tarball <- paste(paste(where, 
                       getPackageVersion(where), 
                       sep = '_'), 
                 'tar.gz', 
                 sep = '.')
system(paste('rm', tarball))

setwd(foo)
rm(foo)

# Now rebuild, re-install. This will throw an
# error if you screwed anything up. 
sourcepath <- file.path(here::here(), where)
system(paste('rm', 
             file.path(sourcepath, 
                       '.DS_Store')))
system(paste('rm', 
             file.path(sourcepath, 
                       'inst', 
                       '.DS_Store')))
system(paste('rm', 
             file.path(sourcepath, 
                       'R', 
                       '.DS_Store')))
system(paste('rm', 
             file.path(sourcepath, 
                       'data', 
                       '.DS_Store')))
system(paste('rm', 
             file.path(sourcepath, 
                       'man', 
                       '.DS_Store')))
system(paste('rm', 
             file.path(sourcepath, 
                       'inst', 
                       'docs', 
                       '.DS_Store')))
system(paste('rm', 
             file.path(sourcepath, 
                       'vignettes', 
                       '.DS_Store')))
system(paste('rm -Rf', 
             file.path(sourcepath, 'inst', 'doc')))

message('')
message('Watch for build notes, warnings here:')
system(paste("R CMD build", where))
# install the newest version:
try(remove.packages(where))
install.packages(tarball, 
                 repos = NULL, 
                 type = 'source')

# to run code examples from the vignette
# you need to unzip the csv files first:
if(what == 'where') {
  system(paste('gunzip', file.path(.libPaths(), where, 'data', 'discreteSimElog.csv.gz')))
  system(paste('gunzip', file.path(.libPaths(), where, 'data', 'cdnowElog.csv.gz')))
}

# Now check the build; but first,
# copy new vignette goods where 
# R CMD check wants to see them:
dir.create(file.path(sourcepath, 'inst', 'doc'))
drop_log <- paste(where, 'walkthrough.log', sep = '-')
system(paste('rm', 
             file.path(sourcepath, 
                       'vignettes', drop_log)))
system(paste('cp -R',
             file.path(sourcepath,
                       'vignettes/'),
             file.path(sourcepath,
                       'inst', 'doc')))

# And drop the .Rcheck folder.
# `R CMD check` will write a fresh one:
system(paste('rm', file.path(here::here(), paste(where, 'Rcheck', sep = '.'))))

message('')
message('Watch for check notes, warnings here:')
checkit <- paste("R CMD check", get(what))
if(what == 'tarball') {
  checkit <- paste("R CMD check --no-manual --ignore-vignettes --as-cran", get(what))
}  
system(checkit)

# you need this to rebuild the help library, 
# but it only works if you're in RStudio already:
# .rs.restartR()
