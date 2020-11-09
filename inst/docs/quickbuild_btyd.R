# Modeled after newbuild.R
# Run this from the command line with Rscript:
# Rscript quickbuild_btyd.R
library("purrr")
library("devtools")

setwd(file.path(here::here(), 'BTYD'))

document()
build()
install()

# some cleanup
sourcepath <- file.path(here::here())
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
system(paste('rm', 
             file.path(sourcepath, 
                       'vignettes', 'BTYD-walkthrough.log')))

message('Restoring inst/doc because check() expects it')
dir.create(file.path('inst', 'doc'))
system(paste('cp -R',
             'vignettes/',
             file.path('inst', 'doc')))

check()
