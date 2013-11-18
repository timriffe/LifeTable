
# Author: triffe
###############################################################################

library(devtools)
load_all("/home/triffe/git/LifeTable/LifeTable")
document("/home/triffe/git/LifeTable/LifeTable")

library(tools)
parent.path <- "/home/triffe/git/LifeTable/"
package.path <- file.path(parent.path,"LifeTable")
Rdfiles     <- list.files(file.path(package.path , "man"))
sapply(Rdfiles, function(xxx, parent.path, package.path){
            htmlname <- gsub("\\.Rd", "\\.html", xxx)
            Rd2HTML(file.path(package.path,"man", xxx),
                    out = file.path(parent.path,"help",htmlname),
                    stylesheet = "R.css")
            #go to folder
        },parent.path = parent.path, package.path = package.path)
system(paste0("cd ", parent.path, " \n git checkout gh-pages \n git add help * \n git commit -m 'update help'"))


