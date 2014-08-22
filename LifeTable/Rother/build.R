
# Author: triffe
###############################################################################
install.packages("svcm")
library(devtools)
load_all("/home/tim/git/LifeTable/LifeTable")
document("/home/tim/git/LifeTable/LifeTable")

library(tools)
parent.path <- "/home/triffe/git/LifeTable"
package.path <- file.path(parent.path,"LifeTable")
Rdfiles     <- list.files(file.path(package.path , "man"))
sapply(Rdfiles, function(xxx, parent.path, package.path){
            htmlname <- gsub("\\.Rd", "\\.html", xxx)
            html.path<-file.path(parent.path,"help",htmlname)
            if (file.exists(html.path)){
                file.remove(html.path)
            }
            Rd2HTML(file.path(package.path,"man", xxx),
                    out = html.path,
                    stylesheet = "R.css")
            #go to folder
        },parent.path = parent.path, package.path = package.path)
system(paste0("cd ", parent.path, " \n git checkout gh-pages \n git add help * \n git commit -m 'update help'"))


# ---------------------------------------------
library(devtools)
install_github("LifeTable", subdir = "LifeTable", username = "timriffe",ref="master")
args(install_github)
library(LifeTable, lib="/home/triffe/R/x86_64-pc-linux-gnu-library/2.13")
data(UKRmales1965)
Nx <- UKRmales1965[, 3]
Dx <- UKRmales1965[, 2]
LT(Nx, Dx, ages = 0:110, axsmooth = TRUE)$ex[1]
source("/home/triffe/git/LifeTable/LifeTable/R/LT.R")

installed.packages()














