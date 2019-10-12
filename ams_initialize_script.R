#clear variables and plots
rm(list=ls()) #clear all variables

#user name  
user = Sys.getenv("USER")

#set seed
set.seed(123456)  

## Ensure that no library from the user home will be loaded such that
## we run with the production packages only
.libPaths(grep("home", .libPaths(), value=TRUE, invert=TRUE))    
.libPaths(grep("usr", .libPaths(), value=TRUE, invert=TRUE))

#add key packages
#library(GGally)
#library(MASS)  
#library(rmarkdown)
library(knitr)
#library(foreign)
#library(caTools)
#library(grid)
#library(gridExtra)
#library(reshape2)
library(stringr)
library(tidyr)
library(dplyr)
library(ggplot2)
library(readr)
library(readxl)
#library(png)

info = Sys.info()
if (info["sysname"] == "Linux" & info["login"] == "steinanf") {
  library(RxODE, lib = "./")
  library(xgxr, lib = "./")
} else {
  library(RxODE)
  library(xgxr)
}


#source additional key functions

#helper fuction for reading parameter files
#it assumes model has been defined above and only model$pin params are kept in the vector
read.param.file = function(filename) {
  d                      = read_excel(filename, 1)
  param.as.double        = suppressWarnings(as.numeric(d$Value))
  names(param.as.double) = d$Parameter
  #param.as.double        = param.as.double[model$pin] #keep only parameters used in ODE
  return(param.as.double)
}

lseq = function(from, to, length.out){
  # Arguments:
  #   from : initial value of the variable
  #   to : teminal value of the variable
  #   length.out : fold number of <to - from>
  # Return :
  #   A vector of length <length.out>
  
  sequence = seq(log(from), log(to), length.out=length.out)
  sequence = exp(sequence)
  return(sequence)
}

drugs = c("Atezolizumab",
          "Bevacizumab_VEGFR1",
          "Bevacizumab_VEGFR2",
          "Pembrolizumab",
          "Ramucirumab",
          "Siltuximab",
          "Tocilizumab")
parameter_files = paste0("parameters/ModelG_",drugs,"_Params.xlsx")
names(parameter_files) = drugs

#flag for labeling figures as draft
draft.flag           = FALSE
print.filenames.flag = TRUE

#set directories
top = paste0(normalizePath("./"),"/")
dirs= list(top.level         = top,
           Rscript.relative  = "./",
           results.relative  = "./results/",
           results.relative.print = "./results/",
           parent_dir = top,
           rscript_dir = "./",
           results_dir = "./results/")

#ggplot settings
theme_set(theme_bw())

#scaling params
scale.mpk2nmol = 70*1e-3/150e3*1e9
scale_mpk2nmol = scale.mpk2nmol
scale.nmol2mpk = 1/scale.mpk2nmol #nM->mg/kg
scale_nmol2mpk = scale.nmol2mpk
scale.mg2nmol  = 1e-3/150e3*1e9