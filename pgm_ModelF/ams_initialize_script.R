#clear variables and plots
  rm(list=ls()) #clear all variables

#user name
  user = Sys.getenv("USER")

#add key packages
  library(ggplot2)
  #library(MASS)
  library(rmarkdown)
  library(knitr)
  library(foreign)
  library(caTools)
  library(grid)
  library(gridExtra)
  library(reshape2)
  library(stringr)
  #library(GGally)
  library(tidyr)
  library(dplyr)
  library(ggrepel)
  library(readxl)
  library(RxODE)

  mapvalues = plyr::mapvalues #only need this one function from plyr

#source additional key functions
  source("ams_graphics_v2.R")
  source("ams_tmdd_helper.R")
  source("AFIRT_calculation.R")
  source("ivsc_4cmtct_shedct.R")

#flag for labeling figures as draft
  draft.flag      = FALSE

#set directories
  top = paste0(normalizePath("../"),"/")
  dirs= list(top =top,
             topz=str_replace(paste0(normalizePath("../"),".zip"),
                              paste0("/view/",user,"_view"),""),
             pgm= paste0(top,"pgm/"),
             data.out   = paste0(top,"data/"),
             results    = paste0(top,"results/"),
             mlxruns    = paste0(top,"mlxruns/"),
             mlxtran    = paste0(top,"mlxtran/"),
             pgm.local  = "pgm/",
             results.local = "results/")

#ggplot settings
  #theme_set(theme_classic() + theme(panel.border=element_rect(color="black",fill=NA)))
  theme_set(theme_bw())
  
#scaling params
  scale.mpk2nmol = 70*1e-3/150e3*1e9
  scale.nmol2mpk = 1/scale.mpk2nmol #nM->mg/kg
  scale.mg2nmol  = 1e-3/150e3*1e9
