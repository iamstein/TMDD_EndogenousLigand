#headers
  library(dplyr)
  library(readr)
  library(readxl)
  library(ggplot2)
  source("../../pgm_ModelG/xgx_packages_functions.R")

#molecular weights  
  MW_IL6   = 25;  #kDa
  MW_sIL6R = 80;  #kDa
  MW_toci  = 148; #kDa
  MW_silt  = 145; #kDa    
  MW_CRP   = 25.1 #kDa
  
#load in the datasets  
  col_names = c("TIME","VALUE")
  d.toci = read_csv("toci_ugml.csv",col_names=col_names) %>%
    mutate(VALUE    = VALUE*1e3/MW_toci,
           VALNAME  = "Tocilizumab",
           VALUNIT  = "nM")
  d.IL6  = read_csv("IL6_pgml.csv",col_names=col_names) %>%
    mutate(VALUE   = VALUE*1e-3/MW_IL6,
           VALNAME = "IL6",
           VALUNIT = "nM")
  d.sIL6R= read_csv("sIL6R_ngml.csv",col_names=col_names) %>%
    mutate(VALUE   = VALUE/MW_sIL6R,
           VALNAME = "sIL6R",
           VALUNIT = "nM")
  d.CRP  = read_csv("CRP_mgL.csv",col_names=col_names) %>%
    mutate(VALUE   = VALUE*1e3/MW_CRP,
           VALNAME = "CRP",
           VALUNIT = "nM")

#merge the datasets
  d = bind_rows(d.toci,d.IL6,d.sIL6R,d.CRP) %>%
    mutate(REF = "Charoin10",
           TIMEUNIT = "days",
           TIME = round(TIME)) %>%
    arrange(-VALUE) %>%
    mutate(VALNAME = factor(VALNAME,levels=unique(VALNAME))) %>%
    arrange(TIME,VALNAME) %>%
    select(TIME,TIMEUNIT,VALNAME,VALUE,VALUNIT,everything())
  
  write_csv(d,"./Task01_Charoin10_Data_1pat.csv")
  
#plot data ----
  g = ggplot(data=d,aes(x=TIME,y=VALUE,color=VALNAME,shape=VALNAME))
  g = g + geom_point()
  g = g + geom_line()
  g = g + scale.y.log10()
  g = g + scale_x_units(t.end = max(d$TIME))
  g = g + labs(y="Concentration (nM)",
               color = "Measurement",
               shape = "Measurement")
  g = g + ggtitle("Charoin10 Data from 1 patient")
  print(g)
  ggsave("./Task01_Charoin10_Data_1pat.png",width=5,height=5)
  
             
    