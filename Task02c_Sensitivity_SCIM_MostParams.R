source("ams_initialize_script.R")
source("SCIM_calculation.R")  
source("ivsc_2cmt_RR_V1.R")
dirs$rscript_name = "Task02c_Sensitivity_SCIM_MostParams.R"
dirs$filename_prefix =  str_extract(dirs$rscript_name,"^Task\\d\\d\\w?_")

model = ivsc_2cmt_RR_v1(target = TRUE)

#Drug list to loop through for finding file names
drugs_list = list("Pembro","VEGFR1","VEGFR2","Atezoli","Ramuc","Siltux","Tociliz") #ADD THIS LINE

dfs = data.frame() #ADD THIS LINE

# Create paths to data files for each drug.
param = NULL #ADD THIS LINE
i = 1 #ADD THIS LINE
#Get params for all the drugs
for (drugs in drugs_list) { #ADD THIS LOOP
  
  filename = list.files(path = "parameters/",pattern = drugs) #change filename line to this
  if (length(filename)>1)
    stop("check and see if you have any temp files open or something.  maybe close excel")
  
  param[[i]] = paste0("parameters/",filename)
  i = i+1
}

# List of parameters of interest.
parameters = c("ksynT",  "keDT", "koff_TL", "kon_TL", "koff_DT","kon_DT")
#parameters = c("ksynT", "ksynL", "keT", "keL", "keDT", "keTL", "koff_TL", "kon_TL", "koff_DT","kon_DT")
#"keL","keTL" "ksynL", "keT", - these causes issues

# Dose time, frequency, compartment, nominal dose
tmax = 52*7 #days
tau  = 21   #days
compartment = 2
dose.mpk  = 100
dose.nmol = dose.mpk*scale.mpk2nmol
isSoluble = FALSE

param.list = list()
all_params <- data.frame() #ADD THIS LINE
for (i in 1:length(drugs_list)){ #loop over all the drugs in the list
  

  # Load parameters.
  param.as.double =  read.param.file(param[i]) #ADD THIS LINE (CHANGED VARIABLE NAME TO PARAM)
  df_param =  as.data.frame(t(param.as.double))
  
  param.list[[i]] = data.frame(t(param.as.double)) %>%
    mutate(drug = drugs_list[[i]]) %>%
    dplyr::select(drug,everything())
  
  # Set range for parameters of interest in SA.
  # Check which parameters are nonzero, not including dose which isn't in df_param.
  nnzero = df_param[parameters[which(parameters != "dose")]] != 0
  nnzero = colnames(nnzero)[which(nnzero)]
  params.to.iterate = data.frame(lapply(df_param[nnzero], function(x) lseq(x*1e-4, x*1e4, 15)))
  
  dfs = list() #Reset the temp list for every drug
  temp_dfs <- data.frame() #Reset the temporary dataframe
  
  # Iterate all of the parameters for a single drug.
  for (j in 1:ncol(params.to.iterate)){
    dfs[[j]] = compare.thy.sim(model = model,
                               param.as.double = param.as.double,
                               dose.nmol = dose.nmol,
                               tmax = tmax,
                               tau  = tau,
                               compartment = compartment,
                               param.to.change = names(params.to.iterate)[j],
                               param.to.change.range = params.to.iterate[[j]],
                               soluble = isSoluble)
    #REMOVE THIS LINE
    #dfs[[j]] = dfs[[j]] %>% mutate(drug = drugs_list[i], isSol = isSoluble)
  }
  #ADD THESE LINES
  temp_dfs <- bind_rows(dfs) #create a temp dataframe for all the data
  temp_dfs$drug <- as.character(drugs_list[i])
  all_params <- rbind(all_params,temp_dfs) #Cat data frame
}


param.table = bind_rows(param.list) %>%
  mutate(Kd_DT  = koff_DT/kon_DT,
         Kss_DT = (koff_DT + keDT)/kon_DT,
         Kd_TL  = koff_TL/kon_TL,
         Kss_TL = (koff_TL + keTL)/kon_TL)
#View(param.table)

drug = param.table$drug
param.tablet = param.table %>%
  dplyr::select(-drug) %>%
  t() %>%
  as.data.frame() %>%
  setNames(drug)
#View(param.tablet)



write.csv(all_params, file = "task02a_sensitivity_all drugs and params_100_updated 04_24.csv")

data.plot = all_params %>%
  dplyr::select(fold.change.param, SCIM_sim, SCIM_thy_keTL_negroot, SCIM_thy_keTL0, AFIR_thy, drug,param) %>%
  gather(key,value,-c(fold.change.param,drug,param))

g <- ggplot(data.plot, aes(x=fold.change.param,y=value,color=key,linetype=key)) + 
  geom_line(size = 1, alpha = .5) +
  facet_grid(drug ~ param,scales = "free_y", switch = "y") + 
  xgx_scale_x_log10() + 
  xgx_scale_y_log10() + 
  scale_color_manual(values = c(SCIM_sim       = "black",
                                SCIM_thy_keTL0 = "blue",
                                SCIM_thy_keTL_negroot = "green",
                                AFIR_thy = "red")) + 
  scale_linetype_manual(values = c(SCIM_sim = "solid",
                                   SCIM_thy_keTL0 = "dotted",
                                   SCIM_thy_keTL_negroot = "dashed",
                                   AFIR_thy = "solid")) + 
  theme(legend.position="top") + 
  ggtitle(paste0("Dose = ", dose.mpk, "mg/kg every ", tau/7, "weeks")) +
  xlab("Fold Change in Parameter")

g = xgx_save(11,11,dirs,paste0("MostSensitivityAnalysis_",dose.mpk,"mpk"),draft.flag)
print(g)

