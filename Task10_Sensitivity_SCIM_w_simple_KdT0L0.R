source("ams_initialize_script.R")
#rxSetIni0(FALSE)
source("SCIM_calculation.R")  
source("ivsc_2cmt_RR_V1.R")
dirs$Rscript.name = "Task10_Sensitivity_SCIM_w_simple_KssT0L0.R"
dirs$output.prefix= str_extract(dirs$Rscript.name,"^Task\\d\\d\\w?_")

model = ivsc_2cmt_RR_KeqT0L0()

# List of parameters of interest.
parameters = c("ksynT", "keDT","koff_TL", "kon_TL", "koff_DT","kon_DT")

# Dose time, frequency, compartment, nominal dose
tmax = 52*7 #days
tau  = 21   #days
compartment = 2
dose.nmol = 100*scale.mpk2nmol

param.list = list()
all_params <- data.frame() #ADD THIS LINE
for (i in 1:length(drugs)){ #loop over all the drugs in the list
  
  drug = drugs[i]
  # Load parameters.
  param.as.double =  read.param.file(parameter_files[i]) #ADD THIS LINE (CHANGED VARIABLE NAME TO PARAM)
  df_param =  as.data.frame(t(param.as.double))
  
  param.list[[i]] = data.frame(t(param.as.double)) %>%
    mutate(drug = drug) %>%
    dplyr::select(drug,everything())
  
  # Set range for parameters of interest in SA.
  # Check which parameters are nonzero, not including dose which isn't in df_param.
  nnzero = df_param[parameters[which(parameters != "dose")]] != 0
  nnzero = colnames(nnzero)[which(nnzero)]
  params.to.iterate = data.frame(lapply(df_param[nnzero], function(x) lseq(x*0.000001, x*1000000, 13)))
  
  dfs = list() #Reset the temp list for every drug
  temp_dfs <- data.frame() #Reset the temporary dataframe
  
  # Iterate all of the parameters for a single drug.
  for (j in 1:ncol(params.to.iterate)){
    param.to.change       = names(params.to.iterate)[j]
    param.to.change.range = params.to.iterate[[j]]
    dfs[[j]] = compare.thy.sim(model = model,
                               param.as.double = param.as.double,
                               dose.nmol = dose.nmol,
                               tmax = tmax,
                               tau  = tau,
                               compartment = compartment,
                               param.to.change = param.to.change,
                               param.to.change.range = param.to.change.range)
  }
  #ADD THESE LINES
  temp_dfs <- bind_rows(dfs) #create a temp dataframe for all the data
  temp_dfs$drug <- as.character(drugs_list[i])
  all_params <- rbind(all_params,temp_dfs) #Cat data frame
}


write.csv(all_params, file = "Task09_Sensitivity_SCIM_w_simple.csv")

#plot results ----
data.plot = all_params %>%
  dplyr::select(fold.change.param, drug, param, 
                SCIM_sim, SCIM_thy, AFIR_thy) %>%
  gather(key,value,-c(fold.change.param,drug,param)) %>%
  mutate(AFIR_SCIM  = ifelse(str_detect(key,"AFIR"),"AFIR","SCIM"),
         theory_sim = ifelse(str_detect(key,"sim"),"sim","thy"),
         approx     = str_extract(key,"_\\w+_"),
         approx     = str_replace(approx,"_",""),
         approx     = ifelse(is.na(approx),"none",approx))

g = ggplot(data.plot, aes(x=fold.change.param,y=value, color = key, linetype = key))
                          #,size=AFIR_SCIM, linetype = theory_sim, color = approx))
g = g + geom_line(size = 1, alpha = .5) 
g = g + facet_grid(drug ~ param,scales = "free_y", switch = "y") 
g = g + scale_x_log10() 
g = g + scale_y_log10()
#g = g + scale_size_manual(values = c(AFIR = 1, SCIM = 2))
# g = g + scale_color_manual   (values = c(SCIM_sim       = "black",
#                                          SCIM_thy_keTL0 = "blue",
#                                          SCIM_thy_keTL_negroot = "green",
#                                          AFIR_thy = "red"))
# g = g + scale_linetype_manual(values = c(SCIM_sim = "solid",
#                                          SCIM_thy_keTL0 = "dotted",
#                                          SCIM_thy_keTL_negroot = "dashed",
#                                          AFIR_thy = "solid"))

print(g)

#check the initial condition and steady state ----
check = all_params %>%
  dplyr::select(fold.change.param, drug, param, TLss_frac_change, TL0_05tau_frac_change) %>%
  gather(key,value,-c(fold.change.param,drug,param))

g = ggplot(check,aes(value))
g = g + geom_histogram()
g = g + facet_wrap(~key, scales = "free_x")
print(g)
