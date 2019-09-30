source("ams_initialize_script.R")
#rxSetIni0(FALSE)
source("SCIM_calculation.R")  
source("ivsc_2cmt_RR_V1.R")
dirs$Rscript.name = "Task07_Check_Repar.R"
dirs$output.prefix= str_extract(dirs$Rscript.name,"^Task\\d\\d\\w?_")

model = ivsc_2cmt_RR_KeqT0L0(target = TRUE)

#Drug list to loop through for finding file names
drugs_list = list("Pembro","VEGFR1","VEGFR2","Atezoli","Ramuc","Siltux","Tociliz") #ADD THIS LINE

dfs = data.frame() #ADD THIS LINE

# Create paths to data files for each drug.
param = NULL #ADD THIS LINE
i = 1 #ADD THIS LINE
#Get params for all the drugs
drugs = drugs_list[1] #ADD THIS LOOP
  
  filename = "parameters/ModelG_Pembro_Params_update_RR.xlsx"


# List of parameters of interest.
parameters = c("T0", "keDT","Kss_TL", "kon_TL", "Kss_DT","kon_DT")


# Dose time, frequency, compartment, nominal dose
tmax = 52*7 #days
tau  = 21   #days
compartment = 2
dose.nmol = 100*scale.mpk2nmol

param.list = list()
all_params <- data.frame() #ADD THIS LINE

i = 1

  # Load parameters.
  param.as.double =  read.param.file(filename) #ADD THIS LINE (CHANGED VARIABLE NAME TO PARAM)
  param.as.double["T0"] = 1
  param.as.double["L0"] = 1
  param.as.double["Kss_DT"] = 1 #with(df_param,(koff_DT+keDT)/kon_DT)
  param.as.double["Kss_TL"] = 1 #with(df_param,(koff_TL+keTL)/kon_TL)
  param.as.double = param.as.double[model$pin]
  
  df_param =  as.data.frame(t(param.as.double))
  
  
  
  param.list[[i]] = data.frame(t(param.as.double)) %>%
    mutate(drug = drugs_list[[i]]) %>%
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



# param.table = bind_rows(param.list) %>%
#   mutate(Kd_DT  = koff_DT/kon_DT,
#          Kss_DT = (koff_DT + keDT)/kon_DT,
#          Kd_TL  = koff_TL/kon_TL,
#          Kss_TL = (koff_TL + keTL)/kon_TL)
# #View(param.table)
# 
# drug = param.table$drug
# param.tablet = param.table %>%
#   dplyr::select(-drug) %>%
#   t() %>%
#   as.data.frame() %>%
#   setNames(drug)
# #View(param.tablet)


write.csv(all_params, file = "Task07_sensitivity_all drugs and params_100_updated 04_24.csv")

data.plot = all_params %>%
  dplyr::select(fold.change.param, SCIM_sim, SCIM_thy_keTL_negroot, SCIM_thy_keTL0, AFIR_thy, drug,param) %>%
  gather(key,value,-c(fold.change.param,drug,param))

g = ggplot(data.plot, aes(x=fold.change.param,y=value,color=key,linetype=key))
g = g + geom_line(size = 1, alpha = .5) 
g = g + facet_grid(drug ~ param,scales = "free_y", switch = "y") 
g = g + scale_x_log10() 
g = g + scale_y_log10()
g = g + scale_color_manual   (values = c(SCIM_sim       = "black",
                                         SCIM_thy_keTL0 = "blue",
                                         SCIM_thy_keTL_negroot = "green",
                                         AFIR_thy = "red"))
g = g + scale_linetype_manual(values = c(SCIM_sim = "solid",
                                         SCIM_thy_keTL0 = "dotted",
                                         SCIM_thy_keTL_negroot = "dashed",
                                         AFIR_thy = "solid"))

print(g)

# Compare simplified SCIM eqns. 
# SCIM_thy_keTL_negroot is the most complex i.e not simplified version of SCIM
# SCIM_sim is the SCIM from the simulation
# 26, 29, and 31 refer to the eqn numbers in the latex doc. for the simplified SCIMs
data.SCIMs = all_params %>%
  dplyr::select(fold.change.param, SCIM_sim, SCIM_thy_keTL_negroot, SCIM_thy_keTL_negroot26, SCIM_thy_keTL_negroot31, drug,param) %>%
  gather(key,value,-c(fold.change.param,drug,param))

g <- ggplot(data.SCIMs, aes(x=fold.change.param,y=value,color=key,linetype=key)) + 
  geom_line(size = 1, alpha = .6) +
  facet_grid(drug ~ param,scales = "free_y", switch = "y") + 
  scale_x_log10() + 
  scale_y_log10() + 
  scale_color_manual(values = c(SCIM_sim = "gray25",
                                SCIM_thy_keTL_negroot = "green3",
                                SCIM_thy_keTL_negroot26 = "dodgerblue4",
                                #SCIM_thy_keTL_negroot29 = "yellow",
                                SCIM_thy_keTL_negroot31 = "darkorange")) + 
  scale_linetype_manual(values = c(SCIM_sim = "solid",
                                   SCIM_thy_keTL_negroot = "dashed",
                                   SCIM_thy_keTL_negroot26 = "dotdash",
                                   #SCIM_thy_keTL_negroot29 = "dotted",
                                   SCIM_thy_keTL_negroot31 = "dashed"))

print(g)

# Compare simplified SCIMs and AFIR. 
data.SCIMs2 = all_params %>%
  dplyr::select(fold.change.param, SCIM_sim, SCIM_thy_keTL_negroot, SCIM_thy_keTL_negroot26, AFIR_thy, drug,param) %>%
  gather(key,value,-c(fold.change.param,drug,param))

g <- ggplot(data.SCIMs2, aes(x=fold.change.param,y=value,color=key,linetype=key)) + 
  geom_line(size = 1, alpha = .6) +
  facet_grid(drug ~ param,scales = "free_y", switch = "y") + 
  scale_x_log10() + 
  scale_y_log10() + 
  scale_color_manual   (values = c(SCIM_sim = "gray25",
                                   SCIM_thy_keTL_negroot = "green3",
                                   SCIM_thy_keTL_negroot26 = "dodgerblue4",
                                   AFIR_thy = "red")) + 
  scale_linetype_manual(values = c(SCIM_sim = "solid",
                                   SCIM_thy_keTL_negroot = "dashed",
                                   SCIM_thy_keTL_negroot26 = "dotdash",
                                   AFIR_thy = "solid"))

print(g)