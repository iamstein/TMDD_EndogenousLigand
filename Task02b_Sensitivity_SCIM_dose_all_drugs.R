source("ams_initialize_script.R")
source("SCIM_calculation.R")  
source("ivsc_2cmt_RR_V1.R")
dirs$rscript_name = "Task02b_Sensitivity_SCIM_dose_all_drugs.R"
dirs$filename_prefix =  str_extract(dirs$rscript_name,"^Task\\d\\d\\w?_")

model = ivsc_2cmt_RR_v1(target = TRUE)
parameters = c("dose")

#Drug list to loop through for finding file names
drugs_list = list("Pembro","VEGFR1","VEGFR2","Atezoli","Ramuc","Siltux","Tociliz") #ADD THIS LINE

dfs = data.frame() #ADD THIS LINE

# Create paths to data files for each drug.
param = NULL #ADD THIS LINE
i = 1 #ADD THIS LINE
#Get params for all the drugs
for (drugs in drugs_list) { #ADD THIS LOOP
  
  filename        = list.files("parameters/",pattern = drugs) #change filename line to this
  param[[i]] = paste0("parameters/",filename)
  i = i+1
}

# Dose time, frequency, compartment, nominal dose
tmax = 364 #days
tau  = 21   #days
compartment = 2
dose.nmol = 1*scale.mpk2nmol
isSoluble = FALSE

# --------------------------------------------------------------------------------
# Iterate over all the drugs.
# --------------------------------------------------------------------------------


all_params <- data.frame() #ADD THIS LINE
for (i in 1:length(drugs_list)){ #loop over all the drugs in the list
# Load parameters.
param.as.double =  read.param.file(param[i])
df_param =  as.data.frame(t(param.as.double))

# Set range for parameters of interest in SA.
# Check which parameters are nonzero, not including dose which isn't in df_param.

param.as.double["dose"] = dose.nmol
params.to.iterate = data.frame(dose = dose.nmol*lseq(0.01,1000,21))

# iterate over doses
  dfs = compare.thy.sim(model           = model,
                        param.as.double = param.as.double,
                        dose.nmol       = dose.nmol,
                        tmax            = tmax,
                        tau             = tau,
                        compartment     = compartment,
                        param.to.change = "dose",
                        param.to.change.range = params.to.iterate$dose,
                        soluble         = FALSE)

  #ADD THESE LINES
  temp_dfs <- bind_rows(dfs) #create a temp dataframe for all the data
  temp_dfs$drug <- as.character(drugs_list[i])
  all_params <- rbind(all_params,temp_dfs) #Cat data frame
  }
write.csv(all_params, "Task02b_sens_SCIM_dose_all_drugs.csv")

# plot results ----
data.plot = all_params %>%
  dplyr::select(param.to.change, SCIM_sim, SCIM_thy_keTL_negroot, SCIM_thy_keTL0, AFIR_thy, drug,param) %>%
  gather(key,value,-c(param.to.change,drug,param))

g <- ggplot(data.plot, aes(x=param.to.change/dose.nmol,y=value,color=key,linetype=key)) + 
  geom_line(size = 1, alpha = .5) +
  facet_wrap(~drug) + 
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
  xlab("Dose (mg/kg) every three weeks")

g = xgx_save(6,6,dirs,"DoseSensitivityAnalysis",draft.flag)
print(g)


