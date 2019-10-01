source("ams_initialize_script.R")
source("SCIM_calculation.R")  
source("ivsc_2cmt_RR_V1.R")
dirs$Rscript.name = "Task02b_Sensitivity_SCIM_dose.R"
dirs$output.prefix= str_extract(dirs$Rscript.name,"^Task\\d\\d\\w?_")

model = ivsc_2cmt_RR_v1()
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
dose.nmol = 10*scale.mpk2nmol

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

params.to.iterate = data.frame(dose = lseq(scale.mpk2nmol*0.1,scale.mpk2nmol*100,100))

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

# modify dfs for plotting to compare theory to simulation
  d.plot = dfs %>%
    gather(metric_full,value,-c(param.to.change,fold.change.param,time_idx,param.to.change1,fold.change.param1,param)) %>%
    mutate(type      = str_replace(metric_full,"^.*\\.",""),
           metric    = str_replace(metric_full,"\\..*$",""))

#plot theory vs simulation
  g = ggplot(d.plot,aes(x=param.to.change,y=value,color=type,linetype=type))
  g = g + geom_line()
  g = g + facet_wrap(~metric)
  g = g + scale.x.log10()
  g = g + scale.y.log10()
  g = g + labs(x="Dose",
               y="Value",
               caption = "Andy isn't sure what D and Dss are.
                          (average value at steady state?  trough?)  
                          Name should be changed to be more descriptive 
                          and to be the same for both sim and thy for easy plotting.
                          Also, TLss should be named either .sim or .thy")
  gg = ggsaveplot(width=6,height=6,dirs,"SCIM_DoseSensitivity_Accuracy")
  grid.arrange(gg)





