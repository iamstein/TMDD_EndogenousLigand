source("ams_initialize_script.R")
source("SCIM_calculation_2019_05_10.R")  
source("ivsc_2cmt_RR_V1.R")
dirs$Rscript.name = "Task02a_Sensitivity_SCIM_k_pars.R"
dirs$output.prefix= str_extract(dirs$Rscript.name,"^Task\\d\\d\\w?_")

model = ivsc_2cmt_RR_v1(target = TRUE)

#Drug list to loop through for finding file names
drugs_list = list("VEGFR2") 

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
parameters = c("koff_TL", "kon_TL")


# Dose time, frequency, compartment, nominal dose
tmax = 52*7 #days
tau  = 21   #days
compartment = 2
dose.nmol = 100*scale.mpk2nmol
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
  params.to.iterate = data.frame(lapply(df_param[nnzero], function(x) lseq(x*0.0001, x*1000, 13)))
  
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


write.csv(all_params, file = "task04.csv")

data.plot <- all_params %>%
  dplyr::select(fold.change.param, 
         TL0_sim,
         TL0_keTL0_thy,
         TL0_negroot_thy,
         TL0_posroot_thy,
         drug, param) %>%
  gather(key, value, -c(fold.change.param,drug,param))

g <- ggplot(data.plot, aes(x=fold.change.param, y=value, color=key, linetype=key)) + 
  facet_grid(drug ~ param, scales = "free_y", switch = "y") + 
  geom_line(size=1, alpha=0.6) +
  scale_x_log10() + 
  scale_y_log10()
  #scale_color_manual(values = c(#SCIM_thy_ketl_pos="orange",
  #  TL0.sim="black",
  #  TL0.thy="red",
  #  TL0_neg.thy="green",
  #  TL0_pos.thy="blue"))
print(g)



#compute the simulation with the final parameter set, just to make sure it's getting to steady state
ev = eventTable(amount.units="nmol", time.units="days")
sample.points = c(seq(-7, tmax, 0.1), 10^(-3:0)) # sample time, increment by 0.1
sample.points = sort(sample.points)
sample.points = unique(sample.points)
ev$add.sampling(sample.points)
#ev$add.dosing(dose=dose.nmol, nbr.doses=floor(tmax/tau)+1, dosing.interval=tau,
#              dosing.to=compartment)

init = model$init(param.as.double)
out  = model$rxode$solve(param.as.double, ev, init)
out  = model$rxout(out)

out.long = out %>%
  dplyr::select(time,T,L,TL) %>%
  gather(key,value,-c(time))
  
pars = as.data.frame(t(param.as.double))
a = with(pars,keTL^2)
b = with(pars,-(keTL) * (ksynT +ksynL) - (((koff_TL+keTL)/kon_TL) * keT *keL))
c = with(pars, ksynL*ksynT)

TL0_pos <- ((-b) + sqrt((b^2)-4*a*c))/(2*a)

TL0_neg <- ((-b) -sqrt((b^2)-4*a*c))/(2*a)



g = ggplot(out.long,aes(x=time,y=value,color=key))
g = g + geom_line()
g = g + scale_y_log10()
g = g + geom_hline(yintercept = TL0_neg,color="purple",linetype="dashed")
print(g)

