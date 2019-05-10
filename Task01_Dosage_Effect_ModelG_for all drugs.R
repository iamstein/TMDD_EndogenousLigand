source("ams_initialize_script.R")
source("ivsc_2cmt_RR_V1.R")
source("SCIM_calculation_2019_05_10.R")
dirs$Rscript.name = "Task01_Doseage_Effect_ModelG.R"
dirs$output.prefix= str_extract(dirs$Rscript.name,"^Task\\d\\d\\w?_")

model  = ivsc_2cmt_RR_v1()

#drugs_list should contain a partial drug name from the parameter files
drugs_list = list("Pembro","VEGFR1","VEGFR2","Atezoli","Ramuc","Siltux","Tociliz") #ADD THIS LINE


dfs = data.frame() #ADD THIS LINE

#Loop through all the drugs
for (drugs in drugs_list) { #ADD THIS LINE
  
  filename        = paste0("parameters/",
                           list.files("parameters/",pattern = drugs)) #change filename line to this
  param           = read.param.file(filename)
  gedf_param        = as.data.frame(t(param))
  
  d = read_excel(filename) %>%
    select(-Order) %>%
    mutate(Value = signif(Value,2))
  kable(d)
  
  #Set dosing parameters
  tmax            = 52*7 #days (maximum time point)
  tau             = 21  #days between doses
  compartment     = 2    #compartment drug is dosed into
  param.as.double =  param
  
  #Make temp dataframe
  temp_df = data.frame() #ADD THIS LINE
  
  #set up dosing
  dose.i = lseq(scale.mpk2nmol*0.01,scale.mpk2nmol*100,10)
  for(dose.nmol in dose.i) {
    
    #set sampling and dosing event (ev) table
    ev = eventTable(amount.units="nmol", time.units="days")
    sample.points = c(seq(-7, tmax, 0.1), 10^(-3:0)) # sample time, increment by 0.1
    sample.points = sort(sample.points)
    sample.points = unique(sample.points)
    ev$add.sampling(sample.points)
    ev$add.dosing(dose            = dose.nmol, 
                  nbr.doses       = floor(tmax/tau)+1, 
                  dosing.interval = tau,
                  dosing.to       = compartment)
    
    #model = do.call(model_name, list()) # model file can contain only one model
    init = model$init(param.as.double)
    out  = model$rxode$solve(param.as.double, ev, init)
    out  = model$rxout(out)
    out = bind_cols(out, dose = rep(dose.nmol, dim(out)[1]))
    temp_df = bind_rows(temp_df, out) #ADD THIS LINE (change dfs to temp_df in this loop)
  }
  #Adds the drug name to data frame
  temp_df$drug = drugs #ADD THIS LINE
  
  #Appends temp_df to dfs to make a dataframe containing all infor for all drugs
  dfs <- rbind(dfs,temp_df) # ADD THIS LINE

}
write.csv(dfs, "all_res.csv")# can use other names for the output file

#modify output dataframe for plotting ----
# dfs.plot = dfs %>%
#   select(time,dose,T,DT,L,TL,D,Ttot,Ltot,Dtot) %>%
#   gather(variable,value,-c(time,dose)) %>%
#   arrange(desc(dose)) %>%
#   mutate(dose.label = paste(signif(dose,2),"nmol"),
#          dose.label = factor(dose.label,levels=unique(dose.label)))
# 
# #plot the results
# g = ggplot(dfs.plot,aes(x=time,y=value,color=dose.label))
# g = g + geom_line()
# g = g + scale.y.log10()
# g = g + facet_wrap(~variable)
# g = g + scale_x_units(units.input = "day",units.output = "month",increment = 1,t.end = tmax/30)
# g = g + labs(y="Concentration",color="Dose")
# 
# gg = ggsaveplot(width=6,height=6,dirs,"DosageEffect")
# grid.arrange(gg)






