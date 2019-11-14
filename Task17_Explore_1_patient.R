source("ams_initialize_script.R")
source("SCIM_calculation.R")  
source("ivsc_2cmt_RR_V1.R")
dirs$rscript_name = "Task17_Explore_1_patient.R"
dirs$filename_prefix= str_extract(dirs$rscript_name,"^Task\\d\\d\\w?_")


data_in = read.csv("results/Task15_2019-11-12_10e3.csv",stringsAsFactors = FALSE)
data_in$id = seq(1,10000)

#simulate a patient where theory and simulation disagree ----
source("SCIM_calculation.R")  
source("ivsc_2cmt_RR_V1.R")
model = ivsc_2cmt_RR_KdT0L0()
library(RxODE)

id_simulate = 1770

param = data %>%
  filter(id==id_simulate)

assumptions = param %>%
  select(contains("assumption")) %>%
  t() 
kable(assumptions)

tmax = 16*7 #days
tau  = param$tau   #days
dose_nmol = param$dose_nmol
compartment = 2
infusion = TRUE

nam   = names(param)
param_as_double = param %>%
  as.numeric() %>%
  setNames(nam)
param_as_double = param_as_double[model$pin]

ev = eventTable(amount.units="nmol", time.units="days")
sample.points = c(seq(0, tmax, 0.1), 10^(-3:0)) # sample time, increment by 0.1
sample.points = sort(sample.points)
sample.points = unique(sample.points)
ev$add.sampling(sample.points)
if (infusion == FALSE) {
  ev$add.dosing(dose=dose_nmol, start.time = tau, nbr.doses=floor(tmax/tau), dosing.interval=tau, dosing.to=compartment)
} else {
  ev$add.dosing(dose=dose_nmol, start.time = tau, nbr.doses=floor(tmax/tau)+1, dosing.interval=tau, dosing.to=compartment, dur = tau)
}  

sim = lumped.parameters.simulation(model, param_as_double, dose_nmol, tmax, tau, compartment, infusion)
thy = lumped.parameters.theory    (       param_as_double, dose_nmol,       tau,              infusion)

sim_rename = sim
nam = names(sim_rename) %>%
  str_replace_all("_sim$","")
names(sim_rename) = nam
sim_rename$type = "sim"

thy_rename = thy
nam = names(thy_rename) %>%
  str_replace_all("_thy$","")
names(thy_rename) = nam
thy_rename$type = "thy"

compare = bind_rows(sim_rename,thy_rename) %>%
  select(type,Dss,T0,L0,TL0,Ttotss,Lss,TLss,AFIR,SCIM)
kable(compare)


init = model$init(param_as_double)
out  = model$rxode$solve(model$repar(param_as_double), ev, init)
out  = model$rxout(out)

out_plot = out %>%
  select(time,D,T,DT,L,TL) %>%
  gather(cmt,value,-time)
out_last = out_plot[(out$time==max(out$time)),]

g = ggplot(out_plot,aes(x=time,y=value, color = cmt, group= cmt))
g = g + geom_line()
g = g + geom_label(data = out_last, aes(label = cmt), show.legend = FALSE, hjust=1)
g = g + geom_vline(xintercept = tau, linetype = "dotted")
g = g + xgx_scale_x_time_units(units_dataset = "days", units_plot = "weeks")
g = g + xgx_scale_y_log10()
g = g + labs(y = "Concentration (nm)", color = "")
g = g + ggtitle(paste0(  "id = ",param$id,
                         "\nAFIR_thy  = ",signif(thy$AFIR_thy,2),
                         "\nAFIR_sim  = ",signif(sim$AFIR_sim,2),
                         "\nSCIM_thy = ",signif(thy$SCIM_adhoc_thy,2),
                         "\nSCIM_sim = ",signif(sim$SCIM_sim,2)))
filepref = paste0("one_parameter_set_",id_simulate)
g = xgx_save(5,5,dirs,filepref,"")
print(g)



