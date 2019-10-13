source("ams_initialize_script.R")
source("SCIM_calculation.R")  
source("ivsc_2cmt_RR_V1.R")
dirs$rscript_name = "Task10c_Sensitivity_SCIM_w_simple_KdT0L0_Atezo_BigT0.R"
dirs$filename_prefix= str_extract(dirs$rscript_name,"^Task\\d\\d\\w?_")

model = ivsc_2cmt_RR_KdT0L0()

#read in parameter ranges to explore
param_minmax.in = readxl::read_excel("parameters/Task10_Param_Ranges.xlsx")
param_minmax = param_minmax.in %>%
  as.data.frame() %>%
  select(Parameter,min,max,units,fixed) %>%
  filter(!(is.na(fixed))) %>%
  filter(fixed==0)
rownames(param_minmax) = param_minmax$Parameter

# Dose time, frequency, compartment, nominal dose
tmax = 52*7 #days
tau  = 21   #days
compartment = 2
n_points = 10
infusion = TRUE

for (dose_original in c(10)) {
result = list()
i_row = 0
for (drug in "Atezolizumab"){ #loop over all the drugs in the list
  param_as_double_original = read.param.file(parameter_files[drug])[model$pin]

  # Iterate all of the parameters for a single drug.
  i_param_name = 0
  for (param_name in "T0"){
    param_values = with(filter(param_minmax,Parameter==param_name),
                        lseq(min,max,n_points))
    
    for (param_value in param_values[n_points]) {
      param_as_double = param_as_double_original
      
      if (param_name == 'dose'){
        dose_nmol = param_value*scale_mpk2nmol
      } else {
        dose_nmol = dose_original*scale_mpk2nmol
        param_as_double[param_name] = param_value
      }
      
      #key lines for computing theory and simulation
      param.as.double = param_as_double
      dose.nmol       = dose_nmol
      sim = lumped.parameters.simulation(model, param_as_double, dose_nmol, tmax, tau, compartment, infusion)
      thy = lumped.parameters.theory    (       param_as_double, dose_nmol,       tau)
      
      #all parameter values for the output table
      par = param_as_double %>% 
        t() %>% 
        as.data.frame() %>%
        mutate(param_name  = param_name,
               param_value = param_value,
               drug        = drug)
    
      #create result table
      i_row = i_row + 1
      result[[i_row]] = bind_cols(sim,thy,par)
        
    }
  }
}
results = bind_rows(result)
}


ev = eventTable(amount.units="nmol", time.units="days")
sample.points = c(seq(0, tmax, 0.1), 10^(-3:0)) # sample time, increment by 0.1
sample.points = sort(sample.points)
sample.points = unique(sample.points)
ev$add.sampling(sample.points)

#add dur  tau for a long infusion
if (infusion == FALSE) {
  ev$add.dosing(dose=dose.nmol, start.time = tau, nbr.doses=floor(tmax/tau), dosing.interval=tau, dosing.to=compartment)
} else {
  ev$add.dosing(dose=dose.nmol, start.time = tau, nbr.doses=floor(tmax/tau)+1, dosing.interval=tau, dosing.to=compartment, dur = tau)
}  

init = model$init(param.as.double)
out  = model$rxode$solve(model$repar(param.as.double), ev, init)
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
g = g + ggtitle(paste0(  "AFIR_thy  = ",signif(results$AFIR_thy,2),
                       "\nAFIR_sim  = ",signif(results$AFIR_sim,2),
                       "\nSCIM_sim = ",signif(results$SCIM_sim,2)))
print(g)
