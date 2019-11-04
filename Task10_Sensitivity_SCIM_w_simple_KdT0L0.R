source("ams_initialize_script.R")
source("SCIM_calculation.R")  
source("ivsc_2cmt_RR_V1.R")
dirs$rscript_name = "Task10_Sensitivity_SCIM_w_simple_KssT0L0.R"
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
infusion = FALSE
n_points = 10



for (dose_original in c(10,100,1000)) 
{


result = list()
i_row = 0
for (drug in drugs){ #loop over all the drugs in the list
  param_as_double_original = read.param.file(parameter_files[drug])[model$pin]

  # Iterate all of the parameters for a single drug.
  i_param_name = 0
  for (param_name in param_minmax$Parameter){
    param_values = with(filter(param_minmax,Parameter==param_name),
                        lseq(min,max,n_points))
    
    for (param_value in param_values) {
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
      sim = lumped.parameters.simulation(model, param_as_double, dose_nmol, tmax, tau, compartment, infusion = infusion)
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
write.csv(results, file = "results/Task10_Sensitivity_SCIM_w_simple_KdT0L0.csv")

#check the initial condition and steady state ----
check = results %>%
  select(param_value, drug, param_name, TLss_frac_change, TL0_05tau_frac_change) %>%
  gather(key,value,-c(param_value,drug,param_name)) %>%
  mutate(value      = abs(value))

g = ggplot(check,aes(value))
g = g + geom_histogram()
g = g + facet_grid(drug~key, switch = "y", scales = "free_x")
g = g + xgx_scale_x_log10()
g = xgx_save(8,4,dirs,paste0("Check_Init_SteadyState_",dose_original,"mpk"),status = "")
print(g)


#plot results ----
data_plot = results %>%
  select(param_value, drug, param_name, TLss_frac_change, TL0_05tau_frac_change,
                SCIM_sim, AFIR_thy, SCIM_adhoc_thy) %>%
  gather(key,value,-c(drug,param_name,param_value,TLss_frac_change, TL0_05tau_frac_change)) %>%
  mutate(AFIR_SCIM  = ifelse(str_detect(key,"AFIR"),"AFIR","SCIM"),
         theory_sim = ifelse(str_detect(key,"sim"),"sim","thy"),
         approx     = str_extract(key,"_\\w+_"),
         approx     = str_replace(approx,"_",""),
         approx     = ifelse(is.na(approx),"none",approx))

threshold = 0.05
data_errss = data_plot %>%
  filter(abs(TLss_frac_change)>=threshold) %>%
  filter(theory_sim == "sim")
print(paste0(nrow(data_errss),": Number of rows with TLss_frac_change > 0.1"))

data_err0 = data_plot %>%
  filter(abs(TL0_05tau_frac_change)>=threshold) %>%
  filter(theory_sim == "sim")
print(paste0(nrow(data_err0),": Number of rows with TL0_05tau_frac_change > 0.1"))

g = ggplot(data_plot, aes(x=param_value,y=value, color = key, linetype = key))
g = g + geom_line(size = 1, alpha = .5) 
g = g + geom_point(data = data_err0,  mapping = aes(x=param_value,y=value), color = "red")
g = g + geom_point(data = data_errss, mapping = aes(x=param_value,y=value), color = "purple")
g = g + facet_grid(drug ~ param_name,scales = "free", switch = "y") 
g = g + xgx_scale_x_log10(breaks = c(1e-4,1e-2,1,100,1e4))#, minor_breaks = 1) 
g = g + xgx_scale_y_log10(breaks = c(1e-4,1e-2,1,100,1e4))#, minor_breaks = 1)
g = g + scale_color_manual(values = c("AFIR_thy" = "red",
                                      "SCIM_sim" = "black",
                                      "SCIM_adhoc_thy" = "green4",
                                      "SCIM_thy" = "blue"))
g = g + scale_linetype_manual(values = c("AFIR_thy" = "dotted",
                                         "SCIM_sim" = "solid",
                                         "SCIM_adhoc_thy" = "longdash",
                                         "SCIM_thy" = "dashed"))
g = g + ggtitle(paste0("Dose = ",dose_original," mg/kg"))
g = xgx_save(17,10,dirs,paste0("reparKdT0L0_",dose_original,"mpk"),status = "")
print(g)

# ----
}
