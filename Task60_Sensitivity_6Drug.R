source("ams_initialize_script.R")
source("SCIM_calculation.R")  
source("ivsc_2cmt_RR_V1.R")
dirs$rscript_name = "Task60_Sensitivity_6Drug.R"
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



for (dose_original in c(10, 100)) 
{
  
  
  cat(paste("dose = ", dose_original, "mg/kg\n"))
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
    select(param_value, drug, param_name,
           SCIM_sim, SCIM_Lfold_adhoc_thy) %>% #Remove AFIR thy
    gather(key,value,-c(drug,param_name,param_value)) %>%
    mutate(drug = factor(drug, levels = unique(drug)),
           key  = plyr::mapvalues(key,
                                  c("SCIM_sim","AFIR_thy","SCIM_Lfold_adhoc_thy"),
                                  c("ASIR simulation","AFIR theory","ASIR theory")))
  
  g = ggplot(data_plot, aes(x=param_value,y=1-value, color = key, linetype = key))
  g = g + geom_line(size = 1, alpha = .5) 
  g = g + facet_grid(drug ~ param_name,scales = "free", switch = "y") 
  g = g + xgx_scale_x_log10(breaks = c(0.01,1,100))#, minor_breaks = 1) 
  breaks = c(0,90,99,99.9,99.99)/100
  labels = paste0(breaks*100,"%")
  g = g + xgx_scale_y_reverselog10(breaks = breaks, labels = labels)
  g = g + labs(x = "Parameter Value",
               y = "Steady State Inhibition Metric (SSIM)",
               caption = "")
  g = g + scale_color_manual(values = c("AFIR theory" = "red",
                                        "ASIR simulation" = "black",
                                        "ASIR theory" = "blue"))
  g = g + scale_linetype_manual(values = c("AFIR theory" = "dashed",
                                           "ASIR simulation" = "solid",
                                           "ASIR theory" = "dotted"))
  g = g + ggtitle(paste0("Dose = ",dose_original," mg/kg"))
  g = g + theme(legend.position = "top")
  filename = paste0("figures/Task60_sensitivity_6drug_",dose_original,"mpk.png")
  ggsave(width = 10, height= 6, filename = filename)
  print(g)
}
