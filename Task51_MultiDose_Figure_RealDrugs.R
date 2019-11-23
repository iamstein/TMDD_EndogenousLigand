source("ams_initialize_script.R")
source("SCIM_calculation.R")  
source("ivsc_2cmt_RR_V1.R")
dirs$rscript_name = "Task51_MultiDose_Figure_RealDrugs.R"
dirs$filename_prefix= str_extract(dirs$rscript_name,"^Task\\d\\d\\w?_")

model = ivsc_2cmt_RR_KdT0L0()

#read in parameter ranges to explore
param_minmax.in = readxl::read_excel("parameters/Task51_Param_Ranges.xlsx")
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

i_row = 0
i_drug = 0
result = list()
for (drug in drugs){ #loop over all the drugs in the list
  i_drug = i_drug + 1
  param.as.double = read.param.file(parameter_files[drug])[model$pin]
  
  for (dose_mpk in 10^(seq(-3,3,by=0.25))) {
    dose.nmol = dose_mpk*scale_mpk2nmol
    sim = lumped.parameters.simulation(model, param.as.double, dose.nmol, tmax, tau, compartment, infusion = infusion)
    thy = lumped.parameters.theory    (       param.as.double, dose.nmol,       tau)
        
    #all parameter values for the output table
    par = param.as.double %>% 
      t() %>% 
      as.data.frame() %>%
      mutate(dose_mpk    = dose_mpk,
             drug        = drug,
             target      = drug_target[i_drug],
             ligand      = drug_ligand[i_drug],
             order       = drug_order[i_drug],
             id          = i_row,
             tmax        = tmax,
             tau         = tau,
             dose_nmol   = dose.nmol)
      
    #create result table
    i_row = i_row + 1
    result[[i_row]] = bind_cols(sim,thy,par)
  }
}
results = bind_rows(result)
write.csv(results, file = "results/Task51_MultiDose_Figure.csv")

#plot results ----
data_plot = results %>%
  filter(drug %in% c("Tocilizumab","Siltuximab","Atezolizumab")) %>%
  select(dose_mpk, drug, target, ligand, order, SCIM_sim, AFIR_thy, SCIM_Lfold_adhoc_thy) %>%
  gather(key,value,-c(drug,dose_mpk,drug,target,ligand,order)) %>%
  arrange(order) %>%
  mutate(drug = factor(drug, levels = unique(drug)),
         target = paste("Target:",target),
         ligand = paste("Ligand:",ligand),
         key    = plyr::mapvalues(key,
                                  c("SCIM_sim","AFIR_thy","SCIM_Lfold_adhoc_thy"),
                                  c("SCIM simulation","AFIR theory","SCIM theory")))

g = ggplot(data_plot, aes(x=dose_mpk,y=1-value, color = key, linetype = key))
g = g + geom_line(size = 1, alpha = .5) 
g = g + facet_wrap(~drug+target+ligand)#, dir = "v", nrow = 2) )
g = g + xgx_scale_x_log10(breaks = 10^seq(-2,20,by=1))#, minor_breaks = 1)
breaks = c(0,90,99,99.9,99.99)/100
labels = paste0(breaks*100,"%")
g = g + xgx_scale_y_reverselog10(breaks = breaks, labels = labels)
#g = g + xgx_scale_y_log10()#, minor_breaks = 1)
g = g + scale_color_manual(values = c("AFIR theory"     = "red",
                                      "SCIM simulation" = "black",
                                      "SCIM theory" = "blue"))
g = g + scale_linetype_manual(values = c("AFIR theory" = "dotted",
                                         "SCIM simulation" = "solid",
                                         "SCIM theory" = "dashed"))
g = g + labs(x = "Dose (mg/kg) every 3 weeks",
             y = "Inhibition Metric",
             caption = "")
g = g + theme(legend.position = "top")

ggsave(width = 6.5, height= 4, filename = "./figures/Task51_DoseRange_Drugs.png")
print(g)
