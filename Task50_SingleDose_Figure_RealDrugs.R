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
tmax = 12*7+21 #days
tau  = 21   #days
compartment = 2
infusion = FALSE

i_row = 0
i_drug = 0
result = list()
for (drug in drugs){ #loop over all the drugs in the list
  i_drug = i_drug + 1
  param.as.double = read.param.file(parameter_files[drug])[model$pin]
  
  for (dose_mpk in 10) { #10^(seq(-3,3,by=0.25))) {
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
    
    par = bind_cols(sim,thy,par)
    out = plot_param(par,model,plot_flag = FALSE, infusion = FALSE)
    
    sim = out$sim %>%
      mutate(dose        = dose_mpk,
             drug        = drug,
             target      = drug_target[i_drug],
             ligand      = drug_ligand[i_drug],
             order       = drug_order[i_drug])
    
    #create result table
    i_row = i_row + 1
    result[[i_row]] = sim
  }
}
results = bind_rows(result)
write.csv(results, file = "results/Task50_SingleDose_Figure.csv")

#plot results ----
data_plot = results %>%
  gather(cmt,value,c(D,T,DT,L,TL)) %>%
  arrange(order) %>%
  mutate(drug = factor(drug, levels = unique(drug)),
         target = paste("Target:",target),
         ligand = paste("Ligand:",ligand)) %>%
  mutate(time = time - 21)

data_last = data_plot %>%
  filter(time == max(time)) %>%
  mutate(time = time + 3)

g = ggplot(data_plot,aes(x=time/7,y=value, color = cmt, group= cmt))
g = g + geom_line()
g = g + geom_text(data = data_last, aes(label = cmt), show.legend = FALSE, hjust=0)
g = g + scale_x_continuous(breaks = seq(-3,100,by=3),
                           lims   = c(-3,tmax/7+3))
g = g + xgx_scale_y_log10()
g = g + facet_wrap(~drug+target+ligand, dir = "v", nrow = 2) 
g = g + labs(y = "Concentration (nM)", color = "")
g = g + ggtitle(paste("Dose:",dose_mpk,"mg/kg every three weeks"))
print(g)
ggsave(width = 6, height= 6, filename = "./figures/Task50_SingleDose_Drugs.png")



