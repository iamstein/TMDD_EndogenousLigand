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
    K = 2
    par = param.as.double %>% 
      t() %>% 
      as.data.frame() %>%
      bind_cols(sim,thy) %>%
      mutate(dose_mpk    = dose_mpk,
             drug        = drug,
             target      = drug_target[i_drug],
             ligand      = drug_ligand[i_drug],
             order       = drug_order[i_drug],
             id          = i_row,
             tmax        = tmax,
             tau         = tau,
             dose_nmol   = dose.nmol,
             koff_DT     = Kd_DT/kon_DT,
             Kss_DT      = Kd_DT + keDT/kon_DT,
             koff_TL     = Kd_TL/kon_TL,
             Kss_TL      = Kd_TL + keDT/kon_TL,
             assumption_SCIM_lt_30    = SCIM_Lfold_adhoc_thy < 0.30,
             assumption_drug_gg_Ttot  = Dss_thy > K*Ttotss_thy,
             assumption_drug_gg_KssDT = Dss_thy > K*Kss_DT,
             assumption_koffDT_gt_keT = koff_DT > keT,
             assumption_Dss_gt_Ccrit = Dss_thy > K*Ccrit_thy,
             assumption_Dss_gg_LssKssDT_KssTL = Dss_thy > K*Kss_DT*Lss_thy/Kss_TL,
             assumption_ODE_tolerance = Dss_thy/TLss_thy < 1e12,
             assumption_all_SCIM =    assumption_SCIM_lt_30 & 
               assumption_drug_gg_Ttot &
               assumption_Dss_gt_Ccrit &
               assumption_ODE_tolerance &
               assumption_Dss_gg_LssKssDT_KssTL)
      
    #create result table
    i_row = i_row + 1
    result[[i_row]] = par
  }
}
results = bind_rows(result)
write.csv(results, file = "results/Task51_MultiDose_Figure.csv")

#plot results ----
data_plot_all = results %>%
  select(dose_mpk, drug, target, ligand, order, SCIM_sim, AFIR_thy, SCIM_Lfold_adhoc_thy, assumption = assumption_drug_gg_Ttot) %>%
  gather(key,value,-c(drug,dose_mpk,drug,target,ligand,order,assumption)) %>%
  arrange(order) %>%
  mutate(drug = factor(drug, levels = unique(drug)),
         target = paste("Target:",target),
         ligand = paste("Ligand:",ligand),
         key    = plyr::mapvalues(key,
                                  c("SCIM_sim","AFIR_thy","SCIM_Lfold_adhoc_thy"),
                                  c("ASIR simulation","AFIR theory","ASIR theory")))

data_plot = data_plot_all %>%  
  filter(drug %in% c("Atezolizumab","Siltuximab","Tocilizumab"))

# data_assumption_false = data_plot %>%
#   filter(assumption == FALSE) %>%
#   filter(key == "ASIR simulation")

g = ggplot(data_plot, aes(x=dose_mpk,y=1-value, color = key, linetype = key))
g = g + geom_line(size = 1, alpha = .5) 
#g = g + geom_point(data = data_assumption_false, color = "red", show.legend = FALSE)
g = g + facet_wrap(~drug+target+ligand)#, dir = "v", nrow = 2) )
g = g + xgx_scale_x_log10(breaks = 10^seq(-2,20,by=1))#, minor_breaks = 1)
breaks = c(0,90,99,99.9,99.99)/100
labels = paste0(breaks*100,"%")
g = g + xgx_scale_y_reverselog10(breaks = breaks, labels = labels)
#g = g + xgx_scale_y_log10()#, minor_breaks = 1)
g = g + scale_color_manual(values = c("AFIR theory"     = "red",
                                      "ASIR simulation" = "black",
                                      "ASIR theory" = "blue"))
g = g + scale_linetype_manual(values = c("AFIR theory" = "dotted",
                                         "ASIR simulation" = "solid",
                                         "ASIR theory" = "dashed"))
g = g + labs(x = "Dose (mg/kg) every 3 weeks",
             y = "Steady State Inhibition Metric\nSSIM = 1 - TLss/TL0",
             caption = "")
g = g + theme(legend.position = "top")

ggsave(width = 6.5, height= 4, filename = "./figures/Task51_DoseRange_Drugs.png")
print(g)

g = g %+% data_plot_all
ggsave(width = 6.5, height= 7, filename = "./figures/Task51_DoseRange_All6_Drugs.png")
print(g)

stop()
p = results %>%
  filter(drug == "Atezolizumab" & round(dose_mpk)==3)
out = plot_param(p,model,infusion = FALSE, plot_flag = FALSE)
g = out$plot
g = g + ggtitle("Atezo parameters are not quite right -- too much target in blood")
print(g)
