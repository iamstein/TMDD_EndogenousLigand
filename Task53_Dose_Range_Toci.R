source("ams_initialize_script.R")
source("SCIM_calculation.R")  
source("ivsc_2cmt_RR_V1.R")
dirs$rscript_name = "Task51_MultiDose_Figure_RealDrugs.R"
dirs$filename_prefix= str_extract(dirs$rscript_name,"^Task\\d\\d\\w?_")

model = ivsc_2cmt_RR_KdT0L0()

# Dose time, frequency, compartment, nominal dose
tmax = 13*7+21 #days
tau  = 21   #days
compartment = 2
n_points = 10

i_row = 0
result_metric = list()
result_sim    = list()
param.as.double = read.param.file(parameter_files["Tocilizumab"])[model$pin]

dose_range     = 10^(seq(-1,2.5,by=.25))
infusion_range = c(TRUE, FALSE)

for (infusion in infusion_range) {
  for (dose_mpk in dose_range) {
    dose.nmol = dose_mpk*scale_mpk2nmol
    sim = lumped.parameters.simulation(model, param.as.double, dose.nmol, tmax, tau, compartment, infusion = infusion)
    thy = lumped.parameters.theory    (       param.as.double, dose.nmol,       tau,              infusion = infusion)
        
    #all parameter values for the output table
    par = param.as.double %>% 
      t() %>% 
      as.data.frame() %>%
      mutate(dose_mpk    = dose_mpk,
             id          = i_row,
             tmax        = tmax,
             tau         = tau,
             dose_nmol   = dose.nmol,
             dosereg  = ifelse(infusion==TRUE, "continual infusion", "every 3 weeks")) %>%
      bind_cols(sim,thy)
    
    out = plot_param(par, model, plot_flag = FALSE, infusion = infusion)
    sim = out$sim %>%
      mutate(dose = dose_mpk,
             TL0  = par$TL0_thy,
             dosereg  = ifelse(infusion==TRUE, "continual infusion", "every 3 weeks"))
    
    #create result table
    i_row = i_row + 1
    result_metric[[i_row]] = par
    result_sim[[i_row]] = sim
  }
}
metrics = bind_rows(result_metric)
sims    = bind_rows(result_sim) %>%
  mutate(SCIM = TL/TL0)


#plot results - curves ----
sims_plot = sims %>%
  select(dose,time,dosereg,D,L,SCIM) %>%
  gather(cmt,value,c(D,L,SCIM)) %>%
  arrange(desc(dose)) %>%
  mutate(time     = time - 21,
         dose     = signif(dose,1)) %>%
  filter(dose %in% c(.1,.3,1,3,10,30,100)) %>%
  mutate(dose     = factor(dose, levels = unique(dose)),
         cmt_name = plyr::mapvalues(cmt,
                                    c("D","L","SCIM"),
                                    c("Tocilizumab (D)", "IL-6 (L)", "SCIM (TLss/TL0)"))) %>%
  mutate(cmt_name = factor(cmt_name,levels = c("Tocilizumab (D)", "IL-6 (L)", "SCIM (TLss/TL0)")))


g = ggplot(sims_plot, aes(x=time,y=value, color = dose, group = dose))
g = g + geom_line() 
g = g + facet_grid(cmt_name~dosereg, switch = "y", scales = "free")
g = g + xgx_scale_x_time_units("day")
g = g + xgx_scale_y_log10()
g = g + labs(y = "SCIM or Concentration (nM)",
             color = "dose mg/kg")
ggsave(width = 7, height= 6, filename = "./figures/Task53_Dose_Range_Toci.png")
print(g)

#plot results - Lratio vs SCIM ----
metrics_plot = metrics %>%
  mutate(Lmax = max(Lss_sim),
         L0   = L0_sim,
         Lss  = Lss_sim,
         Lratio = (Lmax-Lss)/(Lmax-L0))


low  = min(metrics_plot$Lratio)
high = max(metrics_plot$Lratio)

g = ggplot(metrics_plot, aes(x = Lratio, y = SCIM_sim))
g = g + xgx_scale_x_log10()
g = g + xgx_scale_y_log10()
g = g + geom_point()
g = g + geom_line()
g = g + facet_wrap(~dosereg)
g = g + annotate("segment",x = low, y = low, xend = high, yend = high, linetype = "dotted")
g = g + labs(x = "Lratio = [Lmax - Lss(dose)]/[Lmax - L0]",
             y = "SCIM")
ggsave(width = 6.5, height= 3, filename = "./figures/Task53_Lratio_Toci.png")
print(g)

