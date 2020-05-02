source("ams_initialize_script.R")
source("SCIM_calculation.R")  
source("ivsc_2cmt_RR_V1.R")
dirs$rscript_name = "Task51_MultiDose_Figure_RealDrugs.R"
dirs$filename_prefix= str_extract(dirs$rscript_name,"^Task\\d\\d\\w?_")

model = ivsc_2cmt_RR_KdT0L0()

# Dose time, frequency, compartment, nominal dose
tmax = 5*28+28 #days
tau  = 28   #days
compartment = 2
n_points = 10

i_row = 0
result_metric = list()
result_sim    = list()
param.as.double = read.param.file(parameter_files["Tocilizumab"])[model$pin]

dose_range     = 10^(seq(-2,2,by=.25))
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
             dosereg  = ifelse(infusion==TRUE, "continual infusion", "every 4 weeks")) %>%
      bind_cols(sim,thy)
    
    out = plot_param(par, model, plot_flag = FALSE, infusion = infusion)
    sim = out$sim %>%
      mutate(dose = dose_mpk,
             TL0  = par$TL0_thy,
             dosereg  = ifelse(infusion==TRUE, "continual infusion", "every 4 weeks"))
    
    #create result table
    i_row = i_row + 1
    result_metric[[i_row]] = par
    result_sim[[i_row]] = sim
  }
}
metrics = bind_rows(result_metric)
sims    = bind_rows(result_sim) %>%
  mutate(SCIM = TL/TL0,
         SSIM = 1 - SCIM)


#plot results - curves ----
sims_plot = sims %>%
  select(dose,time,dosereg,D,L,SSIM) %>%
  mutate(SSIM = SSIM*100) %>%
  gather(cmt,value,c(D,L,SSIM)) %>%
  arrange(desc(dose)) %>%
  mutate(time     = time - 21,
         dose     = signif(dose,1)) %>%
  filter(dose %in% c(.1,.3,1,3,10,30,100)) %>%
  filter(!(cmt=="SSIM" & value<1e-4)) %>%
  mutate(dose     = factor(dose, levels = unique(dose)),
         cmt_name = plyr::mapvalues(cmt,
                                    c("D","L","SSIM"),
                                    c("Tocilizumab (D)", "IL-6 (L)", "SSIM (1-TLss/TL0)"))) %>%
  mutate(cmt_name = factor(cmt_name,levels = c("Tocilizumab (D)", "IL-6 (L)", "SSIM (1-TLss/TL0)")))


g = ggplot(sims_plot, aes(x=time-7,y=value, color = dose, group = dose))
g = g + geom_line() 
g = g + facet_grid(cmt_name~dosereg, switch = "y", scales = "free")
g = g + xgx_scale_x_time_units("day", limits = c(0,90)-7)
g = g + xgx_scale_y_log10()
g = g + labs(y = "Steady State Inhibition Metric (SSIM, %) or Concentration (nM)",
             color = "dose mg/kg")
ggsave(width = 7, height= 6, filename = "./figures/Task53_Dose_Range_Toci.png")
print(g)

#setup data for Lratio plot ----
metrics_toci = metrics %>%
  mutate(Lmax = max(Lss_thy),
         L0   = L0_sim,
         Lss  = Lss_sim,
         Lratio = 100*(Lss-L0)/(Lmax-L0),
         SSIM   = 100*(1-SCIM_sim)) %>%
  mutate(simulation = "Tocilizumab",
         dosereg = str_replace(dosereg,"every 4 weeks","every 2-4 week dosing")) %>%
  select(Lratio, SSIM, dosereg, simulation)
  

data_in = read.csv("results/Task22_2019-11-23_40e3.csv",stringsAsFactors = FALSE)
metrics_sensitivity = data_in %>%
  mutate(Lmax = Lss_thy,
         Lratio = (Lss_sim - L0)/(Lmax-L0),
         dosereg = ifelse(infusion==1,"continual infusion","every 2-4 week dosing"),
         SSIM = 1-SCIM_sim,
         simulation = "Global\nSensitivity\nAnalysis") %>%
  select(Lratio, SSIM, dosereg, simulation)

metrics_all = bind_rows(metrics_sensitivity, metrics_toci)

# Lratio plot ----
low  = 0
high = 100
g = ggplot(metrics_all, aes(x = Lratio*100, y = SSIM*100, 
                            color = simulation, size = simulation, alpha = simulation, shape = simulation))
#g = g + #xgx_scale_x_log10(limits = c(low,high))
#g = g + xgx_scale_x_reverselog10()#xgx_scale_y_log10(limits = c(low,high))
breaks = seq(0,100,by=25)
labels = paste0(breaks,"%")
limits = c(low, high)
g = g + scale_x_continuous(limits = limits, breaks = breaks, labels = labels)
g = g + scale_y_continuous(limits = limits, breaks = breaks, labels = labels)
g = g + ylim(c(0,100))
g = g + geom_point()
g = g + facet_wrap(~dosereg)
g = g + annotate("segment",x = low, y = low, xend = high, yend = high, color = "blue", linetype = "dashed")
g = g + labs(x = "Lratio = [Lss(dose) - L0]/[Lmax - L0]",
             y = "Steady State Inhibition Metric\nSSIM = 1 - TLss/TL0")
g = g + scale_color_manual(values=c("black","red"))
g = g + scale_alpha_manual(values=c(.5,1))
g = g + scale_size_manual(values=c(1,2.5))
ggsave(width = 6.5, height= 3, filename = "./figures/Task53_Lratio_Toci_Sensitivity.png")
print(g)

