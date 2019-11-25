# Setup and Read Data
source("ams_initialize_script.R")
source("SCIM_calculation.R")  
source("ivsc_2cmt_RR_V1.R")
library(RxODE)
model = ivsc_2cmt_RR_KdT0L0()
dirs$rscript_name = "Task52_Global_Sensitivity_Analysis.R"
dirs$filename_prefix= str_extract(dirs$rscript_name,"^Task\\d\\d\\w?_")

#data_in = read.csv("results/Task19_2019-11-17_20e3.csv",stringsAsFactors = FALSE)
data_in = read.csv("results/Task22_2019-11-25_40e3.csv",stringsAsFactors = FALSE)


# Compute various quantities for comparing AFIR and SCIM, theory and simulation
### Using adhoc theory calculation for SCIM

#put data into categories ----
data    = data_in
 
#check the assumptions of the data ----
K = 2
data = data %>%
  mutate(AFIR_SCIM_pcterr = abs((AFIR_thy             - SCIM_sim)),#/SCIM_sim),
         SCIM_SCIM_pcterr = abs((SCIM_Lfold_adhoc_thy - SCIM_sim))) %>% #/SCIM_sim)) %>%
  mutate(assumption_AFIR_lt_30    = AFIR_thy < 0.30,
         assumption_SCIM_lt_30    = SCIM_Lfold_adhoc_thy < 0.30,
         assumption_drug_gg_Ttot  = Dss_thy > K*Ttotss_thy,
         assumption_drug_gg_KssDT = Dss_thy > K*Kss_DT,
         assumption_koffDT_gt_keT = koff_DT > keT,
#         assumption_koffTL_fast   = koff_TL > 1/30,
         assumption_Dss_gt_Ccrit = Dss_thy > K*Ccrit_thy,
         assumption_Dss_gg_LssKssDT_KssTL = Dss_thy > K*Kss_DT*Lss_thy/Kss_TL,
         assumption_ODE_tolerance = Dss_thy/TLss_thy < 1e12,
#         assumption_T0simple    = T0/(ksynT/keT) > 0.5 & T0/(ksynT/keT) < 2, #the simple formula works for T0
         assumption_L_noaccum    = Lfold_thy <= 1.1, #then SCIM = AFIR
#         assumption_Tss_gt_Lss  = Tss_sim > Lss_sim,
         assumption_all_AFIR    = assumption_AFIR_lt_30 & 
                                  #assumption_drug_gg_Ttot &
                                  #assumption_drug_gg_KssDT &
                                  #assumption_koffDT_gt_keT & 
                                  #assumption_koffTL_fast &           
                                  assumption_Dss_gg_LssKssDT_KssTL &
                                  assumption_Dss_gt_Ccrit &
                                  assumption_ODE_tolerance &
                                  assumption_L_noaccum,
         assumption_all_SCIM =    assumption_SCIM_lt_30 & 
                                  assumption_drug_gg_Ttot &
                                  #assumption_drug_gg_KssDT &
                                  #assumption_koffDT_gt_keT & 
                                  #assumption_koffTL_fast &    
                                  assumption_Dss_gt_Ccrit &
                                  assumption_ODE_tolerance &
                                  assumption_Dss_gg_LssKssDT_KssTL)

# histogram of AFIR_theory and SCIM_sim error ----
data1 = data %>%
  select(pcterr = AFIR_SCIM_pcterr, assumptions_all_true = assumption_all_AFIR, infusion) %>%
  mutate(metric = "AFIR")

data2 = data %>%
  select(pcterr = SCIM_SCIM_pcterr, assumptions_all_true = assumption_all_SCIM, infusion) %>%
  mutate(metric = "ASIR")


data_plot = bind_rows(data1,data2) %>%
  mutate(assumptions_all_true = ifelse(assumptions_all_true,"yes","no"),
         infusion             = ifelse(infusion==1,"continual infusion","every 2-4 week dosing"))
    

g = ggplot(data_plot, aes(pcterr*100, fill = assumptions_all_true))
g = g + geom_histogram()
g = g + facet_grid(infusion~metric, scales = "free_x", switch = "y")
breaks = 10^seq(-3,2,by=1)
g = g + xgx_scale_x_log10(limits = c(2e-4,200), breaks = breaks, labels = paste0(breaks,"%"))
g = g + labs(x    = "Difference between\ntheoretical metric (AFIR or ASIR)\nand ASIR simulation",
             y    = "Number of Simulations",
             fill = "Assumptions\nfor metric\nare all true")
g = g + scale_fill_manual(values = c(yes="black",no="grey80"))
g = g + geom_vline(xintercept = 10, color = "grey30")
print(g)
ggsave(width = 8, height= 4, filename = "./figures/Task52_GlobalSensitivityAnalysis.png")


#explore the case with teh largest error ----
stop()
x = data %>%
  filter(assumption_all_SCIM == TRUE) %>%
  arrange(desc(SCIM_SCIM_pcterr))

i=1
xi = x[i,]
plot_param(xi,model,infusion = xi$infusion)


