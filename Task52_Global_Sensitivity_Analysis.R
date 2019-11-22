# Setup and Read Data
source("ams_initialize_script.R")
source("SCIM_calculation.R")  
source("ivsc_2cmt_RR_V1.R")
library(RxODE)
model = ivsc_2cmt_RR_KdT0L0()
dirs$rscript_name = "Task52_Global_Sensitivity_Analysis.R"
dirs$filename_prefix= str_extract(dirs$rscript_name,"^Task\\d\\d\\w?_")

data_in = read.csv("results/Task19_2019-11-17_20e3.csv",stringsAsFactors = FALSE)


# Compute various quantities for comparing AFIR and SCIM, theory and simulation
### Using adhoc theory calculation for SCIM

#put data into categories ----
data    = data_in
 
#check the assumptions of the data ----
data = data %>%
  mutate(AFIR_SCIM_pcterr = abs((AFIR_thy             - SCIM_sim)/SCIM_sim),
         SCIM_SCIM_pcterr = abs((SCIM_Lfold_adhoc_thy - SCIM_sim)/SCIM_sim),
         SCIM_SCIM_pcterr = ifelse(SCIM_SCIM_pcterr < 0.001, 0.0001, SCIM_SCIM_pcterr)) %>%
  mutate(assumption_AFIR_lt_30    = AFIR_thy < 0.30,
         assumption_SCIM_lt_30    = SCIM_Lfold_adhoc_thy < 0.30,
         assumption_drug_gg_T0    = Dss_thy > 5*Ttotss_thy,
         assumption_drug_gg_KssDT = Dss_thy > 5*Kss_DT,
         assumption_koffDT_gt_keT = koff_DT > keT,
         assumption_koffTL_fast   = koff_TL > 1/30,
         assumption_Cavgss_gt_Ccrit = Cavgss_thy > 10*Ccrit_thy,
         assumption_Cavgss_gg_LssKssDT_KssTL = Cavgss_thy > 5*Kss_DT*Lss_thy/Kss_TL,
         assumption_ODE_tolerance = Cavgss_thy/TLss_thy < 1e12,
#         assumption_T0simple    = T0/(ksynT/keT) > 0.5 & T0/(ksynT/keT) < 2, #the simple formula works for T0
         assumption_L_noaccum    = Lfold_thy <= 1.01, #then SCIM = AFIR
#         assumption_Tss_gt_Lss  = Tss_sim > Lss_sim,
         assumption_all_AFIR    = assumption_AFIR_lt_30 & 
                                  assumption_drug_gg_T0 &
                                  assumption_drug_gg_KssDT &
                                  #assumption_koffDT_gt_keT & 
                                  assumption_koffTL_fast &           
                                  assumption_Cavgss_gg_LssKssDT_KssTL &
                                  assumption_Cavgss_gt_Ccrit &
                                  assumption_ODE_tolerance &
                                  assumption_L_noaccum,
         assumption_all_SCIM =    assumption_SCIM_lt_30 & 
                                  assumption_drug_gg_T0 &
                                  assumption_drug_gg_KssDT &
                                  #assumption_koffDT_gt_keT & 
                                  assumption_koffTL_fast &    
                                  assumption_Cavgss_gt_Ccrit &
                                  assumption_ODE_tolerance &
                                  assumption_Cavgss_gg_LssKssDT_KssTL)

# histogram of AFIR_theory and SCIM_sim error ----
data1 = data %>%
  select(pcterr = AFIR_SCIM_pcterr, assumptions_all_true = assumption_all_AFIR) %>%
  mutate(metric = "AFIR")

data2 = data %>%
  select(pcterr = SCIM_SCIM_pcterr, assumptions_all_true = assumption_all_SCIM) %>%
  mutate(metric = "SCIM")


data_plot = bind_rows(data1,data2) %>%
  mutate(assumptions_all_true = ifelse(assumptions_all_true,"yes","no"))
    

g = ggplot(data_plot, aes(pcterr*100, fill = assumptions_all_true))
g = g + geom_histogram()
g = g + facet_wrap(~metric, scales = "free_x")
breaks = 10^seq(-4,2,by=2)
g = g + xgx_scale_x_log10(limits = c(2e-4,200), breaks = breaks, labels = paste0(breaks,"%"))
g = g + labs(x    = "Percent Error between\ntheoretical metric (AFIR or SCIM)\nand SCIM simulation",
             y    = "Number of Simulations",
             fill = "Assumptions\nfor metric\nare all true")
g = g + scale_fill_manual(values = c(yes="black",no="grey70"))
print(g)
ggsave(width = 6.5, height= 4, filename = "./figures/Task52_GlobalSensitivityAnalysis.png")

