# Setup and Read Data
source("ams_initialize_script.R")
source("SCIM_calculation.R")  
source("ivsc_2cmt_RR_V1.R")
library(RxODE)
model = ivsc_2cmt_RR_KdT0L0()
dirs$rscript_name = "Task52_Global_Sensitivity_Analysis.R"
dirs$filename_prefix= str_extract(dirs$rscript_name,"^Task\\d\\d\\w?_")

#data_in = read.csv("results/Task19_2019-11-17_20e3.csv",stringsAsFactors = FALSE)
data_in = read.csv("results/Task22_2020-10-17_40e3.csv",stringsAsFactors = FALSE)


# Compute various quantities for comparing AFIR and SCIM, theory and simulation
### Using adhoc theory calculation for SCIM

#put data into categories ----
data    = data_in
 
#check the assumptions of the data ----
K = 2
data = data %>%
  mutate(AFIR_SCIM_differr  = abs((AFIR_thy             - SCIM_sim)),
         SCIM_SCIM_differr  = abs((SCIM_Lfold_adhoc_thy - SCIM_sim)), 
         SCIM_SCIM_ratioerr = abs((SCIM_Lfold_adhoc_thy - SCIM_sim))/SCIM_sim) %>%
  mutate(assumption_AFIR_lt_30    = AFIR_thy < 0.30,
         assumption_SCIM_lt_30    = SCIM_Lfold_adhoc_thy < 0.30,
         str_SCIM_lt_30          = ifelse(assumption_SCIM_lt_30,"o SCIM < 30%","x SCIM >= 30%"),
         assumption_drug_gg_LssKssDT_KssTL = Dss_thy > K*Kss_DT*Lss_thy/Kss_TL,
         str_drug_gg_LssKssDT_KssTL = ifelse(assumption_drug_gg_LssKssDT_KssTL, "o D > 2*Lss*KssDT/KssTL","x D < 2*Lss*KssDT/KssTL"),
         assumption_drug_gg_Ttot  = Dss_thy > 2*Ttotss_thy,
         str_drug_gg_Ttot         = ifelse(assumption_drug_gg_Ttot,"o D > 2*Ttot","x D < 2*Ttot"),
         assumption_drug_gg_Ccrit = Dss_thy > 2*Ccrit_thy,
         str_drug_gg_Ccrit         = ifelse(assumption_drug_gg_Ccrit,"o D > 2*Ccrit","x D < 2*Ccrit"),
         assumption_ODE_tolerance = Dss_thy/TLss_thy < 1e12,
         assumption_L_noaccum    = Lfold_thy <= 1.1, #then SCIM = AFIR
         assumption_all_AFIR    = assumption_AFIR_lt_30 & 
                                  assumption_drug_gg_Ttot &
                                  assumption_drug_gg_Ccrit &
                                  assumption_drug_gg_LssKssDT_KssTL &
                                  assumption_ODE_tolerance &
                                  assumption_L_noaccum,
         assumption_all_SCIM =    assumption_SCIM_lt_30 & 
                                  assumption_drug_gg_Ttot &
                                  assumption_drug_gg_Ccrit &
                                  assumption_drug_gg_LssKssDT_KssTL &
                                  assumption_ODE_tolerance) %>%
  mutate(assumption_SCIM_list = paste(assumption_SCIM_lt_30, 
                                     assumption_drug_gg_Ttot,
                                     assumption_drug_gg_Ccrit,
                                     assumption_drug_gg_LssKssDT_KssTL),
         assumption_SCIM_str  = paste(ifelse(assumption_SCIM_lt_30, "SCIM<30%; ", "-"),
                                      ifelse(assumption_drug_gg_Ttot, "D>2*Ttot; ", "-"),
                                      ifelse(assumption_drug_gg_Ccrit, "D>2*Ccrit; ", "-"),
                                      ifelse(assumption_drug_gg_LssKssDT_KssTL, "D>Lss*KssDT/KssTL ", "-"),
                                      sep = "\n")) %>%
  arrange(assumption_SCIM_list) %>%
  mutate(assumption_SCIM_str = factor(assumption_SCIM_str, levels = unique(assumption_SCIM_str)))

# histogram of AFIR_theory and SCIM_sim error ----
data1 = data %>%
  select(differr = AFIR_SCIM_differr, assumptions_all_true = assumption_all_AFIR, infusion) %>%
  mutate(metric = "AFIR")

data2 = data %>%
  select(differr = SCIM_SCIM_differr, assumptions_all_true = assumption_all_SCIM, infusion) %>%
  mutate(metric = "ASIR")

data_plot = bind_rows(data1,data2) %>%
  mutate(assumptions_all_true = ifelse(assumptions_all_true,"yes","no"),
         infusion             = ifelse(infusion==1,"continual infusion","every 2-4 week dosing"))

g = ggplot(data_plot, aes(differr*100, fill = assumptions_all_true))
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

# creat histogram based on different assumptions not being true ----
g = ggplot(filter(data, infusion == 1), aes(100*SCIM_SCIM_differr))
g = g + facet_wrap(~str_drug_gg_Ccrit)
#g = g + facet_grid(str_drug_gg_Ccrit + str_SCIM_lt_30 ~ str_drug_gg_Ttot + str_drug_gg_LssKssDT_KssTL, scales = "free_y", switch =)
#g = g + geom_histogram(aes(y=..density..))
g = g + geom_histogram()
breaks = 10^seq(-3,2,by=1)
g = g + xgx_scale_x_log10(limits = c(2e-4,200), breaks = breaks, labels = paste0(breaks,"%"))
g = g + labs(x    = "Difference between\ntheoretical metric (AFIR or ASIR)\nand ASIR simulation",
             y    = "Number of Simulations",
             fill = "Assumptions \nthat are true")
g = g + geom_vline(xintercept = 10, color = "grey30")
#colors = scales::seq_gradient_pal("blue", "red", "Lab")(seq(0,1,length.out=length(unique(data$assumption_SCIM_str))))
#g = g + scale_fill_manual(values = colors)

print(g)
ggsave(width = 8, height= 4, filename = "./figures/Task52_GlobalSensitivityAnalysis_SSIM_Ccrit_differr.png")


# create colored histogram based on different assumptions not being true ----
g = g + facet_grid()
g = ggplot(filter(data, infusion == 1), aes(100*SCIM_SCIM_ratioerr))
g = g + facet_wrap(~str_drug_gg_Ccrit)
#g = g + facet_grid(str_drug_gg_Ccrit + str_SCIM_lt_30 ~ str_drug_gg_Ttot + str_drug_gg_LssKssDT_KssTL, scales = "free_y", switch =)
#g = g + geom_histogram(aes(y=..density..))
g = g + geom_histogram()
breaks = 10^seq(-3,2,by=1)
g = g + xgx_scale_x_log10(limits = c(2e-4,200), breaks = breaks, labels = paste0(breaks,"%"))
g = g + labs(x    = "Ratio difference between \ntheoretical ASIR and simulation",
             y    = "Number of Simulations",
             fill = "Assumptions \nthat are true")
g = g + geom_vline(xintercept = 10, color = "grey30")
#colors = scales::seq_gradient_pal("blue", "red", "Lab")(seq(0,1,length.out=length(unique(data$assumption_SCIM_str))))
#g = g + scale_fill_manual(values = colors)

print(g)
ggsave(width = 8, height= 4, filename = "./figures/Task52_GlobalSensitivityAnalysis_SSIM_Ccrit_Ratioerr.png")

# create colored histogram based on different assumptions not being true ----
g = g + facet_grid()
g = ggplot(filter(data, infusion == 1), aes(100*SCIM_SCIM_ratioerr))
#g = g + facet_wrap(~str_drug_gg_Ccrit)
g = g + facet_grid(str_drug_gg_Ccrit + str_SCIM_lt_30 ~ str_drug_gg_Ttot + str_drug_gg_LssKssDT_KssTL, scales = "free_y", switch =)
#g = g + geom_histogram(aes(y=..density..))
g = g + geom_histogram()
breaks = 10^seq(-3,2,by=1)
g = g + xgx_scale_x_log10(limits = c(2e-4,200), breaks = breaks, labels = paste0(breaks,"%"))
g = g + labs(x    = "Ratio difference between \ntheoretical ASIR and simulation",
             y    = "Number of Simulations",
             fill = "Assumptions \nthat are true")
g = g + geom_vline(xintercept = 10, color = "grey30")
#colors = scales::seq_gradient_pal("blue", "red", "Lab")(seq(0,1,length.out=length(unique(data$assumption_SCIM_str))))
#g = g + scale_fill_manual(values = colors)

print(g)
ggsave(width = 8, height= 4, filename = "./figures/Task52_GlobalSensitivityAnalysis_SSIM_4assum_Ratioerr.png")


#explore the case with teh largest error ----
x = data %>%
  filter(assumption_all_SCIM == TRUE) %>%
  arrange(desc(SCIM_SCIM_ratioerr))

i=1
xi = x[i,]
summi = plot_param(xi,model,infusion = xi$infusion)

print(summi$compare)



