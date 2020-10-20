# Setup and Read Data
source("ams_initialize_script.R")
source("SCIM_calculation.R")  
source("ivsc_2cmt_RR_V1.R")
library(RxODE)
model = ivsc_2cmt_RR_KdT0L0()
dirs$rscript_name = "Task52_Global_Sensitivity_Analysis.R"
dirs$filename_prefix= str_extract(dirs$rscript_name,"^Task\\d\\d\\w?_")

data_in = read.csv("results/Task22_2020-10-17_40e3.csv",stringsAsFactors = FALSE)


# Compute various quantities for comparing AFIR and SCIM, theory and simulation
### Using adhoc theory calculation for SCIM

#put data into categories ----
data    = data_in
 
#check the assumptions of the data ----
K = 2
data = data %>%
  mutate(AFIR_SCIM_differr  = ((AFIR_thy             - SCIM_sim)),
         SCIM_SCIM_differr  = ((SCIM_Lfold_adhoc_thy - SCIM_sim)), 
         SCIM_SCIM_ratio    = SCIM_Lfold_adhoc_thy/SCIM_sim,
         TLss_ratio_sim     = Tss_sim * Lss_sim / TLss_sim,
         Kss_ratio       = Kss_TL/TLss_ratio_sim) %>%
  mutate(assumption_AFIR_lt_30    = AFIR_thy < 0.30,
         assumption_SCIM_lt_30    = SCIM_Lfold_adhoc_thy < 0.30,
         str_SCIM_lt_30          = ifelse(assumption_SCIM_lt_30,"o SCIM < 30%","x SCIM >= 30%"),
         assumption_drug_gg_LssKssDT_KssTL = Dss_thy > K*Kss_DT*Lss_thy/Kss_TL,
         str_drug_gg_LssKssDT_KssTL = ifelse(assumption_drug_gg_LssKssDT_KssTL, "o D > 2*Lss*KssDT/KssTL","x D < 2*Lss*KssDT/KssTL"),
         assumption_drug_gg_Ttot  = Dss_thy > 2*Ttotss_thy,
         str_drug_gg_Ttot         = ifelse(assumption_drug_gg_Ttot,"o D > 2*Ttot","x D < 2*Ttot"),
         assumption_drug_gg_Ccrit = Dss_thy > 4*Ccrit_thy,
         str_drug_gg_Ccrit         = ifelse(assumption_drug_gg_Ccrit,"o D > 4*Ccrit","x D < 4*Ccrit"),
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
                                      ifelse(assumption_drug_gg_Ccrit, "D>4*Ccrit; ", "-"),
                                      ifelse(assumption_drug_gg_LssKssDT_KssTL, "D>Lss*KssDT/KssTL ", "-"),
                                      sep = "\n")) %>%
  arrange(assumption_SCIM_list) %>%
  mutate(assumption_SCIM_str = factor(assumption_SCIM_str, levels = unique(assumption_SCIM_str)))

data = data %>% 
  mutate(small_error = ifelse(SCIM_SCIM_ratio > 0.75 & SCIM_SCIM_ratio < 1.25, "yes", "no")) %>%
  filter(error_simulation == 0)


data_all_assum_TRUE = data %>% filter(assumption_all_SCIM == TRUE)
data_drug_gg_Ccrit  = data %>% filter(assumption_drug_gg_Ccrit == TRUE)


# ratio-err SSIM - histogram - grid 4x4 ----

g = ggplot(data, aes(SCIM_SCIM_ratio, fill = small_error))
g = g + facet_grid(str_drug_gg_Ccrit + str_SCIM_lt_30 ~ str_drug_gg_Ttot + str_drug_gg_LssKssDT_KssTL, 
                   scales = "free_y", switch = "y")
g = g + geom_histogram(binwidth = 0.1, center = 1)
#breaks = 10^seq(-3,2,by=1)
#g = g + xgx_scale_x_log10()
g = g + labs(x    = "ASIR_theory/ASIR_simulation",
             y    = "Number of Simulations",
             fill = "Ratio between 75-125%")
#g = g + geom_vline(xintercept = 1, color = "grey30")
#colors = scales::seq_gradient_pal("blue", "red", "Lab")(seq(0,1,length.out=length(unique(data$assumption_SCIM_str))))
g = g + scale_fill_manual(values = c(yes = "grey50", no = "red"))

print(g)
g1 = g
ggsave(width = 8, height= 4, filename = "./figures/Task52c_GlobalSensitivityAnalysis_SSIM_4assum_ratio.png")

# focus only when all assumptions met ----
g = g1 %+% data_all_assum_TRUE
g = g + scale_x_continuous()
print(g)

# focus only on when Dss_thy > 5*Ccrit
g = g1 
g = g + scale_x_continuous()
g = g + facet_wrap(~str_drug_gg_Ccrit)
print(g)

#Dss_thy vs Dss_sim based on Ccrit ----
g = ggplot(data, aes(x = Dss_thy/Ccrit_thy, y = Dss_sim/Dss_thy, color = small_error))
g = g + geom_point(alpha = 0.1)
g = g + geom_vline(xintercept = 4)
g = g + xgx_scale_x_log10(limits = c(1,1000))
g = g + xgx_scale_y_log10(limits = c(.1, 1.2), breaks = seq(0.1, 1, by = 0.1))
g = g + labs(x = "Dss_theory/Ccrit",
             y = "Dss: simulation/theory",
             fill = "ASIR theory/sim ratio\nbetween 75-125%")
g = g + scale_color_manual(values = c(yes = "black", no = "red"))
g = g + guides(colour = guide_legend(override.aes = list(alpha=1)))
print(g)


#Kss_TL vs (Tss*Lss/TLss)-sim ----
# this shows that the high SCIM ratio correspond to places where 
g = ggplot(data_all_assum_TRUE, aes(x = Kss_ratio, y = SCIM_SCIM_ratio, color = small_error))
g = g + geom_point(alpha = .3)
g = g + labs(x = "Kss_TL / (Tss * Lss / TLss)",
             y = "ASIR_theory/ASIR_simulation",
             color = "ASIR theory/sim ratio\nbetween 75-125%")
#g = g + xgx_scale_x_log10()
#g = g + xgx_scale_y_log10()
g = g + scale_color_manual(values = c(yes = "black", no = "red"))

print(g)

stop()




# ratio-err - Parallel Coordinates - all assumptions met ----
# param2uniform = function(x) {(log(x) - log(min(x)))/(log(max(x))-log(min(x)))}
# data_plot = data %>%
#   mutate(kratio = (koff_DT+keDT)/(koff_TL+keTL)) %>%
#   mutate_at(vars(AFIR:keDT,dose_mpk,kratio), funs(tf=param2uniform(.))) %>%
#   select(id, infusion, is_soluble, contains("AFIR"),contains("SCIM"), T0_tf:keDT_tf, dose_mpk_tf, kratio_tf, contains("assumption")) %>%
#   gather(param,param_value,-c(id, infusion, is_soluble, keT_keDT_ratio_tf, kratio_tf, contains("AFIR"), contains("SCIM"), contains("assumption"))) %>%
#   mutate(param = str_replace(param,"_tf","")) %>%
#   mutate(ratio_lt_10pct = ifelse(SCIM_SCIM_ratio < 0.1, "yes","no")) 
# 
# #sort by average param value in one category to help with visualization 
# data_summ = data_plot %>%
#   group_by(param) %>%
#   filter(assumption_all_SCIM == TRUE) %>%
#   summarise(x = mean(param_value)) %>%
#   arrange(x) %>%
#   ungroup()
# kable(data_summ)

# data_plot = data_plot %>%
#   mutate(param = factor(param,  levels = data_summ$param))
# 
# 
# g = ggplot(data_plot, aes(x=param,y=param_value, group = id, color = ratio_lt_10pct, alpha = ratio_lt_10pct))
# g = g + geom_line()
# g = g + facet_grid(infusion~is_soluble)
# g = g + theme(axis.text.x = element_text(angle = 45, hjust = 1))
# g = g + labs(x = "Parameter", y = "Parameter Value")
# g = g + guides(colour = guide_legend(override.aes = list(alpha = 1)))
# g = g + scale_color_manual(values = c(yes = "blue", no = "red"))
# g = g + scale_alpha_manual(values = c(yes = .01, no = 0.01))
# #g = xgx_save(7,7,dirs,"Parallel_Coord","")
# print(g)

#explore individual patients ----
x = data %>%
  filter(assumption_all_SCIM == TRUE & !is.na(SCIM_SCIM_ratio)) %>%
  arrange(desc(SCIM_SCIM_ratio)) %>% 
  mutate(koff_DT = Kd_DT*kon_DT, koff_TL = Kd_TL*kon_TL)

xsub = x %>%
  select(SCIM_SCIM_ratio, SCIM_SCIM_differr, Kss_TL, Kss_DT, koff_TL, koff_DT) %>%
  mutate(koff_ratio = koff_DT/koff_TL)

i = 2
xi = x[i,]
summi = plot_param(xi,model,infusion = xi$infusion)

print(summi$param %>% select(Kd_DT, Kd_TL, kon_TL, kon_DT, Ccrit))
print(summi$compare)

# some notable patients

# DRUG CONC ISN'T HIGH ENOUGH TO BE OVER CCRIT
# This was found by setting i = nrow(x)
# id = 11977 - SCIM_sim is higher, because drug conc wasn't quite linear
#            - Dss_thy = 64, and Ccrit = 32, but then Dss_sim = 32
# id = 36357 - when I requires Dss_thy > 5*Ccrit, then there is about a 20% error in drug concentration and 20% over-estimate of SSIM
#
# KSS approx not accurate
# id = 17015 - 140% diff between SCIM theory and sim




