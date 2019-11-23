source("ams_initialize_script.R")
source("SCIM_calculation.R")  
source("ivsc_2cmt_RR_V1.R")
dirs$rscript_name = "Task14c_Parallel_Coordinates_Soluble_2019-10-12_AFIR.R"
dirs$filename_prefix= str_extract(dirs$rscript_name,"^Task\\d\\d\\w?_")


data_in = read.csv("results/Task13_2019-10-12_78.771e3.csv",stringsAsFactors = FALSE)

#put data into categories ----
data    = data_in %>%
  filter(id > 3e4 & id < 4e4) %>%
  mutate(Kss_TL  = Kd_TL + keTL/kon_TL,
         Kss_DT  = Kd_DT + keDT/kon_DT,
         TL0     = T0*L0/Kss_TL,
         koff_TL = Kd_TL*kon_TL,
         koff_DT = Kd_DT*kon_DT,
         ksynT   = T0*keT + keTL*TL0,
         ksynL   = L0*(kon_TL*T0 + keL) - koff_TL*TL0,
         Lss     = ksynL/keL) %>%
  mutate(AFIR_SCIM_sqerr = (AFIR_thy - SCIM_sim)^2) %>%
  mutate(AFIRthy_category = case_when(AFIR_thy <  0.05 ~ "AFIRthy < 5%",
                                      AFIR_thy >  0.30 ~ "AFIRthy > 30%",
                                      AFIR_thy >= 0.05 &  AFIR_thy <= 0.30 ~ "5% <= AFIRthy <= 30%"),
         AFIRsim_category = case_when(AFIR_sim <  0.05 ~ "AFIRsim < 5%",
                                      AFIR_sim >  0.30 ~ "AFIRsim > 30%",
                                      AFIR_sim >= 0.05 &  AFIR_sim  <= 0.30 ~ "5% <= AFIRsim <= 30%"),
         SCIMthy_category = case_when(SCIM_thy <  0.05 ~ "SCIMthy < 5%",
                                      SCIM_thy >  0.30 ~ "SCIMthy > 30%",
                                      SCIM_thy >= 0.05 &  SCIM_thy  <= 0.30 ~ "5% <= SCIMthy <= 30%"),
         SCIMsim_category = case_when(SCIM_sim <  0.05 ~ "SCIMsim < 5%",
                                      SCIM_sim >  0.30 ~ "SCIMsim > 30%",
                                      SCIM_sim >= 0.05 &  SCIM_sim  <= 0.30 ~ "5% <= SCIMsim <= 30%"),
         AFIRthy_AFIRsim_category = paste0(AFIRthy_category, ", ", AFIRsim_category),
         AFIRthy_SCIMsim_category = paste0(AFIRthy_category, ", ", SCIMsim_category),
         AFIRsim_SCIMsim_category = paste0(AFIRsim_category, ", ", SCIMsim_category),
         SCIMthy_SCIMsim_category = paste0(SCIMthy_category, ", ", SCIMsim_category),
         error_category = case_when(AFIR_SCIM_sqerr < 0.1 ~ "low_error",
                                    TRUE ~ "high_error"))
data = data %>%
  arrange(AFIR_thy) %>%
  mutate(AFIRthy_category = factor(AFIRthy_category, levels = unique(AFIRthy_category))) %>%
  arrange(AFIR_sim) %>%
  mutate(AFIRsim_category = factor(AFIRsim_category, levels = unique(AFIRsim_category))) %>%
  arrange(SCIM_sim) %>%
  mutate(SCIMsim_category = factor(SCIMsim_category, levels = unique(SCIMsim_category)))

#check the assumptions of the data ----
data = data %>%
  mutate(Ttotss                 = T0*Tfold,
         koff_DT                = Kd_DT*kon_DT,
         assumption_plateau     = AFIR_thy < 0.30,
         assumption_drug_gg_T0  = Cavgss > 10*Ttotss,
         assumption_drug_gg_Kss = Cavgss > 10*Kss_DT,
         assumption_koffDT_gt_keT = koff_DT > keT,
         assumption_koffTL_fast   = koff_TL > 1/30,
         assumption_CavgssKssTL_gg_LssKssDT = Cavgss*Kss_TL > 10*Kss_DT*Lss,
         assumption_all         = assumption_plateau & 
                                  assumption_drug_gg_T0 &
                                  assumption_drug_gg_Kss &
                                  assumption_koffDT_gt_keT & 
                                  assumption_koffTL_fast &           
                                  assumption_CavgssKssTL_gg_LssKssDT,
         assumption_category    = "",
         assumption_category    = ifelse(assumption_plateau,"","!plateau, "),
         assumption_category    = ifelse(assumption_drug_gg_T0,"","drug !>> T0"))


#put data into error categories ----
threshold = 0.1
data_errss = data_in %>%
  filter(abs(TLss_frac_change)>=threshold) 
print(paste0(nrow(data_errss)," of ", nrow(data_in), " : Number of rows with TLss_frac_change > 0.1"))

data_err0 = data_in %>%
  filter(abs(TL0_05tau_frac_change)>=threshold)
print(paste0(nrow(data_err0)," of ", nrow(data_in), " : Number of rows with TL0_05tau_frac_change > 0.1"))

# error historgram ----
data_quick_summ = data %>%
  select(id,AFIR_thy, SCIM_sim, AFIR_SCIM_sqerr, TLss_frac_change, TL0_05tau_frac_change) %>%
  gather(key,value,-c(id)) %>%
  mutate(category = case_when((value < threshold) ~ "keep_low",
                              ((value >= threshold) & (key %in% c("AFIR_SCIM_sqerr","SCIM_sim"))) ~ "keep_high",
                              ((value >= threshold) & (key %in% c("AFIR_thy"))) ~ "keep_high_AFIR",
                              TRUE ~ "remove_high_error"))

g = ggplot(data_quick_summ, aes(value, fill = category))
g = g + geom_histogram()
g = g + facet_wrap(~key, scales = "free")
g = g + scale_fill_manual(values = c(keep_low = "grey80", 
                                     keep_high = "grey50",
                                     remove_high_error = "red", 
                                     keep_high_AFIR = "blue"))
g = g + xgx_scale_x_log10()
g = g + ggtitle("")
print(g)

#keep only the simulations with no issues 
data_keep = data %>%
  filter(TLss_frac_change < threshold, 
         TL0_05tau_frac_change < threshold)

#put simulations into different categories
data_summary = data_keep %>%
  group_by(AFIRthy_AFIRsim_category) %>%
  count() %>%
  arrange(desc(n))
kable(data_summary)

#1. AFIRsim vs SCIMsim : 3x3 plot colors ----
param2uniform = function(x) {(log(x) - log(min(x)))/(log(max(x))-log(min(x)))}
data_plot = data_keep %>%
  mutate_at(vars(AFIR:kon_TL,dose_mpk), funs(tf=param2uniform(.))) %>%
  select(id,contains("AFIR"),contains("SCIM"), T0_tf:kon_TL_tf, dose_mpk_tf, contains("assumption")) %>%
  gather(param,param_value,-c(id, contains("AFIR"), contains("SCIM"), contains("assumption"))) %>%
  mutate(param = str_replace(param,"_tf",""))

#sort by average param value in one category to help with visualization ----
data_summ = data_plot %>%
  filter(AFIRthy_AFIRsim_category == "AFIRthy < 5%, AFIRsim > 30%") %>%
  group_by(param,AFIRthy_AFIRsim_category) %>%
  summarise(x = mean(param_value)) %>%
  arrange(x) %>%
  ungroup()
print(data_summ)

data_plot = data_plot %>%
  mutate(param = factor(param, 
                        levels = data_summ$param))


g = ggplot(data_plot, aes(x=param,y=param_value, group = id, color = assumption_all))
g = g + geom_line(alpha = 0.01)
g = g + facet_grid(AFIRsim_category~AFIRthy_category,switch = "y")
g = g + theme(axis.text.x = element_text(angle = 45, hjust = 1))
g = g + labs(x = "Parameter", y = "Parameter Value")
g = g + guides(colour = guide_legend(override.aes = list(alpha = 1)))
g = g + scale_color_manual(values = c(`TRUE` = "blue", `FALSE` = "red"))
g = xgx_save(7,7,dirs,"Parallel_Coord_Soluble_3x3_AFIRthy_AFIRsim","")
g1 = g
print(g)

#explore data data where all assumptions are true and still
#AFIRsim > 30% and AFIRthy < 5% ---- on look, there is lots of L0!!!
#focus on this plot
data_new = data_plot %>%
  filter(AFIRsim_category == "AFIRsim > 30%",
         AFIRthy_category == "AFIRthy < 5%",
         assumption_all == TRUE)
g = g1
g = g %+% data_new
g = g + geom_line(alpha = 0.05)
g = xgx_save(5,5,dirs,"Parallel_Coord_Soluble_AFIRthy_lt_5_AFIRsim_ge_30","")
print(g)


stop()
#2. AFIRthy vs AFIRsim : 2 colors ----
data_plot_color = data_plot %>%
  filter(AFIRthy_AFIRsim_category %in% c("AFIRthy < 5%, AFIRsim < 5%","AFIRthy < 5%, AFIRsim > 30%"))

g = ggplot(data_plot_color, aes(x=param, y=param_value, group = id, 
                                color = AFIRthy_AFIRsim_category,
                                alpha = AFIRthy_AFIRsim_category))
g = g + geom_line()
g = g + geom_point()
g = g + theme(axis.text.x = element_text(angle = 45, hjust = 1))
g = g + scale_color_manual(values = c("grey50","red"))
g = g + scale_alpha_manual(values = c(0.01, .1))
g = g + theme(legend.position = "top", legend.direction = "vertical")
g = g + guides(colour = guide_legend(override.aes = list(alpha = 1)))
g = g + labs(x = "Parameter", y = "Parameter Value")
g = xgx_save(4,4,dirs,"Parallel_Coord_Soluble_2cat_AFIRthy_AFIRsim","")
print(g)

#3. SCIMthy vs SCIMsim ----
data_plot_color = data_plot %>%
  filter(SCIMthy_SCIMsim_category %in% c("SCIMthy < 5%, SCIMsim < 5%","SCIMthy < 5%, SCIMsim > 30%"))

g = ggplot(data_plot_color, aes(x=param, y=param_value, group = id, 
                                color = AFIRsim_SCIMsim_category,
                                alpha = AFIRsim_SCIMsim_category))
g = g + geom_line()
g = g + geom_point()
g = g + theme(axis.text.x = element_text(angle = 45, hjust = 1))
g = g + scale_color_manual(values = c("grey50","red"))
g = g + scale_alpha_manual(values = c(0.01, .1))
g = g + theme(legend.position = "top", legend.direction = "vertical")
g = g + guides(colour = guide_legend(override.aes = list(alpha = 1)))
g = g + labs(x = "Parameter", y = "Parameter Value")
g = xgx_save(4,4,dirs,"Parallel_Coord_Soluble_2cat_AFIRsim_SCIMsim","")