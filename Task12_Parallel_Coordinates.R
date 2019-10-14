source("ams_initialize_script.R")
source("SCIM_calculation.R")  
source("ivsc_2cmt_RR_V1.R")
dirs$rscript_name = "Task12_Parallel_Coordinates.R"
dirs$filename_prefix= str_extract(dirs$rscript_name,"^Task\\d\\d\\w?_")


data_in = read.csv("results/Task11_2019-10-09_1e3.csv",stringsAsFactors = FALSE)
data    = data_in %>%
  mutate(AFIR_SCIM_sqerr = (AFIR_thy - SCIM_sim)^2) %>%
  mutate(AFIRthy_category = case_when(AFIR_thy <  0.05 ~ "AFIRthy < 5%",
                                      AFIR_thy >  0.30 ~ "AFIRthy > 30%",
                                      AFIR_thy >= 0.05 &  AFIR_thy <= 0.30 ~ "5% <= AFIRthy <= 30%"),
         SCIMsim_category = case_when(SCIM_sim <  0.05 ~ "SCIMsim < 5%",
                                      SCIM_sim >  0.30 ~ "SCIMsim > 30%",
                                      SCIM_sim >= 0.05 &  SCIM_sim  <= 0.30 ~ "5% <= SCIMthy <= 30%"),
         AFIR_SCIM_category = paste0(AFIRthy_category, ", ", SCIMsim_category),
         error_category = case_when(AFIR_SCIM_sqerr < 0.1 ~ "low_error",
                                    TRUE ~ "high_error"))
data = data %>%
  arrange(AFIR_thy) %>%
  mutate(AFIRthy_category = factor(AFIRthy_category, levels = unique(AFIRthy_category))) %>%
  arrange(SCIM_sim) %>%
  mutate(SCIMsim_category = factor(SCIMsim_category, levels = unique(SCIMsim_category)))

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
                              ((value >= threshold) & (key %in% c("AFIR_thy"))) ~ "remove_high_AFIR",
                              TRUE ~ "remove_high_error"))

g = ggplot(data_quick_summ, aes(value, fill = category))
g = g + geom_histogram()
g = g + facet_wrap(~key, scales = "free")
g = g + scale_fill_manual(values = c(keep_low = "grey80", 
                                     keep_high = "grey50",
                                     remove_high_error = "red", 
                                     remove_high_AFIR = "blue"))
g = g + xgx_scale_x_log10()
g = g + ggtitle("")
print(g)

#keep only the simulations with no issues 
data_keep = data %>%
  filter(TLss_frac_change < threshold, 
         TL0_05tau_frac_change < threshold)

#put simulations into different categories
data_summary = data_keep %>%
  group_by(AFIR_SCIM_category) %>%
  count() %>%
  arrange(desc(n))
kable(data_summary)

#plot results ----

param2uniform = function(x) {(log(x) - log(min(x)))/(log(max(x))-log(min(x)))}
data_plot = data_keep %>%
  mutate()
  mutate_at(vars(dose:kon_TL), funs(tf=param2uniform(.))) %>%
  select(id,AFIRthy_category, SCIMsim_category, AFIR_SCIM_sqerr, AFIR_SCIM_category, dose_tf:kon_TL_tf) %>%
  gather(param,param_value,-c(id, AFIRthy_category, SCIMsim_category, AFIR_SCIM_sqerr,AFIR_SCIM_category)) %>%
  mutate(param = str_replace(param,"_tf","")) %>%
  mutate(param = factor(param, 
                        levels = c("dose","T0","Kd_DT","kon_DT","keT","keDT","L0","Kd_TL","kon_TL","keL","keTL")))

g = ggplot(data_plot, aes(x=param,y=param_value, group = id))
g = g + geom_line(alpha = .1)
#g = g + scale_color_gradient(low = "blue", high = "red", breaks = c(0,.25,.5,.75,1), limits = c(0,1))
g = g + facet_grid(SCIMsim_category~AFIRthy_category,switch = "y")
g = g + theme(axis.text.x = element_text(angle = 45))
#g = g + scale_color_manual(values = c(low = "black", high = "red"))
#g = g + scale_alpha_manual(values = c(low = .01,    high = 0.00))
print(g)

