# ---
# title: "Task16b - Parallel Coordinates - Overiev of simulations"
# output:
#   html_document:
#     toc: true
#     toc_float: true
#     code_folding: hide
# ---

# Setup and Read Data
#```{r, warning=FALSE, message=FALSE}
source("ams_initialize_script.R")
source("SCIM_calculation.R")  
source("ivsc_2cmt_RR_V1.R")
library(RxODE)
model = ivsc_2cmt_RR_KdT0L0()
dirs$rscript_name = "Task20_SimulationProcess_and_Explore.R"
dirs$filename_prefix= str_extract(dirs$rscript_name,"^Task\\d\\d\\w?_")

data_in = read.csv("results/Task19_2019-11-14_20e3.csv",stringsAsFactors = FALSE)
#```

# Compute various quantities for comparing AFIR and SCIM, theory and simulation
### Using adhoc theory calculation for SCIM

#```{r, warning=FALSE, message=FALSE}
#put data into categories ----
data    = data_in %>%
  mutate(SCIM_thy= SCIM_adhoc_thy) %>%
  mutate(Cavgss  = dose_nmol/(CL*tau),
         Kss_TL  = Kd_TL + keTL/kon_TL,
         Kss_DT  = Kd_DT + keDT/kon_DT,
         TL0     = T0*L0/Kss_TL,
         koff_TL = Kd_TL*kon_TL,
         koff_DT = Kd_DT*kon_DT,
         ksynT   = T0*keT + keTL*TL0,
         ksynL   = L0*(kon_TL*T0 + keL) - koff_TL*TL0,
         Lss     = ksynL/keL,
         Lfold   = Lss/L0) %>%
  mutate(AFIR_SCIM_pcterr = abs((AFIR_thy             - SCIM_sim)/SCIM_sim),
         SCIM_SCIM_pcterr = abs((SCIM_Lfold_adhoc_thy - SCIM_sim)/SCIM_sim),
         SCIM_SCIM_pcterr = ifelse(SCIM_SCIM_pcterr < 0.001, 0.0001, SCIM_SCIM_pcterr)) %>%
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
         SCIMthy_SCIMsim_category = paste0(SCIMthy_category, ", ", SCIMsim_category))
data = data %>%
  arrange(AFIR_thy) %>%
  mutate(AFIRthy_category = factor(AFIRthy_category, levels = unique(AFIRthy_category))) %>%
  arrange(AFIR_sim) %>%
  mutate(AFIRsim_category = factor(AFIRsim_category, levels = unique(AFIRsim_category))) %>%
  arrange(SCIM_sim) %>%
  mutate(SCIMsim_category = factor(SCIMsim_category, levels = unique(SCIMsim_category)))

#check the assumptions of the data ----
data = data %>%
  mutate(TLss_thy                 = TL0_thy*SCIM_thy,
         Ttotss                   = T0*Tfold,
         koff_DT                  = Kd_DT*kon_DT,
         V1                       = keD/CL,
         Ccrit                    = ksynT*V1/CL,
         assumption_AFIR_lt_30    = AFIR_thy < 0.30,
         assumption_SCIM_lt_30    = SCIM_Lfold_adhoc_thy < 0.30,
         assumption_drug_gg_T0    = Cavgss > 5*Ttotss,
         assumption_drug_gg_KssDT = Cavgss > 5*Kss_DT,
         assumption_koffDT_gt_keT = koff_DT > keT,
         assumption_koffTL_fast   = koff_TL > 1/30,
         assumption_Cavgss_gt_Ccrit = Cavgss > 5*Ccrit,
         assumption_Cavgss_gg_LssKssDT_KssTL = Cavgss > 5*Kss_DT*Lss/Kss_TL,
#         assumption_T0simple    = T0/(ksynT/keT) > 0.5 & T0/(ksynT/keT) < 2, #the simple formula works for T0
         assumption_L_noaccum   = Lss/L0 <= 1.01, #then SCIM = AFIR
#         assumption_Tss_gt_Lss  = Tss_sim > Lss_sim,
         assumption_all_AFIR    = assumption_AFIR_lt_30 & 
                                  assumption_drug_gg_T0 &
                                  assumption_drug_gg_KssDT &
                                  #assumption_koffDT_gt_keT & 
                                  assumption_koffTL_fast &           
                                  assumption_Cavgss_gg_LssKssDT_KssTL &
                                  assumption_Cavgss_gt_Ccrit &
                                  assumption_L_noaccum,
         assumption_all_SCIM =    assumption_SCIM_lt_30 & 
                                  assumption_drug_gg_T0 &
                                  assumption_drug_gg_KssDT &
                                  #assumption_koffDT_gt_keT & 
                                  assumption_koffTL_fast &    
                                  assumption_Cavgss_gt_Ccrit &
                                  assumption_Cavgss_gg_LssKssDT_KssTL)


data = data %>%
  arrange(id) %>%
  select(id,everything())
filename = paste0("results/",dirs$filename_prefix,"data.csv")
write.csv(data,filename, row.names = FALSE )

assumptions = data %>%
  select(id,AFIR_thy,AFIR_sim,SCIM_simplest_thy,SCIM_adhoc_thy,SCIM_sim,starts_with("assumption")) %>%
  arrange(SCIM_sim)
nam = names(assumptions) %>%
  str_replace("^assumption_","")
names(assumptions) = nam
#View(assumptions)
#```

# Put data into error categories and summarize
#```{r, warning=FALSE, message=FALSE}
threshold = 0.1
data_errss = data_in %>%
  filter(abs(TLss_frac_change)>=threshold) 
print(paste0(nrow(data_errss)," of ", nrow(data_in), " : Number of rows with TLss_frac_change > 0.1"))

data_err0 = data_in %>%
  filter(abs(TL0_05tau_frac_change)>=threshold)
print(paste0(nrow(data_err0)," of ", nrow(data_in), " : Number of rows with TL0_05tau_frac_change > 0.1"))

# error historgram ----
data_quick_summ = data %>%
  select(id,AFIR_thy, SCIM_sim, AFIR_SCIM_pcterr, TLss_frac_change, TL0_05tau_frac_change) %>%
  gather(key,value,-c(id)) %>%
  mutate(category = case_when((value < threshold) ~ "keep_low",
                              ((value >= threshold) & (key %in% c("AFIR_SCIM_pcterr","SCIM_sim"))) ~ "keep_high",
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
  group_by(AFIRthy_SCIMsim_category) %>%
  count() %>%
  arrange(desc(n))
kable(data_summary)
#```

# histogram of AFIR_theory and SCIM_sim error ----
# g = ggplot(data, aes(Lfold))
# g = g + geom_histogram()
# g = g + xgx_scale_x_log10()
# print(g)

g = ggplot(data, aes(AFIR_SCIM_pcterr, fill = assumption_all_AFIR))
g = g + geom_histogram()
g = g + xgx_scale_x_log10()
g = g + scale_fill_manual(values = c(`TRUE`="grey50",`FALSE`="pink"))
g = g + ggtitle("AFIR theory vs SCIM simulation")
print(g)

data_sort_error = data %>% 
  filter(assumption_all_AFIR == TRUE) %>%
  arrange(desc(AFIR_SCIM_pcterr))

#error seems to be due to small accumulation of L.  See plot below (commented)
#plot_param(data_sort_error[1,],model)

# histogram SCIM_theory vs SCIM_sim ----
g = ggplot(data, aes(SCIM_SCIM_pcterr, fill = assumption_all_SCIM))
g = g + geom_histogram()
g = g + xgx_scale_x_log10()
g = g + scale_fill_manual(values = c(`TRUE`="grey50",`FALSE`="pink"))
g = g + labs(x = "Percent Error",
             y = "Number of Simulations")
g = g + ggtitle("SCIM adhoc Lfold theory vs\nSCIM simulation")
print(g)

#patient with biggest error and all assumptions true
data_sort_error = data %>% 
  filter(assumption_all_SCIM == TRUE) %>%
  arrange(desc(SCIM_SCIM_pcterr))

out = plot_param(data_sort_error[1,], model, plot_flag = FALSE)
g = out$plot + labs(subtitle = "assumptions true\nbig diff btw thy and sim")
print(g)

#patient with small error and a false assumption
data_sort_error = data %>% 
  filter(assumption_all_SCIM == FALSE) %>%
  arrange(SCIM_SCIM_pcterr)

id_plot = 7651# data_sort_error$id[1]
out = plot_param(filter(data,id == id_plot), model, plot_flag = FALSE)
g = out$plot + labs(subtitle = "assumptions false\nsmall diff btw thy and sim")
print(g)

assumption_focus = data %>%
  filter(assumption_all_SCIM == FALSE) %>%
  arrange(SCIM_SCIM_pcterr) %>%
  select(id, SCIM_SCIM_pcterr, contains("assumption")) 
#  filter(id %in% c(data_sort_error$id[1:10])) 
nam = names(assumption_focus) %>%
  str_replace("assumption_","")
names(assumption_focus) = nam
View(assumption_focus)
stop()

# AFIRsim vs SCIMsim : 3x3 plot colors ----
#```{r, warning=FALSE, message=FALSE, results = 'asis'}
param2uniform = function(x) {(log(x) - log(min(x)))/(log(max(x))-log(min(x)))}
data_plot = data_keep %>%
  mutate_at(vars(AFIR:kon_TL,dose_mpk), funs(tf=param2uniform(.))) %>%
  select(id,contains("AFIR"),contains("SCIM"), T0_tf:kon_TL_tf, dose_mpk_tf, contains("assumption")) %>%
  gather(param,param_value,-c(id, contains("AFIR"), contains("SCIM"), contains("assumption"))) %>%
  mutate(param = str_replace(param,"_tf",""))

#sort by average param value in one category to help with visualization ----
data_summ = data_plot %>%
  filter(AFIRthy_SCIMsim_category == "AFIRthy < 5%, SCIMsim > 30%") %>%
  group_by(param,AFIRthy_SCIMsim_category) %>%
  summarise(x = mean(param_value)) %>%
  arrange(x) %>%
  ungroup()
kable(data_summ)

data_plot = data_plot %>%
  mutate(param = factor(param, 
                        levels = data_summ$param))

g = ggplot(data_plot, aes(x=param,y=param_value, group = id, color = assumption_all_AFIR, alpha = assumption_all_AFIR))
g = g + geom_line()
g = g + facet_grid(SCIMsim_category~AFIRthy_category,switch = "y")
g = g + theme(axis.text.x = element_text(angle = 45, hjust = 1))
g = g + labs(x = "Parameter", y = "Parameter Value")
g = g + guides(colour = guide_legend(override.aes = list(alpha = 1)))
g = g + scale_color_manual(values = c(`TRUE` = "blue", `FALSE` = "red"))
g = g + scale_alpha_manual(values = c(`TRUE` = .01, `FALSE` = 0.01))
g = xgx_save(7,7,dirs,"Parallel_Coord_Soluble_3x3_AFIRthy_AFIRsim","")
g1 = g
print(g)
#```

# Focus on when assumptions are true and SCIM > 30% while AFIRthy < 5%
#```{r, warning=FALSE, message=FALSE}
#explore data data where all assumptions are true and still
#AFIRsim > 30% and AFIRthy < 5% ---- on look, there is lots of L0!!!
#focus on this plot
data_new = data_plot %>%
  filter(SCIMsim_category == "SCIMsim > 30%",
         AFIRthy_category == "AFIRthy < 5%",
         assumption_all_AFIR == TRUE)

if (nrow(data_new)==0) {
  stop("there are no examples of AFIR_thy<5% and SCIMsim > 30%")
}
  

g = g1
g = g %+% data_new
g = g + geom_line(alpha = 0.05)
g = xgx_save(5,5,dirs,"Parallel_Coord_Soluble_AFIRthy_lt_5_SCIMsim_ge_30","")
g2= g
#print(g)
#```

# Identify patients where all assumptions true and theory vs sim disagree.
#```{r, warning=FALSE}
id = unique(sort(data_new$id))
print("these IDs, even with all the restrictions, AFIR and SCIM still don't match")
print(id)

for (id_plot in id[1:5]) {

  g = g2
  g = g + geom_line(data = filter(data_new,id==id_plot),
                    size = 2,
                    color = "black")
  g = g + ggtitle(paste("id =",id_plot))
  print(g)
  filepref = paste0("Parallel_Coord_Soluble_AFIRthy_lt_5_SCIMsim_ge_30_",id_plot)
  g = xgx_save(5,5,dirs,filepref,"")
  
  #simulate a patient where theory and simulation disagree 
  param = data %>%
    filter(id==id_plot)
  
  assumptions = param %>%
    select(contains("assumption")) %>%
    t() 
  kable(assumptions)
  
  tmax = 365 #days
  tau  = param$tau   #days
  dose_nmol = param$dose_nmol
  compartment = 2
  infusion = TRUE
  
  nam   = names(param)
  param_as_double = param %>%
    as.numeric() %>%
    setNames(nam)
  param_as_double = param_as_double[model$pin]
  
  param_print = param_as_double %>%
    t() %>%
    as.data.frame() %>%
    mutate(CL = signif(keD/V1,2),
           id = id_plot) %>%
    select(id, CL,T0,L0,Kd_DT,Kd_TL,kon_DT,kon_TL,keT,keL,keDT,keTL)
  
  ev = eventTable(amount.units="nmol", time.units="days")
  sample.points = c(seq(0, tmax, 0.1), 10^(-3:0)) # sample time, increment by 0.1
  sample.points = sort(sample.points)
  sample.points = unique(sample.points)
  ev$add.sampling(sample.points)
  if (infusion == FALSE) {
    ev$add.dosing(dose=dose_nmol, start.time = tau, nbr.doses=floor(tmax/tau), dosing.interval=tau, dosing.to=compartment)
  } else {
    ev$add.dosing(dose=dose_nmol, start.time = tau, nbr.doses=floor(tmax/tau)+1, dosing.interval=tau, dosing.to=compartment, dur = tau)
  }  
  
  sim = lumped.parameters.simulation(model, param_as_double, dose_nmol, tmax, tau, compartment, infusion)
  thy = lumped.parameters.theory    (       param_as_double, dose_nmol,       tau,              infusion)
  
  sim_rename = sim
  nam = names(sim_rename) %>%
    str_replace_all("_sim$","")
  names(sim_rename) = nam
  sim_rename$type = "sim"
  
  thy_rename = thy
  nam = names(thy_rename) %>%
    str_replace_all("_thy$","")
  names(thy_rename) = nam
  thy_rename$type = "thy"
  
  compare = bind_rows(sim_rename,thy_rename) %>%
    select(type,Dss,T0,L0,TL0,Ttotss,Lss,TLss,AFIR,SCIM)

  init = model$init(param_as_double)
  out  = model$rxode$solve(model$repar(param_as_double), ev, init)
  out  = model$rxout(out)
  
  out_plot = out %>%
    select(time,D,T,DT,L,TL) %>%
    gather(cmt,value,-time)
  out_last = out_plot[(out$time==max(out$time)),]
  
  g = ggplot(out_plot,aes(x=time,y=value, color = cmt, group= cmt))
  g = g + geom_line()
  g = g + geom_label(data = out_last, aes(label = cmt), show.legend = FALSE, hjust=1)
  g = g + geom_vline(xintercept = tau, linetype = "dotted")
  g = g + xgx_scale_x_time_units(units_dataset = "days", units_plot = "weeks")
  g = g + xgx_scale_y_log10()
  g = g + labs(y = "Concentration (nm)", color = "")
  g = g + ggtitle(paste0(  "id = ",param$id,
                           "\nAFIR_thy  = ",signif(thy$AFIR_thy,2),
                           "\nAFIR_sim  = ",signif(sim$AFIR_sim,2),
                           "\nSCIM_thy = ",signif(thy$SCIM_adhoc_thy,2),
                           "\nSCIM_sim = ",signif(sim$SCIM_sim,2)))
  filepref = paste0("Parallel_Coord_Soluble_AFIRthy_lt_5_SCIMsim_ge_30_",id_plot)
  g = xgx_save(5,5,dirs,filepref,"")
  print(g)
  
  #unfortunately, kable does not work properly inside for loop
  print(t(param_print))
  print(t(compare))
}
#```

