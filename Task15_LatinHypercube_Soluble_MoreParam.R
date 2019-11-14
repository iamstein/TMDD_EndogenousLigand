source("ams_initialize_script.R")
source("SCIM_calculation.R")  
source("ivsc_2cmt_RR_V1.R")
dirs$rscript_name = "Task15_LatinHypercube_Soluble_MoreParam.R"
dirs$filename_prefix= str_extract(dirs$rscript_name,"^Task\\d\\d\\w?_")

library(lhs)

n_samples = 1e4 #number of parameters sets to simulate

tmax        = 7*20 #days (for soluble target, 16 weeks should be long enough)
compartment = 2
infusion    = TRUE

model = ivsc_2cmt_RR_KdT0L0()

#read in parameter ranges to explore
param_minmax.in = readxl::read_excel("parameters/Task15_Param_Ranges_AFIR_Tfold.xlsx")
param_minmax = param_minmax.in %>%
  as.data.frame() %>%
  select(Parameter, min = min_sol,max = max_sol, units, type_sol) %>%
  mutate(fixed = as.numeric(type_sol)) %>%
  filter(!is.na(fixed))
rownames(param_minmax) = param_minmax$Parameter

#get a latin hyper cube of random variables
n_param   = nrow(param_minmax)

x       = lhs::randomLHS(n_samples,n_param)
log_min = matrix(rep(log(param_minmax$min), each = n_samples), nrow = n_samples, ncol = n_param)
log_max = matrix(rep(log(param_minmax$max), each = n_samples), nrow = n_samples, ncol = n_param)

param = exp(log_min + (log_max - log_min)*x)
param[is.na(param)] = 0
colnames(param) = param_minmax$Parameter
param = as.data.frame(param) %>%
  mutate(keTL      = keL/keL_keTL_ratio,
         Kss_DT    = Kd_DT + keDT/kon_DT,
         Kss_TL    = Kd_TL + keTL/kon_TL,
         TL0       = T0*L0/Kss_TL,
         ksynT     = T0*keT + keTL*TL0,
         Ttotss    = ksynT/keDT,
         Tfold     = Ttotss/T0,
         dose_nmol = Kss_DT*Tfold*CL*tau/AFIR, #choose dose to get AFIR_simple in range
         dose_mpk  = dose_nmol*scale.nmol2mpk,
         keTL      = keL/keL_keTL_ratio,
         tmax      = tmax)
cat("instances of zero values\n")
print(summarise_all(param,funs(sum(.==0))))

#look at the dosing ----
g = ggplot(data = param,aes(dose_mpk))
g = g + geom_histogram()
g = g + xgx_scale_x_log10(breaks = 10^seq(-20,20,by=2))
g = g + geom_vline(aes(xintercept=100),color="red")
g = g + ggtitle(paste0("all parameters, N = (", nrow(param), " )"))
g1 = g

param_reduce = filter(param,dose_mpk<=100)
g = ggplot(data = param_reduce,aes(AFIR))
g = g + geom_histogram()
g = g + geom_vline(xintercept=c(0.05,0.30),color="red")
g = g + xgx_scale_x_log10(breaks = 10^seq(-20,20,by=2))
g = g + ggtitle(paste0("reduced set, N = (", nrow(param_reduce), " )"))
g2 = g

gg = gridExtra::arrangeGrob(g1,g2,nrow = 1, ncol = 2)
gridExtra::grid.arrange(gg)

# Run the simulations ----
start_time = Sys.time()
result = list()

#loop through
for (i in 1:n_samples) { 
  if  ((i %% 100) == 1) {
    cat(paste("run ",i-1," of ", n_samples, "-" , Sys.time(), "\n"))
  }
  param.as.double = param[i,] %>%
    as.numeric()
  names(param.as.double) = names(param)
  dose.nmol       = as.numeric(param.as.double["dose_mpk"])*scale_mpk2nmol
  tau             = param.as.double["tau"]
  
  error_flag = 0
  thy = lumped.parameters.theory    (       param.as.double, dose.nmol,       tau,              infusion = infusion)
  sim = lumped.parameters.simulation(model, param.as.double, dose.nmol, tmax, tau, compartment, infusion = infusion)

  #all parameter values for the output table
  par = param.as.double %>% 
    t() %>% 
    as.data.frame() %>%
    mutate(id = i,
           tmax)

  #create result table
  result[[i]] = bind_cols(sim,thy,par) %>%
    select(id,everything())
  
  if (((i %% 1e5) == 0) || (i==n_samples)) {
    filename = paste0("results/",dirs$filename_prefix,Sys.Date(),"_",i/1e3,"e3.csv")
    results_save = result %>%
      bind_rows()
    
    #want to reduce the number of sig digs to save a bit of space, 
    #but need to keep them for the id
    #id_all_sig_digs = results_save$id
    #results_save = results_save %>%
    #  signif(digits = 6) %>%
    #  mutate(id = id_all_sig_digs)
    
    write.csv(results_save,filename,quote = FALSE, row.names = FALSE)
  }
  
  #plot a simulation after every n_sim simulations
  n_sim = 200
  if ( ((i %% n_sim)==1) & (result[[i]]$error_simulation == FALSE) ) {
    plot_param(result[[i]],model)
  }
}
stop_time = Sys.time()
cat("Total time: total_duration\n")
total_duration = (stop_time-start_time)
print(total_duration)

duration_per_run_sec = total_duration/n_samples
cat("Time per run:")
cat(paste0(signif(as.numeric(total_duration/n_samples, units = "secs"),2), " sec\n"))

#check the initial condition and steady state ----
check = results_save %>%
  select(TLss_frac_change, TL0_05tau_frac_change) %>%
  gather() %>%
  mutate(value      = abs(value))

g = ggplot(check,aes(value))
g = g + geom_histogram()
g = g + facet_wrap(~key, switch = "y", scales = "free_x")
g = g + xgx_scale_x_log10()
print(g)
