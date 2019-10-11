source("ams_initialize_script.R")
source("SCIM_calculation.R")  
source("ivsc_2cmt_RR_V1.R")
dirs$rscript_name = "Task11_Latin_Hyper_Cube_Sample.R"
dirs$filename_prefix= str_extract(dirs$rscript_name,"^Task\\d\\d\\w?_")

library(lhs)

model = ivsc_2cmt_RR_KdT0L0()

#read in parameter ranges to explore
param_minmax.in = readxl::read_excel("parameters/Task10_Param_Ranges.xlsx")
param_minmax = param_minmax.in %>%
  as.data.frame() %>%
  select(Parameter,min,max,units, fixed) %>%
  mutate(fixed = as.numeric(fixed)) %>%
  filter(!is.na(fixed))
rownames(param_minmax) = param_minmax$Parameter

#get a latin hyper cube of random variables
n_samples = 1e2
n_param   = nrow(param_minmax)

x       = lhs::randomLHS(n_samples,n_param)
log_min = matrix(rep(log(param_minmax$min), each = n_samples), nrow = n_samples, ncol = n_param)
log_max = matrix(rep(log(param_minmax$max), each = n_samples), nrow = n_samples, ncol = n_param)

param = exp(log_min + (log_max - log_min)*x)
param[is.na(param)] = 0
colnames(param) = param_minmax$Parameter
param = as.data.frame(param)
cat("instances of zero values\n")
print(summarise_all(param,funs(sum(.==0))))

# Dose time, frequency, compartment, nominal dose
tmax = 52*7 #days
tau  = 21   #days
compartment = 2

start_time = Sys.time()
result = list()
for (i in 1:n_samples) {
  if  ((i %% 100) == 1) {
    cat(paste("run ",i-1," of ", n_samples, "-" , Sys.time(), "\n"))
  }
  param.as.double = param[i,] %>%
    as.numeric()
  names(param.as.double) = names(param)
  dose.nmol       = as.numeric(param.as.double["dose"])*scale_mpk2nmol

  sim = lumped.parameters.simulation(model, param.as.double, dose.nmol, tmax, tau, compartment)
  thy = lumped.parameters.theory    (       param.as.double, dose.nmol,       tau)
  
  #all parameter values for the output table
  par = param.as.double %>% 
    t() %>% 
    as.data.frame() %>%
    mutate(id = i)

  #create result table
  result[[i]] = bind_cols(sim,thy,par) %>%
    mutate(Cavgss = dose*scale.mpk2nmol/(CL*tau)) %>%
    select(-c(TL0_sim, T0_sim, L0_sim, Ttotss_sim, L_sim , D_sim,  time_last_dose,
                       T0_thy,         Ttotss_thy, Lss_thy,Dss_thy, 
                       F, ka, CL, Q, V1, V2, keD, k12, k21, Vm, Km, tau)) %>%
    select(id,everything())
    
  if (((i %% 1e5) == 0) || (i==n_samples)) {
    filename = paste0("results/",dirs$filename_prefix,Sys.Date(),"_",i/1e3,"e3.csv")
    results_save = result %>%
      bind_rows() %>%
      signif(digits = 3)
    write.csv(results_save,filename,quote = FALSE, row.names = FALSE)
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
