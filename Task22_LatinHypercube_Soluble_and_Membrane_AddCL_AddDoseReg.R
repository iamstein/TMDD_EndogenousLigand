source("ams_initialize_script.R")
source("SCIM_calculation.R")  
source("ivsc_2cmt_RR_V1.R")
dirs$rscript_name = "Task22_LatinHypercube_Soluble_and_Membrane_AddCL_AddDoseReg.R"
dirs$filename_prefix= str_extract(dirs$rscript_name,"^Task\\d\\d\\w?_")

library(lhs)

n_samples_simulate = 1e4 #number of parameters to simulate for each target type
n_samples_start = n_samples_simulate*10 #number of parameters to start with.  Will select only those at reasonable doses.  

tmax        = 7*52*2 #days (for soluble target, 16 weeks should be long enough)
compartment = 2

model = ivsc_2cmt_RR_KdT0L0()

#read in parameter ranges to explore
param_minmax_in = readxl::read_excel("parameters/Task22_Param_Ranges_AFIR_Tfold.xlsx")

#soluble parameter ranges ----
param_minmax = param_minmax_in %>%
  as.data.frame() %>%
  select(Parameter, min = min_sol,max = max_sol, Units, type_sol) %>%
  mutate(fixed = as.numeric(type_sol)) %>%
  filter(!is.na(fixed))
rownames(param_minmax) = param_minmax$Parameter

n_param   = nrow(param_minmax)

x       = lhs::randomLHS(n_samples_start,n_param)
log_min = matrix(rep(log(param_minmax$min), each = n_samples_start), nrow = n_samples_start, ncol = n_param)
log_max = matrix(rep(log(param_minmax$max), each = n_samples_start), nrow = n_samples_start, ncol = n_param)

param = exp(log_min + (log_max - log_min)*x)
param[is.na(param)] = 0
colnames(param) = param_minmax$Parameter
param_soluble = as.data.frame(param) %>%
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
         tmax      = tmax,
         is_soluble= 1)
cat("instances of zero values for soluble target\n")
print(summarise_all(param_soluble,funs(sum(.==0))))

#membrane-bound parameter ranges ----
param_minmax = param_minmax_in %>%
  as.data.frame() %>%
  select(Parameter, min = min_mem,max = max_mem, Units, type_mem) %>%
  mutate(fixed = as.numeric(type_mem)) %>%
  filter(!is.na(fixed))
rownames(param_minmax) = param_minmax$Parameter

n_param   = nrow(param_minmax)

x       = lhs::randomLHS(n_samples_start,n_param)
log_min = matrix(rep(log(param_minmax$min), each = n_samples_start), nrow = n_samples_start, ncol = n_param)
log_max = matrix(rep(log(param_minmax$max), each = n_samples_start), nrow = n_samples_start, ncol = n_param)

param = exp(log_min + (log_max - log_min)*x)
param[is.na(param)] = 0
colnames(param) = param_minmax$Parameter
param_membrane = as.data.frame(param) %>%
  mutate(keTL      = keL/keL_keTL_ratio,
         keDT      = keT/keT_keDT_ratio,
         Kss_DT    = Kd_DT + keDT/kon_DT,
         Kss_TL    = Kd_TL + keTL/kon_TL,
         TL0       = T0*L0/Kss_TL,
         ksynT     = T0*keT + keTL*TL0,
         Ttotss    = ksynT/keDT,
         Tfold     = Ttotss/T0,
         dose_nmol = Kss_DT*Tfold*CL*tau/AFIR, #choose dose to get AFIR_simple in range
         dose_mpk  = dose_nmol*scale.nmol2mpk,
         keTL      = keL/keL_keTL_ratio,
         tmax      = tmax,
         is_soluble= 0)
cat("instances of zero values for membrane-bound target\n")
print(summarise_all(param_membrane,funs(sum(.==0))))

param_soluble_reduce = param_soluble %>%
  filter(dose_mpk<=100) %>%
  slice(1:n_samples_simulate)

param_membrane_reduce = param_membrane %>%
  filter(dose_mpk<=100) %>%
  slice(1:n_samples_simulate)

param_reduce = bind_rows(param_membrane_reduce, param_soluble_reduce)
param_all = bind_rows(param_membrane, param_soluble)

# look at the dosing ----
g = ggplot(data = param_all, aes(dose_mpk, fill = as.character(is_soluble)))
g = g + geom_histogram()
g = g + xgx_scale_x_log10(breaks = 10^seq(-20,20,by=2))
g = g + geom_vline(aes(xintercept=100),color="black")
g = g + ggtitle(paste0("all parameters, N = (", nrow(param), " )"))
g = g + theme(legend.position = "top")
g1 = g

# create two sets of simulations
# one that will be an infusion and one that will be dosed every 2-4 weeks
param_infusion = param_reduce %>%
  mutate(infusion = TRUE,
         tau      = 21)

param_bolus = param_reduce %>%
  mutate(infusion = FALSE,
         tau      = sample(7*c(2,3,4), nrow(param_reduce), replace = TRUE))

param = bind_rows(param_bolus, param_infusion) %>%
  mutate(keD = CL/V1,
         k12 = Q/V1,
         k21 = Q/V2)

g = ggplot(data = param_reduce,aes(AFIR, fill = as.character(is_soluble)))
g = g + geom_histogram()
g = g + geom_vline(xintercept=c(0.05,0.30),color="black")
g = g + xgx_scale_x_log10(breaks = 10^seq(-20,20,by=2))
g = g + ggtitle(paste0("reduced set, N = (", nrow(param_reduce), " )"))
g = g + theme(legend.position = "top")
g2 = g

gg = gridExtra::arrangeGrob(g1,g2,nrow = 1, ncol = 2)
gridExtra::grid.arrange(gg)



# Run the simulations ----
start_time = Sys.time()
result = list()

#loop through each simulation
n_samples = nrow(param)
for (i in 1:n_samples) { 
  if  ((i %% 100) == 1) {
    cat(paste("run ",i-1," of ", n_samples, "-" , Sys.time(), "\n"))
  }
  param.as.double = param[i,] %>%
    as.numeric()
  names(param.as.double) = names(param)
  dose.nmol       = as.numeric(param.as.double["dose_mpk"])*scale_mpk2nmol
  tau             = param.as.double["tau"]
  
  thy = lumped.parameters.theory    (       param.as.double, dose.nmol,       param$tau[i],              infusion = param$infusion[i])
  sim = lumped.parameters.simulation(model, param.as.double, dose.nmol, tmax, param$tau[i], compartment, infusion = param$infusion[i])

  #identify runs with issues
  if (!is.na(sim$SCIM_sim)) {
    if (sim$TLss_sim < 0 || sim$SCIM_sim < 0) {
      print(paste("there is an issue with run",i)) 
    }
  }
  
  #all parameter values for the output table
  par_in = param.as.double %>% 
    t() %>% 
    as.data.frame()
  name_in = names(par_in)

  #include parameters used in the ODE, but not in the input parameter vector
  par_ode = param.as.double %>%
    model$repar() %>%
    t() %>%
    as.data.frame()
  name_ode = names(par_ode)
  
  par_ode = par_ode[,setdiff(name_ode,name_in)]

  #create result table
  result[[i]] = bind_cols(sim,thy,par_in,par_ode) %>%
    mutate(id   = i, 
           tmax = tmax) %>%
    select(id,everything())
  
  if (((i %% 1e4) == 0) || (i==n_samples)) {
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
    out = plot_param(result[[i]],model, plot_flag = FALSE, infusion = param$infusion[i])
    subtitle = ifelse(result[[i]]$is_soluble==1,"soluble","membrane-bound")
    g = out$plot + labs(subtitle = subtitle)
    print(g)
  }
  
  #if (i==67)
  #  browser()
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
