source("ams_initialize_script.R")
source("SCIM_calculation.R")  
source("ivsc_2cmt_RR_V1.R")
dirs$rscript_name = "Task13b_LHS_Soluble_DaVinci.R"
dirs$filename_prefix= str_extract(dirs$rscript_name,"^Task\\d\\d\\w?_")

library(lhs)

n_samples = 1e6 #number of parameters sets to simulate

model = ivsc_2cmt_RR_KdT0L0()

#read in parameter ranges to explore
param_minmax.in = readxl::read_excel("parameters/Task13_Param_Ranges_AFIR_Tfold.xlsx")
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
  mutate(Kss_DT    = Kd_DT + keDT/kon_DT,
         Tfold     = keT/keDT,
         dose_nmol = Kss_DT*Tfold*CL*tau/AFIR,
         dose_mpk  = dose_nmol*scale.nmol2mpk,
         keTL      = keL/Lfold)
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
tmax = 16*7 #days (for soluble target, 16 weeks should be long enough)
tau  = param$tau[1]   #days
compartment = 2
infusion    = TRUE

start_time = Sys.time()
result = list()

#start at 30076 because there is an error and I want to confirm error checking works
for (i in c(30076,1:n_samples)) { 
  if  ((i %% 100) == 1) {
    cat(paste("run ",i-1," of ", n_samples, "-" , Sys.time(), "\n"))
  }
  param.as.double = param[i,] %>%
    as.numeric()
  names(param.as.double) = names(param)
  dose.nmol       = as.numeric(param.as.double["dose_mpk"])*scale_mpk2nmol
  
  
  error_flag = 0
  thy = lumped.parameters.theory    (       param.as.double, dose.nmol,       tau,              infusion = infusion)
  sim = lumped.parameters.simulation(model, param.as.double, dose.nmol, tmax, tau, compartment, infusion = infusion)

  #all parameter values for the output table
  par = param.as.double %>% 
    t() %>% 
    as.data.frame() %>%
    mutate(id = i)

  #create result table
  result[[i]] = bind_cols(sim,thy,par) %>%
    mutate(Cavgss = dose_nmol*scale.mpk2nmol/(CL*tau)) %>%
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
  
  #plot a simulation after every n_sim simulations
  n_sim = 500
  if ( ((i %% n_sim)==0) & (result[[i]]$error_simulation == FALSE) ) {
    ev = eventTable(amount.units="nmol", time.units="days")
    sample.points = c(seq(0, tmax, 0.1), 10^(-3:0)) # sample time, increment by 0.1
    sample.points = sort(sample.points)
    sample.points = unique(sample.points)
    ev$add.sampling(sample.points)
    
    #add dur  tau for a long infusion
    if (infusion == FALSE) {
      ev$add.dosing(dose=dose.nmol, start.time = tau, nbr.doses=floor(tmax/tau), dosing.interval=tau, dosing.to=compartment)
    } else {
      ev$add.dosing(dose=dose.nmol, start.time = tau, nbr.doses=floor(tmax/tau)+1, dosing.interval=tau, dosing.to=compartment, dur = tau)
    }  
    
    init = model$init(param.as.double)
    out  = model$rxode$solve(model$repar(param.as.double), ev, init)
    out  = model$rxout(out)
    
    out_plot = out %>%
      select(time,D,T,DT,L,TL) %>%
      gather(cmt,value,-time)
    out_last = out[length(out),]
    
    g = ggplot(out_plot,aes(x=time,y=value, color = cmt, group= cmt))
    g = g + geom_line()
    g = g + geom_label(data = out_last, aes(label = cmt))
    g = g + geom_vline(xintercept = tau, linetype = "dotted")
    g = g + xgx_scale_x_time_units(units_dataset = "days", units_plot = "weeks")
    g = g + xgx_scale_y_log10()
    g = g + ggtitle(paste0("run ", i, 
                           "\nAFIR_thy  = ",signif(result[[i]]$AFIR_thy,2),
                           "\nAFIR_sim  = ",signif(result[[i]]$AFIR_sim,2),
                           "\nSCIM_sim = ",signif(result[[i]]$SCIM_sim,2)))
    print(g)
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
