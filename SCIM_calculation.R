# Theory #----------------------------------------------------------------------------------
lumped.parameters.theory = function(param.as.double = param.as.double,
                                    dose.nmol       = dose.nmol,
                                    tau             = tau, 
                                    infusion        = FALSE){
  # Arguments:
  #   params_file_path: full path of the parameters file.
  #   dose.nmol: dosing amout in nmol
  #   tau: dosing interval in days
  # Return:
  #   A data frame of lumped parameters calculated from theory
  
  pars    = as.data.frame(t(model$repar(param.as.double)))
  Lss     = with(pars,ksynL/keL)
  Ttotss  = with(pars,ksynT/keDT)

  Kss_TL = with(pars,(koff_TL + keTL)/kon_TL)
  Kss_DT = with(pars,(koff_DT + keDT)/kon_DT)

  
  #compute Ctrough
  dose = dose.nmol
  
  if (infusion==FALSE) {
    alpha= with(pars, .5*(k12+k21+keD + sqrt( (k12+k21+keD)^2 - 4*keD*k21)))
    beta = with(pars, .5*(k12+k21+keD - sqrt( (k12+k21+keD)^2 - 4*keD*k21)))
    A    = with(pars, dose*(k21-alpha)/(V1*(beta-alpha)))
    B    = with(pars, dose*(k21-beta) /(V1*(alpha-beta)))
    Dss = (A*exp(-alpha*tau)/(1-exp(-alpha*tau)) + B*exp(-beta*tau)/(1-exp(-beta*tau)))
  } else {
    CL  = with(pars, keD*V1)
    Dss = dose.nmol/(CL*tau)
  }
  
  #TL0_pos <- ((-b) + sqrt((b^2)-4*a*c))/(2*a) 
  if (pars$keTL == 0) {
    T0    = with(pars,ksynT/keT)
    L0    = with(pars,ksynL/keL)
    
    TL0   = with(pars, ksynL*T0/(koff_TL/kon_TL*keL)) # membrane bound keTL=0
    
    Tfold = Ttotss/T0
    Lfold = Lss/L0
    
    SCIM         = Tfold*Kss_TL/Lss * 1/(Kss_TL/Lss*(Dss/Kss_DT + 1) + 1)
    SCIM_adhoc   = Tfold*Kss_TL/Lss * 1/(Kss_TL/Lss*(Dss/Kss_DT + 1) + Tfold)
    SCIM_simpler = Tfold*Kss_TL/Lss * 1/(Kss_TL/Lss*(Dss/Kss_DT)     + 1)
    SCIM_simplest= Tfold*Kss_TL/Lss * 1/(Kss_TL/Lss*(Dss/Kss_DT)        )
    SCIM_Lfold       = with(pars, Kss_DT*Tfold*Lfold/Dss)
    SCIM_Lfold_adhoc = with(pars,(Kss_DT*Tfold*Lfold)/(Dss+Kss_DT*Tfold*Lfold))
    
    #SCIM = with(pars, 1/(koff_TL/kon_TL * Dss /(Lss*((koff_DT+keDT)/kon_DT)) + koff_TL/kon_TL/Lss + 1)*(ksynT/(keDT*TL0)))

    #SCIM          = with(pars, 1/(Kss_TL * Dss /(Lss*Kss_DT + Kss_TL/Lss + 1)*Ttotss/TL0))
    #SCIM_simpler  = with(pars, 1/(Kss_TL * Dss /(Lss*Kss_DT + Kss_TL/Lss + 1)*Ttotss/TL0))
    #SCIM_simplest = with(pars, 1/(Kss_TL * Dss /(Lss*Kss_DT + Kss_TL/Lss + 1)*Ttotss/TL0))
    
  } else {
    a = with(pars,keTL^2)
    b = with(pars,-(keTL) * (ksynT +ksynL) - (((koff_TL+keTL)/kon_TL) * keT *keL))
    c = with(pars, ksynL*ksynT)

    TL0   = ((-b) -sqrt((b^2)-4*a*c))/(2*a)
    T0    = with(pars,(ksynT - keTL*TL0)/keT)
    L0    = with(pars,(ksynL - keTL*TL0)/keL)
    
    Tfold = Ttotss/T0 #- this is not quite correct.  See simulation 1424 from Task15/16
    Lfold = Lss/L0
    #Tfold = with(pars,keT/keDT)
    
    SCIM             = with(pars,Ttotss/TL0 * 1/(Kss_TL/Lss*(Dss/Kss_DT + 1) + 1))
    SCIM_adhoc       = with(pars,Ttotss/TL0 * 1/(Kss_TL/Lss*(Dss/Kss_DT + 1) + Ttotss/TL0))
    SCIM_simpler     = with(pars,Ttotss/TL0 * 1/(Kss_TL/Lss*(Dss/Kss_DT)     + 1))
    SCIM_simplest    = with(pars,Ttotss/TL0 * 1/(Kss_TL/Lss*(Dss/Kss_DT)        ))
    SCIM_Lfold       = with(pars, Kss_DT*Tfold*Lfold/Dss)
    SCIM_Lfold_adhoc = with(pars,(Kss_DT*Tfold*Lfold)/(Dss+Kss_DT*Tfold*Lfold))
    
    
    #SCIM          = with(pars,Ttotss/((((Kss_TL*Dss*keL)/(Kss_DT*ksynL))+((Kss_TL*keL)/(ksynL))+1)*TL0))
    #SCIM_thy_simpler  = with(pars,Ttot/((((KssTL*keL/ksynL)*(Dss/KssDT)) +1) *TL0_neg) )
    #SCIM_thy_simplest = with(pars,  Ttot/KssTL * Lss/TL0_neg * KssDT/Dss)
    
    #SCIM_thy = with(pars,Ttotss/TL0 * 1/(Kss_TL/Kss_DT*Dss/Lss + Kss_TL/Lss + 1))
    
  }
  
  
  #SCIM_thy_ketl_pos = with(pars,Ttot/((((KssTL*Dss*keL)/(KssDT*ksynL))+((KssTL*keL)/(ksynL))+1)*TL0_pos))
  #SCIM_thy_ketl_neg_31 = with(pars, (ksynT*ksynL*KssDT) / (KssTL*keL*keDT*Dss*TL0_neg) ) #same as simplest
  
  AFIR        = with(pars,(Kss_DT*Tfold)/(Dss+Tfold*Kss_DT))
  AFIR_simple = with(pars,(Kss_DT*Tfold)/Dss)
  
  kon_TL = with(pars, koff_TL/Kd_TL)
  lumped_parameters_theory = data.frame(
    T0_thy = T0,
    TL0_thy = TL0,
    L0_thy  = L0,
    Ttotss_thy = Ttotss,
    Lss_thy = Lss,
    Tss_thy  = -99,
    TLss_thy = SCIM*TL0,
    Tfold_thy = Ttotss/T0,
    Lfold_thy = Lss/L0,
    TssLss_TLss_thy = with(pars, Kd_TL + keTL/kon_TL),

    Dss_thy = Dss,
    Cavgss_thy = Dss,
    Ccrit_thy  = with(pars, ksynT/keD),
    
    AFIR_thy          = AFIR,
    AFIR_simple_thy   = AFIR_simple,
    
    SCIM_thy          = SCIM,
    SCIM_adhoc_thy    = SCIM_adhoc,
    SCIM_simpler_thy  = SCIM_simpler,
    SCIM_simplest_thy = SCIM_simplest,
    SCIM_Lfold_thy    = SCIM_Lfold,
    SCIM_Lfold_adhoc_thy = SCIM_Lfold_adhoc,
    
    stringsAsFactors = FALSE
  )
  return(lumped_parameters_theory)
}

# Simulation ----------------------------------------------------------------------------------
lumped.parameters.simulation = function(model           = model, 
                                        param.as.double = param.as.double,
                                        dose.nmol       = dose.nmol, 
                                        tmax            = tmax, 
                                        tau             = tau, 
                                        compartment,
                                        infusion        = FALSE){
  
  # Arguments:
  #   model_name: name of the model
  #   params_file_path: full path of the parameters file.
  #   dose.nmol: dosing amount in nmol
  #   tmax: maximum doing period in days
  #   tau: dosing interval in days
  #   compartment: compartment to which dosing is applied
  #   (in model F case, compartment=2)
  #   infusion.  default FALSE.  If True, then dose is a long infusino, throughout tau
  # Return:
  #   A data frame of lumped parameters calculated from simulation
  
  # Run simulation
  #d <- xlsx::read.xlsx(params_file_path, 1)
  #param.as.double = d$Value
  #names(param.as.double) = d$Parameter
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
  
  error_simulation = FALSE
  out = tryCatch({
    model$rxode$solve(model$repar(param.as.double), ev, init, atol = 1e-12, rtol = 1e-12)
  }, error = function(e) {
    data.frame(time = NA, D = NA, T = NA, DT = NA, L = NA, TL = NA)
  })
  
  if (!is.na(out[1,"time"])) {
    out  = model$rxout(out)
    
    # Calculate initial condition
    initial_state = out %>% filter(time==0)
    TL0 = initial_state$TL
    T0  = initial_state$T
    L0  = initial_state$L
    
    #Calculate steady state
    dose_times       = seq(0,tmax,by=tau)
    t_last_dose      = dose_times[length(dose_times)]
    id_last_dose     = which(out$time==t_last_dose)
    id_last_troughss = id_last_dose-1
    
    t_penultimate_dose      = dose_times[length(dose_times)-1]
    id_penultimate_dose     = which(out$time==t_penultimate_dose)
    id_penultimate_troughss = id_penultimate_dose-1
    
    #time_idx = tmax - 0.1
    TLss = out$TL[id_last_troughss]
    Lss  = out$L[id_last_troughss]
    Dss  = out$D[id_last_troughss]
    Tss  = out$T[id_last_troughss]
    
    Ttotss = out$Ttot[id_last_troughss] 
    TLss_prev = out$TL[id_penultimate_troughss]
    
    #first dose is at tau - this is half way to tau.  
    #it's meant for a check, to make sure we start at steady state
    id_05tau  = which.min((out$time-tau/2)^2)
    TL05tau   = out$TL[id_05tau] 
    
    SCIM = TLss/TL0
    AFIR = Tss/T0
  } else {
    error_simulation = TRUE
    TL0 = T0 = L0 = Ttotss = Lss = Dss = AFIR = SCIM = 
      t_last_dose = TLss = Tss = TLss_frac_change = TL0_05tau_frac_change = 
      TLss_prev = TL05tau = NA
  }
  lumped_parameters_sim = data.frame(
    TL0_sim = TL0,
    T0_sim = T0,
    L0_sim = L0,
    TssLss_TLss_sim = Tss*Lss/TLss,
    
    Ttotss_sim = Ttotss,
    Tss_sim = Tss,
    Lss_sim = Lss,
    Dss_sim = Dss,
    TLss_sim = TLss,
    AFIR_sim = AFIR,
    SCIM_sim = SCIM,
    time_last_dose = t_last_dose,
    TLss_frac_change = (TLss-TLss_prev)/TLss, #can be used to check we're at steady state
    TL0_05tau_frac_change = (TL05tau-TL0)/TL0,
    error_simulation = as.numeric(error_simulation),
    stringsAsFactors = FALSE) 
}

# Theory + Simulation: Compare ----------------------------------------------------------------------------------
#  on the user inputted parameter 

# Input: 
# model - model system of ODE's solved with RxODE. In this project, it is 'ivsc_4cmtct_shedct'.
# param.as.double - read parameters from Excel file. read.param.file("file directory").
# dose.nmol - dose in nmol.
# tmax - time of treatment in days
# tau - frequency of administering dose in days
# compartment - compartment where drug is administered
# param.to.change - parameter on which to do SA. This must be a string.
# param.to.change.range - range of parameter on which to do SA. The range must be symmetric in fold change. This must be a vector of odd length.

# Output:
# Data frame of AFIRT vs parameter value

compare.thy.sim = function(model                 = model,
                           param.as.double       = param.as.double,
                           dose.nmol             = dose.nmol,
                           tmax                  = tmax,
                           tau                   = tau,
                           compartment           = compartment,
                           param.to.change       = param.to.change,
                           param.to.change.range = param.to.change.range,
                           infusion              = FALSE) {
  
  # Store the orignal parameter set and parameter to be changed. 
  # This is needed to divide by the baseline value when calculating the fold change.
  param.original = param.as.double[param.to.change]
  dose.original  = dose.nmol
  param.to.change.name = param.to.change
  
  #Iterate through parameters
  df_sim = list()
  df_thy = list()
  i      = 0
  for (param.iter in param.to.change.range){
    if (param.to.change == 'dose'){
      dose.nmol = param.iter
    } else {
      param.as.double[param.to.change] = param.iter
    }
    
    #KEY LINES FOR COMPUTING THEORY AND SIMULATION
    i=i+1
    df_sim[[i]] = lumped.parameters.simulation(model, param.as.double, dose.nmol, tmax, tau, compartment)
    df_thy[[i]] = lumped.parameters.theory    (       param.as.double, dose.nmol,       tau)
  }
  
  #store final results in data.frame    
  df_thy = bind_rows(df_thy) %>% mutate(param.to.change = param.to.change.range)
  df_sim = bind_rows(df_sim) #%>% mutate(param.to.change = param.to.change.range)
  
  # Arrange theory and simulation in single data frame.
  df_compare = bind_cols(df_thy,df_sim)
  df_compare = df_compare %>%
    mutate(param = param.to.change.name) %>%
    mutate_if(is.numeric,signif,6)
  
  if (param.to.change == 'dose'){
    df_compare = df_compare %>% mutate(fold.change.param = param.to.change.range/dose.original)
  } else {
    df_compare = df_compare %>% mutate(fold.change.param = param.to.change.range/param.original)
  }
  
  return(df_compare)
}

# plot results ----
plot_param = function(param = param,
                      model = model,
                      infusion = TRUE,
                      plot_flag = TRUE) {
  
  tmax = param$tmax 
  tau  = param$tau  
  dose_nmol = param$dose_nmol
  compartment = 2

  nam   = names(param)
  options(warn = -1)
  param_as_double = param %>%
    as.numeric() %>%
    setNames(nam)
  options(warn = 0)
  param_as_double = param_as_double[model$pin]
  
  param_print = param_as_double %>%
    t() %>%
    as.data.frame() %>%
    mutate(id = param$id,
           CL = signif(keD/V1,2)) %>%
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
    mutate(Kss_TL = param$Kss_TL) %>%
    select(type, Dss, T0, L0, TL0, Ttotss, Lss, TLss, TssLss_TLss, Kss_TL, AFIR, SCIM)

  init = model$init(param_as_double)
  out  = model$rxode$solve(model$repar(param_as_double), ev, init, atol = 1e-12, rtol = 1e-12)
  out  = model$rxout(out)
  
  out_plot = out %>%
    select(time,D,T,DT,L,TL) %>%
    gather(cmt,value,-time)
  out_last = out_plot[(out$time==max(out$time)),]

  g = ggplot(out_plot,aes(x=time,y=value, color = cmt, group= cmt))
  g = g + geom_line()
  g = g + geom_label(data = out_last, aes(label = cmt), show.legend = FALSE, hjust=1)
  g = g + geom_hline(yintercept = thy$Ccrit_thy, linetype = "dashed", color = "red")
  g = g + annotate(geom = "label", x = 0.9*tmax, y = thy$Ccrit_thy, color = "red", label = "Ccrit")
  g = g + geom_vline(xintercept = tau, linetype = "dotted")
  g = g + xgx_scale_x_time_units(units_dataset = "days", units_plot = "weeks")
  g = g + xgx_scale_y_log10()
  g = g + labs(y = "Concentration (nM)", color = "")
  g = g + ggtitle(paste0(  "id = ",param$id,
                           "\nAFIR_thy                       = ",signif(thy$AFIR_thy,2),
                           "\nAFIR_sim                      = " ,signif(sim$AFIR_sim,2),
                           "\nSCIM_Lfold_adhoc_thy = "          ,signif(thy$SCIM_Lfold_adhoc_thy,2),
                           "\nSCIM_adhoc_thy           = "      ,signif(thy$SCIM_adhoc_thy,2),
                           "\nSCIM_sim                      = " ,signif(sim$SCIM_sim,2)))
  if (plot_flag == TRUE) {
    print(g)
  }
  
  par#unfortunately, kable does not work properly inside for loop
  out = list(
    param   = param_print,
    compare = compare,
    plot    = g,
    sim     = out)
  return(out)
}
