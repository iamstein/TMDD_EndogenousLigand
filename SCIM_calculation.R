# Helper function that returns a range of variable when performing sensitvity analysis  -----------------------------------------------------
read.param.file = function(filename) {
  d                      = read_excel(filename, 1)
  param.as.double        = suppressWarnings(as.numeric(d$Value))
  names(param.as.double) = d$Parameter
  param.as.double        = param.as.double[model$pin] #keep only parameters used in ODE
}

lseq = function(from, to, length.out){
  # Arguments:
  #   from : initial value of the variable
  #   to : teminal value of the variable
  #   length.out : fold number of <to - from>
  # Return :
  #   A vector of length <length.out>
  
  sequence = seq(log(from), log(to), length.out=length.out)
  sequence = exp(sequence)
  return(sequence)
}

# Theory #----------------------------------------------------------------------------------
lumped.parameters.theory = function(param.as.double = param.as.double,
                                    dose.nmol       = dose.nmol,
                                    tau             = tau){
  # Arguments:
  #   params_file_path: full path of the parameters file.
  #   dose.nmol: dosing amout in nmol
  #   tau: dosing interval in days
  # Return:
  #   A data frame of lumped parameters calculated from theory
  
  pars    = as.data.frame(t(param.as.double))
  #CL = with(pars, (keD*V1))
  #TL0 = with(pars, ksynT*ksynL/(koff_TL/kon_TL*keT*(keL+kon_TL) - ksynT*koff_TL))
  T0  = with(pars, ksynT/keT)
  TL0_keDT0 = with(pars, ksynL*T0/(koff_TL/kon_TL*keL))
  #Tacc = with(pars,keT  /keDT)
  Lss  = with(pars,ksynL/keL)
  Ttot = with(pars,ksynT/keDT)
  KssTL = with(pars,(koff_TL+keTL)/kon_TL)
  KssDT = with(pars,(koff_DT+keDT)/kon_DT)
  
  
  #compute Ctrough
  dose = dose.nmol
  alpha= with(pars, .5*(k12+k21+keD + sqrt( (k12+k21+keD)^2 - 4*keD*k21)))
  beta = with(pars, .5*(k12+k21+keD - sqrt( (k12+k21+keD)^2 - 4*keD*k21)))
  A    = with(pars, dose*(k21-alpha)/(V1*(beta-alpha)))
  B    = with(pars, dose*(k21-beta) /(V1*(alpha-beta)))
  Dss = (A*exp(-alpha*tau)/(1-exp(-alpha*tau)) + B*exp(-beta*tau)/(1-exp(-beta*tau)))
  
  
  #AFIR_thy = with(pars, 1/(koff_TL/kon_TL * Dss /(Lss*koff_DT/kon_DT) + koff_TL/kon_TL/Lss + 1)*(ksynT/(keDT*TL0)))
  SCIM_thy_keTL0 = with(pars, 1/(koff_TL/kon_TL * Dss /(Lss*((koff_DT+keDT)/kon_DT)) + koff_TL/kon_TL/Lss + 1)*(ksynT/(keDT*TL0_keDT0)))
  AFIR_thy_simple = with(pars,(KssDT*(Ttot/T0))/Dss)
  
  Tacc     = Ttot/T0
  AFIR_thy = with(pars,(KssDT*Tacc)/(Dss+Tacc))
  
  #For KeTL ~= 0
  a = with(pars,keTL^2)
  b = with(pars,-(keTL) * (ksynT +ksynL) - (((koff_TL+keTL)/kon_TL) * keT *keL))
  c =with(pars, ksynL*ksynT)
  
  TL0_pos <- ((-b) + sqrt((b^2)-4*a*c))/(2*a)
 
  TL0_neg <- ((-b) -sqrt((b^2)-4*a*c))/(2*a)
  
  
  SCIM_thy_ketl_pos = with(pars,Ttot/((((KssTL*Dss*keL)/(KssDT*ksynL))+((KssTL*keL)/(ksynL))+1)*TL0_pos))
  SCIM_thy_ketl_neg = with(pars,Ttot/((((KssTL*Dss*keL)/(KssDT*ksynL))+((KssTL*keL)/(ksynL))+1)*TL0_neg))
  
  SCIM_thy_ketl_neg_26 = with(pars,Ttot/KssTL * Lss/TL0_neg * KssDT/Dss)
  
  SCIM_thy_ketl_neg_29 = with(pars,Ttot/( ( ((KssTL*keL/ksynL)*(Dss/KssDT)) +1) *TL0_neg) )
  
  SCIM_thy_ketl_neg_31 = with(pars, (ksynT*ksynL*KssDT) / (KssTL*keL*keDT*Dss*TL0_neg) )

  lumped_parameters_theory = data.frame(
    TL0_keTL0_thy = TL0_keDT0,
    TL0_negroot_thy = TL0_neg,
    TL0_posroot_thy = TL0_pos,
    T0_thy = T0,
    Ttot_thy = Ttot,
    Lss_thy = Lss,
    Dss_thy = Dss,
    SCIM_thy_keTL0        = SCIM_thy_keTL0,
    # Compare simplified SCIM eqns. 
    # SCIM_thy_keTL_negroot is the most complex i.e not simplified version of SCIM
    # 26, 29, and 31 refer to the eqn numbers in the latex doc. for the simplified SCIMs
    SCIM_thy_keTL_negroot = SCIM_thy_ketl_neg,
    SCIM_thy_keTL_negroot26 = SCIM_thy_ketl_neg_26,
    SCIM_thy_keTL_negroot29 = SCIM_thy_ketl_neg_29,
    SCIM_thy_keTL_negroot31 = SCIM_thy_ketl_neg_31,
    SCIM_thy_keTL_posroot = SCIM_thy_ketl_pos,
    AFIR_thy = AFIR_thy,
    AFIR_thy_simple = AFIR_thy_simple,
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
  sample.points = c(seq(-7, tmax, 0.1), 10^(-3:0)) # sample time, increment by 0.1
  sample.points = sort(sample.points)
  sample.points = unique(sample.points)
  ev$add.sampling(sample.points)
  
  #add dur  tau for a long infusion
  if (infusion == FALSE) {
    ev$add.dosing(dose=dose.nmol, nbr.doses=floor(tmax/tau)+1, dosing.interval=tau, dosing.to=compartment)
  } else {
    ev$add.dosing(dose=dose.nmol, nbr.doses=floor(tmax/tau)+1, dosing.interval=tau, dosing.to=compartment, dur = tau)
  }  
  
  init = model$init(param.as.double)
  out  = model$rxode$solve(param.as.double, ev, init)
  out  = model$rxout(out)
  
  # Calculate initial condition
  initial_state = out %>%
    filter(time==0)
  TL0 = initial_state$TL
  T0 = initial_state$T
  
  #TODO - think more about this code - it's jsut a bit confusing here.  Why - 3?  Why 10*?
  idx_vec = seq((floor(tmax/tau) - 3)*tau, tmax, 0.1)
  out_vec = out[which(round(10*out$time) %in% round(10*idx_vec)),]
  time_idx = out_vec$time[which(out_vec$D == min(out_vec$D))]
  
  #time_idx = tmax - 0.1
  TLss = out$TL[which(round(10*out$time) == round(10*time_idx))]
  SCIM_sim = TLss/TL0
  
  Ttot_sim = out$Ttot[which(round(10*out$time) == round(10*time_idx))]
  L_sim = out$L[which(round(10*out$time) == round(10*time_idx))]
  D_sim = out$D[which(round(10*out$time) == round(10*time_idx))]
  
  lumped_parameters_sim = data.frame(
    TL0_sim = TL0,
    T0_sim = T0,
    Ttot_sim = Ttot_sim,
    L_sim = L_sim,
    D_sim = D_sim,
    SCIM_sim = SCIM_sim,
    time_idx,
    TLss,
    stringsAsFactors = FALSE) #having one named sim will be helpful later on in Task01, Task02, etc.
  
  return(lumped_parameters_sim)
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
  param.as.double.original = param.as.double
  param.to.change.original = param.to.change
  
  
  # SIMULATION: Iterate through parameters
  df_sim = data.frame()
  for (param.iter in param.to.change.range){
    if (param.to.change == 'dose'){
      dose.nmol = param.iter
    } else {
      param.as.double[param.to.change] = param.iter
    }
    #KEY LINE FOR COMPUTED PARAMETERS FROM SIMULATION
    row = lumped.parameters.simulation(model, param.as.double, dose.nmol, tmax, tau, compartment)
    df_sim = rbind(df_sim, row)
  }
  
  if (param.to.change == 'dose'){
    df_sim = df_sim %>% mutate(param.to.change = param.to.change.range,
                               fold.change.param = param.to.change.range/dose.nmol)
  } else {
    df_sim = df_sim %>% mutate(param.to.change = param.to.change.range,
                               fold.change.param = param.to.change.range/param.as.double.original[param.to.change.original])
  }
  
  
  #THEORY: Iterate through parameters
  df_thy = data.frame()
  for (param.iter in param.to.change.range){
    if(param.to.change == 'dose'){
      dose.nmol = param.iter 
    } else {
      param.as.double[param.to.change] = param.iter
    }
    #KEY LINE FOR COMPUTED PARAMETERS FROM THEORY
    row = lumped.parameters.theory(param.as.double, dose.nmol, tau)
    df_thy = rbind(df_thy, row)
  }
  
  if (param.to.change == 'dose'){
    df_thy = df_thy %>% mutate(param.to.change = param.to.change.range,
                               fold.change.param = param.to.change.range/dose.nmol)
  } else {
    df_thy = df_thy %>% mutate(param.to.change = param.to.change.range,
                               fold.change.param = param.to.change.range/param.as.double.original[param.to.change.original])
  }
  
  # Arrange theory and simulation in single data frame.
  df_compare = bind_cols(df_thy,df_sim)
  param      = param.to.change
  df_compare = df_compare %>%
    mutate(param = param) %>%
    mutate_if(is.numeric,signif,6)
  
  return(df_compare)
}
