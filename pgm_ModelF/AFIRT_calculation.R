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

# Theoretical lumped parameters #----------------------------------------------------------------------------------
lumped.parameters.theory = function(param.as.double = param.as.double,
                                    dose.nmol       = dose.nmol,
                                    tau             = tau,
                                    soluble         = FALSE){
    # Arguments:
    #   params_file_path: full path of the parameters file.
    #   dose.nmol: dosing amout in nmol
    #   tau: dosing interval in days
    #   soluble: flag saying whether or not drug is soluble. Default to FALSE.
    # Return:
    #   A data frame of lumped parameters calculated from theory

    pars    = as.data.frame(t(param.as.double))

    # Calculate M3tot.ss, M30, S3tot.ss and S30
    numerator.DM3   = with(pars, k13DM*(VD1/VD3)* ksynM1                     + (keDM1 + kshedDM1 + k13DM)* ksynM3)
    numerator.M3    = with(pars, k13M *(VD1/VD3)* ksynM1                     + (keM1  + kshedM1  + k13M )* ksynM3)
    numerator.DM1   = with(pars, k31DM*(VD3/VD1)* ksynM3                     + (keDM3 + kshedDM3 + k31DM)* ksynM1)
    numerator.M1    = with(pars, k31M *(VD3/VD1)* ksynM3                     + (keM3  + kshedM3  + k31M )* ksynM1)

    denomenator.DM_1or3 = with(pars, (keDM1 + kshedDM1 + k13DM)*(keDM3 + kshedDM3 + k31DM)-k31DM*k13DM)
    denomenator.M_1or3  = with(pars, (keM1  + kshedM1  + k13M )*(keM3  + kshedM3  + k31M )-k31M *k13M )
    
    M1tot.ss = numerator.DM1 / denomenator.DM_1or3
    M3tot.ss = numerator.DM3 / denomenator.DM_1or3
    M10      = numerator.M1  / denomenator.M_1or3
    M30      = numerator.M3  / denomenator.M_1or3   
    
    #note that this aligns with the numerator columns above and can be copied and pasted for comparison
    numerator.DS3  = with(pars, k13DS*(VD1/VD3)*(ksynS1 + kshedDM1*M1tot.ss) + (keDS1             + k13DS)*(ksynS3 + kshedDM3*M3tot.ss))
    numerator.S3   = with(pars, k13S *(VD1/VD3)*(ksynS1 + kshedM1 *M10)      + (keS1              + k13S) *(ksynS3 + kshedM3 *M30))
    
    denomenator.DS3= with(pars, (keDS1            + k13DS)*(keDS3            + k31DS)-k31DS*k13DS)
    denomenator.S3 = with(pars, (keS1             + k13S )*(keS3             + k31S )-k31S *k13S )
    
    S3tot.ss       = numerator.DS3 / denomenator.DS3
    S30            = numerator.S3  / denomenator.S3


    if (!soluble){
      Kssd = with(pars, (koff3 + keDM3 + kshedDM3 + k31DM)/kon3)
      Kss  = with(pars, (koff3 + keDM3 + kshedDM3        )/kon3)
      Kd   = with(pars,  koff3                            /kon3)
      Tacc.tum = M3tot.ss / M30
      
    } else {
      Kssd = with(pars, (koff3 + keDS3 +          + k31DS)/kon3)
      Kss  = with(pars, (koff3 + keDS3                   )/kon3)
      Kd   = with(pars,  koff3                            /kon3)
      Tacc.tum = S3tot.ss / S30
    }

    # Biodistribution coefficient (reference: ModelF_Appendix)
    B = with(pars, (k13D/(keD3 + k31D) * (VD1/VD3)))
    Biodist = B
    
    # Clearance
    CL = with(pars, (keD1*VD1))

    # Average drug concentration in the central compartment
    Cavg1 = dose.nmol/(CL*tau)

    # Compute various AFIRTs
    AFIRT.Kssd = Kssd*Tacc.tum/(B*Cavg1)
    AFIRT.Kss  = Kss *Tacc.tum/(B*Cavg1)
    AFIRT.Kd   = Kd  *Tacc.tum/(B*Cavg1)

    a0 = with(pars, keD1*k21D*k31D)
    a1 = with(pars, keD1*k31D + k21D*k31D + k21D*k13D + keD1*k21D + k31D*k12D)
    a2 = with(pars, keD1 + k12D + k13D + k21D + k31D)
    
    p  = a1 - (a2^2)/3
    q  = 2*(a2^3)/27 - a1*a2/3 + a0 
    r1 = (-(p^3)/27)^0.5
    r2 = 2*(r1^(1/3))
    
    phi = acos(-q/(2*r1))/3
    

    alpha = -(cos(phi)         *r2 - a2/3)
    beta  = -(cos(phi + 2*pi/3)*r2 - a2/3)
    gamma = -(cos(phi + 4*pi/3)*r2 - a2/3)

    V = with(pars, VD1)
    A = with(pars, (1/V) * ((k21D - alpha)/(alpha - beta)) * ((k31D - alpha)/(alpha - gamma)))
    B = with(pars, (1/V) * ((k21D - beta )/(beta - alpha)) * ((k31D - beta )/(beta - gamma )))
    C = with(pars, (1/V) * ((k21D - gamma)/(gamma - beta)) * ((k31D - gamma)/(gamma - alpha)))

    D = dose.nmol
    Cmin1 = D*((A*exp(-alpha*tau))/(1 - exp(-alpha*tau)) + 
              (B*exp(-beta *tau))/(1 - exp(-beta *tau)) + 
              (C*exp(-gamma*tau))/(1 - exp(-gamma*tau)))

    if(!soluble){
        Tfold = M3tot.ss/M30
    }else{
        Tfold = S3tot.ss/S30
    }
     
    TFIRT.Kssd = Kssd*Tfold/(Biodist*Cmin1)
    TFIRT.Kss  = Kss *Tfold/(Biodist*Cmin1)
    TFIRT.Kd   = Kd  *Tfold/(Biodist*Cmin1)

    # Implementation of TEC_50 and AFIRT*
    # Use Kssd as K_eq
    # I will use AFIRT_ for AFIRT*, because * is means multiplication in R 

    TEC_50 = Kssd*Tfold
    AFIRT.Kssd_ = TEC_50/(TEC_50 + Biodist*Cavg1)
    TFIRT.Kssd_ = TEC_50/(TEC_50 + Biodist*Cmin1)
    
    TEC_50      = Kss*Tfold
    AFIRT.Kss_  = TEC_50/(TEC_50 + Biodist*Cavg1)
    TFIRT.Kss_  = TEC_50/(TEC_50 + Biodist*Cmin1)
    
    TEC_50      = Kd*Tfold
    AFIRT.Kd_   = TEC_50/(TEC_50 + Biodist*Cavg1)
    TFIRT.Kd_   = TEC_50/(TEC_50 + Biodist*Cmin1)    

    lumped_parameters_theory = data.frame(type       = "theory",
                                          M30        = M30,
                                          S30        = S30,
                                          M3tot.ss   = M3tot.ss,
                                          S3tot.ss   = S3tot.ss,
                                          Tacc.tum   = Tacc.tum,
                                          B          = Biodist,
                                          Cavg1      = Cavg1,
                                          Cavg3      = Biodist*Cavg1,
                                          Cmin.thy   = Cmin1,
                                          Cmin       = Cmin1,
                                          AFIRT.Kssd = AFIRT.Kssd,
                                          AFIRT.Kss  = AFIRT.Kss,
                                          AFIRT.Kd   = AFIRT.Kd,
                                          AFIRT      = AFIRT.Kssd,
                                          TFIRT.Kssd = TFIRT.Kssd,
                                          TFIRT.Kss  = TFIRT.Kss,
                                          TFIRT.Kd   = TFIRT.Kd,
                                          TFIRT      = TFIRT.Kssd,
                                          TEC_50     = TEC_50,
                                          AFIRT.Kssd_= AFIRT.Kssd_,
                                          TFIRT.Kssd_= TFIRT.Kssd_,
                                          AFIRT.Kss_ = AFIRT.Kss_,
                                          TFIRT.Kss_ = TFIRT.Kss_,
                                          AFIRT.Kd_  = AFIRT.Kd_,                                          
                                          TFIRT.Kd_  = TFIRT.Kd_,                                          
                                          stringsAsFactors = FALSE
                                          )
    return(lumped_parameters_theory)
 }

# Simulated lumped parameters ----------------------------------------------------------------------------------

lumped.parameters.simulation = function(model           = model, 
                                        param.as.double = param.as.double,
                                        dose.nmol       = dose.nmol, 
                                        tmax            = tmax, 
                                        tau             = tau, 
                                        compartment,
                                        soluble         = FALSE){

    # Arguments:
    #   model_name: name of the model
            #   params_file_path: full path of the parameters file.
    #   dose.nmol: dosing amount in nmol
    #   tmax: maximum doing period in days
    #   tau: dosing interval in days
    #   compartment: compartment to which dosing is applied
    #   (in model F case, compartment=2)
    #   soluble: flag saying whether or not drug is soluble. Default to FALSE.
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
    ev$add.dosing(dose=dose.nmol, nbr.doses=floor(tmax/tau)+1, dosing.interval=tau,
                dosing.to=compartment)

    #model = do.call(model_name, list()) # model file can contain only one model
    init = model$init(param.as.double)
    out  = model$rxode$solve(param.as.double, ev, init)
    out  = model$rxout(out)
    out  = out %>%
        mutate(Sfree.pct = S3/init["S3"],
               Mfree.pct = M3/init["M3"],
               dose.nmol = dose.nmol)

    # Calculate lumped parameters
    initial_state = out %>%
        filter(time==0)
    M30 = initial_state$M3
    S30 = initial_state$S3

    ## Assume the system reaches steady state during the last dosing period
    steady_state = out %>%
      filter(time > (floor(tmax/tau)-1)*tau & time <tmax)
    M3tot.ss = mean(steady_state$Mtot3)
    S3tot.ss = mean(steady_state$Stot3)
    
    # Drug concentration
    dose_applied = out %>%
        filter(time > 0)
  
    # Average drug concentration in central compartment
    Cavg1 = mean(steady_state$D1)
    # Minimum drug concentration in central compartment
    Cmin1 = min(steady_state$D1)

    # Average drug concentration in tumor compartment
    Cavg3 = mean(steady_state$D3)
    # Minimum drug concentration in tumor compartment
    Cmin3 = min(steady_state$D3)
  
    # AFIRT and target accumulation
    if (soluble) {
      AFIRT    = mean(steady_state$Sfree.pct)
      Tacc.tum = S3tot.ss / S30
    } else {
      AFIRT    = mean(steady_state$Mfree.pct)
      Tacc.tum = M3tot.ss / M30
    }
    
    # Simulation of TFIRT
    if (soluble) {
        TFIRT = max(steady_state$Sfree.pct)
    } else {
        TFIRT = max(steady_state$Mfree.pct)
    } 


    lumped_parameters_sim = data.frame(type      = "simulation",
                                       M30       = M30,
                                       M3tot.ss  = M3tot.ss,
                                       S30       = S30,
                                       S3tot.ss  = S3tot.ss,
                                       Tacc.tum  = Tacc.tum,
                                       Cavg1     = Cavg1,
                                       Cavg3     = Cavg3,
                                       Cmin.sim  = Cmin1,
                                       Cmin      = Cmin1,
                                       B         = Cavg3/Cavg1,
                                       AFIRT     = AFIRT,
                                       AFIRT.sim = AFIRT,
                                       TFIRT     = TFIRT,
                                       TFIRT.sim = TFIRT,
                                       stringsAsFactors = FALSE) #having one named sim will be helpful later on in Task01, Task02, etc.

    return(lumped_parameters_sim)
}

# Compare Theory to Simulation for sensitivity analysis ----------------------------------------------------------------------------------
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
# soluble - boolean that is true/false if the drug is soluble/insoluble. Need this since soluble and insoluble are treated differently.

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
                           soluble               = FALSE) {
  
  # Store the orignal parameter set and parameter to be changed. 
  # This is needed to divide by the baseline value when calculating the fold change.
  param.as.double.original = param.as.double
  param.to.change.original = param.to.change
  
  # Simulation
  
  df_sim = data.frame()
  
  # Iterate through values in range.
  for (param.iter in param.to.change.range){
    if (param.to.change == 'dose'){
      dose.nmol = param.iter
    } else {
      param.as.double[param.to.change] = param.iter
    }
    row = lumped.parameters.simulation(model, param.as.double, dose.nmol, tmax, tau, compartment, soluble)
    df_sim = rbind(df_sim, row)
  }
  
  if (param.to.change == 'dose'){
    df_sim = df_sim %>% mutate(param.to.change = param.to.change.range,
                               fold.change.param = param.to.change.range/dose.nmol)
  } else {
    df_sim = df_sim %>% mutate(param.to.change = param.to.change.range,
                               fold.change.param = param.to.change.range/param.as.double.original[param.to.change.original])
  }

  # Theory
  
  df_thy = data.frame()
  
  # Iterate through values in range.
  for (param.iter in param.to.change.range){
    if(param.to.change == 'dose'){
      dose.nmol = param.iter 
    } else {
      param.as.double[param.to.change] = param.iter
    }
    row = lumped.parameters.theory(param.as.double, dose.nmol, tau, soluble)
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
  df_compare = bind_rows(df_thy,df_sim)
  param      = param.to.change
  df_compare = df_compare %>%
    mutate(param = param) %>%
    arrange(param.to.change,type) %>%
    mutate_if(is.numeric,signif,2)
  
  return(df_compare)
}
