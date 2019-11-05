ivsc_2cmt_RR_v1 = function() {
  model           = list()
  model$name = 'ivsc_2cmt_DTL'
  
  #COMPARTMENTS
  model$cmtshort  = c('AmtD0','AmtD','AmtD2','T','DT','L','LT')
  #CALCULATE INITIAL CONDITION WITH NO DRUG PRESENT AND ASSUMING STEADY STATE
  model$init      = function(p){
    p    = p %>% t() %>% as.data.frame()

    if (p$keTL == 0) {
      TL0 = with(p,kon_TL*ksynT*ksynL/(koff_TL*keL*keT))
      T0  = with(p,ksynT/keT)
      L0  = with(p,ksynL/keL)
    } else {
      a = with(p,keTL^2)
      b = with(p,-(keTL) * (ksynT + ksynL) - (((koff_TL+keTL)/kon_TL) * keT *keL))
      c = with(p, ksynL*ksynT)
      
      TL0 = ( (-b) -sqrt((b^2)-4*a*c) ) / (2*a)
      T0  = with(p,(ksynT - keTL*TL0)/keT)
      L0 =  with(p,(ksynL - keTL*TL0)/keL)
    }    
    
    init = c(AmtD0=0,
             AmtD=0,
             AmtD2=0,
             T=T0,
             DT=0,
             L=L0,
             TL=TL0)

    return(init)
  }
  
  #PARAMEETRS IN MODEL
  model$pin       = c('F','ka','V1',
                      'k12','k21',
                      'ksynT','ksynL',
                      'keD','keT','keL','keDT','keTL',
                      'Vm','Km',
                      'kon_DT','koff_DT','kon_TL','koff_TL'); #input parameters
  model$pode      = model$pin
  
  model$repar = function(p) {return(p)}
  
  #                        INPUT      BINDING DRUG                BINDING LIGAND        SYNTHESIS  ELIMINATION            DISTRIBUTION
  model$rxode.str = '
  D            = AmtD/V1;
  d/dt(AmtD0)  =       -ka*AmtD0;
  d/dt(AmtD )  = V1*(F *ka*AmtD0 - kon_DT*D*T + koff_DT*DT                                   - keD *D - Vm*D/(D+Km)) - k12*AmtD  + k21*AmtD2;
  d/dt(AmtD2)  =                                                                                                     - k21*AmtD2 + k12*AmtD ; 
  d/dt(T)      =                 - kon_DT*D*T + koff_DT*DT - kon_TL*L*T + koff_TL*TL + ksynT - keT *T ; 
  d/dt(DT)     =                   kon_DT*D*T - koff_DT*DT                                   - keDT*DT; 
  d/dt(L)      =                                           - kon_TL*L*T + koff_TL*TL + ksynL - keL *L; 
  d/dt(TL)     =                                             kon_TL*L*T - koff_TL*TL         - keTL*TL; 
  '
  
  model$rxode     = RxODE(model = model$rxode.str, modName = model$name)
  
  model$rxout     = function(result)    {
    result        = as.data.frame(result)
    result = mutate(result,
                    Ttot = T + DT + TL,
                    Ltot = L + TL,
                    Dtot = D + DT)
  }
  
  return(model)
}


ivsc_2cmt_RR_KssT0L0 = function() {
  model           = ivsc_2cmt_RR_v1()
  model$name      = 'ivsc_2cmt_DTL'
  
  #COMPARTMENTS
  model$cmtshort  = c('AmtD0','AmtD','AmtD2','T','DT','L','LT')
  #CALCULATE INITIAL CONDITION WITH NO DRUG PRESENT AND ASSUMING STEADY STATE
  model$init      = function(p){
    p    = p %>% t() %>% as.data.frame()
    
    TL0 = with(p,T0*L0/Kss_TL)
    
    init = c(AmtD0 = 0,
             AmtD  = 0,
             AmtD2 = 0,
             T     = p$T0,
             DT    = 0,
             L     = p$L0,
             TL    = TL0)
    
    return(init)
  }
  
  #PARAMEETRS IN MODEL
  model$pode       =  c('F','ka','V1', 'k12','k21','ksynT','ksynL','keD','keT','keL','keDT','keTL', 'Vm','Km','kon_DT','koff_DT','kon_TL','koff_TL'); #ode parameters
  model$pin        =  c('F','ka','V1', 'k12','k21','T0'   ,'L0'   ,'keD','keT','keL','keDT','keTL', 'Vm','Km','kon_DT','Kss_DT' ,'kon_TL','Kss_TL'); #input parameters
  
  model$repar = function(p) {
    p     = as.data.frame(as.list(p))
    TL0   = with(p,T0*L0/Kss_TL)
    p     = mutate(p,
                   koff_TL = Kss_TL*kon_TL - keTL,
                   koff_DT = Kss_DT*kon_DT - keDT,
                   ksynT   = T0*keT + keTL*TL0,
                   #ksynL   = L0*(kon_TL*T0 + keL) - koff_TL*TL0) 
                   ksynL   = L0*keL + keTL*TL0)
    return(unlist(p))
  }
  return(model)
}


ivsc_2cmt_RR_KdT0L0 = function() {
  model           = ivsc_2cmt_RR_v1()
  model$name      = 'ivsc_2cmt_DTL'
  
  #COMPARTMENTS
  model$cmtshort  = c('AmtD0','AmtD','AmtD2','T','DT','L','LT')
  #CALCULATE INITIAL CONDITION WITH NO DRUG PRESENT AND ASSUMING STEADY STATE
  model$init      = function(p){
    p       = p %>% t() %>% as.data.frame()
    Kss_TL  = p$Kd_TL + p$keTL/p$kon_TL
    TL0     = with(p, T0*L0/Kss_TL)
    
    init = c(AmtD0 = 0,
             AmtD  = 0,
             AmtD2 = 0,
             T     = p$T0,
             DT    = 0,
             L     = p$L0,
             TL    = TL0)
    
    return(init)
  }
  
  #PARAMETERS IN MODEL
  model$pode       =  c('F','ka','V1', 'k12','k21','ksynT','ksynL','keD','keT','keL','keDT','keTL', 'Vm','Km','kon_DT','koff_DT','kon_TL','koff_TL'); #ode parameters
  model$pin        =  c('F','ka','V1', 'k12','k21','T0'   ,'L0'   ,'keD','keT','keL','keDT','keTL', 'Vm','Km','kon_DT','Kd_DT'  ,'kon_TL','Kd_TL'); #input parameters
  
  model$repar = function(p) {
    p     = as.data.frame(as.list(p))
    Kss_TL= with(p, Kd_TL + keTL/kon_TL)
    TL0   = with(p, T0*L0/Kss_TL)
    p     = mutate(p,
                   koff_TL = Kd_TL*kon_TL,
                   koff_DT = Kd_DT*kon_DT,
                   ksynT   = T0*keT + keTL*TL0,
                   #ksynL   = L0*(kon_TL*T0 + keL) - koff_TL*TL0)
                   ksynL   = L0*keL + keTL*TL0)
    return(unlist(p))
  }
  return(model)
}
