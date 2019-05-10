#target = presence of target
#if TRUE, target is present
#if FALSE, then target concentratinons are zero
ivsc_2cmt_RR_v1 = function(target = TRUE) {
  model           = list()
  model$name = 'ivsc_2cmt_DTL'
  
  #COMPARTMENTS
  model$cmtshort  = c('AmtD0','AmtD','AmtD2','T','DT','L','LT')
  #CALCULATE INITIAL CONDITION WITH NO DRUG PRESENT AND ASSUMING STEADY STATE
  model$init      = function(p){
    init = c(AmtD0=0,AmtD=0,AmtD2=0,T=0,DT=0,L=0,TL=0)
    p    = p %>% t() %>% as.data.frame()
    p
    Ttot0_approx = with(p,ksynT/keT)
    Ltot0_approx = with(p,ksynL/keL)
    
    init["T"]=Ttot0_approx
    init["L"]=Ltot0_approx
    
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
                    Ttot = T+DT+TL,
                    Ltot = L+TL,
                    Dtot = D+DT)
  }
  
  return(model)
}