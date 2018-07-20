#target = presence of target
#if TRUE, target is present
#if FALSE, then target concentratinons are zero
ivsc_4cmtct_shedct = function(target = TRUE) {
  model           = list()
  # model$name      = as.character(sys.calls()[[sys.nframe()]])
  model$name = 'ivsc_4cmtct_shedct'
  
  #COMPARTMENTS
  model$cmtshort  = c('AmtD0','D1','D2','D3','S1','S3','M1','M3','DS1','DS3','DM1','DM3')
  #CALCULATE INITIAL CONDITION WITH NO DRUG PRESENT AND ASSUMING STEADY STATE
  model$init      = function(p){
    init = c(AmtD0=0,D1=0,D2=0,D3=0,S1=0,S3=0,M1=0,M3=0,DS1=0,DS3=0,DM1=0,DM3=0)
    p    = p %>% t() %>% as.data.frame()
    
    
    ksyn = with(p,c(ksynS1,ksynS3,ksynM1,ksynM3))
    K    = with(p,matrix(c( keS1+k13S    , -k31S*VS3/VS1       , -kshedM1           ,  0,
                            -k13S*VS1/VS3 ,  keS3+k31S          ,  0                 , -kshedM3,
                            0            ,  0                  ,  kshedM1+keM1+k13M , -k31M*VD3/VD1,
                            0            ,  0                  , -k13M*VD1/VD3      ,  kshedM3+keM3+k31M), 
                         nrow = 4, byrow=TRUE))
    x    = solve(K,ksyn)
    
    if (target){
      init["S1"] = unlist(x[1])
      init["S3"] = unlist(x[2])
      init["M1"] = unlist(x[3])
      init["M3"] = unlist(x[4]) 
    } else {
      init["S1"] = 0
      init["S3"] = 0
      init["M1"] = 0
      init["M3"] = 0 
    }
    
    return(init)
  }
  
  #PARAMEETRS IN MODEL
  model$pin       = c('F','ka','VD1','VD2','VD3','VS1','VS3','VDS1','VDS3','VM1','VM3','VDM1','VDM3',
                      'k12D','k21D','k13D','k31D','k13S','k31S','k13DS','k31DS','k13M','k31M','k13DM','k31DM',
                      'ksynS1','ksynS3','ksynM1','ksynM3',
                      'keD1','keD3','keS1','keS3','keDS1','keDS3','keM1','keM3','keDM1','keDM3',
                      'kon1','koff1','kon3','koff3',
                      'kshedM3','kshedDM3','kshedM1','kshedDM1'); #input parameters
  model$pode      = model$pin
  
  #                INPUT/SYNTHESIS/SHED DISTRIBUTION (CENTRAL/PERIPH) DISTRIBUTION (CENTRAL/TUMOR)        BINDING           UNBINDING         ELIMINATION                     
  
  if (target){
    model$rxode.str = '
    D1           = AmtD1/VD1                                                                                                                               ;
    d/dt(AmtD0)  =  -ka *AmtD0                                                                                                                             ;
    d/dt(AmtD1)  =(F*ka *AmtD0/VD1 - k12D*D1 + k21D*VD2/VD1*D2      - k13D *D1  + k31D *VD3 /VD1* D3    - kon1*D1*(S1+M1) + koff1*(DS1+DM1) - keD1 *D1)*VD1;
    d/dt(D2)     =                 - k21D*D2 + k12D*VD1/VD2*D1                                                                                             ;
    d/dt(D3)     =                                                  - k31D *D3  + k13D *VD1 /VD3* D1    - kon3*D3*(S3+M3) + koff3*(DS3+DM3) - keD3 *D3     ;
    d/dt(S1)     = ksynS1 + kshedM1*M1                              - k13S *S1  + k31S *VS3 /VS1* S3    - kon1*D1*S1      + koff1*DS1       - keS1 *S1     ;
    d/dt(S3)     = ksynS3 + kshedM3*M3                              - k31S *S3  + k13S *VS1 /VS3* S1    - kon3*D3*S3      + koff3*DS3       - keS3 *S3     ;
    d/dt(M1)     = ksynM1 - kshedM1*M1                              - k13M *M1  + k31M *VM3 /VM1* M3    - kon1*D1*M1      + koff1*DM1       - keM1 *M1     ;
    d/dt(M3)     = ksynM3 - kshedM3*M3                              - k31M *M3  + k13M *VM1 /VM3* M1    - kon3*D3*M3      + koff3*DM3       - keM3 *M3     ;
    d/dt(DS1)    =          kshedDM1*DM1                            - k13DS*DS1 + k31DS*VDS3/VDS1*DS3   + kon1*D1*S1      - koff1*DS1       - keDS1*DS1    ;
    d/dt(DS3)    =          kshedDM3*DM3                            - k31DS*DS3 + k13DS*VDS1/VDS3*DS1   + kon3*D3*S3      - koff3*DS3       - keDS3*DS3    ;
    d/dt(DM1)    =        - kshedDM1*DM1                            - k13DM*DM1 + k31DM*VDM3/VDM1*DM3   + kon1*D1*M1      - koff1*DM1       - keDM1*DM1    ;
    d/dt(DM3)    =        - kshedDM3*DM3                            - k31DM*DM3 + k13DM*VDM1/VDM3*DM1   + kon3*D3*M3      - koff3*DM3       - keDM3*DM3    ;
    '
  } else {
    model$rxode.str = '
    D1           = AmtD1/VD1                                                                                                                               ;
    d/dt(AmtD0)  =  -ka *AmtD0                                                                                                                             ;
    d/dt(AmtD1)  =(F*ka *AmtD0/VD1 - k12D*D1 + k21D*VD2/VD1*D2      - k13D *D1  + k31D *VD3 /VD1* D3    - kon1*D1*(S1+M1) + koff1*(DS1+DM1) - keD1 *D1)*VD1;
    d/dt(D2)     =                 - k21D*D2 + k12D*VD1/VD2*D1                                                                                             ;
    d/dt(D3)     =                                                  - k31D *D3  + k13D *VD1 /VD3* D1    - kon3*D3*(S3+M3) + koff3*(DS3+DM3) - keD3 *D3     ;
    d/dt(S1)     = 0;
    d/dt(S3)     = 0;
    d/dt(M1)     = 0;
    d/dt(M3)     = 0;
    d/dt(DS1)    =          kshedDM1*DM1                            - k13DS*DS1 + k31DS*VDS3/VDS1*DS3   + kon1*D1*S1      - koff1*DS1       - keDS1*DS1    ;
    d/dt(DS3)    =          kshedDM3*DM3                            - k31DS*DS3 + k13DS*VDS1/VDS3*DS1   + kon3*D3*S3      - koff3*DS3       - keDS3*DS3    ;
    d/dt(DM1)    =        - kshedDM1*DM1                            - k13DM*DM1 + k31DM*VDM3/VDM1*DM3   + kon1*D1*M1      - koff1*DM1       - keDM1*DM1    ;
    d/dt(DM3)    =        - kshedDM3*DM3                            - k31DM*DM3 + k13DM*VDM1/VDM3*DM1   + kon3*D3*M3      - koff3*DM3       - keDM3*DM3    ;
    '
  }
  
  model$rxode     = RxODE(model = model$rxode.str, modName = model$name)
  
  model$rxout     = function(result)    {
    result        = as.data.frame(result)
    result = mutate(result,
                    Dtot1 = D1+DS1,
                    Stot1 = S1+DS1,
                    Dtot3 = D3+DS3,
                    Stot3 = S3+DS3,
                    Mtot1 = M1+DM1,
                    Mtot3 = M3+DM3)
  }
  
  return(model)
}
