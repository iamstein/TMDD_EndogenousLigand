d/dt(AmtD0)  =       -ka*AmtD0                                                                                               ;
d/dt(AmtD )  = V1*(F *ka*AmtD0 - kon_DT*D*T + koff_DT*DT                                   - keD *D) - k12*AmtD  + k21*AmtD2);
d/dt(AmtD2)  =                                                                                       - k21*AmtD2 + k12*AmtD  ; 
d/dt(T)       =                - kon_DT*D*T + koff_DT*DT - kon_TL*L*T + koff_TL*LT + ksynT - keT *T ; 
d/dt(DT)     =                   kon_DT*D*T - koff_DT*DT                                   - keDT*DT; 
d/dt(L)      =                                           - kon_TL*L*T + koff_TL*LT + ksynL - keL  *L; 
d/dt(LT)     =                                           + kon_TL*L*T - koff_TL*LT         - keLT*LT; 
