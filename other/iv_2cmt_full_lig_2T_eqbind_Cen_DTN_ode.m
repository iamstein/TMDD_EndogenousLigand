possible initial condition as a simpler place to start
        Ts0       = p.ksynTs/p.keTs;
        Tm0       = p.ksynTm/p.keTm;
        N0        = p.ksynN/p.keN;
        NTs0      = 1e-15; 
        NTm0      = 1e-15;


Please switch notation though to 
	S = soluble target
	M = membrane-bound target


function[dY] = iv_2cmt_full_lig_2T_eqbind_Cen_DTN_ode(~,Y,p)


AmtD   = Y(1); D = AmtD/p.Vc;  %amount of free drug D in central cmt
AmtDp  = Y(2);                 %amount of free drug D in peripheral cmt
Ts     = Y(3);                 %free target
Tm     = Y(4);
DTs    = Y(5);                 %drug-target complex
DTm    = Y(6);
N      = Y(7);                 %amount of native species that can also bind target
NTs    = Y(8);                 %endogenous species - target complex
NTm    = Y(9);

T      =  Ts +  Tm;
DT     = DTs + DTm;
NT     = NTs + NTm;

%            DRUG-TARGET BINDING            TARGET-ENDOGENOUS BINDING      DOSE/SYNTHESIS    ENIMINATION    					    TRANSIT TO PERIPHERAL  
dAmtD  =  (- p.konD*D*T + p.koffD*DT)*p.Vc                                + p.dose_inf      - p.keD *AmtD - p.kcp*AmtD + p.kpc*AmtDp;   %AmtD - free drug central
dAmtDp =                                                                                                    p.kcp*AmtD - p.kpc*AmtDp;   %AmtDp- free drug peripheral
dTs    =   - p.konD*D*Ts+ p.koffD*DTs      - p.konN*N*Ts+ p.koffN*NTs     + p.ksynTs        - p.keTs *Ts;  %T    - free target
dTm    =   - p.konD*D*Tm+ p.koffD*DTm      - p.konN*N*Tm+ p.koffN*NTm     + p.ksynTm        - p.keTm *Tm;  %T    - free target
dDTs   =     p.konD*D*Ts- p.koffD*DTs                                                       - p.keDTs*DTs; %DT   - drug-target complex
dDTm   =     p.konD*D*Tm- p.koffD*DTm                                                       - p.keDTm*DTm; %DT   - drug-target complex
dN     =                                   - p.konN*N*T + p.koffN*NT      + p.ksynN         - p.keN  *N;   %AmtN - free native species (endogenous)
dNTs   =                                   + p.konN*N*Ts- p.koffN*NTs                       - p.keNTs*NTs; %NT   - target-endogenous complex
dNTm   =                                   + p.konN*N*Tm- p.koffN*NTm                       - p.keNTm*NTm; %NT   - target-endogenous complex


dY = [dAmtD dAmtDp dTs dTm dDTs dDTm dN dNTs dNTm]';
1;
