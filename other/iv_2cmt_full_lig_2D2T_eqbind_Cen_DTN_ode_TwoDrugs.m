function[dY] = iv_2cmt_full_lig_2D2T_eqbind_Cen_DTN_ode(t,Y,p)


AmtD1  = Y(1); D1 = AmtD1/p.VcD1;  %amount of free drug D in central cmt
AmtD1p = Y(2);                    %amount of free drug D in peripheral cmt
Ts     = Y(3);                    %free target
Tm     = Y(4);
D1Ts   = Y(5);                    %drug-target complex
D1Tm   = Y(6);
N      = Y(7);                    %amount of native species that can also bind target
NTs    = Y(8);                    %endogenous species - target complex
NTm    = Y(9);
D2N    = Y(10);
AmtD2  = Y(11); D2 = AmtD2/p.VcD2; %amount of free drug D in central cmt
AmtD2p = Y(12);                   %amount of free drug D in peripheral cmt

T      = Ts +  Tm;
D1T    = D1Ts + D1Tm;
NT     = NTs + NTm;

if isfield(p,'ksynN_fun')
    p.ksynN = feval(p.ksynN_fun,t);
end

%             DRUG-TARGET BINDING                 TARGET-ENDOGENOUS BINDING       DOSE/SYNTHESIS     ENIMINATION     %DISTRIBUTION					    TRANSIT TO PERIPHERAL  
dAmtD1  =  (- p.konD1*D1*T + p.koffD1*D1T)*p.VcD1                                                  - p.keD1 *AmtD1    - p.kcpD1*AmtD1 + p.kpcD1*AmtD1p;   %AmtD1 - free drug central
dAmtD1p =                                                                                                               p.kcpD1*AmtD1 - p.kpcD1*AmtD1p;   %AmtD1p- free drug peripheral
dTs     =   - p.konD1*D1*Ts+ p.koffD1*D1Ts      - p.konN*N*Ts + p.koffN*NTs     + p.ksynTs        - p.keTs *Ts;    %T    - free target
dTm     =   - p.konD1*D1*Tm+ p.koffD1*D1Tm      - p.konN*N*Tm + p.koffN*NTm     + p.ksynTm        - p.keTm *Tm;    %T    - free target
dD1Ts   =     p.konD1*D1*Ts- p.koffD1*D1Ts                                                        - p.keD1Ts*D1Ts; %D1T   - drug-target complex
dD1Tm   =     p.konD1*D1*Tm- p.koffD1*D1Tm                                                        - p.keD1Tm*D1Tm; %D1T   - drug-target complex

dNTs    =                                       + p.konN*N*Ts - p.koffN*NTs                       - p.keNTs*NTs;   %NT   - target-endogenous complex
dNTm    =                                       + p.konN*N*Tm - p.koffN*NTm                       - p.keNTm*NTm;   %NT   - target-endogenous complex

dN      =   - p.konD2*D2*N + p.koffD2*D2N       - p.konN*N*T  + p.koffN*NT      + p.ksynN         - p.keN  *N;     %AmtN - free native species (endogenous)
dD2N    =     p.konD2*D2*N - p.koffD2*D2N                                                         - p.keD2N*D2N;
dAmtD2  =  (- p.konD2*D2*N + p.koffD2*D2N)*p.VcD2                                                 - p.keD2 *AmtD2    - p.kcpD2*AmtD2 + p.kpcD2*AmtD2p;    %AmtD2 - free drug central
dAmtD2p =                                                                                                               p.kcpD2*AmtD2 - p.kpcD2*AmtD2p;    %AmtD2p- free drug peripheral

%     1      2       3   4   5     6     7  8    9    10   11     12
dY = [dAmtD1 dAmtD1p dTs dTm dD1Ts dD1Tm dN dNTs dNTm dD2N dAmtD2 dAmtD2p]';
1;
