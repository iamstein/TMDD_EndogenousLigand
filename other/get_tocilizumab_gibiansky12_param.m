function[p q units] = get_tocilizumab_gibiansky12_param()
% From Gibiansky and Frey, J PKPD, 39, 5 (2012)

%from Table 1
    p.CL   = .289*1e3; %ml/d : L/d * 1000ml/L
    p.Vc   = 3.71*1e3; %ml   : L   * 1000ml/L
    p.Q    = .202*1e3; %ml/d : L/d * 1000ml/L
    p.Vp   = 3.45*1e3; %ml   : L   * 1000ml/L
    p.Vm   = 1.41*(p.Vc*1e3); %ug/d -> ug/d
    p.Km   = .367    ; %ug/ml 
    
    p.om2_CL = .076;
    p.om2_Vc = .100;
    p.om2_Q  = .129;
    p.om2_Vp = .296;
    p.om2_Vm = .056;
    
    
%from Table 2
    p.ksyn = 47.8*1e-3;%ug/ml/d : ng/ml/day * 1ug/1e3ng
    p.kdeg = 1.26    ; %1/d
    p.Imax = .939    ; %-   
    p.Kss  = .182    ; %ug/ml
    
    p.om2_kdeg     = .056;
    p.om2_kdeg_Kss = -.128;
    p.om2_Kss      = .895;
        
    %kint = -(p.Imax - 1)*p.kdeg; %Equation 2, solving for kint

%units
    units.Time      = 'days';
    units.Dose      = 'ug';
    units.Conc      = 'ug/ml';
    
    units.Volume    = 'ml';
    units.Clearance = 'ml/day';
    
    units.ksyn      = 'ug/ml/d';
    units.Vmax      = 'ug/d';

    q      = micro2macro_bolus_2cmt(p);
