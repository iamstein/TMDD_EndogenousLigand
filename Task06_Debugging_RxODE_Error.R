library(RxODE)

p = readxl::read_excel("parameters/Param_Crash.xlsx") 
nam = names(p)
p$Vm = 0
p$Km = 1

p = as.numeric(p)
names(p) = nam

dose.nmol = 0
tmax = 364
tau = 21
compartment = 2

ev = eventTable(amount.units="nmol", time.units="days")
sample.points = c(seq(-7, tmax, 0.1), 10^(-3:0)) # sample time, increment by 0.1
sample.points = sort(sample.points)
sample.points = unique(sample.points)
ev$add.sampling(sample.points)
ev$add.dosing(dose=dose.nmol, nbr.doses=floor(tmax/tau)+1, dosing.interval=tau,
              dosing.to=compartment)

init = c(AmtD0=0,AmtD=0,AmtD2=0,T=0,DT=0,L=0,TL=0)
init["T"]=p["ksynT"]/p["keT"]
init["L"]=p["ksynL"]/p["keL"]


equations = '
  D            = AmtD/V1;
  d/dt(AmtD0)  =       -ka*AmtD0;
  d/dt(AmtD )  = V1*(F *ka*AmtD0 - kon_DT*D*T + koff_DT*DT                                   - keD *D - Vm*D/(D+Km)) - k12*AmtD  + k21*AmtD2;
  d/dt(AmtD2)  =                                                                                                     - k21*AmtD2 + k12*AmtD ; 
  d/dt(T)      =                 - kon_DT*D*T + koff_DT*DT - kon_TL*L*T + koff_TL*TL + ksynT - keT *T ; 
  d/dt(DT)     =                   kon_DT*D*T - koff_DT*DT                                   - keDT*DT; 
  d/dt(L)      =                                           - kon_TL*L*T + koff_TL*TL + ksynL - keL *L; 
  d/dt(TL)     =                                             kon_TL*L*T - koff_TL*TL         - keTL*TL; 
  '

model = RxODE(model = equations)

out  = model$solve(p, ev, init)

head(out)
