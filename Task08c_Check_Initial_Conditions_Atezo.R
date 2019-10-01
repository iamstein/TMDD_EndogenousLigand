source("ams_initialize_script.R")
source("SCIM_calculation.R")  
source("ivsc_2cmt_RR_V1.R")
dirs$Rscript.name = "Task08_Check_Initial_Conditions.R"
dirs$output.prefix= str_extract(dirs$Rscript.name,"^Task\\d\\d\\w?_")

#set up the model and dosing
model1  = ivsc_2cmt_RR_v1()
model2  = ivsc_2cmt_RR_KeqT0L0()

# Dose time, frequency, compartment, nominal dose
tmax = 52*7 #days
tau  = 21   #days
compartment = 2
dose.nmol = 100*scale.mpk2nmol

#compute the simulation with the final parameter set, just to make sure it's getting to steady state
ev = eventTable(amount.units="nmol", time.units="days")
sample.points = c(seq(0, tmax, 0.1), 10^(-3:0)) # sample time, increment by 0.1
sample.points = sort(sample.points)
sample.points = unique(sample.points)
ev$add.sampling(sample.points)
#ev$add.dosing(dose=dose.nmol, nbr.doses=floor(tmax/tau)+1, dosing.interval=tau,
#              dosing.to=compartment, dur = tau)

i=0
out_long = list()
drug = drugs[1]
#for (drug in drugs) {
  i=i+1
  filename = paste0("parameters/ModelG_",drug,"_Params.xlsx")
  param.as.double =  read.param.file(filename)
  
  #reparameterize for model with T0, L0
    p    = param.as.double %>% t() %>% as.data.frame()
    if (p$keTL == 0) {
      TL0 = with(p,kon_TL*ksynT*ksynL/(koff_TL*keL*keT))
      T0  = with(p,ksynT/keT)
      L0  = with(p,ksynL/keL)
    } else {
      a = with(p,keTL^2)
      b = with(p,-(keTL) * (ksynT +ksynL) - (((koff_TL+keTL)/kon_TL) * keT *keL))
      c = with(p, ksynL*ksynT)
      
      TL0 = ((-b) -sqrt((b^2)-4*a*c))/(2*a)
      T0  = with(p,(ksynT - keTL*TL0)/keT)
      L0 =  with(p,(ksynL - keTL*TL0)/keL)
    }    
  
    param.as.double["T0"] = T0
    param.as.double["L0"] = L0
    param.as.double["Kss_DT"] = with(p,(koff_DT + keDT)/kon_DT)
    param.as.double["Kss_TL"] = with(p,(koff_TL + keTL)/kon_TL)
    
  #solve the model
    init1 = model1$init(param.as.double[model1$pin])
    init2 = model2$init(param.as.double[model2$pin])
    print(bind_rows(init1,init2))
    
    p11   = param.as.double[model1$pin] #parameter 1
    p21   = model2$repar(param.as.double[model2$pin])[model2$pode] #convert param2 back to param1
    print(bind_rows(p11,p21))
    
    out1  = model1$rxode$solve(model1$repar(param.as.double), ev, init1)
    out1  = model1$rxout(out1)
    
    out2  = model2$rxode$solve(model2$repar(param.as.double), ev, init2)
    out2  = model2$rxout(out2)