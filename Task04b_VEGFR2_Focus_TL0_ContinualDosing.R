source("ams_initialize_script.R")
source("SCIM_calculation.R")  
source("ivsc_2cmt_RR_V1.R")
dirs$Rscript.name = "Task02a_Sensitivity_SCIM_k_pars.R"
dirs$output.prefix= str_extract(dirs$Rscript.name,"^Task\\d\\d\\w?_")

model = ivsc_2cmt_RR_v1()

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

#HERE IS WHERE THE DOSING IS ADDED
#bolus dosing
#ev$add.dosing(dose=dose.nmol, nbr.doses=floor(tmax/tau)+1, dosing.interval=tau,
#              dosing.to=compartment)
#continual dosing
ev$add.dosing(dose=dose.nmol, nbr.doses=floor(tmax/tau)+1, dosing.interval=tau,
              dosing.to=compartment, dur = tau)



param.as.double =  read.param.file("parameters/Bevacizumab_VEGFR2_Params.xlsx") #ADD THIS LINE (CHANGED VARIABLE NAME TO PARAM)
init = model$init(param.as.double)
out  = model$rxode$solve(param.as.double, ev, init)
out  = model$rxout(out)

out.long = out %>%
  dplyr::select(time,T,L,TL) %>%
  gather(key,value,-c(time))
  
pars = as.data.frame(t(param.as.double))
a = with(pars,keTL^2)
b = with(pars,-(keTL) * (ksynT +ksynL) - (((koff_TL+keTL)/kon_TL) * keT *keL))
c = with(pars, ksynL*ksynT)

TL0_pos <- ((-b) + sqrt((b^2)-4*a*c))/(2*a)

TL0_neg <- ((-b) -sqrt((b^2)-4*a*c))/(2*a)



g = ggplot(out.long,aes(x=time,y=value,color=key))
g = g + geom_line()
g = g + scale_y_log10()
g = g + geom_hline(yintercept = TL0_neg,color="purple",linetype="dashed")
print(g)

