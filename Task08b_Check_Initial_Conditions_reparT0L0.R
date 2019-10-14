source("ams_initialize_script.R")
source("SCIM_calculation.R")  
source("ivsc_2cmt_RR_V1.R")
dirs$Rscript.name = "Task08_Check_Initial_Conditions.R"
dirs$output.prefix= str_extract(dirs$Rscript.name,"^Task\\d\\d\\w?_")

#set up the model and dosing
model  = ivsc_2cmt_RR_KssT0L0()
#model1 = ivsc_2cmt_RR_KeqT0L0()

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
for (drug in drugs) {
  i=i+1
  filename = paste0("parameters/ModelG_",drug,"_Params.xlsx")
  param.as.double =  read.param.file(filename) 
  
  #reparameterize for model with T0, L0
    p    = param.as.double %>% t() %>% as.data.frame()
    if (p$keTL == 0) {
      TL0 = with(p,kon_TL*ksynT*ksynL/(koff_TL*keL*keT))
      T0  = with(p,ksynT/keT)
    } else {
      a = with(p,keTL^2)
      b = with(p,-(keTL) * (ksynT +ksynL) - (((koff_TL+keTL)/kon_TL) * keT *keL))
      c = with(p, ksynL*ksynT)
      
      TL0 <- ((-b) -sqrt((b^2)-4*a*c))/(2*a)
      T0 = with(p,(ksynT - keTL*TL0)/keT)
    }    
    L0 = with(p,(ksynL + koff_TL*TL0)/(kon_TL*T0 + keL))
  
    param.as.double["T0"] = T0
    param.as.double["L0"] = L0
    param.as.double["Kss_DT"] = with(p,(koff_DT + keDT)/kon_DT)
    param.as.double["Kss_TL"] = with(p,(koff_TL + keTL)/kon_TL)
    
    param.as.double = param.as.double[model$pin]
  
  #solve the model
    init = model$init(param.as.double)
    out  = model$rxode$solve(model$repar(param.as.double), ev, init)
    out  = model$rxout(out)

    out_long[[drug]] = out %>%
      dplyr::select(time,T,L,TL) %>%
      gather(key,value,-c(time)) %>%
      mutate(drug = drug,
             T0   = param.as.double["T0"],
             L0   = param.as.double["L0"])
}
out_long = bind_rows(out_long)
out1     = out_long %>%
  filter(!duplicated(drug))

g = ggplot(out_long,aes(x=time,y=value,color=key))
g = g + geom_line()
g = g + xgx_scale_y_log10()
g = g + facet_wrap(~drug)
g = g + geom_hline(data = out1, aes(yintercept = T0 ),color="purple",linetype="dashed",alpha=.5,size=2)
g = g + geom_hline(data = out1, aes(yintercept = L0 ),color="purple",linetype="dashed",alpha=.5,size=2)
g = g + ggtitle("Purple line is theory")
print(g)

ggsave("results/Task08b_init_check_reparT0L0.png",width=6,height=6)


