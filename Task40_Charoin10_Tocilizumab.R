source("ams_initialize_script.R")
source("SCIM_calculation.R")  
source("ivsc_2cmt_RR_V1.R")
dirs$rscript_name = "Task40_Charoin10_Tocilizumab.R"
dirs$filename_prefix= str_extract(dirs$rscript_name,"^Task\\d\\d\\w?_")

model = ivsc_2cmt_RR_KdT0L0()

# read in the data ----
data = read_excel("data/charoin10_ugml_nM.xlsx") %>%
  mutate(species_units = paste0(species, " (", valunit, ")")) %>%
  mutate(species = plyr::mapvalues(species,
                                   c("tocilizumab"         ,"IL-6"         ,"sIL-6R"),
                                   c("free tocilizumab (D)","free IL-6 (L)","soluble IL-6R (-)")))

# create plot like in Charoin ----
g = ggplot(data,aes(x=time,y=value, group = species_units, color = species_units))
g = g + geom_line()
g = g + geom_point()
g = g + xgx_scale_x_time_units("day")
g_charoin = g
print(g)

# create plot transformed to common units ----
g = g_charoin
g = g %+% filter(data, species!="CRP")
g = g + aes(y = val_nM, group = species, color = species)
g = g + xgx_scale_y_log10(breaks = 10^c(-10:10))
g = g + labs(y = "Concentration (nM)")
g_data = g
print(g)

# set up and run the simulation ----
tmax = 63+21 #days
compartment = 2
dose_mpk = 10
dose_nmol = dose_mpk*scale.mpk2nmol

ev = eventTable(amount.units="nmol", time.units="days")
sample.points = seq(0, tmax, 0.1) # sample time, increment by 0.1
ev$add.sampling(sample.points)
ev$add.dosing(dose=dose_nmol, start.time = 21, dosing.to=compartment)
 
param_as_double = read.param.file("parameters/ModelG_Tocilizumab_Params_Charoin10.xlsx")[model$pin]

init = model$init(param_as_double)
sim  = model$rxode$solve(model$repar(param_as_double), ev, init, atol = 1e-12, rtol = 1e-12)
sim  = model$rxout(sim)

# final plot ----
sim_plot = sim %>%
  select(time,
         `free tocilizumab (D)` = D, 
         `free IL-6 (L)` = L, 
         `free mIL-6R (T)` = T, 
         `drug-target complex (DT)` = DT,
         `signaling complex (TL)` = TL,) %>%
  gather(species,val_nM,-time) %>%
  mutate(time = time - 14) %>%
  filter(val_nM>0)

sim_last = sim_plot %>%
  filter(time == 56) %>%
  bind_rows(data.frame(
    time = max(sim$time),
    species = "soluble IL-6R (-)",
    val_nM  = 1.8)) %>%
  arrange(desc(val_nM))

sim_plot = sim_plot %>%
  bind_rows(data.frame(
    time = max(sim$time),
    species = "soluble IL-6R (-)",
    val_nM  = 2e0)) %>%
  mutate(species = factor(species, levels = unique(sim_last$species)))

data_plot = data %>%
  filter(species != "CRP") %>%
  mutate(species = factor(species, levels = unique(sim_last$species)))


data_plot$time = data_plot$time - 7
sim_plot$time  = sim_plot$time - 7

sim_last = sim_plot %>%
  filter(time == 63) %>%
  mutate(time = time + 1.5,
         cmt = str_extract(species,"[A-Z]+\\)$"),
         cmt = str_replace(cmt,"\\)",""))


# create plot like in Charoin ----
g = ggplot(sim_plot,aes(x=time,y=val_nM, 
                        group = species, color = species, shape = species, linetype = species))
g = g + geom_line()
g = g + geom_point(data = data_plot)
g = g + geom_text(data = sim_last, aes(label = cmt), show.legend = FALSE, hjust=0)
g = g + xgx_scale_x_time_units("day")
g = g + xgx_scale_y_log10(breaks = 10^(-10:10))
g = g + labs(y = "Concentration (nM)")
g = g + scale_shape_manual(values = c(16,17,46,46,15,46))
g = g + scale_linetype_manual(values = c("solid","blank","solid","solid","solid","solid"))
print(g)
ggsave(width = 6.5, height= 3.5, filename = "./figures/Task40_Tocilizumab.png")

