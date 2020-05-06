source("ams_initialize_script.R")
source("ivsc_2cmt_RR_V1.R")
model = ivsc_2cmt_RR_KdT0L0()

# Create paths to data files for each drug.
param = NULL
for (drug in drugs) {
  param[[drug]] = read_excel(parameter_files[[drug]], 1) %>%
    mutate(Drug = drug,
           Value = as.numeric(Value)) %>%
    select(Order,Drug,ParamType,Molecule,Description,Parameter,Value,Units) %>%
    filter(Parameter %in% model$pin) %>%
    filter(!(Parameter %in% c("F","ka")))
}
d = bind_rows(param) %>%
  select(Drug,Parameter,Value) %>%
  mutate(Value = as.character(signif(Value,2))) %>%
  spread(Drug,Value)

units = param[[1]] %>%
  select(Parameter,Units)
d = left_join(d,units,by="Parameter") %>% 
  select(Parameter,Units,everything()) %>%
  filter(!(Parameter %in% c("Vm","Km")))

#reorder parameters
order = c("V1" ,"k12"   ,"k21"   ,"keD" ,"L0" ,"T0" ,"keT" ,"keDT" ,"keL" ,"keTL" ,"Kd_DT","Kd_TL","kon_DT","kon_TL")
latex = c("V_1","k_{12}","k_{21}","\\keD","L_0","T_0","\\keT","\\keDT","\\keL","\\keTL","\\KdDT","\\KdTL","\\konDT","\\konTL")
d = d %>%
  mutate(order_temp = factor(Parameter, levels = order)) %>%
  arrange(order_temp) %>%
  select(-order_temp) %>%
  mutate(Parameter = plyr::mapvalues(Parameter,order,latex)) %>%
  mutate(Parameter = paste0("$",Parameter,"$"))

file.prefix = paste0("figures/Task54_ParamTable")
#readr::write_csv(d,paste0(file.prefix,".csv")) #to csv file
xtab = xtable::xtable(d) #to latex file
print(xtab,
      file=paste0(file.prefix,".tex"),
      include.rownames=FALSE,
      floating = FALSE,
      sanitize.text.function = function(x){return(x)})
print(kable(d))
