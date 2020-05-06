source("ams_initialize_script.R")
source("ivsc_2cmt_RR_V1.R")
model = ivsc_2cmt_RR_KdT0L0()

# Create paths to data files for each drug.
param = read_excel("parameters/Task22_Param_Ranges_AFIR_Tfold.xlsx") %>%
  filter(type_sol == 0 | type_mem == 0) %>%
  select(Parameter = ParameterTex, Units, min_sol, max_sol, min_mem, max_mem) %>%
  mutate(Parameter = paste0("$",Parameter,"$"),
         min_sol   = as.character(signif(min_sol,2)),
         max_sol   = as.character(signif(max_sol,2)),
         min_mem   = as.character(signif(min_mem,2)),
         max_mem   = as.character(signif(max_mem,2))) %>%
  rename(`Sol. Min.` = min_sol,
         `Sol. Max.` = max_sol,
         `Mem. Min.` = min_mem,
         `Mem. Max.` = max_mem)

file.prefix = paste0("figures/Task55_Param_Range_Sensitivity")
#readr::write_csv(d,paste0(file.prefix,".csv")) #to csv file
xtab = xtable::xtable(param)
print(xtab,
      file=paste0(file.prefix,".tex"),
      include.rownames=FALSE,
      floating = FALSE,
      sanitize.text.function = function(x){return(x)})
print(kable(param))
