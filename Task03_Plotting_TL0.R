source("ams_initialize_script.R")

params <- read.csv(file ="task02a_sensitivity_all drugs and params_100_updated 04_24.csv",header = TRUE)



data.plot <- params %>%
  dplyr::select(fold.change.param, 
         TL0_sim,
         TL0_negroot_thy,
         drug, param) %>%
  gather(key, value, -c(fold.change.param,drug,param))

g <- ggplot(data.plot, aes(x=fold.change.param, y=value, color=key, linetype=key)) + 
  facet_grid(drug ~ param, scales = "free_y", switch = "y") + 
  geom_line(size=1, alpha=0.6) +
  scale_x_log10() + 
  scale_y_log10() +
  scale_color_manual(values = c(#SCIM_thy_ketl_pos="orange",
                               TL0_sim="black",
                               TL0_negroot_thy="green"))
print(g)