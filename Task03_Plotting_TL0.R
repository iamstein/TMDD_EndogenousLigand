params <- read.csv(file ="task02a_sensitivity_all drugs and params_100_updated 04_24.csv",header = TRUE)

library(ggplot2)
library(scales)


data.plot <- all_params %>%
  select(fold.change.param, 
         #AFIR_thy, 
         #SCIM_thy_ketl_neg,
         #SCIM_thy_ketl_pos,
         #SCIM.sim, 
         #SCIM.thy,
         TL0.sim,
         #TL0.thy,
         TL0_neg.thy,
         #TL0_pos.thy,
         drug, param) %>%
  gather(key, value, -c(fold.change.param,drug,param))

g <- ggplot(data.plot, aes(x=fold.change.param, y=value, color=key, linetype=key)) + 
  facet_grid(drug ~ param, scales = "free_y", switch = "y") + 
  geom_line(size=1, alpha=0.6) +
  scale_x_log10() + 
  scale_y_log10() +
  scale_color_manual(values = c(#SCIM_thy_ketl_pos="orange",
                               TL0.sim="black",
                               TL0.thy="red",
                               TL0_neg.thy="green",
                               TL0_pos.thy="blue"))
print(g)