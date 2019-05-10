params <- read.csv(file ="task02a_sensitivity_all drugs and params_100_updated 04_24.csv",header = TRUE)

library(ggplot2)
library(scales)


g <- ggplot(params, aes(x=fold.change.param)) + 
  facet_grid(drug ~ param,scales = "free_y") + 
  geom_line(aes(y = AFIR_thy),color = "green",size = 1) +
  geom_line(aes(y = SCIM_thy_ketl_neg), color = "red",size = 1) + 
  geom_line(aes(y = SCIM.sim),color = "blue",size = 1,linetype = "dashed")

g + scale_x_continuous(trans = 'log10') +annotation_logticks(sides = 'b') + scale_y_continuous(trans = 'log10') +
  ylab('Signaling complex inhibition metric \n (SCIM)')

params_dose <- read.csv(file = "Task02b_sens_SCIM_dose_all_drugs.csv",header =TRUE)

g2 <- ggplot(params_dose, aes(x = param.to.change*scale.nmol2mpk))+ facet_grid(~drug,scale = "free_y") +
  geom_line(aes(y = SCIM.thy), color = "red",size = 1) + geom_line(aes(y=SCIM.sim),color = "blue",size = 1,linetype = "dashed")     

g2 + scale_x_continuous(trans = 'log10') + annotation_logticks(sides = 'b') + scale_y_continuous(trans = 'log10') +
  ylab('Signaling complex inhibition metric \n (SCIM)') + xlab('Dose (mg/kg)')

plot(g2)


all_Res <- read.csv(file = "all_res.csv",header =TRUE)

ggplot(all_Res,aes(x = time, y=D)) + facet_grid(.~drug,scale = 'free_y') + geom_line(aes(y = D),color = 'red')

g <- ggplot(params, aes(x=AFIR_thy,y= SCIM_thy_ketl_neg)) +geom_point() +facet_grid(.~drug,scale = "free_y") + scale_y_continuous(trans = 'log10') +
  scale_x_continuous(trans = 'log10')
plot(g)

g <- ggplot(params, aes(x=AFIR_thy,y= SCIM_thy_ketl_pos)) +geom_point() +facet_grid(.~drug,scale = "free_y") + scale_y_continuous(trans = 'log10') +
  scale_x_continuous(trans = 'log10')
plot(g)