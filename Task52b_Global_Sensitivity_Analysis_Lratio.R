# Setup and Read Data
source("ams_initialize_script.R")
source("SCIM_calculation.R")  
source("ivsc_2cmt_RR_V1.R")
library(RxODE)
model = ivsc_2cmt_RR_KdT0L0()
dirs$rscript_name = "Task52b_Global_Sensitivity_Analysis_Lratio.R"
dirs$filename_prefix= str_extract(dirs$rscript_name,"^Task\\d\\d\\w?_")

data_in = read.csv("results/Task22_2019-11-23_40e3.csv",stringsAsFactors = FALSE)
data    = data_in
 
#comute Lratio and plot
data = data %>%
  mutate(Lmax = Lss_thy,
         Lratio = (Lss_sim - L0)/(Lmax-L0),
         infusion = ifelse(infusion==1,"continual infusion","every 2-4 week dosing"),
         SSIM = 1-SCIM_sim)

#plot result ----
low = 1e-4
high = 1
g = ggplot(data,aes(x = Lratio, y = SSIM))
g = g + xgx_scale_x_log10(limits = c(low,high))
g = g + xgx_scale_y_log10(limits = c(low,high))
g = g + geom_point(alpha = .1)
g = g + facet_wrap(~infusion)
g = g + annotate("segment",x = low, y = low, xend = high, yend = high, color = "blue")
g = g + labs(x = "Lratio = [Lss(dose) - L0]/[Lmax - L0]",
             y = "SSIM")
#ggsave(width = 6, height= 3, filename = "./figures/Task52b_Lratio_Global_Sensitivity.png")
print(g)
