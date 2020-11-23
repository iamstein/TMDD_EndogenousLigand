# this is a parent script that calls all scripts that create all plots and tables for our manuscript

#the function below takes a very long time to run and is commented out.
#it performs the ~40,000 simulations needed for the global sensitivity analysis.  plotting is done by another function
#source("Task22_LatinHypercube_Soluble_and_Membrane_AddCL_AddDoseReg.R") #this function

#create plots
source("Task50_SingleDose_Figure_RealDrugs.R") #plot simulations multiple dosing for 4 drugs studied 
source("Task51_MultiDose_Figure_RealDrugs.R") #plot dose vs ASIR for multiple dosing at steady state for 4 drugs studied
#source("Task52_Global_Sensitivity_Analysis.R") #historgrams of global sensitivity analysis
source("Task52b_Global_Sensitivity_Analysis_Lratio.R") #Lratio scatter plot (theory vs sim) from global sens analysis sims
source("Task52c_Global_Sensitivity_Analysis_PosNeg.R")
source("Task53_Dose_Range_Toci.R") #Simulate Tociliszumab: Drug, Ligand (IL-6), and SSIM
source("Task60_Sensitivity_6Drug.R") #Local sensitivity analysis for 4 drugs

#create tables
source("Task54_Parameter_Table_6Drugs.R") #Create parameter tables used for the 4 drugs studied
source("Task55_Parameter_Range_Sensitivity.R") #Summarize the parameter ranges used in the global sensitivity analysis