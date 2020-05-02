# This is the code used for a publication on the TMDD+L model that includes drug, target and endogenous ligand
This work was started as part of the 2018 Novartis-Academia Quantitative Hackathon

## Folder organization

* data/       - contains data from Charoin presentation for Tocilizumab
* parameters/ - contains parameter files and data.  numbers for the IL-6/IL-6R and PD-1/PD-L1 are what was used in the paper.  Numbers for VEGF/VEGFR are more questionable and we are not confident they are all correct.  We could not find enough data for the VEGF/VEGFR drugs.  Therefore, they are not referenced in the paper.  But we keep our work here in case it is later helpful as a starting point.
* results/    - where results of some scripts were stored.  These results are mostly old and not used in the manuscript
* figures/    - contains figures and tables used in the manuscript

## Code organization

* Task00_Create_All_Plots_And_Tables.R - creates all plots and tables from the manuscript.  It makes use of the scripts below.
* ams_initialize_script.R - called at the top of every R script
* ivsc_2cmt_RR_V1.R - contains the ODE and the initial conditions for the system we have studied.
* SCIM_Calculation - contains the key code for calculating the theoretical and simulated ASIR and SSIM
  * lumped.parameters.simulation() - calls the ODE in ivsc_2cmt_RR_v1 for a particular set of parameters and doses and computes the SCIM parameters
  * lumped.parameters.theory() - performs the theoretical calculation for SCIM, AFIR, etc., using only the model parameters.
  * compare.thy.sim() - puts the theory and the simulation into a very large table, for use in plotting.



