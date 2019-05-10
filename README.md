# TMDD_EndogenousLigand
Novartis Hackathon 2018 - TMDD model with endogenous ligand

## Folder organization
parameters/ - contains parameter files and data
results/    - where results should be stored.

## Code organization
* ivsc_2cmt_RR_V1.R - contains the ODE and the initial conditions
* SCIM_Calculation (or SCIM_calculation_2019_05_10.R) contains a lot of the key code for comparing theory and simulation
  * lumped.parameters.simulation() - calls the ODE in ivsc_2cmt_RR_v1 for a particular set of parameters and doses and computes the SCIM parameters
  * lumped.parameters.theory() - performs the theoretical calculation for SCIM, AFIR, etc., using only the model parameters.
  * compare.thy.sim() - puts the theory and the simulation into a very large table, for use in plotting.
* Task01, Task02, etc.
  * These call compare.thy.sim(), and then should ultimately contain some plotting code.


