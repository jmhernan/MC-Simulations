# MC-Simulations
This repository is dedicated to my dissertation work where I used Monte Carlo Simulations to assess the performance of different 
Propensity Score Matching methods.  The study examined best approaches for comparing groups when randomization is not available (i.e., quasi-experimental designs).  I used Monte Carlo simulation to investigate multilevel data based on real-world parameters and compared commonly used Propensity Score Matching (PSM) methods with non-matching based methods.

The file "SimulationCode.R" is a sample of one of the many MC experiments that were conducted for the study.  
The file "SimFunctions.R" calls the 'getICC' function used for the estimation of the Intra-class correlation for each dataset and the r.cor.matrix function is used to create a data set with a given correlational structure. 

