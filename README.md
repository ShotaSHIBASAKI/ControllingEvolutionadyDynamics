# ControllingEvolutionadyDynamics
Python codes and csv files used in Shibasaki and Mitri 2019 (*bioRxiv*) "Controlling evolutionary dynamics to optimize microbial bioremedation".
You can replicate the simulation, figures, and calculations with the codes. 
There are also some csv files used drowing some figures.

## FigsMainText.py
This code includes the functions to drow the figures in the main text (except for Figs. 1A, 2A and E, 3, and 5A).
For analysis of Markov chian models in Figs. 5B (ergodic one) or A.11A (non-ergodic one), you can also see MarkovChain.c. The two csv files, "optDPMC_ergodics.csv", and "optDPMC_Nonergodic.csv" provide the optimal intruction rates of cooperators in each model given a time step, respectively. 

## Appendix1.py
This code drows Fig. A.1(the phase plain analysis) when there exist only one strategy in a chemostat.

## Appendix2.py
This codes shows the calcualtions in section A.2 (invasion analysis), and drows Figs.2 and 3.

## Appendix3.py
This codes analyzes the 3x3 Jacobian matrices in section A.3 (linear stability of equilibria), and drows Figs. A6 (evolutionary dynamics when two species exist), and A7 (Monte Calro simulation). The summary of the results of Monte Carlo simulation is available in "parameter_sweep_total_case0.csv" (rCo and sCh) or "parameter_sweep_total_case1.csv" (sCo and rCh). The raw data is available on request (as the file sizes are too large).

## Appendix4.py
This codes analyzes the 4x4 Jacobian matrix in section A.4.

## Appendix5.py
This code shows the dynamics when the mutation rates are not zero, as in section A.5.

## Appendix7.py
This code calculates the stationary distribution as in Eq (A.64).
