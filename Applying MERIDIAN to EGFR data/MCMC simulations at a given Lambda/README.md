The scripts in this directory will allow one to sample the signaling network parameters using MCMC. Here, I briefly describe individual scripts

all_derivatives.m: Describes the rates of change of different chemical species in the EGFR/Akt model

binProp.m: assigns histogram bin fractions to data

decide_which_bin_akt.m: Assigns a particular bin value to predicted pAkt levels when given the bin locations

decide_which_bin_egfr.m: Assigns a particular bin value to predicted sEGFR levels when given the bin locations

evaluate_energy.m: Evaluates the "energy" = -log P(\theta) to perform MCMC

extract_preds.m: Extract pAkt and sEGFR predictions from solution to the ODE

initial_conditions.m: Initializes initial conditions

is_within_cutoffs.m: checks that pAkt and sEGFR values are within predefined cutoffs

lmbs_notail.mat: A set of Lambda values

modelPreds.m: Solves differential equations with given initial conditions and Ligand concentrations

monte_carlo_step.m:	Performs one Monte Carlo step

run_simulation.m: Performs several monte carlo steps
