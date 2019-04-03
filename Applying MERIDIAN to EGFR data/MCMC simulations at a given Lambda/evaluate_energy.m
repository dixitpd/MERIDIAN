function [ en ] = evaluate_energy( lambda_akt,lambda_egfr,bins_akt,bins_egfr,dat_akt,dat_egfr )

% This script evaluates the energy used in MCMC sampling

bI_akt  = decide_which_bin_akt( bins_akt,dat_akt );
bI_egfr = decide_which_bin_egfr( bins_egfr,dat_egfr );

en = sum(sum(lambda_akt.*bI_akt)) + sum(sum(lambda_egfr.*bI_egfr));



end

