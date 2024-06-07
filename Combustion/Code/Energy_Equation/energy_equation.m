function [Tcc] = energy_equation(Tcc_s,T0,n_new,sumR)
addpath('./cp/');
[h_fhat_N2,s_hat_N2,a0_N2,a1_N2,a2_N2,a3_N2,a4_N2] = N2_properties;
[h_fhat_CO,s_hat_CO,a0_CO,a1_CO,a2_CO,a3_CO,a4_CO] = CO_properties;
[h_fhat_CO2,s_hat_CO2,a0_CO2,a1_CO2,a2_CO2,a3_CO2,a4_CO2] = CO2_properties;
[h_fhat_H2,s_hat_H2,a0_H2,a1_H2,a2_H2,a3_H2,a4_H2] = H2_properties;
[h_fhat_H2O,s_hat_H2O,a0_H2O,a1_H2O,a2_H2O,a3_H2O,a4_H2O] = H2O_properties;

% Assume adiabatic CC
Qdot_lost = 0;

% Compute Products enthalpies of formation sum
sumPhf = (n_new(1) * h_fhat_CO2) + (n_new(2) * h_fhat_H2O) + (n_new(3) * h_fhat_CO) + (n_new(4) * h_fhat_H2) + (n_new(5) * h_fhat_N2);  

% Compute Products integrals sum
int_CO2 = cp_mean_k(Tcc_s,T0,a0_CO2,a1_CO2,a2_CO2,a3_CO2,a4_CO2);
int_H2O = cp_mean_k(Tcc_s,T0,a0_H2O,a1_H2O,a2_H2O,a3_H2O,a4_H2O);
int_CO  = cp_mean_k(Tcc_s,T0,a0_CO,a1_CO,a2_CO,a3_CO,a4_CO);
int_H2  = cp_mean_k(Tcc_s,T0,a0_H2,a1_H2,a2_H2,a3_H2,a4_H2);
int_N2  = cp_mean_k(Tcc_s,T0,a0_N2,a1_N2,a2_N2,a3_N2,a4_N2);
sumP    = (n_new(1) * int_CO2) + (n_new(2) * int_H2O) + (n_new(3) * int_CO) + (n_new(4) * int_H2) + (n_new(5) * int_N2);

% Compute combustion temperature
Tcc = T0 + ((sumR - sumPhf - Qdot_lost) / sumP);

end