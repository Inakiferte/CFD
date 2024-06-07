function [g_hat] = gibbs_free_energy(h_fhat,s_hat,a0,a1,a2,a3,a4,T)
addpath('./cp/');
T0 = 25 + 273;   % [K]

int1 = (cp_hat_int(T,a0,a1,a2,a3,a4) - cp_hat_int(T0,a0,a1,a2,a3,a4));
int2 = (cp_hat_intT(T,a0,a1,a2,a3,a4) - cp_hat_intT(T0,a0,a1,a2,a3,a4));
g_hat = h_fhat + int1 - T * (s_hat + int2);
end