function [Kp_res] = Kp(g_C,c,g_D,d,g_A,a,g_B,b,T)
% Universal data
R_hat = 8.31447;          % [kJ kmol^-1 K^-1]

deltaG_hat = ((c * g_C) + (d * g_D)) - ((a * g_A) + (b * g_B));
Kp_res = exp(-deltaG_hat / (R_hat * T));

end