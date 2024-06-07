function [cp_hat] = cp_hat_int(T,alpha,beta,gamma,delta,epsilon)
% In this function the integral of cp is done. The used data is for air

% Air data
R_hat = 8.31447;          % [kJ kmol^-1 K^-1]

cp_hat = ((alpha * T) + (beta / 2 * T^2) + (gamma / 3 * T^3) + (delta / 4 * T^4) + (epsilon / 5 * T^5)) * R_hat;

end