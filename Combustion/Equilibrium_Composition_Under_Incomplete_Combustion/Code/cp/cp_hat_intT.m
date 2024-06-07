function [cp_hat] = cp_hat_intT(T,alpha,beta,gamma,delta,epsilon)
% In this function the integral of cp is done. The used data is for air

% Air data
R_hat = 8.31447;          % [kJ kmol^-1 K^-1]

cp_hat = ((alpha * log(T)) + (beta * T) + (gamma / 2 * T^2) + (delta / 3 * T^3) + (epsilon / 4 * T^4)) * R_hat;

end