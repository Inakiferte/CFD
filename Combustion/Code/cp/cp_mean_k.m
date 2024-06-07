function [cp_mean] = cp_mean_k(Tcc,T0,alpha,beta,gamma,delta,epsilon)

cp_mean = (1 / (Tcc - T0)) * (cp_hat_int(Tcc,alpha,beta,gamma,delta,epsilon) - cp_hat_int(T0,alpha,beta,gamma,delta,epsilon));

end