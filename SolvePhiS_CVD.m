function [phi_s_vec_Next] = SolvePhiS_CVD(rho_net_vec_Next,K,sigma,kappa_0,epsilon,PlsPoints)

% Use specific capacitance equation

phi_s_vec_Next = rho_net_vec_Next/K /(1+sigma/kappa_0);

end