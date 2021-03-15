function [phi_s_vec_Current] = phi_sICs_CCC(rho_net_vec_Current,K,sigma,kappa_0,epsilon)

% Use specific capacitance equation

phi_s_vec_Current = rho_net_vec_Current/K/(1+sigma/kappa_0);

end