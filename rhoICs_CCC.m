function [rho_vec_Current,rho_net_vec_Current] = rhoICs_CCC(var_rho,K,sigma,kappa_0,V_start,epsilon,lsPoints)

rho_add = -V_start*K*(1+sigma/kappa_0); 
rho_vec_Current = ones(lsPoints,1)*var_rho + rho_add;
rho_net_vec_Current = rho_vec_Current - var_rho; % creates net charge density vector

end