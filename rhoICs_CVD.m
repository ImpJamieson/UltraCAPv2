function [rho_vec_Current,rho_net_vec_Current] = rhoICs_CVD(rho_0,rho_net_vec_Current)

% Load data

rho_vec_Current = rho_net_vec_Current + rho_0; 

end