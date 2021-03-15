function [c_vec_Current,c_net_vec_Current] = cICs_CCC(epsilon,c_0,rho_0,D_prime,D_0,F,rho_net_vec_Current,PllPoints,PlsPoints)

% SET C IN THE SOLID PHASE REGION

c_net_vec_Current_Solid =  1/2/F*rho_net_vec_Current*(D_prime*(1-epsilon)/rho_0/D_0/epsilon^1.5);

% SET C IN THE LIQUID PHASE REGION

c_net_vec_Current_Liquid = transpose(linspace(c_net_vec_Current_Solid(end),0,PllPoints-PlsPoints));

% CONCATENATE VECTORS

c_net_vec_Current= vertcat(c_net_vec_Current_Solid,c_net_vec_Current_Liquid);

c_vec_Current = c_net_vec_Current + c_0;

end