function [c_vec_Next,c_net_vec_Next ] = SolveC_CVD(C_1,D_star,c_vec_Current,epsilon,c_0,rho_0,D_prime,D_0,F,rho_net_vec_Next,PlsPoints,PllPoints,lsStep)

% SOLVE FOR C IN THE SOLID PHASE REGION

c_net_vec_Next_Solid =  1/2/F*rho_net_vec_Next*(D_prime*(1-epsilon)/rho_0/D_0/epsilon^1.5);

% SOLVE FOR C IN THE LIQUID PHASE REGION
c_net_vec_Next_Liquid = transpose(linspace(c_net_vec_Next_Solid(end),0,PllPoints-PlsPoints));

% CONCATENATE VECTORS

c_net_vec_Next = vertcat(c_net_vec_Next_Solid,c_net_vec_Next_Liquid);

c_vec_Next = c_net_vec_Next + c_0;

end