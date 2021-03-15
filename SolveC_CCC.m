function [c_vec_Next,c_net_vec_Next ] = SolveC_CCC(C_1,D_star,c_vec_Current,epsilon,c_0,rho_0,D_prime,D_0,F,rho_net_vec_Next,PlsPoints,PllPoints,lsStep)

% SOLVE FOR C IN THE SOLID PHASE REGION

c_net_vec_Next_Solid =  1/2/F*rho_net_vec_Next*(D_prime*(1-epsilon)/rho_0/D_0/epsilon^1.5);

% SOLVE FOR C IN THE LIQUID PHASE REGION
c_net_vec_Next_Liquid = transpose(linspace(c_net_vec_Next_Solid(end),0,PllPoints-PlsPoints));

%C_2 = zeros(1,PllPoints-PlsPoints-2);
%C_2(1) = D_star*(c_0+c_net_vec_Next_Solid(end));
%C_2(end) = D_star*c_0;

%c_vec_Next_Liquid = C_1*c_vec_Current(PlsPoints+3:PllPoints) + transpose(C_2);

% CONCATENATE VECTORS

c_net_vec_Next = vertcat(c_net_vec_Next_Solid,c_net_vec_Next_Liquid);

c_vec_Next = c_net_vec_Next + c_0;

end