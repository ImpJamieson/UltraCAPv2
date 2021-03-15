function [phi_l_vec_Current] = phi_lICs_CCC(phi_s_vec_Current,kappa_0,sigma,epsilon,PllPoints,PlsPoints)

% SOLVE FOR PHI_L IN THE SOLID PHASE REGION

phi_l_vec_Current_Solid = - sigma/kappa_0*phi_s_vec_Current;

% SOLVE FOR C IN THE LIQUID PHASE REGION

phi_l_vec_Current_Liquid = transpose(linspace(phi_l_vec_Current_Solid(end),0,PllPoints-PlsPoints));

% CONCATENATE VECTORS

phi_l_vec_Current = vertcat(phi_l_vec_Current_Solid,phi_l_vec_Current_Liquid );

end