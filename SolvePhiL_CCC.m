function [phi_l_vec_Next] = SolvePhiL_CCC(phi_s_vec_Next,kappa_0,sigma,epsilon,PlsPoints,PllPoints)

% SOLVE FOR PHI_L IN THE SOLID PHASE REGION

phi_l_vec_Next_Solid = - sigma/kappa_0*phi_s_vec_Next;

% SOLVE FOR C IN THE LIQUID PHASE REGION

phi_l_vec_Next_Liquid = transpose(linspace(phi_l_vec_Next_Solid(end),0,PllPoints-PlsPoints));

% CONCATENATE VECTORS

phi_l_vec_Next = vertcat(phi_l_vec_Next_Solid,phi_l_vec_Next_Liquid );

end