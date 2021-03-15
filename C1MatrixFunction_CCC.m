function [C_1,D_star] = C1MatrixFunction_CCC(tStep,llStep,D_0,epsilon,PlsPoints,PllPoints)

%% Initialise blank matrix

C_1 = zeros(PllPoints-PlsPoints-2,PllPoints-PlsPoints-2);

%% Compute D-star

D_star = D_0*epsilon^1.5*tStep/llStep/llStep;

%% Populate interior of C_1 matrix

% Inner values

for i = 2:PllPoints-PlsPoints-3
    
    C_1(i,i-1) = D_star; % sets above diagonals
    C_1(i,i) = 1 - 2*D_star; % sets diagonals
    C_1(i,i+1) = D_star; % sets below diagonals
    
end

% Boundary values

        C_1(1,1) = 1 - 2*D_star; % sets top corner
        C_1(1,2) = D_star; % sets top corner plus one to right
        C_1(end,end) = 1 - 2*D_star; % sets below top corner
        C_1(end,end-1) = D_star; % sets bottom corner plus one to left


end