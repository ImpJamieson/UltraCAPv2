%% Clear data for initialisation

close all
clear variables
clc

%% Input parameters

% TIME PARAMETERS

tBeg = 12.65; % time always starts at zero [s]
tEnd = 17.65; % simulation end time, equal to maximum permissible charge duration [s]
tSteps = (tEnd-tBeg)*100; % number of intervals for time to be computed at
tStep = (tEnd-tBeg)/tSteps; % time step between each calculation [s]
tVec = tBeg:tStep:tEnd; % temporal vector
tPoints = length(tVec); % number of time points to be computed

% GEOMETRICAL PARAMETERS

% Solid phase domain

lsBeg = 0; % solid phase length starts at zero [m]
lsEnd = 50e-6; % solid phase end length [m]
lsSteps = 500; % number of intervals for solid phase points to be computed at
lsStep = (lsEnd - lsBeg)/lsSteps; % spatial step between points [m]
lsVec = lsBeg:lsStep:lsEnd; % solid phase spatial vector
lsPoints = length(lsVec); % number of solid phase spatial points to be computed

% Liquid phase domain

llBeg = lsEnd; % liquid phase length starts at zero [m]
llEnd = 62.5e-6; % liquid phase end length [m]
llStep = lsStep;
llVec = llBeg:llStep:llEnd; % liquid phase spatial vector
llPoints = lsPoints + length(llVec); % number of liquid phase spatial points to be computed

% VARYING PHYSICAL PARAMETERS

density = 2100; % density of the conductor [kg/m^3]
molar_mass = 12.01115/1000; % molar mass of the conductor
Free_E_per_molecule = 1;  % free electrons per molecule
D_0 = 3.5e-11; % Diffusion coefficient (free solution) in liquid phase [m^{2}s^{?1}]
c_0 = 930; % base ionic concentration in liquid phase [mol m^{?3}]
T = 298.15; % ambient absolute temperature [K]
I_in = 100; % current supply in [A, C s^-1]
V_t = 0.5; % target voltage at terminal [V]
A_xc = 2.747; % electrode cross sectional area [m^{2}] chosen to match Verbrugge
lambda_s = 0.3e-09; % Stern layer thickness [m]
varepsilon_l = 36.6; % relative permittivity of solvent in liquid phase
a_d = 1.9696e09; % specific surface area [m^{2} m^{-3}]
epsilon = 0.67; % void fraction
V_start_raw = 1.5; % capacitor starting voltage (V)
V_end= 1.4; % discharge/charge voltage (V)
sigma = 0.0521; % electrode effective conductivity [S m^{-1}]

% CONSTANT PHYSICAL PARAMETERS

k_B = 1.380649e-23; % Boltzmann's constant [J K^{-1}]
N_A = 6.02214076e23; % Avogadro constant [mol^{-1}]
e = -1.6e-19; % electron elementary charge [C]
F = abs(N_A*e); % Faraday's constant [C/mol]
varepsilon_0 = 8.85418782e-12; % permittivity of free space [F m^{-1}]
varepsilon = 0.1750; % permittivity of electrode [F m^{-1}]

% CAPACITANCE MODEL CALCULATION

C_D = varepsilon_l*varepsilon_0/lambda_s; % Helmholtz prediction of capacitance
K = a_d*C_D; % effective specific capacitance

% BOUNDARY CONDITION SPECIFICATION

% Choose either 'ConstantCurrent' or 'ConstantVoltage'

connection_type = 'ConstantVoltage';

% DERIVED SIMULATION PARAMETERS

Free_E = density/molar_mass*Free_E_per_molecule*N_A; % free electron density of conductor [m^-3]
V_start= V_start_raw/2; % capacitor starting voltage (V)
j_in = I_in/A_xc; % current supply in [C s^-1 m^-2]
rho_0 = Free_E*e; % base level of charge concentration [C m^-3]
kappa_0 = 2*e^2*N_A/k_B/T*D_0*epsilon^1.5*c_0; % electrolyte conductivity [S m^{-1}] 
var_rho = -K*(1+sigma*(1-epsilon)^(3/2)/kappa_0)*V_end/2; % value of rho at edge under constant voltage conditions
Gamma = tStep/(lsStep)^2; % Fourier number
D_prime = sigma*k_B*T*(1-epsilon)^(1/2)/e; % diffusivity paramter collection of coefficients
Lambda = sigma*tStep*(1-epsilon)^(1/2)/varepsilon; % dimensionless number
Upsilon = sigma*k_B*T*(1-epsilon)^(1/2)/e*Gamma; % collection of values
xi = varepsilon_0/(lsStep)^2; % collection of values
alpha = a_d*C_D/2/abs(e)/N_A/epsilon;
beta = D_0*Gamma/epsilon;
tStep_Critical = varepsilon_0/sigma; % maximum time step for numerical stability
C_total = a_d*C_D*A_xc*lsEnd; % total effective capacitance

rho_t = var_rho+rho_0; % TEST VALUE FOR NOW, NEEDS TO BE CHANGED FOR CONSTANT VOLTAGE
Theta = sigma*k_B*T/e/rho_t*Gamma; % collection of values

% OTHER DERVIED PARAMETERS

R = 1/sigma*lsEnd; % effective resistance (Ohm)
C = varepsilon/lsEnd; % effective capacitance (Farad)
TheoreticalTimeConstant = R*C; % effective time constant (s)
ComputationalTimeConstant = varepsilon/sigma; % computational time constant (s)

%% Preparations to save output data

tLinesNumber = 10; % number of time evolution lines to plot
tInterval = (tPoints-1)/tLinesNumber; % time intervals to plot at

for i = 1:tLinesNumber+1
    
    if i ==1
        
        tEvolutionPoints = [];
        
    end

tSelect = round(i*tInterval - (tInterval-1));
tEvolutionPoints = vertcat(tEvolutionPoints,tSelect);

end

phi_s_mat = zeros(length(tEvolutionPoints),length(lsVec));
rho_mat = zeros(length(tEvolutionPoints),length(lsVec));
ls_mat = zeros(length(tEvolutionPoints),length(lsVec));
phi_l_mat = zeros(length(tEvolutionPoints),length(lsVec)+length(llVec));
c_mat = zeros(length(tEvolutionPoints),length(lsVec)+length(llVec));
ll_mat = zeros(length(tEvolutionPoints),length(lsVec)+length(llVec));

switch connection_type
   case 'ConstantCurrent'
       RlsPoints = lsPoints; % trims number of length points for solid phase rho matrices
       PlsPoints = lsPoints-1; % trims number of length points for solid phase phi matrices
       PllPoints = llPoints-1; % trims number of length points for solid phase phi matrices
       var_phi1 = 0; % sets end value of phi1 to be 0
   case 'ConstantVoltage'
       RlsPoints = lsPoints-1; % trims number of length points for solid phase rho matrices
       PlsPoints = lsPoints-1; % trims number of length points for solid phase phi matrices
       PllPoints = llPoints-1; % trims number of length points for solid phase phi matrices
       var_phi1 = 0; % sets end value of phi1 to be 0
   otherwise
     disp('ERROR: Specified connection type is unknown.')
     disp('Cannot correctly trim matrices.')
     error('An error occured. Simulation halted. See the above messages for information.')
end

%% Initial conditions

% CALL INITIAL CONDITION FUNCTIONS

load CCCValidation.mat  rho_net_vec_Current phi_s_vec_Current phi_l_vec_Current c_net_vec_Current

[rho_vec_Current,rho_net_vec_Current] = rhoICs_CVD(rho_0,rho_net_vec_Current); % creates rho initial conditions
[c_vec_Current,c_net_vec_Current] = cICs_CVD(c_0,c_net_vec_Current); % creates c initial conditions
[phi_s_vec_Current] = phi_sICs_CVD(phi_s_vec_Current); % creates phi_s initial conditions
[phi_l_vec_Current] = phi_lICs_CVD(phi_l_vec_Current); % creates phi_l initial conditions
D_vec_Current = D_prime./rho_vec_Current; % sets initial condition for the diffusivity parameter at all x
Sigma_vec = D_vec_Current*Gamma;
kappa_vec = 2*e^2*N_A/k_B/T*D_0*epsilon^1.5*c_vec_Current; % vector of intitial values of kappa collection of values
varkappa_vec = kappa_vec*Gamma/a_d/C_D; % initial values of varkappa collection of values
v_terminal = 2*phi_s_vec_Current(1); % initial value of v_terminal

err = 1; % sets initial error (dummy value to allow computation to start)
tol = 1e-27; % sets acceptable tolerance limit
j=1; % sets up iteration counter
t=0; % sets up time counting

%% Display information in command window

disp('Time step [s] is:')
disp(tStep)
disp('Dimensionless Lambda is:')
disp(Lambda)
disp('Maximum dimensionless varkappa is:')
disp(max(varkappa_vec))
disp('Dimensionless Beta is:')
disp(beta)

if Lambda > 1
    disp('Lambda stability criterion not satisfied.')
    disp('Reduce time step [s] to below')
    disp(tStep_Critical)
else
    disp('Lambda stability criterion satisfied.')
end 

disp('UltraCAP initialisation complete.')

%% Main computation

% CALL FUNCTIONS TO GENERATE TIME-INVARIANT MATRIX AND VECTOR FUNCTIONS

[P_1_inv] = P1MatrixFunction_CVD(Lambda,PlsPoints);
[P_2] = P2MatrixFunction_CVD(PlsPoints);
[C_1,D_star] = C1MatrixFunction_CVD(tStep,llStep,D_0,epsilon,PlsPoints,PllPoints);

while V_t >= v_terminal  % DUMMY FOR NOW
    
% DETERMINE VALUE OF LAMBDA AT THIS TIME STEP

lambda = rho_net_vec_Current(end);
lambda = 0; % TEST VALUE

% CALL FUNCTIONS TO GENERATE TIME-VARAIANT MATRIX AND VECTOR FUNCTIONS

[R_1] = R1MatrixFunction_CVD(Sigma_vec,Theta,e,epsilon,sigma,k_B,T,j_in,lsStep,RlsPoints,connection_type);
[A_1] = A1VectorFunction_CVD(Sigma_vec,rho_t,Upsilon,Theta,RlsPoints,connection_type);

% SOLVE FOR SOLID PHASE ELECTRIC CHARGE DENSITY AT NEXT TIME LEVEL

[rho_vec_Next,rho_net_vec_Next] = SolveRho_CVD(rho_vec_Current,R_1,A_1,rho_0,rho_t,connection_type);

% SOLVE FOR LIQUID PHASE IONIC CONCENTRATION AT NEXT TIME LEVEL

[c_vec_Next,c_net_vec_Next ] = SolveC_CVD(C_1,D_star,c_vec_Current,epsilon,c_0,rho_0,D_prime,D_0,F,rho_net_vec_Next,PlsPoints,PllPoints,lsStep);
 
 % SOLVE FOR SOLID PHASE ELECTRIC POTENTIAL AT NEXT TIME LEVEL

[phi_s_vec_Next] = SolvePhiS_CVD(rho_net_vec_Next,K,sigma,kappa_0,epsilon,PlsPoints); 
 
 % SOLVE FOR LIQUID PHASE ELECTRIC POTENTIAL AT NEXT TIME LEVEL

[phi_l_vec_Next] = SolvePhiL_CVD(phi_s_vec_Next,kappa_0,sigma,epsilon,PlsPoints,PllPoints);

% RECORD AND SAVE TERMINAL VOLTAGE, CHARGE AND CURRENT INFORMATION

if j == 1
    I_vec = [];
    j_vec = [];
    Q_vec = [];
    rho_edge_vec = [];
end 

j_terminal = D_prime./rho_vec_Next(1)/lsStep*(- 1/2*rho_vec_Next(3) + 2*rho_vec_Next(2) - 3/2*rho_vec_Next(1)); % records estimate of specific current at current collector 
I_terminal = j_terminal*A_xc; % records estimate of current at current collector
rho_terminal = rho_net_vec_Next(1); % records net charge density at current collector

j_vec = horzcat(j_vec,j_terminal);
I_vec = horzcat(I_vec,I_terminal);
rho_edge_vec = horzcat(rho_edge_vec,rho_terminal);


if j == 1
    v_vec = [];
end 

v_terminal = 2*phi_s_vec_Next(1); % records potential at current collector 

v_vec = horzcat(v_vec,v_terminal);

% SAVE OUTPUT DATA AT REGULAR TIME INTERVALS

tMatch = find(j == tEvolutionPoints);

if isempty(tMatch) == 0

phi_s_mat(tMatch,:) = phi_s_vec_Next;
rho_mat(tMatch,:) = rho_net_vec_Next;
ls_mat(tMatch,:) = lsVec;
phi_l_mat(tMatch,:) = phi_l_vec_Next;
c_mat(tMatch,:) = c_vec_Next;
ll_mat(tMatch,:) = horzcat(lsVec,llVec);

end

j = j + 1; % update step count
t = t + tStep; % update time count
err = rho_vec_Current(1) - rho_vec_Next(1); % provisional error calculation
rho_vec_Current = rho_vec_Next;
rho_net_vec_Current = rho_net_vec_Next;
phi_s_vec_Current = phi_s_vec_Next;
phi_l_vec_Current = phi_l_vec_Next;
c_vec_Current = c_vec_Next;
D_vec_Current = D_prime./rho_vec_Current;
Sigma_vec = D_vec_Current*Gamma;

if t > tEnd-tBeg
        disp('Simulation time limit reached.')
        break
end

end

%% Test output data

tVec = tVec(1:length(I_vec)); % trims time vector to correct length

%figure(9)
%plot(tVec,I_vec,'k');
%title('Terminal current')
%xlabel('Time ,s') 
%ylabel('Current, A')

vVecAbs = abs(v_vec);

%figure(10)
%plot(tVec,vVecAbs,'k');
%title('Terminal potential')
%xlabel('Time, s') 
%ylabel('Potential, V')

save CVDValidation.mat  rho_mat phi_s_mat phi_l_mat c_mat ls_mat ll_mat tVec vVecAbs tVec I_vec