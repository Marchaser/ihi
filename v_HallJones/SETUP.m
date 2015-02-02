function params = SETUP()
% Inferring Health Inequality
% Setup environments and parameters
% Author: Wenlan Luo at Georgetown University, luowenlan at gmail.com

%% Load library
run('./../SET_PATH.m')
set_dpopt;
% run CMEX file to preload memory
pp = struct('form','MKLpp','breaks',{{[1 2 3 4]}},...
    'Values',[1 2 3 4],'coefs',[],'order',[4],...
    'Method',[],'ExtrapolationOrder',[],'thread',1,...
    'orient','curvefit');
warmUp_pp = myppual(pp);

%% Set parameters
% demographics
Jw = 45;
Jr = 55;
J = Jw+Jr;

% utility
Dbeta = 1.0200; % discount factor
Dlambda = 1; % consumption-leisure weight
Dzeta = 1; % consumption weight
Drho = 3; % elasticit between consumption and health
Dnu = 3; % risk aversion
Dd = 4.625; % flow of being alive

% production
Dalpha = 0.36;
Ddelta = 0.0833;

% survival
Dpsi1 = 1; % survival probability level
Dpsi2 = 1.61875; % survival probability curvature

% health shock
Dhdelta0 = 1.90735e-09;
Dhdelta1 = 0.0002;
Dhdelta2 = 2.25e-06;

% Dhdelta0 = 0;
% Dhdelta1 = 0;
% Dhdelta2 = 0;

% Dhdelta = 0.0; % health depreciation
Dmu_z0 = -3.1; % level
% Dmu_z0 = -inf; % level
Dmu_z1 = 0.01275; % age slope
Dmu_z2 = 5e-6; % age quadratic
Dvar_z = 0.1744; % variance of disease shock
% Dh1 = (Dpsi2-0.998).^(-1/Dpsi1);
Dcata = 0.00;   
Ddisease = 1; % probability of disease shock

% health technology
Dtheta = 0.25; % curvature
DA_h = 1.69688e-1; % multiplier

% tax
% Dtau_ss = 0.124; % social security tax
% Dtau0_gs = 0.258; % Gouveia and Strauss level
% Dtau1_gs = 0.768; % Gouveia and Strauss curvature
% Dtau2_gs = 1; % Gouveia and Strauss scale
% flat tax
Dtau_ss = 0;
Dtau0_gs = 0.3; % Gouveia and Strauss level
Dtau1_gs = 1e-2; % Gouveia and Strauss curvature
Dtau2_gs = 1; % Gouveia and Strauss scale
GY_ratio = 0.09;

% income process
% NOTE(wenlan): income process comes from Storeletten etc. (2004)
Drho_y = 0.9989; % persistence of wage shock
Dvar_alpha = 0.2105; % variance of skill
Dvar_eta = 0.0161;
Dvar_epsilon = 0.0630;
% Dvar_epsilon = 0.2;
Dbeta1 = 0.015; % Mincerian age slope
Dbeta2 = -0.0001; % Mincerian age quadratic

% insurance
% medicaid
Dgamma_MC = 0.15; % coninsurance
Dchi_MC = 0.017; % deductible
Dgamma_EI = 0.15;
Dchi_EI = 0.017;
Dgamma_MD = 0.15;
Dchi_MD = 0.017;
ninc_test = 0.01; % mean test for medicaid

% safenet
c_low = 1e-12; % consumption at bankruptcy

%% prices
r = 0.0584;
w = 1.0812;
Y = 3;
w_mean = 1.6940;
Dtau0_gs = 0.3440; % Gouveia and Strauss level
p_EI = 0.01*ones(1, Jw); % insurance premium
p_MD = 0;
% Tss = [0.2 0.6];
Tr = 0.0012;
Bequest = 0;

%% adjust for scale
DA = 1;
w = DA * w;
% Tss = DA * Tss;

%% simulation
seed = 0823;
num_of_agents = 1e4;

%% computation
num_of_threads = 8;
num_of_procs = 8;
nslots = getenv('NSLOTS');
if (strcmp(nslots, '') == 0)
    num_of_threads = str2double(nslots)
    num_of_procs = str2double(nslots);
end
% num_of_threads = 1;
tol_x = 1e-6; % first order condition to be satisfied at optimal
tol_con = 1e-12;
% eq
tol_eq = 1e-2;
% update speed
iter_max_eq = 100;
speed_eq_small = 0.1;
speed_eq_big = 0.5;
% boundary value
c_min = 1e-6;
l_min = 1e-6;
n_min = 1e-6;
n_max = 1;
m_min = 1e-6;

% initial values;
KLRatio = 4.7207;
r = Dalpha * KLRatio^(Dalpha-1) - Ddelta;
w = (1-Dalpha) * KLRatio^Dalpha;
w_mean = 1.6940;
Dtau0_gs = 0.3440;
Y = 1.9582;

% read insurance premium
Y2003 = 39682.47;
p_EI_data = csvread('MOMENTS/pEI_by_age.csv');
p_EI = p_EI_data(:,2)' / Y2003;

%% params initial
cali_x0 = [Dbeta Dd Dpsi2 Dhdelta0 Dhdelta1 Dhdelta2 Dmu_z0 Dmu_z1 Dmu_z2 Dvar_z DA_h Dtheta];
cali_lb = [0 0 1 0 0 0 -inf 0 0 0.1 0 0.1];
cali_ub = [inf inf inf inf inf inf inf inf inf inf inf inf];
cali_scale = [0.01 1 0.1 1e-3 1e-4 1e-5 0.1 1e-3 1e-5 5e-2 0.01 0.1];
assert(length(cali_x0) == length(cali_lb));
assert(length(cali_x0) == length(cali_ub));
assert(length(cali_x0) == length(cali_scale));

% EqObjects0 = [ ...
%     Dbeta ...
%     Dd ...
%     Dpsi2 ...
%     Dhdelta0 ...
%     Dhdelta1 ...
%     Dhdelta2 ...
%     Dmu_z0 ...
%     Dmu_z1 ...
%     Dmu_z2 ...
%     Dvar_z ...
%     DA_h ...
%     Dtheta ...
%     KLRatio ...
%     p_EI0 ...
%     p_EI1 ...
%     p_EI2 ...
%     Dtau0_gs ...
%     w_mean ...
%     ];

EqObjects0 = [
%     Dbeta ...
%     Dvar_z ...
    Y ...
    KLRatio ...
    Tr ...
    Dtau0_gs ...
    w_mean ...
%     Dd ...
%     Dpsi2 ...
%     Dhdelta0 ...
%     Dhdelta1 ...
%     Dhdelta2 ...
%     Dmu_z0 ...
%     Dmu_z1 ...
%     Dmu_z2 ...
%     DA_h ...
%     Dtheta ...
    ];

EqObjectsLb = [0 0.1 0 0 0 0 0 0 0 -inf -inf -inf -inf -inf -inf 1e-1 0.1];
EqObjectsLb = -inf * ones(size(EqObjects0));

EqObjectsScale = [1 1 1 1 1 1 1 1 1 1e-5 1e-5 1e-5 1 1 1e-5 1 1];
EqObjectsScale = ones(size(EqObjects0));

%% pack variables
params = v2struct;
end
