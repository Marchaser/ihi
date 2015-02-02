function newparams = GRID(params)
% Inferring Health Inequality
% Construct grids and meshes used for VFI
% Author: Wenlan Luo at Georgetown University, luowenlan at gmail.com

%% unpack input
v2struct(params);
clear params;

%% Grids
% k_grid
k_grid = csvread('KGRID.csv');
k_grid = k_grid(1:12:end) * DA;
k_min = min(k_grid);
k_pts = length(k_grid);
k_max= max(k_grid);

% h_grid
Dh1 = (-log(1-0.999)/Dpsi1)^(1/Dpsi2);
h_grid = [0.001 0.1 0.3 0.5 0.6 0.7 0.8 0.85 0.9 0.925 0.95 0.975 0.99 0.995 0.999]';
h_grid = [0.001 0.1 0.3 0.5 0.6 0.7 0.8 0.85 0.9 0.95 0.99 0.999]';
% h_grid = (Dpsi2-h_grid).^(-1/Dpsi1);
h_grid = (-log(1-h_grid)/Dpsi1).^(1/Dpsi2);
% h_grid
h_pts = length(h_grid);
h_min = min(h_grid);
h_max = max(h_grid);
% h_pts = 11;
% h_min = 1e-3;
% h_max = 1;
% h_grid = linspace(h_min, h_max, h_pts);

% alpha_grid, skill
alpha_pts = 2;
alpha_trans = [0.5 0.5];
alpha_grid = [-Dvar_alpha^0.5, Dvar_alpha^0.5];

% y_grid, earnings process persistent part
y_pts = 5;
y_range = 2;
[y_trans, y_grid] = markovappr(Drho_y, Dvar_eta^0.5, y_range, y_pts);

% epsilon_grid, earnings process, transitory part
epsilon_pts = 2;
epsilon_trans = [0.5 0.5; 0.5 0.5];
epsilon_grid = [-Dvar_epsilon^0.5, Dvar_epsilon^0.5];

% NOTE(wenlan): careful, epsilon goes in memory first
epsilon_y_trans = kron(y_trans, epsilon_trans);

% mu_z, mean of health shock

% mu_z = -3*ones(1,J);

% z_grid, health shock
z_pts = 5;
z_range = 2;
[z_trans, z_grid] = markovappr(0, 1, z_range, z_pts);
% get the first row only since it's iid
z_trans = z_trans(1,:);
% attach non disease and a small catastrophic shock
% z_grid = [-inf inf z_grid];
% z_trans = [(1-Ddisease) (1-Ddisease)*Dcata Ddisease*(1-Dcata)*z_trans];
% z_pts = z_pts+2;
% attach non disease shock
% z_grid = [-inf z_grid];
% z_trans = [(1-Ddisease) Ddisease*z_trans];
% z_pts = z_pts+1;

% health shock
Dhdelta = Dhdelta0 * ones(1,J) + Dhdelta1*[1:J] + Dhdelta2*[1:J].^2; % health depreciation
mu_z = Dmu_z0 + Dmu_z1*[1:J] + Dmu_z2*([1:J].^2);
Dvar_z = Dvar_z*ones(1,J);

% health production
Dtheta = Dtheta * ones(1,J);
DA_h = DA_h * ones(1,J);

% convert to nominal terms
Dchi_MC = 0.017 * Y;
Dchi_EI = 0.017 * Y;
Dchi_MD = 0.017 * Y;
p_EI = p_EI * Y;

%% pack params
newparams = v2struct;
end
