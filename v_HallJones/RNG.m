function newparams = RNG(params)
% Inferring Health Inequality
% Gen random number for alpha shock, epsilon shock, z shock, and death
% shock
% Author: Wenlan Luo at Georgetown University, luowenlan at gmail.com

%% unpack input
v2struct(params);
clear params;

% set seed
rng(seed);

% alpha shock
cum_alpha_trans = cumsum(alpha_trans(1,:), 2);
alpha_cut_point = floor(num_of_agents * cum_alpha_trans);
alpha_cut_point(end) = num_of_agents;  % avoid rounding error
alpha_cut_point = [1 alpha_cut_point];
alpha_shock_idx = zeros(1, num_of_agents);
for i=1:alpha_pts
    alpha_shock_idx(alpha_cut_point(i):alpha_cut_point(i+1)) = i;
end
rand_idx = randperm(num_of_agents);
alpha_shock_idx = alpha_shock_idx(rand_idx);
alpha_shock = alpha_grid(alpha_shock_idx);

% epsilon shock
cum_epsilon_trans = cumsum(epsilon_trans(1,:), 2);
epsilon_cut_point = floor(num_of_agents * cum_epsilon_trans);
epsilon_cut_point(end) = num_of_agents;  % avoid rounding error
epsilon_cut_point = [1 epsilon_cut_point];
epsilon_shock_idx_base = zeros(1, num_of_agents);
for i=1:epsilon_pts
    epsilon_shock_idx_base(epsilon_cut_point(i):epsilon_cut_point(i+1)) = i;
end
epsilon_shock_idx = zeros(num_of_agents, Jw);
epsilon_shock = zeros(num_of_agents, Jw);
for j=1:Jw
    rand_idx = randperm(num_of_agents);
    epsilon_shock_idx(:,j) = epsilon_shock_idx_base(rand_idx);
    epsilon_shock(:,j) = epsilon_grid(epsilon_shock_idx(:,j));
end

% y shock
[y_shock_idx, y_shock] = rng_markov(y_grid, y_trans, num_of_agents, Jw);

% population wage
e = zeros(num_of_agents, Jw);
for j=1:Jw
    e(:,j) = exp(alpha_shock(:) + Dbeta1*j + Dbeta2*j^2 + y_shock(:,j) + epsilon_shock(:,j));
end
e_mean = mean(e');
PhiCoefs = polyfit(e(:,j), e_mean(:), 1);

% z shock
z_shock = zeros(num_of_agents, J);
% align catastrophic shock and disease shock
z_shock_xtile = zeros(1, num_of_agents);
num_of_non_disease = round(num_of_agents*(1-Ddisease));
num_of_cata = round(num_of_agents*Ddisease*Dcata);
num_of_normal = num_of_agents - num_of_non_disease - num_of_cata;
z_shock_xtile(1:num_of_non_disease) = 0;
z_shock_xtile(num_of_non_disease+1:num_of_non_disease+num_of_cata) = 1;
% normal disease
z_shock_xtile(num_of_non_disease+num_of_cata+1:end) = linspace(1/num_of_agents, 1-1/num_of_agents, num_of_normal);
for j=1:J
    z_shock(:,j) = norminv(z_shock_xtile(randperm(num_of_agents)));
end
% z_shock = randn(num_of_agents, J);
% z_shock(rand(num_of_agents,J)<Dcata) = 10;

newparams = v2struct;
end