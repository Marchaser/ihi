% Inferring Health Inequality
% Some test statistics
% Author: Wenlan Luo at Georgetown University, luowenlan at gmail.com

% remove age effect
v0 = load('v0');
v1 = load('v1');


log_c1 = log(v1.c) - repmat(mean(log(v1.c)),v1.num_of_agents,1);
log_c0 = log(v0.c) - repmat(mean(log(v0.c)),v0.num_of_agents,1);