% Inferring Health Inequality
% Entry point
% Author: Wenlan Luo at Georgetown University

% params = SETUP;
% params = MOMENTS(params);
% params = GRID(params);
% params = RNG(params);
% % tic;
% % [vfi_result, vfi_exitflag] = VFI(params);
% % toc;
% % tic;
% % [simulate_result, simulate_exitflag] = SIMULATE(vfi_result, params);
% % toc;
% % AGGREGATE(simulate_result, params);
% EQ(params);

params = SETUP;
params = MOMENTS(params);
CALIBRATE(params);
% EQ3(params);
% EQ2(params);
% CALIBRATE(params);
