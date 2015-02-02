% w/ health
clear;
params = SETUP;
params = GRID(params);
params = RNG(params);
tic;
[vfi_result, vfi_exitflag] = VFI(params);
toc;
tic;
[simulate_result, simulate_exitflag] = SIMULATE(vfi_result, params);
toc;
[aggregate_result, aggregate_exitflag] = AGGREGATE(simulate_result, params);
c0 = aggregate_result.simulate_result.c;

% w/o health
params.mu_z = -inf*ones(1,params.J);
params = GRID(params);
params = RNG(params);
tic;
[vfi_result, vfi_exitflag] = VFI(params);
toc;
tic;
[simulate_result, simulate_exitflag] = SIMULATE(vfi_result, params);
toc;
[aggregate_result, aggregate_exitflag] = AGGREGATE(simulate_result, params);
c1 = aggregate_result.simulate_result.c;

var_log_c0 = var(log(c0));
var_log_c1 = var(log(c1));

figure;
hold on;
plot([21:65],var_log_c1(1:45),'k-','LineWidth',2); % W/O health
plot([21:65],var_log_c0(1:45),'k--','LineWidth',2); % with health
xlabel('Age');
ylabel('Variance of log consumption');
legend('W/O Health','W/ Health','Location','southeast');
set(findall(gcf,'type','text'),'FontSize',16,'FontName','Times New Roman')
print('graph/VarLogCByAge.pdf', '-dpdf');
hold off;
