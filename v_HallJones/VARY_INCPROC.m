clear;

%{
var_alpha = [0.0646 0.0785 0.0881 0.1052 0.1181 0.1277];
var_epsilon = [0.0777	0.0872	0.0968	0.1052	0.1137	0.1179];

var_log_c0 = zeros(size(var_alpha));
var_log_c1 = zeros(size(var_alpha));

params = SETUP;
% overwrite parameters
for i=1:length(var_alpha)
    params.Dvar_alpha = var_alpha(i);
    params.Dvar_epsilon = var_epsilon(i);
    params = GRID(params);
    params = RNG(params);
    tic;
    [vfi_result, vfi_exitflag] = VFI(params);
    toc;
    tic;
    [simulate_result, simulate_exitflag] = SIMULATE(vfi_result, params);
    toc;
    [aggregate_result, aggregate_exitflag] = AGGREGATE(simulate_result, params);
    
    % demean age effect
    log_c = log(aggregate_result.simulate_result.c(:,1:params.J)) - repmat(mean(log(aggregate_result.simulate_result.c(:,1:params.J))),params.num_of_agents,1);
    var_log_c0(i) = var(log_c(:));
%     save('var_log_c0.mat', 'var_log_c0');
end

% no health
params.mu_z = -inf*ones(1,params.J);
for i=1:length(var_alpha)
    params.Dvar_alpha = var_alpha(i);
    params.Dvar_epsilon = var_epsilon(i);
    params = GRID(params);
    params = RNG(params);
    tic;
    [vfi_result, vfi_exitflag] = VFI(params);
    toc;
    tic;
    [simulate_result, simulate_exitflag] = SIMULATE(vfi_result, params);
    toc;
    [aggregate_result, aggregate_exitflag] = AGGREGATE(simulate_result, params);
    
    % demean age effect
    log_c = log(aggregate_result.simulate_result.c(:,1:params.J)) - repmat(mean(log(aggregate_result.simulate_result.c(:,1:params.J))),params.num_of_agents,1);
    var_log_c1(i) = var(log_c(:));
%     save('var_log_c1.mat', 'var_log_c1');
end
%}

datain = load('var_log_c0');
var_log_c0 = datain.var_log_c0;
datain = load('var_log_c1');
var_log_c1 = datain.var_log_c1;

% load moments
var_log_cdata = csvread('MOMENTS/cons_ineq.csv');
var_log_cdata = var_log_cdata(:,2);
delta_var_log_cdata = var_log_cdata - var_log_cdata(1); 

delta_var_log_c0 = var_log_c0 - var_log_c0(1);
delta_var_log_c1 = var_log_c1 - var_log_c1(1);

var_log_c_figure = figure;
hold on;
year = [1980 1985 1990 1995 2000 2003];
plot(year,delta_var_log_c0,'k-','LineWidth',2);
plot(year,delta_var_log_c1,'k--','LineWidth',2);
plot([1980 2003],delta_var_log_cdata([1 end]),'k-.','LineWidth',2);
xlabel('Year');
ylabel('Variance of log consumption, change from 1980 value');
legend('W/O Health','W/ Health','Data','Location','southeast');
set(findall(gcf,'type','text'),'FontSize',16,'FontName','Times New Roman')
print(var_log_c_figure, 'graph/VarLogCByYear.pdf', '-dpdf');
hold off;
