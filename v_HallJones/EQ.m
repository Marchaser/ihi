function [eq_result, exitflag] = EQ(params)
% Inferring Health Inequality
% Compute equilibrium
% Author: Wenlan Luo at Georgetown University, luowenlan at gmail.com

IsCoverged = 0;

function budget = budget_vary_tau0(taxable_income, population, tau0, G, Taxss, Bequest)
    tax = GS(taxable_income, tau0, params.Dtau1_gs, params.Dtau2_gs);
    budget = wmean(tax(:), population(:)) - G + Taxss + Bequest;
end

speed = params.speed;

while(~IsCoverged)
    %
    OldEqObjects = [
        params.w ...
        params.r ...
        params.w_mean ...
        params.Dtau0_gs ... 
        params.Y ...
        params.p_EI ...
        ];
        
    % solve vfi
    tic;
    [vfi_result, vfi_exitflag] = VFI(params);
    toc;
    % simulate
    tic;
    [simulate_result, simulate_exitflag] = SIMULATE(vfi_result, params);
    toc;
    % aggregate
    [aggregate_result, aggregate_exitflag] = AGGREGATE(simulate_result, params);
    
     % tabulates some diagnostic stats
%     fprintf('%-8s%-8s%-8s%-8s\n', 'Object', 'Min', 'Max', 'Mean');
%     fprintf('%-8s%-8.4f%-8.4f%-8.4f\n', 'medical', min(simulate_result.m), max(simulate_result.m), mean(simulate_result.m));
%     [[1:params.J]' min(simulate_result.m)' max(simulate_result.m)']


    % compute implied equilibrium objects
    % new_tau2_gs is computed to balance government budget
    TotalG = aggregate_result.G_MC + aggregate_result.G_MD + aggregate_result.Tss + params.GY_ratio*aggregate_result.Y;
    
    budget_vary_tau0_inline = @(tau0) budget_vary_tau0(simulate_result.taxable_income, simulate_result.population, tau0, TotalG, aggregate_result.Taxss, aggregate_result.Bequest);
    fprintf('test tax:\n');
    budget_min = budget_vary_tau0_inline(0)
    budget_max = budget_vary_tau0_inline(1)
    if (budget_min>0)
        new_tau0_gs = 0;
    elseif (budget_max<0)
        new_tau0_gs = 1;
    else
        new_tau0_gs = fzero(budget_vary_tau0_inline, [0 1]);
    end
    
    % update objects
    NewEqObjects = [
        aggregate_result.new_w ...
        aggregate_result.new_r ...
        aggregate_result.WMean ...
        new_tau0_gs ...
        aggregate_result.Y ...
        aggregate_result.new_p_EI ...
        ];
    [OldEqObjects; NewEqObjects]
    eq_metric = max(abs([
        params.w - aggregate_result.new_w ...
        params.r - aggregate_result.new_r ...
        params.w_mean - aggregate_result.WMean ...
        params.Dtau0_gs - new_tau0_gs ...
        params.Y - aggregate_result.Y ...
        params.p_EI(2:end) - aggregate_result.new_p_EI(2:end) ...
        ]));
    fprintf('Equilibrium metric: %.4f\n', eq_metric);
    display(eq_metric);
    IsCoverged = (eq_metric < params.tol_eq);
    
    % compute data model moments distance
    if (isfield(params, 'MomentsData')==1)
        [aggregate_result.MomentsModel(:)';params.MomentsData(:)']
        moments_metric = sum((aggregate_result.MomentsModel - params.MomentsData).^2 .* params.MomentsWeight);
        fprintf('Moments metric: %.4f\n', moments_metric);
    else
        moments_metric = [];
    end
    
    params.w = params.w*(1-speed) + aggregate_result.new_w*speed;
    params.r = params.r*(1-speed) + aggregate_result.new_r*speed;
    params.w_mean = params.w_mean*(1-speed) + aggregate_result.WMean*speed;
    params.Dtau0_gs = params.Dtau0_gs*(1-speed) + new_tau0_gs*speed;
    params.Y = params.Y*(1-speed) + aggregate_result.Y*speed;
    params.p_EI = params.p_EI*(1-speed) + aggregate_result.new_p_EI*speed;
end

eq_result = v2struct(NewEqObjects, moments_metric);
exitflag = 0;

end

function s = wmean(x, weight)
s = sum(x.*weight)./sum(weight);
end
