function [eq_result, exitflag] = EQ(params)
% Inferring Health Inequality
% Compute equilibrium
% Author: Wenlan Luo at Georgetown University, luowenlan at gmail.com

function budget = budget_vary_tau0(taxable_income, population, tau0, G, Taxss, Bequest)
    tax = GS(taxable_income, tau0, params.Dtau1_gs, params.Dtau2_gs);
    budget = wmean(tax(:), population(:)) - G + Taxss + Bequest;
end

speed = params.speed_eq_small;
iter_max = params.iter_max_eq;
tol = params.tol_eq;
last_delta_direction = [-1 -1 -1 -1 -1];
eq_metric = 1;
iter = 0;

while (eq_metric >= tol && iter < iter_max)
    iter = iter + 1;
    
    %
    OldEqObjects = [
        params.Y
        params.KLRatio
        params.Tr
        params.Dtau0_gs
        params.w_mean
        ]';
    
    params.r = params.Dalpha * params.KLRatio^(params.Dalpha-1) - params.Ddelta;
    params.w = (1-params.Dalpha) * params.KLRatio^params.Dalpha;
        
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
    % new_tau0_gs is computed to balance government budget
    TotalG = aggregate_result.G_MC + aggregate_result.G_MD + aggregate_result.Tss + params.GY_ratio*aggregate_result.Y;
    TaxGap = TotalG - aggregate_result.Bequest;
    new_tau0_gs = TaxGap / (aggregate_result.Tax / params.Dtau0_gs);
    
    % update objects
    NewEqObjects = [
        aggregate_result.Y
        aggregate_result.KLRatio
        aggregate_result.Tr
        new_tau0_gs
        aggregate_result.WMean
        ]';
    [OldEqObjects; NewEqObjects]
    eq_metric = max(abs(OldEqObjects - NewEqObjects));
    fprintf('Equilibrium metric: %.4f\n', eq_metric);
    
    % delta direction
    delta_direction = OldEqObjects ...
        >= NewEqObjects;
    if min(last_delta_direction == delta_direction) == 0
        % overshooting, be cautious
        speed = max(speed-0.1, params.speed_eq_small);
    else
        % safe, update aggresively
        speed = min(speed+0.1, params.speed_eq_big);
    end
    last_delta_direction = delta_direction;
    
    % compute data model moments distance
    if (isfield(params, 'MomentsData')==1)
        [aggregate_result.MomentsModel(:)';params.MomentsData(:)']
        moments_metric = sum((aggregate_result.MomentsModel - params.MomentsData).^2 .* params.MomentsWeight);
        fprintf('Moments metric: %.4f\n', moments_metric);
    else
        moments_metric = [];
    end
    
    
    params.KLRatio = params.KLRatio*(1-speed) + aggregate_result.KLRatio*speed;
    params.w_mean = params.w_mean*(1-speed) + aggregate_result.WMean*speed;
    params.Dtau0_gs = params.Dtau0_gs*(1-speed) + new_tau0_gs*speed;
    params.Y = params.Y*(1-speed) + aggregate_result.Y*speed;
    params.Tr = params.Tr*(1-speed) + aggregate_result.Tr*speed;
end

eq_result = v2struct(NewEqObjects, moments_metric);
exitflag = 0;

end

function s = wmean(x, weight)
s = sum(x.*weight)./sum(weight);
end
