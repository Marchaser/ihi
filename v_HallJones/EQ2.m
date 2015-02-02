% finding equilibrium using fsolve
function [eq_result, exitflag] = EQ2(params)
params.aggregate_display = false;
last_x = [];
eq_objects = [];
    function eq_metric = compute_eq_metric(x)      
        last_x = x;
        
        newparams = params;

        newparams.Dbeta = x(1);
        newparams.Dvar_z = x(2);
        newparams.Y = x(3);
        newparams.KLRatio = x(4);
        newparams.Tr = x(5);
        newparams.Dtau0_gs = x(6);
        newparams.w_mean = x(7);
        
        newparams.r = params.Dalpha * newparams.KLRatio^(params.Dalpha-1) - params.Ddelta;
        newparams.w = (1-params.Dalpha) * newparams.KLRatio^params.Dalpha;
        
        newparams = GRID(newparams);
        newparams = RNG(newparams);
        
%         tic;
        vfi_result = VFI(newparams);
%         toc;
%         tic;
        simulate_result = SIMULATE(vfi_result, newparams);
%         toc;
        aggregate_result = AGGREGATE(simulate_result, newparams);
        
        eq_metric(1) = aggregate_result.KYRatio - newparams.KYRatioData;
        eq_metric(2) = aggregate_result.MCv - newparams.MCvData;
        eq_metric(3) = aggregate_result.Y - newparams.Y;
        eq_metric(4) = aggregate_result.KLRatio - newparams.KLRatio;
        eq_metric(5) = aggregate_result.Tr - newparams.Tr;
        eq_metric(6) = aggregate_result.GYRatio - newparams.GYRatioData;
        eq_metric(7) = aggregate_result.WMean - newparams.w_mean;
        
        eq_objects = [
            aggregate_result.KYRatio ...
            aggregate_result.MCv ...
            aggregate_result.Y ...
            aggregate_result.KLRatio ...
            aggregate_result.Tr ...
            aggregate_result.GYRatio ...
            aggregate_result.WMean ...
            ];
        
%         eq_objects
    end

options = optimoptions('fsolve','Display','iter','DiffMinChange',1e-6,'TolFun',1e-2);
problem.objective = @compute_eq_metric;
problem.x0 = params.EqObjects0;
problem.solver = 'fsolve';
problem.options = options;
fsolve(problem);

% return
eq_result = v2struct(vfi_result, simulate_result, aggregate_result, eq_objects, last_x);
exitflag = 0;
end
