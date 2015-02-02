% finding equilibrium using fsolve
function [eq_result, exitflag] = EQ3(params)
params.aggregate_display = false;
last_x = [];
eq_objects = [];
    function eq_metric = compute_eq_metric(x)      
        last_x = x;

        x = x .* params.EqObjectsScale;
        
        x
        
        x = max(x, params.EqObjectsLb);
        

        newparams = params;

        newparams.Dbeta = x(1);
        newparams.Dvar_z = x(2);
        newparams.Y = x(3);
        newparams.KLRatio = x(4);
        newparams.Tr = x(5);
        newparams.Dtau0_gs = x(6);
        newparams.w_mean = x(7);
        newparams.Dd = x(8);
        newparams.Dpsi2 = x(9);
        newparams.Dhdelta0 = x(10);
        newparams.Dhdelta1 = x(11);
        newparams.Dhdelta2 = x(12);
        newparams.Dmu_z0 = x(13);
        newparams.Dmu_z1 = x(14);
        newparams.Dmu_z2 = x(15);
        newparams.DA_h = x(16);
        newparams.Dtheta = x(17);
        
        newparams.r = params.Dalpha * newparams.KLRatio^(params.Dalpha-1) - params.Ddelta;
        newparams.w = (1-params.Dalpha) * newparams.KLRatio^params.Dalpha;
        
        newparams = GRID(newparams);
        newparams = RNG(newparams);
        
        tic;
        vfi_result = VFI(newparams);
        toc;
        tic;
        simulate_result = SIMULATE(vfi_result, newparams);
        toc;
        aggregate_result = AGGREGATE(simulate_result, newparams);
        
        eq_objects = [
            aggregate_result.KYRatio
            aggregate_result.MCv
            aggregate_result.Y
            aggregate_result.KLRatio
            aggregate_result.Tr
            aggregate_result.GYRatio
            aggregate_result.WMean
            aggregate_result.MMeanTop50(1)
            aggregate_result.MMeanTop50(end)
            aggregate_result.MMeanBot50(1)
            aggregate_result.MMeanBot50(end)
            aggregate_result.PsiTop50(1)
            aggregate_result.PsiTop50(20)
            aggregate_result.PsiTop50(end)
            aggregate_result.PsiBot50(1)
            aggregate_result.PsiBot50(20)
            aggregate_result.PsiBot50(end)
            ]';
        
        target_objects = [
            newparams.KYRatioData
            newparams.MCvData
            newparams.Y
            newparams.KLRatio
            newparams.Tr
            newparams.GYRatioData
            newparams.w_mean
            newparams.MedicalDataTop50(1)
            newparams.MedicalDataTop50(end)
            newparams.MedicalDataBot50(1)
            newparams.MedicalDataBot50(end)
            newparams.PsiDataTop50(1)
            newparams.PsiDataTop50(20)
            newparams.PsiDataTop50(end)
            newparams.PsiDataBot50(1)
            newparams.PsiDataBot50(20)
            newparams.PsiDataBot50(end)
            ]';
        
        eq_metric = eq_objects - target_objects;
        
        [eq_objects; target_objects]
    end

options = optimoptions('fsolve','Display','iter','DiffMinChange',1e-6,'TolFun',1e-2);
problem.objective = @compute_eq_metric;
problem.x0 = params.EqObjects0 ./ params.EqObjectsScale;
problem.solver = 'fsolve';
problem.options = options;
fsolve(problem);

% return
eq_result = v2struct(vfi_result, simulate_result, aggregate_result, eq_objects, last_x);
exitflag = 0;
end
