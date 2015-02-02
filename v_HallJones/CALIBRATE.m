function CALIBRATE(params)
eq_result = [];
WarmUpEqObjects = params.EqObjects0;
LastEqObjects = WarmUpEqObjects;
min_metric_found = 999;
    function metric = compute_moments_metric(cali)
        % scale back
        cali = cali .* params.cali_scale;
        % display parameters
        fprintf('Current parameters:\n');
        fprintf('Dbeta, Dd, Dpsi2: %g, %g, %g\n', cali(1), cali(2), cali(3));
        fprintf('Dhdelta0, Dhdelta1, Dhdelta2: %g, %g, %g\n', cali(4), cali(5), cali(6));
        fprintf('Dmu_z0, Dmu_z1, Dmu_z2: %g, %g, %g\n', cali(7), cali(8), cali(9));
        fprintf('Dvar_z, DA_h, Dtheta: %g, %g, %g\n', cali(10), cali(11), cali(12));
        newparams = SETUP;
        % overwrite moments to params
        newparams.Dbeta = cali(1);
        newparams.Dd = cali(2);
        newparams.Dpsi2 = cali(3);
        newparams.Dhdelta0 = cali(4);
        newparams.Dhdelta1 = cali(5);
        newparams.Dhdelta2 = cali(6);
        newparams.Dmu_z0 = cali(7);
        newparams.Dmu_z1 = cali(8);
        newparams.Dmu_z2 = cali(9);
        newparams.Dvar_z = cali(10);
        newparams.DA_h = cali(11);
        newparams.Dtheta = cali(12);
        
        % overwrite warm_up price
        newparams.Y = WarmUpEqObjects(1);
        newparams.KLRatio = WarmUpEqObjects(2);
        newparams.Tr = WarmUpEqObjects(3);
        newparams.Dtau0_gs = WarmUpEqObjects(4);
        newparams.w_mean = WarmUpEqObjects(5);
        
        newparams = MOMENTS(newparams);
        newparams = GRID(newparams);
        newparams = RNG(newparams);
        
        eq_timer = tic;
        [eq_result, eq_exitflag] = EQ4(newparams);
        fprintf('Time for EQ:\n');
        toc(eq_timer);
        
        metric = eq_result.moments_metric;
        fprintf('Current metric: %g\n', metric);
        
        % store warmup price
        if (metric < min_metric_found)
            LastEqObjects = eq_result.NewEqObjects;
            min_metric_found = metric;
        end
    end

    function [stop,options,optchanged] = outfun(x,optimValues,state)
        stop = false;
        options = [];
        optchanged = [];
        switch state
            case 'iter'
                WarmUpEqObjects = LastEqObjects;
            otherwise
        end
    end

options = psoptimset('MaxMeshSize', 1, 'InitialMeshSize', 1, 'Display', 'iter', 'PollingOrder', 'Consecutive','CompletePoll','On','OutputFcns',@outfun);
patternsearch(@compute_moments_metric, params.cali_x0./params.cali_scale, [], [], [], [], params.cali_lb./params.cali_scale, params.cali_ub./params.cali_scale, [], options);
end
