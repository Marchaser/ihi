function CALIBRATE2(params)
eq_result = [];
WarmUpEqObjects = params.EqObjects0;
LastEqObjects = WarmUpEqObjects;
min_metric_found = 999;

    function [metric] = compute_moments_metric(cali)
        % scale back
        cali = cali .* params.cali_scale;
        % display parameters
        
        newparams = params;
        % overwrite moments to params
        newparams.Dd = cali(1);
        newparams.Dpsi2 = cali(2);
        newparams.Dhdelta0 = cali(3);
        newparams.Dhdelta1 = cali(4);
        newparams.Dhdelta2 = cali(5);
        newparams.Dmu_z0 = cali(6);
        newparams.Dmu_z1 = cali(7);
        newparams.Dmu_z2 = cali(8);
        newparams.DA_h = cali(9);
        newparams.Dtheta = cali(10);
        
        % overwrite warm_up price
        newparams.EqObjects0 = WarmUpEqObjects;
        
        fprintf('Current parameters:\n');
        cali
        WarmUpEqObjects
        
        eq_timer = tic;
        [eq_result, eq_exitflag] = EQ2(newparams);
        fprintf('Time for EQ:\n');
        toc(eq_timer);
        
        metric = sum((eq_result.aggregate_result.MomentsModel - params.MomentsData).^2 .* params.MomentsWeight);
%         metric = eq_result.aggregate_result.MomentsModel - params.MomentsData;
        
        % display something
        [metric]
        [eq_result.last_x]
        [eq_result.eq_objects]
        [eq_result.aggregate_result.MomentsModel; params.MomentsData]

        if (metric < min_metric_found)
            LastEqObjects = eq_result.last_x;
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

rng(0823);
options = psoptimset('MaxMeshSize', 1, 'InitialMeshSize', 1, 'Display', 'iter', 'PollingOrder', 'Consecutive', 'ScaleMesh', 'Off','CompletePoll','On','OutputFcns',@outfun);
patternsearch(@compute_moments_metric, params.cali_x0./params.cali_scale, [], [], [], [], params.cali_lb./params.cali_scale, params.cali_ub./params.cali_scale, [], options);

% options = optimoptions('fmincon','GradObj','off','Display','iter','DiffMinChange',1e-1,'Algorithm','sqp','OutputFcn',@outfun);
% problem.objective=@compute_moments_metric;
% problem.x0 = params.cali_x0 ./ params.cali_scale;
% problem.lb = params.cali_lb ./ params.cali_scale;
% problem.ub = params.cali_ub ./ params.cali_scale;
% problem.solver = 'fmincon';
% problem.options = options;
%
% fmincon(problem);
end

