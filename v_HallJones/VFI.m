function [result, exitflag] = VFI(params)
% Inferring Health Inequality
% Value function iterations solving agent's problem
% Author: Wenlan Luo at Georgetown University, luowenlan at gmail.com

%% unpack params
v2struct(params);
clear params;

%% old's problem
% construct old's mesh
[z_mesh, epsilon_mesh, y_mesh, alpha_mesh, k_mesh, h_mesh] = ...
    ndgrid(z_grid, epsilon_grid, y_grid, alpha_grid, k_grid, h_grid);
num_of_problems = numel(z_mesh);
% set the solve options
problem_old.options.NumThreads = num_of_threads;
problem_old.options.NumProcs = num_of_procs;
problem_old.options.Algorithm = 'donlp2';
problem_old.options.TolX = tol_x;
problem_old.options.TolCon = tol_con;

% vector function index, dimension z goes in memory first
vec_idx = repmat(1:epsilon_pts*y_pts*alpha_pts,z_pts,1);
vec_idx = vec_idx(:);
vec_idx = repmat(vec_idx, 1, k_pts*h_pts);
problem_old.idx = vec_idx(:)';
% last period wage
wage = w*exp(epsilon_mesh + y_mesh + alpha_mesh + Dbeta1*Jw + Dbeta2*Jw^2);
% Tss as a function of last period wage
Tss = TSS(wage, PhiCoefs(2), PhiCoefs(1), w_mean);
% upper bound
% budget = (1+r)*(k_mesh(:)+Tr) - GS(r*(k_mesh(:)+Tr),Dtau0_gs,Dtau1_gs,Dtau2_gs) + Tss(:);
budget = (1+r)*(k_mesh(:)+Tr) - r*(k_mesh(:)+Tr)*Dtau0_gs + Tss(:);

% initial values
init = repmat([k_min;m_min],1,num_of_problems);
x = init;

% initiate storage space
value_old = zeros([epsilon_pts y_pts alpha_pts k_pts h_pts Jr+1]);
policy_old = zeros([2 z_pts epsilon_pts y_pts alpha_pts k_pts h_pts Jr]);

% do iterations
for j=J:-1:Jw+1
    % construct spline
    v0 = reshape(value_old(:,:,:,:,:,j-Jw+1), [], k_pts, h_pts);
    ppform = struct('form','MKLpp','breaks',{{k_grid,h_grid}},...
        'Values',v0,'coefs',[],'order',[4 4],...
        'Method',[],'ExtrapolationOrder',[],'thread',num_of_threads,...
        'orient','curvefit');
    ppform = tensor_pchip({k_grid,h_grid}, v0);
    problem_old.pp = myppual(ppform);
    
    depreciation = exp(z_mesh*Dvar_z(j)^0.5 + mu_z(j));
    hl_after_depreciation = max(h_mesh*(1-Dhdelta(j))-depreciation, 0);
%     hl_after_depreciation = h_mesh*(1-Dhdelta)-depreciation;
    
    problem_old.shared_data = [
        Dbeta; Dlambda; Dzeta; Drho; Dnu; Dd;
        Dpsi1; Dpsi2;
        DA_h(j); Dtheta(j);
        Dtau_ss; Dtau0_gs; Dtau1_gs; Dtau2_gs;
        ];
    
    data = [
        budget(:)';
        h_mesh(:)';
        hl_after_depreciation(:)';
        ];
    problem_old.x0 = x;
%     problem_old.lb = repmat([k_min;m_min],1,num_of_problems);
    problem_old.lb = [
        k_min*ones(1,num_of_problems);
%         max((max(h_min - hl_after_depreciation(:)',0)/DA_h(j)).^(1/Dtheta(j)), m_min);
        m_min*ones(1,num_of_problems);
        ];
    problem_old.ub = [
        k_max*ones(1,num_of_problems);
        max(((h_mesh(:)*(1-Dhdelta(j))-hl_after_depreciation(:))/DA_h(j)).^(1/Dtheta(j)),m_min)';
%         1e20*ones(1,num_of_problems);
        ];
    problem_old.clb = c_min*ones(1,num_of_problems);
    problem_old.cub = 1e20*ones(1,num_of_problems);
    problem_old.data = [data; Dgamma_MC*ones(1,num_of_problems); Dchi_MC*ones(1,num_of_problems)];
    [x,v,opt_exitflag] = dpopt_utility_old(problem_old);
    policy_old(:,:,:,:,:,:,:,j-Jw) = reshape(x, [2 z_pts epsilon_pts y_pts alpha_pts k_pts h_pts]);
    
    % integrate over disease shock
    value_old(:,:,:,:,:,j-Jw) = reshape(z_trans * reshape(v, z_pts, []), [epsilon_pts y_pts alpha_pts k_pts h_pts]);
end


%% young's problem
% construct young's mesh
[z_mesh, epsilon_mesh, y_mesh, alpha_mesh, k_mesh, h_mesh] = ...
    ndgrid(z_grid, epsilon_grid, y_grid, alpha_grid, k_grid, h_grid);
num_of_problems = numel(z_mesh);
% set the solve options
problem_young.options.NumThreads = num_of_threads;
problem_young.options.NumProcs = num_of_procs;
problem_young.options.Algorithm = 'donlp2';
problem_young.options.TolX = tol_x;
problem_young.options.TolCon = tol_con;

% shared vector function index, dimension z goes in memory first
vec_idx = repmat(1:epsilon_pts*y_pts*alpha_pts,z_pts,1);
vec_idx = vec_idx(:);
vec_idx = repmat(vec_idx, 1, k_pts*h_pts);
problem_young.idx = vec_idx(:)';
% upper bound
budget = (1+r)*(k_mesh(:) + Tr);

% shared data
kinc = r*(k_mesh(:) + Tr);

% initial values
init = [
    reshape(policy_old(1,:,:,:,:,:,:,1), 1, []);
    reshape(policy_old(2,:,:,:,:,:,:,1), 1, []);
    n_max*ones(1,num_of_problems);
    ];
init = [
    k_min*ones(1,num_of_problems);
    m_min*ones(1,num_of_problems);
    n_max*ones(1,num_of_problems);
    ];
x_EI = init;
x_MD = init;
x_NI = init;

% initiate storage space
value_young = zeros([epsilon_pts y_pts alpha_pts k_pts h_pts Jw+1]);
Ev_EI = value_young;
Ev_MD = value_young;
Ev_NI = value_young;
policy_young = zeros([3 z_pts epsilon_pts y_pts alpha_pts k_pts h_pts Jw]);
insurance_policy_young = zeros([epsilon_pts y_pts alpha_pts k_pts h_pts Jw]);
% last period value function is given by value_old
value_young(:,:,:,:,:,Jw+1) = value_old(:,:,:,:,:,1);

% do iterations
for j=Jw:-1:1
    % construct spline
    v0 = reshape(value_young(:,:,:,:,:,j+1), [], k_pts, h_pts);
    ppform = struct('form','MKLpp','breaks',{{k_grid,h_grid}},...
        'Values',v0,'coefs',[],'order',[4 4],...
        'Method',[],'ExtrapolationOrder',[],'thread',num_of_threads,...
        'orient','curvefit');
    ppform = tensor_pchip({k_grid,h_grid}, v0);
    problem_young.pp = myppual(ppform);
    
    depreciation = exp(z_mesh*Dvar_z(j)^0.5 + mu_z(j));
    hl_after_depreciation = max(h_mesh*(1-Dhdelta(j))-depreciation, 0);
%     hl_after_depreciation = h_mesh*(1-Dhdelta)-depreciation;
    
    wage = w*exp(alpha_mesh + Dbeta1*j + Dbeta2*j^2 + epsilon_mesh);
    budget_full_work = budget(:) + wage(:)*(1-Dtau_ss)*n_max - GS(wage(:)*(1-Dtau_ss)*n_max+kinc(:), Dtau0_gs, Dtau1_gs, Dtau2_gs);
    % shared_data
    problem_young.shared_data = [
        Dbeta; Dlambda; Dzeta; Drho; Dnu; Dd;
        Dpsi1; Dpsi2;
        DA_h(j); Dtheta(j);
        Dtau_ss; Dtau0_gs; Dtau1_gs; Dtau2_gs;
        ];
    
	% individual data
    data = [
        budget(:)';
        h_mesh(:)';
        hl_after_depreciation(:)';
        wage(:)';
        kinc(:)';
        ];
    problem_young.lb = repmat([k_min;m_min;n_min],1,num_of_problems);
    problem_young .lb = [
        k_min*ones(1,num_of_problems);
%         max((max(h_min - hl_after_depreciation(:)',0)/DA_h(j)).^(1/Dtheta(j)), m_min);
        m_min*ones(1,num_of_problems);
        n_min*ones(1,num_of_problems);
        ];
    problem_young.ub = [
        k_max*ones(1, num_of_problems);
        n_max*ones(1, num_of_problems);
        max(((h_mesh(:)*(1-Dhdelta(j))-hl_after_depreciation(:))/DA_h(j)).^(1/Dtheta(j)),m_min)';
%         1e20*ones(1, num_of_problems);
        ];
    
    % first one is consumption, second one is leisure
    problem_young.clb = [
        c_min*ones(1,num_of_problems);
        l_min*ones(1,num_of_problems);
        ];
    problem_young.cub = 1e20*ones(2,num_of_problems);
    
    % solve problem EI
    problem_young.x0 = x_EI;
    problem_young.ub(3,:) = n_max; % n is not bounded
    problem_young.data = [data; Dgamma_EI*ones(1,num_of_problems); Dchi_EI*ones(1,num_of_problems)];
    problem_young.data(1,:) = budget(:)' - p_EI(j); % pay premium
    % make sure it's feasible
    is_feasible = (budget_full_work(:)' - p_EI(j) >= k_min+c_min+m_min);
    problem_young.data(1,is_feasible==0) = k_min + c_min + m_min;
    
    % tic;
    [x_EI,v_EI,opt_exitflag] = dpopt_utility_young(problem_young);
    % toc;
    % adjust v to -inf if it's not feasible
    v_EI(is_feasible==0) = -1e20; 
    
    %{
    % max budget under means test income
    max_ninc_under_test = min(wage(:)'*n_max, ninc_test);
    budget_under_test = budget(:)' + max_ninc_under_test*(1-Dtau_ss) - GS(max_ninc_under_test*(1-Dtau_ss)+kinc(:)', Dtau0_gs, Dtau1_gs, Dtau2_gs);
    % solve problem MD
    problem_young.x0 = x_MD;
    problem_young.ub(3,:) = ninc_test ./ wage(:)'; % means test of labor income
    problem_young.data = [data; Dgamma_MD*ones(1,num_of_problems); Dchi_MD*ones(1,num_of_problems)];
    problem_young.data(1,:) = budget(:)' - p_MD;
    % make sure it's feasible
    is_feasible = (budget_under_test(:)' - p_MD >= k_min+c_min+m_min);    
    problem_young.data(1,is_feasible==0) = k_min + c_min + m_min;
    
    % tic;
    [x_MD,v_MD,opt_exitflag] = dpopt_utility_young(problem_young);
    % toc;
    % adjust v to -inf if it's not feasible
    v_MD(is_feasible==0) = -1e20;
    %}
    x_MD = x_EI;
    v_MD = v_EI;
    
    % solve problem NI
    problem_young.x0 = x_NI;
    problem_young.ub(3,:) = n_max;
    problem_young.data = [data; ones(1,num_of_problems); zeros(1,num_of_problems)];
    problem_young.data(1,:) = budget(:)';
    % always feasible
    % tic;
    [x_NI,v_NI,opt_exitflag] = dpopt_utility_young(problem_young);
    % toc;
    
    % integrate over disease shock
    Ev_EI(:,:,:,:,:,j) = reshape(z_trans * reshape(v_EI, z_pts, []), [epsilon_pts y_pts alpha_pts k_pts h_pts]);
    Ev_MD(:,:,:,:,:,j) = reshape(z_trans * reshape(v_MD, z_pts, []), [epsilon_pts y_pts alpha_pts k_pts h_pts]);
    Ev_NI(:,:,:,:,:,j) = reshape(z_trans * reshape(v_NI, z_pts, []), [epsilon_pts y_pts alpha_pts k_pts h_pts]);
    
    v1 = max(max(Ev_EI(:,:,:,:,:,j),Ev_MD(:,:,:,:,:,j)),Ev_NI(:,:,:,:,:,j));
    
    % integrate over income shock
    v2 = epsilon_y_trans * reshape(v1, epsilon_pts*y_pts, []);
    value_young(:,:,:,:,:,j) = reshape(v2, [epsilon_pts y_pts alpha_pts k_pts h_pts]);
end

result.value_old = value_old;
result.value_young = value_young;
result.Ev_EI = Ev_EI;
result.Ev_MD = Ev_MD;
result.Ev_NI = Ev_NI;
exitflag = 0;
end

function tempplot
mesh(squeeze(k_mesh(1,1,1,:,:)), squeeze(h_mesh(1,1,1,:,:)), squeeze(value_old(1,1,:,:,j+1)));
mesh(squeeze(k_mesh(1,1,1,1,:,:)), squeeze(h_mesh(1,1,1,1,:,:)), squeeze(value_old(1,1,1,:,:,j-Jw)));
mesh(squeeze(k_mesh(1,1,1,:,:)), squeeze(h_mesh(1,1,1,:,:)), squeeze(value_old(1,1,:,:,Jr-1)));
mesh(squeeze(k_mesh(1,1,1,:,:)), squeeze(h_mesh(1,1,1,:,:)), squeeze(policy_old(1,1,1,1,:,:,1)));
mesh(squeeze(k_mesh(1,1,1,:,:)), squeeze(h_mesh(1,1,1,:,:)), squeeze(policy_old(2,1,1,1,:,:,1)));
mesh(squeeze(k_mesh(1,1,1,:,:)), squeeze(h_mesh(1,1,1,:,:)), squeeze(policy_old(3,1,1,1,:,:,1)));
mesh(squeeze(k_mesh(1,1,1,:,:)), squeeze(h_mesh(1,1,1,:,:)), squeeze(insurance_policy_young(7,2,:,:,j)));
mesh(squeeze(k_mesh(1,1,1,:,:)), squeeze(h_mesh(1,1,1,:,:)), squeeze(value_young(1,1,:,:,Jw+1)));
mesh(squeeze(k_mesh(1,1,1,:,:)), squeeze(h_mesh(1,1,1,:,:)), squeeze(value_young(1,1,:,:,Jw)));
mesh(squeeze(k_mesh(1,1,1,:,:)), squeeze(h_mesh(1,1,1,:,:)), squeeze(value_young(1,1,:,:,j+1)));
mesh(squeeze(k_mesh(1,1,1,1,:,:)), squeeze(h_mesh(1,1,1,1,:,:)), squeeze(value_young(2,3,1,:,:,j+1)));

v_low_young_temp = reshape(v_low_young, [z_pts epsilon_pts alpha_pts k_pts h_pts]);
mesh(squeeze(k_mesh(1,1,1,:,:)), squeeze(h_mesh(1,1,1,:,:)), squeeze(v_low_young_temp(1,1,1,:,:)));

v_EI_temp = max(reshape(v_NI, [z_pts epsilon_pts alpha_pts k_pts h_pts]), 0);
mesh(squeeze(k_mesh(1,1,1,:,:)), squeeze(h_mesh(1,1,1,:,:)), squeeze(v_EI_temp(2,1,1,:,:)));
p_temp = max(reshape(x_EI, [3 z_pts epsilon_pts y_pts alpha_pts k_pts h_pts]), 0);
mesh(squeeze(k_mesh(1,1,1,1,:,:)), squeeze(h_mesh(1,1,1,1,:,:)), squeeze(p_temp(3,2,2,5,2,:,:)));
end
