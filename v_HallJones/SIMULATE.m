function [output, exitflag] = SIMULATE(vfi_result, params)
% Inferring Health Inequality
% Simulate the economy using value function solved
% Author: Wenlan Luo at Georgetown University, luowenlan at gmail.com

% unpack structure
v2struct(vfi_result);
v2struct(params);

%% initialize storage
clear Dist Stats;

k = k_min*ones(num_of_agents, J+1);
k(:,1) = Bequest;
m = zeros(num_of_agents, J);
ms = zeros(num_of_agents, J); % out_of_pocket 
hl = Dh1 * ones(num_of_agents, J+1);
hl_after_depreciation = zeros(num_of_agents, J);
c = zeros(num_of_agents, J);

% Decision made by young only
% wage
wage = zeros(num_of_agents, J);
% labor supply
e = zeros(num_of_agents, J); % labor efficiency
n = zeros(num_of_agents, J);
l = zeros(num_of_agents, J);
% insurance
insurance = zeros(num_of_agents, J);

% value
v = zeros(num_of_agents, J);

% ninc kinc
kinc = zeros(num_of_agents, J);
ninc = zeros(num_of_agents, J);
taxable_income = zeros(num_of_agents, J);

% tax
tax = zeros(num_of_agents, J);

% social security transfer
Tss = zeros(num_of_agents, J);

% death shock and alive status
survival = ones(num_of_agents, J);
population = ones(num_of_agents, J+1);

%% simulate young
% construct young's interpolation problem
pp_Ev_EI = tensor_pchip({k_grid,h_grid}, reshape(permute(Ev_EI, [6 1 2 3 4 5]), [], k_pts, h_pts));
pp_Ev_MD = tensor_pchip({k_grid,h_grid}, reshape(permute(Ev_MD, [6 1 2 3 4 5]), [], k_pts, h_pts));
pp_Ev_NI = tensor_pchip({k_grid,h_grid}, reshape(permute(Ev_NI, [6 1 2 3 4 5]), [], k_pts, h_pts));
pp_value = tensor_pchip({k_grid,h_grid}, reshape(permute(value_young, [6 1 2 3 4 5]), [], k_pts, h_pts));
% pp_Ev_EI = struct('form','MKLpp','breaks',{{k_grid,h_grid}},...
%         'Values',[],'coefs',[],'order',[4 4],...
%         'Method',[],'ExtrapolationOrder',[],'thread',num_of_threads,...
%         'orient','curvefit');
% pp_Ev_MD = pp_Ev_EI;
% pp_Ev_NI = pp_Ev_EI;
% pp_value = pp_Ev_EI;
% pp_Ev_EI.Values = reshape(permute(Ev_EI, [6 1 2 3 4 5]), [], k_pts, h_pts);
% pp_Ev_MD.Values = reshape(permute(Ev_MD, [6 1 2 3 4 5]), [], k_pts, h_pts);
% pp_Ev_NI.Values = reshape(permute(Ev_NI, [6 1 2 3 4 5]), [], k_pts, h_pts);
% pp_value.Values = reshape(permute(value_young, [6 1 2 3 4 5]), [], k_pts, h_pts);

% convert to mklpp
pp_Ev_EI = myppual(pp_Ev_EI);
pp_Ev_MD = myppual(pp_Ev_MD);
pp_Ev_NI = myppual(pp_Ev_NI);
pp_value = myppual(pp_value);
% attach num of threads
pp_Ev_EI.thread = num_of_threads;
pp_Ev_MD.thread = num_of_threads;
pp_Ev_NI.thread = num_of_threads;
pp_value.thread = num_of_threads;
value_vec_size = [Jw+1 epsilon_pts y_pts alpha_pts];

problem_young.options.Algorithm = 'donlp2';
% problem_young.options.Algorithm = 'npopt';
problem_young.options.NumThreads = num_of_threads;
problem_young.options.NumProcs = num_of_procs;
problem_young.options.TolX = tol_x;
problem_young.options.TolCon = tol_con;
problem_young.pp = pp_value;

% problem_young.hessian = zeros(50, num_of_agents);

init = [
    k_min*ones(1, num_of_agents);
    m_min*ones(1, num_of_agents);
    n_max*ones(1, num_of_agents);
    ];
x = init;
% x_below_chi = x;
% x_above_chi = x;

for j=1:1:Jw
    % wage shock is realized
    e(:,j) = exp(alpha_shock(:) + Dbeta1*j + Dbeta2*j^2 + y_shock(:,j) + epsilon_shock(:,j));
    wage(:,j) = w*e(:,j);
    
    % NOTE(wenlan): compute value interpolation vector function index
    % use current period value for insurance decision!!!!!!
    % memory order, epsilon, y, alpha
    Ev_idx = sub2ind(value_vec_size, j*ones(1, num_of_agents), epsilon_shock_idx(:,j)', y_shock_idx(:,j)', alpha_shock_idx(:)');
    % use future value for consumption saving!
    value_idx = sub2ind(value_vec_size, (j+1)*ones(1, num_of_agents), epsilon_shock_idx(:,j)', y_shock_idx(:,j)', alpha_shock_idx(:)');
    problem_young.idx = value_idx(:)';
    
    % make insurance decision
    % interpolate insurance value
    Ev_EI = myppual(pp_Ev_EI, [k(:,j)';hl(:,j)'], [], Ev_idx(:)', []);
    Ev_MD = myppual(pp_Ev_MD, [k(:,j)';hl(:,j)'], [], Ev_idx(:)', []);
    Ev_NI = myppual(pp_Ev_NI, [k(:,j)';hl(:,j)'], [], Ev_idx(:)', []);    
    % compare value and make insurance decision
    insurance(:,j) = [ ...
        1 * (Ev_EI>=Ev_MD & Ev_EI>=Ev_NI) + ...
        2 * (Ev_MD>Ev_EI & Ev_MD>Ev_NI) + ...
        3 * (Ev_NI>Ev_EI & Ev_NI>=Ev_MD)]';

    % disease shock is realized
    depreciation = exp(z_shock(:,j)*Dvar_z(j)^0.5 + mu_z(j));
    hl_after_depreciation(:,j) = max(hl(:,j)*(1-Dhdelta(j))-depreciation, 0);
%     hl_after_depreciation(:,j) = hl(:,j)*(1-Dhdelta)-depreciation;
    
    problem_young.shared_data = [
        Dbeta; Dlambda; Dzeta; Drho; Dnu; Dd;
        Dpsi1; Dpsi2;
        DA_h(j); Dtheta(j);
        Dtau_ss; Dtau0_gs; Dtau1_gs; Dtau2_gs;
        ];
    
    % kinc
    budget = (k(:,j)+Tr) * (1+r);
    kinc(:,j) = (k(:,j)+Tr) * r;
    
    % attach data and construct problem
    % some shared structure
    problem_young.lb = repmat([k_min m_min n_min]', 1, num_of_agents);
    problem_young.lb = [
        k_min*ones(1, num_of_agents);
%         max((max(h_min - hl_after_depreciation(:,j)',0)/DA_h(j)).^(1/Dtheta(j)), m_min);
        m_min*ones(1, num_of_agents);
        n_min*ones(1, num_of_agents);
        ];
    problem_young.ub = [
        k_max*ones(1, num_of_agents);
%         1e20*ones(1, num_of_agents);
%         max((depreciation(:)/DA_h(j)).^(1/Dtheta(j)),m_min)';  
        max(((hl(:,j)*(1-Dhdelta(j))-hl_after_depreciation(:,j))/DA_h(j)).^(1/Dtheta(j)),m_min)';  
        n_max*ones(1, num_of_agents);
        ];
    % first one is consumption, second is leisure
    problem_young.clb = [
        c_min*ones(1, num_of_agents);
        l_min*ones(1, num_of_agents);
        ];
    problem_young.cub = 1e20*ones(2, num_of_agents);
    data = [...
        budget(:)';
        hl(:,j)';
        hl_after_depreciation(:,j)';
        wage(:,j)';
        kinc(:,j)';
        ];
    % carefully attach insruance policy
    problem_young.data = [
        data;
        [(insurance(:,j)==1)*Dgamma_EI + (insurance(:,j)==2)*Dgamma_MD + (insurance(:,j)==3)*1]';
        [(insurance(:,j)==1)*Dchi_EI + (insurance(:,j)==2)*Dchi_MD + (insurance(:,j)==3)*0]';
        ];
    % now budget should take away insurance premium
    problem_young.data(1,:) = [budget(:) - ...
        (insurance(:,j)==1)*p_EI(j) - ...
        (insurance(:,j)==2)*p_MD]';
    % labor supply should be bounded to means test level for MD
    problem_young.ub(3,:) = [ ...
        (insurance(:,j)==2)*ninc_test ./ wage(:,j) + ...
        (insurance(:,j)~=2)*n_max ...
        ]';
    problem_young.x0 = init;
    problem_young.x0 = x;

    % a final check for is_feasible, since interpolated insurance value may
    % be not reliable
    max_ninc_under_test = min(wage(:,j)*n_max, ninc_test);
    budget_full_work = budget(:) + wage(:,j)*(1-Dtau_ss)*n_max - GS(wage(:,j)*(1-Dtau_ss)*n_max+kinc(:,j), Dtau0_gs, Dtau1_gs, Dtau2_gs);
    budget_under_test = budget(:) + max_ninc_under_test*(1-Dtau_ss) - GS(max_ninc_under_test*(1-Dtau_ss)+kinc(:,j), Dtau0_gs, Dtau1_gs, Dtau2_gs);
    is_feasible = ...
        (insurance(:,j)==1).*(budget_full_work - p_EI(j)) + ...
        (insurance(:,j)==2).*(budget_under_test - p_MD) + ...
        (insurance(:,j)==3).*(budget_full_work) >= ...
        k_min+c_min+m_min;
    % make sure it's feasible
    problem_young.data(1,is_feasible==0) = k_min+c_min+m_min;
    [x,v_j,opt_exitflag] = dpopt_utility_young(problem_young);
    % Sometimes feasibility search fails, in those cases, solve problem
    % with initial values always feasible
    fail_index = (x<problem_young.lb | x>problem_young.ub);
    fail_index = max(fail_index);
    if (max(fail_index(:))==1)
        resolve = problem_young;
        resolve.x0 = init(:,fail_index); % use the always feasible initial value
        resolve.lb = resolve.lb(:,fail_index);
        resolve.ub = resolve.ub(:,fail_index);
        resolve.clb = resolve.clb(:,fail_index);
        resolve.cub = resolve.cub(:,fail_index);
        resolve.data = resolve.data(:,fail_index);
        resolve.idx = resolve.idx(:,fail_index);
        [resolve_x, resolve_v_j] = dpopt_utility_young(resolve);
        x(:,fail_index) = resolve_x;
        v_j(:,fail_index) = resolve_v_j;
    end
    % if still fails, then enforce it to lower bound to force feasibility
    x = min(max(x, problem_young.lb), problem_young.ub);
    
    % penalize if it's not feasible
    v_j(is_feasible==0) = -1e20;
    v(:,j) = v_j;
    
    % law of motion
    k(:,j+1) = x(1,:)';
    m(:,j) = x(2,:)';
    
    n(:,j) = x(3,:)';
    hl(:,j+1) = hl_after_depreciation(:,j) + DA_h(j).*(m(:,j).^Dtheta(j));
    
    % tax
    ninc(:,j) = n(:,j).*wage(:,j)*(1-Dtau_ss);
    taxable_income(:,j) = kinc(:,j)+ninc(:,j);
    tax(:,j) = GS(kinc(:,j)+ninc(:,j), Dtau0_gs, Dtau1_gs, Dtau2_gs);
    
    % take away insurance reimbursement
    ms(:,j) = m(:,j) - ...
        (1- ([(insurance(:,j)==1)*Dgamma_EI + (insurance(:,j)==2)*Dgamma_MD + (insurance(:,j)==3)*1])) .* ...
        max(m(:,j) - ([(insurance(:,j)==1)*Dchi_EI + (insurance(:,j)==2)*Dchi_MD + (insurance(:,j)==3)*0]), 0);
    c(:,j) = problem_young.data(1,:)'+ninc(:,j)-tax(:,j)-k(:,j+1)-ms(:,j);
    l(:,j) = 1-n(:,j);
    
    % uncomment this to shut down labor leisure choice
    l(:,j) = 1;
    
    % death shock
%     survival(:,j) = Dpsi2-hl(:,j+1).^(-Dpsi1);
    survival(:,j) = 1 - exp(-Dpsi1*hl(:,j+1).^Dpsi2);
    population(:,j+1) = population(:,j) .* survival(:,j);
end

%% simulate old
% construct old's interpolation problem
pp_value = tensor_pchip({k_grid,h_grid}, reshape(permute(value_old, [6 1 2 3 4 5]), [], k_pts, h_pts));
% pp_value = struct('form','MKLpp','breaks',{{k_grid,h_grid}},...
%         'Values',[],'coefs',[],'order',[4 4],...
%         'Method',[],'ExtrapolationOrder',[],'thread',num_of_threads,...
%         'orient','curvefit');
% pp_value.Values = reshape(permute(value_old, [6 1 2 3 4 5]), [], k_pts, h_pts);
% convert to mklpp
pp_value = myppual(pp_value);
% attach num of threads
pp_value.thread = num_of_threads;
value_vec_size = [Jr+1 epsilon_pts y_pts alpha_pts];

problem_old.options.Algorithm = 'donlp2';
problem_old.options.NumThreads = num_of_threads;
problem_old.options.NumProcs = num_of_procs;
problem_old.options.TolX = tol_x;
problem_old.options.TolCon = tol_con;
problem_old.pp = pp_value;

% use young's solution as initial value
init = init(1:2,:);
% x = x(1:2,:);
x = init;
for j=(Jw+1):1:J
    % NOTE(wenlan): for olds, Tss is determined by last period wage
    % compute value interpolation vector function index
    % use future value for consumption saving!
    value_idx = sub2ind(value_vec_size, (j-Jw+1)*ones(1, num_of_agents), epsilon_shock_idx(:,Jw)', y_shock_idx(:,Jw)', alpha_shock_idx(:)');
    problem_old.idx = value_idx(:)';
    
    % disease shock is realized
    depreciation = exp(z_shock(:,j)*Dvar_z(j)^0.5 + mu_z(j));
    hl_after_depreciation(:,j) = max(hl(:,j)*(1-Dhdelta(j))-depreciation, 0);
%     hl_after_depreciation(:,j) = hl(:,j)*(1-Dhdelta)-depreciation;
    
    problem_old.shared_data = [
        Dbeta; Dlambda; Dzeta; Drho; Dnu; Dd;
        Dpsi1; Dpsi2;
        DA_h(j); Dtheta(j);
        Dtau_ss; Dtau0_gs; Dtau1_gs; Dtau2_gs;
        ];

    % Tss
    Tss(:,j) = TSS(wage(:,Jw), PhiCoefs(2), PhiCoefs(1), w_mean);
    
    % kinc
    kinc(:,j) = (k(:,j)+Tr) * r ;
%     budget = Tss(:,j) + (k(:,j)+Tr) * (1+r) - GS(kinc(:,j),Dtau0_gs,Dtau1_gs,Dtau2_gs);
    budget = Tss(:,j) + (k(:,j)+Tr) * (1+r) - r*(kinc(:,j)+Tr)*Dtau0_gs;
    
    % attach data and construct problem
    % some shared structure
    problem_old.lb = repmat([k_min m_min]', 1, num_of_agents);
    problem_young.lb = [
        k_min*ones(1, num_of_agents);
%         max((max(h_min - hl_after_depreciation(:,j)',0)/DA_h(j)).^(1/Dtheta(j)), m_min);
        m_min*ones(1, num_of_agents);
        ];
    problem_old.ub = [
        k_max*ones(1,num_of_agents);
%         1e20*ones(1, num_of_agents);
%         max((depreciation(:)/DA_h(j)).^(1/Dtheta(j)),m_min)'; 
        max(((hl(:,j)*(1-Dhdelta(j))-hl_after_depreciation(:,j))/DA_h(j)).^(1/Dtheta(j)),m_min)';
        ];
    % first inequality constraint is consumption constraint
    problem_old.clb = c_min*ones(1, num_of_agents);
    problem_old.cub = 1e20*ones(1, num_of_agents);
    data = [...
        budget(:)';
        hl(:,j)';
        hl_after_depreciation(:,j)';
        ];
    % carefully attach insruance policy
    problem_old.data = [
        data;
        Dgamma_MD*ones(1,num_of_agents);
        Dchi_MD*ones(1,num_of_agents);
        ];
    % now medical expenditure should be bounded above
    problem_old.x0 = init;
    problem_old.x0 = x;
    % always feasible
    [x,v_j,opt_exitflag] = dpopt_utility_old(problem_old);
    % Sometimes feasibility search fails, in those cases, solve problem
    % with initial values always feasible
    fail_index = (x<problem_old.lb | x>problem_old.ub);
    fail_index = max(fail_index);
    if (max(fail_index(:))==1)
        resolve = problem_old;
        resolve.x0 = init(:,fail_index); % use the always feasible initial value
        resolve.lb = resolve.lb(:,fail_index);
        resolve.ub = resolve.ub(:,fail_index);
        resolve.clb = resolve.clb(:,fail_index);
        resolve.cub = resolve.cub(:,fail_index);
        resolve.data = resolve.data(:,fail_index);
        resolve.idx = resolve.idx(:,fail_index);
        [resolve_x, resolve_v_j] = dpopt_utility_old(resolve);
        x(:,fail_index) = resolve_x;
        v_j(:,fail_index) = resolve_v_j;
    end
    x = min(max(x, problem_old.lb), problem_old.ub);
    
    % bankruptcy decision
    v(:,j) = v_j;
    
    % law of motion
    k(:,j+1) = x(1,:);
    m(:,j) = x(2,:);
    hl(:,j+1) = hl_after_depreciation(:,j) + DA_h(j).*(m(:,j).^Dtheta(j));
    
    % tax
    ninc(:,j) = Tss(:,j);
    taxable_income(:,j) = kinc(:,j);
    tax(:,j) = r*(kinc(:,j)+Tr)*Dtau0_gs;

    % everyone has Medicare insurance
    ms(:,j) = m(:,j) - (1-Dgamma_MC)*max(m(:,j)-Dchi_MC, 0);
    c(:,j) = budget(:)-tax(:,j)-k(:,j+1)-ms(:,j);
    l(:,j) = 1;
    
    % death shock
%     survival(:,j) = Dpsi2-hl(:,j+1).^(-Dpsi1);
    survival(:,j) = 1 - exp(-Dpsi1*hl(:,j+1).^Dpsi2);
    population(:,j+1) = population(:,j) .* survival(:,j);
end

k(:,J+1) = [];
h = hl(:,2:J+1);
hl(:,J+1) = [];
population(:,J+1) = [];

output = v2struct(k,m,ms,hl,hl_after_depreciation,h,c,wage,e,n,l,insurance,kinc,ninc,taxable_income,tax,Tss,population,survival,v);
exitflag = 0;
end
