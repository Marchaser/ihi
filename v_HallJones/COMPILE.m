function COMPILE
% Inferring Health Inequality
% Compile dpopt
% Author: Wenlan Luo at Georgetown University, luowenlan at gmail.com

options.MaxDim = 10;
options.MaxData = 50;
options.MaxVec = 1;
% options.debug = 'on';
% options.openmp = 'off';
options.donlp2 = 'on';
% options.snopt = 'on';
make_dpopt('utility_old',options);
make_dpopt('utility_young', options);
end