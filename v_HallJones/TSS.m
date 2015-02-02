% Inferring Health Inequality
% Social security as a function of last period wage
% Author: Wenlan Luo at Georgetown University
function tss = Tss(W, k0, k1, W_bar)
% Y is the last period wage
% predict average wage from last_period_wage
average_wage =  k0 + k1*W;
% normalize average_wage by population_average_wage
W_tilde = W/W_bar;
% compute social security transfer using brackets
if W_tilde<=0.3
    tss = 0.9*W_tilde;
elseif W_tilde<=2
    tss = 0.27 + 0.32*(W_tilde-0.3);
elseif W_tilde<=4.1
    tss = 0.81 + 0.15*(W_tilde-2);
else
    tss = 1.13;
end

% convert back to nominal value
tss = W_bar * tss;
end