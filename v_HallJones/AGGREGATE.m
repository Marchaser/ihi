function [aggregate_result, exitflag] = AGGREGATE(simulate_result, params)
v2struct(params);
v2struct(simulate_result);

% survival probability
Population = sum(population);
Psi = Population(:,2:J) ./ Population(:,1:J-1);
GrossSurvival = Population / num_of_agents;

% medical expenditure
MMean = wmean(ms(:), population(:));
MVar = wvar(ms(:), population(:));
MCv = MVar^0.5 / MMean;
% population_cross_sectional = population(:);
% ms_cross_sectional = ms(:);
% positive_medical_index = ms_cross_sectional>1e-5;
% MVar = wvar(log(ms_cross_sectional(), population);

% income
WMean = wmean(reshape(wage(:,1:Jw),[],1), reshape(population(:,1:Jw),[],1));
WVar = wvar(log(wage), population);
% IncVar = nanvar(log(kinc+ninc));

% saving
KMean = wmean(k, population);

% consumption
CMean = wmean(c, population);
CVar = wvar(log(c), population);

% labor
NMean = wmean(n, population);
NVar = wvar(log(n), population);

% per capita aggregate
K = wmean(k(:), population(:));
N = wmean(reshape(NMean(1:Jw),[],1), reshape(Population(1:Jw),[],1));
L = wmean(n(:).*e(:), population(:));
KLRatio = K/L;
Tax = wmean(tax(:), population(:));
Tss = wmean(Tss(:), population(:));
Taxss = L*w*Dtau_ss;
Bequest = sum(reshape((1-survival(:,1:J-1)).*k(:,2:J).*population(:,1:J-1),1,[])) / sum(population(:));

% implied price
Y = DA*K^Dalpha*L^(1-Dalpha);
new_w = DA*(1-Dalpha)*Y/L;
new_r = Dalpha*Y/K - Ddelta;

% government spending
% medicaid
G_MC = sum(sum(max(m(:,Jw+1:J)-Dchi_MC,0).*population(:,Jw+1:J))) / sum(Population) * (1-Dgamma_MC);
G_MD = sum(sum(max(m(:,1:Jw)-Dchi_MD,0).*(insurance(:,1:Jw)==2).*population(:,1:Jw))) / sum(Population) * (1-Dgamma_MD);

% private insurance
% new_p_EI = sum(max(m(:,1:Jw)-Dchi_EI,0).*(insurance(:,1:Jw)==1).*population(:,1:Jw)) ...
%     ./ sum((insurance(:,1:Jw)==1).*population(:,1:Jw)) ...
%     * (1-Dgamma_EI);
% new_p_EI(isnan(new_p_EI)) = 0;
% new_p_EI_Coefs = polyfit([2:Jw],new_p_EI(2:Jw),2);
% new_p_EI0 = new_p_EI_Coefs(3);
% new_p_EI1 = new_p_EI_Coefs(2);
% new_p_EI2 = new_p_EI_Coefs(1);
% profit of private insurance
Tr = sum(sum( ...
    (repmat(p_EI,num_of_agents,1) - (1-Dgamma_EI)*max(m(:,1:Jw)-Dchi_EI,0)) ...
    .*(insurance(:,1:Jw)==1).*population(:,1:Jw))) / sum(population(:));

% by income statistics
% medical expenditure
inc = ninc + kinc;
median_inc = zeros(1,J);
MMeanTop50 = zeros(1,J);
MMeanBot50 = zeros(1,J);
PsiTop50 = zeros(1,J);
PsiBot50 = zeros(1,J);
for j=1:J-1
    median_inc(j) = weightedMedian(inc(:,j),population(:,j)/Population(j));
    top50 = inc(:,j) >= median_inc(j);
    bot50 = inc(:,j) <= median_inc(j);
    MMeanTop50(j) = wmean(ms(top50==1,j),population(top50==1,j));
    MMeanBot50(j) = wmean(ms(bot50==1,j),population(bot50==1,j));
%     IncMeanTop50(j) = wmean(inc(top50==1,j), population(top50==1,j));
%     IncMeanBot50(j) = wmean(inc(bot50==1,j), population(bot50==1,j));
%     MMeanTop50(j) = wmean(ms(alpha_shock_idx(:)==2,j),population(alpha_shock_idx(:)==2,j));
%     MMeanBot50(j) = wmean(ms(alpha_shock_idx(:)==1,j),population(alpha_shock_idx(:)==1,j));
    PsiTop50(j) = sum(population(top50==1,j+1))/sum(population(top50==1,j));
    PsiBot50(j) = sum(population(bot50==1,j+1))/sum(population(bot50==1,j));
%     PsiTop50(j) = wmean(survival(alpha_shock_idx(:)==2,j), population(alpha_shock_idx(:)==2,j));
%     PsiBot50(j) = wmean(survival(alpha_shock_idx(:)==1,j), population(alpha_shock_idx(:)==1,j));
end

% convert to be consistent with data
MMeanTop50 = MMeanTop50(1:60);
MMeanBot50 = MMeanBot50(1:60);
% normalize expenditure to per capita income
MMeanTop50 = MMeanTop50 / Y;
MMeanBot50 = MMeanBot50 / Y;
% compute diff
MMeanDiff50 = MMeanTop50 - MMeanBot50;

% convert to be consistent with data
PsiTop50 = PsiTop50(30:70);
PsiBot50 = PsiBot50(30:70);
% compute diff
PsiDiff50 = PsiTop50 - PsiBot50;

% KYRatio
KYRatio = K/Y;

%
G_cost = G_MC + G_MD + Tss;
G_revenue = Tax + Taxss + Bequest;
GYRatio = (G_revenue - G_cost) / Y;

% construct model moments
MomentsModel = [KYRatio MCv GYRatio PsiTop50(:)' PsiBot50(:)' MMeanTop50(:)' MMeanBot50(:)'];
%{
% welfare decomposition
% as of age 0
full_cov = cov([c l hl]);
beta_cum = Dbeta.^[0:J-1];
LMean = wmean(l, population);
HMean = wmean(h, population);
psi_mean = 1-exp(-Dpsi1*HMean.^Dpsi2);
psi_prime_mean = Dpsi1*Dpsi2*HMean.^(Dpsi2-1).*exp(-Dpsi1*HMean.^Dpsi2);
psi_cum = cumprod(psi_mean);
[u,du_dc,du_dl,du_dh] = U(CMean,LMean,HMean,Dlambda,Dzeta,Drho,Dnu);
u = u+Dd;
Uc = beta_cum.*psi_cum.*du_dc;
Ul = beta_cum.*psi_cum.*du_dl;
Uh = zeros(1, J);
for j=1:J
    Uh(j) = beta_cum(j)*psi_cum(j)*du_dh(j);
    for i=j+1:J
        Uh(j) = Uh(j) + beta_cum(i)*psi_cum(i)/psi_mean(j)*psi_prime_mean(j)*u(i);
    end
end
var_U_c = Uc * full_cov(1:J,1:J) * Uc';
var_U_l = Ul * full_cov(J+1:2*J, J+1:2*J) * Ul';
var_U_h = Uh * full_cov(2*J+1:3*J, 2*J+1:3*J) * Uh';
var_U_total = [Uc Ul Uh] * full_cov * [Uc Ul Uh]';
var_U_h
var_U_h / var_U_total

% as of age 65
full_cov_Jr = cov([c(:,Jw+1:J) l(:,Jw+1:J) hl(:,Jw+1:J)]);
var_U_c_Jr = Uc(Jw+1:J) * full_cov_Jr(1:Jr,1:Jr) * Uc(Jw+1:J)';
var_U_l_Jr = Ul(Jw+1:J) * full_cov_Jr(Jr+1:2*Jr,Jr+1:2*Jr) * Ul(Jw+1:J)';
var_U_h_Jr = Uh(Jw+1:J) * full_cov_Jr(2*Jr+1:3*Jr,2*Jr+1:3*Jr) * Uh(Jw+1:J)';
var_U_total_Jr = [Uc(Jw+1:J) Ul(Jw+1:J) Uh(Jw+1:J)] * full_cov_Jr * [Uc(Jw+1:J) Ul(Jw+1:J) Uh(Jw+1:J)]';
var_U_h_Jr
var_U_h_Jr / var_U_total_Jr

% figure;
% plot(PsiTop50(1:75),'k-');
% hold on;
% plot(PsiBot50(1:75),'r.');
% hold off;
%}

%{
% NOTE(wenlan): These staff should be moved to CALIBRATION
PsiFigure = figure;
hold on;
plot([50:90],PsiTop50,'k-','LineWidth',2);
plot([50:90],PsiBot50,'r-','LineWidth',2);
plot(PsiDataAge,PsiDataTop50,'k--','LineWidth',2);
plot(PsiDataAge,PsiDataBot50,'r--','LineWidth',2);
xlabel('Age');
ylabel('Conditional survival probability');
set(findall(gcf,'type','text'),'FontSize',16,'FontName','Times New Roman')
hold off;
print(PsiFigure, 'graph/PsiModelVsData.pdf', '-dpdf');

% normalize to 2003 per capita income
% Y2003 = 39682.47; % from google
% MMeanTop50 = MMeanTop50 / Y * Y2003;
% MMeanBot50 = MMeanBot50 / Y * Y2003;
MedicalFigure = figure;
hold on;
plot([21:80],MMeanTop50,'k-','LineWidth',2);
plot([21:80],MMeanBot50,'r-','LineWidth',2);
plot(MedicalDataAge,MedicalDataTop50,'k--','LineWidth',2);
plot(MedicalDataAge,MedicalDataBot50,'r--','LineWidth',2);
xlabel('Age');
ylabel('2003 $');
set(findall(gcf,'type','text'),'FontSize',16,'FontName','Times New Roman')
hold off;
print(MedicalFigure, 'graph/MedicalModelVsData.pdf', '-dpdf');
%}

% figure;
% plot(IncMeanTop50,'k-');
% hold on;
% plot(IncMeanBot50,'r.');
% hold off;
% save('v0');
aggregate_result = v2struct(K,N,L,Tax,Taxss,Tss,Tr,Y,new_w,new_r,WMean,G_MC,G_MD,MCv,Bequest,MMeanTop50,MMeanBot50,PsiTop50,PsiBot50,MomentsModel,GYRatio,KLRatio,KYRatio);
exitflag = 0;
end


function s = wmean(x, weight)
s = sum(x.*weight)./sum(weight);
end

function v = wvar(x, weight)
s = sum(x.*weight) ./ sum(weight);
v = sum(((x-repmat(s,size(x,1),1)).^2).*weight) ./ sum(weight);
end

function wMed = weightedMedian(D,W)
D = D(:)';
W = W(:)';

% sort D;
[sortD, iD] = sort(D);
% sort M;
sortW = W(iD);
% cumsum M;
sortW = sortW / sum(sortW);
cumsumW = cumsum(W);
% find closest to 0.5 cumsum;
diff_to_50 = abs(cumsumW - 0.5);
[~,median_idx] = min(diff_to_50);

wMed = sortD(median_idx);
end
