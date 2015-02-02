function params = MOMENTS(params)
% Inferring Health Inequality
% Read in moments
% Author: Wenlan Luo at Georgetown University, luowenlan at gmail.com

v2struct(params);
clear params;

% survival probability
PsiData = csvread('MOMENTS/psi_age_inc.csv');
PsiDataAge = PsiData(:,1);
PsiDataBot50 = PsiData(:,2);
PsiDataTop50 = PsiData(:,3);
% convert to diff
PsiDataDiff50 = PsiDataTop50 - PsiDataBot50;

% medical expenditure
MedicalData = csvread('MOMENTS/medical_by_age_inc.csv');
MedicalDataAge = MedicalData(:,1);
MedicalDataBot50 = MedicalData(:,2);
MedicalDataTop50 = MedicalData(:,3);
% normalize by per capita income
MedicalDataBot50 = MedicalDataBot50/Y2003;
MedicalDataTop50 = MedicalDataTop50/Y2003;
% convert to difference
MedicalDataDiff50 = MedicalDataTop50 - MedicalDataBot50;

% KY ratio
KYRatioData = 2.7;

% cv of medical expenditure
MCvData = 2.1340;

% mandatory expenditure as a fraction of GDP
GYRatioData = 0.09;

% formulate data moments
MomentsData = [KYRatioData MCvData GYRatioData PsiDataTop50(:)' PsiDataBot50(:)' MedicalDataTop50(:)' MedicalDataBot50(:)'];
MomentsWeight = ones(size(MomentsData));
% MomentsWeight(1:3) = 100;

params = v2struct();

