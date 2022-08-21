%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%%% 03_OoS_Test_Ibov
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
clc
warning off

addpath(genpath('/Volumes/CTRUCIOS_SD/UNICAMP/Ongoing Research/Resampling Portfolios/ResamplingBRPortfolios/Results')) % MFE Toolbox

R60 = importdata('Rport_60_ibov.csv');
R60 = R60.data;
R120 = importdata('Rport_120_ibov.csv');
R120 = R120.data;

OoS_returns = R120;

[T N] = size(OoS_returns);
k = N/3;
b = floor(T/20);
for i = 1:(k-1)
    rng(i+1234);
    Rtwo = OoS_returns(:, [1, 1 + i]);
    optimal_b = optimalblrobustVariance(Rtwo,1, 1000);
    [~, ~, ~, ~, pvalues_table(i, 1)] = robustVariance(Rtwo, 1, 5000, optimal_b);
    [pvalues_table(i, 2), ~, ~, ~] = bootInference(Rtwo, optimal_b, 5000);
    Rtwo = OoS_returns(:, [k + 1, k + 1 + i]);
    optimal_b = optimalblrobustVariance(Rtwo,1, 1000);
    [~, ~, ~, ~, pvalues_table(i, 3)] = robustVariance(Rtwo, 1, 5000, optimal_b);
    [pvalues_table(i, 4), ~, ~, ~] = bootInference(Rtwo, optimal_b, 5000);
    Rtwo = OoS_returns(:, [2*k + 1, 2*k + 1 + i]);
    optimal_b = optimalblrobustVariance(Rtwo,1, 1000);
    [~, ~, ~, ~, pvalues_table(i, 5)] = robustVariance(Rtwo, 1, 5000, optimal_b);
    [pvalues_table(i, 6), ~, ~, ~] = bootInference(Rtwo, optimal_b, 5000);
end
pvalues_table 