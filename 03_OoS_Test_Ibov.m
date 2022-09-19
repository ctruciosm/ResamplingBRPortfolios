%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%%% Out-of-Sample tests for SD and SR 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
clc
warning off
addpath(genpath('/Volumes/CTRUCIOS_SD/UNICAMP/Ongoing Research/Resampling Portfolios/ResamplingBRPortfolios/Results')) % MFE Toolbox
R60 = importdata('Rport_60_ibov.csv');
R60 = R60.data;
R120 = importdata('Rport_120_ibov.csv');
R120 = R120.data;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                 T = 60                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
OoS_returns = R60;
[T N] = size(OoS_returns);
k = N/3;
for i = 1:(k-1)
    rng(i+1234);
    Rtwo = OoS_returns(:, [1, 1 + i]);
    optimal_b = optimalblrobustVariance(Rtwo,1, 1000, 199, 'G', [1 3 6 10 15]);
    [~, log_diff(i, 1), hac_sd(i, 1), ~, pvalues_table(i, 1)] = robustVariance_onesided(Rtwo, 1, 5000, optimal_b, 'G');
    [pvalues_table(i, 2), ~, ~, ~] = bootInference(Rtwo, optimal_b, 5000, 'G');
    Rtwo = OoS_returns(:, [k + 1, k + 1 + i]);
    optimal_b = optimalblrobustVariance(Rtwo,1, 1000, 199, 'G', [1 3 6 10 15]);
    [~, log_diff(i, 2), hac_sd(i, 2), ~, pvalues_table(i, 3)] = robustVariance_onesided(Rtwo, 1, 5000, optimal_b, 'G');
    [pvalues_table(i, 4), ~, ~, ~] = bootInference(Rtwo, optimal_b, 5000, 'G');
    Rtwo = OoS_returns(:, [2*k + 1, 2*k + 1 + i]);
    optimal_b = optimalblrobustVariance(Rtwo,1, 1000, 199, 'G', [1 3 6 10 15]);
    [~, log_diff(i, 3), hac_sd(i, 3), ~, pvalues_table(i, 5)] = robustVariance_onesided(Rtwo, 1, 5000, optimal_b, 'G');
    [pvalues_table(i, 6), ~, ~, ~] = bootInference(Rtwo, optimal_b, 5000, 'G');
end
pvalues60 = pvalues_table;
log_diff60 = log_diff;
hac_sd60 = hac_sd;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                 T = 120                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
OoS_returns = R120;
[T N] = size(OoS_returns);
k = N/3;
for i = 1:(k-1)
    rng(i+1234);
    Rtwo = OoS_returns(:, [1, 1 + i]);
    optimal_b = optimalblrobustVariance(Rtwo,1, 1000, 199, 'G', [1 3 6 10 15]);
    [~, log_diff(i, 1), hac_sd(i, 1), ~, pvalues_table(i, 1)] = robustVariance_onesided(Rtwo, 1, 5000, optimal_b, 'G');
    [pvalues_table(i, 2), ~, ~, ~] = bootInference(Rtwo, optimal_b, 5000, 'G');
    Rtwo = OoS_returns(:, [k + 1, k + 1 + i]);
    optimal_b = optimalblrobustVariance(Rtwo,1, 1000, 199, 'G', [1 3 6 10 15]);
    [~, log_diff(i, 2), hac_sd(i, 2), ~, pvalues_table(i, 3)] = robustVariance_onesided(Rtwo, 1, 5000, optimal_b, 'G');
    [pvalues_table(i, 4), ~, ~, ~] = bootInference(Rtwo, optimal_b, 5000, 'G');
    Rtwo = OoS_returns(:, [2*k + 1, 2*k + 1 + i]);
    optimal_b = optimalblrobustVariance(Rtwo,1, 1000, 199, 'G', [1 3 6 10 15]);
    [~, log_diff(i, 3), hac_sd(i, 3), ~, pvalues_table(i, 5)] = robustVariance_onesided(Rtwo, 1, 5000, optimal_b, 'G');
    [pvalues_table(i, 6), ~, ~, ~] = bootInference(Rtwo, optimal_b, 5000, 'G');
end
pvalues_table 
log_diff
hac_sd


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Tangency Portfolios pairwise comparison %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
OoS_returns = R60;
[T N] = size(OoS_returns);
K = N/3
for i = 1:(K-1)
    for j = (i+1):K
        rng(i + j +1234);
        optimal_b = optimalblrobustVariance(OoS_returns(:, [K + i, K + j]),1, 1000, 199, 'G', [1 3 6 10 15]);
        [~, ~, ~, ~, pvalues_table_sd_tp60(i, j)] = robustVariance_onesided(OoS_returns(:, [K + i, K + j]), 1, 5000, optimal_b, 'G');
        [~, ~, ~, ~, pvalues_table_sd_tp60(j, i)] = robustVariance_onesided(OoS_returns(:, [K + j, K + i]), 1, 5000, optimal_b, 'G');
    end
end
OoS_returns = R120;
[T N] = size(OoS_returns);
K = N/3
for i = 1:(K-1)
    for j = (i+1):K
        rng(i + j +1234);
        optimal_b = optimalblrobustVariance(OoS_returns(:, [K + i, K + j]),1, 1000, 199, 'G', [1 3 6 10 15]);
        [~, ~, ~, ~, pvalues_table_sd_tp120(i, j)] = robustVariance_onesided(OoS_returns(:, [K + i, K + j]), 1, 5000, optimal_b, 'G');
        [~, ~, ~, ~, pvalues_table_sd_tp120(j, i)] = robustVariance_onesided(OoS_returns(:, [K + j, K + i]), 1, 5000, optimal_b, 'G');
    end
end

