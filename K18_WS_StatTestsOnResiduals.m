% Author:      Daniele Zannoni
% Name:        K18_WS_dependency.m
% Description: Plot k/18 as a function of WS
% Date:        Last revision 07/12/2020

%% Run script
% DetermineBestKineticXIAO2017_18Oand2H_minimizationNOSPRAY_HF
%DetermineBestKineticNOONE2012_18Oand2H_minimizationNOSPRAY_HF
% trick?
%kinetic_2_1 = 0.88*kinetic_18_16;
K_WS_dependency_V4 % run usual script first
kinetic_18_16 = BK18_KP;
kinetic_2_1 = BK2_KP;
clf
%% bin estimations
WSbinner = 0:.5:12;
WSbinnerVAls = WSbinner+.25;
WSbinnerVAls(end) = [];
kinetic_18_16_WS = zeros(1, length(WSbinner)-1);
kinetic_2_1_WS = zeros(1, length(WSbinner)-1);
for j = 1 : length(WSbinner) - 1
    indexesOI =  WS>=WSbinner(j) & WS<WSbinner(j+1);
    kinetic_18_16_WS(j) = mean(kinetic_18_16(indexesOI));
    SE_kinetic_18_16_WS(j) = std(kinetic_18_16(indexesOI))/sqrt(sum(indexesOI));
    kinetic_2_1_WS(j) = mean(kinetic_2_1(indexesOI));
    SE_kinetic_2_1_WS(j) = std(kinetic_2_1(indexesOI))/sqrt(sum(indexesOI));
end
%% Calculate k value with MJ79 parametrization
% Estimate kinetic fractionation from model, only smooth regime
SimWS_LOWER = 0.1:0.1:6;
SimWS_UPPER = 6.1:0.1:15;
ustar = MJ79_ustar(SimWS_LOWER, 10);
k_smooth_LOWER = MJ79_k(ustar, 1000, 18, 'smooth');
ustar = MJ79_ustar(SimWS_UPPER, 10);
k_smooth_UPPER = MJ79_k(ustar, 1000, 18, 'smooth');

% Estimate kinetic fractionation from model, only smooth regime
ustar = MJ79_ustar(SimWS_LOWER, 10);
k_rough_LOWER = MJ79_k(ustar, 1000, 18, 'rough');
ustar = MJ79_ustar(SimWS_UPPER, 10);
k_rough_UPPER = MJ79_k(ustar, 1000, 18, 'rough');

% Estimate kinetic fractionation from model, smooth and rough regime
ustar = MJ79_ustar(SimWS_LOWER, 10);
k_SandR_LOWER = MJ79_k(ustar, 1000, 18, 'smooth');
ustar = MJ79_ustar(SimWS_UPPER, 10);
k_SandR_UPPER = MJ79_k(ustar, 1000, 18, 'rough');

%% Calculate correlation
clc
indexer = WSbinnerVAls > .5 & WSbinnerVAls < 10;
[rho, pval] = corr(WSbinnerVAls(indexer)', kinetic_18_16_WS_KP(indexer)');
fprintf('Correlation between k18 and WS: %.2f (pval =%.5f)\n', rho, pval)
[rho, pval] = corr(WSbinnerVAls(indexer)', kinetic_2_1_WS_KP(indexer)');
fprintf('Correlation between k2 and WS: %.2f (pval =%.5f)\n', rho, pval)

%% Calculate k18 for smooth regime
WS_Threshold = 6;
ustar = MJ79_ustar(WS, 10);
k18_smooth_OBS = MJ79_k(ustar, 1000, 18, 'smooth');
k2_smooth_OBS = MJ79_k(ustar, 1000, 2, 'smooth');

alphalevel = 0.01;

for j = 1 : length(WSbinner) - 1
    indexesOI =  WS>=WSbinner(j) & WS<WSbinner(j+1);
    kinetic_18_16_MJ79(j) = mean(k18_smooth_OBS(indexesOI));
    SE_kinetic_18_16_MJ79(j) = std(k18_smooth_OBS(indexesOI))/sqrt(sum(indexesOI));
    kinetic_2_1_MJ79(j) = mean(k2_smooth_OBS(indexesOI));
    SE_kinetic_2_1_MJ79(j) = std(k2_smooth_OBS(indexesOI))/sqrt(sum(indexesOI));
end
[h, p] = kstest(kinetic_18_16_MJ79(WSbinnerVAls<WS_Threshold) - kinetic_18_16_WS(WSbinnerVAls<WS_Threshold));
if h
    fprintf('Kolmogorov-Smirnov test: 18-O residuals not normally distributed with pval: %.4f\n', p)
end
[h, p] = kstest(kinetic_2_1_MJ79(WSbinnerVAls<WS_Threshold) - kinetic_2_1_WS(WSbinnerVAls<WS_Threshold));
if h
    fprintf('Kolmogorov-Smirnov test: 2-H residuals not normally distributed with pval: %.4f\n', p)
end
% [h, p] = lillietest(kinetic_18_16_MJ79(WSbinnerVAls<WS_Threshold) - kinetic_18_16_WS(WSbinnerVAls<WS_Threshold));
% if h
%     fprintf('Lilliefors test: 18-O residuals not normally distributed with pva: %.4f\n', p)
% end
%  jbtest, goodness-of-fit test of whether sample data have the skewness and kurtosis matching a normal distribution.
[h, p] = swtest(kinetic_18_16_MJ79(WSbinnerVAls<WS_Threshold) - kinetic_18_16_WS(WSbinnerVAls<WS_Threshold));
if h
    fprintf('Shapiro-Wilk test: 18-O residuals not normally distributed with pval: %.4f\n', p)
end
[h, p] = swtest(kinetic_2_1_MJ79(WSbinnerVAls<WS_Threshold) - kinetic_2_1_WS(WSbinnerVAls<WS_Threshold));
if h
    fprintf('Shapiro-Wilk test: 2-H residuals not normally distributed with pval: %.4f\n', p)
end

[h, p] = ttest(kinetic_18_16_MJ79(WSbinnerVAls<WS_Threshold) - kinetic_18_16_WS(WSbinnerVAls<WS_Threshold));
if h
    fprintf('T-test: 18-O residuals not normally distributed with pval: %.4f\n', p)
end
[h, p] = ttest(kinetic_2_1_MJ79(WSbinnerVAls<WS_Threshold) - kinetic_2_1_WS(WSbinnerVAls<WS_Threshold));
if h
    fprintf('T-test: 2-H residuals not normally distributed with pval: %.4f\n', p)
end

fprintf('Mean difference for 18-O: %.3f‰\n', ...
        nanmean(kinetic_18_16_MJ79(WSbinnerVAls<WS_Threshold) - kinetic_18_16_WS(WSbinnerVAls<WS_Threshold)))
fprintf('Mean difference for 2-H: %.3f‰\n', ...
        nanmean(kinetic_2_1_MJ79(WSbinnerVAls<WS_Threshold) - kinetic_2_1_WS(WSbinnerVAls<WS_Threshold)))    

%% Calculate k18 for rough regime
clc
WS_Threshold = 6;
ustar = MJ79_ustar(WS, 10);
k18_smooth_OBS = MJ79_k(ustar, 1000, 18, 'rough');
k2_smooth_OBS = MJ79_k(ustar, 1000, 2, 'rough');

alphalevel = 0.01;

for j = 1 : length(WSbinner) - 1
    indexesOI =  WS>=WSbinner(j) & WS<WSbinner(j+1);
    kinetic_18_16_MJ79(j) = mean(k18_smooth_OBS(indexesOI));
    SE_kinetic_18_16_MJ79(j) = std(k18_smooth_OBS(indexesOI))/sqrt(sum(indexesOI));
    kinetic_2_1_MJ79(j) = mean(k2_smooth_OBS(indexesOI));
    SE_kinetic_2_1_MJ79(j) = std(k2_smooth_OBS(indexesOI))/sqrt(sum(indexesOI));
end
[h, p] = kstest(kinetic_18_16_MJ79(WSbinnerVAls>WS_Threshold) - kinetic_18_16_WS(WSbinnerVAls>WS_Threshold));
if h
    fprintf('Kolmogorov-Smirnov test: 18-O residuals not normally distributed with pval: %.4f\n', p)
end
[h, p] = kstest(kinetic_2_1_MJ79(WSbinnerVAls>WS_Threshold) - kinetic_2_1_WS(WSbinnerVAls>WS_Threshold));
if h
    fprintf('Kolmogorov-Smirnov test: 2-H residuals not normally distributed with pval: %.4f\n', p)
end
% [h, p] = lillietest(kinetic_18_16_MJ79(WSbinnerVAls<WS_Threshold) - kinetic_18_16_WS(WSbinnerVAls<WS_Threshold));
% if h
%     fprintf('Lilliefors test: 18-O residuals not normally distributed with pva: %.4f\n', p)
% end
%  jbtest, goodness-of-fit test of whether sample data have the skewness and kurtosis matching a normal distribution.
[h, p] = swtest(kinetic_18_16_MJ79(WSbinnerVAls>WS_Threshold) - kinetic_18_16_WS(WSbinnerVAls>WS_Threshold));
if h
    fprintf('Shapiro-Wilk test: 18-O residuals not normally distributed with pval: %.4f\n', p)
end
[h, p] = swtest(kinetic_2_1_MJ79(WSbinnerVAls>WS_Threshold) - kinetic_2_1_WS(WSbinnerVAls>WS_Threshold));
if h
    fprintf('Shapiro-Wilk test: 2-H residuals not normally distributed with pval: %.4f\n', p)
end

[h, p] = ttest(kinetic_18_16_MJ79(WSbinnerVAls>WS_Threshold) - kinetic_18_16_WS(WSbinnerVAls>WS_Threshold));
if h
    fprintf('T-test: 18-O residuals not normally distributed with pval: %.4f\n', p)
end
[h, p] = ttest(kinetic_2_1_MJ79(WSbinnerVAls>WS_Threshold) - kinetic_2_1_WS(WSbinnerVAls>WS_Threshold));
if h
    fprintf('T-test: 2-H residuals not normally distributed with pval: %.4f\n', p)
end

fprintf('Mean difference for 18-O: %.3f‰\n', ...
        nanmean(kinetic_18_16_MJ79(WSbinnerVAls(2:20)>WS_Threshold) - kinetic_18_16_WS(WSbinnerVAls(2:20)>WS_Threshold)))
fprintf('Mean difference for 2-H: %.3f‰\n', ...
        nanmean(kinetic_2_1_MJ79(WSbinnerVAls(2:20)>WS_Threshold) - kinetic_2_1_WS(WSbinnerVAls(2:20)>WS_Threshold)))    

goodindexes = WSbinnerVAls>6 & WSbinnerVAls<10;
fprintf('Mean difference for 18-O in 6 - 10m/s range: %.3f‰\n', ...
        nanmean(kinetic_18_16_MJ79(goodindexes) - kinetic_18_16_WS(goodindexes)))
fprintf('Mean difference for 2-H in 6 - 10m/s range: %.3f‰\n', ...
        nanmean(kinetic_2_1_MJ79(goodindexes) - kinetic_2_1_WS(goodindexes)))    
fprintf('Reltive difference for 18O is %.2f %%\n',...
    abs(100*nanmean(kinetic_18_16_MJ79(goodindexes) - kinetic_18_16_WS(goodindexes))./nanmean(kinetic_18_16_MJ79(goodindexes))))
%% Calculate k18 for smooth and rough regime
clc
WS_Threshold = 9000;
k_smoothandrough_OBS = nan(length(k18_smooth_OBS), 1);

ustar = MJ79_ustar(WS<6, 10);
k18_smooth_OBS = MJ79_k(ustar, 1000, 18, 'smooth');
k2_smooth_OBS = MJ79_k(ustar, 1000, 2, 'smooth');
ustar = MJ79_ustar(WS>=6, 10);
k18_rough_OBS = MJ79_k(ustar, 1000, 18, 'rough');
k2_rough_OBS = MJ79_k(ustar, 1000, 2, 'rough');
k18_smoothandrough_OBS = k18_smooth_OBS + k18_rough_OBS;
k2_smoothandrough_OBS = k2_smooth_OBS + k2_rough_OBS;

for j = 1 : length(WSbinner) - 1
    indexesOI =  WS>=WSbinner(j) & WS<WSbinner(j+1);
    kinetic_18_16_MJ79(j) = mean(k18_smoothandrough_OBS(indexesOI));
    SE_kinetic_18_16_MJ79(j) = std(k18_smoothandrough_OBS(indexesOI))/sqrt(sum(indexesOI));
    kinetic_2_1_MJ79(j) = mean(k2_smoothandrough_OBS(indexesOI));
    SE_kinetic_2_1_MJ79(j) = std(k2_smoothandrough_OBS(indexesOI))/sqrt(sum(indexesOI));
end

[h, p] = kstest(kinetic_18_16_MJ79(WSbinnerVAls<WS_Threshold) - kinetic_18_16_WS(WSbinnerVAls<WS_Threshold));
if h
    fprintf('Kolmogorov-Smirnov test: 18-O residuals not normally distributed with pval: %.4f\n', p)
end
[h, p] = kstest(kinetic_2_1_MJ79(WSbinnerVAls<WS_Threshold) - kinetic_2_1_WS(WSbinnerVAls<WS_Threshold));
if h
    fprintf('Kolmogorov-Smirnov test: 2-H residuals not normally distributed with pval: %.4f\n', p)
end
% [h, p] = lillietest(kinetic_18_16_MJ79(WSbinnerVAls<WS_Threshold) - kinetic_18_16_WS(WSbinnerVAls<WS_Threshold));
% if h
%     fprintf('Lilliefors test: 18-O residuals not normally distributed with pva: %.4f\n', p)
% end
%  jbtest, goodness-of-fit test of whether sample data have the skewness and kurtosis matching a normal distribution.
[h, p] = swtest(kinetic_18_16_MJ79(WSbinnerVAls<WS_Threshold) - kinetic_18_16_WS(WSbinnerVAls<WS_Threshold));
if h
    fprintf('Shapiro-Wilk test: 18-O residuals not normally distributed with pval: %.4f\n', p)
end
[h, p] = swtest(kinetic_2_1_MJ79(WSbinnerVAls<WS_Threshold) - kinetic_2_1_WS(WSbinnerVAls<WS_Threshold));
if h
    fprintf('Shapiro-Wilk test: 2-H residuals not normally distributed with pval: %.4f\n', p)
end

[h, p] = ttest(kinetic_18_16_MJ79(WSbinnerVAls<WS_Threshold) - kinetic_18_16_WS(WSbinnerVAls<WS_Threshold));
if h
    fprintf('T-test: 18-O residuals not normally distributed with pval: %.4f\n', p)
end
[h, p] = ttest(kinetic_2_1_MJ79(WSbinnerVAls<WS_Threshold) - kinetic_2_1_WS(WSbinnerVAls<WS_Threshold));
if h
    fprintf('T-test: 2-H residuals not normally distributed with pval: %.4f\n', p)
end

fprintf('Mean difference for 18-O: %.3f‰\n', ...
        nanmean(kinetic_18_16_MJ79(WSbinnerVAls<WS_Threshold) - kinetic_18_16_WS(WSbinnerVAls<WS_Threshold)))
fprintf('Mean difference for 2-H: %.3f‰\n', ...
        nanmean(kinetic_2_1_MJ79(WSbinnerVAls<WS_Threshold) - kinetic_2_1_WS(WSbinnerVAls<WS_Threshold)))    