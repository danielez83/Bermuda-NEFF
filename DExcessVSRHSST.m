% Author:      Daniele Zannoni
% Name:        DExcessVSRHSST.m
% Description: Display relatiopship between d-excess and RHsst. Estimate
%              d-excess of water vapor at equilibrium (RHSST = 100)
%              1) Data is subsetted following the same approach for estimating
%              kinetic fractionation values, e.g. wind sector, non-rainy
%              days etc.
%              2) All computations are made using Top inlet values because
%              they are more representative of the study area and less
%              influenced by local evaporation.
%              3) The d-excess of water vapor in equilibrium with ocean
%              water is reported following measurements of Steen-Larsen
%              (personal communication), Benetti et al. 2014, Voelker et
%              al. 2015, LeGrande and Smidth 2006.
%
% Date:        Last revision 30/11/2020
tic
%% Load and filter data
LoadAndFilter

d18O_ocean = [1.06 1.12];
dD_ocean = [6.2 8.29];

%% Plot d-excess vs RHSST line
% Fit linear model (Top Inlet)
clc
dexcessRHSSTmodel = fitlm(TopInletRHSST/100, TopInletSYNC.Dexess);
xvals = (0:.01:1.2)';
yvals = predict(dexcessRHSSTmodel, xvals);
fprintf('Using all observations at top inlet: \n')
fprintf('d-excess at h = 100%% is:  %.2f ‰\n', yvals(xvals == 100));
fprintf('Slope of d vs RHSST is:        %.4f±%.4f \n', dexcessRHSSTmodel.Coefficients.Estimate(2), dexcessRHSSTmodel.Coefficients.SE(2));
fprintf('Intercept of d vs RHSST is:    %.4f±%.4f \n', dexcessRHSSTmodel.Coefficients.Estimate(1), dexcessRHSSTmodel.Coefficients.SE(1));
fprintf('Adj. R sq of d vs RHSST is:    %.2f \n', dexcessRHSSTmodel.Rsquared.Adjusted);

TopDexcess = TopdD - 8*Topd18O;
% Plot all datapoints
% scatter(MeteoDataSYNC.RHSST, TopInletSYNC.Dexess, 64, [.75 .75 .75], '+')
scatter(TopInletRHSST/100, TopInletSYNC.Dexess, 64, [.75 .75 .75], '+')
% Plot only datapoints used for the linear model
hold on
scatter(RHSSTtop/100, TopDexcess, 64, PBLH, 'filled', 'LineWidth', 1.5, 'MarkerEdgeColor', [0 0 0])
% Fit linear model (Top Inlet)
dexcessRHSSTmodelTOP = fitlm(RHSSTtop/100, TopDexcess);
% dexcessRHSSTmodelTOP = fitlm(TopInletRHSST, TopInletSYNC.Dexess);
xvals = (0:.01:1.20)';
yvals = predict(dexcessRHSSTmodelTOP, xvals);
fprintf('Using top inlet data: \n')
fprintf('d-excess at h = 1 is:  %.2f ‰\n', yvals(xvals == 100));
fprintf('Slope of d vs RHSST is:        %.4f±%.4f \n', dexcessRHSSTmodelTOP.Coefficients.Estimate(2), dexcessRHSSTmodelTOP.Coefficients.SE(2));
fprintf('Intercept of d vs RHSST is:    %.4f±%.4f \n', dexcessRHSSTmodelTOP.Coefficients.Estimate(1), dexcessRHSSTmodelTOP.Coefficients.SE(1));
fprintf('Adj. R sq of d vs RHSST is:    %.2f \n', dexcessRHSSTmodelTOP.Rsquared.Adjusted);
plot(xvals, yvals, 'LineWidth', 1.5, 'Color', [1 0 0])
%plot(xvals, yvals, 'LineWidth', 1.5, 'Color', [1 0 0])
% Plot 95% prediction boundaries
[ypred,yci] = predict(dexcessRHSSTmodelTOP, (0:120)','Alpha',0.05, 'Prediction' , 'observation');
%[ypred,yci] = predict(dexcessRHSSTmodel, (0:100)','Prediction' , 'observation');
plot(yci,'DisplayName','95% P.I', 'LineStyle', '-.', 'Color', [1 0 0])

text(.55, 35, sprintf('d-excess = (%.2f ± %.2f)*h + (%.0f ± %.0f)‰', ...
dexcessRHSSTmodelTOP.Coefficients.Estimate(2), dexcessRHSSTmodelTOP.Coefficients.SE(2), ...
dexcessRHSSTmodelTOP.Coefficients.Estimate(1), dexcessRHSSTmodelTOP.Coefficients.SE(1)), ...
'Color', [0 0 1])


%% D-excess of ocean

% Max d-excess of equilibrium vapor accounting for SST variability (1 STD)
d18O_ocean_VAPEQ_max = d18_V(max(d18O_ocean), nanmean(MeteoDataSYNC.SST)+273.15);%-nanstd(MeteoDataSYNC.SST)+273.15); % most enriched value with lower temperature
dD_ocean_VAPEQ_max = d2_V(max(dD_ocean), nanmean(MeteoDataSYNC.SST)+273.15);%-nanstd(MeteoDataSYNC.SST)+273.15);
dex_ocean_VAPEQ_max = dD_ocean_VAPEQ_max - 8.*d18O_ocean_VAPEQ_max;
% % Max d-excess of equilibrium vapor accounting for SST variability (1 STD)
d18O_ocean_VAPEQ_min = d18_V(min(d18O_ocean), nanmean(MeteoDataSYNC.SST)+273.15);%+nanstd(MeteoDataSYNC.SST)+273.15); % most depleted value with higher temperature
dD_ocean_VAPEQ_min = d2_V(min(dD_ocean), nanmean(MeteoDataSYNC.SST)+273.15);%+nanstd(MeteoDataSYNC.SST)+273.15);
dex_ocean_VAPEQ_min = dD_ocean_VAPEQ_min - 8.*d18O_ocean_VAPEQ_min;

% Mean d-excess
d18O_ocean_VAPEQ = d18_V(mean(SeaWater_d18O), nanmean(MeteoDataSYNC.SST)+273.15);
dD_ocean_VAPEQ = d2_V(mean(SeaWater_dD), nanmean(MeteoDataSYNC.SST)+273.15);
dex_ocean_VAPEQ = dD_ocean_VAPEQ - 8.*d18O_ocean_VAPEQ;
hold on
% scatter(100, mean(dex_ocean_VAPEQ), 32, [1 0 0]);
% errorbar(100, mean(dex_ocean_VAPEQ),std(dex_ocean_VAPEQ), 'Color', [1 0 0])
% line([0 120], [mean(dex_ocean_VAPEQ)+std(dex_ocean_VAPEQ) mean(dex_ocean_VAPEQ)+std(dex_ocean_VAPEQ)], 'LineStyle', '--', 'LineWidth', 1.5, 'Color', [.5 .5 .5])
% line([0 120], [mean(dex_ocean_VAPEQ)-std(dex_ocean_VAPEQ) mean(dex_ocean_VAPEQ)-std(dex_ocean_VAPEQ)], 'LineStyle', '--', 'LineWidth', 1.5, 'Color', [.5 .5 .5])
% Plot equilibrium line --------------
% x = [0 120 120 0];
%y = [dex_ocean_VAPEQ_max dex_ocean_VAPEQ_max ...
%     dex_ocean_VAPEQ_min dex_ocean_VAPEQ_min];
% patch(x,y, [0 0 1], 'EdgeColor', [0 0 0], 'FaceAlpha',.3)
scatter(1, dex_ocean_VAPEQ, 32, [0.0745 0.6235 1.0000], 'LineWidth', 2);
errorbar(1, dex_ocean_VAPEQ,(dex_ocean_VAPEQ_max - dex_ocean_VAPEQ_min), 'Color', [0.0745 0.6235 1.0000], 'LineWidth', 2)
%% Set graphix
colorbar
xlim([min(RHSST/100)-0.10 1.05])
ylim([-5 40])
ax = gca;
%ax.XLabel.String = 'RH_{SST} (%)';
ax.XLabel.String = 'h (-)';
ax.XLabel.FontSize = 12;
ax.XLabel.FontWeight = 'bold';
ax.YLabel.String = 'd-excess (‰)';
ax.YLabel.FontSize = 12;
ax.YLabel.FontWeight = 'bold';
box on

