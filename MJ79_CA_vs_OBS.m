LoadAndFilter
load('asd_jet.mat')

%% Simulate Vapor from measurements using the with Closure Assumption
Rsw18   = ((SeaWater_d18O/1000)+1)*R18SMOW;
Rsw2    = ((SeaWater_dD/1000)+1)*R2SMOW;

% Kinetic fractionation values
% k18 = 6.06;
% k2 = 20.6; 
% k2 = k18*0.88; 
% 04.05.2021 Estimates
% k18 = 5.1034;
% k2 = 16.8385;
% k2 = 0.5;
k18 = 5.2;
% k2 = 14.61;
% k2 = 3.615;
k2 = 4.3;

% Test estimates
%k18 = 7.5;
%k2 = 0.88*k18; 

% Whole data
for j = 1: length(BottomInletSYNC.d18O)
    % Test RH CALCULATION --------------------__----------------------
    e_s = (exp(77.3450 + 0.0057 .* (MeteoDataSYNC.T(j)) - (7235./(MeteoDataSYNC.T(j))))./(MeteoDataSYNC.T(j)).^8.2);
    e_a = TopInletSYNC.H2O(j)*1e-6*(MeteoDataSYNC.SLP(j)*100);
    currRH = 100*e_a/e_s;
    % End of test RH CALCULATION --------------------------------------
    R18_Vapor_CGCA = MJ79_CGCA(MeteoDataSYNC.RH(j), MeteoDataSYNC.T(j)+273.15, MeteoDataSYNC.SST(j)+273.15, Rsw18, k18, 18);
    d18V_CGCA(j, 1) = ((R18_Vapor_CGCA/R18SMOW)-1)*1e3;
    R2_Vapor_CGCA = MJ79_CGCA(MeteoDataSYNC.RH(j), MeteoDataSYNC.T(j)+273.15, MeteoDataSYNC.SST(j)+273.15, Rsw2, k2, 2);
    d2V_CGCA(j, 1) = ((R2_Vapor_CGCA/R2SMOW)-1)*1e3;
    dexV_CGCA(j, 1) = d2V_CGCA(j, 1) - 8*d18V_CGCA(j, 1);
end

% Subset data
for j = 1: length(RH)
    % Test RH CALCULATION --------------------__----------------------
    e_s = (exp(77.3450 + 0.0057 .* (AirT(j)) - (7235./(AirT(j))))./(AirT(j)).^8.2);
    e_a = Topq(j)*1e-6*(SLP(j)*100);
    currRH = 100*e_a/e_s;
    % End of test RH CALCULATION --------------------------------------    
    R18_Vapor_CGCA = MJ79_CGCA(currRH, AirT(j), SST(j), Rsw18, k18, 18);
    d18V_CGCA_subset(j, 1) = ((R18_Vapor_CGCA/R18SMOW)-1)*1e3;
    R2_Vapor_CGCA = MJ79_CGCA(currRH, AirT(j), SST(j), Rsw2, k2, 2);
    d2V_CGCA_subset(j, 1) = ((R2_Vapor_CGCA/R2SMOW)-1)*1e3;
    dexV_CGCA_subset(j, 1) = d2V_CGCA_subset(j, 1) - 8*d18V_CGCA_subset(j, 1);
end

r = Topq./Bottomq;

%% Plot d-excess vs RHSST
clf

hold on

% RHSST Top
e_s = (exp(77.3450 + 0.0057 .* (MeteoDataSYNC.SST+273.15) - (7235./(MeteoDataSYNC.SST+273.15)))./(MeteoDataSYNC.SST+273.15).^8.2);
e_a = TopInletSYNC.H2O*1e-6.*(MeteoDataSYNC.SLP*100);
RHSST_new_TOP = 100.*e_a./e_s;


scatter(RHSST_new_TOP/100, TopInletSYNC.Dexess, 32, ERA5DataSYNC.BLH, 'filled', 'LineWidth', 1.5, 'MarkerEdgeColor', [0 0 0])

MJ79bestfitMDL = fitlm(RHSST_new_TOP/100, dexV_CGCA);
fprintf('MJ79 slope: %.2f ± %.2f\n', MJ79bestfitMDL.Coefficients.Estimate(2), MJ79bestfitMDL.Coefficients.SE(2))
fprintf('MJ79 intercept: %.2f ± %.2f\n', MJ79bestfitMDL.Coefficients.Estimate(1), MJ79bestfitMDL.Coefficients.SE(1))
OBSbestfitMDL = fitlm(RHSST_new_TOP/100, TopInletSYNC.Dexess);
fprintf('OBS slope: %.2f ± %.2f\n', OBSbestfitMDL.Coefficients.Estimate(2), OBSbestfitMDL.Coefficients.SE(2))
fprintf('OBS intercept: %.2f ± %.2f\n', OBSbestfitMDL.Coefficients.Estimate(1), OBSbestfitMDL.Coefficients.SE(1))

fprintf('Average difference between MJ79 and observed d-excess %.2f\n', nanmean(dexV_CGCA - TopInletSYNC.Dexess))
asderfit = fitlm(dexV_CGCA, TopInletSYNC.Dexess);


hold on
plot(0:.01:1, predict(MJ79bestfitMDL, (0:.01:1)'),  'LineStyle', '--', 'LineWidth', 1.5, 'Color', [1 0 0])
hold off


% Graphics
colormap(asd_jet)
colorbar
grid on
box on
% xlim([-20 -5])
ax = gca;
ax.YAxis.FontSize = 13;
ax.XAxis.FontSize = 13;
ax.XLim = [0.3 1.05];
ylabel('d-excess (‰)')
xlabel('h (-)')
set(gcf, 'Color', [1 1 1]);

% Now build a regression model with BLH and CA
fprintf('------\n')
new_mod = fitlm(horzcat((dexV_CGCA), (ERA5DataSYNC.BLH)), (TopInletSYNC.Dexess));
new_mod_normal = fitlm(horzcat(normalize(dexV_CGCA), normalize(ERA5DataSYNC.BLH)), normalize(TopInletSYNC.Dexess));
predictions = predict(new_mod, horzcat((dexV_CGCA), (ERA5DataSYNC.BLH)));

fprintf('d-excess variability explained by h (CA): %.2f %%\n', ...
        100*new_mod_normal.Rsquared.Ordinary*new_mod_normal.Coefficients.Estimate(2))
fprintf('d-excess variability explained by PBLH: %.2f %%\n', ...
        100*new_mod_normal.Rsquared.Ordinary*new_mod_normal.Coefficients.Estimate(3))

