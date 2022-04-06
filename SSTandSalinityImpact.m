%%
LoadAndFilter

load('../Matlab data/SST_SYNC_OSTIA1x1deg.mat')
load('../Matlab data/Salinity_SYNC.mat')

%% Colors and symbols
HogReefCOL = 'r';
CrescentReefCOL = 'g';
StGeorgeCOL = 'b';
BATSoceanCOL = [0 0.4471 0.7412];
BATSoceanSYM = '^';
MkrSyze = 16;

%% Make data subset
SST_OSTIA           = SST_SYNC.SST_OSTIA(GoodIndexes);
SST_Cresceng        = SST_SYNC.SST_Crescent(GoodIndexes);
SST_Hog             = SST_SYNC.SST_HOG(GoodIndexes);
SST_StGeorge        = SST_SYNC.SST_StGeorge(GoodIndexes);

Salinity_BATS       = Sal_DataSYNC.goodSALBATS(GoodIndexes);
Salinity_Cresceng   = Sal_DataSYNC.goodSALCrescent(GoodIndexes);
Salinity_Hog        = Sal_DataSYNC.goodSALHOG(GoodIndexes);

%% Salinity to isotopic ratio conversion function
% Sal2d18O    = @(sal) 0.39*sal-13.1; % From BATS s-delta relationship
% Sal2dD      = @(sal) 8.64*sal-307.6; % From BATS s-delta relationship
Sal2d18O    = @(sal) 0.32*sal-10.5; % From Benetti et al. (2017) s-delta relationship
Sal2dD      = @(sal) 2.04*sal-65.8; % From Benetti et al. (2017) s-delta relationship
% Sal2d18O    = @(sal) 0.3605*sal-11.849; % From Benetti et al. (2017) + BATS 2012 data s-delta relationship, see OceanIsotopicCompostion.xlsx Tab Benetti et al 2017
% Sal2dD      = @(sal) 2.139*sal-69.539; % From Benetti et al. (2017) + BATS 2012 data  s-delta relationship, see OceanIsotopicCompostion.xlsx Tab Benetti et al 2017

%% Tempearure for d18O and dD, deviation from OSTIA
x_values_deltas = -2:.01:+2;
% histogram(d18_V(SeaWater_d18O, SST_Hog + 273.15) - d18_V(SeaWater_d18O, SST_OSTIA + 273.15), 'Normalization', 'pdf' )
% histogram(d2_V(SeaWater_dD, SST_Hog + 273.15) - d2_V(SeaWater_dD, SST_OSTIA + 273.15), 'Normalization', 'pdf' )
pd = fitdist(d18_V(SeaWater_d18O, SST_Hog + 273.15) - d18_V(SeaWater_d18O, SST_OSTIA + 273.15),'Kernel', 'Kernel','normal', 'Bandwidth', .05 );
pdfHOG = pdf(pd,x_values_deltas);
yyaxis left
plot(x_values_deltas, pdfHOG, HogReefCOL, 'LineStyle', '--', 'LineWidth', 1.5)
hold on
yyaxis right
pd = fitdist(d2_V(SeaWater_dD, SST_Hog + 273.15) - d2_V(SeaWater_dD, SST_OSTIA + 273.15),'Kernel', 'Kernel','normal', 'Bandwidth', .5 );
pdfHOG = pdf(pd,x_values_deltas);
plot(x_values_deltas, pdfHOG, HogReefCOL, 'LineStyle', '-', 'LineWidth', 1.5)

% histogram(d18_V(SeaWater_d18O, SST_Cresceng + 273.15) - d18_V(SeaWater_d18O, SST_OSTIA + 273.15), 'Normalization', 'pdf' )
pd = fitdist(d18_V(SeaWater_d18O, SST_Cresceng + 273.15) - d18_V(SeaWater_d18O, SST_OSTIA + 273.15),'Kernel', 'Kernel','normal', 'Bandwidth', .05 );
pdfCRESCENT = pdf(pd,x_values_deltas);
yyaxis left
plot(x_values_deltas, pdfCRESCENT, CrescentReefCOL, 'LineStyle', '--', 'LineWidth', 1.5)
pd = fitdist(d2_V(SeaWater_dD, SST_Cresceng + 273.15) - d2_V(SeaWater_dD, SST_OSTIA + 273.15),'Kernel', 'Kernel','normal', 'Bandwidth', .5 );
pdfCRESCENT = pdf(pd,x_values_deltas);
yyaxis right
plot(x_values_deltas, pdfCRESCENT, CrescentReefCOL, 'LineStyle', '-', 'LineWidth', 1.5)


% histogram(d18_V(SeaWater_d18O, SST_StGeorge + 273.15) - d18_V(SeaWater_d18O, SST_OSTIA + 273.15), 'Normalization', 'pdf' )
pd = fitdist(d18_V(SeaWater_d18O, SST_StGeorge + 273.15) - d18_V(SeaWater_d18O, SST_OSTIA + 273.15),'Kernel', 'Kernel','normal', 'Bandwidth', .05 );
pdfSTGEORGE = pdf(pd,x_values_deltas);
yyaxis left
plot(x_values_deltas, pdfSTGEORGE, StGeorgeCOL, 'LineStyle', '--', 'LineWidth', 1.5)
pd = fitdist(d2_V(SeaWater_dD, SST_StGeorge + 273.15) - d2_V(SeaWater_dD, SST_OSTIA + 273.15),'Kernel', 'Kernel','normal', 'Bandwidth', .5 );
pdfSTGEORGE = pdf(pd,x_values_deltas);
yyaxis right
plot(x_values_deltas, pdfSTGEORGE, StGeorgeCOL, 'LineStyle', '-', 'LineWidth', 1.5)

hold off
title('Deviation from OSTIA equilibrium water vapor')
% Shared axis option
ax = gca;
ax.YAxis(1).Label.String = "PDF \delta^{18}O";
ax.YAxis(1).Color = [0 0 0];
ax.YAxis(2).Label.String = "PDF \deltaD";
ax.YAxis(2).Color = [0 0 0];
ax.XAxis(1).Label.String = '\delta* values (‰)';
legend(["Hog Reef \delta^{18}","Crescent Reef \delta^{18}","St. Georgef \delta^{18}",...
        "Hog Reef \deltaD","Crescent Reef \deltaD","St. Georgef \deltaD"])
grid on
% %% Tempearure for dD
% histogram(d2_V(SeaWater_dD, SST_Hog + 273.15) - d2_V(SeaWater_dD, SST_OSTIA + 273.15), 'Normalization', 'pdf' )
% hold on
% histogram(d2_V(SeaWater_dD, SST_Cresceng + 273.15) - d2_V(SeaWater_dD, SST_OSTIA + 273.15), 'Normalization', 'pdf' )
% histogram(d2_V(SeaWater_dD, SST_StGeorge + 273.15) - d2_V(SeaWater_dD, SST_OSTIA + 273.15), 'Normalization', 'pdf' )
% hold off
fprintf('Mean deviation using Hog instead of OSTIA is d18O = %.2f‰, dD = %.2f‰\n', ...
        mean(d18_V(SeaWater_d18O, SST_Hog + 273.15) - d18_V(SeaWater_d18O, SST_OSTIA + 273.15)), ...
        mean(d2_V(SeaWater_dD, SST_Hog + 273.15) - d2_V(SeaWater_dD, SST_OSTIA + 273.15)));
fprintf('Mean deviation using Crescent instead of OSTIA is d18O = %.2f‰, dD = %.2f‰\n', ...
        mean(d18_V(SeaWater_d18O, SST_Cresceng + 273.15) - d18_V(SeaWater_d18O, SST_OSTIA + 273.15)), ...
        mean(d2_V(SeaWater_dD, SST_Cresceng + 273.15) - d2_V(SeaWater_dD, SST_OSTIA + 273.15)));    
fprintf('Mean deviation using St.George instead of OSTIA is d18O = %.2f‰, dD = %.2f‰\n', ...
        mean(d18_V(SeaWater_d18O, SST_StGeorge + 273.15) - d18_V(SeaWater_d18O, SST_OSTIA + 273.15)), ...
        mean(d2_V(SeaWater_dD, SST_StGeorge + 273.15) - d2_V(SeaWater_dD, SST_OSTIA + 273.15)));        
%% Salinity for d18O and dD
x_values_deltas = -2:.01:+2;
% histogram(Sal2d18O(Salinity_Hog) - Sal2d18O(Salinity_BATS), 'Normalization', 'pdf' )
pd = fitdist(Sal2d18O(Salinity_Hog) - Sal2d18O(Salinity_BATS),'Kernel', 'Kernel','normal', 'Bandwidth', .05 );
pdfHOG_sal = pdf(pd,x_values_deltas);
yyaxis left
plot(x_values_deltas, pdfHOG_sal, HogReefCOL, 'LineStyle', '--', 'LineWidth', 1.5)
hold on
yyaxis right
pd = fitdist(Sal2dD(Salinity_Hog) - Sal2dD(Salinity_BATS),'Kernel', 'Kernel','normal', 'Bandwidth', .5 );
pdfHOG_sal = pdf(pd,x_values_deltas);
plot(x_values_deltas, pdfHOG_sal, HogReefCOL, 'LineStyle', '-', 'LineWidth', 1.5)


x_values_deltas = -2:.01:+2;
% histogram(Sal2d18O(Salinity_Hog) - Sal2d18O(Salinity_BATS), 'Normalization', 'pdf' )
pd = fitdist(Sal2d18O(Salinity_Cresceng) - Sal2d18O(Salinity_BATS),'Kernel', 'Kernel','normal', 'Bandwidth', .05 );
pdfCRESCENT_sal = pdf(pd,x_values_deltas);
yyaxis left
plot(x_values_deltas, pdfCRESCENT_sal, CrescentReefCOL, 'LineStyle', '--', 'LineWidth', 1.5)
hold on
yyaxis right
pd = fitdist(Sal2dD(Salinity_Cresceng) - Sal2dD(Salinity_BATS),'Kernel', 'Kernel','normal', 'Bandwidth', .5 );
pdfCRESCENT_sal = pdf(pd,x_values_deltas);
plot(x_values_deltas, pdfCRESCENT_sal, CrescentReefCOL, 'LineStyle', '-', 'LineWidth', 1.5)
title('Ocean composition, deviation from BATS')
% Shared axis option
ax = gca;
ax.YAxis(1).Label.String = "PDF \delta^{18}O";
ax.YAxis(1).Color = [0 0 0];
ax.YAxis(2).Label.String = "PDF \deltaD";
ax.YAxis(2).Color = [0 0 0];
ax.XAxis(1).Label.String = '\delta values (‰)';
legend(["Hog Reef \delta^{18}","Crescent Reef \delta^{18}",...
        "Hog Reef \deltaD","Crescent Reef \deltaD"])
grid on

fprintf('Deviation of ocean composition using Hog instead of mean value is d18O = %.2f‰, dD = %.2f‰\n', ...
        nanmean(Sal2d18O(Salinity_Hog) - Sal2d18O(Salinity_BATS)), ...
        nanmean(Sal2dD(Salinity_Hog) - Sal2dD(Salinity_BATS)));
fprintf('Deviation of ocean composition using Crescent instead of mean value is d18O = %.2f‰, dD = %.2f‰\n', ...
        nanmean(Sal2d18O(Salinity_Cresceng) - Sal2d18O(Salinity_BATS)), ...
        nanmean(Sal2dD(Salinity_Cresceng) - Sal2dD(Salinity_BATS)));
% %% Salinity for dD
% histogram(Sal2dD(Salinity_Hog) - Sal2dD(Salinity_BATS), 'Normalization', 'pdf' )
% hold on
% histogram(Sal2dD(Salinity_Cresceng) - Sal2dD(Salinity_BATS), 'Normalization', 'pdf' )
% hold off
%%
