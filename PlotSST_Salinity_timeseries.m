% Author:      Daniele Zannoni
% Name:        PlotSST_Salinity_timeseries.m
% Description: Plot k/18 as a function of WS
% Date:        Last revision 12/05/2021


%% Import data
LoadAndFilter
data_folder             = '../GIT Data/';
NOAA_BATS_SYNC          = readtable(strcat(data_folder, 'NOAA_BATS_SYNC.txt'));
HogReef                 = readtable(strcat(data_folder, 'HogReef_ASYNC.txt'));
Crescent                = readtable(strcat(data_folder, 'Crescent_ASYNC.txt'));
BATS                    = readtable(strcat(data_folder, 'BATS_ASYNC.txt'));


% Geolocations
HogReefPOS = [32.46, -64.83]; % https://www.ncei.noaa.gov/access/ocean-carbon-data-system/oceans/Moorings/Hog_Reef.html
CrescentReefPOS = [32.40, -64.79]; % https://www.ncei.noaa.gov/access/ocean-carbon-data-system/oceans/Moorings/Crescent_Reef.html
BATSoceanPOS = [31.69 -64.17]; % approximated from here http://bats.bios.edu/about/
StGeorgePOS = [32.38, -64.67]; % St.George harbor


%% Isotopic values
%SeaWater_d18O          = 1.26;       % [‰], from excel file "OceanIsotopicComposition.xlsx"
%SeaWater_dD            = 8.24;       % [‰], from excel file "OceanIsotopicComposition.xlsx"

%% Colors and symbols
HogReefCOL = 'r';
CrescentReefCOL = 'g';
StGeorgeCOL = 'b';
BATSoceanCOL = 'y'; %[0 0.4471 0.7412];
BATSoceanSYM = '^';
MkrSyze = 16;

%% Plot positions
BoxCoordinates = [31.5 32.5; -65 -64]; % Latitude limits; Longitude limits

h = worldmap(BoxCoordinates(1,:) ,BoxCoordinates(2,:));
% Change colors
p = findobj(h,'type','patch'); % Find background
%set(p,'FaceColor',[0.0745 0.6235 1.0000]); % Change background to white

% Download e.g. from here https://data.humdata.org/dataset/cod-ab-bmu
coast = shaperead('../../BermudaIsland/nw036zp7611.shp','UseGeoCoords',true);%,'RecordNumbers',2);


geoshow(coast, 'FaceColor', [.5 .5 .5])

geoshow(HogReefPOS(1), HogReefPOS(2), 'DisplayType', 'point','Marker','o', 'MarkerEdgeColor','k','MarkerFaceColor', HogReefCOL, 'MarkerSize', MkrSyze)
geoshow(CrescentReefPOS(1), CrescentReefPOS(2), 'DisplayType', 'point','Marker','o', 'MarkerEdgeColor','k','MarkerFaceColor', CrescentReefCOL, 'MarkerSize', MkrSyze)
geoshow(StGeorgePOS(1), StGeorgePOS(2), 'DisplayType', 'point','Marker','o', 'MarkerEdgeColor','k','MarkerFaceColor',StGeorgeCOL, 'MarkerSize', MkrSyze)
geoshow(BATSoceanPOS(1), BATSoceanPOS(2), 'DisplayType', 'point','Marker',BATSoceanSYM, 'MarkerEdgeColor','k','MarkerFaceColor',BATSoceanCOL, 'MarkerSize', MkrSyze)
legend(["OSTIA";"Hog Reef";"Crescent Reef";"St.George Harnour"; "BATS"]);

%% Plot Timeseries
subplot(2,1,1)
    ylimits = [10 40];
    for j = 1 : length(ObsDates)
        line([ObsDates(j) ObsDates(j)], ylimits, 'Color', [.5 .5 .5], 'LineWidth', 1, 'LineStyle', ':')
        hold on
    end
    plot(HogReef.Time, HogReef.SST_C, 'Color', HogReefCOL, 'LineWidth', 1.5)
    plot(Crescent.Time, Crescent.SST_C, 'Color', CrescentReefCOL, 'LineWidth', 1.5)
    plot(NOAA_BATS_SYNC.TimeVector, NOAA_BATS_SYNC.SST_StGeorge, 'Color', StGeorgeCOL, 'LineWidth', 2)
    scatter(BATS.Time, BATS.SST_C, 'Marker', BATSoceanSYM, 'MarkerEdgeColor', 'k', 'MarkerFaceColor',BATSoceanCOL, 'LineWidth', 1, 'SizeData', 92)
    % OSTIA data
    plot(MeteoDataSYNC.Time, MeteoDataSYNC.SST, 'Color', [0 0 0], 'LineWidth', 3)
    ylabel('SST (ºC)')
    xlim([datetime('01/05/2013','InputFormat','dd/MM/uuuu'); ...
          datetime('01/02/2014','InputFormat','dd/MM/uuuu')]);  
    ylim([18 32])
    hold off
    grid on
    box on
    
    
subplot(2,1,2)
    ylimits = [10 40];
    for j = 1 : length(ObsDates)
        line([ObsDates(j) ObsDates(j)], ylimits, 'Color', [.5 .5 .5], 'LineWidth', 1, 'LineStyle', ':')
        hold on
    end
    plot(HogReef.Time, HogReef.SAL_psu, 'Color', HogReefCOL, 'LineWidth', 1.5)
    plot(Crescent.Time, Crescent.SAL_psu, 'Color', CrescentReefCOL, 'LineWidth', 1.5)
    %plot(BermudaWS_DataSYNC.TimeVector, BermudaWS_DataSYNC.goodSST, 'Color', StGeorgeCOL)
    scatter(BATS.Time, BATS.SAL_psu, 'Marker', BATSoceanSYM, 'MarkerEdgeColor', 'k', 'MarkerFaceColor',BATSoceanCOL, 'LineWidth', 1, 'SizeData', 92)
    ylabel('Salinity (PSU)')
    xlim([datetime('01/05/2013','InputFormat','dd/MM/uuuu'); ...
        datetime('01/02/2014','InputFormat','dd/MM/uuuu')]);
    ylim([35.4 37.2])
    hold off
    grid on
    box on

% %% Difference between OSTIA SST and other sites
% % Resample and aggregate timeseries on daily basis
% TTOSTIA = table2timetable(MeteoDataSYNC);
% TTOSTIA = retime(TTOSTIA,'daily','mean');
% dateIndx1 = TTOSTIA.Time(1);
% dateIndx2 = TTOSTIA.Time(end);
% 
% OBSdays = table(ObsDates, ones(length(ObsDates),1));
% TTobs = table2timetable(OBSdays);
% TTobs = retime(TTobs,'daily','mean');
% % Add few days on top to match array size
% datesremaining = (TTOSTIA.Time(1):TTobs.ObsDates(1)-days(1))';
% faketable = table2timetable( ...
%             table(datesremaining, nan(length(datesremaining), 1), 'VariableNames',{'ObsDates';'Var2'})...
%             );
% TTobs = [faketable; TTobs];
% % Add few days on bottom to match array size
% datesremaining = (TTobs.ObsDates(end)+days(1):TTOSTIA.Time(end))';
% faketable = table2timetable( ...
%             table(datesremaining, nan(length(datesremaining), 1), 'VariableNames',{'ObsDates';'Var2'})...
%             );
% TTobs = [TTobs; faketable];
% 
% TTHog = table2timetable(HogReef);
% TTHog = retime(TTHog,'daily','mean');
% TTHog = TTHog(find(TTHog.Time == dateIndx1, 1):find(TTHog.Time == dateIndx2, 1), :); % trim edges
% 
% TTCrescent = table2timetable(Crescent);
% TTCrescent = retime(TTCrescent,'daily','mean');
% TTCrescent = TTCrescent(find(TTCrescent.Time == dateIndx1, 1):find(TTCrescent.Time == dateIndx2, 1), :); % trim edges
% 
% TTStGeorge = table2timetable(NOAA_BATS_SYNC);
% TTStGeorge = retime(TTStGeorge,'daily','mean');
% TTStGeorge = TTStGeorge(find(TTStGeorge.TimeVector == dateIndx1, 1):find(TTStGeorge.TimeVector == dateIndx2, 1), :); % trim edges
% 
% TTBATS = table2timetable(table(BATS.Time, BATS.SST_C, BATS.SAL_psu, 'VariableNames',{'Time';'SST_C';'SAL_psu'}));
% TTBATS = retime(TTBATS,'daily','mean');
% % Add few days to match array size
% datesremaining = (TTBATS.Time(end)+days(1):TTOSTIA.Time(end))';
% faketable = table2timetable( ...
%             table(datesremaining, nan(length(datesremaining), 1), nan(length(datesremaining), 1), 'VariableNames',{'Time';'SST_C';'SAL_psu'})...
%             );
% TTBATS = [TTBATS; faketable];
% TTBATS = TTBATS(find(TTBATS.Time == dateIndx1, 1):find(TTBATS.Time == dateIndx2, 1), :); % trim edges
% 
% %% Plot Timeseries
% goodIndexes = ~isnan(TTobs.Var2);
% % SST dataset - SST OSTIA
% %plot(1.5:5.5, 11:15, '.-', 'LineWidth',10, 'DisplayName',' 1.0', 'Color',h1a.Color)
% h1b = plot(TTHog.Time, TTHog.SST_C - TTOSTIA.SST, 'Color', HogReefCOL, 'LineWidth', 3); h1b.Color(4)=0.10;  % 75% opaque
% hold on
% scatter(TTHog.Time(goodIndexes), TTHog.SST_C(goodIndexes) - TTOSTIA.SST(goodIndexes), 'MarkerEdgeColor', 'k', 'MarkerFaceColor',HogReefCOL)
% 
% h2b = plot(TTCrescent.Time, TTCrescent.SST_C - TTOSTIA.SST, 'Color', CrescentReefCOL, 'LineWidth', 3); h2b.Color(4)=0.10;  % 75% opaque
% scatter(TTCrescent.Time(goodIndexes), TTCrescent.SST_C(goodIndexes) - TTOSTIA.SST(goodIndexes), 'MarkerEdgeColor', 'k', 'MarkerFaceColor',CrescentReefCOL)
% 
% h3b = plot(TTStGeorge.TimeVector, TTStGeorge.SST_StGeorge - TTOSTIA.SST, 'Color', StGeorgeCOL, 'LineWidth', 2); h3b.Color(4)=0.10;  % 75% opaque
% scatter(TTStGeorge.TimeVector(goodIndexes), TTStGeorge.SST_StGeorge(goodIndexes) - TTOSTIA.SST(goodIndexes), 'MarkerEdgeColor', 'k', 'MarkerFaceColor',StGeorgeCOL)
% 
% scatter(TTBATS.Time, TTBATS.SST_C - TTOSTIA.SST, 'Marker', BATSoceanSYM,  'MarkerEdgeColor', 'k', 'MarkerFaceColor',BATSoceanCOL)
% % OSTIA data
% %plot(MeteoDataSYNC.Time, MeteoDataSYNC.SST, 'Color', [0 0 0], 'LineWidth', 3)
% ylabel('SST (ºC)')
% xlim([datetime('01/05/2013','InputFormat','dd/MM/uuuu'); ...
%       datetime('01/02/2014','InputFormat','dd/MM/uuuu')]);  
% hold off
% grid on
% box on
% 
% %% Isotopic composition of equilibrium vapor
% EqStGeorged18O = d18_V(SeaWater_d18O, TTStGeorge.SST_StGeorge(goodIndexes)+273.15);
% EqHogReefd18O = d18_V(SeaWater_d18O, TTHog.SST_C(goodIndexes)+273.15);
% EqCrescentd18O = d18_V(SeaWater_d18O, TTCrescent.SST_C(goodIndexes)+273.15);
% EqOSTIAd18O = d18_V(SeaWater_d18O, TTOSTIA.SST(goodIndexes)+273.15);
% 
% EqStGeorgedD = d2_V(SeaWater_dD, TTStGeorge.SST_StGeorge(goodIndexes)+273.15);
% EqHogReefdD = d2_V(SeaWater_dD, TTHog.SST_C(goodIndexes)+273.15);
% EqCrescentdD = d2_V(SeaWater_dD, TTCrescent.SST_C(goodIndexes)+273.15);
% EqOSTIAdD = d2_V(SeaWater_dD, TTOSTIA.SST(goodIndexes)+273.15);
% 
% %% Estimate distribution of differences between eq. vap with OSTIA and
% x_values = -1:0.01:1;
% % other measurement stations - d18O
% % Hog Reef - OSTIA
% pd = fitdist(EqHogReefd18O-EqOSTIAd18O,'Kernel', 'Kernel','normal', 'Bandwidth', .1);
% PDF = pdf(pd,x_values);
% plot(x_values, PDF, 'Color', HogReefCOL, 'LineWidth', 1.5)
% hold on
% % Crescent Reef - OSTIA
% pd = fitdist(EqCrescentd18O-EqOSTIAd18O,'Kernel', 'Kernel','normal', 'Bandwidth', .1);
% PDF = pdf(pd,x_values);
% plot(x_values, PDF, 'Color', CrescentReefCOL, 'LineWidth', 1.5)
% % StGeorge Reef - OSTIA
% pd = fitdist(EqStGeorged18O-EqOSTIAd18O,'Kernel', 'Kernel','normal', 'Bandwidth', .1);
% PDF = pdf(pd,x_values);
% plot(x_values, PDF, 'Color', StGeorgeCOL, 'LineWidth', 1.5)
% % Deviation from OSTIA (0)
% line([0 0], [0 7], 'LineStyle', '--', 'Color', [0 0 0])
% hold off
% xlabel('\delta^{18}O_{eq}* - \delta^{18}O_{eq} OSTIA (‰)')
% ylabel('PDF')
% xlim([-.5 .5])
% ylim([0 5])
% xticks(-.5:.2:.5)
% box on
% grid on
% legend({'Hog Reef';'Crescent Reef';'St.George'})
% %% Estimate distribution of differences between eq. vap with OSTIA and
% x_values = -5:0.01:5;
% % other measurement stations - dD
% % Hog Reef - OSTIA
% pd = fitdist(EqHogReefdD-EqOSTIAdD,'Kernel', 'Kernel','normal', 'Bandwidth', 1);
% PDF = pdf(pd,x_values);
% plot(x_values, PDF, 'Color', HogReefCOL, 'LineWidth', 1.5)
% hold on
% % Crescent Reef - OSTIA
% pd = fitdist(EqCrescentdD-EqOSTIAdD,'Kernel', 'Kernel','normal', 'Bandwidth', 1);
% PDF = pdf(pd,x_values);
% plot(x_values, PDF, 'Color', CrescentReefCOL, 'LineWidth', 1.5)
% % StGeorge Reef - OSTIA
% pd = fitdist(EqStGeorgedD-EqOSTIAdD,'Kernel', 'Kernel','normal', 'Bandwidth', 1);
% PDF = pdf(pd,x_values);
% plot(x_values, PDF, 'Color', StGeorgeCOL, 'LineWidth', 1.5)
% % Deviation from OSTIA (0)
% line([0 0], [0 .7], 'LineStyle', '--', 'Color', [0 0 0])
% xlabel('\deltaD_{eq}* - \deltaD_{eq} OSTIA (‰)')
% ylabel('PDF')
% xlim([-5 5])
% ylim([0 .5])
% xticks(-5:1:5)
% box on
% grid on
% legend({'Hog Reef';'Crescent Reef';'St.George'})
% hold off
% 
% %% Calculate variaiblity of surface composition by mean of salinity data
% SlopeSd18O = 0.36; % ‰/PSU
% SlopeSdD =2.1; % ‰/PSU
% 
% % Salinity Difference between Crescent - Hog Reef
% SalDiff_Crescent_Hog = TTCrescent.SAL_psu - TTHog.SAL_psu;
% SalDiff_BATS_Hog = TTBATS.SAL_psu - TTHog.SAL_psu;
% 
% %% d18O
% x_values = -10:0.01:10;
% pd = fitdist(SalDiff_Crescent_Hog(goodIndexes)*SlopeSd18O,'Kernel', 'Kernel','normal', 'Bandwidth', .1);
% PDF = pdf(pd,x_values);
% hold on
% %histogram(SalDiff_Crescent_Hog(goodIndexes)*SlopeSd18O, 'Normalization', 'pdf')
% plot(x_values, PDF, 'Color', [0 0 0], 'LineWidth', 1.5)
% 
% % dD
% %yyaxis right
% pd = fitdist(SalDiff_Crescent_Hog(goodIndexes)*SlopeSdD,'Kernel', 'Kernel','normal', 'Bandwidth', .2);
% PDF = pdf(pd,x_values);
% plot(x_values, PDF, 'Color', [0 0 0], 'LineStyle','--', 'LineWidth', 1.5)
% %histogram(SalDiff_Crescent_Hog(goodIndexes)*SlopeSdD, 'Normalization', 'pdf')
% line([0 0], [0 4], 'LineStyle', '--', 'Color', [0 0 0])
% 
% xlabel('\delta_{Ocean} CR - \delta_{Ocean} HR (‰)')
% ylabel('PDF')
% xlim([-.75 .75])
% ylim([0 4])
% xticks(-.75:.25:.75)
% box on
% grid on
% hold off
% legend({'\delta^{18}O', '\deltaD'})

