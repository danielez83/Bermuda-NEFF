% Author:      Daniele Zannoni
% Name:        DisplayTimeSeries.m
% Description: Plot timeseries of several variables
% Date:        Last revision 07/12/2020

LoadAndFilter
load('traj_analysis/centroid_pos_WATERSIP_latweight.mat')
WATERSIPFLAG = 1; % Show results for watersip
centroid_WATERSIP = centroid;

%load('../Matlab data/WVsourcesCentroid.mat')

%% Graphical configuration
xlimits = [datetime('2013-06-15') datetime('2014-01-5') ];
n = 6;
%% Plot d18O
% pos1 = [0.06 0.82 0.88 1/n-0.02];
pos1 = [0.1 0.82 0.80 1/n-0.02];
subplot('Position',pos1)

ylimits = [-21 -7];
for j = 1 : length(ObsDates)
    line([ObsDates(j) ObsDates(j)], ylimits, 'Color', [.5 .5 .5], 'LineWidth', .5, 'LineStyle', ':')
    hold on
end
% plot(BottomInletSYNC.Date, BottomInletSYNC.d18O, 'Color', [0 0 0], 'LineWidth', 1)
plot(TopInletSYNC.Date, TopInletSYNC.d18O, 'Color', [0 0 0], 'LineWidth', 1)
% scatter(ObsDates, Bottomd18O, 24, [1 0 0], 'filled');%, '.', 'Color', [1 0 0], 'LineWidth', 2)
hold off
box on
grid on
ax = gca;
set(ax,'XTickLabel',[]);
ax.XLim = xlimits;
ax.YLim = ylimits;
ax.YAxis.FontSize = 11;
ax.YAxis.Label.String = '\delta^{18}O (‰)';

%% Plot d excess
% pos2 = [0.06 0.82-1/n+0.02, 0.88 1/n-0.02];
pos2 = [0.1 0.82-1/n+0.02, 0.80 1/n-0.02];
subplot('Position',pos2)

ylimits = [-5 40];
yyaxis right
for j = 1 : length(ObsDates)
    line([ObsDates(j) ObsDates(j)], ylimits, 'Color', [.5 .5 .5], 'LineWidth', .5, 'LineStyle', ':')
    hold on
end
% plot(BottomInletSYNC.Date, BottomInletSYNC.Dexess, 'Color', [0 0 0], 'LineWidth', 1)
plot(TopInletSYNC.Date, TopInletSYNC.Dexess, 'Color', [0 0 0], 'LineWidth', 1)
% scatter(ObsDates, Bottomd18O, 24, [1 0 0], 'filled');%, '.', 'Color', [1 0 0], 'LineWidth', 2)
hold off
box on
grid on
ax = gca;
set(ax,'XTickLabel',[]);
yyaxis left
ax = gca;
set(ax,'YTickLabel',[]);
set(ax,'XTickLabel',[]);
yyaxis right
ax = gca;
ax.XLim = xlimits;
ax.YLim = ylimits;
ax.YAxis(2).FontSize = 11;
ax.YAxis(2).Label.String = 'd-excess (‰)';
ax.YAxis(2).Label.Color = [0 0 0];
ax.YAxis(2).Color = [0 0 0];

%% Plot q
% pos3 = [0.06 0.6933-1/n, 0.88 1/n-0.02];
pos3 = [0.1 0.6933-1/n, 0.80 1/n-0.02];
subplot('Position',pos3)

ylimits = [5000 40000];
for j = 1 : length(ObsDates)
    line([ObsDates(j) ObsDates(j)], ylimits, 'Color', [.5 .5 .5], 'LineWidth', .5, 'LineStyle', ':')
    hold on
end
% plot(BottomInletSYNC.Date, BottomInletSYNC.H2O, 'Color', [0 0 0], 'LineWidth', 1)
plot(TopInletSYNC.Date, TopInletSYNC.H2O, 'Color', [0 0 0], 'LineWidth', 1)
% scatter(ObsDates, Bottomd18O, 24, [1 0 0], 'filled');%, '.', 'Color', [1 0 0], 'LineWidth', 2)
hold off
box on
grid on
ax = gca;
set(ax,'XTickLabel',[]);
ax.XLim = xlimits;
ax.YLim = ylimits;
ax.YAxis.FontSize = 11;
ax.YAxis.Label.String = 'H_{2}O (ppmv)';
ax.YAxis.Label.Color = [0 0 0];
ax.YAxis.Color = [0 0 0];

%% Plot WS
% pos4 = [0.06 0.5466-1/n, 0.88 1/n-0.02];
pos4 = [0.1 0.5466-1/n, 0.80 1/n-0.02];
subplot('Position',pos4)

ylimits = [0 15];
yyaxis left
for j = 1 : length(ObsDates)
    line([ObsDates(j) ObsDates(j)], [0 100], 'Color', [.5 .5 .5], 'LineWidth', .5, 'LineStyle', ':')
    hold on
end
plot(MeteoDataSYNC.Time, MeteoDataSYNC.WS, 'Color', [0 0 1], 'LineWidth', 1)
ax = gca;
set(ax,'XTickLabel',[]);
ax.XLim = xlimits;
ax.YLim = ylimits;
ax.YColor = [0 0 1];
ax.YAxis(1).FontSize = 11;
ax.YAxis(1).Label.String = 'WS (m s^{-1})';


yyaxis right
plot(MeteoDataSYNC.Time, MeteoDataSYNC.SST, 'Color', [1.0000 0.0745 0.6510], 'LineWidth', 1)
ax = gca;
ax.YAxis(2).FontSize = 11;
ax.YAxis(2).Label.String = 'SST (°C)';
ax.YAxis(2).Label.Color = [1.0000 0.0745 0.6510];
ax.YAxis(2).Color = [1.0000 0.0745 0.6510];

box on
grid on

%% Calculate centroid dates
% centroid_dates = [datetime(2013, 6, 1); ...
%                   datetime(2013, 7, 1); ...
%                   datetime(2013, 8, 1); ...
%                   datetime(2013, 9, 1); ...
%                   datetime(2013, 10, 1); ...  
%                   datetime(2013, 11, 1); ...  
%                   datetime(2013, 12, 1)];
centroid_dates = datetime(2013, centroid_months, 15);
%% Plot centroid latitude
% pos5 = [0.06 0.3999-1/n, 0.88 1/n-0.02];
pos5 = [0.1 0.3999-1/n, 0.80 1/n-0.02];
subplot('Position',pos5)

%yyaxis left  

for j = 1 : length(ObsDates)
    line([ObsDates(j) ObsDates(j)], [-90 90], 'Color', [.5 .5 .5], 'LineWidth', .5, 'LineStyle', ':')
    hold on
end

if WATERSIPFLAG == 1
    plot(centroid_dates, centroid_WATERSIP(:, 1), 'Color', [1 0 0], 'LineWidth', 1)
    ylimits = [27 37];
else
    plot(centroid_time, centroid(:, 1), 'Color', [1 0 0], 'LineWidth', 1)
    ylimits = [18 45];
end
% scatter(ObsDates, Bottomd18O, 24, [1 0 0], 'filled');%, '.', 'Color', [1 0 0], 'LineWidth', 2)
hold off
box on
ax = gca;
set(ax,'XTickLabel',[]);
ax.XLim = xlimits;
ax.YLim = ylimits;
ax.YAxis.FontSize = 11;
ax.YAxis.Label.String = 'Centr. Lat (d)';
ax.YAxis.Label.Color = [1 0 0];
ax.YAxis.Color = [1 0 0];

% yyaxis right
% ylimits = [-90 -40];
% %plot(centroid_time, centroid(:, 2), 'Color', [1 0 0], 'LineWidth', 1)
% % scatter(ObsDates, Bottomd18O, 24, [1 0 0], 'filled');%, '.', 'Color', [1 0 0], 'LineWidth', 2)
% hold off
% box on
% ax = gca;
% % set(ax,'YTickLabel',[]);
% yyaxis right
% ax = gca;
% ax.YLim = ylimits;
% ax.YAxis(2).FontSize = 13;
% ax.YAxis(2).Label.String = 'Centr. Lon (d)';
% ax.YAxis(2).Label.Color = [1 0 0];
% ax.YAxis(2).Color = [1 0 0];

%% Plot centroid longitude
% pos6 = [0.06 0.2532-1/n, 0.88 1/n-0.02];
pos6 = [0.1 0.2532-1/n, 0.80 1/n-0.02];

subplot('Position', pos6)


yyaxis right
for j = 1 : length(ObsDates)
    line([ObsDates(j) ObsDates(j)], [-90 90], 'Color', [.5 .5 .5], 'LineWidth', .5, 'LineStyle', ':')
    hold on
end

if WATERSIPFLAG == 1
    plot(centroid_dates, centroid_WATERSIP(:, 2), 'Color', [1 0 0], 'LineWidth', 1)
    ylimits = [-68 -58];
else
    plot(centroid_time, centroid(:, 2), 'Color', [1 0 0], 'LineWidth', 1)
    ylimits = [-100 -40];
end
% scatter(ObsDates, Bottomd18O, 24, [1 0 0], 'filled');%, '.', 'Color', [1 0 0], 'LineWidth', 2)
hold off
box on
yyaxis left

ax = gca;
set(ax,'YTickLabel',[]);
yyaxis right
ax = gca;
ax.XLim = xlimits;
ax.YLim = ylimits;
ax.YAxis(2).FontSize = 11;
ax.YAxis(2).Label.String = 'Centr. Lon (d)';
ax.YAxis(2).Label.Color = [1 0 0];
ax.YAxis(2).Color = [1 0 0];

ax.XAxis(1).FontSize = 13;


% %% Plot P-E
% load('../Matlab data/ERA5_Evap_Prec_TS.mat')
% PRECminusEVAP = timetable(datetimevector, (p*3600)+(e*3600));
% TT1 = retime(PRECminusEVAP,'daily','mean');
% 
% pos6 = [0.06 0.2532-1/n, 0.88 1/n-0.02];
% subplot('Position',pos6)
% 
% ylimits = [-1 4];
% % yyaxis right
% for j = 1 : length(ObsDates)
%     line([ObsDates(j) ObsDates(j)], ylimits, 'Color', [.5 .5 .5], 'LineWidth', .5, 'LineStyle', ':')
%     hold on
% end
% plot(TT1.datetimevector, TT1.Var1, 'Color', [0 0 0], 'LineWidth', 1)
% % scatter(ObsDates, Bottomd18O, 24, [1 0 0], 'filled');%, '.', 'Color', [1 0 0], 'LineWidth', 2)
% hold off
% box on
% grid on
% ax = gca;
% set(ax,'YTickLabel',[]);
% ax = gca;
% ax.XLim = xlimits;
% ax.YLim = ylimits;
% ax.YAxis.FontSize = 13;
% ax.YAxis.Label.String = 'P - E (mm day^{-1})';
% ax.YAxis.Label.Color = [0 0 0];
% ax.YAxis.Color = [0 0 0];
% 
% ax.XAxis.FontSize = 13;
%% Final settings
set(gcf, 'Color', [1 1 1])
%set(gcf, 'Position', [212 173 830 621])
set(gcf, 'Position', [465   176   602   621])