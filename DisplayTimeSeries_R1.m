% Author:      Daniele Zannoni
% Name:        DisplayTimeSeries.m
% Description: Plot timeseries of several variables
% Date:        Last revision 07/12/2020

LoadAndFilter
%load('traj_analysis/centroid_pos_WATERSIP_latweight.mat')
%WATERSIPFLAG = 1; % Show results for watersip
%centroid_WATERSIP = centroid;

%load('../Matlab data/WVsourcesCentroid.mat')

%% Graphical configuration
xlimits = [datetime('2013-06-15') datetime('2014-01-5') ];
%n = 6; % INITIAL SUBMISSION
n = 5; % R1
%% Plot d18O
%pos1 = [0.1 0.82 0.80 1/n-0.02]; % INITIAL SUBMISSION
pos1 = [0.1 0.82 0.80 1/n-0.03]; % R1
subplot('Position',pos1)

ylimits = [-21 -7];
% plot(BottomInletSYNC.Date, BottomInletSYNC.d18O, 'Color', [0 0 0], 'LineWidth', 1)
plot(TopInletSYNC.Date, TopInletSYNC.d18O, 'Color', [0 0 0], 'LineWidth', 1)
hold on
% Assign NaN to observation not used
d18O_masked = BottomInletSYNC.d18O;
d18O_masked(~GoodIndexes) = nan;
% Plot data used for flux calculation
plot(TopInletSYNC.Date, d18O_masked, 'Color', [1 0 0], 'LineWidth', 2)
hold off

box on
grid on
grid minor
ax = gca;
set(ax,'XTickLabel',[]);
ax.XLim = xlimits;
ax.YLim = ylimits;
ax.YAxis.FontSize = 11;
ax.YAxis.Label.String = '\delta^{18}O (‰)';

%% Plot d excess
% pos2 = [0.1 0.82-1/n+0.02, 0.80 1/n-0.02]; % INITIAL SUBMISSION
pos2 = [0.1 0.83-1/n, 0.80 1/n-0.03]; %R1
subplot('Position',pos2)

ylimits = [0 35];
yyaxis right

plot(TopInletSYNC.Date, TopInletSYNC.Dexess, 'Color', [0 0 0], 'LineWidth', 1)
hold on
% Assign NaN to observation not used
d_masked = TopInletSYNC.Dexess;
d_masked(~GoodIndexes) = nan;
% Plot data used for flux calculation
plot(TopInletSYNC.Date, d_masked, 'Color', [1 0 0], 'LineWidth', 2)
hold off

box on
ax = gca;
set(ax,'XTickLabel',[]);
ytl = get(gca, 'YTick'); 
yyaxis left
ylim(ylimits)
grid on
grid minor
ax = gca;
set(ax,'YTickLabel',[]);
set(ax,'XTickLabel',[]);
yyaxis right
ytr = get(gca, 'YTick');
set(gca, 'YTick',ytr, 'YTickLabel',ytr)  
ax = gca;
ax.XLim = xlimits;
ax.YLim = ylimits;
ax.YAxis(2).FontSize = 11;
ax.YAxis(2).Label.String = 'd-excess (‰)';
ax.YAxis(2).Label.Color = [0 0 0];
ax.YAxis(2).Color = [0 0 0];
ax.YAxis(1).Color = [0 0 0];


%% Plot q
pos3 = [0.1 0.635-1/n, 0.80 1/n-0.03];
subplot('Position',pos3)

ylimits = [5000 40000];

% plot(BottomInletSYNC.Date, BottomInletSYNC.H2O, 'Color', [0 0 0], 'LineWidth', 1)
plot(TopInletSYNC.Date, TopInletSYNC.H2O, 'Color', [0 0 0], 'LineWidth', 1)
hold on
% Assign NaN to observation not used
q_masked = TopInletSYNC.H2O;
q_masked(~GoodIndexes) = nan;
% Plot data used for flux calculation
plot(TopInletSYNC.Date, q_masked, 'Color', [1 0 0], 'LineWidth', 2)

hold off
box on
ax = gca;
set(ax,'XTickLabel',[]);
ax.XLim = xlimits;
ax.YLim = ylimits;
ax.YAxis.FontSize = 11;
ax.YAxis.Label.String = 'H_{2}O (ppmv)';
ax.YAxis.Label.Color = [0 0 0];
ax.YAxis.Color = [0 0 0];
grid on
grid minor

%% Plot WS
pos4 = [0.1 0.44-1/n, 0.80 1/n-0.03];
subplot('Position',pos4)

ylimits = [0 15];
yyaxis left
plot(MeteoDataSYNC.Time, MeteoDataSYNC.WS, 'Color', [0 0 0], 'LineWidth', 1)
hold on
% Assign NaN to observation not used
WS_masked = MeteoDataSYNC.WS;
WS_masked(~GoodIndexes) = nan;
% Plot data used for flux calculation
plot(TopInletSYNC.Date, WS_masked, 'Color', [1 0 0], 'LineWidth', 2)
hold off

ax = gca;
set(ax,'XTickLabel',[]);
ax.XLim = xlimits;
ax.YLim = ylimits;
ax.YColor = [0 0 0];
ax.YAxis(1).FontSize = 11;
ax.YAxis(1).Label.String = 'WS (m s^{-1})';


yyaxis right
%plot(MeteoDataSYNC.Time, MeteoDataSYNC.SST, 'Color', [1.0000 0.0745 0.6510], 'LineWidth', 1)
plot(MeteoDataSYNC.Time, MeteoDataSYNC.SST, 'Color', [0 0 1], 'LineWidth', 1)
%plot(MeteoDataSYNC.Time, MeteoDataSYNC.SST, 'Color', [0, 0, 0], 'LineWidth', 1, 'LineStyle', '-.')
hold on
% Assign NaN to observation not used
SST_masked = MeteoDataSYNC.SST;
SST_masked(~GoodIndexes) = nan;
% Plot data used for flux calculation
scatter(TopInletSYNC.Date, SST_masked, 10, 'MarkerFaceColor', 'r')
hold off
ax = gca;
ax.YAxis(2).FontSize = 11;
ax.YAxis(2).Label.String = 'SST (°C)';
%ax.YAxis(2).Label.Color = [1.0000 0.0745 0.6510];
%ax.YAxis(2).Color = [1.0000 0.0745 0.6510];
ax.YAxis(2).Label.Color = [0 0 1];
ax.YAxis(2).Color = [0 0 1];

box on
grid on
grid minor

% ------------ R1
%% Plot hs
pos5 = [0.1 0.24-1/n, 0.80 1/n-0.03];
subplot('Position',pos5)
ylimits = [.3 1];

box on
%ax = gca;
%set(ax,'XTickLabel',[]);

plot(TopInletSYNC.Date, TopInletRHSST/100, 'Color', [0 0 0], 'LineWidth', 1)

hold on
% Assign NaN to observation not used
h_masked = TopInletRHSST/100;
h_masked(~GoodIndexes) = nan;
% Plot data used for flux calculation
plot(TopInletSYNC.Date, h_masked, 'Color', [1 0 0], 'LineWidth', 2)
hold off

grid on 
grid minor

box on
ylabel('\it{h_{s}} (-)', 'FontSize', 11);
ylim(ylimits)
xlim(xlimits)
ax = gca;
ax.YAxis(1).FontSize = 11;
%ax.XAxis(1).FontSize = 13;
%ax.XAxis(1).Limits = xlimits;

%% Final settings
set(gcf, 'Color', [1 1 1])
%set(gcf, 'Position', [212 173 830 621])
set(gcf, 'Position', [465   176   602   621])