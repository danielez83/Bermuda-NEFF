%% Configuration and load
% Load simulation summary output
data_folder                 = '../GIT Data/'; % %% Where the data is located?
simulation_thres_data       = readtable(strcat(data_folder, 'Threshold_simulations.txt'));

col_names                   = string(simulation_thres_data.Properties.VariableNames);
case_name                   = string(simulation_thres_data.Threshold);
num                         = horzcat(simulation_thres_data.FGFD18OE, simulation_thres_data.err, ...
                                      simulation_thres_data.FGFDDE, simulation_thres_data.err_1, ...
                                      simulation_thres_data.KPFD18OE, simulation_thres_data.err_2, ...
                                      simulation_thres_data.KPFDDE, simulation_thres_data.err_3,...
                                      simulation_thres_data.TopD18O, ...
                                      simulation_thres_data.TopDD, ...
                                      simulation_thres_data.k18, simulation_thres_data.err_4, ...
                                      simulation_thres_data.k2, simulation_thres_data.err_5);
F18O_col        = 5;
F18O_col_err    = 6;
k18O_col        = 11;
k18O_col_err    = 12;

F2H_col        = 7;
F2H_col_err    = 8;
k2H_col        = 13;
k2H_col_err    = 14;

%% Prepare subplots
f = figure;
f.Position = [440 70 560 727];

%% Prepare plot for FO-18
pos1 = [0.08 0.55 0.40 0.35];


s = errorbar(1:length(case_name)-1, num(1:end-1,F18O_col), num(1:end-1,F18O_col_err),'o');
s.MarkerSize = 10;
s.MarkerEdgeColor = [0 0 0];
s.LineWidth = 1;
s.MarkerFaceColor = [1 1 1];
s.DisplayName = 'd18OF obs.';
s.Color = [0,0,0];

% Plot all thresholds on
v = [0 num(end,F18O_col)-num(end,F18O_col_err);...
    10 num(end,F18O_col)-num(end,F18O_col_err);...
    10 num(end,F18O_col)+num(end,F18O_col_err);...
    0 num(end,F18O_col)+num(end,F18O_col_err)];
f = [1 2 3 4];
patch('Faces',f,'Vertices',v,'FaceColor',[0.5,0.5,0.5], 'FaceAlpha', 0.3)

ax = gca;

ax.XLim = [0,6];
ax.YLim = [-8,0];
ax.XTick = [1:5];
ax.XTickLabel = [];
ax.YLabel.String = '\delta{18}O_{E} (‰)';
grid on

%% Prepare plot for FH-2
%subplot(2,2,3)
pos2 = [0.08 0.15 0.40 0.35];
subplot('Position', pos2)

s = errorbar(1:length(case_name)-1, num(1:end-1,F2H_col), num(1:end-1,F2H_col_err),'o');
s.MarkerSize = 10;
s.MarkerEdgeColor = [0 0 0];
s.LineWidth = 1;
s.MarkerFaceColor = [1 1 1];
s.DisplayName = 'd18OF obs.';
s.Color = [0,0,0];

% Plot all thresholds on
v = [0 num(end,F2H_col)-num(end,F2H_col_err);...
    10 num(end,F2H_col)-num(end,F2H_col_err);...
    10 num(end,F2H_col)+num(end,F2H_col_err);...
    0 num(end,F2H_col)+num(end,F2H_col_err)];
f = [1 2 3 4];
patch('Faces',f,'Vertices',v,'FaceColor',[0.5,0.5,0.5], 'FaceAlpha', 0.3)

ax = gca;

ax.XLim = [0,6];
ax.YLim = [-60,0];
ax.YLabel.String = '\deltaD_{E} (‰)';
ax.XTick = [1:5];
ax.XTickLabel{1} = 'All OFF';
ax.XTickLabel{2} = 'Time';
ax.XTickLabel{3} = 'WD';
ax.XTickLabel{4} = 'iso';
ax.XTickLabel{5} = 'w';
grid on

%% Prepare plot for kO-18
%subplot(2,2,2)
pos3 = [0.55 0.55 0.40 0.35];
subplot('Position', pos3)

s = errorbar(1:length(case_name)-1, num(1:end-1,k18O_col), num(1:end-1,k18O_col_err),'o');
s.MarkerSize = 10;
s.MarkerEdgeColor = [0 0 0];
s.LineWidth = 1;
s.MarkerFaceColor = [1 1 1];
s.DisplayName = 'k18.';
s.Color = [0,0,0];

% Plot all thresholds on
v = [0 num(end,k18O_col)-num(end,k18O_col_err);...
    10 num(end,k18O_col)-num(end,k18O_col_err);...
    10 num(end,k18O_col)+num(end,k18O_col_err);...
    0 num(end,k18O_col)+num(end,k18O_col_err)];
f = [1 2 3 4];
patch('Faces',f,'Vertices',v,'FaceColor',[0.5,0.5,0.5], 'FaceAlpha', 0.3)

ax = gca;

ax.XLim = [0,6];
ax.YLim = [4,8];
ax.YLabel.String = 'k_{18} (‰)';
ax.XTick = [1:5];
ax.XTickLabel = [];
grid on

%% Prepare plot for kH-2
%subplot(2,2,4)
pos4 = [0.55 0.15 0.40 0.35];
subplot('Position', pos4)

s = errorbar(1:length(case_name)-1, num(1:end-1,k2H_col), num(1:end-1,k2H_col_err),'o');
s.MarkerSize = 10;
s.MarkerEdgeColor = [0 0 0];
s.LineWidth = 1;
s.MarkerFaceColor = [1 1 1];
s.DisplayName = 'k18.';
s.Color = [0,0,0];

% Plot all thresholds on
v = [0 num(end,k2H_col)-num(end,k2H_col_err);...
    10 num(end,k2H_col)-num(end,k2H_col_err);...
    10 num(end,k2H_col)+num(end,k2H_col_err);...
    0 num(end,k2H_col)+num(end,k2H_col_err)];
f = [1 2 3 4];
patch('Faces',f,'Vertices',v,'FaceColor',[0.5,0.5,0.5], 'FaceAlpha', 0.3)

ax = gca;

ax.XLim = [0,6];
ax.YLim = [-2,20];
ax.YLabel.String = 'k_{2} (‰)';
ax.XTick = [1:5];
ax.XTickLabel{1} = 'All OFF';
ax.XTickLabel{2} = 'Time';
ax.XTickLabel{3} = 'WD';
ax.XTickLabel{4} = 'iso';
ax.XTickLabel{5} = 'w';
grid on
