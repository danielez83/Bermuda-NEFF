% Load simulation summary output
[num,txt,raw] = xlsread('../Excel/Threshold_Impact_on_Flux_and_K.xlsx');

col_names = txt(1,2:end);
case_name = txt(2:end,1);
num = num(:, 4:end);

F18O_col        = 6;
F18O_col_err    = 7;
k18O_col        = 12;
k18O_col_err    = 13;

F2H_col        = 8;
F2H_col_err    = 9;
k2H_col        = 14;
k2H_col_err    = 15;

%% Prepare subplots
f = figure;
f.Position = [440 70 560 727];

%% Prepare plot for FO-18
pos1 = [0.08 0.55 0.40 0.35];
subplot('Position', pos1)
%subplot(2,2,1)

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
%ax.XTickLabel{1} = 'All OFF';
%ax.XTickLabel{2} = 'Time';
%ax.XTickLabel{3} = 'WD';
%ax.XTickLabel{4} = 'iso';
%ax.XTickLabel{5} = 'w';
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
%ax.XTickLabel{1} = 'All OFF';
%ax.XTickLabel{2} = 'Time';
%ax.XTickLabel{3} = 'WD';
%ax.XTickLabel{4} = 'iso';
%ax.XTickLabel{5} = 'w';
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
