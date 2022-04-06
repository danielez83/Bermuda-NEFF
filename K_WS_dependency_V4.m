% Author:      Daniele Zannoni
% Name:        K18_WS_dependency.m
% Description: Plot k/18 as a function of WS
% Date:        Last revision 30/03/2023
%% Configuration
Show_no_tres_q_iso  = 0; % Display kinetic fractionation values vs WS when isotope and humidity thresholds are removed?
Show_k2             = 1; % Display results for k2
Show_k18_bestfit    = 1; % Display k18 best fit?
%% Data loade
tic
LoadAndFilter
if Show_no_tres_q_iso == 1
    load('../Matlab Data/WS_dep_no_q_iso_thresholds.mat')
    %WSbinnerVAls_no_thres = WSbinnerVAls_no_thres(2:end);
    %kinetic_18_16_WS_KP_no_thres = kinetic_18_16_WS_KP_no_thres(2:end);
    %kinetic_2_1_WS_KP_no_thres = kinetic_2_1_WS_KP_no_thres(2:end);
    %nVAls_no_thres = nVAls_no_thres(2:end);
end

%% Calculation of flux composition
% Preallocate memory for speed-up
d18OE_observed_FG   = nan(length(Topd18O), 1);      
dDE_observed_FG     = d18OE_observed_FG;
dexE_observed_FG    = d18OE_observed_FG;
d18OE_observed_KP   = d18OE_observed_FG;
dDE_observed_KP     = d18OE_observed_FG;
d18OE_CG            = nan(length(Topd18O), length(kVals));
dDE_CG              = d18OE_CG;
dexE_CG             = d18OE_CG;

for j = 1:length(Topd18O)
    % check if there are nans, skip calculation in case
    if ~isnan(sum([RH(j) AirT(j) Topd18O(j) SST(j) Bottomq(j) Topq(j) Bottomd18O(j) ...
                   TopdD(j) BottomdD(j)]))
        % Ocean Isotopic composition
        if UseConstantOceanDelta == 0
            SeaWater_d18O = SeaWater_d18O_COPY(j);
            SeaWater_dD = SeaWater_dD_COPY(j);
        end
        % ---------------------------------------------------------------------------------------------
        % FLUX GRADIENT CALCULATION -------------------------------------------------------------------
        % ---------------------------------------------------------------------------------------------       
        % Convert delta values to molar concentrations
        % Top inlet
        Top_x_H2O       = Topq(j)*1e-6;
        Top_x_18        = ((Topd18O(j)*1e-3 + 1)*R18SMOW) * Top_x_H2O;
        Top_x_2         = ((TopdD(j)*1e-3 + 1)*R2SMOW) * Top_x_H2O;
        % Bottom inlet
        Bottom_x_H2O    = Bottomq(j)*1e-6;
        Bottom_x_18     = ((Bottomd18O(j)*1e-3 + 1)*R18SMOW) * Bottom_x_H2O;
        Bottom_x_2      = ((BottomdD(j)*1e-3 + 1)*R2SMOW) * Bottom_x_H2O;
        % Calculate flux with eq.2
        RE_18 = (Top_x_18 - Bottom_x_18)/(Top_x_H2O - Bottom_x_H2O);
        RE_2  = (Top_x_2 - Bottom_x_2)/(Top_x_H2O - Bottom_x_H2O);
        % Convert to delta value
        d18OE_observed_FG(j) = (RE_18/R18SMOW - 1)*1e3;      
        dDE_observed_FG(j) = (RE_2/R2SMOW - 1)*1e3;
        dexE_observed_FG(j) = dDE_observed_FG(j) - 8 * d18OE_observed_FG(j); 
        % This is the same
        % mdl = fitlm([Top_x_H2O; Bottom_x_H2O], [Top_x_18; Bottom_x_18]); % Paragraph [21]
        % d18OE_obsered_COPY(j) = ((mdl.Coefficients.Estimate(2)/R18SMOW)-1)*1e3;      
        % ---------------------------------------------------------------------------------------------
        % Keeling Plot CALCULATION --------------------------------------------------------------------
        % --------------------------------------------------------------------------------------------- 
        mdl = fitlm([1/Topq(j); 1/Bottomq(j)], [Topd18O(j); Bottomd18O(j)]);
        d18OE_observed_KP(j) = mdl.Coefficients.Estimate(1);
        mdl = fitlm([1/Topq(j); 1/Bottomq(j)], [TopdD(j); BottomdD(j)]);
        dDE_observed_KP(j) = mdl.Coefficients.Estimate(1);
        % ---------------------------------------------------------------------------------------------
        % Craig-Gordon CALCULATION --------------------------------------------------------------------
        % --------------------------------------------------------------------------------------------- 
        switch inletForCG
            case 'top'
            e_s = (exp(77.3450 + 0.0057 .* (AirT(j)) - (7235./(AirT(j))))./(AirT(j)).^8.2);
            e_a = Topq(j)*1e-6*(SLP(j)*100);
            currRH = 100*e_a/e_s;
            d18OE_CG(j,:) = CG_dE_18_MJ79(currRH, AirT(j), Topd18O(j), SST(j), SeaWater_d18O, kVals);
            dDE_CG(j,:) = CG_dE_2_MJ79(currRH, AirT(j), TopdD(j), SST(j), SeaWater_dD, kVals);
            dexE_CG(j,:) = dDE_CG(j,:) - 8* d18OE_CG(j,:);
            case 'bottom'
            e_s = (exp(77.3450 + 0.0057 .* (AirT(j)) - (7235./(AirT(j))))./(AirT(j)).^8.2);
            e_a = Bottomq(j)*1e-6*(SLP(j)*100);
            currRH = 100*e_a/e_s;
            d18OE_CG(j,:) = CG_dE_18_MJ79(currRH, AirT(j), Bottomd18O(j), SST(j), SeaWater_d18O, kVals);
            dDE_CG(j,:) = CG_dE_2_MJ79(currRH, AirT(j), BottomdD(j), SST(j), SeaWater_dD, kVals);
            dexE_CG(j,:) = dDE_CG(j,:) - 8* d18OE_CG(j,:);
        end
                
    end
end
%% --------------------------------------------------------------------------------------------
% Find best k-value using FG ------------------------------------------------------------------
% ---------------------------------------------------------------------------------------------
% Step 1
DF18 = (d18OE_CG - repmat(d18OE_observed_FG, 1, size(d18OE_CG, 2))).^2;
DF2 = (dDE_CG - repmat(dDE_observed_FG, 1, size(dDE_CG, 2))).^2;
% Step 2
MIN18 = min(DF18, [], 2);
MIN2 = min(DF2, [], 2);
% Step 3
BM18 = DF18 == repmat(MIN18, 1, size(DF18, 2));
BM2 = DF2 == repmat(MIN2, 1, size(DF2, 2));
% Step 4
BV18 = sum(d18OE_CG.*BM18, 2);
BV2 = sum(dDE_CG.*BM2, 2);
% Step 5
BK18_FG = sum(repmat(kVals', size(BM18, 1), 1).*BM18, 2);
BK2_FG = sum(repmat(kVals', size(BM2, 1), 1).*BM2, 2);
%% --------------------------------------------------------------------------------------------
% Find best k-value using KP ------------------------------------------------------------------
% ---------------------------------------------------------------------------------------------
% Step 1
DF18 = (d18OE_CG - repmat(d18OE_observed_KP, 1, size(d18OE_CG, 2))).^2;
DF2 = (dDE_CG - repmat(dDE_observed_KP, 1, size(dDE_CG, 2))).^2;
% Step 2
MIN18 = min(DF18, [], 2);
MIN2 = min(DF2, [], 2);
% Step 3
BM18 = DF18 == repmat(MIN18, 1, size(DF18, 2));
BM2 = DF2 == repmat(MIN2, 1, size(DF2, 2));
% Step 4
BV18 = sum(d18OE_CG.*BM18, 2);
BV2 = sum(dDE_CG.*BM2, 2);
% Step 5
BK18_KP = sum(repmat(kVals', size(BM18, 1), 1).*BM18, 2);
BK2_KP = sum(repmat(kVals', size(BM2, 1), 1).*BM2, 2);

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

%%
BK18_KP_copy = BK18_KP; 
BK2_KP_copy = BK2_KP;
BK18_KP_copy(BK18_KP==max(kVals) | BK18_KP==min(kVals)) = nan;
BK2_KP_copy(BK2_KP==max(kVals) | BK2_KP==min(kVals)) = nan;
WSbinner = 0:.5:12;
% WSbinner = 0:2:12;
WSbinnerVAls = WSbinner+.25;
% WSbinnerVAls = WSbinner+1;
% WSbinner = 0:.2:12;
% WSbinnerVAls = WSbinner+.1;
WSbinnerVAls(end) = [];
kinetic_18_16_WS_KP = zeros(1, length(WSbinner)-1);
kinetic_2_1_WS_KP = zeros(1, length(WSbinner)-1);
nVAls = WSbinnerVAls;
for j = 1 : length(WSbinner) - 1
    indexesOI =  WS>=WSbinner(j) & WS<WSbinner(j+1);
    nVAls(j) = sum(indexesOI);
    kinetic_18_16_WS_KP(j) = nanmean(BK18_KP_copy(indexesOI));
    SE_kinetic_18_16_WS_KP(j) = nanstd(BK18_KP_copy(indexesOI))/sqrt(sum(indexesOI));  
    kinetic_2_1_WS_KP(j) = nanmean(BK2_KP_copy(indexesOI));
    SE_kinetic_2_1_WS_KP(j) = nanstd(BK2_KP_copy(indexesOI))/sqrt(sum(indexesOI));
    PBLasd(j) =  nanmean(PBLH(indexesOI));
    SSTasd(j) =  nanmean(SST(indexesOI));
    RHSSTasd(j) =  nanmean(RHSSTtop(indexesOI));
end
%% Plot estimated k-2
if Show_k2 == 1
    clf

    hold on
    s = errorbar(WSbinnerVAls,kinetic_2_1_WS_KP, SE_kinetic_2_1_WS_KP,'o');
    s.MarkerSize = 14;
    s.MarkerEdgeColor = [0 0 0];
    s.LineWidth = 2;
    s.MarkerFaceColor = [0.0745 0.6235 1.0000];
    hold off


    ax = gca;
    ax.XLabel.FontSize = 18;
    ax.XAxis.FontSize = 16;
    % ax.XLabel.FontWeight = 'bold';
    ax.XLabel.String = 'Wind Speed @10m (m s^{-1})';

    ax.YLabel.FontSize = 18;
    ax.YAxis.FontSize = 16;
    %ax.YLabel.FontWeight = 'bold';
    ax.YLabel.String = 'k_{2}  (‰)';

    xlim([0 15])
    ylim([-30 30])
    box on
    grid on

    s.MarkerEdgeColor = [0 0 0];
    s.LineWidth = 1;
    hold on
    % Plot linear model
    fitMDL = fitlm(WSbinnerVAls(1:21), kinetic_2_1_WS_KP(1:21));  % range 0 - 10 m/s
    plot(WSbinnerVAls, feval(fitMDL, WSbinnerVAls), 'LineWidth', 1.5, 'LineStyle', '-', 'Color', [0 0 0])
    [rho, pval] = corr(WSbinnerVAls(2:21)', kinetic_2_1_WS_KP(2:21)');
    fprintf('In wind speed range 0.5 10 m/s, correlation betwen k2 and WS is %.5f, p-value is %.5f\n', rho, pval)

    % Only smooth regime
    plot(SimWS_LOWER, k_smooth_LOWER*.88, 'LineWidth', 2, 'LineStyle', '--', 'Color', [0.4667    0.6745    0.1882]);
    % plot(SimWS_UPPER, k_smooth_UPPER*.88, 'LineWidth', 2, 'LineStyle', '--', 'Color', [0 0 0])

    % Only rough regime
    %plot(SimWS_LOWER, k_rough_LOWER*.88, 'LineWidth', 2, 'LineStyle', '--', 'Color', [0 0 0])
    plot(SimWS_UPPER, k_rough_UPPER*.88, 'LineWidth', 2, 'LineStyle', '--', 'Color', [0.3020    0.7451    0.9333]);

    % Plot linear model
    % plot(xfitMDL, yfitMDL)
    hold off


    ax = gca;
    %ax.XLabel.FontSize = 14;
    % ax.XLabel.FontWeight = 'bold';
    ax.XLabel.String = 'Wind Speed @10m (m s^{-1})';

    %ax.YLabel.FontSize = 14;
    %ax.YLabel.FontWeight = 'bold';
    ax.YLabel.String = 'k_{2} (‰)';

    xlim([0 15])

    ax.XTick = 0:3:15;

    ylim([-10 20])
    box on
    grid on
end


%% Store different values
if exist('thresholds', 'var')
    thresholds      = vertcat(thresholds, [q_diff_threshold, d18O_diff_threshol, dD_diff_threshol]);
    k18_data        = vertcat(k18_data, kinetic_18_16_WS_KP);
    k18_data_errors = vertcat(k18_data_errors, SE_kinetic_18_16_WS_KP);
else
    WS_main_data    = WSbinnerVAls;
    thresholds      = [q_diff_threshold, d18O_diff_threshol, dD_diff_threshol];
    k18_data        = kinetic_18_16_WS_KP;
    k18_data_errors = SE_kinetic_18_16_WS_KP;
end

%% Prepare subplots
f = figure;
f.Position = [440 70 560 727];
%% Plot with different thtresholds
%Plot estimated k-18
clf
%subplot(3,1,1)
pos1 = [0.1 0.5 0.80 0.45];
subplot('Position',pos1)

if Show_no_tres_q_iso == 1
    s1 = scatter(WSbinnerVAls_no_thres, kinetic_18_16_WS_KP_no_thres, 150, 'x');
    s1.MarkerEdgeColor = [.7, .7, .7];
    s1.LineWidth = 2;
    s1.DisplayName = 'Binned k obs. (no q-iso thres.)';
    % Color 20-24 in red...few datapoints
    hold on 
    s2 = scatter(WSbinnerVAls_no_thres(21:end), kinetic_18_16_WS_KP_no_thres(21:end), 150, 'x');
    s2.MarkerEdgeColor = [.7, .7, .7];
    s2.LineWidth = 2;
    s2.Annotation.LegendInformation.IconDisplayStyle = 'off';
end

hold on
Mysymbols = ["o"; "s"; "d"; "^"; "p"; "h"];
Mycolors = [0 0 0; 0 0 1; 0 1 0; 1 0 0; .5 .5 .5];
for j = 1:size(k18_data, 1)
    s = errorbar(WSbinnerVAls,k18_data(j,:),k18_data_errors(j, :), Mysymbols(j),'Color', 'k');
    s.MarkerSize = 10;
    s.MarkerEdgeColor = Mycolors(j, :);
    s.DisplayName = 'Binned k obs.';
end
for j = 21:24 % Color 20-24 in red...few datapoints
    hold on
    s3 = errorbar(WSbinnerVAls(j),k18_data(j), k18_data_errors(j),'o','Color', 'r');
    s3.MarkerSize = 10;
    s3.Annotation.LegendInformation.IconDisplayStyle = 'off';
end

hold off

s.MarkerEdgeColor = [0 0 0];
s.LineWidth = 1;
hold on
% Plot linear model
fitMDL = fitlm(WSbinnerVAls(2:21), kinetic_18_16_WS_KP(2:21));  % range 0 - 10 m/s
slope       = [fitMDL.Coefficients.Estimate(2), fitMDL.Coefficients.SE(2)];
intercept   = [fitMDL.Coefficients.Estimate(1), fitMDL.Coefficients.SE(1)];
Rsq         = fitMDL.Rsquared.Ordinary;
fprintf('Slope of fit:          %.2f ± %.2f\n', slope(1), slope(2))
fprintf('Intercept of fit:      %.2f ± %.2f\n', intercept(1), intercept(2))
fprintf('Rsquared of fit:       %.2f\n', Rsq) 
if Show_k18_bestfit == 1
    fit_plot = plot(WSbinnerVAls, feval(fitMDL, WSbinnerVAls), 'LineWidth', 1.5, 'LineStyle', '-', 'Color', [0 0 0]);
    fit_plot.DisplayName = 'Linear fit (1 - 10 m s^{-1})';
end

% Only smooth regime
smooth_plot = plot(SimWS_LOWER, k_smooth_LOWER, 'LineWidth', 2, 'LineStyle', '--', 'Color', [0.4667    0.6745    0.1882]);
%plot(SimWS_UPPER, k_smooth_UPPER, 'LineWidth', 2, 'LineStyle', '--', 'Color', [0 0 0])
smooth_plot.DisplayName = 'Smooth reg. param.';

% Only rough regime
%plot(SimWS_LOWER, k_rough_LOWER, 'LineWidth', 2, 'LineStyle', '--', 'Color', [0 0 0])
rough_plot = plot(SimWS_UPPER, k_rough_UPPER, 'LineWidth', 2, 'LineStyle', '--', 'Color', [0.3020    0.7451    0.9333]);
rough_plot.DisplayName = 'Rough reg. param.';

hold off

ax = gca;
ax.XLabel.FontSize = 12;
ax.XAxis.FontSize = 12;
ax.XLabel.FontWeight = 'bold';
ax.XAxis.FontWeight = 'bold';
%ax.XLabel.String = 'Wind Speed  @10m (m s^{-1})';
ax.XTick = 0:3:15;
ax.XTickLabel = [];

ax.YLabel.FontSize = 12;
ax.YAxis.FontSize = 12;
ax.YLabel.FontWeight = 'bold';
ax.YAxis.FontWeight = 'bold';
ax.YLabel.String = 'k_{18} (‰)';

xlim([0 12.5])
ylim([0 10])
box on
grid on

legend

%% Plot fractionation factor ratio
%subplot(3,1,2)
pos2 = [0.1 0.3 0.80 0.18];
subplot('Position',pos2)
if Show_no_tres_q_iso == 1
    plt = plot(WSbinnerVAls_no_thres(1:21), kinetic_2_1_WS_KP_no_thres(1:21)./kinetic_18_16_WS_KP_no_thres(1:21),...
        'LineWidth', 2, 'Color', [0.8 0.8 0.8]);
    plt.DisplayName = 'Binned k obs. (no q-iso thres.)';
end
hold on
plt2 = plot(WSbinnerVAls(1:21), kinetic_2_1_WS_KP(1:21)./kinetic_18_16_WS_KP(1:21),...
        'LineWidth', 2, 'Color', [0 0 0]);
plt2.DisplayName = 'Binned k obs.';

plt3 = plot(WSbinnerVAls(21:end), kinetic_2_1_WS_KP(21:end)./kinetic_18_16_WS_KP(21:end),...
        'LineWidth', 2, 'LineStyle', '--', 'Color', [1 0 0]);
plt3.Annotation.LegendInformation.IconDisplayStyle = 'off';

% Show text
annotation('textarrow', [.78, .78] , [.35, 0.30], ...
    'String', sprintf('%.2f', kinetic_2_1_WS_KP(22)./kinetic_18_16_WS_KP(22)),...
    'Color', 'r');

ax = gca;
ax.XLabel.FontSize = 12;
ax.XAxis.FontSize = 12;
ax.XLabel.FontWeight = 'bold';
ax.XAxis.FontWeight = 'bold';
%ax.XLabel.String = 'Wind Speed  @10m (m s^{-1})';
ax.XTick = 0:3:15;
ax.XTickLabel = [];

ax.YLabel.FontSize = 12;
ax.YAxis.FontSize = 12;
ax.YLabel.FontWeight = 'bold';
ax.YAxis.FontWeight = 'bold';
ax.YLabel.String = 'k_{2}/k_{18}';

xlim([0 12.5])
ylim([-1, 2.2])

lineplt = line([0 10], [nanmean(kinetic_2_1_WS_KP(1:21)./kinetic_18_16_WS_KP(1:21)) nanmean(kinetic_2_1_WS_KP(1:21)./kinetic_18_16_WS_KP(1:21))], ...
        'LineWidth', 2, 'LineStyle', '--', 'Color', [0 0 0]);
lineplt.Annotation.LegendInformation.IconDisplayStyle = 'off';
box on
grid on
legend('Location', 'southwest')
legend

mean_ratio = nanmean(kinetic_2_1_WS_KP(1:21)./kinetic_18_16_WS_KP(1:21));
std_error_ratio = nanstd(kinetic_2_1_WS_KP(1:21)./kinetic_18_16_WS_KP(1:21))/sqrt(length(kinetic_2_1_WS_KP(1:21)));
fprintf('k2/k18 ratio --> Mean = %.1f, Std. error = %.1f\n',mean_ratio, std_error_ratio);

if Show_no_tres_q_iso == 1
    mean_ratio = nanmean(kinetic_2_1_WS_KP_no_thres(1:21)./kinetic_18_16_WS_KP_no_thres(1:21));
    std_error_ratio = nanstd(kinetic_2_1_WS_KP_no_thres(1:21)./kinetic_18_16_WS_KP_no_thres(1:21))/sqrt(length(kinetic_2_1_WS_KP_no_thres(1:21)));
    fprintf('k2/k18 ratio --> Mean = %.1f, Std. error = %.1f (NO q or iso threshold)\n',mean_ratio, std_error_ratio);
end

%% Plot number of observations
%subplot(3,1,3)
pos3 = [0.1 0.1 0.80 0.18];
subplot('Position',pos3)

if Show_no_tres_q_iso == 1
    h2 = bar(WSbinnerVAls_no_thres, nVAls_no_thres');
    h2.FaceColor = [0.8 0.8 0.8];
    h2.FaceAlpha = .3;
    h2.EdgeColor = [0, 0, 0];
    h2.DisplayName = 'Binned k obs. (no q-iso thres.)';
    hold on
    h3 = bar(WSbinnerVAls_no_thres(22:end), nVAls_no_thres(22:end)');
    h3.FaceColor = [1 0 0];
    %h3.FaceAlpha = .3;
    h3.EdgeColor = [1, 0, 0];
    h3.Annotation.LegendInformation.IconDisplayStyle = 'off';
end

h = bar(WSbinnerVAls, nVAls');
h.FaceColor = [0 0 0];
h.EdgeColor = [0, 0, 0];
h.DisplayName = 'Binned k obs';
hold on
h4 = bar(WSbinnerVAls(21:24), nVAls(21:24)');
h4.FaceColor = [1 0 0];
%h3.FaceAlpha = .3;
h4.EdgeColor = [1, 0, 0];
h4.Annotation.LegendInformation.IconDisplayStyle = 'off';


grid('on')
xlim([0 12.5])
%ylim([0,150]);
ax = gca;
ax.XLabel.FontSize = 12;
ax.XAxis.FontSize = 12;
ax.XLabel.FontWeight = 'bold';
ax.XAxis.FontWeight = 'bold';
ax.XLabel.String = 'Wind Speed  @10m (m s^{-1})';
ax.XTick = 0:3:15;

ax.YLabel.FontSize = 12;
ax.YAxis.FontSize = 12;
ax.YLabel.FontWeight = 'bold';
ax.YAxis.FontWeight = 'bold';
ax.YLabel.String = 'Counts';

legend



%% Utils
k_smooth_LOWER(find(SimWS_LOWER==6))
k_rough_UPPER(find(SimWS_UPPER==6.1))

