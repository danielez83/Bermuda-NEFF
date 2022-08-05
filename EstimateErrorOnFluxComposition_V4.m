% Author:      Daniele Zannoni
% Name:        EstimateErrorOnFluxComposition_V3.m
% Description: 
% Date:        Last revision 29/03/2022
tic
%% Configuration
% Bootstrapping parameters
nIterations             = 1e2;     % Repetition for bootstrapping
rng('default');                    % Set default random  number generator to default seed (0), required to reproduce results
s                       = rng;

% Plot regression lines for each day?
PlotRegressions         = 0;

% Show Flux Gradient in addition to Keeling Plot results?
ShowFG                  = 0;

% Load and filter data
LoadAndFilter

%% Subset data day by day
% this way you will have more than two points to fit a line and you can
% calculate regression error

DayOfInterest               = (datetime(2013, 06, 20):days(1):datetime(2013, 12, 29))';

% Preallocate memory
% FGValues                    = NaN(length(Topq), 2); % Isotopic composition of flux in delta notation, d18O and dD
FGValues                    = NaN(length(DayOfInterest), 2); % Isotopic composition of flux in delta notation, d18O and dD
epsilonValues_FG            = FGValues;
epsilonValues_delta_FG      = FGValues;
sigmaET_R_FG                = FGValues;
sigmaET_delta_FG            = FGValues;

% KPValues                    = NaN(length(Topq), 2); % Isotopic composition of flux in delta notation, d18O and dD
KPValues                    = NaN(length(DayOfInterest), 2); % Isotopic composition of flux in delta notation, d18O and dD
epsilonValues_KP            = KPValues;
epsilonValues_delta_KP      = KPValues;
sigmaET_R_KP                = KPValues;
sigmaET_delta_KP            = KPValues;

% Observations subset
TopVals_Averages            = NaN(length(DayOfInterest), 3);
TopVals_STDs                = NaN(length(DayOfInterest), 3);
WSs                         = NaN(length(DayOfInterest), 3);
SSTs                        = NaN(length(DayOfInterest), 1);
SSTdifferences              = NaN(length(DayOfInterest), 1);
PBLHs                       = NaN(length(DayOfInterest), 1);
RHSSTs                      = NaN(length(DayOfInterest), 1);
TopVals_dexcess             = NaN(length(DayOfInterest), 1);
Ns                          = NaN(length(DayOfInterest), 1);


% kinetic fractionation values
MinkVal                 = -60;
MaxkVal                 = 60;
kVals = (MinkVal:0.01:MaxkVal)';

% R squared of regressions
RegresssionStatsTable = table(DayOfInterest, nan(length(DayOfInterest),1), nan(length(DayOfInterest),1),  nan(length(DayOfInterest),1),  nan(length(DayOfInterest),1),  nan(length(DayOfInterest),1), ...
                      'VariableNames', ["Date","NumObservations", "FGd18O","FGdD", "KPd18O", "KPdD"]);

%% Dataset is small, uncertainty cannot be approximated with eq. 25 in Good et 2012
%load('../Matlab data/SST_SYNC_OSTIA1x1deg.mat')
%local_SST_mean = nanmean(horzcat(SST_SYNC.SST_Crescent(GoodIndexes), SST_SYNC.SST_HOG(GoodIndexes), SST_SYNC.SST_StGeorge(GoodIndexes)), 2);

for j = 1 : length(DayOfInterest)
    GoodOBS = month(ObsDates) == month(DayOfInterest(j)) & day(ObsDates) == day(DayOfInterest(j));
    if sum(GoodOBS) > 2
        N = 2*sum(GoodOBS);
        % Subset data
        % mixratio, d18o, dD
        TopVals_deltas          = [Topq(GoodOBS) Topd18O(GoodOBS) TopdD(GoodOBS)];
        BottomVals_deltas       = [Bottomq(GoodOBS) Bottomd18O(GoodOBS) BottomdD(GoodOBS)];
        TopVals_Averages(j, :)  = mean(TopVals_deltas);
        TopVals_STDs(j, :)      = std(TopVals_deltas);
        % wind speed
        WSs(j,1)                = mean(WS(GoodOBS));
        % SSTs
        SSTs(j,1)                = mean(SST(GoodOBS));
        
        %SSTdifferences      
        %SSTdifferences(j,1)     = mean(SST(GoodOBS) - local_SST_mean(GoodOBS)-273.15);
        
        % BLH
        PBLHs(j,1)              = mean(PBLH(GoodOBS));
        % BLH
        RHSSTs(j,1)             = mean(RHSST(GoodOBS));
        % Number of obervations for fitting FG and KP
        Ns(j, 1)                = sum(GoodOBS);
        % Ocean Isotopic composition
        if UseConstantOceanDelta == 0
            SeaWater_d18O = mean(SeaWater_d18O_COPY(GoodOBS));
            SeaWater_dD = mean(SeaWater_dD_COPY(GoodOBS));
        end
        
        
        % Flux and uncertainty estimation following Good et al 2012
        % ---------------------------------------------------------------------------------------------
        % FLUX GRADIENT CALCULATION -------------------------------------------------------------------
        % ---------------------------------------------------------------------------------------------
        % Convert to molar ratio
        Top_X_H2O           = TopVals_deltas(:,1)*1e-6;
        Bottom_X_H2O        = BottomVals_deltas(:,1)*1e-6;
        % Oxygen 18 vap ratio, Deuterium vap ratio
        TopVals_X           = ((TopVals_deltas(:,2:3)*1e-3 + 1).*[repmat(R18SMOW, N/2, 1), repmat(R2SMOW, N/2, 1)])...
                              .* Top_X_H2O;
        BottomVals_X        = ((BottomVals_deltas(:,2:3)*1e-3 + 1).*[repmat(R18SMOW, N/2, 1), repmat(R2SMOW, N/2, 1)])...
                              .* Bottom_X_H2O;
        % ----------------------------------------
        % Calculate regression model for Oxygen 18
        % ----------------------------------------
        mdl                             = fitlm(vertcat(Top_X_H2O, Bottom_X_H2O), vertcat(TopVals_X(:,1), BottomVals_X(:,1))); % Paragraph [21]
        epsilonValues_FG(j, 1) = sqrt((1/(N-2))*sum((vertcat(TopVals_X(:,1), BottomVals_X(:,1))...
                                                  - repmat(mdl.Coefficients.Estimate(1), N, 1)...
                                                  - repmat(mdl.Coefficients.Estimate(2), N, 1) .* vertcat(Top_X_H2O, Bottom_X_H2O)...
                                                  ).^2 ...
                                                 )...
                                   );% See eq. 14 in Good et al 2012
        % Equations 18 and 19, see Paragraph [23]
        epsilonValues_delta_FG(j, 1)    = epsilonValues_FG(j, 1)/(R18SMOW*mean(vertcat(TopVals_X(:,1), BottomVals_X(:,1))));
        sigmaET_R_FG(j, 1)              = epsilonValues_FG(j, 1)*sqrt(N/(N*(N-1)*std(vertcat(TopVals_X(:,1), BottomVals_X(:,1)))^2)); % Uncertainty on flux (18O) composition as molar ratio
        % coeff_of_var                    = std(vertcat(Top_X_H2O, Bottom_X_H2O))/mean(vertcat(Top_X_H2O, Bottom_X_H2O));
        coeff_of_var                    = std(vertcat(TopVals_deltas(:,1), BottomVals_deltas(:,1)))/mean(vertcat(TopVals_deltas(:,1), BottomVals_deltas(:,1)));
        sigmaET_delta_FG(j, 1)          = epsilonValues_delta_FG(j, 1)*(1/(coeff_of_var*sqrt(N)))*sqrt(N/(N-1)); % Uncertainty on flux (18O) composition in delta notation
        FGValues(j, 1)                  = ((mdl.Coefficients.Estimate(2)/R18SMOW)-1)*1e3; % See eq.6 Good et al 2012, Flux (18O) isotopic composition in delta notation
        % Save regression R
        RegresssionStatsTable.FGd18O(j) = mdl.Rsquared.Adjusted;
        RegresssionStatsTable.NumObservations(j) = mdl.NumObservations;
        if PlotRegressions
            subplot(2,2,1)
            plot(mdl)
            title(sprintf('FG d18O_E for %s', DayOfInterest(j)))
            xlabel('H_{2}^{16}O')
            ylabel('H_{2}^{18}O')
            text(.2,.9,sprintf('H_{2}^{18}O/H_{2}^{16}O = %.4f', mdl.Coefficients.Estimate(2)),'Units','normalized')
            text(.2,.8,sprintf('d^{18}O_E = %.2f', FGValues(j, 1)),'Units','normalized')
            text(.2,.7,sprintf('R^2 = %.2f', mdl.Rsquared.Adjusted),'Units','normalized')
            legend off
        end
        % ----------------------------------------
        % Calculate regression model for Deuterium
         % ----------------------------------------
        mdl                             = fitlm(vertcat(Top_X_H2O,Bottom_X_H2O), vertcat(TopVals_X(:,2), BottomVals_X(:,2)));
        epsilonValues_FG(j, 2) = sqrt((1/(N-2))*sum((vertcat(TopVals_X(:,2), BottomVals_X(:,2))...
                                                  - repmat(mdl.Coefficients.Estimate(1), N, 1)...
                                                  - repmat(mdl.Coefficients.Estimate(2), N, 1) .* vertcat(Top_X_H2O, Bottom_X_H2O)...
                                                  ).^2 ...
                                                 )...
                                   );% See eq. 14 in Good et al 2012
        % Equations 18 and 19, see Paragraph [23]
        epsilonValues_delta_FG(j, 2)    = epsilonValues_FG(j, 2)/(R2SMOW*mean(vertcat(TopVals_X(:,2), BottomVals_X(:,2))));
        sigmaET_R_FG(j, 2)              = epsilonValues_FG(j, 2)*sqrt(N/(N*(N-1)*std(vertcat(TopVals_X(:,2), BottomVals_X(:,2)))^2)); % Uncertainty on flux (D) composition as molar ratio
        % coeff_of_var                    = std(vertcat(Top_X_H2O, Bottom_X_H2O))/mean(vertcat(Top_X_H2O, Bottom_X_H2O));
        coeff_of_var                    = std(vertcat(TopVals_deltas(:,1), BottomVals_deltas(:,1)))/mean(vertcat(TopVals_deltas(:,1), BottomVals_deltas(:,1)));
        sigmaET_delta_FG(j, 2)          = epsilonValues_delta_FG(j, 2)*(1/(coeff_of_var*sqrt(N)))*sqrt(N/(N-1)); % Uncertainty on flux (D) composition in delta notation
        FGValues(j, 2)                  = ((mdl.Coefficients.Estimate(2)/R2SMOW)-1)*1e3; % See eq.6 Good et al 2012, Flux (D) isotopic composition in delta notation
        % Save regression R
        RegresssionStatsTable.FGdD(j) = mdl.Rsquared.Adjusted;        
        if PlotRegressions
            subplot(2,2,2)
            plot(mdl)
            title(sprintf('FG dD_E for %s', DayOfInterest(j)))
            xlabel('H_{2}O')
            ylabel('HDO')
            text(.2,.9,sprintf('HD^{16}O/H_{2}^{16}O = %.4f', mdl.Coefficients.Estimate(2)),'Units','normalized')
            text(.2,.8,sprintf('dD_E = %.2f', FGValues(j, 2)),'Units','normalized')
            text(.2,.7,sprintf('R^2 = %.2f', mdl.Rsquared.Adjusted),'Units','normalized')
            legend off
        end
        % ---------------------------------------------------------------------------------------------
        % KEELING PLOT CALCULATION --------------------------------------------------------------------
        % ---------------------------------------------------------------------------------------------
        H2Oinv                          = vertcat(1./TopVals_deltas(:,1), 1./BottomVals_deltas(:,1));
        % ----------------------------------------
        % Calculate regression model for Oxygen 18
        % ----------------------------------------
        mdl                             = fitlm(H2Oinv, vertcat(TopVals_deltas(:,2), BottomVals_deltas(:,2)));
        epsilonValues_delta_KP(j, 1)    = sqrt((1/(N-2))*sum((vertcat(TopVals_deltas(:,2), BottomVals_deltas(:,2))...
                                                             - repmat(mdl.Coefficients.Estimate(1), N, 1)...
                                                             - repmat(mdl.Coefficients.Estimate(2), N, 1) .* H2Oinv...
                                                             ).^2 ...
                                                             )...
                                              );% See eq. 14 in Good et al 2012
        coeff_of_var                    = std(vertcat(TopVals_deltas(:,1), BottomVals_deltas(:,1)))/mean(vertcat(TopVals_deltas(:,1), BottomVals_deltas(:,1)));
        sigmaET_delta_KP(j, 1)          = epsilonValues_delta_KP(j, 1)*(1/(coeff_of_var*sqrt(N)))*sqrt(coeff_of_var^2 + (N/(N-1)));
        KPValues(j, 1)                  = mdl.Coefficients.Estimate(1);
        % Save regression R
        RegresssionStatsTable.KPd18O(j) = mdl.Rsquared.Adjusted;        
        if PlotRegressions
            subplot(2,2,3)
            plot(mdl)
            title(sprintf('KP d18O_E for %s', DayOfInterest(j)))
            xlabel('1/(H2O) (1/ppm)')
            ylabel('\delta^{18}O_E')
            text(.5,.9,sprintf('d^{18}O_E = %.2f', mdl.Coefficients.Estimate(1)),'Units','normalized')
            %text(.5,.8,sprintf('d^{18}O_E = %.2f', FGValues(j, 1)),'Units','normalized')
            text(.5,.8,sprintf('R^2 = %.2f', mdl.Rsquared.Adjusted),'Units','normalized')
            legend off
        end       
        % ----------------------------------------
        % Calculate regression model for Deuterium
         % ----------------------------------------
        mdl                             = fitlm(H2Oinv, vertcat(TopVals_deltas(:,3), BottomVals_deltas(:,3)));
        epsilonValues_delta_KP(j, 2)    = sqrt((1/(N-2))*sum((vertcat(TopVals_deltas(:,3), BottomVals_deltas(:,3))...
                                                             - repmat(mdl.Coefficients.Estimate(1), N, 1)...
                                                             - repmat(mdl.Coefficients.Estimate(2), N, 1) .* H2Oinv...
                                                             ).^2 ...
                                                             )...
                                              );% See eq. 14 in Good et al 2012
        coeff_of_var                    = std(vertcat(TopVals_deltas(:,1), BottomVals_deltas(:,1)))/mean(vertcat(TopVals_deltas(:,1), BottomVals_deltas(:,1)));
        sigmaET_delta_KP(j, 2)          = epsilonValues_delta_KP(j, 2)*(1/(coeff_of_var*sqrt(N)))*sqrt(coeff_of_var^2 + (N/(N-1)));
        KPValues(j, 2)                  = mdl.Coefficients.Estimate(1);
        % Save regression R
        RegresssionStatsTable.KPdD(j) = mdl.Rsquared.Adjusted;        
        if PlotRegressions
            fprintf('Iteration #%d, ', j)
            fprintf('%s (N = %d)\n', DayOfInterest(j), 2*Ns(j))
            subplot(2,2,4)
            plot(mdl)
            title(sprintf('KP dD_E for %s', DayOfInterest(j)))
            xlabel('1/(H2O) (1/ppm)')
            ylabel('\deltaD_E')
            text(.5,.9,sprintf('dD_E = %.2f', mdl.Coefficients.Estimate(1)),'Units','normalized')
            %text(.5,.8,sprintf('d^{18}O_E = %.2f', FGValues(j, 1)),'Units','normalized')
            text(.5,.8,sprintf('R^2 = %.2f', mdl.Rsquared.Adjusted),'Units','normalized')
            legend off
            pause
        end      
        % -----------------------------------------------------------------
        % Calculate flux with Craig-Gordon --> chose inlet in cofiguration
        % -----------------------------------------------------------------
        AirT_       = mean(AirT(GoodOBS));
        SLP_        = mean(SLP(GoodOBS));
        SST_        = mean(SST(GoodOBS));
        switch inletForCG 
            case 'top' %use top data
                Topd18O_    = mean(Topd18O(GoodOBS));
                TopdD_      = mean(TopdD(GoodOBS));
                Topq_       = mean(Topq(GoodOBS));                
                e_s = (exp(77.3450 + 0.0057 .* (AirT_) - (7235./(AirT_)))./(AirT_).^8.2);
                e_a = Topq_*1e-6*(SLP_*100);
                currRH = 100*e_a/e_s;
                d18OE_CG(j,:) = CG_dE_18_MJ79(currRH, AirT_, Topd18O_, SST_, SeaWater_d18O, kVals);
                dDE_CG(j,:) = CG_dE_2_MJ79(currRH, AirT_, TopdD_, SST_, SeaWater_dD, kVals);
                dexE_CG(j,:) = dDE_CG(j,:) - 8* d18OE_CG(j,:);
            case 'bottom'
                Bottomd18O_    = mean(Bottomd18O(GoodOBS));
                BottomdD_      = mean(BottomdD(GoodOBS));
                Bottomq_       = mean(Bottomq(GoodOBS));   
                e_s = (exp(77.3450 + 0.0057 .* (AirT_) - (7235./(AirT_)))./(AirT_).^8.2);
                e_a = Bottomq_*1e-6*(SLP_*100);
                currRH = 100*e_a/e_s;
                d18OE_CG(j,:) = CG_dE_18_MJ79(currRH, AirT_, Bottomd18O_, SST_, SeaWater_d18O, kVals);
                dDE_CG(j,:) = CG_dE_2_MJ79(currRH, AirT_, BottomdD_, SST_, SeaWater_dD, kVals);
                dexE_CG(j,:) = dDE_CG(j,:) - 8* d18OE_CG(j,:); 
            otherwise % just use top data as default selection
                Topd18O_    = mean(Topd18O(GoodOBS));
                TopdD_      = mean(TopdD(GoodOBS));
                Topq_       = mean(Topq(GoodOBS));                     
                e_s = (exp(77.3450 + 0.0057 .* (AirT_) - (7235./(AirT_)))./(AirT_).^8.2);
                e_a = Topq_*1e-6*(SLP_*100);
                currRH = 100*e_a/e_s;
                d18OE_CG(j,:) = CG_dE_18_MJ79(currRH, AirT_, Topd18O_, SST_, SeaWater_d18O, kVals);
                dDE_CG(j,:) = CG_dE_2_MJ79(currRH, AirT_, TopdD_, SST_, SeaWater_dD, kVals);
                dexE_CG(j,:) = dDE_CG(j,:) - 8* d18OE_CG(j,:);
        end
    end
end

%% Load flux isotopic composition
%load('../Matlab data/FluxComposition Data/BootstrappingData_Noone.mat', 'd18OE_CG_Top_sample', 'dDE_CG_Top_sample')
% d18OE_MixModel  = d18OE_CG_Top_sample;
d18OE_MixModel  = KPValues(:,1);
% dDE_MixModel   = dDE_CG_Top_sample;
dDE_MixModel  = KPValues(:,2);
%load('../Matlab data/FluxComposition Data/BootstrappingData_FLUX.mat', 'd18OE_CG_Top_sample', 'dDE_CG_Top_sample', 'Topd18O', 'TopdD');
% d18OE_FLUX      = d18OE_CG_Top_sample;
d18OE_FLUX      = FGValues(:,1);
% dDE_FLUX        = dDE_CG_Top_sample;
dDE_FLUX        = FGValues(:,2);

%% Plot configuration

display_histograms      = 0;
display_kernelPDF       = 1;
display_flux            = 1;
xlimits18               = [0 10];
xlimits2                = [0 100];

%% Plot error distributions for d18O
clf
% Flux boxplot
if display_flux
    box on
    % boxplot(horzcat(d18OE_MixModel, d18OE_FLUX, [Topd18O; nan(length(d18OE_MixModel) - length(Topd18O), 1)]),...
    boxplot(horzcat(d18OE_MixModel, d18OE_FLUX, TopVals_Averages(:,2)),...
        'Orientation', 'horizontal',...
        'labels', ["","", ""],...% 'labels', ["KP","FG", "d18Top"],...
        'Colors', [0 0 0]);
    xlabel('\delta^{18}O (‰)')
    ax = gca;
    ax.YAxis.FontSize = 13;
    ax.XAxis.FontSize = 13;
    set(gcf, 'Color', [1 1 1]);
    xlim([-25 5])
    % Boxplot properties
    a = get(get(gca,'children'),'children');   % Get the handles of all the objects
    t = get(a,'tag');   % List the names of all the objects 
    box1 = a(7);   % The 7th object is the first box
    outliers1 = a(1);
    box2 = a(8);   % The 8th object is the first box
    outliers2 = a(2);
    box3 = a(9);   % The 9th object is the first box
    outliers3 = a(3);
    % set(box1, 'Color', 'g');   % Set the color of the first box to green
    box1.LineWidth = 2;
    outliers1.MarkerEdgeColor = [0 0 0];
    box2.LineWidth = 2;
    box2.Color = [1 0 0];
    outliers2.MarkerEdgeColor = [1 0 0];
    box3.LineWidth = 2;
    box3.Color = [0 0 1];
    outliers3.MarkerEdgeColor = [0 0 1];    
end
fprintf('Mean difference between FG and KP 18-O day-by-day flux fis %.3f‰\n', nanmean((FGValues(:,1) - KPValues(:,1))))
fprintf('Mean difference between FG and KP 18-O flux with bootstrapping is %.3f‰\n', nanmean((d18OE_MixModel(:,1) - d18OE_FLUX(:,1))))

if display_histograms
    axes('Position',[.22 .27 .35 .35])
    hold on
    d18Hist = histogram(sigmaET_delta_FG(:,1),'Normalization', 'pdf', 'NumBins', 25);
    d18Hist.EdgeColor = [1 0 0];
    d18Hist.FaceColor = [1 0 0];
    d18Hist.EdgeAlpha = .5;%.2;
    d18Hist.FaceAlpha = .2;

    d18Hist = histogram(sigmaET_delta_KP(:,1),'Normalization', 'pdf', 'NumBins', 25);
    d18Hist.EdgeColor = [0 0 1];
    d18Hist.FaceColor = [0 0 1];
    d18Hist.EdgeAlpha = .8;
    d18Hist.FaceAlpha = .2;
    hold off
end

if display_kernelPDF
    axes('Position',[.10 .28 .35 .35])
    hold on
    x_values = xlimits18(1):.1:xlimits18(2);
    pd = fitdist(sigmaET_delta_FG(:,1),'Kernel', 'Kernel', 'normal', 'Bandwidth', 1);
    y = pdf(pd, x_values);    
    plot(x_values, y, 'LineWidth',2, 'Color', [1 0 0])
    
    pd = fitdist(sigmaET_delta_KP(:,1),'Kernel', 'Kernel', 'normal', 'Bandwidth', 1);
    x_values = xlimits18(1):.1:xlimits18(2);
    y = pdf(pd, x_values);
    hold on
    plot(x_values, y, 'LineWidth',2, 'Color', [0 0 1])
    hold off
end

%legend({"\sigma_{\delta_{FG}}", "\sigma_{\delta_{KP}}", "\sigma_{\delta_{FG}}", "\sigma_{\delta_{KP}}"}, 'Location', 'southwest')

% Plot additional features
xlabel('\sigma_{\delta^{18}O_{E}} (‰)')
ylabel('PDF')
xlim(xlimits18)
box on
ax = gca;
ax.YAxis.FontSize = 13;
ax.XAxis.FontSize = 13;
ax.XAxis.TickValues = 0:2:10;
set(gcf, 'Color', [1 1 1]);
grid on

set(gcf, 'Position',  [100, 100, 500, 400])

%% Plot error distributions for dD
clf
% Flux boxplot
if display_flux
    box on
    % boxplot(horzcat(dDE_MixModel, dDE_FLUX, [TopdD; nan(length(dDE_MixModel) - length(TopdD), 1)]),...
    boxplot(horzcat(dDE_MixModel, dDE_FLUX, TopVals_Averages(:,3)),...
        'Orientation', 'horizontal',...
        'labels', ["","", ""],...%'labels', ["KP","FG", "dDTop"],...
        'Colors', [0 0 0]);
    xlabel('\deltaD (‰)')
    ax = gca;
    ax.YAxis.FontSize = 13;
    ax.XAxis.FontSize = 13;
    set(gcf, 'Color', [1 1 1]);
    xlim([-150 5])
    
        % Boxplot properties
    a = get(get(gca,'children'),'children');   % Get the handles of all the objects
    t = get(a,'tag');   % List the names of all the objects 
    box1 = a(7);   % The 7th object is the first box
    outliers1 = a(1);
    box2 = a(8);   % The 8th object is the first box
    outliers2 = a(2);
    box3 = a(9);   % The 9th object is the first box
    outliers3 = a(3);
    % set(box1, 'Color', 'g');   % Set the color of the first box to green
    box1.LineWidth = 2;
    outliers1.MarkerEdgeColor = [0 0 0];
    box2.LineWidth = 2;
    box2.Color = [1 0 0];
    outliers2.MarkerEdgeColor = [1 0 0];
    box3.LineWidth = 2;
    box3.Color = [0 0 1];
    outliers3.MarkerEdgeColor = [0 0 1];  
end
fprintf('Mean difference between FG and KP 2-H day-by-day flux fis %.3f‰\n', nanmean((FGValues(:,1) - KPValues(:,1))))
fprintf('Mean difference between FG and KP 2-H flux with bootstrapping is %.3f‰\n', nanmean((d18OE_MixModel(:,1) - d18OE_FLUX(:,1))))

if display_histograms
    axes('Position',[.22 .3 .35 .35])
    hold on
    dDHist = histogram(sigmaET_delta_FG(:,2),'Normalization', 'pdf', 'NumBins', 25);
    dDHist.EdgeColor = [1 0 0];
    dDHist.FaceColor = [1 0 0];
    dDHist.EdgeAlpha = .8;
    dDHist.FaceAlpha = .2;

    dDHist = histogram(sigmaET_delta_KP(:,2),'Normalization', 'pdf', 'NumBins', 25);
    dDHist.EdgeColor = [0 0 1];
    dDHist.FaceColor = [0 0 1];
    dDHist.EdgeAlpha = .8;
    dDHist.FaceAlpha = .2;
    hold off
end

if display_kernelPDF
    axes('Position',[.10 .28 .35 .35])
    hold on
    x_values = xlimits2(1):.1:xlimits2(2);
    pd = fitdist(sigmaET_delta_FG(:,2),'Kernel', 'Kernel', 'normal', 'Bandwidth', 10);
    y = pdf(pd, x_values);    
    plot(x_values, y, 'LineWidth',2, 'Color', [1 0 0])
    
    pd = fitdist(sigmaET_delta_KP(:,2),'Kernel', 'Kernel', 'normal', 'Bandwidth', 10);
    x_values = xlimits2(1):.1:xlimits2(2);
    y = pdf(pd, x_values);
    hold on
    plot(x_values, y, 'LineWidth',2, 'Color', [0 0 1])
    hold off
end

%legend({"\sigma_{\delta_{FG}}", "\sigma_{\delta_{KP}}", "\sigma_{\delta_{FG}}", "\sigma_{\delta_{KP}}"}, 'Location', 'southwest')

% Plot additional features
xlabel('\sigma_{\deltaD_{E}} (‰)')
ylabel('PDF')
xlim(xlimits2)
box on
ax = gca;
ax.YAxis.FontSize = 13;
ax.XAxis.FontSize = 13;
ax.XAxis.TickValues = 0:20:100;
set(gcf, 'Color', [1 1 1]);
grid on

set(gcf, 'Position',  [100, 100, 500, 400])

fprintf('Mean difference between FG and KP 2-H day-by-day flux is %.3f‰\n', nanmean((FGValues(:,2) - KPValues(:,2))))
fprintf('Mean difference between FG and KP 2-H flux with bootstrapping is %.3f‰\n', nanmean((dDE_MixModel(:,1) - dDE_FLUX(:,1))))

%% Save error distributions
Error_distribution_FG = table(sigmaET_delta_FG(~isnan(sigmaET_delta_FG(:,1)),1), sigmaET_delta_FG(~isnan(sigmaET_delta_FG(:,1)),2),'VariableNames',["d18O";"dD"]);
Error_distribution_KP = table(sigmaET_delta_KP(~isnan(sigmaET_delta_KP(:,1)),1), sigmaET_delta_KP(~isnan(sigmaET_delta_KP(:,1)),2),'VariableNames',["d18O";"dD"]);

%% Bootstrapping for KP
% Reduce observations by available flux
d18OE_CG_reduced            = d18OE_CG(~isnan(d18OE_MixModel(:,1)), :);
dDE_CG_reduced              = dDE_CG(~isnan(d18OE_MixModel(:,1)), :);
d18OE_MixModel_reduced      = d18OE_MixModel(~isnan(d18OE_MixModel(:,1)), :);
dDE_MixModel_reduced        = dDE_MixModel(~isnan(d18OE_MixModel(:,1)), :);
sigmaET_delta_KP_reduced    = sigmaET_delta_KP(~isnan(d18OE_MixModel(:,1)), :);
WSs_reduced                 = WSs(~isnan(d18OE_MixModel(:,1)), :);

% nIterations = 20;
% Generate matrix of random indeces with the same size of original sample
rndIndexes = 1 + round((size(d18OE_CG_reduced, 1)-1)*rand([size(d18OE_CG_reduced, 1) nIterations]));

% Preallocate  memory for storing the
% distrbutions of kinetic fractionation values and other parameters
kinetic_18_16_sample_KP = nan(nIterations, 1);
kinetic_2_1_sample_KP = nan(nIterations, 1);
%kinetic_2_1_sample_expected = nan(nIterations, 1);
%d18OE_CG_Top_sample = kinetic_18_16_sample_KP;
%dDE_CG_Top_sample = kinetic_18_16_sample_KP;
RMSEvals_d18O_KP = nan(size(d18OE_CG_reduced, 2), 1)';
RMSEvals_dD_KP = RMSEvals_d18O_KP;
RMSE_18_16_sample_KP = nan(nIterations, 1);
RMSE_2_1_sample_KP = nan(nIterations, 1);

% Calculate weight vector based on error: vector of RSD normalized to 1
w18_16  = 1./abs((sigmaET_delta_KP_reduced(:, 1)./d18OE_MixModel_reduced));
w18_16  = w18_16./sum(w18_16);
w2_1  = 1./abs((sigmaET_delta_KP_reduced(:, 2)./dDE_MixModel_reduced));
w2_1  = w2_1./sum(w2_1);

%f = waitbar(0,'Please wait...');
% Use CG flux estimated from top inlet data
for k = 1 : nIterations
    % Step 1
    DF18 = sqrt((d18OE_CG_reduced(rndIndexes(:,k),:) - repmat(d18OE_MixModel_reduced(rndIndexes(:,k),:), 1, size(d18OE_CG_reduced(rndIndexes(:,k),:), 2))).^2);
    DF2 = sqrt((dDE_CG_reduced(rndIndexes(:,k),:) - repmat(dDE_MixModel_reduced(rndIndexes(:,k),:), 1, size(dDE_CG_reduced(rndIndexes(:,k),:), 2))).^2);
    % DF18 = abs((d18OE_CG_reduced(rndIndexes(:,k),:) - repmat(d18OE_MixModel_reduced(rndIndexes(:,k),:), 1, size(d18OE_CG_reduced(rndIndexes(:,k),:), 2))));
    % DF2 = abs((dDE_CG_reduced(rndIndexes(:,k),:) - repmat(dDE_MixModel_reduced(rndIndexes(:,k),:), 1, size(dDE_CG_reduced(rndIndexes(:,k),:), 2))));
    % Step 2
    MIN18 = min(DF18, [], 2);
    MIN2 = min(DF2, [], 2);
    % Step 3
    BM18 = DF18 == repmat(MIN18, 1, size(DF18, 2));
    BM2 = DF2 == repmat(MIN2, 1, size(DF2, 2));
    % Step 4
    BV18 = sum(d18OE_CG_reduced(rndIndexes(:,k),:).*BM18, 2);
    BV2 = sum(dDE_CG_reduced(rndIndexes(:,k),:).*BM2, 2);
    % Step 5
    BK18_KP = sum(repmat(kVals', size(BM18, 1), 1).*BM18, 2);
    BK2_KP = sum(repmat(kVals', size(BM2, 1), 1).*BM2, 2);
    % parameters.pptx"
    % Find minimum RMSE, associate kinetic fractionation parameter
    kinetic_18_16_sample_KP(k) = sum(w18_16.*BK18_KP);
    kinetic_2_1_sample_KP(k) = sum(w2_1.*BK2_KP);
    % Save RMSE for further analysis
    RMSE_18_16_sample_KP(k) = min(RMSEvals_d18O_KP);
    RMSE_2_1_sample_KP(k) = min(RMSEvals_dD_KP);
    clc
    fprintf('%.2f %%\n', 100*k/nIterations)
    %waitbar(k/nIterations, f, sprintf('Bootstrapping mean k-distribution using KP'));
end
%close(f)

%% Bootstrapping for FG
% Reduce observations by available flux
d18OE_FLUX_reduced      = d18OE_FLUX(~isnan(d18OE_FLUX(:,1)), :);
dDE_FLUX_reduced        = dDE_MixModel(~isnan(d18OE_FLUX(:,1)), :);
sigmaET_delta_FG_reduced    = sigmaET_delta_FG(~isnan(d18OE_FLUX(:,1)), :);

%nIterations = 20;
% Generate matrix of random indeces with the same size of original sample
rndIndexes = 1 + round((size(d18OE_CG_reduced, 1)-1)*rand([size(d18OE_CG_reduced, 1) nIterations]));

% Preallocate  memory for storing the
% distrbutions of kinetic fractionation values and other parameters
kinetic_18_16_sample_FG = nan(nIterations, 1);
kinetic_2_1_sample_FG = nan(nIterations, 1);
%kinetic_2_1_sample_expected = nan(nIterations, 1);
%d18OE_CG_Top_sample = kinetic_18_16_sample_KP;
%dDE_CG_Top_sample = kinetic_18_16_sample_KP;
RMSEvals_d18O_FG = nan(size(d18OE_CG_reduced, 2), 1)';
RMSEvals_dD_FG = RMSEvals_d18O_FG;
RMSE_18_16_sample_FG = nan(nIterations, 1);
RMSE_2_1_sample_FG = nan(nIterations, 1);

% Calculate weight vector based on error: vector of RSD normalized to 1
w18_16  = 1./abs((sigmaET_delta_FG_reduced(:, 1)./d18OE_FLUX_reduced));
w18_16  = w18_16./sum(w18_16);
w2_1  = 1./abs((sigmaET_delta_FG_reduced(:, 2)./dDE_FLUX_reduced));
w2_1  = w2_1./sum(w2_1);

% Use CG flux estimated from top inlet data
f = waitbar(0,'Please wait...');
for k = 1 : nIterations
   % Step 1
    DF18 = sqrt((d18OE_CG_reduced(rndIndexes(:,k),:) - repmat(d18OE_FLUX_reduced(rndIndexes(:,k),:), 1, size(d18OE_CG_reduced(rndIndexes(:,k),:), 2))).^2);
    DF2 = sqrt((dDE_CG_reduced(rndIndexes(:,k),:) - repmat(dDE_FLUX_reduced(rndIndexes(:,k),:), 1, size(dDE_CG_reduced(rndIndexes(:,k),:), 2))).^2);
    % DF18 = abs((d18OE_CG_reduced(rndIndexes(:,k),:) - repmat(d18OE_FLUX_reduced(rndIndexes(:,k),:), 1, size(d18OE_CG_reduced(rndIndexes(:,k),:), 2))));
    % DF2 = abs((dDE_CG_reduced(rndIndexes(:,k),:) - repmat(dDE_FLUX_reduced(rndIndexes(:,k),:), 1, size(dDE_CG_reduced(rndIndexes(:,k),:), 2))));
    % Step 2
    MIN18 = min(DF18, [], 2);
    MIN2 = min(DF2, [], 2);
    % Step 3
    BM18 = DF18 == repmat(MIN18, 1, size(DF18, 2));
    BM2 = DF2 == repmat(MIN2, 1, size(DF2, 2));
    % Step 4
    BV18 = sum(d18OE_CG_reduced(rndIndexes(:,k),:).*BM18, 2);
    BV2 = sum(dDE_CG_reduced(rndIndexes(:,k),:).*BM2, 2);
    % Step 5
    BK18_FG = sum(repmat(kVals', size(BM18, 1), 1).*BM18, 2);
    BK2_FG = sum(repmat(kVals', size(BM2, 1), 1).*BM2, 2);
    % parameters.pptx"
    % Find minimum RMSE, associate kinetic fractionation parameter
    kinetic_18_16_sample_FG(k) = sum(w18_16.*BK18_FG);
    kinetic_2_1_sample_FG(k) = sum(w2_1.*BK2_FG);
    % Save RMSE for further analysis
    RMSE_18_16_sample_FG(k) = min(RMSEvals_d18O_FG);
    RMSE_2_1_sample_FG(k) = min(RMSEvals_dD_FG);
    clc
    fprintf('%.2f %%\n', 100*k/nIterations)
    waitbar(k/nIterations, f, sprintf('Bootstrapping mean k-distribution using FG'));
end
close(f)

%% Descriptive statistics
clc
fprintf('Flux stats\n')
fprintf('       Mean    Median          IQR             sigmadE\n')
fprintf('------------------------------------------------\n')
IQR = icdf(fitdist(d18OE_FLUX_reduced,'Normal'),[0.25,0.75]);
SdF = mean(sigmaET_delta_FG_reduced);
fprintf('FGd18O  %.2f    %.2f     %.2f ; %.2f       %.2f \n', mean(d18OE_FLUX_reduced), median(d18OE_FLUX_reduced), IQR(1), IQR(2), SdF(1));
IQR = icdf(fitdist(dDE_FLUX_reduced,'Normal'),[0.25,0.75]);
SdF = mean(sigmaET_delta_FG_reduced);
fprintf('FGdD    %.2f   %.2f   %.2f ; %.2f        %.2f \n', mean(dDE_FLUX_reduced), median(dDE_FLUX_reduced), IQR(1), IQR(2), SdF(2));
fprintf('\n')
IQR = icdf(fitdist(d18OE_MixModel_reduced,'Normal'),[0.25,0.75]);
SdF = mean(sigmaET_delta_KP_reduced);
fprintf('KPd18O  %.2f    %.2f     %.2f ; %.2f       %.2f \n', mean(d18OE_MixModel_reduced), median(d18OE_MixModel_reduced), IQR(1), IQR(2), SdF(1));
IQR = icdf(fitdist(dDE_MixModel_reduced,'Normal'),[0.25,0.75]);
SdF = mean(sigmaET_delta_KP_reduced);
fprintf('KPdD    %.2f   %.2f   %.2f ; %.2f        %.2f \n', mean(dDE_MixModel_reduced), median(dDE_MixModel_reduced), IQR(1), IQR(2), SdF(2));

fprintf('\n')
IQR = icdf(fitdist(Topd18O,'Normal'),[0.25,0.75]);
%SdF = mean(sigmaET_delta_KP_reduced);
fprintf('Topd18O %.2f    %.2f     %.2f ; %.2f\n', mean(Topd18O), median(Topd18O), IQR(1), IQR(2));
IQR = icdf(fitdist(TopdD,'Normal'),[0.25,0.75]);
%SdF = mean(sigmaET_delta_KP_reduced);
fprintf('TopdD   %.2f   %.2f   %.2f ; %.2f\n', mean(TopdD), median(TopdD), IQR(1), IQR(2));