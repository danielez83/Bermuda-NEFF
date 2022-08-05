%% Configuration
k18                         = 5.2
k2                          = 4.3
% VSMOW ratio
R18SMOW                     = 2005.2e-6;
R2SMOW                      = 155.76e-6;

%% Include functions and script
if ~contains(string(path), 'Functions')
    addpath(genpath('Functions')) % Include Functions to search path
end



%% Load d18O ocean composition
% Download this data from here (LeGrande and Schmidt 2006)
% https://pubs.giss.nasa.gov/abs/le09100s.html
filename                    = '../../NetCDF data/LeGrandeCalculated_d18O_v1_1.nc';
%finfo = ncinfo(filename);
d18O_Ocean_Surface          = ncread(filename, 'd18o');
d18O_Ocean_Surface          = squeeze(d18O_Ocean_Surface(:,:,1));
d18O_Ocean_Surface_lat      = ncread(filename, 'lat');
d18O_Ocean_Surface_lon      = ncread(filename, 'lon');
%% Load water vapor dataset 1 (ACTIV CRUISE)
% Download this data from here (Benetti et al. 2017)
% http://cds-espri.ipsl.fr/isowvdataatlantic/
filename                    = '../Matlab Data/ACT_15_Min.mat';
d18O_to_dD                  = @(d18O) 6.57.*d18O-0.67; % Subpolar Gyre d18O vs dD relationship following Benetti et al. (2017)


load(filename)
d18O_data_1                 = d18O_VSMOW;
dD_data_1                   = dD_VSMOW;
latitude_data_1             = Longitude;
longitude_data_1            = Latitude;
SST_data_1                  = SST_OSTIA_Larger;
RH_SST_data_1               = RHS_10m_OSTIA;
WS_data_1                   = Wind_Speed_ERA;
dexcess_data_1             = dD_data_1 - 8* d18O_data_1;
% Preallocate memory
d18O_Ocean_data_1           = nan(length(latitude_data_1), 1);
dD_Ocean_data_1             = d18O_Ocean_data_1;
d18O_CA_data_1              = d18O_Ocean_data_1;
dD_CA_data_1                = d18O_Ocean_data_1;

start_idx   = 1;
stop_idx    = 8482;  %ok for ACTIV CRUISE.... length(latitude_data_1)
for j = start_idx:stop_idx
    if ~isnan(d18O_data_1(j)) && ~isnan(SST_data_1(j))
        % Calculate isotopic composition of ocean
        lat_idx = find(abs(latitude_data_1(j) - d18O_Ocean_Surface_lat) == min(abs(latitude_data_1(j) - d18O_Ocean_Surface_lat)));
        lon_idx = find(abs(longitude_data_1(j) - d18O_Ocean_Surface_lon) == min(abs(longitude_data_1(j) - d18O_Ocean_Surface_lon)));
        d18O_Ocean_data_1(j) = nanmean(nanmean(d18O_Ocean_Surface(lon_idx-1:lon_idx+1, lat_idx-1:lat_idx+1)));
        if d18O_Ocean_data_1(j) > + 1.5 || d18O_Ocean_data_1(j) < - 2
            d18O_Ocean_data_1(j) = NaN;
        else
            dD_Ocean_data_1(j) = d18O_to_dD(d18O_Ocean_data_1(j));
        end
        % Convert current isotopes ratios from delta to R
        Rsw18   = ((d18O_Ocean_data_1(j)/1000)+1)*R18SMOW;
        Rsw2    = ((dD_Ocean_data_1(j)/1000)+1)*R2SMOW;
        R18_Vapor_CGCA = MJ79_CGCA(RH_SST_data_1(j), SST_data_1(j) + 273.15, SST_data_1(j) + 273.15, Rsw18, k18, 18);
        d18O_CA_data_1(j, 1) = ((R18_Vapor_CGCA/R18SMOW)-1)*1e3;
        R2_Vapor_CGCA = MJ79_CGCA(RH_SST_data_1(j), SST_data_1(j) + 273.15, SST_data_1(j) + 273.15, Rsw2, k2, 2);
        dD_CA_data_1(j, 1) = ((R2_Vapor_CGCA/R2SMOW)-1)*1e3;
    end
end
dexcess_CA_data_1 = dD_CA_data_1 - 8*d18O_CA_data_1;
%% Load water vapor dataset 2 (RARA CRUISE)
% Download this data from here (Benetti et al. 2017)
% http://cds-espri.ipsl.fr/isowvdataatlantic/
filename                    = '../Matlab Data/RAR_15_Min.mat';
d18O_to_dD                  = @(d18O) ((5.80+6.36)/2).*d18O + ((0.91+0.50)/2); % Mean of Western and Eastern NASTG d18O vs dD relationship following Benetti et al. (2017)


load(filename)
d18O_data_2                 = d18O_VSMOW;
dD_data_2                   = dD_VSMOW;
latitude_data_2             = Longitude;
longitude_data_2            = Latitude;
SST_data_2                  = SST_TSG;
RH_SST_data_2               = RHS_10m;
WS_data_2                   = Wind_Speed;
dexcess_data_2             = dD_data_2 - 8* d18O_data_2;
% Preallocate memory
d18O_Ocean_data_2           = nan(length(latitude_data_2), 1);
dD_Ocean_data_2             = d18O_Ocean_data_2;
d18O_CA_data_2              = d18O_Ocean_data_2;
dD_CA_data_2                = d18O_Ocean_data_2;

start_idx   = 1;
stop_idx    = length(latitude_data_2);
for j = start_idx:stop_idx
    if ~isnan(d18O_data_2(j)) && ~isnan(SST_data_2(j))
        % Calculate isotopic composition of ocean
        lat_idx = find(abs(latitude_data_2(j) - d18O_Ocean_Surface_lat) == min(abs(latitude_data_2(j) - d18O_Ocean_Surface_lat)));
        lon_idx = find(abs(longitude_data_2(j) - d18O_Ocean_Surface_lon) == min(abs(longitude_data_2(j) - d18O_Ocean_Surface_lon)));
        d18O_Ocean_data_2(j) = nanmean(nanmean(d18O_Ocean_Surface(lon_idx-1:lon_idx+1, lat_idx-1:lat_idx+1)));
        if d18O_Ocean_data_2(j) > + 1.5 || d18O_Ocean_data_2(j) < - 2
            d18O_Ocean_data_2(j) = NaN;
        else
            dD_Ocean_data_2(j) = d18O_to_dD(d18O_Ocean_data_2(j));
        end
        % Convert current isotopes ratios from delta to R
        Rsw18   = ((d18O_Ocean_data_2(j)/1000)+1)*R18SMOW;
        Rsw2    = ((dD_Ocean_data_2(j)/1000)+1)*R2SMOW;
        R18_Vapor_CGCA = MJ79_CGCA(RH_SST_data_2(j), SST_data_2(j) + 273.15, SST_data_2(j) + 273.15, Rsw18, k18, 18);
        d18O_CA_data_2(j, 1) = ((R18_Vapor_CGCA/R18SMOW)-1)*1e3;
        R2_Vapor_CGCA = MJ79_CGCA(RH_SST_data_2(j), SST_data_2(j) + 273.15, SST_data_2(j) + 273.15, Rsw2, k2, 2);
        dD_CA_data_2(j, 1) = ((R2_Vapor_CGCA/R2SMOW)-1)*1e3;
    end
end
dexcess_CA_data_2 = dD_CA_data_2 - 8*d18O_CA_data_2;
%% Load water vapor dataset 3 (STRASSE CRUISE)
% Download this data from here (Benetti et al. 2017)
% http://cds-espri.ipsl.fr/isowvdataatlantic/
filename                    = '../Matlab Data/STR_15_Min.mat';
d18O_to_dD                  = @(d18O) 5.80.*d18O + 0.91; % Eastern NASTG d18O vs dD relationship following Benetti et al. (2017)


load(filename)
d18O_data_3                 = d18O_VSMOW;
dD_data_3                   = dD_VSMOW;
latitude_data_3             = Longitude;
longitude_data_3            = Latitude;
SST_data_3                  = SST_TSG;
RH_SST_data_3               = RHS_10m;
WS_data_3                   = Wind_Speed;
dexcess_data_3             = dD_data_3 - 8* d18O_data_3;
% Preallocate memory
d18O_Ocean_data_3           = nan(length(latitude_data_3), 1);
dD_Ocean_data_3             = d18O_Ocean_data_3;
d18O_CA_data_3              = d18O_Ocean_data_3;
dD_CA_data_3                = d18O_Ocean_data_3;

start_idx   = 1;
stop_idx    = length(latitude_data_3);
for j = start_idx:stop_idx
    if ~isnan(d18O_data_3(j)) && ~isnan(SST_data_3(j))
        % Calculate isotopic composition of ocean
        lat_idx = find(abs(latitude_data_3(j) - d18O_Ocean_Surface_lat) == min(abs(latitude_data_3(j) - d18O_Ocean_Surface_lat)));
        lon_idx = find(abs(longitude_data_3(j) - d18O_Ocean_Surface_lon) == min(abs(longitude_data_3(j) - d18O_Ocean_Surface_lon)));
        d18O_Ocean_data_3(j) = nanmean(nanmean(d18O_Ocean_Surface(lon_idx-1:lon_idx+1, lat_idx-1:lat_idx+1)));
        if d18O_Ocean_data_3(j) > + 1.5 || d18O_Ocean_data_3(j) < - 2
            d18O_Ocean_data_3(j) = NaN;
        else
            dD_Ocean_data_3(j) = d18O_to_dD(d18O_Ocean_data_3(j));
        end
        % Convert current isotopes ratios from delta to R
        Rsw18   = ((d18O_Ocean_data_3(j)/1000)+1)*R18SMOW;
        Rsw2    = ((dD_Ocean_data_3(j)/1000)+1)*R2SMOW;
        R18_Vapor_CGCA = MJ79_CGCA(RH_SST_data_3(j), SST_data_3(j) + 273.15, SST_data_3(j) + 273.15, Rsw18, k18, 18);
        d18O_CA_data_3(j, 1) = ((R18_Vapor_CGCA/R18SMOW)-1)*1e3;
        R2_Vapor_CGCA = MJ79_CGCA(RH_SST_data_3(j), SST_data_3(j) + 273.15, SST_data_3(j) + 273.15, Rsw2, k2, 2);
        dD_CA_data_3(j, 1) = ((R2_Vapor_CGCA/R2SMOW)-1)*1e3;
    end
end
dexcess_CA_data_3 = dD_CA_data_3 - 8*d18O_CA_data_3;

%% Create regression models
% Names
name_data_1         = 'ACTIV';
name_data_2         = 'RARA';
name_data_3         = 'STRASSE';
% ACTIV
obs_d_vs_h_data_1   = fitlm(RH_SST_data_1/100, dexcess_data_1);
CA_d_vs_h_data_1    = fitlm(RH_SST_data_1/100, dexcess_CA_data_1); 
MAE_data_1          = nanmean(abs(dexcess_data_1 - dexcess_CA_data_1));
RMSE_data_1         = fitlm(dexcess_data_1, dexcess_CA_data_1);
RMSE_data_1         = RMSE_data_1.RMSE;
% RARA
obs_d_vs_h_data_2   = fitlm(RH_SST_data_2/100, dexcess_data_2);
CA_d_vs_h_data_2    = fitlm(RH_SST_data_2/100, dexcess_CA_data_2); 
MAE_data_2          = nanmean(abs(dexcess_data_2 - dexcess_CA_data_2));
RMSE_data_2         = fitlm(dexcess_data_2, dexcess_CA_data_2);
RMSE_data_2         = RMSE_data_2.RMSE;
% STRASSE
obs_d_vs_h_data_3   = fitlm(RH_SST_data_3/100, dexcess_data_3);
CA_d_vs_h_data_3    = fitlm(RH_SST_data_3/100, dexcess_CA_data_3);
MAE_data_3          = nanmean(abs(dexcess_data_3 - dexcess_CA_data_3));
RMSE_data_3         = fitlm(dexcess_data_3, dexcess_CA_data_3);
RMSE_data_3         = RMSE_data_3.RMSE;
% Print table for manuscript
fprintf('-------------------------------------------------------------------------------------------------\n')
fprintf('Cruise     Slope_obs     Intercept_obs     R2       Slope_CA     Intercept_CA     MAE     RMSE\n')
fprintf('-------------------------------------------------------------------------------------------------\n')
fprintf('%s \t %.2f ± %.2f \t %.2f ± %.2f \t %.2f \t %.2f ± %.2f \t %.2f ± %.2f \t %.2f \t %.2f\n',...
        name_data_1, ...
        obs_d_vs_h_data_1.Coefficients.Estimate(2), obs_d_vs_h_data_1.Coefficients.SE(2), ...
        obs_d_vs_h_data_1.Coefficients.Estimate(1), obs_d_vs_h_data_1.Coefficients.SE(1), ...
        obs_d_vs_h_data_1.Rsquared.Ordinary, ...
        CA_d_vs_h_data_1.Coefficients.Estimate(2), CA_d_vs_h_data_1.Coefficients.SE(2), ...
        CA_d_vs_h_data_1.Coefficients.Estimate(1), CA_d_vs_h_data_1.Coefficients.SE(1), ...
        MAE_data_1, ...
        RMSE_data_1)
    
fprintf('%s \t %.2f ± %.2f \t %.2f ± %.2f \t %.2f \t %.2f ± %.2f \t %.2f ± %.2f \t %.2f \t %.2f\n',...
        name_data_2, ...
        obs_d_vs_h_data_2.Coefficients.Estimate(2), obs_d_vs_h_data_2.Coefficients.SE(2), ...
        obs_d_vs_h_data_2.Coefficients.Estimate(1), obs_d_vs_h_data_2.Coefficients.SE(1), ...
        obs_d_vs_h_data_2.Rsquared.Ordinary, ...
        CA_d_vs_h_data_2.Coefficients.Estimate(2), CA_d_vs_h_data_2.Coefficients.SE(2), ...
        CA_d_vs_h_data_2.Coefficients.Estimate(1), CA_d_vs_h_data_2.Coefficients.SE(1), ...
        MAE_data_2, ...
        RMSE_data_2)
        
fprintf('%s  %.2f ± %.2f \t %.2f ± %.2f \t %.2f \t %.2f ± %.2f \t %.2f ± %.2f \t %.2f \t %.2f\n',...
        name_data_3, ...
        obs_d_vs_h_data_3.Coefficients.Estimate(2), obs_d_vs_h_data_3.Coefficients.SE(2), ...
        obs_d_vs_h_data_3.Coefficients.Estimate(1), obs_d_vs_h_data_3.Coefficients.SE(1), ...
        obs_d_vs_h_data_3.Rsquared.Ordinary, ...
        CA_d_vs_h_data_3.Coefficients.Estimate(2), CA_d_vs_h_data_3.Coefficients.SE(2), ...
        CA_d_vs_h_data_3.Coefficients.Estimate(1), CA_d_vs_h_data_3.Coefficients.SE(1), ...
        MAE_data_3, ...
        RMSE_data_3)
    
fprintf('%s  %.2f ± %.2f \t %.2f ± %.2f \t %.2f \t %.2f ± %.2f \t %.2f ± %.2f \t %.2f \t %.2f\n',...
        'This s.', ...
        -48.34, 0.23, ...
        47.91, 0.16, ...
        0.83, ...
        -36.22, 0.06, ...
        35.48, 0.04, ...
        3.82, ...
        4.46)
    
%% Create plot for manuscript
ylimits = [-10, 30];
xlimits = [.4, 1.1];
figure(1)
f = gcf;
f.Position = [118 205 450 514];

subplot(3,1,1)
    scatter(RH_SST_data_1/100, dexcess_data_1, 64, ...
            'MarkerEdgeColor',[0 0 0],...
            'MarkerFaceColor',[1 1 1],...
            'LineWidth',1.5, ...
            'MarkerFaceAlpha', 0.3)
    hold on
    plot(RH_SST_data_1/100, obs_d_vs_h_data_1.Fitted, 'LineWidth', 2, 'Color', [1,1,1])
    plot(RH_SST_data_1/100, CA_d_vs_h_data_1.Fitted, 'LineWidth', 2, 'Color', [1,0,0])
    hold off
    ylim(ylimits)
    ylabel('d-excess (‰)')
    xlim(xlimits)
    set(gca, 'XTickLabel', [])
    text(1, 20, name_data_1)
    box on
    grid on

subplot(3,1,2)
    scatter(RH_SST_data_2/100, dexcess_data_2, 64, ...
            'MarkerEdgeColor',[0 0 0],...
            'MarkerFaceColor',[1 1 1],...
            'LineWidth',1.5, ...
            'MarkerFaceAlpha', 0.3)
    hold on
    plot(RH_SST_data_2/100, obs_d_vs_h_data_2.Fitted, 'LineWidth', 2, 'Color', [1,1,1])
    plot(RH_SST_data_2/100, CA_d_vs_h_data_2.Fitted, 'LineWidth', 2, 'Color', [1,0,0])
    hold off
    ylim(ylimits)
    ylabel('d-excess (‰)')
    xlim(xlimits)
    set(gca, 'XTickLabel', [])
    text(1, 20, name_data_2)
    box on
    grid on
    
subplot(3,1,3)
    scatter(RH_SST_data_3/100, dexcess_data_3, 64, ...
            'MarkerEdgeColor',[0 0 0],...
            'MarkerFaceColor',[1 1 1],...
            'LineWidth',1.5, ...
            'MarkerFaceAlpha', 0.3)
    hold on
    plot(RH_SST_data_3/100, obs_d_vs_h_data_3.Fitted, 'LineWidth', 2, 'Color', [1,1,1])
    plot(RH_SST_data_3/100, CA_d_vs_h_data_3.Fitted, 'LineWidth', 2, 'Color', [1,0,0]);
    
    hold off
    ylim(ylimits)
    ylabel('d-excess (‰)')
    xlim(xlimits)
    xlabel('\it{h} (-)')
    text(1, 20, name_data_3)
    box on
    grid on


