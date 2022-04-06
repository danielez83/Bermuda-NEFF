% Author:      Daniele Zannoni
% Name:        LoadAndFilterData.m
% Description: This script load all the required data to run the other
%              scripts.
% Date:        Last revision 08/07/2021

%% Include functions and script
if ~contains(string(path), 'Functions')
    addpath(genpath('Functions')) % Include Functions to search path
end    

%% Filter Configuration (activate/deactivate filters)
ShowFilterImpact    = 1; % 1 ON, 0 OFF
FilterWD            = 1; % 1 ON, 0 OFF
FilterWS            = 1; % 1 ON, 0 OFF
FilterP             = 1; % 1 ON, 0 OFF
FilterRH            = 1; % 1 ON, 0 OFF
Filterq             = 1; % 1 ON, 0 OFF
Filterd18O          = 1; % 1 ON, 0 OFF
FilterdD            = 1; % 1 ON, 0 OFF
DayTimeOffset       = 2; % number of hours after sunrise and before sunset (2 is OK), USE 9999 value to exclude filter

% Wind -------------------------------------------------------------------
GoodWindSector          = [180 340]; % clockwise direction
WS_Lower_Threshold      = .1;       % lower wind thresholsd
WS_Upper_Threshold      = 15;        % higher wind threshold
% Precipitation -----------------------------------------------------------
hRainOffset             = 2;        % remove datapoint after each rain for n 
                                    % hours (e.g. if set to 2 and rain occured
                                    % at 18:00 the followind data point are
                                    % removed 18:00, 19:00, 20:00
                                    % If  > 0 data number reduced too
                                    % much...
% Relative Humidity -------------------------------------------------------
RH_Upper_Threshold      = 100;       % Subset data lower thant this RH threshold
RH_Lower_Threshold      = 0;       % Subset data higer thant this RH threshold

% Vertical profile structure  ---------------------------------------------
q_diff_threshold       = 100;      % Minimum PPMV difference between 
                                    % top and bottom inlet
d18O_diff_threshol     = 0.1;        % Minimum d18O difference (‰) between top and bottom inlet
                                    % Better if this value is two time the
                                    % d18O measurement uncertainty (top -
                                    % bottom), value from L2120-i datasheet
                                    
dD_diff_threshol       = 1;         % Minimum dD difference (‰) between top and bottom inlet
                                    % Better if this value is two time the
                                    % dD measurement uncertainty (top -
                                    % bottom), value from L2120-i datasheet

% Ocean isotopic composition  --------------------------------------------
UseConstantOceanDelta   = 1;                    % Turn on/off ocean isotopic composition variability
UseSSToffset            = 0;                    % Turn on/off ocean isotopic composition variability
SeaWater_d18O           = 1.0900;               % [‰], estimated from salinity, use this as a CONSTANT VALUE
SeaWater_dD             = 7.2450;               % [‰], estimated from salinity, use this as a  CONSTANT VALUE
Sal2d18O                = @(sal) 0.32*sal-10.5; % From Benetti et al. (2017) s-delta relationship
Sal2dD                  = @(sal) 2.04*sal-65.8; % From Benetti et al. (2017) s-delta relationship
SeaWater_d18O_offset    = -0.06;
SeaWater_dD_offset      = -0.38;
% Sal2d18O                = @(sal) 0.39*sal-13.1; % From BATS s-delta relationship
% Sal2dD                  = @(sal) 8.64*sal-307.6; % From BATS s-delta relationship
% Sal2d18O    = @(sal) 0.3605*sal-11.849; % From Benetti et al. (2017) + BATS 2012 data s-delta relationship, see OceanIsotopicCompostion.xlsx Tab Benetti et al 2017
% Sal2dD      = @(sal) 2.139*sal-69.539; % From Benetti et al. (2017) + BATS 2012 data  s-delta relationship, see OceanIsotopicCompostion.xlsx Tab Benetti et al 2017

% SST configuration   ----------------------------------------------------
UseOSTIASST             = 1;                   % Set to 1 to use OSTIA SST
                                               % Set to 0 to use mean SST ftom Hog Reef, Crescent Reef and St.George 
                                               
% Offset for debugging (remember to set all values here to 0 after debugging!) ---------------------------  
PicarroOffsetEnable = 1; % This switch ON/OFF all the following offsets
TopH2O_Offset       = 0; % 
Topd18O_Offset      = 0; %
TopdD_Offset        = 0; %
BottomH2O_Offset    = 0; %
Bottomd18O_Offset   = +0.07; % old offset is +0.07
BottomdD_Offset     = +0.75; % old offset is +0.77

% Other configuration  --------------------------------------------
GMWL_Slope              = 8;        % Has no impact for the determination of
                                    % kinetic_2 fractionation factor
KernelBandwidth        = 1;         % Kernel bandiwth to fit kinetic factors distribution                                    
% SeaWater_d18O          = 1.3516 + d18ObiasOW;    % [‰]
% SeaWater_dD            = 8.5881 + dDbiasOW;      % [‰]
MinThetaN              = -3;
MaxThetaN              = 3;

% kinetic fractionation values
MinkVal                 = -60;
MaxkVal                 = 60;
kVals = (MinkVal:0.01:MaxkVal)';

% Molecular diffusivities from Merlivat 1978
Di1816 = 0.9727; % From Horita 2008 is Di1816 = 0.9727;
Di21 = 0.9757; % From Horita 2008 is Di21   = 0.9757;
D18_D16 = (1 - Di1816) * 1000; % [‰]
D2_D1 = (1 - Di21) * 1000; % [‰]

% VSMOW ratio
R18SMOW     = 2005.2e-6;
R2SMOW     = 155.76e-6;
%R2SMOW     = 2*155.76e-6;


% CG model, use top or bottom data? (Top is the deault choice)
inletForCG      = 'top'; % 'top' or 'bottom'
% WV analysis, use top or bottom data? (Top is the deault choice)
inlet           = 'top'; % 'top' or 'bottom'


if ~exist('SensitivityStudySST', 'var') && ~exist('SensitivityStudyOceanR', 'var')
    tic
end
%% On the fly functions

%% Load data
% ----------------- Picarro data --------------------
load('../Matlab data/BottomTopData_SYNC_HF_5minAVG.mat')
if PicarroOffsetEnable == 0 
    TopH2O_Offset       = 0;
    Topd18O_Offset      = 0;
    TopdD_Offset        = 0; 
    BottomH2O_Offset    = 0;
    Bottomd18O_Offset   = 0; 
    BottomdD_Offset     = 0; 
end
% ----------------- ERA5 Data --------------------
load('../Matlab data/ERA5DataSYNC_HF.mat');
load('../Matlab data/ERA5_OceanDataSYNC.mat')
% ----------------- Meteorological and SST observations --------------------
load('../Matlab data/MeteoDataSYNC_GridCellSST_1x1_HF.mat')       % Do not comment, SST is from OSTIA
% load('../Matlab data/BermudaMeteoDataSYNC_Local_HF.mat')          % Comment/Uncomment if you want to use SST measured St. George
load('../Matlab data/BermudaMeteoDataSYNC_Local_AVGSST_HogAndCrescent_HF.mat')     % Comment/Uncomment if you want to use and averaged SST from Crescent, Hog and St. George
% Change SST value from local to OSTIA
if UseOSTIASST == 0
    MeteoDataSYNC.SST = BermudaWS_DataSYNC.goodSST;                  
end
    
% ----------------- Salinity observations--------------------
load('../Matlab data/Salinity_SYNC.mat')


%% Sensitivity study
% SensitivityStudyOceanR = 0; % Set to 1 if you want to perform a sensitivity study for ocean isotopic composition variability
if exist('SensitivityStudySST', 'var') % perform a sensitivity study for SST variability
    MeteoDataSYNC.RHSST     = NewRHSST;
    MeteoDataSYNC.SST       = NewSST;
end
if exist('SensitivityStudyOceanR', 'var') % perform a sensitivity study for SST variability
    SeaWater_d18O           = NewSeaWater_d18O;
    SeaWater_dD             = NewSeaWater_dD;
end
%% Subset data by wind direction without rain
if hRainOffset > 0
    % Build a logic vector with n + shift dimension
    BooleanIndexes = zeros(hRainOffset+2, length(MeteoDataSYNC.AccRain) + hRainOffset);
    BooleanIndexes(1, 1:end - hRainOffset) = MeteoDataSYNC.AccRain ~= 0;
    for j = 1 : hRainOffset
        BooleanIndexes(1+j, 1 + j: (length(BooleanIndexes) - hRainOffset + j)) = MeteoDataSYNC.AccRain ~= 0;
    end
    BooleanIndexes(end, :) = sum(BooleanIndexes(1:end-1, :));
    BooleanIndexes(end, :) = BooleanIndexes(end, :)>0;
    BooleanIndexes = BooleanIndexes(end, 1:length(BooleanIndexes)-hRainOffset);
    % Subset data
    %OceanSectorIndexes = MeteoDataSYNC.WD > GoodWindSector(1) & MeteoDataSYNC.WD < GoodWindSector(2) & BooleanIndexes' ~= 0;
    NotRain = ~(BooleanIndexes' ~= 0);
else
    %OceanSectorIndexes = MeteoDataSYNC.WD > GoodWindSector(1) & MeteoDataSYNC.WD < GoodWindSector(2) & MeteoDataSYNC.AccRain == 0;
    NotRain = MeteoDataSYNC.AccRain == 0;
end

%% Subset data by ocean sector
if FilterWD == 1
    OceanSectorIndexes = MeteoDataSYNC.WD > GoodWindSector(1) & MeteoDataSYNC.WD < GoodWindSector(2);
else
    OceanSectorIndexes = ones(length(MeteoDataSYNC.WD), 1);
end
%% Only data points that show decreasing q with height
if Filterq == 1
    Dec_q_with_heightIndex = TopInletSYNC.H2O < BottomInletSYNC.H2O - q_diff_threshold;
else
    Dec_q_with_heightIndex = ones(length(TopInletSYNC.H2O), 1);
end

%% Olny data points decreasing/increasing deltas with height larger than instrumental uncertainty
% Dec_deltas_with_heightIndex = (TopInletSYNC.d18O < BottomInletSYNC.d18O - d18O_diff_threshol) & ...
%                               (TopInletSYNC.dD < BottomInletSYNC.dD - dD_diff_threshol) ;
if Filterd18O & FilterdD == 1
    Dec_deltas_with_heightIndex = (abs(TopInletSYNC.d18O - BottomInletSYNC.d18O) > d18O_diff_threshol) & ...
                                  (abs(TopInletSYNC.dD - BottomInletSYNC.dD) > dD_diff_threshol) ;
else
    Dec_deltas_with_heightIndex = ones(length(TopInletSYNC.H2O), 1);
end
%% Only data points during quiscient conditions (wind speed less than 1 m/s)
if FilterWS == 1 
    WS_VeryLow = MeteoDataSYNC.WS < WS_Upper_Threshold & MeteoDataSYNC.WS > WS_Lower_Threshold;
else
    WS_VeryLow = ones(length(MeteoDataSYNC.WS), 1);
end

%% Only data points within ranged RH 
RH_Threshold = MeteoDataSYNC.RH < RH_Upper_Threshold & MeteoDataSYNC.RH > RH_Lower_Threshold;

%% Only daytime observations
if DayTimeOffset ~= 9999
    UTCzone = -4;
    [SRISE,SSET,NOON] = sunrise(32.26, -64.88, 0, UTCzone, BottomInletSYNC.Date);
    SRISE = datetime(datestr(SRISE));
    SSET = datetime(datestr(SSET));
    DayTime_Threshols = (hour(MeteoDataSYNC.Time) + UTCzone) > (hour(SRISE) + DayTimeOffset) & (hour(MeteoDataSYNC.Time) + UTCzone) < (hour(SSET) - DayTimeOffset);
else
    DayTime_Threshols = ones(length(BottomInletSYNC.Date), 1);
end

%% Show percentage impact of filtering
if ShowFilterImpact == 1
    rel_imp = 100 * (1 - sum(DayTime_Threshols)/size(BermudaWS_DataSYNC, 1));
    fprintf('Time threshold impact: %.0f\n', rel_imp);
    rel_imp = 100 * (1 - sum(DayTime_Threshols & OceanSectorIndexes)/size(BermudaWS_DataSYNC, 1));
    fprintf('Wind sector threshold impact: %.0f\n', rel_imp);
    rel_imp = 100 * (1 - sum(DayTime_Threshols & OceanSectorIndexes ...
                             & Dec_deltas_with_heightIndex)/size(BermudaWS_DataSYNC, 1));
    fprintf('Isotope diff. threshold impact: %.0f\n', rel_imp);
    rel_imp = 100 * (1 - sum(DayTime_Threshols & OceanSectorIndexes ...
                             & Dec_deltas_with_heightIndex & Dec_q_with_heightIndex )/size(BermudaWS_DataSYNC, 1));
    fprintf('Humidity diff. threshold impact: %.0f\n', rel_imp);
    rel_imp = 100 * (1 - sum(DayTime_Threshols & OceanSectorIndexes ...
                             & Dec_deltas_with_heightIndex & Dec_q_with_heightIndex ...
                             & NotRain)/size(BermudaWS_DataSYNC, 1));
    fprintf('No prec in the previous %.0f hours, impact: %.0f\n', hRainOffset, rel_imp);
end
%% Susbet data
GoodIndexes = NotRain & OceanSectorIndexes & Dec_q_with_heightIndex & WS_VeryLow & Dec_deltas_with_heightIndex & RH_Threshold & DayTime_Threshols;
if ShowFilterImpact == 1
    fprintf('Dataset reduced from %.0f to %.0f \n', size(BermudaWS_DataSYNC, 1), sum(GoodIndexes))
end
%GoodIndexes = OceanSectorIndexes & Dec_q_with_heightIndex & WS_VeryLow & Dec_deltas_with_heightIndex & RH_Threshold & DayTime_Threshols & BLHthreshold;
GoodIndexes = GoodIndexes(1:end);
% Picarro data
Topq        = TopInletSYNC.H2O(GoodIndexes)     + TopH2O_Offset;
Bottomq     = BottomInletSYNC.H2O(GoodIndexes)  + BottomH2O_Offset;
Topd18O     = TopInletSYNC.d18O(GoodIndexes)    + Topd18O_Offset;
Bottomd18O  = BottomInletSYNC.d18O(GoodIndexes) + Bottomd18O_Offset;
TopdD       = TopInletSYNC.dD(GoodIndexes)      + TopdD_Offset;
BottomdD    = BottomInletSYNC.dD(GoodIndexes)   + BottomdD_Offset;
ObsDates    = BottomInletSYNC.Date(GoodIndexes);
% Meteo data
AirT        = MeteoDataSYNC.T(GoodIndexes) + 273.15; % Convert celsius to K
SST         = MeteoDataSYNC.SST(GoodIndexes) + 273.15; % Convert celsius to K
RH          = MeteoDataSYNC.RH(GoodIndexes);
RHSST       = MeteoDataSYNC.RHSST(GoodIndexes);
% --------------- WIND DATA -----------------
WS          = MeteoDataSYNC.WS(GoodIndexes); WS = WScorr(WS, 10,50, 0.0002); % Implement Wind speed profile correction
% WS          = ERA5DataSYNC.WS10mERaA5(GoodIndexes);
WD          = MeteoDataSYNC.WD(GoodIndexes);
% -------------------------------------------
SLP         = MeteoDataSYNC.SLP(GoodIndexes);
% Additional ERA5 data
CBH         = ERA5DataSYNC.CBH(GoodIndexes);
PBLH        = ERA5DataSYNC.BLH(GoodIndexes);
TCC         = ERA5DataSYNC.TCC(GoodIndexes);
U10         = ERA5DataSYNC.WS10mERaA5(GoodIndexes);
% SWH         = ERA5_OceanDataSYNC.p140121(GoodIndexes); % Significant wave height
% MWP         = ERA5_OceanDataSYNC.mwp(GoodIndexes); % Mean wave period

%% RHSST Top Inlet
e_s = (exp(77.3450 + 0.0057 .* (MeteoDataSYNC.SST+273.15) - (7235./(MeteoDataSYNC.SST+273.15)))./(MeteoDataSYNC.SST+273.15).^8.2);
e_a = TopInletSYNC.H2O.*1e-6.*(MeteoDataSYNC.SLP*100);
TopInletRHSST =  100*e_a./e_s;
RHSSTtop = TopInletRHSST(GoodIndexes);

%% RHSST Bottom Inlet
e_s = (exp(77.3450 + 0.0057 .* (MeteoDataSYNC.SST+273.15) - (7235./(MeteoDataSYNC.SST+273.15)))./(MeteoDataSYNC.SST+273.15).^8.2);
e_a = BottomInletSYNC.H2O.*1e-6.*(MeteoDataSYNC.SLP*100);
BottomInletRHSST =  100*e_a./e_s;
RHSSTbottom = BottomInletRHSST(GoodIndexes);

%% Local Data
SSTreef = BermudaWS_DataSYNC.goodSST+273.15;
e_s = (exp(77.3450 + 0.0057 .* (SSTreef) - (7235./(SSTreef)))./(SSTreef).^8.2);
e_a = BottomInletSYNC.H2O.*1e-6.*(MeteoDataSYNC.SLP*100);
RHSSTreef = 100.*e_a./e_s;
RHSSTreefsubset = RHSSTreef(GoodIndexes);

%% Ocean isotopic composition
if UseConstantOceanDelta == 0
    % Estimate mean salinity in the area
    SeaWater_d18O_COPY   = Sal2d18O(nanmean(horzcat(Sal_DataSYNC.goodSALBATS, Sal_DataSYNC.goodSALCrescent, Sal_DataSYNC.goodSALHOG), 2));
    SeaWater_d18O_COPY   = SeaWater_d18O_COPY(GoodIndexes);
    SeaWater_dD_COPY     = Sal2dD(nanmean(horzcat(Sal_DataSYNC.goodSALBATS, Sal_DataSYNC.goodSALCrescent, Sal_DataSYNC.goodSALHOG), 2));
    SeaWater_dD_COPY     = SeaWater_dD_COPY(GoodIndexes);
end

if UseSSToffset == 1
    SeaWater_d18O = SeaWater_d18O + SeaWater_d18O_offset;
    SeaWater_dD = SeaWater_dD  + SeaWater_dD_offset;
end
