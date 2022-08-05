% Author:      Daniele Zannoni
% Name:        ImportBermudaDEM.m
% Description: Import Bermuda coastal DEM
%              
% Date:        Last revision 11/10/2020

tic
%% Include functions and script
if ~contains(string(path), 'Functions')
    addpath(genpath('Functions')) % Include Functions to search path
end  
%% Configuration
SiteCoordinates = [32.2644427 -64.8941763 29]; % Latitude Longitude Z (m ASL)
% Variables that do not require timeseries extraction
%noTSVarNames = {'time', 'latitude', 'longitude'};
%dt = minutes(30);

%% Custom colormap
size_cmap = 25;
rchan = linspace(0.8750, 0, size_cmap);
gchan = linspace(1, 0, size_cmap);
bchan = linspace(1, 0.5, size_cmap);
cmap_blue = [rchan', gchan', bchan'];
%% Load NetCDF file
% https://www.ncei.noaa.gov/metadata/geoportal/rest/metadata/item/gov.noaa.ngdc.mgg.dem:5010/html
filename = '../../NetCDF data/bermuda_3_msl_2013.nc';
finfo = ncinfo(filename);

% Get variable name info
for j = 1:size(finfo.Variables, 2)
    VarName = finfo.Variables(j).Name;
    eval(sprintf('%s = ncread("%s", "%s");', VarName, filename, VarName));
    % Preallocate memory for time series
    %if sum(strcmp(VarName,noTSVarNames)) == 0
    %    eval(sprintf('%s_TS = zeros(length(%s), 1);', VarName, VarName));
    %end
end

%% Plot positions
BoxCoordinates = [31.6 32.6; -65.5 -64.1]; % Latitude limits; Longitude limits

h = worldmap(BoxCoordinates(1,:) ,BoxCoordinates(2,:));
% Change colors
p = findobj(h,'type','patch'); % Find background
%set(p,'FaceColor',[0.0745 0.6235 1.0000]); % Change background to white

coast = shaperead('../../BermudaIsland/nw036zp7611.shp','UseGeoCoords',true);%,'RecordNumbers',2);


%% Plot points
HogReefPOS = [32.46, -64.83]; % https://www.ncei.noaa.gov/access/ocean-carbon-data-system/oceans/Moorings/Hog_Reef.html
CrescentReefPOS = [32.40, -64.79]; % https://www.ncei.noaa.gov/access/ocean-carbon-data-system/oceans/Moorings/Crescent_Reef.html
BATSoceanPOS = [31.69 -64.17]; % approximated from here http://bats.bios.edu/about/
StGeorgePOS = [32.38, -64.67]; % St.George harbor
HogReefCOL = 'r';
CrescentReefCOL = 'g';
StGeorgeCOL = 'b';
BATSoceanCOL = 'y';%;[0 0.4471 0.7412];
BATSoceanSYM = '^';
MkrSyze = 16;
geoshow(HogReefPOS(1), HogReefPOS(2), 'DisplayType', 'point','Marker','o', 'MarkerEdgeColor','k','MarkerFaceColor', HogReefCOL, 'MarkerSize', MkrSyze)
geoshow(CrescentReefPOS(1), CrescentReefPOS(2), 'DisplayType', 'point','Marker','o', 'MarkerEdgeColor','k','MarkerFaceColor', CrescentReefCOL, 'MarkerSize', MkrSyze)
geoshow(StGeorgePOS(1), StGeorgePOS(2), 'DisplayType', 'point','Marker','o', 'MarkerEdgeColor','k','MarkerFaceColor',StGeorgeCOL, 'MarkerSize', MkrSyze)
geoshow(BATSoceanPOS(1), BATSoceanPOS(2), 'DisplayType', 'point','Marker',BATSoceanSYM, 'MarkerEdgeColor','k','MarkerFaceColor',BATSoceanCOL, 'MarkerSize', ...
    MkrSyze + 2, "LineWidth", 2)
legend(["OSTIA";"Hog Reef";"Crescent Reef";"St.George Harnour"; "BATS"]);

%% Convert time vector into datetime vector
Band1corr = Band1;
Band1corr(Band1corr>=0) = -0.1;
Band1corr = -1*Band1corr;
R = georefcells([lat(1) lat(end)], [lon(1) lon(end)], [size(Band1, 2) size(Band1, 1)],'ColumnsStartFrom','south');
geoshow(Band1corr', R, 'DisplayType', 'texture')
geoshow(coast, 'FaceColor', [1 1 1], 'EdgeColor', [0, 0, 0], 'LineWidth', 1.5)
%colormap(flip(haxby(50)))
colormap(cmap_blue)
set(gca,'ColorScale','log')
colorbar

geoshow(HogReefPOS(1), HogReefPOS(2), 'DisplayType', 'point','Marker','o', 'MarkerEdgeColor','k','MarkerFaceColor', HogReefCOL, 'MarkerSize', MkrSyze)
geoshow(CrescentReefPOS(1), CrescentReefPOS(2), 'DisplayType', 'point','Marker','o', 'MarkerEdgeColor','k','MarkerFaceColor', CrescentReefCOL, 'MarkerSize', MkrSyze)
geoshow(StGeorgePOS(1), StGeorgePOS(2), 'DisplayType', 'point','Marker','o', 'MarkerEdgeColor','k','MarkerFaceColor',StGeorgeCOL, 'MarkerSize', MkrSyze)
geoshow(BATSoceanPOS(1), BATSoceanPOS(2), 'DisplayType', 'point','Marker',BATSoceanSYM, 'MarkerEdgeColor','k','MarkerFaceColor',BATSoceanCOL, 'MarkerSize', MkrSyze)

a = legend;
legend(a.String(1:5))