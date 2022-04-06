% Estimate wind speed at height z with log law
% https://websites.pmc.ucsc.edu/~jnoble/wind/extrap/
% z0 is 0.0002 for open sea
% WScorr, corrected wind speed value [m/s]
% WSref, measured wind speed [m/s]
% z, height of WScorr [m]
% zref, height of measured WSref [m]
% z0, roughness lenght [m]
function WScorr = WScorr(WSref, z, zref, z0)
    WScorr = WSref.*(log(z./z0)./log(zref./z0));
end