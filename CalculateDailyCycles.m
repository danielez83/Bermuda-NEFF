%% Daily cycles
hour_vect = hour(MeteoDataSYNC.Time);

hours_vect  = 0:23;
T_vect      = zeros(1,24);
RH_vect     = zeros(1,24);
RHsst_vect  = zeros(1,24);
WS_vect     = zeros(1,24);

for i = 1:24
    T_vect(i) = nanmean(MeteoDataSYNC.T(hour(MeteoDataSYNC.Time) == hours_vect(i)));
    RHsst_vect(i) = nanmean(MeteoDataSYNC.RHSST(hour(MeteoDataSYNC.Time) == hours_vect(i)));
    RH_vect(i) = nanmean(MeteoDataSYNC.RH(hour(MeteoDataSYNC.Time) == hours_vect(i)));
    WS_vect(i) = nanmean(MeteoDataSYNC.WS(hour(MeteoDataSYNC.Time) == hours_vect(i)));
end
%%
plot(hours_vect, WS_vect)