% Estimate isotopic composition of water vapor flux with
% the Craig-Gordon model (for 18O) using the Merlivat and Jouzel 1979
% formualtion (eq. 9)
% Computation is made for sea water (activity 0.98), see Horita 2008 for
% additional informations
% RH: 1 - 100 [%]
% Ait temperature [k]
% d18O of atmosphere [‰] 
% Sea/Water temperature [k]
% d18O of water [‰]
% kVal: kinetic fractionation value following MJ79 [‰]  
% - RH is not normalized to SST, it will be normalized to accounting for salinity,
% - kinetic fract C_k as following Merlivat (1978)

%------------------------------------------------------------------------------------
% Notes on saturation water vapor pressure. The formula used in this
% algorithm is based on the Engineering Toolbox one
% (http://www.engineeringtoolbox.com/water-vapor-saturation-pressure-air-d_689.html).
% Reference for this formula can be taken from:
% - Vladilo et al. 2013, "The habitable zone of Earth-like planets with different
%   levels of atmospheric pressure"
% - Rao et al. 2008, "Convective condensation of vapor in the presence of a
%   non-condensable gas of high concentration in laminar flow in a vertical pipe"
% - Saraireh & Thorpe 2010, "Modelling of Heat and Mass Transfer Involvin
%   vapour ..."
%------------------------------------------------------------------------------------

% % MJ79 VERSION
function CG_dE_18_MJ79 = CG_dE_18_MJ79(RH, T_air, d18_atmos, T_water, d18_water, kVal)

    ImplActivity = 'yes'; % 'yes'/'no' if you want to implement activity correction of humidity, default is yes
    % For activity chose one of the following values:
    % 0.98 Sea water
    % 1 Freshwater
    activity    = 0.98;
    % activity = 1;
    
    % Convert all permil values into decimal values
    d18_atmos   = d18_atmos/1000;
    d18_water   = d18_water/1000;
    kVal        = kVal/1000;
    
    % normalize humidity at the surface 
    switch ImplActivity
        case 'yes' % Following Horita et al., 2008
            e_a = (RH./100)*(exp(77.3450 + 0.0057 .* (T_air) - (7235./(T_air)))./(T_air).^8.2);
            e_s = (exp(77.3450 + 0.0057 .* (T_water) - (7235./(T_water)))./(T_water).^8.2);
            one_minus_h = (activity*e_s-e_a)/(activity*e_s);
            h = 1 - (activity*e_s-e_a)/(activity*e_s);
        case 'no'
            e_a = (RH./100)*(exp(77.3450 + 0.0057 .* (T_air) - (7235./(T_air)))./(T_air).^8.2);
            e_s = (exp(77.3450 + 0.0057 .* (T_water) - (7235./(T_water)))./(T_water).^8.2);
            one_minus_h = 1 - e_a/e_s;
            h = e_a/e_s;
        otherwise
    end
    
    alpha18_VL = 1./alpha18_LV(T_water);

    CG_dE_18_MJ79 = (((1-kVal).*((1+d18_water).*alpha18_VL - h.*(1+d18_atmos))./one_minus_h)-1).*1e3;
    
end