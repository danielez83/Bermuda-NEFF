% Estimate isotopic ratio of atmospheric water vapor with
% the Craig-Gordon model and closure assumption
% Computation is made for sea water (activity 0.98), see Horita 2008 for
% additional informations
% RH: 1 - 100 [%]
% Ait temperature [k]
% Sea/Water temperature [k]
% Isotopic ratio of water (R)
% kinetic fractionation factoR (e.g. 6.1 ‰)
% Isotopic specie of water of interest (2, 18)

function MJ79_CGCA = MJ79_CGCA(RH, T_air, T_water, Rsw,  k, isotope)

    activity = 0.98;
    % normalize humidity at the surface 
    e_a = (RH/100)*(exp(77.3450 + 0.0057 * (T_air) - (7235/(T_air)))/(T_air)^8.2);
    e_s = (exp(77.3450 + 0.0057 * (T_water) - (7235/(T_water)))/(T_water)^8.2);
    one_minus_h = (activity*e_s-e_a)/(activity*e_s);
    h = 1 - (activity*e_s-e_a)/(activity*e_s);
    
    % Equilibrium fractionation factor
    switch isotope
        case 2
            alpha_VL = 1./alpha2_LV(T_water);
        case 18
            alpha_VL = 1./alpha18_LV(T_water);
    end
    
    a_k = (k*1e-3)+1;
    
    % Calculate Isotopic ratio
    MJ79_CGCA = (Rsw.*alpha_VL)./(h + a_k.*one_minus_h);
    
end