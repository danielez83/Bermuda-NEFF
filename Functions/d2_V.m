% Evaluate dD from rain to vapour IN DELTA PER MIL, T in kelvin
function d2_V = d2_V(d2L, T)
    % Kendall pg59
    d2_V = ((1000+d2L)./alpha2_LV(T))-1000;
end