% Evaluate d18O from vapor to liquid water IN DELTA PER MIL
% T in kelvin
function d2_L = d2_L(d2_V, T)
    d2_L = ((d2_V + 1000).*alpha2_LV(T) - 1000);
end