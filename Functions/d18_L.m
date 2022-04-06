% Evaluate d18O from vapor to liquid water IN DELTA PER MIL
% T in kelvin
function d18_L = d18_L(d18_V, T)
% Kendall pg59
    d18_L = ((d18_V + 1000).*alpha18_LV(T) - 1000);
end