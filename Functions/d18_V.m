% Evaluate d18O from rain to vapour IN DELTA PER MIL, T in kelvin
function d18_V = d18_V(d18L, T)
    d18_V = ((1000+d18L)./alpha18_LV(T))-1000;
end