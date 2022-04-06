%Fractionation factor of 18O Ice or Liquid/Vapor
% Ice-Liquid point at 268.15K (-5°C) 
% T [k]
% See http://www.hydrochemistry.eu/exmpls/istp.html
function alpha18_LV = alpha18_LV(T)
    %if T >= 268.15 %Majoube 1971, Temperature in Kelvin
%         alpha18_LV = exp((1137./(T.^2))-(0.4156./T)-0.0020667);
        % Horita and Wesolowski 1994
        alpha18_LV = exp((-7.685 + 6.7123.*((1e3)./T) - 1.6664.*((1e6)./T.^2) + 0.35041.*((1e9)./T.^3))./1e3);
    %else 
    %    alpha18_LV = exp((1137./(T.^2))-(0.4156./T)-0.0020667); % do not differentiate
        %alpha18_LV = exp((8312.58./(T.^2))-(49.192./T)+0.0831); % Ellehoj et al., 2013
    %end
end

% function alpha18_LV = alpha18_LV(T) % Horita and Wesolowksiy 1994, V/L. Temperature in Kelvin
%     alpha18_LV = exp((0.35041e-6*T^(-3))-(1.6664e3*T^(-2))+(6.7123*T^(-1))-7.685e-3);
% end