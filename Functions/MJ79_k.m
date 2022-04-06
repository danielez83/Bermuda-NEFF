  
% Kinetic fractionation factor for Global Closure Assumption (MJ79)
% result = k [‰]
% u_star, Friction velocity [cm/s]
% z, height [cm]
% isotope, 2 for Deuterium, 18 for Oxygen-18
% regime, 'smooth' or 'rough'

function MJ79_k = MJ79_k(u_star, z, isotope, regime)

    % Assuming wind measured at 10 m
    % u_star = 0:50; % cm/s
    % z = 10; % cm ?

    % Calculate of surface roughness
    z_0 = (u_star.^2) ./ (81.1.*981); %  acceleration of gravity [cm/s^2]

    % Constants needed
    chi = 0.4; % Von Karman constant, adimensional
    nu =  0.1568; % cm^2/s (1.568e-5 m^2/s) Kinematic viscosity of dry air at 300 K (http://www.engineeringtoolbox.com/dry-air-properties-d_973.html)
    D = 0.282; % cm^2/s molecular diffusion of water in air at 298 K

    % Molecular diffusions
    switch isotope
        case 2
            e_D = 0.0251; % 25.1‰ for Deuterioum (following Merlivat and Jouzel 1979), Di/D = 0.9755 (Merlivat 1978 --> [1/(Di/D)] - 1 = 1.0251-1 = 25.1‰
        case 18
            e_D = 0.0285; % 28.5‰ for Oxygen 18 (following Merlivat and Jouzel 1979), Di/D = 0.9723 (Merlivat 1978 --> [1/(Di/D)] - 1 = 1.0285-1 = 25.1‰
    end

    % Compute k by wind regime
    switch regime
        case 'smooth' % smooth regime
            n = 2/3;
            R_e = u_star.*z_0./nu; % Reynold number...not needed here but calculated for completness
            r_T_over_r_M = ((1/chi) * log((u_star .* z)/(30 .* nu)))/(13.6 .* (nu ./ D).^(2/3));
            MJ79_k = 1000*((1 + e_D).^n - 1)./((1+e_D).^n + r_T_over_r_M);
            %MJ79_k = k_smooth;
        case 'rough' % rough regime
            n = 1/2;
            R_e = u_star.*z_0/nu; % Reynold number
            r_T_over_r_M = (((1/chi) .* log(z ./ z_0)) - 5) ./ (7.3 .* (R_e.^(1./4)) * (nu ./ D).^(1/2));
            MJ79_k = 1000.*((1 + e_D).^n - 1)./((1+e_D).^n + r_T_over_r_M);
        otherwise
                disp('No wind regime selected');
    end
end