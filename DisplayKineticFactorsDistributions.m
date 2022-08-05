%% Configuration
Show_FG             = 0; % 1 ON, 0 OFF


%%
% Molecular diffusivities from Merlivat 1978
Di1816 = 0.9727; % From Horita 2008 is Di1816 = 0.9727;
Di21 = 0.9757; % From Horita 2008 is Di21   = 0.9757;
% Di1816 = 0.9943;
% Di21 = 0.9880;
D18_D16 = (1 - Di1816) * 1000; % [‰]
D2_D1 = (1 - Di21) * 1000; % [‰]


if ~exist('kinetic_18_16_sample_FG', 'var')
load('../Matlab data/SImulations/SImulationOutputAVG5min_VaporOffset_30032022.mat') % Fixed Ocean composition, SST from OSTIA, Isotopic data recalibrated using 5min averaging window instead of 10.5 (default in HC data)
% load('../Matlabdata/SImulationsSImulationOutputAVG5min_VaporOffset_VaryinSST_18072022.mat') % Same as above but using day-by-day SST difference between OSTIA and St. George (not average offset)
%load('../Matlab data/BottomTopData_SYNC_HF_5minAVG_1e3_repetitions.mat')    % FIxed Ocean composition, SST from OSTIA, Isotopic data recalibrated using 7.5min averaging window instead of 10.5 (default in HC data)
    disp('Load saved data')
    pause(1)
else
    disp('Use current values')
    pause(1)
end


x_values_18_16 = -30:.01:30;
x_values_2_1 = -30:.01:30;
% Flux gradient
pd = fitdist(kinetic_18_16_sample_FG,'Kernel', 'Kernel','normal');
FGPDF18 = pdf(pd,x_values_18_16);
FGk18 = kinetic_18_16_sample_FG;
pd = fitdist(kinetic_2_1_sample_FG,'Kernel', 'Kernel','normal');
FGPDF2 = pdf(pd,x_values_2_1);
FGk2   = kinetic_2_1_sample_FG;

% Keeling Plot
pd = fitdist(kinetic_18_16_sample_KP,'Kernel', 'Kernel','normal');
KPPDF18 = pdf(pd,x_values_18_16);
KPk18 = kinetic_18_16_sample_KP;
pd = fitdist(kinetic_2_1_sample_KP,'Kernel', 'Kernel','normal');
KPPDF2 = pdf(pd,x_values_2_1);
KPk2   = kinetic_2_1_sample_KP;

plot_lines = 1;
plot_MJ79_full_area = 1;
%% Plot distributions
% https://se.mathworks.com/matlabcentral/fileexchange/42905-break-x-axis
clf
subplot(2,1,1)
    if Show_FG == 1
        h1 = histogram(FGk18);
        h1.Normalization = 'pdf';
        h1.NumBins = 20;
        h1.EdgeColor = [0 0 0];
        % h1.FaceColor = [0.9294 0.6941 0.1255];
        h1.FaceColor = [1 0 0];
        h1.EdgeAlpha = 1;
        h1.FaceAlpha = .1;
    end
    hold on
    h2 = histogram(KPk18);
    h2.Normalization = 'pdf';
    h2.NumBins = 20;
    h2.EdgeColor = [0 0 0];
    % h1.FaceColor = [0.9294 0.6941 0.1255];
    h2.FaceColor = [0.5 0.5 .5];
    h2.EdgeAlpha = 1;
    h2.FaceAlpha = .1;
    
    % plot(x_values, y, 'LineWidth',2, 'Color', [0.9294 0.6941 0.1255])
    if Show_FG == 1
        plot(x_values_18_16, FGPDF18, 'LineWidth', 2, 'Color', [1 0.0 0.0])
    end
    plot(x_values_18_16, KPPDF18, 'LineWidth', 2, 'Color', [0 0 0])
    % Graphics
    grid on
    box on
    xlim([0 29])
    ylim([0 0.8])
    ax = gca;
    ax.YAxis.FontSize = 13;
    ax.XAxis.FontSize = 13;
    ax.XAxis.FontWeight = 'bold';
    ylabel('PDF')
    %xlabel('\epsilon_{k-18} (‰)')
    xlabel('k_{18} (‰)')
    set(gcf, 'Color', [1 1 1]);
    max_y = max(FGPDF18);
    if plot_lines == 1
        % Plot diffusivity value from literature (Horita 2008)
        % H. Craig, L.I. Gordon, Y. Horibe. Isotopic exchange effects in the evaporation of water. 1. 
        % Low-temperature experimetnal results. J. Geophys. Res., 68, 5079–5087 (1963).
        %line([5.7 5.7], [0 max(max_y)+0.1], 'LineStyle', '--', 'LineWidth', 2, 'Color', [0 0 0])

        %﻿L. Merlivat. Molecular diffusivities of (H2O)-O-16, HD16O, and (H2O)-O-18 in 
        % gases. J. Chem. Phys., 69, 2864–2871 (1978).
        line([D18_D16 D18_D16], [0 max(max_y)+0.1], 'LineStyle', '--', 'LineWidth', 2,  'Color', [0 0 0])

        % C.D. Cappa, M.B. Hendricks, D.J. DePaolo, R.C. Cohen. Isotopic fractionation of 
        % water during evaporation. J. Geophys. Res. Atmos., 108, 4525–4535 (2003).
        % @ 20°C
        %line([11.8 11.8], [0 max(max_y)+0.1], 'LineStyle', '--', 'LineWidth', 2,  'Color', [0 0 0])

        % Pfahl and Wenli 2009, ﻿Lagrangian simulations of stable isotopes in water vapor: An evaluation of 
        % nonequilibrium fractionation in the Craig-Gordon model
        % Good for Mediterranean Sea
        line([7.5 7.5], [0 max(max_y)+0.1], 'LineStyle', '--', 'LineWidth', 2,  'Color', [0 0 0])

        % Merlivat and Jouzel 1979, ﻿Global Climatic Interpretation of the Deuterium-Oxygen 18 Relationship for Precipitation
        % For smooth regime, reported in Benetti et al., 2018 ﻿A Framework to Study Mixing Processes in 
        % the Marine Boundary Layer Using Water Vapor Isotope Measurements
        % line([6 6], [0 max(max_y)+0.01], 'LineStyle', '--', 'LineWidth', 2,  'Color', [0 0 0])

        % Uemura et al 2010,﻿Triple isotope composition of oxygen in atmospheric water vapor
        line([(1-(1/1.0083))*1e3 (1-(1/1.0083))*1e3], [0 max(max_y)+0.1], 'LineStyle', '--', 'LineWidth', 2,  'Color', [0 0 0])
    end
    
    if plot_MJ79_full_area
        %line([2.59 2.59], [0 max(max_y)+1], 'LineStyle', '-', 'LineWidth', .5,  'Color', [0.7 0.7 0.7])
        line([6.6453 6.6453], [0 max(max_y)+1], 'LineStyle', '-', 'LineWidth', .5,  'Color', [0.7 0.7 0.7])
        line([5.6741 5.6741], [0 max(max_y)+1], 'LineStyle', '-', 'LineWidth', .5,  'Color', [0.7 0.7 0.7])
        line([2.5937 2.5937], [0 max(max_y)+1], 'LineStyle', '-', 'LineWidth', .5,  'Color', [0.7 0.7 0.7])
        line([4.4985 4.4985], [0 max(max_y)+1], 'LineStyle', '-', 'LineWidth', .5,  'Color', [0.7 0.7 0.7])
        %line([8.52 8.52], [0 max(max_y)+1], 'LineStyle', '-', 'LineWidth', .5,  'Color', [0.7 0.7 0.7])
        line([D18_D16 D18_D16], [0 max(max_y)+0.1], 'LineStyle', '--', 'LineWidth', 2,  'Color', [0 0 0])
    end

    % text(D18_D16, max(max_y) + 0.01, sprintf('%.1f‰', D18_D16), 'Color', [0 0 0],'FontSize',14);
    hold off
    ax = gca;
    ax.XAxis.FontWeight = 'bold';
    ax.YAxis.FontWeight = 'bold';
    
    h = breakxaxis([10 27], 0.01);
% https://se.mathworks.com/matlabcentral/answers/159417-how-to-calculate-the-confidence-interval  
    
subplot(2,1,2)
    if Show_FG == 1
        h1 = histogram(FGk2);
        h1.Normalization = 'pdf';
        h1.NumBins = 20;
        h1.EdgeColor = [0 0 0];
        % h1.FaceColor = [0.9294 0.6941 0.1255];
        h1.FaceColor = [1 0 0];
        h1.EdgeAlpha = 1;
        h1.FaceAlpha = .1;
    end
    hold on
    h2 = histogram(KPk2);
    h2.Normalization = 'pdf';
    h2.NumBins = 20;
    h2.EdgeColor = [0 0 0];
    % h1.FaceColor = [0.9294 0.6941 0.1255];
    h2.FaceColor = [0.5 0.5 .5];
    h2.EdgeAlpha = 1;
    h2.FaceAlpha = .1;
    % plot(x_values, y, 'LineWidth',2, 'Color', [0.9294 0.6941 0.1255])
    if Show_FG == 1
        plot(x_values_2_1, FGPDF2, 'LineWidth', 2, 'Color', [1 0.0 0.0])
    end
    plot(x_values_2_1, KPPDF2, 'LineWidth', 2, 'Color', [0 0 0])
    ylabel('PDF')
    % Graphics
    max_y = max(FGPDF2);
    % h2.FaceColor = [1.0000 0.4118 0.1608];
    % Expected kinetic fractionation values with same thetaN estimated for
    % d18O
    if plot_lines == 1
        % Plot diffusivity value from literature (Horita 2008)
        % H. Craig, L.I. Gordon, Y. Horibe. Isotopic exchange effects in the evaporation of water. 1. 
        % Low-temperature experimetnal results. J. Geophys. Res., 68, 5079–5087 (1963).
        %line([12 12], [0 max(max_y)+0.02], 'LineStyle', '--', 'LineWidth', 2, 'Color', [0 0 0])

        %﻿L. Merlivat. Molecular diffusivities of (H2O)-O-16, HD16O, and (H2O)-O-18 in 
        % gases. J. Chem. Phys., 69, 2864–2871 (1978).
        line([D2_D1 D2_D1], [0 max(max_y)+0.02], 'LineStyle', '--', 'LineWidth', 2,  'Color', [0 0 0])

        % C.D. Cappa, M.B. Hendricks, D.J. DePaolo, R.C. Cohen. Isotopic fractionation of 
        % water during evaporation. J. Geophys. Res. Atmos., 108, 4525–4535 (2003).
        % @ 20°C
        %line([13.1 13.1], [0 max(max_y)+0.02], 'LineStyle', '--', 'LineWidth', 2,  'Color', [0 0 0])

        % Pfahl and Wenli 2009, ﻿Lagrangian simulations of stable isotopes in water vapor: An evaluation of 
        % nonequilibrium fractionation in the Craig-Gordon model
        % Good for Mediterranean Sea
        line([3.9 3.9], [0 max(max_y)+0.02], 'LineStyle', '--', 'LineWidth', 2,  'Color', [0 0 0])

        % Merlivat and Jouzel 1979, ﻿Global Climatic Interpretation of the Deuterium-Oxygen 18 Relationship for Precipitation
        % For smooth regime, reported in Benetti et al., 2018 ﻿A Framework to Study Mixing Processes in 
        % the Marine Boundary Layer Using Water Vapor Isotope Measurements
        %line([5.3 5.3], [0 max(max_y)+0.01], 'LineStyle', '--', 'LineWidth', 2,  'Color', [0 0 0])
    end
    
    if plot_MJ79_full_area
        %line([2.59 2.59], [0 max(max_y)+1], 'LineStyle', '-', 'LineWidth', .5,  'Color', [0.7 0.7 0.7])
        line([6.6453 6.6453].*.88, [0 max(max_y)+1], 'LineStyle', '-', 'LineWidth', .5,  'Color', [0.7 0.7 0.7])
        line([5.6741 5.6741].*.88, [0 max(max_y)+1], 'LineStyle', '-', 'LineWidth', .5,  'Color', [0.7 0.7 0.7])
        line([2.5937 2.5937].*.88, [0 max(max_y)+1], 'LineStyle', '-', 'LineWidth', .5,  'Color', [0.7 0.7 0.7])
        line([4.4985 4.4985].*.88, [0 max(max_y)+1], 'LineStyle', '-', 'LineWidth', .5,  'Color', [0.7 0.7 0.7])
        %line([8.52 8.52], [0 max(max_y)+1], 'LineStyle', '-', 'LineWidth', .5,  'Color', [0.7 0.7 0.7])
        line([D2_D1 D2_D1], [0 max(max_y)+0.02], 'LineStyle', '--', 'LineWidth', 2,  'Color', [0 0 0])
    end
    % text(D18_D16, max(max_y) + 0.01, sprintf('%.1f‰', D18_D16), 'Color', [0 0 0],'FontSize',14);
    hold off
    % Graphics
    grid on
    box on
    xlim([-5 30])    
    ylim([0 0.15])
    ax = gca;
    ax.XAxis.FontWeight = 'bold';
    ax.YAxis.FontWeight = 'bold';
    hold off;
    xlabel('k_{2} (‰)')
    set(gcf, 'Color', [1 1 1]);
    h = breakxaxis([12 21], 0.01);
    ax = gca;
    ax.XAxis.FontWeight = 'bold';
% last configuration
set(gcf, 'Position', [465   354   480   443]);
    
%% Print statistics KP
mean18value = mean(kinetic_18_16_sample_KP);
std18value = std(kinetic_18_16_sample_KP);
mean2value = mean(kinetic_2_1_sample_KP);
std2value = std(kinetic_2_1_sample_KP);
% Find 95 CI
low18value = x_values_18_16(find(cumsum(KPPDF18)>=2.5, 1)); 
up18value = x_values_18_16(find(cumsum(KPPDF18)>=97.5, 1));  
low2value = x_values_2_1(find(cumsum(KPPDF2)>=2.5, 1));
up2value = x_values_2_1(find(cumsum(KPPDF2)>=97.5, 1)); 
clc
fprintf('KP stats\n')
fprintf('    Mean    SD      CI5     CI95\n')
fprintf('---------------------------------\n')
fprintf('k18 %.2f    %.2f    %.2f    %.2f\n', mean18value, std18value, low18value, up18value);
fprintf('k2  %.2f    %.2f    %.2f    %.2f\n', mean2value, std2value, low2value, up2value);
fprintf('\n')
fprintf('\n')
% Find 95 CI
%value = x_values_18_16(find(cumsum(FGPDF18)>=2.5, 1));
%value = x_values_18_16(find(cumsum(FGPDF18)>=97.5, 1)):  
%value = x_values_2_1(find(cumsum(FGPDF2)>=2.5, 1));
%value = x_values_2_1(find(cumsum(FGPDF2)>=97.5, 1)): 
%% Print statistics FG
mean18value = mean(kinetic_18_16_sample_FG);
std18value = std(kinetic_18_16_sample_FG);
mean2value = mean(kinetic_2_1_sample_FG);
std2value = std(kinetic_2_1_sample_FG);
% Find 95 CI
low18value = x_values_18_16(find(cumsum(FGPDF18)>=2.5, 1)); 
up18value = x_values_18_16(find(cumsum(FGPDF18)>=97.5, 1));  
low2value = x_values_2_1(find(cumsum(FGPDF2)>=2.5, 1));
up2value = x_values_2_1(find(cumsum(FGPDF2)>=97.5, 1)); 
fprintf('FG stats\n')
fprintf('    Mean    SD      CI5     CI95\n')
fprintf('---------------------------------\n')
fprintf('k18 %.2f    %.2f    %.2f    %.2f\n', mean18value, std18value, low18value, up18value);
fprintf('k2  %.2f    %.2f    %.2f    %.2f\n', mean2value, std2value, low2value, up2value);
fprintf('\n')
fprintf('\n')