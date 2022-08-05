# Bermuda-NEFF
The code can be used to reproduce the results in Zannoni et al. (2022), currently under review at Journal of Geophysical Research: Atmospheres.
A preprint of the manuscript is available here: https://doi.org/10.1002/essoar.10511947.1


Description the Matlab scripts.

## LoadAndFilter
Load and filter all the data required to run the code for estimating the isotopic composition of evaporation flux and non-equilibrium fractionation.

## ImportBermudaDEM
Plot the study area (Figure 1 in the manuscript). Some external data is needed to reproduce the plot (links provided as comments).

## DisplayTimeSeries_R1.m
Show timeseries od $\delta^{18}$O, d-excess, wind speed, SST and $h_{s}$ (Figure 2 in the manuscript).

## EstimateErrorOnFluxComposition_V4
Run this script to calculate the isotopic composition of the evaporation flux and distribution of non-equilibrium fractionation factors (Figure 3 in the manuscript). 
* By default the random resampling is set to 100 repetition. to reproduce the result of the manuscript keep the same seed and change the 
*nIterations* from 10e2 to 10e3.
* By default the script calculate the isotopic composition of the flux with Keeling Plot (KP) and with Flux Gradient (FG) methods.

## K_WS_dependency_V4
This script calculates the non-equilibrium fractionation factors for each pair of water vapor observations (top and bottom inlet) with the Keeling Plot. Afterwards, fractionation factors are binned by wind speed classes in order to fit a linear model between k18 and wind speed (Figure 4 in the manuscript). 
Run *K18_WS_StatTestsOnResiduals* to perform statistical tests on residuals (MJ79 parametrization vs observations).

## PlotSST_Salinity_timeseries
Run this script to reproduce the time series of SST and Salinity shown in Figure 5 in the manuscript.

## SSTandSalinityImpact
Run this script to reproduce the distributions of SST and Salinity difference (converted to isotopic composition of water vapor and ocean water), as shown in Figure 5 in the manuscript.

## Test_kvalues_over_ocean
Run this script to reproduce Table 4 in the Manuscript and Figure S7 in supplementary material. External datasets are required to run the script. Link provided as comment followiin the pubblication of Benetti et al. (2017).

