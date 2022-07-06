function cvalue = todd_constants(cname)

% Return the value of a specified physical constant for use in another
% script or function, eliminating the need to define values of the
% constants in each new script and thus making it easier for analyses to be
% physically consistent across scripts.  Values will not necessarily be the
% same as those from Tiffany's physics.mat file.

cpair = 1004; %Heat capacity of dry air, WH p. 467, J/kg/K, 1/13/17
R     = 287; %Gas constant of dry air, WH p. 467, J/kg/K, 1/13/17
Rv    = 461; %Gas constant of water vapor, WH p. 467, J/kg/K, 3/23/17
Lv    = 2.5008e6; %Latent heat of vaporization, probably actually used in
    %ECHAM6, J/kg, 1/13/17
Lf    = 3.337e5; %Latent heat of fusion, probably actually used in ECHAM6,
    %J/kg, 1/13/17
g     = 9.81; %Gravity, nominal, m/s^2, 1/13/17
a     = 6371000; %Earth radius, nominal, m, 1/13/17
Tref  = 273.16; %Reference temperaure for Clausius-Clapeyron, Frierson et
    %al. 2007 p. 1682, K, 3/23/17
esref = 610.78; %Saturation vapor pressure at Tref for Clausius-Clapeyron,
    %Frierson et al. 2007 p. 1682, Pa, 3/23/17

cvalue = eval(cname);