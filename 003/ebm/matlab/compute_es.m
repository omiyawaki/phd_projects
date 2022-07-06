function es = compute_es(T)

% This function computes the saturation vapor pressure (es, units of Pa) as
% a function of temperature (T, units of K).  The formula used is equation
% (1) of Frierson et al. 2007, which should be standard.

% 1) LOAD CONSTANTS
esref = todd_constants('esref'); %Reference saturation vapor pressure, Pa
Lv    = todd_constants('Lv'); %Latent heat of vaporization, J/kg
Rv    = todd_constants('Rv'); %Gas constant of water vapor, J/kg/K
Tref  = todd_constants('Tref'); %Reference temperature, K

% 2) COMPUTE SATURATION VAPOR PRESSURE
es    = esref*exp(-(Lv/Rv)*(bsxfun(@rdivide,1,T)-(1/Tref)));
