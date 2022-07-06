function [dm_dy, lath, m, latf, H, Nim, Nout, nc] = EBM_steady(N,Do,m_i)

% This function is an energy balance model very similar to that used by
% Hwang and Frierson 2010.  Key differences are the following:
%
% 1) The model explicitly computes a steady-state solution, rather than
%        iterating until convergence is achieved.
% 2) The net MSE generation rate with which the MSE flux divergence is in
%        balance is directly specified in its entirety, rather than being a
%        function of temperature.  Thus there is no need for explicit
%        reference to temperature or any other specific subcomponent of
%        moist static energy anywhere in the code.  However, to get the
%        assumed steady state the global mean MSE generation rate must be
%        zero--this is imposed by making a globally uniform adjustment to
%        the input N before using it to drive the model.  This also implies
%        that the model's global mean MSE will remain fixed at whatever
%        value it was initialized at.  Since global mean MSE is an input
%        parameter of this model, rather than an output, the MSE gradient
%        is arguably more physically meaningful and the input global mean
%        MSE must be chosen with care if the full MSE profile is to be
%        physically interpreted.
% 3) Diffusivity is allowed to vary with latitude.
% 4) Numerics are more explicitly finite-volume, in order to obtain the
%        desired conservation properties.
%
% INPUTS (REQUIRED)
% N     Zonal mean MSE generation rate in W/m^2.  Defined on N evenly-
%           spaced full latitudes stretching from pole to pole--given this
%           assumption, the grid boxes have meridional extent 180/N deg and
%           the most southerly box has its center at -90+90/N deg N.
% Do    Zonal mean diffusivity in kg/s.  **This is the EBM definition of
%           diffusivity, not the D = -H/G diagnosed diffusivity.  The two
%           diffusivities are related by Do = D/cos(lat).**  Must be
%           specified on N-1 half latitudes evenly spaced between the N
%           full latitudes.  Can alternatively just pass a scalar, which
%           will be automatically expanded.
%
% INPUT (OPTIONAL)
% m_i   Initial MSE profile, units J/kg.  Can be defined on the N full
%           latitudes, or by passing a scalar to specify the global mean
%           MSE.  Default value is 0.
%
% OUTPUTS
% dm_dy MSE gradient, units J/kg/m.  Defined on the half latitudes.
% lath  Half latitudes, deg N.
% m     MSE, units J/kg.  Defined on the full latitudes.
% latf  Full latitudes, deg N.
% H     MSE transport, W.  Defined on the half latitudes.
% Nim   Global mean MSE generation rate associated with the input N
%           profile, units W/m^2.  Scalar.
% Nout  MSE generation rate profile with the global mean removed to enable
%           appropriate forcing of the model.
% nc    Degree of non-closure in transport calculations.  Two-element row
%           vector, left (right) element is implied flux at north (south)
%           pole from south-to-north (north-to-south) transport integral.
%           Should be very small due to explicit adjustment of N to have
%           zero global mean.  Units are W/m^2.

% 1) VALIDATE INPUTS
if ~(isequal(numel(N),numel(Do)+1) || isequal(numel(Do),1))
    error('Do has incorrect length')
elseif ~exist('m_i','var')
    m_i = 0;
elseif ~(isequal(numel(N),numel(m_i)) || isequal(numel(m_i),1))
    error('m_i has incorrect length')
end

% 2) REMOVE GLOBAL MEAN MSE GENERATION RATE
nlats = numel(N); %Number of full latitudes
lath  = linspace(-90,90,nlats+1)'; %Create half latitudes
latf  = 0.5*(lath(1:(end-1))+lath(2:end)); %Create full latitudes
dsinl = sind(lath(2:end))-sind(lath(1:(end-1)));
Nim   = 0.5*sum(N(:).*dsinl(:)); %Global mean MSE generation rate, W/m^2
Nout  = N-Nim; %Now has zero global mean for further use

% 3) COMPUTE MSE TRANSPORTS
a      = todd_constants('a'); %Earth radius, m
diff_H = 2*pi*a*a*(Nout(:).*dsinl(:)); %Difference between meridional MSE
    %transports at successive half latitudes, W
s2n    = [0; cumsum(diff_H,'forward')]; %South to north
n2s    = [-cumsum(diff_H,'reverse'); 0];%North to south
nc     = [s2n(end) n2s(1)]/(4*pi*a*a); %Store degree of non-closure
H      = 0.5*(s2n+n2s); %Average across both integration directions, W
H      = H(2:(end-1)); %Do not return transports at the poles

% 4) COMPUTE MSE GRADIENTS
lath   = lath(2:(end-1)); %Drop poles
dm_dy  = -H(:)./(2*pi*a*(Do(:).*cosd(lath(:)))); %Employs definition of
    %diffusivity to compute the gradient, J/kg/m
    
% 5) COMPUTE MSE PROFILE
dm     = ((pi*a)/nlats)*dm_dy; %MSE increment between successive boxes
m      = [0; cumsum(dm,'forward')]; %MSE increment relative to the value in
    %the most southerly grid box, J/kg
mm     = 0.5*sum(m(:).*dsinl(:)); %Global mean of the MSE increment profile
m_im   = 0.5*sum(m_i(:).*dsinl(:)); %Global mean of initial MSE profile
m      = m-mm+m_im; %Adjust MSE increment profile so that it has the same
    %mean as the initial MSE profile, i.e. could have been produced by
    %rearranging the initial condition, J/kg.
