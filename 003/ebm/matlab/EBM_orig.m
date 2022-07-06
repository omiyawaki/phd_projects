function [To, lat, qo, mo, trn, lath, nl] = EBM_orig(SW,a,b,D,hs)

% This function is intended as an easily runnable form of (basically) the
% Hwang and Frierson EBM supplied by Tiffany in a 7/8/17 email.  Some
% physical constants differ slightly from their Hwang and Frierson values
% and others are now set as inputs to the function, but the model structure
% and solution algorithm are exactly the same.
%
% REQUIRED INPUT:
% SW    Specified components of the MSE generation rate--i.e., everything
%           other than clear sky OLR.  One way to break this down is as
%           (net downward shortwave at the TOA + longwave CRF - net
%           downward energy flux at the surface).  Units are W/m^2, must be
%           a column vector.  The specified values are assumed to be
%           defined at the centers of gridboxes of uniform meridional size
%           and the first element is assumed to be closest to the south
%           pole.
%
% OPTIONAL INPUTS:
% a     Coefficient in aT-b, the linearized formula for clear sky OLR.  A
%           scalar, default value is 2.07 W/m^2/K.
% b     Coefficient in aT-b, a scalar.  Default value is 332.4 W/m^2.
% D     Diffusivity, a scalar.  Default value is 1.06e6 m^2/s.
% hs    Relative humidity, a scalar.  Default is 0.8 (nondimensional).
%
% OUTPUTS:
% To    Final temperature profile, K
% lat   Assumed latitudes, deg
% qo    Final specific humidity profile, dimensionless
% mo    Final MSE profile, J/kg
% trn   Energy transport on the half latitudes, W
% lath  Assumed half latitudes, deg
% nl    Number of iterations required for model to converge

% 1) SET DEFAULT VALUES
if ~exist('a','var')
    a = 2.07;
end
if ~exist('b','var')
    b = 332.4;
end
if ~exist('D','var')
    D = 1.06e6;
end
if ~exist('hs','var')
    hs = 0.8;
end

% 2) SET NON-ADJUSTABLE PARAMETERS
lat_number = numel(SW);
lat_grid = pi/lat_number; %Latitude grid spacing in radians
lat      = -0.5*(lat_grid+pi)+lat_grid*(1:lat_number)';%1:lat_number makes 
    %a vector of the correct length and spacing 1.  Multiplying by lat_grid
    %changes spacing to lat_grid.  Subtracting pi/2 starts the array one
    %gridbox north of the south pole.  Subtracting lat_grid/2 starts the
    %array half a gridbox north of the south pole, as needed.
To    = 200*ones(lat_number,1); % initial condition for temperature, K 
del_t = 10^-4; %time step, will assume it's years (in other words,
    %0.876 days)
Cp    = todd_constants('cpair');%heat capacity of air at constant p, J/kg/K
r     = todd_constants('a'); %Earth radius, m
Lv    = todd_constants('Lv'); %Latent heat of vaporization, J/kg
Rd    = todd_constants('R'); %Gas constant for dry air, J/kg/K
Rv    = todd_constants('Rv'); %Gas constant of water vapor, J/kg/K
g     = todd_constants('g'); %Gravity, nominal, m/s^2
Ps    = 1e4*g; %Surface pressure, Pa
Do    = Ps*D/g; %Normalized diffusion coefficient, kg/s
Ta    = todd_constants('Tref'); %C-C reference temp, K
eo    = todd_constants('esref'); %C-C reference vapor pressure, Pa
C     = 1; %ocean heat capacity [Actually this parameter doesn't
    %really have physical meaning except in the case where we are emulating
    %a slab ocean model with no heat transport.  However this doesn't
    %matter if we are only interested in the model's steady-state behavior.
    %Given our assumption about the units of del_t, this has units of
    %W yr/m^2/K. 1 W yr/m^2/K = 3.1536e7 J/m^2/K.
A = zeros(lat_number,lat_number); %Square array w/ each dimension having
    %the length of the latitude vector.  Ultimately we wind up with an
    %array shaped like this:
for k=1:lat_number-1
    %iterate over all but last latitude
    A(k,k)   = -1;
    A(k,k+1) = 1;
end
B = zeros(lat_number,lat_number); %Same dimensions as A
for k=2:lat_number
    %iterate over all but the first latitude
    B(k,k)   = 1;
    B(k,k-1) = -1;
end
lat_ph = lat+lat_grid/2; %For now let's assume that lat_grid is the spacing
lat_nh = lat-lat_grid/2; %of the latitude grid points in radians and that
    %this is uniform.  Thus lat_ph (lat_nh) can be considered a vector of
    %northern (southern) edges of grid cells (assuming latitude is in
    %radians N).

% 3) RUN EBM
for n = 1:1000000
    %iterate over a very large number of steps
    es  = eo*exp((-Lv/Rv)*(1./To-1/Ta)); %Compute saturation vapor pressure
        %given current temperature To
    mo  = Cp*To+Lv*hs*Rd*es/(Rv*Ps); %Moist static energy cp*T+L*q,
        %assuming that the specific humidity and mixing ratios are
        %approximately equal and that the surface pressure is created
        %entirely by dry air.  z = 0 because we are evaluating this at
        %ground level.
    OLR = a*To-b; %outgoing longwave, W/m^2

    % general form of temperature tendency due to diffusion (horizontal
    % transport)
    Diff = Do*r^(-2).*(cos(lat)).^(-1).*(lat_grid)^(-2).*...
        (cos(lat_ph).*(A*mo)-cos(lat_nh).*(B*mo));

    % for south boundary (first grid point)
    Diff(1) = Do*r^(-2)*(cos(lat(1)))^(-1)*(lat_grid)^(-2)*...
        (cos(lat(2))*(mo(3)-mo(2))-cos(lat(1))*(mo(2)-mo(1)));

    % for north boundary (last grid point)
    Diff(lat_number) = Do*r^(-2)*(cos(lat(lat_number)))^(-1)*...
        (lat_grid)^(-2)*...
        (cos(lat(lat_number))*(mo(lat_number)-mo(lat_number-1))-...
        cos(lat(lat_number-1))*(mo(lat_number-1)-mo(lat_number-2)));

    Tn = To+(del_t*(SW-OLR+Diff))/C; %Technically this functional form is
        %appropriate only to the case where column heat capacity is
        %dominated by the ocean in the slab ocean/no ocean heat transport
        %case.  If this is not so (e.g., when seeking to emulate an AGCM),
        %the moist static energy capacity dm/dT that relates the moist
        %static energy change rate dm/dt to the temperature change rate
        %dT/dt has a well-defined (albeit in principle
        %temperature-dependent) value which will not generally equal the
        %specified C.  However, heat capacity is not relevant to the
        %equilibrium behavior of the model and thus this equation can be
        %considered a heuristic adjustment that moves the system closer to
        %equilibrium--warmer if MSE is being generated/converged into the
        %latitude band, cooler if MSE is being destroyed/diverged from the
        %latitude band.
    To = Tn; %updating surface temperature

    % break the loop if equilibrium is reached
    W = (max(abs(SW-OLR+Diff))<0.001 && n>1); %Set to 1 if the magnitude of
        %the MSE change rate at any given latitude is <0.001 W/m^2 and this
        %is not the first passage through the loop, 0 otherwise.
    if W==1
        break; %Stop iterating if convergence has been achieved
    end
end

% 4) CREATE OUTPUT VARIABLES
lat  = (180/pi)*lat; %Unit change to deg
es   = eo*exp((-Lv/Rv)*(1./To-1/Ta)); %Compute saturation vapor pressure
        %given current temperature To
qo   = hs*Rd*es/(Rv*Ps); %Compute specific humidity
mo   = Cp*To+Lv*qo; %Compute MSE
lath = 0.5*(lat(1:(end-1))+lat(2:end)); %Half lats at which transport is
    %defined
trn  = -2*pi*Do*cosd(lath).*(diff(mo,1,1)/lat_grid); %Compute transports on
    %the half latitudes, units of W
nl   = n;
