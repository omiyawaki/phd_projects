function [To, qo, mo, H, lath, nl, dt_c, last_loop, full_iter, ...
    n_attempts] = EBM_flux(Ns,lat,a,b,Do,ptot,rh,zg,Too,dt_c,max_s,cnv_t)

% This function is an energy balance model very similar to that used by
% Hwang and Frierson 2010.  Key differences are the following:
%
% 1) Overall numerical algorithm has been carefully formulated in flux form
%        to improve the conservation properties of the model, in contrast
%        to HF2010's finite difference formulation.
% 2) The boundary condition (no flux through the poles) makes more sense
%        than the HF2010 boundary condition which I don't fully understand.
% 3) Diffusivity can vary with latitude and is input in units of kg/s.
% 4) Relative humidity can vary with latitude.
% 5) A more accurate formula for computing specific humidity from relative
%        humidity and temperature is used.
% 6) A geopotential height contribution to the mean state moist static
%        energy can be specified.
% 7) Default first guess temperature profile is 286.65+13.5*cosd(2*lat),
%        inspired by the Neale and Hoskins 2001 Qobs profile.  Can
%        alternatively be specified by the user as an input.
% 8) Automated selection of the timestep used in the iterative process,
%        based on making sure that temperature values never exceed the
%        0-400K range and restarting the iterative process with a smaller
%        timestep if they do.
%
% REQUIRED INPUTS:
% Ns    Specified components of the MSE generation rate--i.e., everything
%           other than clear sky OLR.  One way to break this down is as
%           (net downward shortwave at the TOA + longwave CRF - net
%           downward energy flux at the surface).  Units are W/m^2, must be
%           a column vector.  The specified values are assumed to be
%           defined on the full latitude grid.
% lat   Column vector of full latitudes, units deg N.
%
% OPTIONAL INPUTS:
% Passage of [] as the value of an optional input will result in the
% selection of its default value.
% a     Coefficient in aT-b, the linearized formula for clear sky OLR.  A
%           scalar, default value is 2.07 W/m^2/K.
% b     Coefficient in aT-b, a scalar.  Default value is 332.4 W/m^2.
% Do    Diffusivity, a scalar or column vector.  If the column vector has
%           the same number of elements as lat, Do is defined on the full
%           latitudes and will be interpolated to the half latitudes via a
%           spline.  If the column vector has one fewer element than lat,
%           it is assumed to be defined directly on the half latitudes.
%           Default value is 1.06e10 kg/s.  **This is the EBM definition of
%           diffusivity, not the D = -H/G diagnosed diffusivity.  The two
%           diffusivities are related by Do = D/cos(lat).**
% ptot  Total air pressure to use in the conversion from relative to
%           specific humidity.  Default value is 92500 Pa.
% rh    Relative humidity, a scalar or column vector defined on the full
%           latitudes.  Default value is 0.8.
% zg    Specified geopotential height contribution to the mean state moist
%           static energy.  A scalar or column vector defined on the full
%           latitudes, default value is 0 m.
% Too   First guess temperature profile used to initialize the iterative
%           process, default value is 286.65+13.5*cosd(2*lat) K.
% dt_c  First guess of the ratio of timestep to assumed system heat
%           capacity, which controls the temperature change made in
%           response to a given total heating rate.  If iteration does not
%           appear to be converging (i.e., departs from the 0-400K
%           temperature range), dt_c is halved and the iterative process is
%           restarted.  This process of halving is repeated until
%           convergence is achieved or the maximum number of steps for a
%           given value of dt_c is exceeded.  Default value is 2^(-11) K
%           m^2/W.
% max_s First guess of the maximum allowed number of iteration steps, is
%           doubled for each necessary halving of dt_c (i.e., the time
%           for which the integration is run doesn't change the steps are
%           just made smaller).  If max_s is reached before convergence to
%           within the target accuracy is achieved, the iteration process
%           will quit anyway.  Default value is 2^22.
% cnv_t Tolerance to within which net heating rate must converge everywhere
%           in the domain to consider the model to have reached a steady
%           state.  Default value is 0.001 W/m^2.
%
% OUTPUTS:
% To          Final temperature profile, K
% qo          Final specific humidity profile, dimensionless
% mo          Final MSE profile, J/kg
% H           Energy transport on the half latitudes, W
% lath        Assumed half latitudes, deg
% nl          Number of iterations required for model to converge, will be
%                 negative if iteration process stopped before convergence
%                 was achieved
% dt_c        Ratio of timestep to assumed system heat capacity actually
%                 used for successful iteration, K m^2/W
% last_loop   Length of time in seconds required to run the final
%                 (successful) iteration attempt
% full_iter   Length of time in seconds required to run all of the
%                 iteration attempts
% n_attempts  Number of iteration attempts made (i.e., number of dt_c/max_s
%                 combinations used before convergence was achieved or
%                 max_s steps were taken without achieving convergence)

% 1) SET DEFAULT VALUES
if (~exist('a','var')) || isempty(a)
    a  = 2.07; %W/m^2/K
end
if (~exist('b','var')) || isempty(b)
    b  = 332.4;%W/m^2
end
if (~exist('Do','var')) || isempty(Do)
    Do = 1.06e10*ones(numel(lat)-1,1); %kg/s, half latitudes
end
if (~exist('ptot','var')) || isempty(ptot)
    ptot = 92500; %Pa
end
if (~exist('rh','var')) || isempty(rh)
    rh = 0.8*ones(size(lat)); %dimensionless, full latitudes
end
if (~exist('zg','var')) || isempty(zg)
    zg = zeros(size(lat)); %m, full latitudes
end
if (~exist('Too','var')) || isempty(Too)
    Too = 286.65+13.5*cosd(2*lat);%27C at the equator and 0C at the poles--
        %same range as Neale and Hoskins 2001 Qobs, but much smoother
end
if (~exist('dt_c','var')) || isempty(dt_c)
    dt_c = 2^(-11); %Proportionality constant controlling the relationship
        %between the TOA energy imbalance and the adjustment to the
        %low-level temperature, units K m^2/W.  Larger value should speed
        %integration, unless it causes the model to become unstable or to
        %oscillate around the true solution.  If oscillatory behavior is
        %detected, reduce this value to enhance stability.  This is
        %basically a ratio of a timestep to a heat capacity.
end
if (~exist('max_s','var')) || isempty(max_s)
    max_s = 2^22; %Perform up to 2^22 iteration steps while seeking steady
        %solution to the model
end
if (~exist('cnv_t','var')) || isempty(cnv_t)
    cnv_t = 0.001; %Tolerance to within which net heating rate must
        %converge everywhere in the domain to consider the model to have
        %reached a steady state, units W/m^2.
end
g  = todd_constants('g'); %m/s^2
zg = g*zg; %Convert geopotential height units to J/kg

% 2) EXPAND DIFFUSIVITY TO VECTOR AS NEEDED
lath = 0.5*(lat(1:(end-1))+lat(2:end));
if lath(1)<0
    %south to north
    lath = [-90; lath; 90];
else
    lath = [90; lath; -90];
end
if isscalar(Do)
    %Diffusivity is globally uniform
    Do = Do*ones(numel(lat)-1,1); %Diffusivity should be defined on the N-1
        %half latitudes that are not the poles
elseif isequal(size(Do),size(lat))
    %Diffusivity is defined on N full latitudes
    Do = interp1(lat,Do,lath(2:(end-1)),'spline'); %Now on N-1 half lats
elseif ~isequal(numel(Do),numel(lat)-1)
    %Diffusivity is not defined on N-1 half latitudes
    error('Unexpected number of Do latitudes!')
end

% 3) INITIALIZE EBM VARIABLES
cp    = todd_constants('cpair');%heat capacity of air at constant p, J/kg/K
r     = todd_constants('a'); %Earth radius, m
Lv    = todd_constants('Lv'); %Latent heat of vaporization, J/kg
H     = zeros(size(lath)); %MSE transport divided by 2*pi on N+1 half lats
dlat  = (pi/180)*diff(lat,1,1); %Full latitude spacings on N-1 half lats
awt   = sind(lath(2:end))-sind(lath(1:(end-1))); %Differences of successive
    %sines of half latitudes--this quantity is defined on N full latitudes
    %and is proportional to the surface area of each full latitude strip
coslath = cosd(lath(2:(end-1))); %Cosines of N-1 half latitudes
H_pf    = -((coslath.*Do)./dlat);%Transform MSE differences to transports
    %divided by 2*pi, N-1 half lats
inv_dA  = 1./(r*r*awt); %This quantity is proportional to the inverse of
    %the surface area of each full latitude strip.  N full latitudes, units
    %m^{-2}.
    
% 4) RUN EBM
is_converged = false; %Iteration is not finished
n_attempts   = 0; %Number of attempts at finding iterative solution
full_iter_s  = tic; %Time the entire iteration process
while ~is_converged
    %Experiment w/ different values of dt_c until stable solution is found
    To          = Too; %Initialize model run with temperature profile Too
    n_attempts  = n_attempts+1;
    last_loop_s = tic; %Time the final loop of the iteration process
for ii = 1:max_s
    %Take up to max_s iteration steps
    qo      = relative_to_specific(rh,To,ptot); %Compute specific humidity
    mo      = cp*To+Lv*qo+zg; %Compute low-level MSE on N full lats
    H(2:(end-1)) = H_pf.*diff(mo,1,1); %Compute MSE transports divided by
        %2*pi from differences in adjacent MSE gridboxes, N-1 half lats
    Nh  = -inv_dA.*diff(H,1,1); %Heating rate of each full gridbox due to
        %transport, W/m^2
    Ncs = b-a*To; %Heating rate of each full gridbox due to clear sky OLR,
        %W/m^2
    Nt  = Ns+Ncs+Nh; %Total heating rate of each full gridbox, W/m^2
    if max(abs(Nt))<cnv_t
        %Max. heating rate associated w/ this MSE profile is <cnv_t W/m^2
        last_loop = toc(last_loop_s); %Read out the timers
        full_iter = toc(full_iter_s);
        is_converged = true;
        break; %When this happens, we consider the iteration converged
            %Break exits the for loop, while resetting is_converged causes
            %exit from while loop
    end
    To           = To+dt_c*Nt; %Update temperatures
    if any(To>400) || any(To<0)
        %Check to see if temperatures have become unreasonable
        dt_c  = dt_c/2; %Shorten timestep
        max_s = 2*max_s;%Boost maximum allowed number of iterations
        break; %Leaves the for loop, reinitializes temperature profile
    end
end
    if ii==max_s
        %Finished going through for loop, but still did not converge
        last_loop    = toc(last_loop_s); %Read out the timers
        full_iter    = toc(full_iter_s);
        is_converged = true; %Not really, but we want to terminate run
        ii = -ii; %Flip the sign of ii (and by extension nl) to indicate
            %the lack of convergence
    end
end

% 5) CREATE VARIABLES FOR RETURN
H    = 2*pi*H; %Convert to actual transport in units of W on N+1 half lats
H    = H(2:(end-1)); %Reduce to N-1 half lats
lath = lath(2:(end-1));
nl   = ii; %Number of passages through the loop required to converge
