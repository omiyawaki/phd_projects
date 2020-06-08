clc; close all; clear variables;

addpath(genpath('/project2/tas1/miyawaki/matlab'));

% parameters
lat_interp = 'era'; % which latitudinal grid to interpolate to: don (donohoe, coarse) or era (native ERA-Interim, fine)
lat_std = -90:0.25:90; % define standard latitude grid for 'std' interpolation
ep = 0.3; % threshold for RCE definition. RCE is defined as where abs(R1) < ep
ga = 0.7; % threshold for RAE definition. RAE is defined as where R1 > ga

% read data from Aaron Donohoe's mass-corrected atmospheric energy transport calculations
% from ERA-Interim climatology from year 2000 through 2012.
don = load('/project2/tas1/miyawaki/projects/002/data/radiation_dynamics_climatology.mat');

% read radiation climatology from ERA-Interim, 2000 - 2012.
load('/project2/tas1/miyawaki/projects/002/data/radiation_climatology.mat');

if strcmp(lat_interp, 'don')
    lat = don.lat;
    % interpolate raw ERA-Interim data onto Donohoe lat x lon grid
    for fn = {'ssr', 'str', 'tsr', 'ttr', 'sshf', 'slhf'}
        % interpolate to Donohoe lon
        don.(fn{1}) = interp1(interim.lon, interim.(fn{1}), don.lon);
        % interpolate to Donohoe lat
        don.(fn{1}) = permute(don.(fn{1}), [2 1 3]);
        don.(fn{1}) = interp1(interim.lat, don.(fn{1}), don.lat);
        % order dimensions to be consistent with Donohoe data (mon, lat, lon)
        don.(fn{1}) = permute(don.(fn{1}), [3 1 2]);
    end
elseif strcmp(lat_interp, 'era')
    lat = interim.lat;
    % interpolate don TETEN and TEDIV data onto ERA lat x lon grid
    for fn = {'TETEN', 'TEDIV'}
        % order dimensions as (lon, lat, mon)
        don.(fn{1}) = permute(don.(fn{1}), [3 2 1]);
        % interpolate to ERA lon
        don.(fn{1}) = interp1(don.lon, don.(fn{1}), interim.lon, 'spline');
        % interpolate to ERA lat
        don.(fn{1}) = permute(don.(fn{1}), [2 1 3]);
        don.(fn{1}) = interp1(don.lat, don.(fn{1}), interim.lat, 'spline');
        % order dimensions to be consistent with Donohoe data (mon, lat, lon)
        don.(fn{1}) = permute(don.(fn{1}), [3 1 2]);
    end
    % load ERA-Interim data
    for fn = {'ssr', 'str', 'tsr', 'ttr', 'sshf', 'slhf'}
        % order dimensions to be consistent with Donohoe data (mon, lat, lon)
        don.(fn{1}) = permute(interim.(fn{1}), [3 2 1]);
    end
elseif strcmp(lat_interp, 'std') % interpolate all to fine standard grid
    lat = lat_std;
    % interpolate raw ERA-Interim data onto std lat x ERA lon grid
    for fn = {'ssr', 'str', 'tsr', 'ttr', 'sshf', 'slhf'}
        % interpolate to std lat
        don.(fn{1}) = permute(interim.(fn{1}), [2 1 3]);
        don.(fn{1}) = interp1(interim.lat, don.(fn{1}), lat_std, 'splin');
        % order dimensions to be consistent with Donohoe data (mon, lat, lon)
        don.(fn{1}) = permute(don.(fn{1}), [3 1 2]);
    end
    % interpolate don TETEN and TEDIV data onto ERA lat x lon grid
    for fn = {'TETEN', 'TEDIV'}
        % order dimensions as (lon, lat, mon)
        don.(fn{1}) = permute(don.(fn{1}), [3 2 1]);
        % interpolate to ERA lon
        don.(fn{1}) = interp1(don.lon, don.(fn{1}), interim.lon, 'spline');
        % interpolate to ERA lat
        don.(fn{1}) = permute(don.(fn{1}), [2 1 3]);
        don.(fn{1}) = interp1(don.lat, don.(fn{1}), lat_std, 'spline');
        % order dimensions to be consistent with Donohoe data (mon, lat, lon)
        don.(fn{1}) = permute(don.(fn{1}), [3 1 2]);
    end
end

% compute net radiative cooling from radiative fluxes
don.ra = don.tsr - don.ssr + don.ttr - don.str;

% compute surface turbulent fluxes directly from ERA-Interim data
% multiply by negative to define flux from surface to atmosphere as positive
don.tf = -( don.sshf + don.slhf );

% compute surface turbulent fluxes as residual of radiative cooling, storage, and atmospheric heat transport
don.tf_res = don.TETEN + don.TEDIV - don.ra;

% compute MSE storage as residual of radiative cooling, turbulent fluxes, and heat transport
don.TETEN_res = don.ra + don.tf - don.TEDIV;

% compute non-dimensional number R1
% R1 measures the importance of atmospheric heat transport
% defined using Donohoe TETEN closure
don.r1_d = (don.TETEN + don.TEDIV) ./ don.ra;
% defined using ERA-Interim tf closure
don.r1_e = (don.TETEN_res + don.TEDIV) ./ don.ra;

% compute non-dimensional number R2
% R2 measures the importance of surface turbulent fluxes
% defined using Donohoe TETEN closure
don.r2_d = don.tf_res ./ don.ra;
% defined using ERA-Interim tf closure
don.r2_e = don.tf ./ don.ra;

% identify locations of RCE using threshold epsilon (ep)
% take the zonal average of r1
don.r1_d_zon = nanmean(don.r1_d, 3);
don.r1_e_zon = nanmean(don.r1_e, 3);
don.rcae_d = zeros(size(don.r1_d_zon));
don.rcae_e = zeros(size(don.r1_e_zon));
don.rcae_d(abs(don.r1_d_zon) < ep) = 1;
don.rcae_e(abs(don.r1_e_zon) < ep) = 1;
% identify locations of RAE using threshold gamma (ga)
don.rcae_d(don.r1_d_zon > ga) = -1;
don.rcae_e(don.r1_e_zon > ga) = -1;

% save data into mat file
save(['/project2/tas1/miyawaki/projects/002/data/processed_data_' lat_interp '.mat'], 'don', 'lat');
