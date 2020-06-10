clc; close all; clear variables;

addpath(genpath('/project2/tas1/miyawaki/matlab'));

%% set parameters
par.lat_interp = 'std'; % which latitudinal grid to interpolate to: don (donohoe, coarse), era (native ERA-Interim, fine), or std (custom, very fine)
par.lat_std = transpose(-90:0.25:90); % define standard latitude grid for 'std' interpolation
par.ep = 0.1; % threshold for RCE definition. RCE is defined as where abs(R1) < ep
par.ga = 1-par.ep; % threshold for RAE definition. RAE is defined as where R1 > ga

%% call functions
proc_rcae(par)
proc_temp(par)
filt_temp(par)

%% define functions
function proc_rcae(par)
    % read data from Aaron Donohoe's mass-corrected atmospheric energy transport calculations
    % from ERA-Interim climatology from year 2000 through 2012.
    don = load('/project2/tas1/miyawaki/projects/002/data/read/radiation_dynamics_climatology.mat');

    % read ERA-Interim grid data
    load('/project2/tas1/miyawaki/projects/002/data/read/era_grid.mat');
    % read radiation climatology from ERA-Interim, 2000 - 2012 (rad)
    load('/project2/tas1/miyawaki/projects/002/data/read/radiation_climatology.mat');
    % surface turbulent fluxes (stf)
    load('/project2/tas1/miyawaki/projects/002/data/read/turbfluxes_climatology.mat');

    if strcmp(par.lat_interp, 'don')
        lat = don.lat;
        % interpolate raw ERA-Interim data onto Donohoe lat x lon grid
        for fn = rad_vars
            % interpolate to Donohoe lon
            don.(fn{1}) = interp1(lon_era, rad.(fn{1}), don.lon);
            % interpolate to Donohoe lat
            don.(fn{1}) = permute(don.(fn{1}), [2 1 3]);
            don.(fn{1}) = interp1(lat_era, don.(fn{1}), don.lat);
            % order dimensions to be consistent with Donohoe data (mon, lat, lon)
            don.(fn{1}) = permute(don.(fn{1}), [3 1 2]);
        end
        for fn = stf_vars
            % interpolate to Donohoe lon
            don.(fn{1}) = interp1(lon_era, stf.(fn{1}), don.lon);
            % interpolate to Donohoe lat
            don.(fn{1}) = permute(don.(fn{1}), [2 1 3]);
            don.(fn{1}) = interp1(lat_era, don.(fn{1}), don.lat);
            % order dimensions to be consistent with Donohoe data (mon, lat, lon)
            don.(fn{1}) = permute(don.(fn{1}), [3 1 2]);
        end
    elseif strcmp(par.lat_interp, 'era')
        lat = lat_era;
        % interpolate don TETEN and TEDIV data onto ERA lat x lon grid
        for fn = {'TETEN', 'TEDIV'}
            % order dimensions as (lon, lat, mon)
            don.(fn{1}) = permute(don.(fn{1}), [3 2 1]);
            % interpolate to ERA lon
            don.(fn{1}) = interp1(don.lon, don.(fn{1}), lon_era, 'spline');
            % interpolate to ERA lat
            don.(fn{1}) = permute(don.(fn{1}), [2 1 3]);
            don.(fn{1}) = interp1(don.lat, don.(fn{1}), lat_era, 'spline', nan);
            % order dimensions to be consistent with Donohoe data (mon, lat, lon)
            don.(fn{1}) = permute(don.(fn{1}), [3 1 2]);
        end
        % load ERA-Interim data
        for fn = rad_vars
            % order dimensions to be consistent with Donohoe data (mon, lat, lon)
            don.(fn{1}) = permute(rad.(fn{1}), [3 2 1]);
        end
        for fn = stf_vars
            % order dimensions to be consistent with Donohoe data (mon, lat, lon)
            don.(fn{1}) = permute(stf.(fn{1}), [3 2 1]);
        end
    elseif strcmp(par.lat_interp, 'std') % interpolate all to fine standard grid
        lat = par.lat_std;
        % interpolate raw ERA-Interim data onto std lat x ERA lon grid
        for fn = rad_vars
            % interpolate to std lat
            don.(fn{1}) = permute(rad.(fn{1}), [2 1 3]);
            don.(fn{1}) = interp1(lat_era, don.(fn{1}), par.lat_std, 'spline', nan);
            % order dimensions to be consistent with Donohoe data (mon, lat, lon)
            don.(fn{1}) = permute(don.(fn{1}), [3 1 2]);
        end
        for fn = stf_vars
            % interpolate to std lat
            don.(fn{1}) = permute(stf.(fn{1}), [2 1 3]);
            don.(fn{1}) = interp1(lat_era, don.(fn{1}), par.lat_std, 'spline', nan);
            % order dimensions to be consistent with Donohoe data (mon, lat, lon)
            don.(fn{1}) = permute(don.(fn{1}), [3 1 2]);
        end
        % interpolate don TETEN and TEDIV data onto ERA lat x lon grid
        for fn = {'TETEN', 'TEDIV'}
            % order dimensions as (lon, lat, mon)
            don.(fn{1}) = permute(don.(fn{1}), [3 2 1]);
            % interpolate to ERA lon
            don.(fn{1}) = interp1(don.lon, don.(fn{1}), lon_era, 'spline');
            % interpolate to ERA lat
            don.(fn{1}) = permute(don.(fn{1}), [2 1 3]);
            don.(fn{1}) = interp1(don.lat, don.(fn{1}), par.lat_std, 'spline', nan);
            % order dimensions to be consistent with Donohoe data (mon, lat, lon)
            don.(fn{1}) = permute(don.(fn{1}), [3 1 2]);
        end
    end

    % compute net radiative cooling from radiative fluxes
    don.ra = don.tsr - don.ssr + don.ttr - don.str;

    % compute surface turbulent fluxes directly from ERA-Interim data
    % multiply by negative to define flux from surface to atmosphere as positive
    don.stf = -( don.sshf + don.slhf );

    % compute surface turbulent fluxes as residual of radiative cooling, storage, and atmospheric heat transport
    don.stf_res = don.TETEN + don.TEDIV - don.ra;

    % compute MSE storage as residual of radiative cooling, turbulent fluxes, and heat transport
    don.TETEN_res = don.ra + don.stf - don.TEDIV;

    % compute non-dimensional number R1
    % R1 measures the importance of atmospheric heat transport
    % defined using Donohoe TETEN closure
    don.r1_d = (don.TETEN + don.TEDIV) ./ don.ra;
    % defined using ERA-Interim stf closure
    don.r1_e = (don.TETEN_res + don.TEDIV) ./ don.ra;

    % compute non-dimensional number R2
    % R2 measures the importance of surface turbulent fluxes
    % defined using Donohoe TETEN closure
    don.r2_d = don.stf_res ./ don.ra;
    % defined using ERA-Interim stf closure
    don.r2_e = don.stf ./ don.ra;

    % take zonal averages
    for fn = {'ra', 'stf', 'stf_res', 'sshf', 'slhf', 'TETEN', 'TETEN_res', 'TEDIV', 'r1_d', 'r1_e', 'r2_d', 'r2_e'}
        donz.(fn{1}) = nanmean(don.(fn{1}), 3);
    end

    % identify locations of RCE using threshold epsilon (ep)
    donz.rcae_d = zeros(size(donz.r1_d));
    donz.rcae_e = zeros(size(donz.r1_e));
    donz.rcae_d(abs(donz.r1_d) < par.ep) = 1;
    donz.rcae_e(abs(donz.r1_e) < par.ep) = 1;
    % identify locations of RAE using threshold gamma (ga)
    donz.rcae_d(donz.r1_d > par.ga) = -1;
    donz.rcae_e(donz.r1_e > par.ga) = -1;

    % save data into mat file
    foldername = sprintf('/project2/tas1/miyawaki/projects/002/data/proc/eps_%g/', par.ep);
    filename = sprintf('processed_rad_%s.mat', par.lat_interp);
    printname = [foldername filename];
    if ~exist(foldername, 'dir')
        mkdir(foldername)
    end
    save(printname, 'donz', 'lat');

end
function proc_temp(par)
    % read data from Aaron Donohoe's mass-corrected atmospheric energy transport calculations
    % from ERA-Interim climatology from year 2000 through 2012.
    don = load('/project2/tas1/miyawaki/projects/002/data/read/radiation_dynamics_climatology.mat');
    % read ERA-Interim grid data
    load('/project2/tas1/miyawaki/projects/002/data/read/era_grid.mat');
    % load temperature
    load('/project2/tas1/miyawaki/projects/002/data/read/temp_climatology.mat');

    if strcmp(par.lat_interp, 'don')
        lat = don.lat;
        % interpolate raw ERA-Interim data onto Donohoe lat grid
        for fn = vars_3d
            % interpolate to Vertzohoe lat
            vertz.(fn{1}) = interp1(lat_era, vert.(fn{1}), lat);
            % order dimensions to be consistent with Vertzohoe data (mon, lat, plev)
            vertz.(fn{1}) = permute(vertz.(fn{1}), [3 1 2]);
        end
    elseif strcmp(par.lat_interp, 'era')
        lat = lat_era;
        % load ERA-Interim data
        for fn = vars_3d
            % order dimensions to be consistent with Donohoe data (mon, lat, plev)
            vertz.(fn{1}) = permute(vert.(fn{1}), [3 1 2]);
        end
    elseif strcmp(par.lat_interp, 'std') % interpolate all to fine standard grid
        lat = par.lat_std;
        % interpolate raw ERA-Interim data onto std lat x ERA lon grid
        for fn = vars_3d
            % interpolate to std lat
            vertz.(fn{1}) = interp1(lat_era, vert.(fn{1}), par.lat_std, 'spline', nan);
            % order dimensions to be consistent with Donohoe data (mon, lat, plev)
            vertz.(fn{1}) = permute(vertz.(fn{1}), [3 1 2]);
        end
    end

    % save data into mat file
    foldername = sprintf('/project2/tas1/miyawaki/projects/002/data/proc/eps_%g/', par.ep);
    filename = sprintf('processed_vert_%s.mat', par.lat_interp);
    printname = [foldername filename];
    if ~exist(foldername, 'dir')
        mkdir(foldername)
    end
    save(printname, 'vertz', 'lat');

end
function filt_temp(par)
    % load rcae data
    load(sprintf('/project2/tas1/miyawaki/projects/002/data/proc/eps_%g/processed_rad_%s.mat', par.ep, par.lat_interp));
    % load temp
    load(sprintf('/project2/tas1/miyawaki/projects/002/data/proc/eps_%g/processed_vert_%s', par.ep, par.lat_interp));

    % create empty arrays to store filtering array
    rce_d_filt = nan(size(donz.rcae_d));
    rce_e_filt = nan(size(donz.rcae_e));
    rae_d_filt = nan(size(donz.rcae_d));
    rae_e_filt = nan(size(donz.rcae_e));
    % set RCE=1, elsewhere nan
    rce_d_filt(donz.rcae_d==1)=1;
    rce_e_filt(donz.rcae_e==1)=1;
    % set RAE=1, elsewhere nan
    rae_d_filt(donz.rcae_d==-1)=1;
    rae_e_filt(donz.rcae_e==-1)=1;

    % filter temperature
    t_rce_d_filt = rce_d_filt .* vertz.t;
    t_rce_e_filt = rce_e_filt .* vertz.t;
    t_rae_d_filt = rae_d_filt .* vertz.t;
    t_rae_e_filt = rae_e_filt .* vertz.t;

    % averaged temperature profile
    % take time average first
    rce_d_filt_t = squeeze(nanmean(rce_d_filt, 1));
    rce_e_filt_t = squeeze(nanmean(rce_e_filt, 1));
    rae_d_filt_t = squeeze(nanmean(rae_d_filt, 1));
    rae_e_filt_t = squeeze(nanmean(rae_e_filt, 1));
    rce_d_t_t = squeeze(nanmean(t_rce_d_filt, 1));
    rce_e_t_t = squeeze(nanmean(t_rce_e_filt, 1));
    rae_d_t_t = squeeze(nanmean(t_rae_d_filt, 1));
    rae_e_t_t = squeeze(nanmean(t_rae_e_filt, 1));
    % take cosine-weighted zonal average
    vert_filt.rce_d = nansum(cosd(lat).*rce_d_t_t, 1) / nansum(cosd(lat).*rce_d_filt_t');
    vert_filt.rce_e = nansum(cosd(lat).*rce_e_t_t, 1) / nansum(cosd(lat).*rce_e_filt_t');
    vert_filt.rae_d = nansum(cosd(lat).*rae_d_t_t, 1) / nansum(cosd(lat).*rae_d_filt_t');
    vert_filt.rae_e = nansum(cosd(lat).*rae_e_t_t, 1) / nansum(cosd(lat).*rae_e_filt_t');

    % save filtered data
    foldername = sprintf('/project2/tas1/miyawaki/projects/002/data/proc/eps_%g/', par.ep);
    filename = sprintf('processed_vert_filt_%s.mat', par.lat_interp);
    printname = [foldername filename];
    if ~exist(foldername, 'dir')
        mkdir(foldername)
    end
    save(printname, 'vert_filt');
end
