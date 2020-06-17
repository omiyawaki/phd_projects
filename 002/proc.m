clc; close all; clear variables;

addpath(genpath('/project2/tas1/miyawaki/matlab'));

%% set parameters
par.lat_interp = 'std'; % which latitudinal grid to interpolate to: don (donohoe, coarse), era (native ERA-Interim, fine), or std (custom, very fine)
par.lat_std = transpose(-90:0.25:90); % define standard latitude grid for 'std' interpolation
par.ep_swp = 0.3; % threshold for RCE definition. RCE is defined as where abs(R1) < ep
par.ep_cp = 0.5; % threshold for RCE definition. RCE is defined as where abs(R1) < ep
par.cpd = 1005.7; par.Rd = 287; par.g = 9.81; par.L = 2.501e6; par.a = 6357e3;

%% call functions
% proc_era_vh(par) % calculate integrated energy transport (in units of power) using Donohoe MSE flux divergence
% proc_era_net_fluxes(par) % calculate global TOA and surface energy flux imbalance
% proc_era5_net_fluxes(par) % calculate global TOA and surface energy flux imbalance
for i=1:length(par.ep_swp); par.ep = par.ep_swp(i); par.ga = 1-par.ep;
    % proc_era_rcae(par) % calculate energy fluxes in the vertically-integrated MSE budget using ERA-Interim data
    % proc_era_temp(par) % extract raw temperature data from ERA-Interim
    % filt_era_temp(par) % calculate temperature profiles over RCE and RAE using ERA-Interim data
    % proc_era5_fluxes(par) % calculate energy fluxes in the vertically-integrated MSE budget using ERA-Interim data
    % proc_era5_rcae(par) % calculate energy fluxes in the vertically-integrated MSE budget using ERA-Interim data
    % proc_era5_temp(par) % extract raw temperature data from ERA5
    % proc_era5_srfc(par) % extract raw surface data from ERA5
    filt_era5_temp(par) % calculate temperature profiles over RCE and RAE using ERA5 data
    % proc_mpi_rcae(par) % calculate energy fluxes in the vertically-integrated MSE budget using MPI-ESM-LR data
end

%% define functions
function proc_era_rcae(par)
    % ead data from Aaron Donohoe's mass-corrected atmospheric energy transport calculations
    % from ERA-Interim climatology from year 2000 through 2012.
    don = load('/project2/tas1/miyawaki/projects/002/data/read/radiation_dynamics_climatology.mat');

    % read ERA-Interim grid data
    load('/project2/tas1/miyawaki/projects/002/data/read/era-interim/grid.mat');
    % read radiation climatology from ERA-Interim, 2000 - 2012 (rad)
    load('/project2/tas1/miyawaki/projects/002/data/read/era-interim/radiation_climatology.mat');
    % surface turbulent fluxes (stf)
    load('/project2/tas1/miyawaki/projects/002/data/read/era-interim/turbfluxes_climatology.mat');

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
    don.r1_teten = (don.TETEN + don.TEDIV) ./ don.ra;
    % defined using ERA-Interim stf closure
    don.r1_stf = (don.TETEN_res + don.TEDIV) ./ don.ra;

    % compute non-dimensional number R2
    % R2 measures the importance of surface turbulent fluxes
    % defined using Donohoe TETEN closure
    don.r2_teten = don.stf_res ./ don.ra;
    % defined using ERA-Interim stf closure
    don.r2_stf = don.stf ./ don.ra;

    % take zonal averages
    for fn = {'ra', 'stf', 'stf_res', 'sshf', 'slhf', 'TETEN', 'TETEN_res', 'TEDIV', 'r1_teten', 'r1_stf', 'r2_teten', 'r2_stf'}
        fluxez.(fn{1}) = nanmean(don.(fn{1}), 3);
    end

    % identify locations of RCE using threshold epsilon (ep)
    rcae.teten = zeros(size(fluxez.r1_teten));
    rcae.stf = zeros(size(fluxez.r1_stf));
    rcae.teten(abs(fluxez.r1_teten) < par.ep) = 1;
    rcae.stf(abs(fluxez.r1_stf) < par.ep) = 1;
    % identify locations of RAE using threshold gamma (ga)
    rcae.teten(fluxez.r1_teten > par.ga) = -1;
    rcae.stf(fluxez.r1_stf > par.ga) = -1;

    % save energy flux data into mat file
    foldername = sprintf('/project2/tas1/miyawaki/projects/002/data/proc/era/%s/', par.lat_interp);
    printname = [foldername 'fluxes.mat'];
    if ~exist(foldername, 'dir')
        mkdir(foldername)
    end
    save(printname, 'fluxez', 'lat');
    % save rcae data
    foldername = sprintf('/project2/tas1/miyawaki/projects/002/data/proc/era/%s/eps_%g/', par.lat_interp, par.ep);
    printname = [foldername 'rcae.mat'];
    if ~exist(foldername, 'dir')
        mkdir(foldername)
    end
    save(printname, 'rcae', 'lat');

end
function proc_era_net_fluxes(par)
% calculates the global TOA energy imbalance using ERA-Interim data
    % read ERA-Interim grid data
    load('/project2/tas1/miyawaki/projects/002/data/read/era-interim/grid.mat');
    % read radiation climatology from ERA-Interim, 2000 - 2012 (rad)
    load('/project2/tas1/miyawaki/projects/002/data/read/era-interim/radiation_climatology.mat');
    % surface turbulent fluxes (stf)
    load('/project2/tas1/miyawaki/projects/002/data/read/era-interim/turbfluxes_climatology.mat');

    net_toa_raw = rad.tsr + rad.ttr; % compute net radiative fluxes at TOA, positive down
    net_toa_tz = squeeze(nanmean(nanmean( net_toa_raw, 1 ), 3))'; % take zonal and time averages and transpose to have same dimensions as latitude grid
    net_toa = nansum(cosd(lat_era).*net_toa_tz) / nansum(cosd(lat_era));
    disp( sprintf('The net radiative imbalance at TOA is %g Wm^-2.', net_toa) );

    net_sfc_raw = rad.ssr + rad.str + stf.sshf + stf.slhf; % compute net radiative fluxes at surface, positive down
    net_sfc_tz = squeeze(nanmean(nanmean( net_sfc_raw, 1 ), 3))'; % take zonal and time averages and transpose to have same dimensions as latitude grid
    net_sfc = nansum(cosd(lat_era).*net_sfc_tz) / nansum(cosd(lat_era));
    disp( sprintf('The net radiative imbalance at the surface is %g Wm^-2.', net_sfc) );
end
function proc_era_vh(par)
% calculates the power transported by the atmosphere at every latitude
    % read data from Aaron Donohoe's mass-corrected atmospheric energy transport calculations
    % from ERA-Interim climatology from year 2000 through 2012.
    don = load('/project2/tas1/miyawaki/projects/002/data/read/radiation_dynamics_climatology.mat');

    tediv_t = squeeze(nanmean(don.TEDIV, 1)); % take time average
    tediv_tz = trapz(deg2rad(don.lon), tediv_t, 2); % zonally integrate
    vh = cumtrapz(deg2rad(don.lat), par.a^2*cosd(don.lat).*tediv_tz); % cumulatively integrate in latitude

    if strcmp(par.lat_interp, 'don')
        lat = don.lat;
    elseif strcmp(par.lat_interp, 'era')
        load('/project2/tas1/miyawaki/projects/002/data/read/era-interim/grid.mat'); % read ERA-Interim grid data
        lat = lat_era;
        vh = interp1(don.lat, vh, lat_era, 'spline', nan); % interpolate to ERA lat
    elseif strcmp(par.lat_interp, 'std') % interpolate all to fine standard grid
        lat = par.lat_std;
        vh = interp1(don.lat, vh, par.lat_std, 'spline', nan); % interpolate to standard lat
    end

    % save data into mat file
    foldername = sprintf('/project2/tas1/miyawaki/projects/002/data/proc/era/%s/', par.lat_interp);
    printname = [foldername 'vh.mat'];
    if ~exist(foldername, 'dir')
        mkdir(foldername)
    end
    save(printname, 'vh', 'lat');

end
function proc_era_temp(par)
    % read data from Aaron Donohoe's mass-corrected atmospheric energy transport calculations
    % from ERA-Interim climatology from year 2000 through 2012.
    don = load('/project2/tas1/miyawaki/projects/002/data/read/era-interim/radiation_dynamics_climatology.mat');
    % read ERA-Interim grid data
    load('/project2/tas1/miyawaki/projects/002/data/read/era-interim/grid.mat');
    % load temperature
    load('/project2/tas1/miyawaki/projects/002/data/read/era-interim/temp_climatology.mat');

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
    foldername = sprintf('/project2/tas1/miyawaki/projects/002/data/proc/era/%s/', par.lat_interp);
    printname = [foldername 'vert.mat'];
    if ~exist(foldername, 'dir')
        mkdir(foldername)
    end
    save(printname, 'vertz', 'lat');

end
function filt_era_temp(par)
    % load fluxes
    load(sprintf('/project2/tas1/miyawaki/projects/002/data/proc/era/%s/fluxes.mat', par.lat_interp));
    % load rcae data
    load(sprintf('/project2/tas1/miyawaki/projects/002/data/proc/era/%s/eps_%g/rcae.mat', par.lat_interp, par.ep));
    % load temp
    load(sprintf('/project2/tas1/miyawaki/projects/002/data/proc/era/%s/vert', par.lat_interp));

    % create empty arrays to store filtering array
    rce_teten_filt = nan(size(rcae.teten));
    rce_stf_filt = nan(size(rcae.stf));
    rae_teten_filt = nan(size(rcae.teten));
    rae_stf_filt = nan(size(rcae.stf));
    % set RCE=1, elsewhere nan
    rce_teten_filt(rcae.teten==1)=1;
    rce_stf_filt(rcae.stf==1)=1;
    % set RAE=1, elsewhere nan
    rae_teten_filt(rcae.teten==-1)=1;
    rae_stf_filt(rcae.stf==-1)=1;

    % filter temperature
    t_rce_teten_filt = rce_teten_filt .* vertz.t;
    t_rce_stf_filt = rce_stf_filt .* vertz.t;
    t_rae_teten_filt = rae_teten_filt .* vertz.t;
    t_rae_stf_filt = rae_stf_filt .* vertz.t;

    % averaged temperature profile
    % take time average first
    rce_teten_filt_t = squeeze(nanmean(rce_teten_filt, 1));
    rce_stf_filt_t = squeeze(nanmean(rce_stf_filt, 1));
    rae_teten_filt_t = squeeze(nanmean(rae_teten_filt, 1));
    rae_stf_filt_t = squeeze(nanmean(rae_stf_filt, 1));
    rce_teten_t_t = squeeze(nanmean(t_rce_teten_filt, 1));
    rce_stf_t_t = squeeze(nanmean(t_rce_stf_filt, 1));
    rae_teten_t_t = squeeze(nanmean(t_rae_teten_filt, 1));
    rae_stf_t_t = squeeze(nanmean(t_rae_stf_filt, 1));
    % take cosine-weighted zonal average
    vert_filt.rce_teten = nansum(cosd(lat).*rce_teten_t_t, 1) / nansum(cosd(lat).*rce_teten_filt_t');
    vert_filt.rce_stf = nansum(cosd(lat).*rce_stf_t_t, 1) / nansum(cosd(lat).*rce_stf_filt_t');
    vert_filt.rae_teten = nansum(cosd(lat).*rae_teten_t_t, 1) / nansum(cosd(lat).*rae_teten_filt_t');
    vert_filt.rae_stf = nansum(cosd(lat).*rae_stf_t_t, 1) / nansum(cosd(lat).*rae_stf_filt_t');

    % save filtered data
    foldername = sprintf('/project2/tas1/miyawaki/projects/002/data/proc/era/%s/eps_%g/', par.lat_interp, par.ep);
    printname = [foldername 'vert_filt'];
    if ~exist(foldername, 'dir')
        mkdir(foldername)
    end
    save(printname, 'vert_filt');
end
function proc_era5_fluxes(par)
    % read ERA5 grid data
    load('/project2/tas1/miyawaki/projects/002/data/read/era5/grid.mat');
    % read radiation climatology from ERA5
    load('/project2/tas1/miyawaki/projects/002/data/read/era5/radiation_climatology.mat');
    % read hydrology climatology from ERA5
    load('/project2/tas1/miyawaki/projects/002/data/read/era5/pe_climatology.mat');
    % surface turbulent fluxes (stf)
    load('/project2/tas1/miyawaki/projects/002/data/read/era5/turbfluxes_climatology.mat');

    if strcmp(par.lat_interp, 'era5')
        lat = lat_era;
        % load ERA5 data
        for fn = rad_vars
            % order dimensions to be consistent with Donohoe data (mon, lat, lon)
            era5.(fn{1}) = permute(rad.(fn{1}), [3 2 1]);
        end
        for fn = pe_vars
            % order dimensions to be consistent with Donohoe data (mon, lat, lon)
            era5.(fn{1}) = permute(pe.(fn{1}), [3 2 1]);
        end
        for fn = stf_vars
            % order dimensions to be consistent with Donohoe data (mon, lat, lon)
            era5.(fn{1}) = permute(stf.(fn{1}), [3 2 1]);
        end
    elseif strcmp(par.lat_interp, 'std') % interpolate all to fine standard grid
        lat = par.lat_std;
        % interpolate raw ERA5 data onto std lat x ERA5 lon grid
        for fn = rad_vars
            % interpolate to std lat
            era5.(fn{1}) = permute(rad.(fn{1}), [2 1 3]);
            era5.(fn{1}) = interp1(lat_era, era5.(fn{1}), par.lat_std, 'spline', nan);
            % order dimensions to be consistent with Era5ohoe data (mon, lat, lon)
            era5.(fn{1}) = permute(era5.(fn{1}), [3 1 2]);
        end
        for fn = pe_vars
            % interpolate to std lat
            era5.(fn{1}) = permute(pe.(fn{1}), [2 1 3]);
            era5.(fn{1}) = interp1(lat_era, era5.(fn{1}), par.lat_std, 'spline', nan);
            % order dimensions to be consistent with Era5ohoe data (mon, lat, lon)
            era5.(fn{1}) = permute(era5.(fn{1}), [3 1 2]);
        end
        for fn = stf_vars
            % interpolate to std lat
            era5.(fn{1}) = permute(stf.(fn{1}), [2 1 3]);
            era5.(fn{1}) = interp1(lat_era, era5.(fn{1}), par.lat_std, 'spline', nan);
            % order dimensions to be consistent with Era5ohoe data (mon, lat, lon)
            era5.(fn{1}) = permute(era5.(fn{1}), [3 1 2]);
        end
    end

    % compute net radiative cooling from radiative fluxes
    era5.ra = era5.tsr - era5.ssr + era5.ttr - era5.str;

    % compute surface turbulent fluxes directly from ERA5 data
    % multiply by negative to define flux from surface to atmosphere as positive
    era5.stf = -( era5.sshf + era5.slhf );

    % compute MSE storage and flux divergence as residual of radiative cooling and surface turbulent fluxes
    era5.res = era5.ra + era5.stf;
    % compute northward MSE transport using the residual data
    tediv_t = squeeze(nanmean(era5.res, 1)); % take time average
    tediv_tz = trapz(deg2rad(lon_era), tediv_t, 2); % zonally integrate
    vh = cumtrapz(deg2rad(lat), par.a^2*cosd(lat).*tediv_tz); % cumulatively integrate in latitude

    % compute non-dimensional number R1
    % R1 measures the importance of atmospheric heat transport
    era5.r1 = (era5.res) ./ era5.ra;

    % compute non-dimensional number R2
    % R2 measures the importance of surface turbulent fluxes
    era5.r2 = era5.stf ./ era5.ra;

    % take zonal averages
    for fn = {'ra', 'stf', 'sshf', 'slhf', 'res', 'r1', 'r2', 'cp', 'lsp', 'e'}
        fluxez.(fn{1}) = nanmean(era5.(fn{1}), 3);
    end

    % save energy flux data into mat file
    foldername = sprintf('/project2/tas1/miyawaki/projects/002/data/proc/era5/%s/', par.lat_interp);
    printname = [foldername 'fluxes.mat'];
    if ~exist(foldername, 'dir')
        mkdir(foldername)
    end
    save(printname, 'fluxez', 'vh', 'lat');

end
function proc_era5_rcae(par)
    % read ERA5 zonally averaged fluxes
    filename = sprintf('/project2/tas1/miyawaki/projects/002/data/proc/era5/%s/fluxes.mat', par.lat_interp);
    if ~exist(filename); error('Data does not exist. Please run proc_era5_fluxes.m first.'); else
        load(filename);
    end

    % identify locations of RCE using threshold epsilon (ep)
    rcae.def = zeros(size(fluxez.r1));
    rcae.def(abs(fluxez.r1) < par.ep) = 1;
    % identify locations of RAE using threshold gamma (ga)
    rcae.def(fluxez.r1 > par.ga) = -1;
    % add additional criteria for RCE that P-E>0
    rcae.pe = zeros(size(fluxez.r1));
    rcae.pe(abs(fluxez.r1) < par.ep & (fluxez.cp+fluxez.lsp+fluxez.e>0)) = 1; % note that evap is defined negative into atmosphere
    rcae.pe(fluxez.r1 > par.ga) = -1;
    % add additional criteria for RCE that (P_ls - E)<<1
    rcae.cp = zeros(size(fluxez.r1));
    rcae.cp(abs(fluxez.r1) < par.ep & abs(fluxez.lsp./fluxez.cp<par.ep_cp)) = 1; % note that evap is defined negative into atmosphere
    rcae.cp(fluxez.r1 > par.ga) = -1;

    % save rcae data
    foldername = sprintf('/project2/tas1/miyawaki/projects/002/data/proc/era5/%s/eps_%g/', par.lat_interp, par.ep);
    printname = [foldername 'rcae.mat'];
    if ~exist(foldername, 'dir')
        mkdir(foldername)
    end
    save(printname, 'rcae', 'lat');

end
function proc_era5_net_fluxes(par)
% calculates the global TOA energy imbalance using ERA-Interim data
    % read ERA5 grid data
    load('/project2/tas1/miyawaki/projects/002/data/read/era5/grid.mat');
    % read radiation climatology from ERA5, 2000 - 2012 (rad)
    load('/project2/tas1/miyawaki/projects/002/data/read/era5/radiation_climatology.mat');
    % surface turbulent fluxes (stf)
    load('/project2/tas1/miyawaki/projects/002/data/read/era5/turbfluxes_climatology.mat');

    net_toa_raw = rad.tsr + rad.ttr; % compute net radiative fluxes at TOA, positive down
    net_toa_tz = squeeze(nanmean(nanmean( net_toa_raw, 1 ), 3))'; % take zonal and time avera5ges and transpose to have same dimensions as latitude grid
    net_toa = nansum(cosd(lat_era).*net_toa_tz) / nansum(cosd(lat_era));
    disp( sprintf('The net radiative imbalance at TOA is %g Wm^-2.', net_toa) );

    net_sfc_raw = rad.ssr + rad.str + stf.sshf + stf.slhf; % compute net radiative fluxes at surface, positive down
    net_sfc_tz = squeeze(nanmean(nanmean( net_sfc_raw, 1 ), 3))'; % take zonal and time avera5ges and transpose to have same dimensions as latitude grid
    net_sfc = nansum(cosd(lat_era).*net_sfc_tz) / nansum(cosd(lat_era));
    disp( sprintf('The net radiative imbalance at the surface is %g Wm^-2.', net_sfc) );
end
function proc_era5_temp(par)
    load('/project2/tas1/miyawaki/projects/002/data/read/era5/grid.mat'); % read ERA5 grid data
    load('/project2/tas1/miyawaki/projects/002/data/read/era5/temp_climatology.mat'); % load temperature

    if strcmp(par.lat_interp, 'era5')
        lat = lat_era;
        % load ERA5 data
        for fn = vars_3d
            % order dimensions to be consistent with Donohoe data (mon, lat, plev)
            vertz.(fn{1}) = permute(vert.(fn{1}), [3 1 2]);
        end
    elseif strcmp(par.lat_interp, 'std') % interpolate all to fine standard grid
        lat = par.lat_std;
        % interpolate raw ERA5 data onto std lat x ERA5 lon grid
        for fn = vars_3d
            % interpolate to std lat
            vertz.(fn{1}) = interp1(lat_era, vert.(fn{1}), par.lat_std, 'spline', nan);
            % order dimensions to be consistent with Donohoe data (mon, lat, plev)
            vertz.(fn{1}) = permute(vertz.(fn{1}), [3 1 2]);
        end
    end

    % save data into mat file
    foldername = sprintf('/project2/tas1/miyawaki/projects/002/data/proc/era5/%s/', par.lat_interp);
    printname = [foldername 'vert.mat'];
    if ~exist(foldername, 'dir')
        mkdir(foldername)
    end
    save(printname, 'vertz', 'lat');

end
function proc_era5_srfc(par)
    load('/project2/tas1/miyawaki/projects/002/data/read/era5/grid.mat'); % read ERA5 grid data
    load('/project2/tas1/miyawaki/projects/002/data/read/era5/srfc_climatology.mat'); % load surface data

    if strcmp(par.lat_interp, 'era5')
        lat = lat_era;
        % load ERA5 data
        for fn = vars_srfc
            srfcz.(fn{1}) = permute(srfc.(fn{1}), [2 1]); % order dimensions to be consistent with Donohoe data (mon, lat, plev)
        end
    elseif strcmp(par.lat_interp, 'std') % interpolate all to fine standard grid
        lat = par.lat_std;
        % interpolate raw ERA5 data onto std lat x ERA5 lon grid
        for fn = vars_srfc
            srfcz.(fn{1}) = interp1(lat_era, srfc.(fn{1}), par.lat_std, 'spline', nan); % interpolate to std lat
            srfcz.(fn{1}) = permute(srfcz.(fn{1}), [2 1]); % order dimensions to be consistent with Donohoe data (mon, lat, plev)
        end
    end

    % save data into mat file
    foldername = sprintf('/project2/tas1/miyawaki/projects/002/data/proc/era5/%s/', par.lat_interp);
    printname = [foldername 'srfc.mat'];
    if ~exist(foldername, 'dir')
        mkdir(foldername)
    end
    save(printname, 'srfcz', 'lat');

end
function filt_era5_temp(par)
    load('/project2/tas1/miyawaki/projects/002/data/read/era5/grid.mat'); % read ERA5 grid data
    load(sprintf('/project2/tas1/miyawaki/projects/002/data/proc/era5/%s/fluxes.mat', par.lat_interp)); % load fluxes
    load(sprintf('/project2/tas1/miyawaki/projects/002/data/proc/era5/%s/eps_%g/rcae.mat', par.lat_interp, par.ep)); % load rcae data
    load(sprintf('/project2/tas1/miyawaki/projects/002/data/proc/era5/%s/vert.mat', par.lat_interp)); % load temp
    load(sprintf('/project2/tas1/miyawaki/projects/002/data/proc/era5/%s/srfc.mat', par.lat_interp)); % load surface data

    % create surface mask
    ps_3d = repmat(srfcz.sp, [1 1 size(vertz.t, 3)]);
    pa_3d = double(permute(repmat(plev_era, [1 size(srfcz.sp)]), [2 3 1]))*10^2; % convert from hPa to Pa
    surface_mask = nan(size(vertz.t));
    surface_mask(pa_3d < ps_3d) = 1;

    % filter all 3D data with surface mask
    for fn = fieldnames(vertz)'
        vertz.(fn{1}) = vertz.(fn{1}).*surface_mask;
    end
    pa_3d = pa_3d.*surface_mask;

    for fn = fieldnames(rcae)'
        rce_filt.(fn{1}) = nan(size(rcae.(fn{1}))); % create empty arrays to store filtering array
        rae_filt.(fn{1}) = nan(size(rcae.(fn{1})));
        rce_filt.(fn{1})(rcae.(fn{1})==1)=1; % set RCE=1, elsewhere nan
        rae_filt.(fn{1})(rcae.(fn{1})==-1)=1; % set RAE=1, elsewhere nan
        t_rce.(fn{1}) = rce_filt.(fn{1}) .* vertz.t; % filter temperatures in RCE only
        t_rae.(fn{1}) = rae_filt.(fn{1}) .* vertz.t; % filter temperatures in RAE only
        ps_rce.(fn{1}) = rce_filt.(fn{1}) .* srfcz.sp; % filter surface pressures in RCE only
        ps_rae.(fn{1}) = rae_filt.(fn{1}) .* srfcz.sp; % filter surface pressures in RAE only
        rce_filt_t.(fn{1}) = squeeze(nanmean(rce_filt.(fn{1}), 1)); % take time mean of RCE filter
        rae_filt_t.(fn{1}) = squeeze(nanmean(rae_filt.(fn{1}), 1)); % take time mean of RAE filter
        t_rce_t.(fn{1}) = squeeze(nanmean(t_rce.(fn{1}), 1)); % take time mean of temperature
        t_rae_t.(fn{1}) = squeeze(nanmean(t_rae.(fn{1}), 1)); % take time mean of temperature
        vert_filt.rce.(fn{1}) = nansum(cosd(lat).*t_rce_t.(fn{1})) / nansum(cosd(lat).*rce_filt_t.(fn{1})'); % area-weighted meridional average
        vert_filt.rce.tp.(fn{1}) = nansum(cosd(lat(abs(lat)<30)).*t_rce_t.(fn{1})(abs(lat)<30,:)) / nansum(cosd(lat(abs(lat)<30)).*rce_filt_t.(fn{1})(abs(lat)<30)'); % area-weighted meridional average
        vert_filt.rce.nh.(fn{1}) = nansum(cosd(lat(lat>30)).*t_rce_t.(fn{1})(lat>30,:)) / nansum(cosd(lat(lat>30)).*rce_filt_t.(fn{1})(lat>30)'); % area-weighted meridional average
        vert_filt.rce.sh.(fn{1}) = nansum(cosd(lat(lat<-30)).*t_rce_t.(fn{1})(lat<-30,:)) / nansum(cosd(lat(lat<-30)).*rce_filt_t.(fn{1})(lat<-30)'); % area-weighted meridional average
        vert_filt.rae.(fn{1}) = nansum(cosd(lat).*t_rae_t.(fn{1})) / nansum(cosd(lat).*rae_filt_t.(fn{1})'); % area-weighted meridional average
        vert_filt.rae.nh.(fn{1}) = nansum(cosd(lat(lat>0)).*t_rae_t.(fn{1})(lat>0,:)) / nansum(cosd(lat(lat>0)).*rae_filt_t.(fn{1})(lat>0)'); % area-weighted meridional average
        vert_filt.rae.sh.(fn{1}) = nansum(cosd(lat(lat<0)).*t_rae_t.(fn{1})(lat<0,:)) / nansum(cosd(lat(lat<0)).*rae_filt_t.(fn{1})(lat<0)'); % area-weighted meridional average
        vert_filt.rce.(fn{1})(plev_era*100 > min(min(ps_rce.(fn{1})))) = nan;
        vert_filt.rce.tp.(fn{1})(plev_era*100 > min(min(ps_rce.(fn{1})(:,abs(lat)<30)))) = nan;
        vert_filt.rce.nh.(fn{1})(plev_era*100 > min(min(ps_rce.(fn{1})(:,lat>30)))) = nan;
        vert_filt.rce.sh.(fn{1})(plev_era*100 > min(min(ps_rce.(fn{1})(:,lat<-30)))) = nan;
        vert_filt.rae.(fn{1})(plev_era*100 > min(min(ps_rae.(fn{1})))) = nan;
        vert_filt.rae.nh.(fn{1})(plev_era*100 > min(min(ps_rae.(fn{1})(:,lat>0)))) = nan;
        vert_filt.rae.sh.(fn{1})(plev_era*100 > min(min(ps_rae.(fn{1})(:,lat<0)))) = nan;
    end

    % save filtered data
    foldername = sprintf('/project2/tas1/miyawaki/projects/002/data/proc/era5/%s/eps_%g/', par.lat_interp, par.ep);
    printname = [foldername 'vert_filt'];
    if ~exist(foldername, 'dir')
        mkdir(foldername)
    end
    save(printname, 'vert_filt');
end
function proc_mpi_rcae(par)
    load('/project2/tas1/miyawaki/projects/002/data/read/mpi/2d_climatology'); % read 2D mpi data
    load('/project2/tas1/miyawaki/projects/002/data/read/mpi/3d_climatology'); % read 3D mpi data

    lat = mpi_3d.lat;

    % create surface mask
    ps_4d = permute(repmat(mpi_2d.ps, [1 1 1 size(mpi_3d.ta, 3)]), [1 2 4 3]);
    pa_4d = permute(repmat(mpi_3d.plev, [1 size(mpi_2d.ps)]), [2 3 1 4]);
    surface_mask = nan(size(mpi_3d.ta));
    surface_mask(pa_4d < ps_4d) = 1;

    % filter all 3D data with surface mask
    for fn = fieldnames(mpi_3d)'
        if ~strcmp(fn{1}, {'lon', 'lat', 'plev'})
            mpi_3d.(fn{1}) = mpi_3d.(fn{1}).*surface_mask;
        end
    end
    pa_4d = pa_4d.*surface_mask;

    % calculate x component of MSE flux divergence
    dlon = mpi_3d.lon(2)-mpi_3d.lon(1);
    dx = par.a*cosd(mpi_3d.lat)*deg2rad(dlon);
    dx_4d = repmat(dx, [1 size(mpi_3d.ta,1) size(mpi_3d.ta,3) size(mpi_3d.ta,4)]);
    dx_4d = permute(dx_4d, [2 1 3 4]);
    lon_half = (mpi_3d.lon(2:end) + mpi_3d.lon(1:end-1))/2;
    lon_half(end+1) = (mpi_3d.lon(end)+360)/2;
    duh_half = (mpi_3d.uh(2:end,:,:,:) - mpi_3d.uh(1:end-1,:,:,:))./dx_4d(1:end-1,:,:,:);
    duh_half(end+1,:,:,:) = (mpi_3d.uh(1,:,:,:) - mpi_3d.uh(end,:,:,:))./dx_4d(1,:,:,:);
    duh = interp1(lon_half, duh_half, mpi_3d.lon);

    % calculate y component of MSE flux divergence
    dlat = mpi_3d.lat(2) - mpi_3d.lat(1);
    dy = par.a*deg2rad(dlat);
    lat_half = (mpi_3d.lat(2:end) + mpi_3d.lat(1:end-1))/2;
    cos_4d = repmat( cosd(mpi_3d.lat), [1 size(mpi_3d.ta,1), size(mpi_3d.ta,3), size(mpi_3d.ta,4)] );
    cos_4d = permute(cos_4d, [2 1 3 4]);
    dvh_half = (cos_4d(:,2:end,:,:).*mpi_3d.vh(:,2:end,:,:) - cos_4d(:,1:end-1,:,:).*mpi_3d.vh(:,1:end-1,:,:))/dy;
    dvh_half = permute(dvh_half, [2 1 3 4]);
    dvh = interp1(lat_half, dvh_half, mpi_3d.lat);
    dvh = permute(dvh, [2 1 3 4]);

    % vertically integrate
    % fluxes.TEDIV = squeeze( trapz(mpi_3d.plev, -1/par.g*(duh + dvh), 3) ); % total horizontal MSE flux divergence
    for i = 1:length(mpi_3d.lon); long = mpi_3d.lon(i);
        for j = 1:length(mpi_3d.lat); lati = mpi_3d.lat(j);
            fluxes.TEDIV(i,j,:) = squeeze( trapz(squeeze(pa_4d(i,j,:,1)), -1/par.g*(duh(i,j,:,:) + dvh(i,j,:,:)), 3) ); % total horizontal MSE flux divergence
        end
    end
    h = squeeze( trapz(mpi_3d.plev, -1/par.g*(mpi_3d.h), 3) ); % MSE

    % MSE tendency
    dhdt_half = nan(size(h));
    dhdt_half(:,:,2:end) = (h(:,:,2:end) - h(:,:,1:end-1)) / (30.5*86400);
    dhdt_half(:,:,1) = (h(:,:,1) - h(:,:,end)) / (30.5*86400);
    dhdt_half = permute(dhdt_half, [3 1 2]);
    fluxes.TETEN = interp1(0.5:11.5, dhdt_half, 1:12);
    fluxes.TETEN = permute(fluxes.TETEN, [2 3 1]);

    load('/project2/tas1/miyawaki/projects/002/data/read/mpi/2d_climatology'); % read 2D mpi data
    fluxes.slhf = -mpi_2d.hfls; fluxes.sshf = -mpi_2d.hfss; % rename and change sign of turbulent fluxes to be consistent with era
    fluxes.ra = mpi_2d.rsdt - mpi_2d.rsut + mpi_2d.rsus - mpi_2d.rsds + mpi_2d.rlus - mpi_2d.rlds - mpi_2d.rlut; % calculate atmospheric radiative cooling
    fluxes.stf = fluxes.slhf + fluxes.sshf;

    % infer energy fluxes as residuals
    fluxes.TETEN_res = fluxes.ra + fluxes.slhf + fluxes.sshf - fluxes.TEDIV;
    fluxes.stf_res = fluxes.TETEN + fluxes.TEDIV - fluxes.ra;

    fluxes.r1 = (fluxes.TETEN + fluxes.TEDIV)./fluxes.ra; % calculate nondimensional number R1 disregarding MSE budget closure
    fluxes.r2 = fluxes.stf./fluxes.ra; % calculate nondimensional number R2 disregarding MSE budget closure
    fluxes.r1_teten = (fluxes.TETEN + fluxes.TEDIV)./fluxes.ra; % calculate nondimensional number R1 using MSE tendency as closure
    fluxes.r2_teten = fluxes.stf_res./fluxes.ra; % calculate nondimensional number R2 using MSE tendency as closure
    fluxes.r1_stf = (fluxes.TETEN_res + fluxes.TEDIV)./fluxes.ra; % calculate nondimensional number R1 using turbulent fluxes as closure
    fluxes.r2_stf = fluxes.stf./fluxes.ra; % calculate nondimensional number R2 using turbulent fluxes as closure

    % take zonal averages
    for fn = fieldnames(fluxes)'
        fluxez.(fn{1}) = squeeze(nanmean(fluxes.(fn{1}), 1));
        fluxez.(fn{1}) = permute(fluxez.(fn{1}), [2 1]); % rearrange to mon x lat
    end

    % identify locations of RCE using threshold epsilon (ep)
    rcae.none = zeros(size(fluxez.r1));
    rcae.teten = zeros(size(fluxez.r1_teten));
    rcae.stf = zeros(size(fluxez.r1_stf));
    rcae.none(abs(fluxez.r1) < par.ep) = 1;
    rcae.teten(abs(fluxez.r1_teten) < par.ep) = 1;
    rcae.stf(abs(fluxez.r1_stf) < par.ep) = 1;
    % identify locations of RAE using threshold gamma (ga)
    rcae.none(abs(fluxez.r2) < par.ep) = -1;
    rcae.teten(abs(fluxez.r2_teten) < par.ep) = -1;
    rcae.stf(abs(fluxez.r2_stf) < par.ep) = -1;

    % save energy flux data into mat file
    foldername = sprintf('/project2/tas1/miyawaki/projects/002/data/proc/mpi/%s/', par.lat_interp);
    printname = [foldername 'fluxes.mat'];
    if ~exist(foldername, 'dir')
        mkdir(foldername)
    end
    save(printname, 'fluxez', 'lat');
    % save rcae data
    foldername = sprintf('/project2/tas1/miyawaki/projects/002/data/proc/mpi/%s/eps_%g/', par.lat_interp, par.ep);
    printname = [foldername 'rcae.mat'];
    if ~exist(foldername, 'dir')
        mkdir(foldername)
    end
    save(printname, 'rcae', 'lat');

end
