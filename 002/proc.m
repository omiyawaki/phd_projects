clc; close all; clear variables;

addpath(genpath('/project2/tas1/miyawaki/matlab'));

%% set parameters
par.lat_interp = 'std'; % which latitudinal grid to interpolate to: don (donohoe, coarse), era (native ERA-Interim, fine), or std (custom, very fine)
par.lat_std = transpose(-90:0.25:90); % define standard latitude grid for 'std' interpolation
par.ep_swp = 0.1:0.05:0.35; % threshold for RCE definition. RCE is defined as where abs(R1) < ep
par.cpd = 1005.7; par.Rd = 287; par.g = 9.81; par.L = 2.501e6;

%% call functions
for i=1:length(par.ep_swp); par.ep = par.ep_swp(i); par.ga = 1-par.ep;
    % proc_era_rcae(par)
    proc_era_temp(par)
    filt_era_temp(par)
    % proc_mpi_rcae(par)
end

%% define functions
function proc_era_rcae(par)
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
function proc_era_temp(par)
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

function proc_mpi_rcae(par)
    load('/project2/tas1/miyawaki/projects/002/data/read/mpi_3d_climatology'); % read 3D mpi data

    lat = mpi_3d.lat;

    % calculate x component of MSE flux divergence
    lon_4d = repmat(mpi_3d.lon, [1 size(mpi_3d.ta, 2) size(mpi_3d.ta, 3) size(mpi_3d.ta, 4)]);
    lon_half = (mpi_3d.lon(2:end) + mpi_3d.lon(1:end-1))/2;
    lon_half(end+1) = (mpi_3d.lon(end)+360)/2;
    duh_half = (mpi_3d.uh(2:end,:,:,:) - mpi_3d.uh(1:end-1,:,:,:))./(lon_4d(2:end,:,:,:) - lon_4d(1:end-1,:,:,:));
    duh_half(end+1,:,:,:) = (mpi_3d.uh(1,:,:,:) - mpi_3d.uh(end,:,:,:))/(360 - mpi_3d.lon(end));
    duh = interp1(lon_half, duh_half, mpi_3d.lon);

    % calculate y component of MSE flux divergence
    lat_4d = permute(repmat(mpi_3d.lat, [1 size(mpi_3d.ta, 1) size(mpi_3d.ta, 3) size(mpi_3d.ta, 4)]), [2 1 3 4]);
    lat_half = (mpi_3d.lat(2:end) + mpi_3d.lat(1:end-1))/2;
    dvh_half = (mpi_3d.vh(:,2:end,:,:) - mpi_3d.vh(:,1:end-1,:,:))./(lat_4d(:,2:end,:,:) - lat_4d(:,1:end-1,:,:));
    dvh_half = permute(dvh_half, [2 1 3 4]);
    dvh = interp1(lat_half, dvh_half, mpi_3d.lat);
    dvh = permute(dvh, [2 1 3 4]);

    % vertically integrate
    fluxes.TEDIV = squeeze( trapz(mpi_3d.plev, -1/par.g*(duh + dvh), 3) ); % total horizontal MSE flux divergence
    h = squeeze( trapz(mpi_3d.plev, -1/par.g*(mpi_3d.h), 3) ); % MSE

    % MSE tendency
    dhdt_half = nan(size(h));
    dhdt_half(:,:,2:end) = (h(:,:,2:end) - h(:,:,1:end-1)) / (30.5*86400);
    dhdt_half(:,:,1) = (h(:,:,1) - h(:,:,end)) / (30.5*86400);
    dhdt_half = permute(dhdt_half, [3 1 2]);
    fluxes.TETEN = interp1(0.5:11.5, dhdt_half, 1:12);
    fluxes.TETEN = permute(fluxes.TETEN, [2 3 1]);

    load('/project2/tas1/miyawaki/projects/002/data/read/mpi_2d_climatology'); % read 2D mpi data
    fluxes.slhf = mpi_2d.hfls; fluxes.sshf = mpi_2d.hfss; % rename turbulent fluxes to be consistent with era
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
