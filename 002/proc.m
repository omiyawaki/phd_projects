clc; close all; clear variables;

addpath(genpath('/project2/tas1/miyawaki/matlab'));

%% set parameters
par.lat_interp = 'std'; % which latitudinal grid to interpolate to: don (donohoe, coarse), era (native ERA-Interim, fine), or std (custom, very fine)
par.lat_std = transpose(-90:0.25:90); % define standard latitude grid for 'std' interpolation
par.ep_swp = 0.3; % threshold for RCE definition. RCE is defined as where abs(R1) < ep
par.ep_cp = 0.5; % threshold for RCE definition. RCE is defined as where abs(R1) < ep
par.cpd = 1005.7; par.Rd = 287; par.g = 9.81; par.L = 2.501e6; par.a = 6357e3;
gcm_info

%% call functions
type = 'era5'; % data type to run analysis on
choose_proc(type, par)
for k=1:length(par.gcm_models); par.model = par.gcm_models{k};
    type = 'gcm';
    choose_proc(type, par)
end

for i=1:length(par.ep_swp); par.ep = par.ep_swp(i); par.ga = 1-par.ep;
    % choose_proc_ep(type, par)
    for k = 1:length(par.gcm_models); par.model = par.gcm_models{k};
        type = 'gcm';
        % choose_proc_ep(type, par)
    end
end

%% define functions
function choose_proc(type, par)
    % proc_fluxes(type, par) % calculate energy fluxes in the vertically-integrated MSE budget using ERA-Interim data
    proc_net_fluxes(type, par) % calculate net energy fluxes at TOA and surface
end
function proc_fluxes(type, par)
    if strcmp(type, 'era5')
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s', type);
    elseif strcmp(type, 'gcm')
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s', type, par.model);
    end

    load(sprintf('%s/grid.mat', prefix)) % read grid data
    load(sprintf('%s/rad.mat', prefix)) % read radiation data
    load(sprintf('%s/pe.mat', prefix)) % read hydrology data
    load(sprintf('%s/stf.mat', prefix)) % read surface turbulent flux data

    lat = par.lat_std;
    % interpolate onto std lat x lon grid
    for fn = rad_vars % interpolate to std lat
        intp.(fn{1}) = permute(rad.(fn{1}), [2 1 3]);
        intp.(fn{1}) = interp1(grid.dim2.lat, intp.(fn{1}), par.lat_std, 'spline');
        intp.(fn{1}) = permute(intp.(fn{1}), [2 1 3]);
    end
    for fn = pe_vars % interpolate to std lat
        intp.(fn{1}) = permute(pe.(fn{1}), [2 1 3]);
        intp.(fn{1}) = interp1(grid.dim2.lat, intp.(fn{1}), par.lat_std, 'spline');
        intp.(fn{1}) = permute(intp.(fn{1}), [2 1 3]);
    end
    for fn = stf_vars % interpolate to std lat
        intp.(fn{1}) = permute(stf.(fn{1}), [2 1 3]);
        intp.(fn{1}) = interp1(grid.dim2.lat, intp.(fn{1}), par.lat_std, 'spline');
        intp.(fn{1}) = permute(intp.(fn{1}), [2 1 3]);
    end

    if strcmp(type, 'era5') | strcmp(type, 'erai')
        intp.ra = intp.tsr - intp.ssr + intp.ttr - intp.str; % compute net radiative cooling from radiative fluxes
        % compute surface turbulent fluxes directly from INTP data
        % multiply by negative to define flux from surface to atmosphere as positive
        intp.stf = -( intp.sshf + intp.slhf );
    elseif strcmp(type, 'gcm')
        intp.ra = intp.rsdt - intp.rsut + intp.rsus - intp.rsds + intp.rlus - intp.rlds - intp.rlut; % calculate atmospheric radiative cooling
        intp.stf = intp.hfls + intp.hfss;
    end

    intp.res = intp.ra + intp.stf; % infer MSE tendency and flux divergence as residuals
    % compute northward MSE transport using the residual data
    tediv_t = squeeze(nanmean(intp.res, 3)); % take time average
    tediv_tz = trapz(deg2rad(grid.dim2.lon), tediv_t, 1); % zonally integrate
    vh = cumtrapz(deg2rad(lat), par.a^2*cosd(lat).*tediv_tz'); % cumulatively integrate in latitude
    intp.r1 = (intp.res)./intp.ra; % calculate nondimensional number R1 disregarding MSE budget closure
    intp.r2 = intp.stf./intp.ra; % calculate nondimensional number R2 disregarding MSE budget closure

    if strcmp(type, 'era5') | strcmp(type, 'erai')
        % take zonal averages
        for fn = {'ra', 'stf', 'sshf', 'slhf', 'res', 'r1', 'r2', 'cp', 'lsp', 'e'}
            fluxes.(fn{1}) = squeeze(nanmean(intp.(fn{1}), 1));
        end
        foldername = sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/', type, par.lat_interp);
    elseif strcmp(type, 'gcm')
        % take zonal averages
        for fn = {'ra', 'stf', 'hfls', 'hfss', 'res', 'r1', 'r2', 'prc', 'pr', 'evspsbl'}
            fluxes.(fn{1}) = squeeze(nanmean(intp.(fn{1}), 1));
        end
        foldername = sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s/', type, par.model, par.lat_interp);
    end

    % save energy flux data into mat file
    printname = [foldername 'fluxes.mat'];
    if ~exist(foldername, 'dir')
        mkdir(foldername)
    end
    save(printname, 'fluxes', 'vh', 'lat');

end
function proc_net_fluxes(type, par)
% calculates the global TOA energy imbalance using ERA-Interim data
    if strcmp(type, 'era5') | strcmp(type, 'erai')
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s', type);
    elseif strcmp(type, 'gcm')
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s', type, par.model);
    end
    load(sprintf('%s/grid.mat', prefix)); % read grid data
    load(sprintf('%s/rad.mat', prefix)); % read radiation data
    load(sprintf('%s/stf.mat', prefix)); % read surface turbulent flux data

    if strcmp(type, 'era5') | strcmp(type, 'erai')
        net_toa_raw = rad.tsr + rad.ttr; % compute net radiative fluxes at TOA, positive down
    elseif strcmp(type, 'gcm')
        net_toa_raw = - rad.rsut + rad.rsdt - rad.rlut;
    end

    net_toa_tz = squeeze(nanmean(nanmean( net_toa_raw, 1 ), 3))'; % take zonal and time avera5ges and transpose to have same dimensions as latitude grid
    net_toa = nansum(cosd(grid.dim2.lat).*net_toa_tz) / nansum(cosd(grid.dim2.lat));
    if strcmp(type, 'era5') | strcmp(type, 'erai')
        disp( sprintf('The net radiative imbalance in %s at TOA is %g Wm^-2.', type, net_toa) );
    elseif strcmp(type, 'gcm')
        disp( sprintf('The net radiative imbalance in %s at TOA is %g Wm^-2.', par.model, net_toa) );
    end

    if strcmp(type, 'era5') | strcmp(type, 'erai')
        net_sfc_raw = rad.ssr + rad.str + stf.sshf + stf.slhf; % compute net radiative fluxes at surface, positive down
    elseif strcmp(type, 'gcm')
        net_sfc_raw = - rad.rsus + rad.rsds - rad.rlus + rad.rlds - stf.hfss - stf.hfls; % compute net radiative fluxes at surface, positive down
    end
    net_sfc_tz = squeeze(nanmean(nanmean( net_sfc_raw, 1 ), 3))'; % take zonal and time avera5ges and transpose to have same dimensions as latitude grid
    net_sfc = nansum(cosd(grid.dim2.lat).*net_sfc_tz) / nansum(cosd(grid.dim2.lat));
    if strcmp(type, 'era5') | strcmp(type, 'erai')
        disp( sprintf('The net radiative imbalance in %s at the surface is %g Wm^-2.', type, net_sfc) );
    elseif strcmp(type, 'gcm')
        disp( sprintf('The net radiative imbalance in %s at the surface is %g Wm^-2.', par.model, net_sfc) );
    end
end
function proc_era5_srfc(par)
    if strcmp(type, 'era5') | strcmp(type, 'erai')
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s', type);
    elseif strcmp(type, 'gcm')
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s', type, par.model);
    end
    load(sprintf('%s/grid.mat', type)); % read ERA5 grid data
    load(sprintf('%s/srfc.mat', type)); % read ERA5 grid data

    if strcmp(par.lat_interp, 'era5')
        lat = grid.lat;
        % load ERA5 data
        for fn = vars_srfc
            srfcz.(fn{1}) = permute(srfc.(fn{1}), [2 1]); % order dimensions to be consistent with Donohoe data (mon, lat, plev)
        end
    elseif strcmp(par.lat_interp, 'std') % interpolate all to fine standard grid
        lat = par.lat_std;
        % interpolate raw ERA5 data onto std lat x ERA5 lon grid
        for fn = vars_srfc
            srfcz.(fn{1}) = interp1(grid.lat, srfc.(fn{1}), par.lat_std, 'spline', nan); % interpolate to std lat
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

function choose_proc_ep(type, par)
    proc_rcae(type, par) % calculate RCE and RAE regimes
end
function proc_rcae(type, par)
    if strcmp(type, 'era5') | strcmp(type, 'erai')
        filename = sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/fluxes.mat', type, par.lat_interp); % read ERA5 zonally averaged fluxes
        foldername = sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/eps_%g/', type, par.lat_interp, par.ep);
    elseif strcmp(type, 'gcm')
        filename = sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s/fluxes.mat', type, par.model, par.lat_interp); % read gcm zonally averaged fluxes
        foldername = sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s/eps_%g/', type, par.model, par.lat_interp, par.ep);
    end
    if ~exist(filename); error('Data does not exist. Please run proc_fluxes.m first.'); else
        load(filename);
    end

    % identify locations of RCE using threshold epsilon (ep)
    rcae.def = zeros(size(fluxes.r1));
    rcae.def(abs(fluxes.r1) < par.ep) = 1;
    % identify locations of RAE using threshold gamma (ga)
    rcae.def(fluxes.r1 > par.ga) = -1;
    % add additional criteria for RCE that P-E>0
    rcae.pe = zeros(size(fluxes.r1));
    if strcmp(type, 'era5') | strcmp(type, 'erai')
        rcae.pe(abs(fluxes.r1) < par.ep & (fluxes.cp+fluxes.lsp+fluxes.e>0)) = 1; % note that evap is defined negative into atmosphere
    elseif strcmp(type, 'gcm')
        rcae.pe(abs(fluxes.r1) < par.ep & (fluxes.pr-fluxes.evspsbl>0)) = 1;
    end
    rcae.pe(fluxes.r1 > par.ga) = -1;
    % add additional criteria for RCE that (P_ls - E)<<1
    rcae.cp = zeros(size(fluxes.r1));
    if strcmp(type, 'era5') | strcmp(type, 'erai')
        rcae.cp(abs(fluxes.r1) < par.ep & abs(fluxes.lsp./fluxes.cp<par.ep_cp)) = 1; % note that evap is defined negative into atmosphere
    elseif strcmp(type, 'gcm')
        rcae.cp(abs(fluxes.r1) < par.ep & abs((fluxes.pr-fluxes.prc)./fluxes.prc<par.ep_cp)) = 1; % note that evap is defined negative into atmosphere
    end
    rcae.cp(fluxes.r1 > par.ga) = -1;

    % save rcae data
    printname = [foldername 'rcae.mat'];
    if ~exist(foldername, 'dir')
        mkdir(foldername)
    end
    save(printname, 'rcae', 'lat');

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
            vertz.(fn{1}) = interp1(grid.lat, vert.(fn{1}), lat);
            % order dimensions to be consistent with Vertzohoe data (mon, lat, plev)
            vertz.(fn{1}) = permute(vertz.(fn{1}), [3 1 2]);
        end
    elseif strcmp(par.lat_interp, 'era')
        lat = grid.lat;
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
            vertz.(fn{1}) = interp1(grid.lat, vert.(fn{1}), par.lat_std, 'spline', nan);
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

function proc_era5_temp(par)
    load('/project2/tas1/miyawaki/projects/002/data/read/era5/grid.mat'); % read ERA5 grid data
    temp = squeeze(ncread(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/%s/temp/%s_temp_%s.ymonmean.nc', type, type, par.era5.yr_span), 'temp')); % load temp

    if strcmp(par.lat_interp, 'era5')
        lat = grid.lat;
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
            vertz.(fn{1}) = interp1(grid.lat, vert.(fn{1}), par.lat_std, 'spline', nan);
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
function filt_era5_temp(par)
    load('/project2/tas1/miyawaki/projects/002/data/read/era5/grid.mat'); % read ERA5 grid data
    load(sprintf('/project2/tas1/miyawaki/projects/002/data/proc/era5/%s/fluxes.mat', par.lat_interp)); % load fluxes
    load(sprintf('/project2/tas1/miyawaki/projects/002/data/proc/era5/%s/eps_%g/rcae.mat', par.lat_interp, par.ep)); % load rcae data
    load(sprintf('/project2/tas1/miyawaki/projects/002/data/proc/era5/%s/vert.mat', par.lat_interp)); % load temp
    load(sprintf('/project2/tas1/miyawaki/projects/002/data/proc/era5/%s/srfc.mat', par.lat_interp)); % load surface data

    % create surface mask
    ps_3d = repmat(srfcz.sp, [1 1 size(vertz.t, 3)]);
    pa_3d = double(permute(repmat(grid.plev, [1 size(srfcz.sp)]), [2 3 1]))*10^2; % convert from hPa to Pa
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
        vert_filt.rce.(fn{1})(grid.plev*100 > min(min(ps_rce.(fn{1})))) = nan;
        vert_filt.rce.tp.(fn{1})(grid.plev*100 > min(min(ps_rce.(fn{1})(:,abs(lat)<30)))) = nan;
        vert_filt.rce.nh.(fn{1})(grid.plev*100 > min(min(ps_rce.(fn{1})(:,lat>30)))) = nan;
        vert_filt.rce.sh.(fn{1})(grid.plev*100 > min(min(ps_rce.(fn{1})(:,lat<-30)))) = nan;
        vert_filt.rae.(fn{1})(grid.plev*100 > min(min(ps_rae.(fn{1})))) = nan;
        vert_filt.rae.nh.(fn{1})(grid.plev*100 > min(min(ps_rae.(fn{1})(:,lat>0)))) = nan;
        vert_filt.rae.sh.(fn{1})(grid.plev*100 > min(min(ps_rae.(fn{1})(:,lat<0)))) = nan;
    end

    % save filtered data
    foldername = sprintf('/project2/tas1/miyawaki/projects/002/data/proc/era5/%s/eps_%g/', par.lat_interp, par.ep);
    printname = [foldername 'vert_filt'];
    if ~exist(foldername, 'dir')
        mkdir(foldername)
    end
    save(printname, 'vert_filt');
end

function proc_gcm_temp(par, model)

    % create surface mask
    ps_4d = permute(repmat(gcm_2d.ps, [1 1 1 size(gcm_3d.ta, 3)]), [1 2 4 3]);
    pa_4d = permute(repmat(gcm_3d.plev, [1 size(gcm_2d.ps)]), [2 3 1 4]);
    surface_mask = nan(size(gcm_3d.ta));
    surface_mask(pa_4d < ps_4d) = 1;

    % filter all 3D data with surface mask
    for fn = fieldnames(gcm_3d)'
        if ~strcmp(fn{1}, {'lon', 'lat', 'plev'})
            gcm_3d.(fn{1}) = gcm_3d.(fn{1}).*surface_mask;
        end
    end
    pa_4d = pa_4d.*surface_mask;

end
