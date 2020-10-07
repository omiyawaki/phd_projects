clc; close all; clear variables;

addpath(genpath('/project2/tas1/miyawaki/matlab'));

%% set parameters
% par.erai.yr_span = '2000_2012'; % spanning years for ERA-Interim
par.erai.yr_span = '2000_2018'; % spanning years for ERA-Interim
par.era5.yr_span = '1979_2019'; % spanning years for ERA5
par.gcm.clim = 'piControl'; % choose either piControl or abrupt4xCO2
par.echam.clim = '20170908'; % choose from 20170908 (snowball) or 20170915_2 (modern)
par.lat_interp = 'native'; % which latitudinal grid to interpolate to: native (no interpolation), don (donohoe, coarse), era (native ERA-Interim, fine), or std (custom, very fine)
par.lat_std = transpose(-90:0.25:90); % define standard latitude grid for 'std' interpolation
par.ep_swp = 0.1; %[0.25 0.3 0.35]; % threshold for RCE definition. RCE is defined as where abs(R1) < ep
par.ga_swp = 0.9; % optional threshold for RAE. If undefined, the default value is 1-par.ep
par.ep_cp = 0.5; % additional flag for RCE definition using convective precipitation. RCE is defined as where lsp/cp < ep_cp
par.ma_type = 'reversible'; % choose the type of moist adiabat: reversible, pseudo, or std
par.frz = 0; % consider latent heat of fusion in moist adiabat?
par.pa_span = [1000 10]*100; % pressure range for calculating moist adiabat
par.pa = linspace(1000,10,100)*1e2; % high resolution vertical grid to interpolate to
par.si = 1e-5*par.pa; % high resolution vertical grid to interpolate to
par.dpa = -10; % pressure increment for integrating dry adiabat section of moist adiabat (matters for how accurately the LCL is computed)
par.z_span = [0 25]*10^3; % height range for calculating moist adiabat
par.dz = 10; % pressure increment for integrating dry adiabat section of moist adiabat (matters for how accurately the LCL is computed)
par.do_surf = 0; % whether or not to calculate temperature profile in pressure grids including ts and ps data
par.si_eval = [0.8 0.85 0.9]; % sigma level for evaluating inversion strength (T(si_eval) - T(surface))
par.si_bl = 0.85; % sigma level to separate vertical average for close to moist adiabatic and stable surface stratifications
par.si_up = 0.4; % sigma level for upper boundary of vertical average for close to moist adiabatic
% par.era.fw = {'mse', 'dse', 'db13', 'db13s', 'db13t', 'div', 'divt', 'div79'};
% par.era.fw = {'div79', 'mse', 'dse', 'db13', 'db13s', 'db13t', 'div', 'divt'};
par.era.fw = {'div00'};
par.gcm.fw = {'mse', 'dse'};
par.echam.fw = {'mse', 'dse'};
par.cpd = 1005.7; par.cpv = 1870; par.cpl = 4186; par.cpi = 2108; par.Rd = 287; par.Rv = 461; par.g = 9.81; par.L = 2.501e6; par.a = 6357e3; par.eps = 0.622; % common constants, all in SI units for brevity
gcm_info

%% call functions
% comp_flux(par)
% ceres_flux(par)
% choose_disp(par)

type = 'echam_ml'; % data type to run analysis on
% choose_proc(type, par)
for k=1:length(par.gcm_models); par.model = par.gcm_models{k};
    type = 'gcm';
    choose_proc(type, par)
end

for i=1:length(par.ep_swp); par.ep = par.ep_swp(i); par.ga = par.ga_swp(i);
    type = 'erai';
    % choose_proc_ep(type, par)
    for k = 1:length(par.gcm_models); par.model = par.gcm_models{k};
        type = 'gcm';
        % choose_proc_ep(type, par)
    end
end

%% define functions
function choose_proc(type, par)
    % save_mask(type, par) % save land and ocean masks once (faster than creating mask every time I need it)
    % proc_flux(type, par) % calculate energy fluxes in the vertically-integrated MSE budget using ERA-Interim data
    % proc_net_flux(type, par) % calculate net energy fluxes at TOA and surface
    % proc_temp_mon_lat(type, par) % calculate mon x lat temperature profiles
    % proc_ma_mon_lat(type, par) % calculate mon x lat moist adiabats
    % proc_inv_str_mon_lat(type, par) % calculate mon x lat inversion strength
    % proc_ga_diff_si_mon_lat(type, par) % calculate mon x lat gamma percentage difference
    % proc_ga_bl_diff_si_mon_lat(type, par) % calculate mon x lat gamma percentage difference
    % proc_ga_malr_diff_si_mon_lat(type, par) % calculate mon x lat gamma percentage difference
    proc_ga_dalr_bl_diff_si_mon_lat(type, par) % calculate mon x lat gamma percentage difference
    % proc_ga_trop_malr_diff_si_mon_lat(type, par) % calculate mon x lat gamma percentage difference
    % make_tai(type, par) % calculate moist adiabat in lon x lat x mon
    % make_dtdzsi(type, par) % calculate model lapse rate and interpolate to sigma coordinates
    % make_malrsi(type, par) % calculate moist adiabatic lapse rate of model temperature sigma coordinates
end % select which functions to run at a time
function save_mask(type, par)
    if strcmp(type, 'era5') | strcmp(type, 'erai')
        foldername = sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/', type, par.lat_interp);
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s', type);
    elseif strcmp(type, 'gcm')
        foldername = sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s/%s/', type, par.model, par.gcm.clim, par.lat_interp);
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s/%s', type, par.model, par.gcm.clim);
    elseif strcmp(type, 'echam')
        foldername = sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s/', type, par.echam.clim, par.lat_interp);
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s', type, par.echam.clim);
    elseif strcmp(type, 'echam_ml')
        foldername = sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/', type, par.lat_interp);
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s', type);
    end

    load(sprintf('%s/grid.mat', prefix)); % read grid data

    if strcmp(par.lat_interp, 'std')
        lat = par.lat_std;
    else
        lat = grid.dim3.lat;
    end

    mask.land = remove_land(lat, grid.dim3.lon, 12);
    mask.ocean = remove_ocean(lat, grid.dim3.lon, 12);

    % save masks
    printname = [foldername 'masks'];
    if ~exist(foldername, 'dir')
        mkdir(foldername)
    end
    save(printname, 'mask');

end % compute masks once and save it for later use
function proc_flux(type, par)
    if any(strcmp(type, {'era5', 'erai'}))
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s', type);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s', type, par.lat_interp);
        % file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/%s/w500/%s_w500_%s.ymonmean.nc', type, type, par.(type).yr_span));
        % fullpath=sprintf('%s/%s', file.folder, file.name);
        % w500 = ncread(fullpath, 'w');
        % file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/%s/vas/%s_vas_%s.ymonmean.nc', type, type, par.(type).yr_span));
        % fullpath=sprintf('%s/%s', file.folder, file.name);
        % vas = ncread(fullpath, 'v10');
        don = load(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/don/radiation_dynamics_climatology')); % read donohoe data
        prefix_ceres=sprintf('/project2/tas1/miyawaki/projects/002/data/read/ceres'); % prefix for CERES data
        load(sprintf('%s/div.mat', prefix)) % read divergence data
        load(sprintf('%s/tend.mat', prefix)) % read tendency data
        if any(strcmp(type, {'erai'})) & strcmp(par.erai.yr_span, '1979_2018')
            don79 = load(sprintf('%s/dondiv79.mat', prefix)); % read Donohoe data 1979--2018
        elseif any(strcmp(type, {'erai'})) & strcmp(par.erai.yr_span, '2000_2018')
            don79 = load(sprintf('%s/dondiv00.mat', prefix)); % read Donohoe data 1979--2018
        end
    elseif strcmp(type, 'gcm')
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s/%s', type, par.model, par.gcm.clim);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s/%s', type, par.model, par.gcm.clim, par.lat_interp);
        % file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/gcm/%s/%s_Amon_%s_%s_r1i1p1_*.ymonmean.nc', par.model, 'w500', par.model, par.gcm.clim)); % load w500
        % fullpath=sprintf('%s/%s', file.folder, file.name);
        % w500 = squeeze(ncread(fullpath, 'wap'));
        % file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/gcm/%s/%s_Amon_%s_%s_r1i1p1_*.ymonmean.nc', par.model, 'vas', par.model, par.gcm.clim)); % load vas
        % fullpath=sprintf('%s/%s', file.folder, file.name);
        % vas = squeeze(ncread(fullpath, 'vas'));
    elseif strcmp(type, 'echam')
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s', type, par.echam.clim);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s', type, par.echam.clim, par.lat_interp);
    end

    load(sprintf('%s/grid.mat', prefix)) % read grid data
    load(sprintf('%s/rad.mat', prefix)) % read radiation data
    load(sprintf('%s/pe.mat', prefix)) % read hydrology data
    load(sprintf('%s/stf.mat', prefix)) % read surface turbulent flux data
    load(sprintf('%s/masks.mat', prefix_proc)); % load land and ocean masks

    if strcmp(par.lat_interp, 'std')
        lat = par.lat_std;
    else
        lat = grid.dim3.lat;
    end

    % interpolate onto std lat x lon grid
    if any(strcmp(type, {'era5', 'erai'}))
        ceres = struct();
        tmp = load(sprintf('%s/grid.mat', prefix_ceres)); ceres.grid = tmp.grid; % read grid data
        tmp = load(sprintf('%s/rad.mat', prefix_ceres)); ceres.rad = tmp.rad; % read CERES radiation data
        for fn = {'TEDIV', 'TETEN'}; fname = fn{1};
            flux.(fname) = permute(don.(fname), [2 3 1]); % bring lat forward
            flux.(fname) = interp1(don.lat, flux.(fname), lat, 'linear');
            flux.(fname) = permute(flux.(fname), [2 1 3]); % bring lon forward
            flux.(fname) = interp1(don.lon, flux.(fname), grid.dim2.lon, 'linear');
        end
        flux.don79div = interp1(don79.lat, don79.dondiv, lat, 'linear'); % interpolate to standard lat
        flux.don79div = repmat(flux.don79div, [1 1 length(grid.dim2.lon)]); % repeat to longitude because don79 data is already zonally averaged
        flux.don79div = permute(flux.don79div, [3 1 2]); % bring lon to 1st
        for fn = fieldnames(ceres.rad)'; fname = fn{1}; % interpolate to std lat and ERA lon
            ceres.(fname) = interp1(ceres.grid.dim2.lon, ceres.rad.(fname), grid.dim2.lon, 'linear');
            ceres.(fname) = permute(ceres.(fname), [2 1 3]); % bring lat forward
            ceres.(fname) = interp1(ceres.grid.dim2.lat, ceres.(fname), lat, 'linear');
            ceres.(fname) = permute(ceres.(fname), [2 1 3]); % bring lon forward
        end
        ceres.tsr = ceres.tsdr - ceres.tsur; ceres.net = ceres.tsr - ceres.ttr;
        ceres.str = -ceres.str;
        ceres.ra = ceres.tsr - ceres.ssr + ceres.str - ceres.ttr;
    end
    for fn = rad_vars; fname = fn{1}; % interpolate to std lat
        flux.(fname) = permute(rad.(fname), [2 1 3]);
        flux.(fname) = interp1(grid.dim2.lat, flux.(fname), lat, 'linear');
        flux.(fname) = permute(flux.(fname), [2 1 3]);
    end
    for fn = pe_vars; fname = fn{1}; % interpolate to std lat
        flux.(fname) = permute(pe.(fname), [2 1 3]);
        flux.(fname) = interp1(grid.dim2.lat, flux.(fname), lat, 'linear');
        flux.(fname) = permute(flux.(fname), [2 1 3]);
    end
    for fn = stf_vars; fname = fn{1}; % interpolate to std lat
        flux.(fname) = permute(stf.(fname), [2 1 3]);
        flux.(fname) = interp1(grid.dim2.lat, flux.(fname), lat, 'linear');
        flux.(fname) = permute(flux.(fname), [2 1 3]);
    end
    if any(strcmp(type, {'era5', 'erai'}))
        for fn = tend_vars_txt; fname = fn{1}; % interpolate to std lat
            flux.(fname) = permute(tend.(fname), [2 1 3]);
            flux.(fname) = interp1(grid.dim2.lat, flux.(fname), lat, 'linear');
            flux.(fname) = permute(flux.(fname), [2 1 3]);
        end; clear tend
        for fn = div_vars_txt; fname = fn{1}; % interpolate to std lat
            flux.(fname) = permute(div.(fname), [2 1 3]);
            flux.(fname) = interp1(grid.dim2.lat, flux.(fname), lat, 'linear');
            flux.(fname) = permute(flux.(fname), [2 1 3]);
        end; clear div
    end
    % flux.w500 = permute(w500, [2 1 3]);
    % flux.w500 = interp1(grid.dim3.lat, flux.w500, lat, 'linear');
    % flux.w500 = permute(flux.w500, [2 1 3]);
    % flux.vas = permute(vas, [2 1 3]);
    % flux.vas = interp1(grid.dim3.lat, flux.vas, lat, 'linear');
    % flux.vas = permute(flux.vas, [2 1 3]);

    if strcmp(type, 'era5') | strcmp(type, 'erai')
        % compute surface turbulent fluxes directly from INTP data
        % multiply by negative to define flux from surface to atmosphere as positive
        flux.stf.mse = -( flux.sshf + flux.slhf ); flux.stf.mse2 = flux.stf.mse;
        flux.stf.dse = par.L*(flux.cp+flux.lsp) - flux.sshf;
    elseif strcmp(type, 'gcm')
        flux.stf.mse = flux.hfls + flux.hfss; flux.stf.mse2 = flux.stf.mse;
        flux.stf.dse = par.L*flux.pr + flux.hfss;
    elseif strcmp(type, 'echam')
        flux.stf.mse = -(flux.ahfl + flux.ahfs); flux.stf.mse2 = flux.stf.mse;
        flux.stf.dse = par.L*(flux.aprc+flux.aprl) - flux.ahfs;
    end

    if any(strcmp(type, {'era5', 'erai'})); f_vec = par.era.fw;
    elseif strcmp(type, 'gcm'); f_vec = par.gcm.fw;
    elseif strcmp(type, 'echam'); f_vec = par.echam.fw; end
    for f = f_vec; fw = f{1};
        if any(strcmp(type, {'era5', 'erai'}));
            flux.lw = flux.ttr - flux.str; flux.sw = flux.tsr-flux.ssr; % compute net shortwave and longwave flux through atmosphere
            if any(strcmp(fw, {'mse', 'dse'})); flux.ra.(fw) = flux.tsr - flux.ssr + flux.ttr - flux.str; % compute net radiative cooling from radiative fluxes
            elseif contains(fw, 'db13') | contains(fw, 'div'); flux.ra.(fw) = ceres.ra; end % use radiative cooling from CERES data
        elseif strcmp(type, 'gcm');
            flux.lw = flux.rlus - flux.rlds - flux.rlut; flux.sw = flux.rsdt - flux.rsut + flux.rsus - flux.rsds;
            flux.ra.(fw) = flux.rsdt - flux.rsut + flux.rsus - flux.rsds + flux.rlus - flux.rlds - flux.rlut;
        elseif strcmp(type, 'echam');
            flux.lw = flux.trad0 - flux.trads; flux.sw = flux.srad0 - flux.srads;
            flux.ra.(fw) = flux.lw + flux.sw;
        end % calculate atmospheric radiative cooling
        if any(strcmp(fw, {'mse', 'dse'}))
            flux.res.(fw) = flux.ra.(fw) + flux.stf.(fw); % infer MSE tendency and flux divergence as residuals
        elseif any(strcmp(fw, {'mse2'}))
            flux.res.(fw) = flux.lw + flux.stf.(fw);
        elseif strcmp(fw, 'db13')
            flux.res.(fw) = flux.TEDIV + flux.TETEN; % use MSE tendency and flux divergence from DB13
            flux.stf.(fw) = flux.res.(fw) - flux.ra.(fw); % infer the surface turbulent fluxes
        elseif strcmp(fw, 'db13s')
            flux.res.(fw) = flux.TEDIV; % use MSE flux divergence from DB13, ignore MSE tendency term
            flux.stf.(fw) = flux.res.(fw) - flux.ra.(fw); % infer the surface turbulent fluxes
        elseif strcmp(fw, 'db13t')
            flux.res.(fw) = flux.TEDIV + flux.tend; % use MSE flux divergence from DB13, ignore MSE tendency term
            flux.stf.(fw) = flux.res.(fw) - flux.ra.(fw); % infer the surface turbulent fluxes
        elseif strcmp(fw, 'div')
            flux.res.(fw) = flux.divt + flux.divg + flux.divq*par.L; % use MSE tendency and flux divergence from ERA5 output
            flux.stf.(fw) = flux.res.(fw) - flux.ra.(fw); % infer the surface turbulent fluxes
        elseif strcmp(fw, 'divt')
            flux.res.(fw) = flux.divt + flux.divg + flux.divq*par.L + flux.tend; % use MSE tendency and flux divergence from ERA5 output
            flux.stf.(fw) = flux.res.(fw) - flux.ra.(fw); % infer the surface turbulent fluxes
        elseif strcmp(fw, 'div79')
            flux.res.(fw) = flux.don79div + flux.tend; % use MSE tendency from ERA output and MSE flux divergence from Donohoe MSE flux transport data
            flux.stf.(fw) = flux.res.(fw) - flux.ra.(fw); % infer the surface turbulent fluxes
        elseif strcmp(fw, 'div00')
            flux.res.(fw) = flux.don79div + flux.tend; % use MSE tendency from ERA output and MSE flux divergence from Donohoe MSE flux transport data
            flux.stf.(fw) = flux.res.(fw) - ceres.ra; % infer the surface turbulent fluxes
        end
        if strcmp(fw, 'mse2')
            flux.r1.(fw) = (flux.res.(fw))./flux.lw; % calculate nondimensional number R1 disregarding MSE budget closure
            flux.r2.(fw) = flux.stf.(fw)./flux.lw; % calculate nondimensional number R2 disregarding MSE budget closure
        else
            flux.r1.(fw) = (flux.res.(fw))./flux.ra.(fw); % calculate nondimensional number R1 disregarding MSE budget closure
            flux.r2.(fw) = flux.stf.(fw)./flux.ra.(fw); % calculate nondimensional number R2 disregarding MSE budget closure
        end
        if any(strcmp(type, {'era5', 'erai'}));
            flux.ftoa.(fw) = flux.tsr + flux.ttr; flux.fsfc.(fw) = -flux.ssr - flux.str + flux.stf.(fw);
        elseif strcmp(type, 'gcm');
            flux.ftoa.(fw) = flux.rsdt - flux.rsut - flux.rlut;
            flux.fsfc.(fw) = flux.rsus - flux.rsds + flux.rlus - flux.rlds + flux.stf.(fw);
        elseif strcmp(type, 'echam');
            flux.ftoa.(fw) = flux.trad0 + flux.srad0;
            flux.fsfc.(fw) = -flux.trads - flux.srads + flux.stf.(fw);
        end
    end

    if strcmp(type, 'era5') | strcmp(type, 'erai')
        for fn = {'sshf', 'slhf', 'cp', 'lsp', 'e', 'lw', 'sw', 'tend', 'divt', 'divg', 'divq', 'TETEN', 'TEDIV', 'don79div'}; fname = fn{1};
            for l = {'lo', 'l', 'o'}; land = l{1};
                if strcmp(land, 'lo'); flux_n.(land).(fname) = flux.(fname);
                elseif strcmp(land, 'l'); flux_n.(land).(fname) = flux.(fname) .*mask.ocean;
                elseif strcmp(land, 'o'); flux_n.(land).(fname) = flux.(fname) .*mask.land;
                end

                % take zonal averages
                flux_z.(land).(fname) = squeeze(nanmean(flux_n.(land).(fname), 1));

                % take time averages
                for t = {'ann', 'djf', 'jja', 'mam', 'son'}; time = t{1};
                    if strcmp(time, 'ann')
                        flux_t.(land).(time).(fname) = nanmean(flux_n.(land).(fname), 3);
                    elseif strcmp(time, 'djf')
                        flux_shift.(land).(fname) = circshift(flux_n.(land).(fname), 1, 3);
                        flux_t.(land).(time).(fname) = nanmean(flux_shift.(land).(fname)(:,:,1:3), 3);
                    elseif strcmp(time, 'jja')
                        flux_t.(land).(time).(fname) = nanmean(flux_n.(land).(fname)(:,:,6:8), 3);
                    elseif strcmp(time, 'mam')
                        flux_t.(land).(time).(fname) = nanmean(flux_n.(land).(fname)(:,:,3:5), 3);
                    elseif strcmp(time, 'son')
                        flux_t.(land).(time).(fname) = nanmean(flux_n.(land).(fname)(:,:,9:11), 3);
                    end
                    flux_zt.(land).(time).(fname) = squeeze(nanmean(flux_t.(land).(time).(fname), 1));
                end
            end
        end
        foldername = sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/', type, par.lat_interp);
    elseif strcmp(type, 'gcm')
        for fn = {'hfls', 'hfss', 'prc', 'pr', 'evspsbl', 'lw', 'sw'}; fname = fn{1};
            for l = {'lo', 'l', 'o'}; land = l{1};
                if strcmp(land, 'lo'); flux_n.(land).(fname) = flux.(fname);
                elseif strcmp(land, 'l'); flux_n.(land).(fname) = flux.(fname) .*mask.ocean;
                elseif strcmp(land, 'o'); flux_n.(land).(fname) = flux.(fname) .*mask.land;
                end

                % take zonal averages
                flux_z.(land).(fname) = squeeze(nanmean(flux_n.(land).(fname), 1));

                % take time averages
                for t = {'ann', 'djf', 'jja', 'mam', 'son'}; time = t{1};
                    if strcmp(time, 'ann')
                        flux_t.(land).(time).(fname) = nanmean(flux_n.(land).(fname), 3);
                    elseif strcmp(time, 'djf')
                        flux_shift.(land).(fname) = circshift(flux_n.(land).(fname), 1, 3);
                        flux_t.(land).(time).(fname) = nanmean(flux_shift.(land).(fname)(:,:,1:3), 3);
                    elseif strcmp(time, 'jja')
                        flux_t.(land).(time).(fname) = nanmean(flux_n.(land).(fname)(:,:,6:8), 3);
                    elseif strcmp(time, 'mam')
                        flux_t.(land).(time).(fname) = nanmean(flux_n.(land).(fname)(:,:,3:5), 3);
                    elseif strcmp(time, 'son')
                        flux_t.(land).(time).(fname) = nanmean(flux_n.(land).(fname)(:,:,9:11), 3);
                    end
                    flux_zt.(land).(time).(fname) = squeeze(nanmean(flux_t.(land).(time).(fname), 1));
                end
            end
        end
        foldername = sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s/%s/', type, par.model, par.gcm.clim, par.lat_interp);
    elseif strcmp(type, 'echam')
        for fn = {'ahfl', 'ahfs', 'aprc', 'aprl', 'evap', 'lw', 'sw'}; fname = fn{1};
            for l = {'lo', 'l', 'o'}; land = l{1};
                if strcmp(land, 'lo'); flux_n.(land).(fname) = flux.(fname);
                elseif strcmp(land, 'l'); flux_n.(land).(fname) = flux.(fname) .*mask.ocean;
                elseif strcmp(land, 'o'); flux_n.(land).(fname) = flux.(fname) .*mask.land;
                end

                % take zonal averages
                flux_z.(land).(fname) = squeeze(nanmean(flux_n.(land).(fname), 1));

                % take time averages
                for t = {'ann', 'djf', 'jja', 'mam', 'son'}; time = t{1};
                    if strcmp(time, 'ann')
                        flux_t.(land).(time).(fname) = nanmean(flux_n.(land).(fname), 3);
                    elseif strcmp(time, 'djf')
                        flux_shift.(land).(fname) = circshift(flux_n.(land).(fname), 1, 3);
                        flux_t.(land).(time).(fname) = nanmean(flux_shift.(land).(fname)(:,:,1:3), 3);
                    elseif strcmp(time, 'jja')
                        flux_t.(land).(time).(fname) = nanmean(flux_n.(land).(fname)(:,:,6:8), 3);
                    elseif strcmp(time, 'mam')
                        flux_t.(land).(time).(fname) = nanmean(flux_n.(land).(fname)(:,:,3:5), 3);
                    elseif strcmp(time, 'son')
                        flux_t.(land).(time).(fname) = nanmean(flux_n.(land).(fname)(:,:,9:11), 3);
                    end
                    flux_zt.(land).(time).(fname) = squeeze(nanmean(flux_t.(land).(time).(fname), 1));
                end
            end
        end
        foldername = sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s/', type, par.echam.clim, par.lat_interp);
    end
    for fn = {'ra', 'stf', 'res', 'r1', 'r2', 'ftoa', 'fsfc'}; fname = fn{1};
        if any(strcmp(type, {'era5', 'erai'})); f_vec = par.era.fw;
        elseif strcmp(type, 'gcm'); f_vec = par.gcm.fw;
        elseif strcmp(type, 'echam'); f_vec = par.echam.fw; end
        for f = f_vec; fw = f{1};
            for l = {'lo', 'l', 'o'}; land = l{1};
                if strcmp(land, 'lo'); flux_n.(land).(fname).(fw) = flux.(fname).(fw);
                elseif strcmp(land, 'l'); flux_n.(land).(fname).(fw) = flux.(fname).(fw) .*mask.ocean;
                elseif strcmp(land, 'o'); flux_n.(land).(fname).(fw) = flux.(fname).(fw) .*mask.land;
                end

                if strcmp(fname, 'res')
                    % compute northward MSE transport using the residual data
                    tediv_0 = fillmissing(flux_n.(land).res.(fw), 'constant', 0); % replace nans with 0s so missing data doesn't influence transport
                    % tediv_z = squeeze(trapz(deg2rad(grid.dim2.lon), tediv_0, 1)); % zonally integrate
                    tediv_z = squeeze(nanmean(tediv_0, 1)); % zonally average
                    vh_mon.(land).(fw) = cumtrapz(deg2rad(lat), 2*pi*par.a^2*cosd(lat).*tediv_z, 1); % cumulatively integrate
                    for t = {'ann', 'djf', 'jja', 'mam', 'son'}; time = t{1};
                        if strcmp(time, 'ann')
                            tediv_t = squeeze(nanmean(flux_n.(land).res.(fw), 3)); % take time average
                        elseif strcmp(time, 'djf')
                            flux_n_shift.(land).res = circshift(flux_n.(land).res.(fw), 1, 3); % shift month by 1 to allow for djf average
                            tediv_t = squeeze(nanmean(flux_n_shift.(land).res(:,:,1:3), 3)); % take time average
                        elseif strcmp(time, 'jja')
                            tediv_t = squeeze(nanmean(flux_n.(land).res.(fw)(:,:,6:8), 3)); % take time average
                        elseif strcmp(time, 'mam')
                            tediv_t = squeeze(nanmean(flux_n.(land).res.(fw)(:,:,3:5), 3)); % take time average
                        elseif strcmp(time, 'son')
                            tediv_t = squeeze(nanmean(flux_n.(land).res.(fw)(:,:,9:11), 3)); % take time average
                        end
                        % tediv_tz = trapz(deg2rad(grid.dim2.lon), tediv_t, 1); % zonally integrate
                        tediv_tz = squeeze(nanmean(tediv_t, 1)); % zonally integrate
                        tediv_tz = fillmissing(tediv_tz, 'constant', 0);
                        vh.(land).(time).(fw) = cumtrapz(deg2rad(lat), 2*pi*par.a^2*cosd(lat).*tediv_tz',1); % cumulatively integrate in latitude
                    end
                end

                % take zonal average
                flux_z.(land).(fname).(fw) = squeeze(nanmean(flux_n.(land).(fname).(fw), 1));

                % take time averages
                for t = {'ann', 'djf', 'jja', 'mam', 'son'}; time = t{1};
                    if strcmp(time, 'ann')
                        flux_t.(land).(time).(fname).(fw) = nanmean(flux_n.(land).(fname).(fw), 3);
                    elseif strcmp(time, 'djf')
                        flux_shift.(land).(fname).(fw) = circshift(flux_n.(land).(fname).(fw), 1, 3);
                        flux_t.(land).(time).(fname).(fw) = nanmean(flux_shift.(land).(fname).(fw)(:,:,1:3), 3);
                    elseif strcmp(time, 'jja')
                        flux_t.(land).(time).(fname).(fw) = nanmean(flux_n.(land).(fname).(fw)(:,:,6:8), 3);
                    elseif strcmp(time, 'mam')
                        flux_t.(land).(time).(fname).(fw) = nanmean(flux_n.(land).(fname).(fw)(:,:,3:5), 3);
                    elseif strcmp(time, 'son')
                        flux_t.(land).(time).(fname).(fw) = nanmean(flux_n.(land).(fname).(fw)(:,:,9:11), 3);
                    end
                    flux_zt.(land).(time).(fname).(fw) = squeeze(nanmean(flux_t.(land).(time).(fname).(fw), 1));
                end
            end
        end
    end

    % save energy flux data into mat file
    if ~exist(foldername, 'dir')
        mkdir(foldername)
    end
    for v = {'flux', 'flux_t', 'flux_z', 'flux_zt', 'vh', 'vh_mon'}; varname = v{1}; % removed flux and flux_t because it takes forever
        save(sprintf('%s%s', foldername, varname), varname, 'lat', '-v7.3');
    end

end % process radiative fluxes into one struct
function proc_net_flux(type, par)
% calculates the global TOA energy imbalance using ERA-Interim data
    if strcmp(type, 'era5') | strcmp(type, 'erai')
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s', type);
    elseif strcmp(type, 'gcm')
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s/%s', type, par.model, par.gcm.clim);
    elseif strcmp(type, 'echam')
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s', type, par.echam.clim);
    end
    load(sprintf('%s/grid.mat', prefix)); % read grid data
    load(sprintf('%s/rad.mat', prefix)); % read radiation data
    load(sprintf('%s/stf.mat', prefix)); % read surface turbulent flux data

    if strcmp(type, 'era5') | strcmp(type, 'erai')
        net_toa_raw = rad.tsr + rad.ttr; % compute net radiative fluxes at TOA, positive down
    elseif strcmp(type, 'gcm')
        net_toa_raw = - rad.rsut + rad.rsdt - rad.rlut;
    elseif strcmp(type, 'echam')
        net_toa_raw = rad.trad0 + rad.srad0;
    end

    net_toa_tz = squeeze(nanmean(nanmean( net_toa_raw, 1 ), 3))'; % take zonal and time avera5ges and transpose to have same dimensions as latitude grid
    net_toa = nansum(cosd(grid.dim2.lat).*net_toa_tz) / nansum(cosd(grid.dim2.lat));
    if strcmp(type, 'era5') | strcmp(type, 'erai')
        disp( sprintf('The net radiative imbalance in %s at TOA is %g Wm^-2.', type, net_toa) );
    elseif strcmp(type, 'gcm')
        disp( sprintf('The net radiative imbalance in %s at TOA is %g Wm^-2.', par.model, net_toa) );
    elseif strcmp(type, 'echam')
        disp( sprintf('The net radiative imbalance in ECHAM at TOA is %g Wm^-2.', net_toa) );
    end

    if strcmp(type, 'era5') | strcmp(type, 'erai')
        net_sfc_raw = rad.ssr + rad.str + stf.sshf + stf.slhf; % compute net radiative fluxes at surface, positive down
    elseif strcmp(type, 'gcm')
        net_sfc_raw = - rad.rsus + rad.rsds - rad.rlus + rad.rlds - stf.hfss - stf.hfls; % compute net radiative fluxes at surface, positive down
    elseif strcmp(type, 'echam')
        net_sfc_raw = rad.srads + rad.trads + stf.ahfl + stf.ahfs; % compute net radiative fluxes at surface, positive down
    end
    net_sfc_tz = squeeze(nanmean(nanmean( net_sfc_raw, 1 ), 3))'; % take zonal and time avera5ges and transpose to have same dimensions as latitude grid
    net_sfc = nansum(cosd(grid.dim2.lat).*net_sfc_tz) / nansum(cosd(grid.dim2.lat));
    if strcmp(type, 'era5') | strcmp(type, 'erai')
        disp( sprintf('The net radiative imbalance in %s at the surface is %g Wm^-2.', type, net_sfc) );
    elseif strcmp(type, 'gcm')
        disp( sprintf('The net radiative imbalance in %s at the surface is %g Wm^-2.', par.model, net_sfc) );
    elseif strcmp(type, 'echam')
        disp( sprintf('The net radiative imbalance in ECHAM at the surface is %g Wm^-2.', net_sfc) );
    end
end
function proc_temp_mon_lat(type, par)
    if strcmp(type, 'era5') | strcmp(type, 'erai')
        foldername = sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/', type, par.lat_interp);
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s', type);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s', type);
        file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/%s/temp/%s_temp_%s.ymonmean.nc', type, type, par.(type).yr_span));
        fullpath=sprintf('%s/%s', file.folder, file.name);
        ta_orig = ncread(fullpath, 't');
    elseif strcmp(type, 'gcm')
        foldername = sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s/%s/', type, par.model, par.gcm.clim, par.lat_interp);
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s/%s', type, par.model, par.gcm.clim);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s', type, par.model, par.gcm.clim);
        var = 'ta';
        file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/gcm/%s/%s_Amon_%s_%s_r1i1p1_*.ymonmean.nc', par.model, var, par.model, par.gcm.clim));
        fullpath=sprintf('%s/%s', file.folder, file.name);
        ta_orig = ncread(fullpath, var);
    elseif strcmp(type, 'echam')
        foldername = sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s/', type, par.echam.clim, par.lat_interp);
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s', type, par.echam.clim);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s', type, par.echam.clim);
        var = 't';
        file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/echam/ATM_rjg_%s_*.ymonmean.nc', par.echam.clim));
        fullpath=sprintf('%s/%s', file.folder, file.name);
        ta_orig = double(ncread(fullpath, var));
    elseif strcmp(type, 'echam_ml')
        foldername = sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/', type, par.lat_interp);
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s', type);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s', type);
        var = 't';
        file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/echam_ml/ATM_*.ymonmean.nc'));
        fullpath=sprintf('%s/%s', file.folder, file.name);
        ta_orig = double(ncread(fullpath, var));
    end
    load(sprintf('%s/grid.mat', prefix)); % read grid data
    load(sprintf('%s/tempz.mat', prefix)); taz_orig = tempz; clear tempz; % read temp in z coordinates
    load(sprintf('%s/tempsi.mat', prefix));tasi_orig = tempsi; clear tempsi; % read temp in si coordinates
    load(sprintf('%s/srfc.mat', prefix)); % load surface data
    load(sprintf('%s/%s/masks.mat', prefix_proc, par.lat_interp)); % load land and ocean masks

    lat = par.lat_std;

    % interpolate ta to standard lat grid
    ta_orig = permute(ta_orig, [2 1 3 4]);
    ta_orig = interp1(grid.dim3.lat, ta_orig, lat);
    ta_orig = permute(ta_orig, [2 1 3 4]);
    taz_orig = permute(taz_orig, [2 1 3 4]);
    taz_orig = interp1(grid.dim3.lat, taz_orig, lat);
    taz_orig = permute(taz_orig, [2 1 3 4]);
    tasi_orig = permute(tasi_orig, [2 1 3 4]);
    tasi_orig = interp1(grid.dim3.lat, tasi_orig, lat);
    tasi_orig = permute(tasi_orig, [2 1 3 4]);

    % create surface mask
    if strcmp(type, 'era5') | strcmp(type, 'erai')
        ps = permute(srfc.sp, [2 1 3]); % bring lat to front to interpolate
        ps = interp1(grid.dim2.lat, ps, lat, 'linear', 'extrap'); % interpolate to standard grid
        ps = permute(ps, [2 1 3]); % reorder to original
        ps_vert = repmat(ps, [1 1 1 size(ta_orig, 3)]); % dims (lon x lat x time x plev)
        ps_vert = permute(ps_vert, [1 2 4 3]); % dims (lon x lat x plev x time)
        pa = double(permute(repmat(grid.dim3.plev, [1 size(ps)]), [2 3 1 4]));
        zs = permute(srfc.zs, [2 1 3]); % bring lat to front to interpolate
        zs = interp1(grid.dim2.lat, zs, lat); % interpolate to standard grid
        zs = permute(zs, [2 1 3]); % reorder to original
        zs_vert = repmat(zs, [1 1 1 size(taz_orig, 3)]); % dims (lon x lat x time x plev)
        zs_vert = permute(zs_vert, [1 2 4 3]); % dims (lon x lat x plev x time)
        za = permute(repmat(grid.dim3.z, [1 size(zs)]), [2 3 1 4]);
        ts = permute(srfc.t2m, [2 1 3]); % bring lat to front to interpolate
        ts = interp1(grid.dim2.lat, ts, lat); % interpolate to standard grid
        ts = permute(ts, [2 1 3]); % reorder to original
        ts_vert = repmat(ts, [1 1 1 size(taz_orig, 3)]); % dims (lon x lat x time x plev)
        ts_vert = permute(ts_vert, [1 2 4 3]); % dims (lon x lat x plev x time)
    elseif strcmp(type, 'gcm')
        ps = permute(srfc.ps, [2 1 3]); % bring lat to front to interpolate
        ps = interp1(grid.dim2.lat, ps, lat, 'linear', 'extrap'); % interpolate to standard grid
        ps = permute(ps, [2 1 3]); % reorder to original
        ps_vert = repmat(ps, [1 1 1 size(ta_orig, 3)]); % dims (lon x lat x time x plev)
        ps_vert = permute(ps_vert, [1 2 4 3]); % dims (lon x lat x plev x time)
        pa = permute(repmat(grid.dim3.plev, [1 size(ps)]), [2 3 1 4]);
        zs = permute(srfc.zs, [2 1 3]); % bring lat to front to interpolate
        zs = interp1(grid.dim2.lat, zs, lat); % interpolate to standard grid
        zs = permute(zs, [2 1 3]); % reorder to original
        zs_vert = repmat(zs, [1 1 1 size(taz_orig, 3)]); % dims (lon x lat x time x plev)
        zs_vert = permute(zs_vert, [1 2 4 3]); % dims (lon x lat x plev x time)
        za = permute(repmat(grid.dim3.z, [1 size(zs)]), [2 3 1 4]);
        ts = permute(srfc.tas, [2 1 3]); % bring lat to front to interpolate
        ts = interp1(grid.dim2.lat, ts, lat); % interpolate to standard grid
        ts = permute(ts, [2 1 3]); % reorder to original
        ts_vert = repmat(ts, [1 1 1 size(taz_orig, 3)]); % dims (lon x lat x time x plev)
        ts_vert = permute(ts_vert, [1 2 4 3]); % dims (lon x lat x plev x time)
    elseif strcmp(type, 'echam') | strcmp(type, 'echam_ml')
        ps = permute(srfc.aps, [2 1 3]); % bring lat to front to interpolate
        ps = interp1(grid.dim2.lat, ps, lat, 'linear', 'extrap'); % interpolate to standard grid
        ps = permute(ps, [2 1 3]); % reorder to original
        ps_vert = repmat(ps, [1 1 1 size(ta_orig, 3)]); % dims (lon x lat x time x plev)
        ps_vert = permute(ps_vert, [1 2 4 3]); % dims (lon x lat x plev x time)
        pa = permute(repmat(grid.dim3.plev, [1 size(ps)]), [2 3 1 4]);
        zs = permute(srfc.zs, [2 1 3]); % bring lat to front to interpolate
        zs = interp1(grid.dim2.lat, zs, lat); % interpolate to standard grid
        zs = permute(zs, [2 1 3]); % reorder to original
        zs_vert = repmat(zs, [1 1 1 size(taz_orig, 3)]); % dims (lon x lat x time x plev)
        zs_vert = permute(zs_vert, [1 2 4 3]); % dims (lon x lat x plev x time)
        za = permute(repmat(grid.dim3.z, [1 size(zs)]), [2 3 1 4]);
        ts = permute(srfc.temp2, [2 1 3]); % bring lat to front to interpolate
        ts = interp1(grid.dim2.lat, ts, lat); % interpolate to standard grid
        ts = permute(ts, [2 1 3]); % reorder to original
        ts_vert = repmat(ts, [1 1 1 size(taz_orig, 3)]); % dims (lon x lat x time x plev)
        ts_vert = permute(ts_vert, [1 2 4 3]); % dims (lon x lat x plev x time)
    end

    if strcmp(type, 'echam_ml')
        sm = ones(size(ta_orig));
        smz = ones(size(taz_orig));
    else
        sm = nan(size(ta_orig));
        sm(pa < ps_vert) = 1;
        smz = nan(size(taz_orig));
        smz(za > zs_vert) = 1;
    end

    ta_sm.lo = ta_orig.*sm; % filter ta with surface mask
    taz_sm.lo = taz_orig.*smz; % filter taz with surface mask
    tasi_sm.lo = tasi_orig; % surface is already masked in standard sigma coordinates

    if par.do_surf
        % add tsurf data and interpolate to higher resolution vertical grid
        [pa_plus ta_plus] = deal(nan([size(pa,1), size(pa,2) size(pa,3)+1 size(pa,4)])); % create empty grid with one extra vertical level
        pa_plus(:,:,1:end-1,:) = pa; % populate with standard pressure grid
        ta_plus(:,:,1:end-1,:) = ta_sm.lo; % populate with standard atmospheric temperature
        pa_plus(:,:,end,:) = ps_vert(:,:,1,:); % add surface pressure data into standard pressure grid
        ta_plus(:,:,end,:) = ts_vert(:,:,1,:); % add surface temperature data
        pa_plus = permute(pa_plus, [3 1 2 4]); % bring plev dimension to front
        ta_plus = permute(ta_plus, [3 1 2 4]); % bring plev dimension to front
        [pa_plus sort_index] = sort(pa_plus, 1, 'descend'); % sort added surface pressure such that pressure decreases monotonically
        tai_sm.lo = nan(length(par.pa), size(pa, 1), size(pa, 2), size(pa, 4));
        pb = CmdLineProgressBar("Sorting and interpolating temperature to new standard grid...");
        for ilon=1:size(pa_plus,2)
            pb.print(ilon, size(pa_plus,2));
            for ilat=1:size(pa_plus,3)
                for time=1:size(pa_plus,4)
                    ta_plus(:,ilon,ilat,time) = ta_plus(sort_index(:,ilon,ilat,time),ilon,ilat,time); % sort temperature (has to be in loop because sort_index works for vector calls only)
                    tai_sm.lo(:,ilon,ilat,time) = interp1(pa_plus(:,ilon,ilat,time), ta_plus(:,ilon,ilat,time), par.pa, 'linear'); % interpolate to higher resolution vertical grid
                end
            end
        end
        clear pa_plus ta_plus; % clear unneeded variables
        tai_sm.lo = permute(tai_sm.lo, [2 3 1 4]); % bring plev back to third dimension
    end

    % Land/ocean filter 2D variables
    tsurf_sm.lo = ts;
    psurf_sm.lo = ps;
    zsurf_sm.lo = zs;
    tsurf_sm.l = ts.*mask.ocean; %filter out ocean
    psurf_sm.l = ps.*mask.ocean; %filter out ocean
    zsurf_sm.l = zs.*mask.ocean; %filter out ocean
    tsurf_sm.o = ts.*mask.land; %filter out land
    psurf_sm.o = ps.*mask.land; %filter out land
    zsurf_sm.o = zs.*mask.land; %filter out land

    % Land/ocean filter 3D variables
    mask.land_vert = repmat(mask.land, [1 1 1 size(ta_orig, 3)]); % expand land mask to vertical dim
    mask.land_vert = permute(mask.land, [1 2 4 3]); % place vertical dim where it belongs
    mask.ocean_vert = repmat(mask.ocean, [1 1 1 size(ta_orig, 3)]); % expand ocean mask to vertical dim
    mask.ocean_vert = permute(mask.ocean, [1 2 4 3]); % place vertical dim where it belongs
    mask.land_verti = repmat(mask.land, [1 1 1 length(par.pa)]); % expand land mask to vertiical dim
    mask.land_verti = permute(mask.land, [1 2 4 3]); % place vertiical dim where it belongs
    mask.ocean_verti = repmat(mask.ocean, [1 1 1 length(par.pa)]); % expand ocean mask to vertiical dim
    mask.ocean_verti = permute(mask.ocean, [1 2 4 3]); % place vertiical dim where it belongs

    ta_sm.l = ta_sm.lo.*mask.ocean_vert; % filter ta with surface mask
    ta_sm.o = ta_sm.lo.*mask.land_vert; % filter ta with surface mask
    ta_sm.lo = permute(ta_sm.lo, [1 2 4 3]); % bring plev to last dimension
    ta_sm.l = permute(ta_sm.l, [1 2 4 3]); % bring plev to last dimension
    ta_sm.o = permute(ta_sm.o, [1 2 4 3]); % bring plev to last dimension
    taz_sm.l = taz_sm.lo.*mask.ocean_vert; % filter taz with surface mask
    taz_sm.o = taz_sm.lo.*mask.land_vert; % filter taz with surface mask
    taz_sm.lo = permute(taz_sm.lo, [1 2 4 3]); % bring plev to last dimension
    taz_sm.l = permute(taz_sm.l, [1 2 4 3]); % bring plev to last dimension
    taz_sm.o = permute(taz_sm.o, [1 2 4 3]); % bring plev to last dimension
    tasi_sm.l = tasi_sm.lo.*mask.ocean_vert; % filter tasi with surface mask
    tasi_sm.o = tasi_sm.lo.*mask.land_vert; % filter tasi with surface mask
    tasi_sm.lo = permute(tasi_sm.lo, [1 2 4 3]); % bring plev to last dimension
    tasi_sm.l = permute(tasi_sm.l, [1 2 4 3]); % bring plev to last dimension
    tasi_sm.o = permute(tasi_sm.o, [1 2 4 3]); % bring plev to last dimension
    if par.do_surf
        tai_sm.l = tai_sm.lo.*mask.ocean_verti; % filter tai with surface mask
        tai_sm.o = tai_sm.lo.*mask.land_verti; % filter tai with surface mask
        tai_sm.lo = permute(tai_sm.lo, [1 2 4 3]); % bring plev to last dimension
        tai_sm.l = permute(tai_sm.l, [1 2 4 3]); % bring plev to last dimension
        tai_sm.o = permute(tai_sm.o, [1 2 4 3]); % bring plev to last dimension
    end

    mask_t.land = nanmean(mask.land, 3);
    mask_t.ocean = nanmean(mask.ocean, 3);

    for l = {'lo', 'l', 'o'}; land = l{1}; % over land, over ocean, or both
        ta.(land)= squeeze(nanmean(ta_sm.(land), 1)); % zonal average
        taz.(land)= squeeze(nanmean(taz_sm.(land), 1)); % zonal average
        tasi.(land)= squeeze(nanmean(tasi_sm.(land), 1)); % zonal average
        if par.do_surf; tai.(land)= squeeze(nanmean(tai_sm.(land), 1)); end % zonal average
        tsurf.(land)= squeeze(nanmean(tsurf_sm.(land), 1)); % zonal average
        psurf.(land)= squeeze(nanmean(psurf_sm.(land), 1)); % zonal average
        zsurf.(land)= squeeze(nanmean(zsurf_sm.(land), 1)); % zonal average
    end

    for l = {'lo', 'l', 'o'}; land = l{1}; % over land, over ocean, or both
        % take time averages
        for t = {'ann', 'djf', 'jja', 'mam', 'son'}; time = t{1};
            if strcmp(time, 'ann')
                ta_t.(land).(time) = squeeze(nanmean(ta_sm.(land), 3));
                taz_t.(land).(time) = squeeze(nanmean(taz_sm.(land), 3));
                tasi_t.(land).(time) = squeeze(nanmean(tasi_sm.(land), 3));
                if par.do_surf; tai_t.(land).(time) = squeeze(nanmean(tai_sm.(land), 3)); end
                tsurf_t.(land).(time) = squeeze(nanmean(tsurf_sm.(land), 3));
                psurf_t.(land).(time) = squeeze(nanmean(psurf_sm.(land), 3));
                zsurf_t.(land).(time) = squeeze(nanmean(zsurf_sm.(land), 3));
            elseif strcmp(time, 'djf')
                ta_shift.(land) = circshift(ta_sm.(land), 1, 3);
                taz_shift.(land) = circshift(taz_sm.(land), 1, 3);
                tasi_shift.(land) = circshift(tasi_sm.(land), 1, 3);
                if par.do_surf; tai_shift.(land) = circshift(tai_sm.(land), 1, 3); end
                tsurf_shift.(land) = circshift(tsurf_sm.(land), 1, 3);
                psurf_shift.(land) = circshift(psurf_sm.(land), 1, 3);
                zsurf_shift.(land) = circshift(zsurf_sm.(land), 1, 3);
                ta_t.(land).(time) = squeeze(nanmean(ta_shift.(land)(:,:,1:3,:), 3));
                taz_t.(land).(time) = squeeze(nanmean(taz_shift.(land)(:,:,1:3,:), 3));
                tasi_t.(land).(time) = squeeze(nanmean(tasi_shift.(land)(:,:,1:3,:), 3));
                if par.do_surf; tai_t.(land).(time) = squeeze(nanmean(tai_shift.(land)(:,:,1:3,:), 3)); end
                tsurf_t.(land).(time) = squeeze(nanmean(tsurf_shift.(land)(:,:,1:3,:), 3));
                psurf_t.(land).(time) = squeeze(nanmean(psurf_shift.(land)(:,:,1:3,:), 3));
                zsurf_t.(land).(time) = squeeze(nanmean(zsurf_shift.(land)(:,:,1:3,:), 3));
            elseif strcmp(time, 'jja')
                ta_t.(land).(time) = squeeze(nanmean(ta_sm.(land)(:,:,6:8,:), 3));
                taz_t.(land).(time) = squeeze(nanmean(taz_sm.(land)(:,:,6:8,:), 3));
                tasi_t.(land).(time) = squeeze(nanmean(tasi_sm.(land)(:,:,6:8,:), 3));
                if par.do_surf; tai_t.(land).(time) = squeeze(nanmean(tai_sm.(land)(:,:,6:8,:), 3)); end
                tsurf_t.(land).(time) = squeeze(nanmean(tsurf_sm.(land)(:,:,6:8,:), 3));
                psurf_t.(land).(time) = squeeze(nanmean(psurf_sm.(land)(:,:,6:8,:), 3));
                zsurf_t.(land).(time) = squeeze(nanmean(zsurf_sm.(land)(:,:,6:8,:), 3));
            elseif strcmp(time, 'mam')
                ta_t.(land).(time) = squeeze(nanmean(ta_sm.(land)(:,:,3:5,:), 3));
                taz_t.(land).(time) = squeeze(nanmean(taz_sm.(land)(:,:,3:5,:), 3));
                tasi_t.(land).(time) = squeeze(nanmean(tasi_sm.(land)(:,:,3:5,:), 3));
                if par.do_surf; tai_t.(land).(time) = squeeze(nanmean(tai_sm.(land)(:,:,3:5,:), 3)); end
                tsurf_t.(land).(time) = squeeze(nanmean(tsurf_sm.(land)(:,:,3:5,:), 3));
                psurf_t.(land).(time) = squeeze(nanmean(psurf_sm.(land)(:,:,3:5,:), 3));
                zsurf_t.(land).(time) = squeeze(nanmean(zsurf_sm.(land)(:,:,3:5,:), 3));
            elseif strcmp(time, 'son')
                ta_t.(land).(time) = squeeze(nanmean(ta_sm.(land)(:,:,9:11,:), 3));
                taz_t.(land).(time) = squeeze(nanmean(taz_sm.(land)(:,:,9:11,:), 3));
                tasi_t.(land).(time) = squeeze(nanmean(tasi_sm.(land)(:,:,9:11,:), 3));
                if par.do_surf; tai_t.(land).(time) = squeeze(nanmean(tai_sm.(land)(:,:,9:11,:), 3)); end
                tsurf_t.(land).(time) = squeeze(nanmean(tsurf_sm.(land)(:,:,9:11,:), 3));
                psurf_t.(land).(time) = squeeze(nanmean(psurf_sm.(land)(:,:,9:11,:), 3));
                zsurf_t.(land).(time) = squeeze(nanmean(zsurf_sm.(land)(:,:,9:11,:), 3));
            end
        end
    end

    % save filtered data
    printname = [foldername 'ta_mon_lat'];
    if ~exist(foldername, 'dir')
        mkdir(foldername)
    end
    if par.do_surf; save(printname, 'ta', 'taz', 'tasi', 'tai', 'lat', 'tsurf', 'psurf', 'zsurf');
    else save(printname, 'ta', 'taz', 'tasi', 'lat', 'tsurf', 'psurf', 'zsurf', '-v7.3'); end

    printname = [foldername 'ta_lon_lat'];
    if ~exist(foldername, 'dir')
        mkdir(foldername)
    end
    if par.do_surf; save(printname, 'ta_t', 'taz_t', 'tasi_t', 'tai_t', 'lat', 'tsurf_t', 'psurf_t', 'zsurf_t');
    else save(printname, 'ta_t', 'taz_t', 'tasi_t', 'lat', 'tsurf_t', 'psurf_t', 'zsurf_t', '-v7.3'); end
end % compute mon x lat temperature field
function proc_ma_mon_lat(type, par)
% calculate moist adiabats
    if strcmp(type, 'era5') | strcmp(type, 'erai')
        foldername = sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/', type, par.lat_interp);
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s', type);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s', type);
    elseif strcmp(type, 'gcm')
        foldername = sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s/%s/', type, par.model, par.gcm.clim, par.lat_interp);
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s/%s', type, par.model, par.gcm.clim);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s', type, par.model, par.gcm.clim);
    elseif strcmp(type, 'echam')
        foldername = sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s/', type, par.echam.clim, par.lat_interp);
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s', type, par.echam.clim);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s', type, par.echam.clim);
    end
    load(sprintf('%s/grid.mat', prefix)); % read grid data
    load(sprintf('%s/srfc.mat', prefix)); % read surface variable data
    load(sprintf('%s/%s/masks.mat', prefix_proc, par.lat_interp)); % load land and ocean masks

    lat = par.lat_std;
    for fn = fieldnames(srfc)'
        srfc.(fn{1}) = permute(srfc.(fn{1}), [2 1 3]); % bring lat to front
        srfc.(fn{1}) = interp1(grid.dim2.lat, srfc.(fn{1}), lat); % interpolate to standard grid
        srfc.(fn{1}) = permute(srfc.(fn{1}), [2 1 3]); % reorder to original dims
        for l = {'lo', 'l', 'o'}; land = l{1};
            if strcmp(land, 'lo'); srfc_n.(fn{1}).(land) = srfc.(fn{1});
            elseif strcmp(land, 'l'); srfc_n.(fn{1}).(land) = srfc.(fn{1}).*mask.ocean; % filter out ocean
            elseif strcmp(land, 'o'); srfc_n.(fn{1}).(land) = srfc.(fn{1}).*mask.land; % filter out land
            end
        end
    end

    for l = {'lo', 'l', 'o'}; land = l{1};
        for fn_var = fieldnames(srfc)'
            ma.(land).(fn_var{1}) = squeeze(nanmean(srfc_n.(fn_var{1}).(land), 1)); % zonal average
            ma_t.(land).(fn_var{1}) = squeeze(nanmean(srfc_n.(fn_var{1}).(land), 3)); % time average

        end % end srfc variables loop

        if strcmp(type, 'era5') | strcmp(type, 'erai')
            pb = CmdLineProgressBar("Calculating moist adiabats...");
            for ilat = 1:length(lat);
                pb.print(ilat, length(lat)); % output progress of moist adiabat calculation
                for imon = 1:12;
                    for fn_var = fieldnames(srfc)'
                        ima.(fn_var{1}) = ma.(land).(fn_var{1})(ilat, imon);
                    end
                    ma.(land).ta(ilat, imon, :) = calc_ma_dew(ima, grid.dim3.plev, par, type); % compute moist adiabat with RH
                    maz.(land).ta(ilat, imon, :) = calc_maz_dew(ima, grid.dim3.z, par, type); % compute moist adiabat with RH
                end
                for ilon = 1:length(grid.dim3.lon);
                    for fn_var = fieldnames(srfc)'
                        ima_t.(fn_var{1}) = ma_t.(land).(fn_var{1})(ilon, ilat);
                    end
                    ma_t.(land).ta(ilat, imon, :) = calc_ma_dew(ima_t, grid.dim3.plev, par, type); % compute moist adiabat with RH
                    maz_t.(land).ta(ilat, imon, :) = calc_maz_dew(ima_t, grid.dim3.z, par, type); % compute moist adiabat with RH
                end
            end
        elseif strcmp(type, 'gcm')
            pb = CmdLineProgressBar("Calculating moist adiabats...");
            for ilat = 1:length(lat);
                pb.print(ilat, length(lat)); % output progress of moist adiabat calculation
                for imon = 1:12;
                    for fn_var = fieldnames(srfc)'
                        ima.(fn_var{1}) = ma.(land).(fn_var{1})(ilat, imon);
                    end
                    ma.(land).ta(ilat, imon, :) = calc_ma_hurs(ima, grid.dim3.plev, par); % compute moist adiabat with RH
                    maz.(land).ta(ilat, imon, :) = calc_maz_hurs(ima, grid.dim3.z, par); % compute moist adiabat with RH
                end
                for ilon = 1:length(grid.dim3.lon);
                    for fn_var = fieldnames(srfc)'
                        ima_t.(fn_var{1}) = ma_t.(land).(fn_var{1})(ilon, ilat);
                    end
                    ma_t.(land).ta(ilat, imon, :) = calc_ma_hurs(ima_t, grid.dim3.plev, par); % compute moist adiabat with RH
                    maz_t.(land).ta(ilat, imon, :) = calc_maz_hurs(ima_t, grid.dim3.z, par); % compute moist adiabat with RH
                end
            end
        elseif strcmp(type, 'echam')
            pb = CmdLineProgressBar("Calculating moist adiabats...");
            for ilat = 1:length(lat);
                pb.print(ilat, length(lat)); % output progress of moist adiabat calculation
                for imon = 1:12;
                    for fn_var = fieldnames(srfc)'
                        ima.(fn_var{1}) = ma.(land).(fn_var{1})(ilat, imon);
                    end
                    ma.(land).ta(ilat, imon, :) = calc_ma_dew(ima, grid.dim3.plev, par, type); % compute moist adiabat with RH
                    maz.(land).ta(ilat, imon, :) = calc_maz_dew(ima, grid.dim3.z, par, type); % compute moist adiabat with RH
                end
                for ilon = 1:length(grid.dim3.lon);
                    for fn_var = fieldnames(srfc)'
                        ima_t.(fn_var{1}) = ma_t.(land).(fn_var{1})(ilon, ilat);
                    end
                    ma_t.(land).ta(ilat, imon, :) = calc_ma_dew(ima_t, grid.dim3.plev, par, type); % compute moist adiabat with RH
                    maz_t.(land).ta(ilat, imon, :) = calc_maz_dew(ima_t, grid.dim3.z, par, type); % compute moist adiabat with RH
                end
            end
        end
    end % end land option loop

    % save data into mat file
    printname = [foldername 'ma_mon_lat.mat'];
    if ~exist(foldername, 'dir')
        mkdir(foldername)
    end
    save(printname, 'ma', 'maz');

    printname = [foldername 'ma_lon_lat.mat'];
    if ~exist(foldername, 'dir')
        mkdir(foldername)
    end
    save(printname, 'ma_t', 'maz_t');

end % compute mon x lat moist adiabat field
function proc_inv_str_mon_lat(type, par)
    if strcmp(type, 'era5') | strcmp(type, 'erai')
        foldername = sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/', type, par.lat_interp);
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s', type);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s', type);
        file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/%s/temp/%s_temp_%s.ymonmean.nc', type, type, par.(type).yr_span));
        fullpath=sprintf('%s/%s', file.folder, file.name);
        ta_orig = ncread(fullpath, 't');
    elseif strcmp(type, 'gcm')
        foldername = sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s/%s/', type, par.model, par.gcm.clim, par.lat_interp);
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s/%s', type, par.model, par.gcm.clim);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s', type, par.model, par.gcm.clim);
        var = 'ta';
        file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/gcm/%s/%s_Amon_%s_%s_r1i1p1_*.ymonmean.nc', par.model, var, par.model, par.gcm.clim));
        fullpath=sprintf('%s/%s', file.folder, file.name);
        ta_orig = ncread(fullpath, var);
    elseif strcmp(type, 'echam')
        foldername = sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s/', type, par.echam.clim, par.lat_interp);
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s', type, par.echam.clim);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s', type, par.echam.clim);
        var = 't';
        file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/echam/ATM_rjg_%s_*.ymonmean.nc', par.echam.clim));
        fullpath=sprintf('%s/%s', file.folder, file.name);
        ta_orig = double(ncread(fullpath, var));
    end
    load(sprintf('%s/grid.mat', prefix)); % read grid data
    load(sprintf('%s/tempsi.mat', prefix));tasi_orig = tempsi; clear tempsi; % read temp in si coordinates
    load(sprintf('%s/%s/masks.mat', prefix_proc, par.lat_interp)); % load land and ocean masks

    lat = par.lat_std;

    % interpolate ta to standard lat grid
    tasi_orig = permute(tasi_orig, [2 1 3 4]);
    tasi_orig = interp1(grid.dim3.lat, tasi_orig, lat);
    tasi_orig = permute(tasi_orig, [2 1 3 4]);

    tasi_orig = permute(tasi_orig, [3 1 2 4]); % bring sigma forward

    for isig = 1:length(par.si_eval); sig = par.si_eval(isig);
        inv_orig(:,:,:,isig) = squeeze(interp1(grid.dim3.si, tasi_orig, sig) - tasi_orig(1,:,:,:)); % evaluate difference between surface and si_eval
    end

    % Land/ocean filter 2D variables
    mask.land_vert = repmat(mask.land, [1 1 1 length(par.si_eval)]); % expand land mask to si_eval dim
    mask.ocean_vert = repmat(mask.ocean, [1 1 1 length(par.si_eval)]); % expand ocean mask to si_eval dim

    inv0.lo = inv_orig;
    inv0.l = inv0.lo.*mask.ocean_vert; % filter inv0 with surface mask
    inv0.o = inv0.lo.*mask.land_vert; % filter inv0 with surface mask

    for l = {'lo', 'l', 'o'}; land = l{1}; % over land, over ocean, or both
        inv.(land)= squeeze(nanmean(inv0.(land), 1)); % zonal average
        % take time averages
        for t = {'ann', 'djf', 'jja', 'mam', 'son'}; time = t{1};
            if strcmp(time, 'ann')
                inv_t.(land).(time) = squeeze(nanmean(inv0.(land), 3));
            elseif strcmp(time, 'djf')
                inv_shift.(land) = circshift(inv0.(land), 1, 3);
                inv_t.(land).(time) = squeeze(nanmean(inv_shift.(land)(:,:,1:3,:), 3));
            elseif strcmp(time, 'jja')
                inv_t.(land).(time) = squeeze(nanmean(inv0.(land)(:,:,6:8,:), 3));
            elseif strcmp(time, 'mam')
                inv_t.(land).(time) = squeeze(nanmean(inv0.(land)(:,:,3:5,:), 3));
            elseif strcmp(time, 'son')
                inv_t.(land).(time) = squeeze(nanmean(inv0.(land)(:,:,9:11,:), 3));
            end
        end
    end

    % save filtered data
    printname = [foldername 'inv_mon_lat'];
    if ~exist(foldername, 'dir')
        mkdir(foldername)
    end
    save(printname, 'inv', 'lat');

    printname = [foldername 'inv_lon_lat'];
    if ~exist(foldername, 'dir')
        mkdir(foldername)
    end
    save(printname, 'inv_t', 'lat');
end % compute mon x lat inversion strength field with land/ocean masking
function proc_ga_diff_mon_lat(type, par)
    if strcmp(type, 'era5') | strcmp(type, 'erai')
        foldername = sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/', type, par.lat_interp);
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s', type);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s', type);
        file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/%s/temp/%s_temp_%s.ymonmean.nc', type, type, par.(type).yr_span));
        fullpath=sprintf('%s/%s', file.folder, file.name);
        ta_orig = ncread(fullpath, 't');
    elseif strcmp(type, 'gcm')
        foldername = sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s/%s/', type, par.model, par.gcm.clim, par.lat_interp);
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s/%s', type, par.model, par.gcm.clim);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s', type, par.model, par.gcm.clim);
        var = 'ta';
        file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/gcm/%s/%s_Amon_%s_%s_r1i1p1_*.ymonmean.nc', par.model, var, par.model, par.gcm.clim));
        fullpath=sprintf('%s/%s', file.folder, file.name);
        ta_orig = ncread(fullpath, var);
    elseif strcmp(type, 'echam')
        foldername = sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s/', type, par.echam.clim, par.lat_interp);
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s', type, par.echam.clim);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s', type, par.echam.clim);
        var = 't';
        file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/echam/ATM_rjg_%s_*.ymonmean.nc', par.echam.clim));
        fullpath=sprintf('%s/%s', file.folder, file.name);
        ta_orig = double(ncread(fullpath, var));
    end
    load(sprintf('%s/grid.mat', prefix)); % read grid data
    load(sprintf('%s/dtdz.mat', prefix)); % read temp in si coordinates
    load(sprintf('%s/dtmdz.mat', prefix)); % read temp in si coordinates
    load(sprintf('%s/%s/masks.mat', prefix_proc, par.lat_interp)); % load land and ocean masks

    lat = par.lat_std;

    ga_diff_orig = (dtmdz - dtdz)./dtmdz * 1e2; % percentage difference of moist adiabatic vs GCM lapse rates
    clear dtmdz dtdz;
    ga_diff_orig = permute(ga_diff_orig,[2 1 3 4]); % bring lat to 1st
    ga_diff_orig = interp1(grid.dim3.lat, ga_diff_orig, lat); % interpolate to standard lat
    ga_diff_orig = permute(ga_diff_orig, [3 2 1 4]); % bring height front
    ga_diff_orig = interp1(par.pa, ga_diff_orig, 1e2*linspace(1000,200,100)); % prepare to average between 1000-200 hPa
    ga_diff_orig = squeeze(nanmean(ga_diff_orig,1)); % take vertical average

    ga_diff0.lo = ga_diff_orig;
    ga_diff0.l = ga_diff0.lo.*mask.ocean; % filter ga_diff0 with surface mask
    ga_diff0.o = ga_diff0.lo.*mask.land; % filter ga_diff0 with surface mask

    for l = {'lo', 'l', 'o'}; land = l{1}; % over land, over ocean, or both
        ga_diff.(land)= squeeze(nanmean(ga_diff0.(land), 1)); % zonal average
        % take time averages
        for t = {'ann', 'djf', 'jja', 'mam', 'son'}; time = t{1};
            if strcmp(time, 'ann')
                ga_diff_t.(land).(time) = squeeze(nanmean(ga_diff0.(land), 3));
            elseif strcmp(time, 'djf')
                ga_diff_shift.(land) = circshift(ga_diff0.(land), 1, 3);
                ga_diff_t.(land).(time) = squeeze(nanmean(ga_diff_shift.(land)(:,:,1:3,:), 3));
            elseif strcmp(time, 'jja')
                ga_diff_t.(land).(time) = squeeze(nanmean(ga_diff0.(land)(:,:,6:8,:), 3));
            elseif strcmp(time, 'mam')
                ga_diff_t.(land).(time) = squeeze(nanmean(ga_diff0.(land)(:,:,3:5,:), 3));
            elseif strcmp(time, 'son')
                ga_diff_t.(land).(time) = squeeze(nanmean(ga_diff0.(land)(:,:,9:11,:), 3));
            end
        end
    end

    % save filtered data
    printname = [foldername 'ga_diff_mon_lat'];
    if ~exist(foldername, 'dir')
        mkdir(foldername)
    end
    save(printname, 'ga_diff', 'lat');

    printname = [foldername 'ga_diff_lon_lat'];
    if ~exist(foldername, 'dir')
        mkdir(foldername)
    end
    save(printname, 'ga_diff_t', 'lat');
end % compute mon x lat gamma percentage difference field with land/ocean masking
function proc_ga_diff_si_mon_lat(type, par)
    if strcmp(type, 'era5') | strcmp(type, 'erai')
        foldername = sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/', type, par.lat_interp);
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s', type);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s', type);
        file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/%s/temp/%s_temp_%s.ymonmean.nc', type, type, par.(type).yr_span));
        fullpath=sprintf('%s/%s', file.folder, file.name);
        ta_orig = ncread(fullpath, 't');
    elseif strcmp(type, 'gcm')
        foldername = sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s/%s/', type, par.model, par.gcm.clim, par.lat_interp);
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s/%s', type, par.model, par.gcm.clim);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s', type, par.model, par.gcm.clim);
        var = 'ta';
        file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/gcm/%s/%s_Amon_%s_%s_r1i1p1_*.ymonmean.nc', par.model, var, par.model, par.gcm.clim));
        fullpath=sprintf('%s/%s', file.folder, file.name);
        ta_orig = ncread(fullpath, var);
    elseif strcmp(type, 'echam')
        foldername = sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s/', type, par.echam.clim, par.lat_interp);
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s', type, par.echam.clim);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s', type, par.echam.clim);
        var = 't';
        file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/echam/ATM_rjg_%s_*.ymonmean.nc', par.echam.clim));
        fullpath=sprintf('%s/%s', file.folder, file.name);
        ta_orig = double(ncread(fullpath, var));
    end
    load(sprintf('%s/grid.mat', prefix)); % read grid data
    load(sprintf('%s/dtdzzsi.mat', prefix)); % read temp in si coordinates
    load(sprintf('%s/dtmdzzsi.mat', prefix)); % read temp in si coordinates
    load(sprintf('%s/%s/masks.mat', prefix_proc, par.lat_interp)); % load land and ocean masks

    lat = par.lat_std;

    ga_diff_orig = (dtmdzzsi - dtdzzsi)./dtmdzzsi * 1e2; % percentage difference of moist adiabatic vs GCM lapse rates
    clear dtmdz dtdz;
    ga_diff_orig = permute(ga_diff_orig,[2 1 3 4]); % bring lat to 1st
    ga_diff_orig = interp1(grid.dim3.lat, ga_diff_orig, lat); % interpolate to standard lat
    ga_diff_orig = permute(ga_diff_orig, [3 2 1 4]); % bring height front
    ga_diff_orig = interp1(grid.dim3.si, ga_diff_orig, linspace(par.si_bl,par.si_up,100)); % prepare to average between 1000-200 hPa
    ga_diff_orig = squeeze(nanmean(ga_diff_orig,1)); % take vertical average

    ga_diff0.lo = ga_diff_orig;
    ga_diff0.l = ga_diff0.lo.*mask.ocean; % filter ga_diff0 with surface mask
    ga_diff0.o = ga_diff0.lo.*mask.land; % filter ga_diff0 with surface mask

    for l = {'lo', 'l', 'o'}; land = l{1}; % over land, over ocean, or both
        ga_diff.(land)= squeeze(nanmean(ga_diff0.(land), 1)); % zonal average
        % take time averages
        for t = {'ann', 'djf', 'jja', 'mam', 'son'}; time = t{1};
            if strcmp(time, 'ann')
                ga_diff_t.(land).(time) = squeeze(nanmean(ga_diff0.(land), 3));
            elseif strcmp(time, 'djf')
                ga_diff_shift.(land) = circshift(ga_diff0.(land), 1, 3);
                ga_diff_t.(land).(time) = squeeze(nanmean(ga_diff_shift.(land)(:,:,1:3,:), 3));
            elseif strcmp(time, 'jja')
                ga_diff_t.(land).(time) = squeeze(nanmean(ga_diff0.(land)(:,:,6:8,:), 3));
            elseif strcmp(time, 'mam')
                ga_diff_t.(land).(time) = squeeze(nanmean(ga_diff0.(land)(:,:,3:5,:), 3));
            elseif strcmp(time, 'son')
                ga_diff_t.(land).(time) = squeeze(nanmean(ga_diff0.(land)(:,:,9:11,:), 3));
            end
        end
    end

    % save filtered data
    printname = [foldername 'ga_diff_si_mon_lat'];
    if ~exist(foldername, 'dir')
        mkdir(foldername)
    end
    save(printname, 'ga_diff', 'lat');

    printname = [foldername 'ga_diff_si_lon_lat'];
    if ~exist(foldername, 'dir')
        mkdir(foldername)
    end
    save(printname, 'ga_diff_t', 'lat');
end % compute mon x lat gamma percentage difference field with land/ocean masking
function proc_ga_bl_diff_si_mon_lat(type, par)
    if strcmp(type, 'era5') | strcmp(type, 'erai')
        foldername = sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/', type, par.lat_interp);
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s', type);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s', type);
        file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/%s/temp/%s_temp_%s.ymonmean.nc', type, type, par.(type).yr_span));
        fullpath=sprintf('%s/%s', file.folder, file.name);
        ta_orig = ncread(fullpath, 't');
    elseif strcmp(type, 'gcm')
        foldername = sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s/%s/', type, par.model, par.gcm.clim, par.lat_interp);
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s/%s', type, par.model, par.gcm.clim);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s', type, par.model, par.gcm.clim);
        var = 'ta';
        file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/gcm/%s/%s_Amon_%s_%s_r1i1p1_*.ymonmean.nc', par.model, var, par.model, par.gcm.clim));
        fullpath=sprintf('%s/%s', file.folder, file.name);
        ta_orig = ncread(fullpath, var);
    elseif strcmp(type, 'echam')
        foldername = sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s/', type, par.echam.clim, par.lat_interp);
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s', type, par.echam.clim);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s', type, par.echam.clim);
        var = 't';
        file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/echam/ATM_rjg_%s_*.ymonmean.nc', par.echam.clim));
        fullpath=sprintf('%s/%s', file.folder, file.name);
        ta_orig = double(ncread(fullpath, var));
    end
    load(sprintf('%s/grid.mat', prefix)); % read grid data
    load(sprintf('%s/dtdzzsi.mat', prefix)); % read temp in si coordinates
    load(sprintf('%s/dtmdzzsi.mat', prefix)); % read temp in si coordinates
    load(sprintf('%s/%s/masks.mat', prefix_proc, par.lat_interp)); % load land and ocean masks

    lat = par.lat_std;

    ga_bl_diff_orig = (dtmdzzsi - dtdzzsi)./dtmdzzsi * 1e2; % percentage difference of moist adiabatic vs GCM lapse rates
    clear dtmdz dtdz;
    ga_bl_diff_orig = permute(ga_bl_diff_orig,[2 1 3 4]); % bring lat to 1st
    ga_bl_diff_orig = interp1(grid.dim3.lat, ga_bl_diff_orig, lat); % interpolate to standard lat
    ga_bl_diff_orig = permute(ga_bl_diff_orig, [3 2 1 4]); % bring height front
    ga_bl_diff_orig = interp1(grid.dim3.si, ga_bl_diff_orig, linspace(1,par.si_bl,100)); % prepare to average between 1000-200 hPa
    ga_bl_diff_orig = squeeze(nanmean(ga_bl_diff_orig,1)); % take vertical average

    ga_bl_diff0.lo = ga_bl_diff_orig;
    ga_bl_diff0.l = ga_bl_diff0.lo.*mask.ocean; % filter ga_bl_diff0 with surface mask
    ga_bl_diff0.o = ga_bl_diff0.lo.*mask.land; % filter ga_bl_diff0 with surface mask

    for l = {'lo', 'l', 'o'}; land = l{1}; % over land, over ocean, or both
        ga_bl_diff.(land)= squeeze(nanmean(ga_bl_diff0.(land), 1)); % zonal average
        % take time averages
        for t = {'ann', 'djf', 'jja', 'mam', 'son'}; time = t{1};
            if strcmp(time, 'ann')
                ga_bl_diff_t.(land).(time) = squeeze(nanmean(ga_bl_diff0.(land), 3));
            elseif strcmp(time, 'djf')
                ga_bl_diff_shift.(land) = circshift(ga_bl_diff0.(land), 1, 3);
                ga_bl_diff_t.(land).(time) = squeeze(nanmean(ga_bl_diff_shift.(land)(:,:,1:3,:), 3));
            elseif strcmp(time, 'jja')
                ga_bl_diff_t.(land).(time) = squeeze(nanmean(ga_bl_diff0.(land)(:,:,6:8,:), 3));
            elseif strcmp(time, 'mam')
                ga_bl_diff_t.(land).(time) = squeeze(nanmean(ga_bl_diff0.(land)(:,:,3:5,:), 3));
            elseif strcmp(time, 'son')
                ga_bl_diff_t.(land).(time) = squeeze(nanmean(ga_bl_diff0.(land)(:,:,9:11,:), 3));
            end
        end
    end

    % save filtered data
    printname = [foldername 'ga_bl_diff_si_mon_lat'];
    if ~exist(foldername, 'dir')
        mkdir(foldername)
    end
    save(printname, 'ga_bl_diff', 'lat');

    printname = [foldername 'ga_bl_diff_si_lon_lat'];
    if ~exist(foldername, 'dir')
        mkdir(foldername)
    end
    save(printname, 'ga_bl_diff_t', 'lat');
end % compute mon x lat gamma percentage difference field with land/ocean masking
function proc_ga_malr_diff_si_mon_lat(type, par)
    if strcmp(type, 'era5') | strcmp(type, 'erai')
        foldername = sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/', type, par.lat_interp);
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s', type);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s', type);
        file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/%s/temp/%s_temp_%s.ymonmean.nc', type, type, par.(type).yr_span));
        fullpath=sprintf('%s/%s', file.folder, file.name);
        ta_orig = ncread(fullpath, 't');
    elseif strcmp(type, 'gcm')
        foldername = sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s/%s/', type, par.model, par.gcm.clim, par.lat_interp);
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s/%s', type, par.model, par.gcm.clim);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s', type, par.model, par.gcm.clim);
        var = 'ta';
        file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/gcm/%s/%s_Amon_%s_%s_r1i1p1_*.ymonmean.nc', par.model, var, par.model, par.gcm.clim));
        fullpath=sprintf('%s/%s', file.folder, file.name);
        ta_orig = ncread(fullpath, var);
    elseif strcmp(type, 'echam')
        foldername = sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s/', type, par.echam.clim, par.lat_interp);
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s', type, par.echam.clim);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s', type, par.echam.clim);
        var = 't';
        file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/echam/ATM_rjg_%s_*.ymonmean.nc', par.echam.clim));
        fullpath=sprintf('%s/%s', file.folder, file.name);
        ta_orig = double(ncread(fullpath, var));
    end
    load(sprintf('%s/grid.mat', prefix)); % read grid data
    % load(sprintf('%s/dtdzzsi.mat', prefix)); % read temp in si coordinates
    load(sprintf('%s/dtdzsi.mat', prefix)); dtdzzsi = dtdzsi; clear dtdzsi; % read temp in si coordinates
    % load(sprintf('%s/malrzsi.mat', prefix)); % read temp in si coordinates
    load(sprintf('%s/malrsi.mat', prefix)); dtmdzzsi = dtmdzsi; clear dtmdzsi; % read temp in si coordinates
    load(sprintf('%s/%s/masks.mat', prefix_proc, par.lat_interp)); % load land and ocean masks

    lat = par.lat_std;

    ga_malr_diff_orig = (dtmdzzsi - dtdzzsi)./dtmdzzsi * 1e2; % percentage difference of moist adiabatic vs GCM lapse rates
    clear dtmdz dtdz;
    ga_malr_diff_orig = permute(ga_malr_diff_orig,[2 1 3 4]); % bring lat to 1st
    ga_malr_diff_orig = interp1(grid.dim3.lat, ga_malr_diff_orig, lat); % interpolate to standard lat
    ga_malr_diff_orig = permute(ga_malr_diff_orig, [3 2 1 4]); % bring height front
    ga_malr_diff_orig = interp1(grid.dim3.si, ga_malr_diff_orig, linspace(par.si_bl,par.si_up,100)); % prepare to average between 1000-200 hPa
    ga_malr_diff_orig = squeeze(nanmean(ga_malr_diff_orig,1)); % take vertical average

    ga_malr_diff0.lo = ga_malr_diff_orig;
    ga_malr_diff0.l = ga_malr_diff0.lo.*mask.ocean; % filter ga_malr_diff0 with surface mask
    ga_malr_diff0.o = ga_malr_diff0.lo.*mask.land; % filter ga_malr_diff0 with surface mask

    for l = {'lo', 'l', 'o'}; land = l{1}; % over land, over ocean, or both
        ga_malr_diff.(land)= squeeze(nanmean(ga_malr_diff0.(land), 1)); % zonal average
        % take time averages
        for t = {'ann', 'djf', 'jja', 'mam', 'son'}; time = t{1};
            if strcmp(time, 'ann')
                ga_malr_diff_t.(land).(time) = squeeze(nanmean(ga_malr_diff0.(land), 3));
            elseif strcmp(time, 'djf')
                ga_malr_diff_shift.(land) = circshift(ga_malr_diff0.(land), 1, 3);
                ga_malr_diff_t.(land).(time) = squeeze(nanmean(ga_malr_diff_shift.(land)(:,:,1:3,:), 3));
            elseif strcmp(time, 'jja')
                ga_malr_diff_t.(land).(time) = squeeze(nanmean(ga_malr_diff0.(land)(:,:,6:8,:), 3));
            elseif strcmp(time, 'mam')
                ga_malr_diff_t.(land).(time) = squeeze(nanmean(ga_malr_diff0.(land)(:,:,3:5,:), 3));
            elseif strcmp(time, 'son')
                ga_malr_diff_t.(land).(time) = squeeze(nanmean(ga_malr_diff0.(land)(:,:,9:11,:), 3));
            end
        end
    end

    % save filtered data
    printname = [foldername 'ga_malr_diff_si_mon_lat'];
    if ~exist(foldername, 'dir')
        mkdir(foldername)
    end
    save(printname, 'ga_malr_diff', 'lat');

    printname = [foldername 'ga_malr_diff_si_lon_lat'];
    if ~exist(foldername, 'dir')
        mkdir(foldername)
    end
    save(printname, 'ga_malr_diff_t', 'lat');
end % compute mon x lat gamma percentage difference field with land/ocean masking
function proc_ga_dalr_diff_si_mon_lat(type, par)
    if strcmp(type, 'era5') | strcmp(type, 'erai')
        foldername = sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/', type, par.lat_interp);
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s', type);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s', type);
        file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/%s/temp/%s_temp_%s.ymonmean.nc', type, type, par.(type).yr_span));
        fullpath=sprintf('%s/%s', file.folder, file.name);
        ta_orig = ncread(fullpath, 't');
    elseif strcmp(type, 'gcm')
        foldername = sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s/%s/', type, par.model, par.gcm.clim, par.lat_interp);
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s/%s', type, par.model, par.gcm.clim);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s', type, par.model, par.gcm.clim);
        var = 'ta';
        file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/gcm/%s/%s_Amon_%s_%s_r1i1p1_*.ymonmean.nc', par.model, var, par.model, par.gcm.clim));
        fullpath=sprintf('%s/%s', file.folder, file.name);
        ta_orig = ncread(fullpath, var);
    elseif strcmp(type, 'echam')
        foldername = sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s/', type, par.echam.clim, par.lat_interp);
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s', type, par.echam.clim);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s', type, par.echam.clim);
        var = 't';
        file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/echam/ATM_rjg_%s_*.ymonmean.nc', par.echam.clim));
        fullpath=sprintf('%s/%s', file.folder, file.name);
        ta_orig = double(ncread(fullpath, var));
    end
    load(sprintf('%s/grid.mat', prefix)); % read grid data
    load(sprintf('%s/dtdzzsi.mat', prefix)); % read temp in si coordinates
    load(sprintf('%s/%s/masks.mat', prefix_proc, par.lat_interp)); % load land and ocean masks

    lat = par.lat_std;

    dtmdzzsi = 1e3*par.g/par.cpd*ones(size(dtdzzsi));

    ga_dalr_diff_orig = (dtmdzzsi - dtdzzsi)./dtmdzzsi * 1e2; % percentage difference of moist adiabatic vs GCM lapse rates
    clear dtmdz dtdz;
    ga_dalr_diff_orig = permute(ga_dalr_diff_orig,[2 1 3 4]); % bring lat to 1st
    ga_dalr_diff_orig = interp1(grid.dim3.lat, ga_dalr_diff_orig, lat); % interpolate to standard lat
    ga_dalr_diff_orig = permute(ga_dalr_diff_orig, [3 2 1 4]); % bring height front
    ga_dalr_diff_orig = interp1(grid.dim3.si, ga_dalr_diff_orig, linspace(par.si_bl,par.si_up,100)); % prepare to average between 1000-200 hPa
    ga_dalr_diff_orig = squeeze(nanmean(ga_dalr_diff_orig,1)); % take vertical average

    ga_dalr_diff0.lo = ga_dalr_diff_orig;
    ga_dalr_diff0.l = ga_dalr_diff0.lo.*mask.ocean; % filter ga_dalr_diff0 with surface mask
    ga_dalr_diff0.o = ga_dalr_diff0.lo.*mask.land; % filter ga_dalr_diff0 with surface mask

    for l = {'lo', 'l', 'o'}; land = l{1}; % over land, over ocean, or both
        ga_dalr_diff.(land)= squeeze(nanmean(ga_dalr_diff0.(land), 1)); % zonal average
        % take time averages
        for t = {'ann', 'djf', 'jja', 'mam', 'son'}; time = t{1};
            if strcmp(time, 'ann')
                ga_dalr_diff_t.(land).(time) = squeeze(nanmean(ga_dalr_diff0.(land), 3));
            elseif strcmp(time, 'djf')
                ga_dalr_diff_shift.(land) = circshift(ga_dalr_diff0.(land), 1, 3);
                ga_dalr_diff_t.(land).(time) = squeeze(nanmean(ga_dalr_diff_shift.(land)(:,:,1:3,:), 3));
            elseif strcmp(time, 'jja')
                ga_dalr_diff_t.(land).(time) = squeeze(nanmean(ga_dalr_diff0.(land)(:,:,6:8,:), 3));
            elseif strcmp(time, 'mam')
                ga_dalr_diff_t.(land).(time) = squeeze(nanmean(ga_dalr_diff0.(land)(:,:,3:5,:), 3));
            elseif strcmp(time, 'son')
                ga_dalr_diff_t.(land).(time) = squeeze(nanmean(ga_dalr_diff0.(land)(:,:,9:11,:), 3));
            end
        end
    end

    % save filtered data
    printname = [foldername 'ga_dalr_diff_si_mon_lat'];
    if ~exist(foldername, 'dir')
        mkdir(foldername)
    end
    save(printname, 'ga_dalr_diff', 'lat');

    printname = [foldername 'ga_dalr_diff_si_lon_lat'];
    if ~exist(foldername, 'dir')
        mkdir(foldername)
    end
    save(printname, 'ga_dalr_diff_t', 'lat');
end % compute mon x lat gamma percentage difference field with land/ocean masking
function proc_ga_dalr_bl_diff_si_mon_lat(type, par)
    if strcmp(type, 'era5') | strcmp(type, 'erai')
        foldername = sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/', type, par.lat_interp);
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s', type);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s', type);
        file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/%s/temp/%s_temp_%s.ymonmean.nc', type, type, par.(type).yr_span));
        fullpath=sprintf('%s/%s', file.folder, file.name);
        ta_orig = ncread(fullpath, 't');
    elseif strcmp(type, 'gcm')
        foldername = sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s/%s/', type, par.model, par.gcm.clim, par.lat_interp);
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s/%s', type, par.model, par.gcm.clim);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s', type, par.model, par.gcm.clim);
        var = 'ta';
        file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/gcm/%s/%s_Amon_%s_%s_r1i1p1_*.ymonmean.nc', par.model, var, par.model, par.gcm.clim));
        fullpath=sprintf('%s/%s', file.folder, file.name);
        ta_orig = ncread(fullpath, var);
    elseif strcmp(type, 'echam')
        foldername = sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s/', type, par.echam.clim, par.lat_interp);
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s', type, par.echam.clim);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s', type, par.echam.clim);
        var = 't';
        file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/echam/ATM_rjg_%s_*.ymonmean.nc', par.echam.clim));
        fullpath=sprintf('%s/%s', file.folder, file.name);
        ta_orig = double(ncread(fullpath, var));
    end
    load(sprintf('%s/grid.mat', prefix)); % read grid data
    load(sprintf('%s/dtdzsi.mat', prefix)); dtdzzsi = dtdzsi; clear dtdzsi; % read temp in si coordinates
    % load(sprintf('%s/dtdzzsi.mat', prefix)); % read temp in si coordinates
    load(sprintf('%s/%s/masks.mat', prefix_proc, par.lat_interp)); % load land and ocean masks

    dtmdzzsi = 1e3*par.g/par.cpd*ones(size(dtdzzsi));

    lat = par.lat_std;

    ga_dalr_bl_diff_orig = (dtmdzzsi - dtdzzsi)./dtmdzzsi * 1e2; % percentage difference of moist adiabatic vs GCM lapse rates
    clear dtmdz dtdz;
    ga_dalr_bl_diff_orig = permute(ga_dalr_bl_diff_orig,[2 1 3 4]); % bring lat to 1st
    ga_dalr_bl_diff_orig = interp1(grid.dim3.lat, ga_dalr_bl_diff_orig, lat); % interpolate to standard lat
    ga_dalr_bl_diff_orig = permute(ga_dalr_bl_diff_orig, [3 2 1 4]); % bring height front
    ga_dalr_bl_diff_orig = interp1(grid.dim3.si, ga_dalr_bl_diff_orig, linspace(1,par.si_bl,100)); % prepare to average between 1000-200 hPa
    ga_dalr_bl_diff_orig = squeeze(nanmean(ga_dalr_bl_diff_orig,1)); % take vertical average

    ga_dalr_bl_diff0.lo = ga_dalr_bl_diff_orig;
    ga_dalr_bl_diff0.l = ga_dalr_bl_diff0.lo.*mask.ocean; % filter ga_dalr_bl_diff0 with surface mask
    ga_dalr_bl_diff0.o = ga_dalr_bl_diff0.lo.*mask.land; % filter ga_dalr_bl_diff0 with surface mask

    for l = {'lo', 'l', 'o'}; land = l{1}; % over land, over ocean, or both
        ga_dalr_bl_diff.(land)= squeeze(nanmean(ga_dalr_bl_diff0.(land), 1)); % zonal average
        % take time averages
        for t = {'ann', 'djf', 'jja', 'mam', 'son'}; time = t{1};
            if strcmp(time, 'ann')
                ga_dalr_bl_diff_t.(land).(time) = squeeze(nanmean(ga_dalr_bl_diff0.(land), 3));
            elseif strcmp(time, 'djf')
                ga_dalr_bl_diff_shift.(land) = circshift(ga_dalr_bl_diff0.(land), 1, 3);
                ga_dalr_bl_diff_t.(land).(time) = squeeze(nanmean(ga_dalr_bl_diff_shift.(land)(:,:,1:3,:), 3));
            elseif strcmp(time, 'jja')
                ga_dalr_bl_diff_t.(land).(time) = squeeze(nanmean(ga_dalr_bl_diff0.(land)(:,:,6:8,:), 3));
            elseif strcmp(time, 'mam')
                ga_dalr_bl_diff_t.(land).(time) = squeeze(nanmean(ga_dalr_bl_diff0.(land)(:,:,3:5,:), 3));
            elseif strcmp(time, 'son')
                ga_dalr_bl_diff_t.(land).(time) = squeeze(nanmean(ga_dalr_bl_diff0.(land)(:,:,9:11,:), 3));
            end
        end
    end

    % save filtered data
    printname = [foldername 'ga_dalr_bl_diff_si_mon_lat'];
    if ~exist(foldername, 'dir')
        mkdir(foldername)
    end
    save(printname, 'ga_dalr_bl_diff', 'lat');

    printname = [foldername 'ga_dalr_bl_diff_si_lon_lat'];
    if ~exist(foldername, 'dir')
        mkdir(foldername)
    end
    save(printname, 'ga_dalr_bl_diff_t', 'lat');
end % compute mon x lat gamma percentage difference field with land/ocean masking
function proc_ga_trop_malr_diff_si_mon_lat(type, par)
    if strcmp(type, 'era5') | strcmp(type, 'erai')
        foldername = sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/', type, par.lat_interp);
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s', type);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s', type);
        file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/%s/temp/%s_temp_%s.ymonmean.nc', type, type, par.(type).yr_span));
        fullpath=sprintf('%s/%s', file.folder, file.name);
        ta_orig = ncread(fullpath, 't');
    elseif strcmp(type, 'gcm')
        foldername = sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s/%s/', type, par.model, par.gcm.clim, par.lat_interp);
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s/%s', type, par.model, par.gcm.clim);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s', type, par.model, par.gcm.clim);
        var = 'ta';
        file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/gcm/%s/%s_Amon_%s_%s_r1i1p1_*.ymonmean.nc', par.model, var, par.model, par.gcm.clim));
        fullpath=sprintf('%s/%s', file.folder, file.name);
        ta_orig = ncread(fullpath, var);
    elseif strcmp(type, 'echam')
        foldername = sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s/', type, par.echam.clim, par.lat_interp);
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s', type, par.echam.clim);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s', type, par.echam.clim);
        var = 't';
        file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/echam/ATM_rjg_%s_*.ymonmean.nc', par.echam.clim));
        fullpath=sprintf('%s/%s', file.folder, file.name);
        ta_orig = double(ncread(fullpath, var));
    end
    load(sprintf('%s/grid.mat', prefix)); % read grid data
    load(sprintf('%s/dtdzzsi.mat', prefix)); % read temp in si coordinates
    load(sprintf('%s/malrzsi.mat', prefix)); % read temp in si coordinates
    load(sprintf('%s/%s/masks.mat', prefix_proc, par.lat_interp)); % load land and ocean masks

    lat = par.lat_std;

    ga_malr_diff_orig = (dtmdzzsi - dtdzzsi)./dtmdzzsi * 1e2; % percentage difference of moist adiabatic vs GCM lapse rates
    clear dtmdz dtdz;
    ga_malr_diff_orig = permute(ga_malr_diff_orig,[2 1 3 4]); % bring lat to 1st
    ga_malr_diff_orig = interp1(grid.dim3.lat, ga_malr_diff_orig, lat); % interpolate to standard lat
    ga_malr_diff_orig = permute(ga_malr_diff_orig, [3 2 1 4]); % bring height front
    ga_malr_diff_orig = interp1(grid.dim3.si, ga_malr_diff_orig, linspace(par.si_bl,par.si_up,100)); % prepare to average between 1000-200 hPa
    ga_malr_diff_orig = squeeze(nanmean(ga_malr_diff_orig,1)); % take vertical average

    ga_malr_diff0.lo = ga_malr_diff_orig;
    ga_malr_diff0.l = ga_malr_diff0.lo.*mask.ocean; % filter ga_malr_diff0 with surface mask
    ga_malr_diff0.o = ga_malr_diff0.lo.*mask.land; % filter ga_malr_diff0 with surface mask

    for l = {'lo', 'l', 'o'}; land = l{1}; % over land, over ocean, or both
        ga_malr_diff.(land)= squeeze(nanmean(ga_malr_diff0.(land), 1)); % zonal average
        % take time averages
        for t = {'ann', 'djf', 'jja', 'mam', 'son'}; time = t{1};
            if strcmp(time, 'ann')
                ga_malr_diff_t.(land).(time) = squeeze(nanmean(ga_malr_diff0.(land), 3));
            elseif strcmp(time, 'djf')
                ga_malr_diff_shift.(land) = circshift(ga_malr_diff0.(land), 1, 3);
                ga_malr_diff_t.(land).(time) = squeeze(nanmean(ga_malr_diff_shift.(land)(:,:,1:3,:), 3));
            elseif strcmp(time, 'jja')
                ga_malr_diff_t.(land).(time) = squeeze(nanmean(ga_malr_diff0.(land)(:,:,6:8,:), 3));
            elseif strcmp(time, 'mam')
                ga_malr_diff_t.(land).(time) = squeeze(nanmean(ga_malr_diff0.(land)(:,:,3:5,:), 3));
            elseif strcmp(time, 'son')
                ga_malr_diff_t.(land).(time) = squeeze(nanmean(ga_malr_diff0.(land)(:,:,9:11,:), 3));
            end
        end
    end

    % save filtered data
    printname = [foldername 'ga_malr_diff_si_mon_lat'];
    if ~exist(foldername, 'dir')
        mkdir(foldername)
    end
    save(printname, 'ga_malr_diff', 'lat');

    printname = [foldername 'ga_malr_diff_si_lon_lat'];
    if ~exist(foldername, 'dir')
        mkdir(foldername)
    end
    save(printname, 'ga_malr_diff_t', 'lat');
end % compute mon x lat gamma percentage difference field with land/ocean masking
function make_tai(type, par) % add surface data to temperature and interpolate to hi-res grid
    if strcmp(type, 'era5') | strcmp(type, 'erai')
        foldername = sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/', type, par.lat_interp);
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s', type);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s', type);
        file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/%s/temp/%s_temp_%s.ymonmean.nc', type, type, par.(type).yr_span));
        fullpath=sprintf('%s/%s', file.folder, file.name);
        temp = ncread(fullpath, 't');
        file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/%s/zg/%s_zg_%s.ymonmean.nc', type, type, par.(type).yr_span));
        fullpath=sprintf('%s/%s', file.folder, file.name);
        zg = ncread(fullpath, 'z'); zg = zg/par.g; % convert from geopotential to height in m
    elseif strcmp(type, 'gcm')
        foldername = sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s/%s/', type, par.model, par.gcm.clim, par.lat_interp);
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s/%s', type, par.model, par.gcm.clim);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s', type, par.model, par.gcm.clim);
        var = 'ta';
        file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/gcm/%s/%s_Amon_%s_%s_r1i1p1_*.ymonmean.nc', par.model, var, par.model, par.gcm.clim));
        fullpath=sprintf('%s/%s', file.folder, file.name);
        temp = ncread(fullpath, var);
        var = 'zg';
        file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/gcm/%s/%s_Amon_%s_%s_r1i1p1_*.ymonmean.nc', par.model, var, par.model, par.gcm.clim));
        fullpath=sprintf('%s/%s', file.folder, file.name);
        zg = ncread(fullpath, var);
    elseif strcmp(type, 'echam')
        foldername = sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s/', type, par.echam.clim, par.lat_interp);
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s', type, par.echam.clim);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s', type, par.echam.clim);
        var = 't';
        file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/echam/ATM_rjg_%s_*.ymonmean.nc', par.echam.clim));
        fullpath=sprintf('%s/%s', file.folder, file.name);
        temp = double(ncread(fullpath, var));
        var = 'geopoth';
        file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/echam/ATM_rjg_%s_*.ymonmean.nc', par.echam.clim));
        fullpath=sprintf('%s/%s', file.folder, file.name);
        zg = double(ncread(fullpath, var));
    end
    load(sprintf('%s/grid.mat', prefix)); % read grid data
    load(sprintf('%s/srfc.mat', prefix)); % read surface variable data

    % create surface mask
    if strcmp(type, 'era5') | strcmp(type, 'erai')
        ps_vert = repmat(srfc.sp, [1 1 1 size(temp, 3)]); % dims (lon x lat x time x plev)
        ps_vert = permute(ps_vert, [1 2 4 3]); % dims (lon x lat x plev x time)
        pa = double(permute(repmat(grid.dim3.plev, [1 size(srfc.sp)]), [2 3 1 4]));
        ts_vert = repmat(srfc.t2m, [1 1 1 size(temp, 3)]); % dims (lon x lat x time x plev)
        ts_vert = permute(ts_vert, [1 2 4 3]); % dims (lon x lat x plev x time)
        zs_vert = repmat(srfc.zs, [1 1 1 size(temp, 3)]); % dims (lon x lat x time x plev)
        zs_vert = permute(zs_vert, [1 2 4 3]); % dims (lon x lat x plev x time)
    elseif strcmp(type, 'gcm')
        ps_vert = repmat(srfc.ps, [1 1 1 size(temp, 3)]); % dims (lon x lat x time x plev)
        ps_vert = permute(ps_vert, [1 2 4 3]); % dims (lon x lat x plev x time)
        pa = permute(repmat(grid.dim3.plev, [1 size(srfc.ps)]), [2 3 1 4]);
        ts_vert = repmat(srfc.tas, [1 1 1 size(temp, 3)]); % dims (lon x lat x time x plev)
        ts_vert = permute(ts_vert, [1 2 4 3]); % dims (lon x lat x plev x time)
        zs_vert = repmat(srfc.zs, [1 1 1 size(temp, 3)]); % dims (lon x lat x time x plev)
        zs_vert = permute(zs_vert, [1 2 4 3]); % dims (lon x lat x plev x time)
    elseif strcmp(type, 'echam')
        ps_vert = repmat(srfc.aps, [1 1 1 size(temp, 3)]); % dims (lon x lat x time x plev)
        ps_vert = permute(ps_vert, [1 2 4 3]); % dims (lon x lat x plev x time)
        pa = permute(repmat(grid.dim3.plev, [1 size(srfc.aps)]), [2 3 1 4]);
        ts_vert = repmat(srfc.temp2, [1 1 1 size(temp, 3)]); % dims (lon x lat x time x plev)
        ts_vert = permute(ts_vert, [1 2 4 3]); % dims (lon x lat x plev x time)
        zs_vert = repmat(srfc.zs, [1 1 1 size(temp, 3)]); % dims (lon x lat x time x plev)
        zs_vert = permute(zs_vert, [1 2 4 3]); % dims (lon x lat x plev x time)
    end
    surface_mask = nan(size(temp));
    surface_mask(pa < ps_vert) = 1;

    temp_sm.lo = temp.*surface_mask; % filter temp with surface mask
    zg_sm.lo = zg.*surface_mask; % filter zg with surface mask

    % add tsurf data and interpolate to higher resolution vertical grid
    [pa_plus ta_plus zg_plus] = deal(nan([size(pa,1), size(pa,2) size(pa,3)+1 size(pa,4)])); % create empty grid with one extra vertical level
    pa_plus(:,:,1:end-1,:) = pa; % populate with standard pressure grid
    ta_plus(:,:,1:end-1,:) = temp_sm.lo; % populate with standard atmospheric temperature
    zg_plus(:,:,1:end-1,:) = zg_sm.lo; % populate with szgndard atmospheric temperature
    pa_plus(:,:,end,:) = ps_vert(:,:,1,:); % add surface pressure data into standard pressure grid
    ta_plus(:,:,end,:) = ts_vert(:,:,1,:); % add surface temperature data
    zg_plus(:,:,end,:) = zs_vert(:,:,1,:); % add surface temperature data
    pa_plus = permute(pa_plus, [3 1 2 4]); % bring plev dimension to front
    ta_plus = permute(ta_plus, [3 1 2 4]); % bring plev dimension to front
    zg_plus = permute(zg_plus, [3 1 2 4]); % bring plev dimension to front
    [pa_plus sort_index] = sort(pa_plus, 1, 'descend'); % sort added surface pressure such that pressure decreases monotonically
    tai_sm.lo = nan(length(par.pa), size(pa, 1), size(pa, 2), size(pa, 4));
    zgi_sm.lo = nan(length(par.pa), size(pa, 1), size(pa, 2), size(pa, 4));
    pb = CmdLineProgressBar("Sorting and interpolating temperature to new standard grid...");
    for ilon=1:size(pa_plus,2)
        pb.print(ilon, size(pa_plus,2));
        for ilat=1:size(pa_plus,3)
            for time=1:size(pa_plus,4)
                ta_plus(:,ilon,ilat,time) = ta_plus(sort_index(:,ilon,ilat,time),ilon,ilat,time); % sort temperature (has to be in loop because sort_index works for vector calls only)
                zg_plus(:,ilon,ilat,time) = zg_plus(sort_index(:,ilon,ilat,time),ilon,ilat,time); % sort temperature (has to be in loop because sort_index works for vector calls only)

                tapanan = squeeze(pa_plus(:,ilon,ilat,time));
                tananfi = isnan(squeeze(pa_plus(:,ilon,ilat,time))) | isnan(squeeze(ta_plus(:,ilon,ilat,time)));
                tanan = squeeze(ta_plus(:,ilon,ilat,time));
                zgpanan = squeeze(pa_plus(:,ilon,ilat,time));
                zgnanfi = isnan(squeeze(pa_plus(:,ilon,ilat,time))) | isnan(squeeze(zg_plus(:,ilon,ilat,time)));
                zgnan = squeeze(zg_plus(:,ilon,ilat,time));

                tapanan(tananfi) = [];
                tanan(tananfi) = [];
                zgpanan(zgnanfi) = [];
                zgnan(zgnanfi) = [];

                tai_sm.lo(:,ilon,ilat,time) = interp1(tapanan, tanan, par.pa, 'linear', nan); % interpolate to higher resolution vertical grid
                zgi_sm.lo(:,ilon,ilat,time) = interp1(zgpanan, zgnan, par.pa, 'linear', nan); % interpolate to higher resolution vertical grid
            end
        end
    end
    clear pa_plus ta_plus zg_plus; % clear unneeded variables
    tai_sm.lo = permute(tai_sm.lo, [2 3 1 4]); % bring plev back to third dimension
    tai = tai_sm.lo;
    zgi_sm.lo = permute(zgi_sm.lo, [2 3 1 4]); % bring plev back to third dimension
    zgi = zgi_sm.lo;

    if any(strcmp(type, {'era5', 'erai'})); newdir=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s', type);
    elseif strcmp(type, 'gcm'); newdir=sprintf('/project2/tas1/miyawaki/projects/002/data/read/gcm/%s/%s', par.model, par.gcm.clim);
    elseif strcmp(type, 'echam'); newdir=sprintf('/project2/tas1/miyawaki/projects/002/data/read/echam/%s', par.echam.clim); end;
    if ~exist(newdir, 'dir'); mkdir(newdir); end
    filename='tai.mat';
    save(sprintf('%s/%s', newdir, filename), 'tai', 'zgi', '-v7.3');

end
function make_dtdzsi(type, par) % compute model lapse rate in lat x plev x mon (add surface data after converting to sigma)
% calculate moist adiabats
    if strcmp(type, 'era5') | strcmp(type, 'erai')
        foldername = sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/', type, par.lat_interp);
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s', type);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s', type);
        file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/%s/temp/%s_temp_%s.ymonmean.nc', type, type, par.(type).yr_span));
        fullpath=sprintf('%s/%s', file.folder, file.name);
        temp = ncread(fullpath, 't');
        file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/%s/zg/%s_zg_%s.ymonmean.nc', type, type, par.(type).yr_span));
        fullpath=sprintf('%s/%s', file.folder, file.name);
        zg = ncread(fullpath, 'z'); zg = zg/par.g; % convert from geopotential to height in m
    elseif strcmp(type, 'gcm')
        foldername = sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s/%s/', type, par.model, par.gcm.clim, par.lat_interp);
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s/%s', type, par.model, par.gcm.clim);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s', type, par.model, par.gcm.clim);
        var = 'ta';
        file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/gcm/%s/%s_Amon_%s_%s_r1i1p1_*.ymonmean.nc', par.model, var, par.model, par.gcm.clim));
        fullpath=sprintf('%s/%s', file.folder, file.name);
        temp = ncread(fullpath, var);
        var = 'zg';
        file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/gcm/%s/%s_Amon_%s_%s_r1i1p1_*.ymonmean.nc', par.model, var, par.model, par.gcm.clim));
        fullpath=sprintf('%s/%s', file.folder, file.name);
        zg = ncread(fullpath, var);
    elseif strcmp(type, 'echam')
        foldername = sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s/', type, par.echam.clim, par.lat_interp);
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s', type, par.echam.clim);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s', type, par.echam.clim);
        var = 't';
        file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/echam/ATM_rjg_%s_*.ymonmean.nc', par.echam.clim));
        fullpath=sprintf('%s/%s', file.folder, file.name);
        temp = double(ncread(fullpath, var));
        var = 'geopoth';
        file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/echam/ATM_rjg_%s_*.ymonmean.nc', par.echam.clim));
        fullpath=sprintf('%s/%s', file.folder, file.name);
        zg = double(ncread(fullpath, var));
    end
    load(sprintf('%s/grid.mat', prefix)); % read grid data
    load(sprintf('%s/srfc.mat', prefix)); % read surface variable data
    % load(sprintf('%s/orog.mat', prefix)); % read surface orography

    % create surface mask
    if strcmp(type, 'era5') | strcmp(type, 'erai')
        ps_vert = repmat(srfc.sp, [1 1 1 size(temp, 3)]); % dims (lon x lat x time x plev)
        ps_vert = permute(ps_vert, [1 2 4 3]); % dims (lon x lat x plev x time)
        pa = double(permute(repmat(grid.dim3.plev, [1 size(srfc.sp)]), [2 3 1 4]));
    elseif strcmp(type, 'gcm')
        ps_vert = repmat(srfc.ps, [1 1 1 size(temp, 3)]); % dims (lon x lat x time x plev)
        ps_vert = permute(ps_vert, [1 2 4 3]); % dims (lon x lat x plev x time)
        pa = permute(repmat(grid.dim3.plev, [1 size(srfc.ps)]), [2 3 1 4]);
    elseif strcmp(type, 'echam')
        ps_vert = repmat(srfc.aps, [1 1 1 size(temp, 3)]); % dims (lon x lat x time x plev)
        ps_vert = permute(ps_vert, [1 2 4 3]); % dims (lon x lat x plev x time)
        pa = permute(repmat(grid.dim3.plev, [1 size(srfc.aps)]), [2 3 1 4]);
    end
    surface_mask = nan(size(pa));
    surface_mask(pa < ps_vert) = 1;

    ta_sm = temp .* surface_mask;
    zg_sm = zg .* surface_mask;

    ps_vert = permute(ps_vert, [3 1 2 4]);
    pa = permute(pa, [3 1 2 4]);
    ta_sm = permute(ta_sm, [3 1 2 4]);
    zg_sm = permute(zg_sm, [3 1 2 4]);
    % convert to sigma coord.
    pb = CmdLineProgressBar("Sorting and interpolating dtdz to new standard grid...");
    for lo=1:size(pa,2)
        pb.print(lo, size(pa,2));
        for la=1:size(pa,3)
            for mo=1:size(pa,4)
                ta_si(:,lo,la,mo) = interp1(pa(:,lo,la,mo)./ps_vert(1,lo,la,mo), ta_sm(:,lo,la,mo), grid.dim3.si, 'linear');
                zg_si(:,lo,la,mo) = interp1(pa(:,lo,la,mo)./ps_vert(1,lo,la,mo), zg_sm(:,lo,la,mo), grid.dim3.si, 'linear');

                % add surface data
                ta_si(1,lo,la,mo) = srfc.tas(lo,la,mo);
                zg_si(1,lo,la,mo) = srfc.zs(lo,la,mo);

                % only keep nonnan data and redo interpolation
                notnan = find(~isnan(squeeze(ta_si(:,lo,la,mo))) & ~isnan(squeeze(zg_si(:,lo,la,mo))));
                ta_si(:,lo,la,mo) = interp1(grid.dim3.si(notnan), ta_si(notnan,lo,la,mo), grid.dim3.si);
                zg_si(:,lo,la,mo) = interp1(grid.dim3.si(notnan), zg_si(notnan,lo,la,mo), grid.dim3.si);

            end
        end
    end
    clear ta_sm zg_sm

    % % interpolate to higher resolution grid
    % ta_si = interp1(grid.dim3.si, ta_si, )

    % calculate lapse rate before taking zonal average
    dtdzsi = -1e3*(ta_si(2:end,:,:,:)-ta_si(1:end-1,:,:,:))./(zg_si(2:end,:,:,:)-zg_si(1:end-1,:,:,:)); % lapse rate in K/km
    dtdzsi = interp1(1/2*(grid.dim3.si(2:end)+grid.dim3.si(1:end-1)), dtdzsi, grid.dim3.si);
    dtdzsi(1,:,:,:) = -1e3*(ta_si(2,:,:,:)-ta_si(1,:,:,:))./(zg_si(2,:,:,:)-zg_si(1,:,:,:));
    dtdzsi = permute(dtdzsi, [2 3 1 4]);

    if any(strcmp(type, {'era5', 'erai'})); newdir=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s', type);
    elseif strcmp(type, 'gcm'); newdir=sprintf('/project2/tas1/miyawaki/projects/002/data/read/gcm/%s/%s', par.model, par.gcm.clim);
    elseif strcmp(type, 'echam'); newdir=sprintf('/project2/tas1/miyawaki/projects/002/data/read/echam/%s', par.echam.clim); end;
    if ~exist(newdir, 'dir'); mkdir(newdir); end
    filename='dtdzsi.mat';
    save(sprintf('%s/%s', newdir, filename), 'dtdzsi', '-v7.3');

end
function make_dtdzzsi_old(type, par) % compute model lapse rate in lat x plev x mon (add surface data before converting to sigma)
% calculate moist adiabats
    if strcmp(type, 'era5') | strcmp(type, 'erai')
        foldername = sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/', type, par.lat_interp);
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s', type);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s', type);
        file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/%s/temp/%s_temp_%s.ymonmean.nc', type, type, par.(type).yr_span));
        fullpath=sprintf('%s/%s', file.folder, file.name);
        temp = ncread(fullpath, 't');
        file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/%s/zg/%s_zg_%s.ymonmean.nc', type, type, par.(type).yr_span));
        fullpath=sprintf('%s/%s', file.folder, file.name);
        zg = ncread(fullpath, 'z'); zg = zg/par.g; % convert from geopotential to height in m
    elseif strcmp(type, 'gcm')
        foldername = sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s/%s/', type, par.model, par.gcm.clim, par.lat_interp);
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s/%s', type, par.model, par.gcm.clim);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s', type, par.model, par.gcm.clim);
        var = 'ta';
        file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/gcm/%s/%s_Amon_%s_%s_r1i1p1_*.ymonmean.nc', par.model, var, par.model, par.gcm.clim));
        fullpath=sprintf('%s/%s', file.folder, file.name);
        temp = ncread(fullpath, var);
        var = 'zg';
        file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/gcm/%s/%s_Amon_%s_%s_r1i1p1_*.ymonmean.nc', par.model, var, par.model, par.gcm.clim));
        fullpath=sprintf('%s/%s', file.folder, file.name);
        zg = ncread(fullpath, var);
    elseif strcmp(type, 'echam')
        foldername = sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s/', type, par.echam.clim, par.lat_interp);
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s', type, par.echam.clim);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s', type, par.echam.clim);
        var = 't';
        file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/echam/ATM_rjg_%s_*.ymonmean.nc', par.echam.clim));
        fullpath=sprintf('%s/%s', file.folder, file.name);
        temp = double(ncread(fullpath, var));
        var = 'geopoth';
        file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/echam/ATM_rjg_%s_*.ymonmean.nc', par.echam.clim));
        fullpath=sprintf('%s/%s', file.folder, file.name);
        zg = double(ncread(fullpath, var));
    end
    load(sprintf('%s/grid.mat', prefix)); % read grid data
    load(sprintf('%s/srfc.mat', prefix)); % read surface variable data
    % load(sprintf('%s/orog.mat', prefix)); % read surface orography

    % create surface mask
    if strcmp(type, 'era5') | strcmp(type, 'erai')
        ps_vert = repmat(srfc.sp, [1 1 1 size(temp, 3)]); % dims (lon x lat x time x plev)
        ps_vert = permute(ps_vert, [1 2 4 3]); % dims (lon x lat x plev x time)
        pa = double(permute(repmat(grid.dim3.plev, [1 size(srfc.sp)]), [2 3 1 4]));
    elseif strcmp(type, 'gcm')
        ps_vert = repmat(srfc.ps, [1 1 1 size(temp, 3)]); % dims (lon x lat x time x plev)
        ps_vert = permute(ps_vert, [1 2 4 3]); % dims (lon x lat x plev x time)
        pa = permute(repmat(grid.dim3.plev, [1 size(srfc.ps)]), [2 3 1 4]);
    elseif strcmp(type, 'echam')
        ps_vert = repmat(srfc.aps, [1 1 1 size(temp, 3)]); % dims (lon x lat x time x plev)
        ps_vert = permute(ps_vert, [1 2 4 3]); % dims (lon x lat x plev x time)
        pa = permute(repmat(grid.dim3.plev, [1 size(srfc.aps)]), [2 3 1 4]);
    end
    surface_mask = nan(size(pa));
    surface_mask(pa < ps_vert) = 1;

    ta_sm = temp .* surface_mask;
    zg_sm = zg .* surface_mask;

    % add tsurf data and interpolate to higher resolution vertical grid
    [pa_plus ta_plus zg_plus] = deal(nan([size(pa,1), size(pa,2) size(pa,3)+1 size(pa,4)])); % create empty grid with one extra vertical level
    pa_plus(:,:,1:end-1,:) = pa; % populate with standard pressure grid
    ta_plus(:,:,1:end-1,:) = ta_sm; % populate with standard atmospheric temperature
    zg_plus(:,:,1:end-1,:) = zg_sm; % populate with szgndard atmospheric temperature
    pa_plus(:,:,end,:) = ps_vert(:,:,1,:); % add surface pressure data into standard pressure grid
    if any(strcmp(type, {'era5', 'erai'})); ta_plus(:,:,end,:) = srfc.t2m(:,:,:); % add surface temperature data
    elseif strcmp(type, 'gcm'); ta_plus(:,:,end,:) = srfc.tas(:,:,:); % add surface temperature data
    elseif strcmp(type, 'echam'); ta_plus(:,:,end,:) = srfc.temp2(:,:,:); end % add surface temperature data
    % zg_plus(:,:,end,:) = srfc.zs(:,:,:); % add surface height data
    % zg_plus(:,:,end,:) = repmat(orog, [1 1 12]); % add surface height data
    zg_plus(:,:,end,:) = nan(size(srfc.zs)); % add surface height data
    ps_vert = permute(ps_vert, [3 1 2 4]); % bring plev dimension to front
    pa_plus = permute(pa_plus, [3 1 2 4]); % bring plev dimension to front
    ta_plus = permute(ta_plus, [3 1 2 4]); % bring plev dimension to front
    zg_plus = permute(zg_plus, [3 1 2 4]); % bring plev dimension to front
    [pa_plus sort_index] = sort(pa_plus, 1, 'descend'); % sort added surface pressure such that pressure decreases monotonically
    pb = CmdLineProgressBar("Sorting temperature with surface data added...");
    for lo=1:size(pa_plus,2)
        pb.print(lo, size(pa_plus,2));
        for la=1:size(pa_plus,3)
            for mo=1:size(pa_plus,4)
                ta_plus(:,lo,la,mo) = ta_plus(sort_index(:,lo,la,mo),lo,la,mo); % sort temperature (has to be in loop because sort_index works for vector calls only)
                zg_plus(:,lo,la,mo) = zg_plus(sort_index(:,lo,la,mo),lo,la,mo); % sort temperature (has to be in loop because sort_index works for vector calls only)
            end
        end
    end

    % convert to sigma coord.
    pb = CmdLineProgressBar("Sorting and interpolating dtdzz to new standard grid...");
    for lo=1:size(pa_plus,2)
        pb.print(lo, size(pa_plus,2));
        for la=1:size(pa_plus,3)
            for mo=1:size(pa_plus,4)
                ta_si(:,lo,la,mo) = interp1(pa_plus(:,lo,la,mo)./ps_vert(1,lo,la,mo), ta_plus(:,lo,la,mo), grid.dim3.si, 'linear');
                zg_si(:,lo,la,mo) = interp1(pa_plus(:,lo,la,mo)./ps_vert(1,lo,la,mo), zg_plus(:,lo,la,mo), grid.dim3.si, 'linear');
            end
        end
    end
    clear ta_plus zg_plus

    % take zonal average before computing lapse rate
    ta_si_z = squeeze(nanmean(ta_si,2)); clear ta_si;
    zg_si_z = squeeze(nanmean(zg_si,2)); clear zg_si;

    dtdzzsi = -1e3*(ta_si_z(2:end,:,:)-ta_si_z(1:end-1,:,:))./(zg_si_z(2:end,:,:)-zg_si_z(1:end-1,:,:)); % lapse rate in K/km
    dtdzzsi = interp1(1/2*(grid.dim3.si(2:end)+grid.dim3.si(1:end-1)), dtdzzsi, grid.dim3.si);
    dtdzzsi = repmat(dtdzzsi, [1 1 1 length(grid.dim3.lon)]); % repeat in longitude

    dtdzzsi = permute(dtdzzsi, [4 2 1 3]); % bring height to 3rd

    if any(strcmp(type, {'era5', 'erai'})); newdir=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s', type);
    elseif strcmp(type, 'gcm'); newdir=sprintf('/project2/tas1/miyawaki/projects/002/data/read/gcm/%s/%s', par.model, par.gcm.clim);
    elseif strcmp(type, 'echam'); newdir=sprintf('/project2/tas1/miyawaki/projects/002/data/read/echam/%s', par.echam.clim); end;
    if ~exist(newdir, 'dir'); mkdir(newdir); end
    filename='dtdzzsi.mat';
    save(sprintf('%s/%s', newdir, filename), 'dtdzzsi', '-v7.3');

end
function make_dtdzzsi(type, par) % compute model lapse rate in lat x plev x mon (add surface data after converting to sigma)
% calculate moist adiabats
    if strcmp(type, 'era5') | strcmp(type, 'erai')
        foldername = sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/', type, par.lat_interp);
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s', type);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s', type);
        file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/%s/temp/%s_temp_%s.ymonmean.nc', type, type, par.(type).yr_span));
        fullpath=sprintf('%s/%s', file.folder, file.name);
        temp = ncread(fullpath, 't');
        file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/%s/zg/%s_zg_%s.ymonmean.nc', type, type, par.(type).yr_span));
        fullpath=sprintf('%s/%s', file.folder, file.name);
        zg = ncread(fullpath, 'z'); zg = zg/par.g; % convert from geopotential to height in m
    elseif strcmp(type, 'gcm')
        foldername = sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s/%s/', type, par.model, par.gcm.clim, par.lat_interp);
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s/%s', type, par.model, par.gcm.clim);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s', type, par.model, par.gcm.clim);
        var = 'ta';
        file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/gcm/%s/%s_Amon_%s_%s_r1i1p1_*.ymonmean.nc', par.model, var, par.model, par.gcm.clim));
        fullpath=sprintf('%s/%s', file.folder, file.name);
        temp = ncread(fullpath, var);
        var = 'zg';
        file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/gcm/%s/%s_Amon_%s_%s_r1i1p1_*.ymonmean.nc', par.model, var, par.model, par.gcm.clim));
        fullpath=sprintf('%s/%s', file.folder, file.name);
        zg = ncread(fullpath, var);
    elseif strcmp(type, 'echam')
        foldername = sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s/', type, par.echam.clim, par.lat_interp);
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s', type, par.echam.clim);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s', type, par.echam.clim);
        var = 't';
        file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/echam/ATM_rjg_%s_*.ymonmean.nc', par.echam.clim));
        fullpath=sprintf('%s/%s', file.folder, file.name);
        temp = double(ncread(fullpath, var));
        var = 'geopoth';
        file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/echam/ATM_rjg_%s_*.ymonmean.nc', par.echam.clim));
        fullpath=sprintf('%s/%s', file.folder, file.name);
        zg = double(ncread(fullpath, var));
    end
    load(sprintf('%s/grid.mat', prefix)); % read grid data
    load(sprintf('%s/srfc.mat', prefix)); % read surface variable data
    % load(sprintf('%s/orog.mat', prefix)); % read surface orography

    % create surface mask
    if strcmp(type, 'era5') | strcmp(type, 'erai')
        ps_vert = repmat(srfc.sp, [1 1 1 size(temp, 3)]); % dims (lon x lat x time x plev)
        ps_vert = permute(ps_vert, [1 2 4 3]); % dims (lon x lat x plev x time)
        pa = double(permute(repmat(grid.dim3.plev, [1 size(srfc.sp)]), [2 3 1 4]));
    elseif strcmp(type, 'gcm')
        ps_vert = repmat(srfc.ps, [1 1 1 size(temp, 3)]); % dims (lon x lat x time x plev)
        ps_vert = permute(ps_vert, [1 2 4 3]); % dims (lon x lat x plev x time)
        pa = permute(repmat(grid.dim3.plev, [1 size(srfc.ps)]), [2 3 1 4]);
    elseif strcmp(type, 'echam')
        ps_vert = repmat(srfc.aps, [1 1 1 size(temp, 3)]); % dims (lon x lat x time x plev)
        ps_vert = permute(ps_vert, [1 2 4 3]); % dims (lon x lat x plev x time)
        pa = permute(repmat(grid.dim3.plev, [1 size(srfc.aps)]), [2 3 1 4]);
    end
    surface_mask = nan(size(pa));
    surface_mask(pa < ps_vert) = 1;

    ta_sm = temp .* surface_mask;
    zg_sm = zg .* surface_mask;

    ps_vert = permute(ps_vert, [3 1 2 4]);
    pa = permute(pa, [3 1 2 4]);
    ta_sm = permute(ta_sm, [3 1 2 4]);
    zg_sm = permute(zg_sm, [3 1 2 4]);
    % convert to sigma coord.
    pb = CmdLineProgressBar("Sorting and interpolating dtdzz to new standard grid...");
    for lo=1:size(pa,2)
        pb.print(lo, size(pa,2));
        for la=1:size(pa,3)
            for mo=1:size(pa,4)
                ta_si(:,lo,la,mo) = interp1(pa(:,lo,la,mo)./ps_vert(1,lo,la,mo), ta_sm(:,lo,la,mo), grid.dim3.si, 'linear');
                zg_si(:,lo,la,mo) = interp1(pa(:,lo,la,mo)./ps_vert(1,lo,la,mo), zg_sm(:,lo,la,mo), grid.dim3.si, 'linear');
            end
        end
    end
    clear ta_sm zg_sm

    % add surface data
    ta_si(1,:,:,:) = srfc.tas;
    zg_si(1,:,:,:) = srfc.zs;

    % take zonal average before computing lapse rate
    ta_si_z = squeeze(nanmean(ta_si,2)); clear ta_si;
    zg_si_z = squeeze(nanmean(zg_si,2)); clear zg_si;

    dtdzzsi = -1e3*(ta_si_z(2:end,:,:)-ta_si_z(1:end-1,:,:))./(zg_si_z(2:end,:,:)-zg_si_z(1:end-1,:,:)); % lapse rate in K/km
    dtdzzsi = interp1(1/2*(grid.dim3.si(2:end)+grid.dim3.si(1:end-1)), dtdzzsi, grid.dim3.si);
    dtdzzsi = repmat(dtdzzsi, [1 1 1 length(grid.dim3.lon)]); % repeat in longitude
    dtdzzsi = permute(dtdzzsi, [4 2 1 3]); % bring height to 3rd

    if any(strcmp(type, {'era5', 'erai'})); newdir=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s', type);
    elseif strcmp(type, 'gcm'); newdir=sprintf('/project2/tas1/miyawaki/projects/002/data/read/gcm/%s/%s', par.model, par.gcm.clim);
    elseif strcmp(type, 'echam'); newdir=sprintf('/project2/tas1/miyawaki/projects/002/data/read/echam/%s', par.echam.clim); end;
    if ~exist(newdir, 'dir'); mkdir(newdir); end
    filename='dtdzzsi.mat';
    save(sprintf('%s/%s', newdir, filename), 'dtdzzsi', '-v7.3');

end
function make_malr(type, par) % compute moist adiabatic lapse rate at every lat, lon, time
% calculate moist adiabats
    if strcmp(type, 'era5') | strcmp(type, 'erai')
        foldername = sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/', type, par.lat_interp);
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s', type);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s', type);
    elseif strcmp(type, 'gcm')
        foldername = sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s/%s/', type, par.model, par.gcm.clim, par.lat_interp);
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s/%s', type, par.model, par.gcm.clim);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s', type, par.model, par.gcm.clim);
    elseif strcmp(type, 'echam')
        foldername = sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s/', type, par.echam.clim, par.lat_interp);
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s', type, par.echam.clim);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s', type, par.echam.clim);
    end
    load(sprintf('%s/grid.mat', prefix)); % read grid data
    load(sprintf('%s/tai.mat', prefix)); clear zgi; % read interpolated temperature

    pa = repmat(par.pa', [1 length(grid.dim3.lon) length(grid.dim3.lat) 12]);
    pa = permute(pa, [2 3 1 4]);
    pa = pa/1e2; % convert to mb

    dtmdz = nan([length(grid.dim3.lon) length(grid.dim3.lat) length(par.pa) 12]);

    % calculate MALR following equations 3-5 in Stone and Carlson (1979)
    dalr = 9.8; % K/km
    eps = 0.622;
    R = 0.287; % J/g/K
    cp = 1.005; % J/g/K
    es0 = 6.11; % mb
    T0 = 273; % K
    L = 2510 - 2.38*(tai-T0);
    es = es0*exp(eps*L/R.*(1/T0-1./tai));
    desdT = eps*L.*es./(R*tai.^2);

    dtmdz = dalr * (1+eps*L.*es./(pa.*R.*tai))./(1+(eps.*L./(cp.*pa)).*(desdT));

    if any(strcmp(type, {'era5', 'erai'})); newdir=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s', type);
    elseif strcmp(type, 'gcm'); newdir=sprintf('/project2/tas1/miyawaki/projects/002/data/read/gcm/%s/%s', par.model, par.gcm.clim);
    elseif strcmp(type, 'echam'); newdir=sprintf('/project2/tas1/miyawaki/projects/002/data/read/echam/%s', par.echam.clim); end;
    if ~exist(newdir, 'dir'); mkdir(newdir); end
    filename='malr.mat';
    save(sprintf('%s/%s', newdir, filename), 'dtmdz', '-v7.3');

end
function make_malrz(type, par) % compute moist adiabatic lapse rate at every lat, time
% calculate moist adiabats
    if strcmp(type, 'era5') | strcmp(type, 'erai')
        foldername = sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/', type, par.lat_interp);
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s', type);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s', type);
    elseif strcmp(type, 'gcm')
        foldername = sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s/%s/', type, par.model, par.gcm.clim, par.lat_interp);
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s/%s', type, par.model, par.gcm.clim);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s', type, par.model, par.gcm.clim);
    elseif strcmp(type, 'echam')
        foldername = sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s/', type, par.echam.clim, par.lat_interp);
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s', type, par.echam.clim);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s', type, par.echam.clim);
    end
    load(sprintf('%s/grid.mat', prefix)); % read grid data
    load(sprintf('%s/tai.mat', prefix)); clear zgi; % read interpolated temperature

    tai = squeeze(nanmean(tai, 1)); % zonal average

    pa = repmat(par.pa', [1 length(grid.dim3.lat) 12]);
    pa = permute(pa, [2 1 3]);
    pa = pa/1e2; % convert to mb

    dtmdzz = nan([length(grid.dim3.lat) length(par.pa) 12]);

    % calculate MALR following equations 3-5 in Stone and Carlson (1979)
    dalr = 9.8; % K/km
    eps = 0.622;
    R = 0.287; % J/g/K
    cp = 1.005; % J/g/K
    es0 = 6.11; % mb
    T0 = 273; % K
    L = 2510 - 2.38*(tai-T0);
    es = es0*exp(eps*L/R.*(1/T0-1./tai));
    desdT = eps*L.*es./(R*tai.^2);

    dtmdzz = dalr * (1+eps*L.*es./(pa.*R.*tai))./(1+(eps.*L./(cp.*pa)).*(desdT));

    dtmdzz = repmat(dtmdzz, [1 1 1 length(grid.dim3.lon)]); % repeat to longitude
    dtmdzz = permute(dtmdzz, [4 1 2 3]); % bring lon to 1st

    if any(strcmp(type, {'era5', 'erai'})); newdir=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s', type);
    elseif strcmp(type, 'gcm'); newdir=sprintf('/project2/tas1/miyawaki/projects/002/data/read/gcm/%s/%s', par.model, par.gcm.clim);
    elseif strcmp(type, 'echam'); newdir=sprintf('/project2/tas1/miyawaki/projects/002/data/read/echam/%s', par.echam.clim); end;
    if ~exist(newdir, 'dir'); mkdir(newdir); end
    filename='malrz.mat';
    save(sprintf('%s/%s', newdir, filename), 'dtmdzz', '-v7.3');

end
function make_malrsi(type, par) % compute moist adiabatic lapse rate at every lat, time
% calculate moist adiabats
    if strcmp(type, 'era5') | strcmp(type, 'erai')
        foldername = sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/', type, par.lat_interp);
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s', type);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s', type);
    elseif strcmp(type, 'gcm')
        foldername = sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s/%s/', type, par.model, par.gcm.clim, par.lat_interp);
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s/%s', type, par.model, par.gcm.clim);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s', type, par.model, par.gcm.clim);
    elseif strcmp(type, 'echam')
        foldername = sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s/', type, par.echam.clim, par.lat_interp);
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s', type, par.echam.clim);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s', type, par.echam.clim);
    end
    load(sprintf('%s/grid.mat', prefix)); % read grid data
    load(sprintf('%s/srfc.mat', prefix)); % read surface variable data
    load(sprintf('%s/tai.mat', prefix)); clear zgi; % read interpolated temperature

    pa = repmat(par.pa', [1 length(grid.dim3.lon) length(grid.dim3.lat) 12]);
    pa = permute(pa, [2 3 1 4]);
    pa = pa/1e2; % convert to mb

    dtmdz = nan([length(grid.dim3.lon) length(grid.dim3.lat) length(par.pa) 12]);

    % calculate MALR following equations 3-5 in Stone and Carlson (1979)
    dalr = 9.8; % K/km
    eps = 0.622;
    R = 0.287; % J/g/K
    cp = 1.005; % J/g/K
    es0 = 6.11; % mb
    T0 = 273; % K
    L = 2510 - 2.38*(tai-T0);
    es = es0*exp(eps*L/R.*(1/T0-1./tai));
    desdT = eps*L.*es./(R*tai.^2);

    dtmdz = dalr * (1+eps*L.*es./(pa.*R.*tai))./(1+(eps.*L./(cp.*pa)).*(desdT));

    clear L es desdT;
    % 2 m dtmdz
    L = 2510 - 2.38*(srfc.tas-T0);
    es = es0*exp(eps*L/R.*(1/T0-1./srfc.tas));
    desdT = eps*L.*es./(R*srfc.tas.^2);
    dtasmdz = dalr * (1+eps*L.*es./(1e-2*srfc.ps.*R.*srfc.tas))./(1+(eps.*L./(cp.*1e-2*srfc.ps)).*(desdT));

    if strcmp(type, 'era5') | strcmp(type, 'erai')
        ps_vert = repmat(srfc.sp, [1 1 1 size(dtmdz, 3)]); % dims (lon x lat x time x plev)
        ps_vert = permute(ps_vert, [1 2 4 3]); % dims (lon x lat x plev x time)
        pa = double(permute(repmat(par.pa', [1 size(srfc.sp)]), [2 3 1 4]));
    elseif strcmp(type, 'gcm')
        ps_vert = repmat(srfc.ps, [1 1 1 size(dtmdz, 3)]); % dims (lon x lat x time x plev)
        ps_vert = permute(ps_vert, [1 2 4 3]); % dims (lon x lat x plev x time)
        pa = permute(repmat(par.pa', [1 size(srfc.ps)]), [2 3 1 4]);
    elseif strcmp(type, 'echam')
        ps_vert = repmat(srfc.aps, [1 1 1 size(dtmdz, 3)]); % dims (lon x lat x time x plev)
        ps_vert = permute(ps_vert, [1 2 4 3]); % dims (lon x lat x plev x time)
        pa = permute(repmat(par.pa', [1 size(srfc.aps)]), [2 3 1 4]);
    end

    pa = permute(pa, [3 1 2 4]); % bring height to 1st
    dtmdz = permute(dtmdz, [3 1 2 4]); % bring height to 1st
    pb = CmdLineProgressBar("Sorting and interpolating dtmdz to new standard grid...");
    for lo=1:size(pa,2)
        pb.print(lo, size(pa,2));
        for la=1:size(pa,3)
            for mo=1:size(pa,4)
                dtmdzsi(:,lo,la,mo) = interp1(pa(:,lo,la,mo)/ps_vert(lo,la,1,mo), dtmdz(:,lo,la,mo), grid.dim3.si);
            end
        end
    end
    dtmdzsi(1,:,:,:) = dtasmdz;
    dtmdzsi = permute(dtmdzsi, [2 3 1 4]); % bring height to 3rd

    if any(strcmp(type, {'era5', 'erai'})); newdir=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s', type);
    elseif strcmp(type, 'gcm'); newdir=sprintf('/project2/tas1/miyawaki/projects/002/data/read/gcm/%s/%s', par.model, par.gcm.clim);
    elseif strcmp(type, 'echam'); newdir=sprintf('/project2/tas1/miyawaki/projects/002/data/read/echam/%s', par.echam.clim); end;
    if ~exist(newdir, 'dir'); mkdir(newdir); end
    filename='malrsi.mat';
    save(sprintf('%s/%s', newdir, filename), 'dtmdzsi', '-v7.3');

end
function make_malrzsi(type, par) % compute moist adiabatic lapse rate at every lat, time
% calculate moist adiabats
    if strcmp(type, 'era5') | strcmp(type, 'erai')
        foldername = sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/', type, par.lat_interp);
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s', type);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s', type);
    elseif strcmp(type, 'gcm')
        foldername = sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s/%s/', type, par.model, par.gcm.clim, par.lat_interp);
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s/%s', type, par.model, par.gcm.clim);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s', type, par.model, par.gcm.clim);
    elseif strcmp(type, 'echam')
        foldername = sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s/', type, par.echam.clim, par.lat_interp);
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s', type, par.echam.clim);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s', type, par.echam.clim);
    end
    load(sprintf('%s/grid.mat', prefix)); % read grid data
    load(sprintf('%s/srfc.mat', prefix)); % read surface variable data
    load(sprintf('%s/tai.mat', prefix)); clear zgi; % read interpolated temperature

    tai = squeeze(nanmean(tai, 1)); % zonal average

    pa = repmat(par.pa', [1 length(grid.dim3.lat) 12]);
    pa = permute(pa, [2 1 3]);
    pa = pa/1e2; % convert to mb

    dtmdzz = nan([length(grid.dim3.lat) length(par.pa) 12]);

    % calculate MALR following equations 3-5 in Stone and Carlson (1979)
    dalr = 9.8; % K/km
    eps = 0.622;
    R = 0.287; % J/g/K
    cp = 1.005; % J/g/K
    es0 = 6.11; % mb
    T0 = 273; % K
    L = 2510 - 2.38*(tai-T0);
    es = es0*exp(eps*L/R.*(1/T0-1./tai));
    desdT = eps*L.*es./(R*tai.^2);

    dtmdzz = dalr * (1+eps*L.*es./(pa.*R.*tai))./(1+(eps.*L./(cp.*pa)).*(desdT));

    dtmdzz = repmat(dtmdzz, [1 1 1 length(grid.dim3.lon)]); % repeat to longitude
    dtmdzz = permute(dtmdzz, [4 1 2 3]); % bring lon to 1st

    if strcmp(type, 'era5') | strcmp(type, 'erai')
        ps_vert = repmat(srfc.sp, [1 1 1 size(dtmdzz, 3)]); % dims (lon x lat x time x plev)
        ps_vert = permute(ps_vert, [1 2 4 3]); % dims (lon x lat x plev x time)
        pa = double(permute(repmat(par.pa', [1 size(srfc.sp)]), [2 3 1 4]));
    elseif strcmp(type, 'gcm')
        ps_vert = repmat(srfc.ps, [1 1 1 size(dtmdzz, 3)]); % dims (lon x lat x time x plev)
        ps_vert = permute(ps_vert, [1 2 4 3]); % dims (lon x lat x plev x time)
        pa = permute(repmat(par.pa', [1 size(srfc.ps)]), [2 3 1 4]);
    elseif strcmp(type, 'echam')
        ps_vert = repmat(srfc.aps, [1 1 1 size(dtmdzz, 3)]); % dims (lon x lat x time x plev)
        ps_vert = permute(ps_vert, [1 2 4 3]); % dims (lon x lat x plev x time)
        pa = permute(repmat(par.pa', [1 size(srfc.aps)]), [2 3 1 4]);
    end

    pa = permute(pa, [3 1 2 4]); % bring height to 1st
    dtmdzz = permute(dtmdzz, [3 1 2 4]); % bring height to 1st
    pb = CmdLineProgressBar("Sorting and interpolating dtmdzz to new standard grid...");
    for lo=1:size(pa,2)
        pb.print(lo, size(pa,2));
        for la=1:size(pa,3)
            for mo=1:size(pa,4)
                dtmdzzsi(:,lo,la,mo) = interp1(pa(:,lo,la,mo)/ps_vert(lo,la,1,mo), dtmdzz(:,lo,la,mo), grid.dim3.si);
            end
        end
    end
    dtmdzzsi = permute(dtmdzzsi, [2 3 1 4]); % bring height to 3rd

    if any(strcmp(type, {'era5', 'erai'})); newdir=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s', type);
    elseif strcmp(type, 'gcm'); newdir=sprintf('/project2/tas1/miyawaki/projects/002/data/read/gcm/%s/%s', par.model, par.gcm.clim);
    elseif strcmp(type, 'echam'); newdir=sprintf('/project2/tas1/miyawaki/projects/002/data/read/echam/%s', par.echam.clim); end;
    if ~exist(newdir, 'dir'); mkdir(newdir); end
    filename='malrzsi.mat';
    save(sprintf('%s/%s', newdir, filename), 'dtmdzzsi', '-v7.3');

end
function make_dtmdz(type, par) % compute moist adiabat at every lat, lon, time
% calculate moist adiabats
    if strcmp(type, 'era5') | strcmp(type, 'erai')
        foldername = sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/', type, par.lat_interp);
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s', type);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s', type);
    elseif strcmp(type, 'gcm')
        foldername = sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s/%s/', type, par.model, par.gcm.clim, par.lat_interp);
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s/%s', type, par.model, par.gcm.clim);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s', type, par.model, par.gcm.clim);
    elseif strcmp(type, 'echam')
        foldername = sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s/', type, par.echam.clim, par.lat_interp);
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s', type, par.echam.clim);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s', type, par.echam.clim);
    end
    load(sprintf('%s/grid.mat', prefix)); % read grid data
    load(sprintf('%s/srfc.mat', prefix)); % read surface variable data

    dtmdz = nan([length(grid.dim3.lon) length(grid.dim3.lat) length(par.pa) 12]);

    if strcmp(type, 'era5') | strcmp(type, 'erai')
        pb = CmdLineProgressBar("Calculating moist adiabats...");
        for ilon = 1:length(grid.dim3.lon);
            pb.print(ilon, length(grid.dim3.lon)); % output progress of moist adiabat calculation
            for ilat = 1:length(grid.dim3.lat);
                for imon = 1:12;
                    for fn_var = fieldnames(srfc)'
                        ima.(fn_var{1}) = srfc.(fn_var{1})(ilon, ilat, imon);
                    end
                    dtmdz(ilon, ilat, :, imon) = calc_maz_dew_pa(ima, par.pa, par, type); % compute moist adiabat with RH
                end
            end
        end
    elseif strcmp(type, 'gcm')
        pb = CmdLineProgressBar("Calculating moist adiabats...");
        for ilon = 1:length(grid.dim3.lon);
            pb.print(ilon, length(grid.dim3.lon)); % output progress of moist adiabat calculation
            for ilat = 1:length(grid.dim3.lat);
                for imon = 1:12;
                    for fn_var = fieldnames(srfc)'
                        ima.(fn_var{1}) = srfc.(fn_var{1})(ilon, ilat, imon);
                    end
                    dtmdz(ilon, ilat, :, imon) = calc_maz_hurs_pa(ima, par.pa, par); % compute moist adiabat with RH
                end
            end
        end
    end

    if any(strcmp(type, {'era5', 'erai'})); newdir=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s', type);
    elseif strcmp(type, 'gcm'); newdir=sprintf('/project2/tas1/miyawaki/projects/002/data/read/gcm/%s/%s', par.model, par.gcm.clim);
    elseif strcmp(type, 'echam'); newdir=sprintf('/project2/tas1/miyawaki/projects/002/data/read/echam/%s', par.echam.clim); end;
    if ~exist(newdir, 'dir'); mkdir(newdir); end
    filename='dtmdz.mat';
    save(sprintf('%s/%s', newdir, filename), 'dtmdz', '-v7.3');

end
function make_dtmdzz(type, par) % compute moist adiabat at every lon, time
% calculate moist adiabats
    if strcmp(type, 'era5') | strcmp(type, 'erai')
        foldername = sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/', type, par.lat_interp);
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s', type);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s', type);
    elseif strcmp(type, 'gcm')
        foldername = sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s/%s', type, par.model, par.gcm.clim, par.lat_interp);
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s/%s', type, par.model, par.gcm.clim);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s', type, par.model, par.gcm.clim);
    elseif strcmp(type, 'echam')
        foldername = sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s/', type, par.echam.clim, par.lat_interp);
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s', type, par.echam.clim);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s', type, par.echam.clim);
    end
    load(sprintf('%s/grid.mat', prefix)); % read grid data
    load(sprintf('%s/srfc.mat', prefix)); % read surface variable data

    dtmdzz = nan([length(grid.dim3.lat) length(par.pa) 12]);

    for fn_var = fieldnames(srfc)'
        srfc.(fn_var{1}) = squeeze(nanmean(srfc.(fn_var{1}),1)); % zonal average
    end

    if strcmp(type, 'era5') | strcmp(type, 'erai')
        pb = CmdLineProgressBar("Calculating moist adiabats...");
        for ilat = 1:length(grid.dim3.lat);
            pb.print(ilat, length(grid.dim3.lat)); % output progress of moist adiabat calculation
            for imon = 1:12;
                for fn_var = fieldnames(srfc)'
                    ima.(fn_var{1}) = srfc.(fn_var{1})(ilat, imon);
                end
                dtmdzz(ilat, :, imon) = calc_maz_dew_pa(ima, par.pa, par, type); % compute moist adiabat with RH
            end
        end
    elseif strcmp(type, 'gcm')
        pb = CmdLineProgressBar("Calculating moist adiabats...");
        for ilat = 1:length(grid.dim3.lat);
            pb.print(ilat, length(grid.dim3.lat)); % output progress of moist adiabat calculation
            for imon = 1:12;
                for fn_var = fieldnames(srfc)'
                    ima.(fn_var{1}) = srfc.(fn_var{1})(ilat, imon);
                end
                dtmdzz(ilat, :, imon) = calc_maz_hurs_pa(ima, par.pa, par); % compute moist adiabat with RH
            end
        end
    end

    dtmdzz = repmat(dtmdzz, [1 1 1 length(grid.dim3.lon)]); % repeat in longitude
    dtmdzz = permute(dtmdzz, [4 1 2 3]); % bring lon to 1st

    if any(strcmp(type, {'era5', 'erai'})); newdir=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s', type);
    elseif strcmp(type, 'gcm'); newdir=sprintf('/project2/tas1/miyawaki/projects/002/data/read/gcm/%s/%s', par.model, par.gcm.clim);
    elseif strcmp(type, 'echam'); newdir=sprintf('/project2/tas1/miyawaki/projects/002/data/read/echam/%s', par.echam.clim); end;
    if ~exist(newdir, 'dir'); mkdir(newdir); end
    filename='dtmdzz.mat';
    save(sprintf('%s/%s', newdir, filename), 'dtmdzz', '-v7.3');

end
function make_dtmdzzsi(type, par) % compute moist adiabat at every lon, time
% calculate moist adiabats
    if strcmp(type, 'era5') | strcmp(type, 'erai')
        foldername = sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/', type, par.lat_interp);
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s', type);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s', type);
    elseif strcmp(type, 'gcm')
        foldername = sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s/%s', type, par.model, par.gcm.clim, par.lat_interp);
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s/%s', type, par.model, par.gcm.clim);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s', type, par.model, par.gcm.clim);
    elseif strcmp(type, 'echam')
        foldername = sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s/', type, par.echam.clim, par.lat_interp);
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s', type, par.echam.clim);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s', type, par.echam.clim);
    end
    load(sprintf('%s/grid.mat', prefix)); % read grid data
    load(sprintf('%s/srfc.mat', prefix)); % read surface variable data

    dtmdzzsi = nan([length(grid.dim3.lat) length(par.pa) 12]);

    for fn_var = fieldnames(srfc)'
        srfc.(fn_var{1}) = squeeze(nanmean(srfc.(fn_var{1}),1)); % zonal average
    end

    if strcmp(type, 'era5') | strcmp(type, 'erai')
        pb = CmdLineProgressBar("Calculating moist adiabats...");
        for ilat = 1:length(grid.dim3.lat);
            pb.print(ilat, length(grid.dim3.lat)); % output progress of moist adiabat calculation
            for imon = 1:12;
                for fn_var = fieldnames(srfc)'
                    ima.(fn_var{1}) = srfc.(fn_var{1})(ilat, imon);
                end
                dtmdzzsi(ilat, :, imon) = calc_maz_dew_si(ima, par.pa, par, type, grid); % compute moist adiabat with RH
            end
        end
    elseif strcmp(type, 'gcm')
        pb = CmdLineProgressBar("Calculating moist adiabats...");
        for ilat = 1:length(grid.dim3.lat);
            pb.print(ilat, length(grid.dim3.lat)); % output progress of moist adiabat calculation
            for imon = 1:12;
                for fn_var = fieldnames(srfc)'
                    ima.(fn_var{1}) = srfc.(fn_var{1})(ilat, imon);
                end
                dtmdzzsi(ilat, :, imon) = calc_maz_hurs_si(ima, par.pa, par, grid); % compute moist adiabat with RH
            end
        end
    elseif strcmp(type, 'echam')
        pb = CmdLineProgressBar("Calculating moist adiabats...");
        for ilat = 1:length(grid.dim3.lat);
            pb.print(ilat, length(grid.dim3.lat)); % output progress of moist adiabat calculation
            for imon = 1:12;
                for fn_var = fieldnames(srfc)'
                    ima.(fn_var{1}) = srfc.(fn_var{1})(ilat, imon);
                end
                dtmdzzsi(ilat, :, imon) = calc_maz_dew_si(ima, par.pa, par, type, grid); % compute moist adiabat with RH
            end
        end
    end

    dtmdzzsi = repmat(dtmdzzsi, [1 1 1 length(grid.dim3.lon)]); % repeat in longitude
    dtmdzzsi = permute(dtmdzzsi, [4 1 2 3]); % bring lon to 1st

    if any(strcmp(type, {'era5', 'erai'})); newdir=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s', type);
    elseif strcmp(type, 'gcm'); newdir=sprintf('/project2/tas1/miyawaki/projects/002/data/read/gcm/%s/%s', par.model, par.gcm.clim);
    elseif strcmp(type, 'echam'); newdir=sprintf('/project2/tas1/miyawaki/projects/002/data/read/echam/%s', par.echam.clim); end;
    if ~exist(newdir, 'dir'); mkdir(newdir); end
    filename='dtmdzzsi.mat';
    save(sprintf('%s/%s', newdir, filename), 'dtmdzzsi', '-v7.3');

end
function make_dtmdz_lat_plev(type, par) % compute moist adiabat as function of lat
% calculate moist adiabats
    if strcmp(type, 'era5') | strcmp(type, 'erai')
        foldername = sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/', type, par.lat_interp);
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s', type);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s', type);
    elseif strcmp(type, 'gcm')
        foldername = sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s/%s/', type, par.model, par.gcm.clim, par.lat_interp);
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s/%s', type, par.model, par.gcm.clim);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s', type, par.model, par.gcm.clim);
    elseif strcmp(type, 'echam')
        foldername = sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s/', type, par.echam.clim, par.lat_interp);
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s', type, par.echam.clim);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s', type, par.echam.clim);
    end
    load(sprintf('%s/grid.mat', prefix)); % read grid data
    load(sprintf('%s/srfc.mat', prefix)); % read surface variable data

    dtmdz_zt = nan([length(grid.dim3.lat) length(par.pa)]);

    for fn_var = fieldnames(srfc)'
        srfc.(fn_var{1}) = squeeze(nanmean(nanmean(srfc.(fn_var{1}),1),3)); % zonal and time mean
    end

    if strcmp(type, 'era5') | strcmp(type, 'erai')
        pb = CmdLineProgressBar("Calculating moist adiabats...");
        for ilat = 1:length(grid.dim3.lat);
            pb.print(ilat, length(grid.dim3.lat)); % output progress of moist adiabat calculation
            for fn_var = fieldnames(srfc)'
                ima.(fn_var{1}) = srfc.(fn_var{1})(ilat);
            end
            dtmdz_zt(ilat, :) = calc_maz_dew_pa(ima, par.pa, par, type); % compute moist adiabat with RH
        end
    elseif strcmp(type, 'gcm')
        pb = CmdLineProgressBar("Calculating moist adiabats...");
        for ilat = 1:length(grid.dim3.lat);
            pb.print(ilat, length(grid.dim3.lat)); % output progress of moist adiabat calculation
            for fn_var = fieldnames(srfc)'
                ima.(fn_var{1}) = srfc.(fn_var{1})(ilat);
            end
            dtmdz_zt(ilat, :) = calc_maz_hurs_pa(ima, par.pa, par); % compute moist adiabat with RH
        end
    elseif strcmp(type, 'gcm')
        pb = CmdLineProgressBar("Calculating moist adiabats...");
        for ilat = 1:length(grid.dim3.lat);
            pb.print(ilat, length(grid.dim3.lat)); % output progress of moist adiabat calculation
            for fn_var = fieldnames(srfc)'
                ima.(fn_var{1}) = srfc.(fn_var{1})(ilat);
            end
            dtmdz_zt(ilat, :) = calc_maz_dew_pa(ima, par.pa, par, type); % compute moist adiabat with RH
        end
    end

    if any(strcmp(type, {'era5', 'erai'})); newdir=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s', type);
    elseif strcmp(type, 'gcm'); newdir=sprintf('/project2/tas1/miyawaki/projects/002/data/read/gcm/%s/%s', par.model, par.gcm.clim);
    elseif strcmp(type, 'echam'); newdir=sprintf('/project2/tas1/miyawaki/projects/002/data/read/echam/%s', par.echam.clim); end;
    if ~exist(newdir, 'dir'); mkdir(newdir); end
    filename='dtmdz_zt.mat';
    save(sprintf('%s/%s', newdir, filename), 'dtmdz_zt', '-v7.3');

end

function choose_proc_ep(type, par)
    proc_rcae(type, par) % calculate RCE and RAE regimes
    proc_rcae_alt(type, par) % calculate RCE and RAE regimes (all divergence allowed for RCE)
    % proc_temp(type, par) % calculate RCE and RAE temperature profiles
    % proc_ma(type, par) % calculate moist adiabats corresponding to RCE profiles
end % select which ep-functions to run at a time
function proc_rcae(type, par)
    for f = {'flux', 'flux_t', 'flux_z'}; ftype = f{1};
        if strcmp(type, 'era5') | strcmp(type, 'erai')
            filename = sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s.mat', type, par.lat_interp, ftype); % read ERA5 zonally averaged flux
            foldername = sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/eps_%g_ga_%g/', type, par.lat_interp, par.ep, par.ga);
            prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s', type);
            prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s', type, par.lat_interp);
        elseif strcmp(type, 'gcm')
            filename = sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s/%s/%s.mat', type, par.model, par.gcm.clim, par.lat_interp, ftype); % read gcm zonally averaged flux
            foldername = sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s/%s/eps_%g_ga_%g/', type, par.model, par.gcm.clim, par.lat_interp, par.ep, par.ga);
            prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s/%s', type, par.model, par.gcm.clim);
            prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s/%s', type, par.model, par.gcm.clim, par.lat_interp);
        elseif strcmp(type, 'echam')
            filename = sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s/%s.mat', type, par.echam.clim, par.lat_interp, ftype); % read echam zonally averaged flux
            foldername = sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s/eps_%g_ga_%g/', type, par.echam.clim, par.lat_interp, par.ep, par.ga);
            prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s', type, par.echam.clim);
            prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s', type, par.echam.clim, par.lat_interp);
        end
        if ~exist(filename); error(sprintf('Data does not exist. Please run proc_%s.m first.', ftype)); else
            load(filename);
        end

        load(sprintf('%s/masks.mat', prefix_proc)); % load land and ocean masks
        load(sprintf('%s/vh.mat', prefix_proc)); % load atmospheric heat transport
        load(sprintf('%s/vh_mon.mat', prefix_proc)); % load atmospheric heat transport

        % identify locations of RCE and RAE
        if strcmp(ftype, 'flux')
            rcae = def_rcae(type, flux, vh_mon, par); % lon x lat x time structure
            rcae_rc = def_rcae_recomp_r1(type, flux, vh_mon, par); % lon x lat x time structure
            printname = [foldername 'rcae.mat'];
            printname_rc = [foldername 'rcae_rc.mat'];
        elseif strcmp(ftype, 'flux_t')
            for l = {'lo', 'l', 'o'}; land = l{1};
                % lon x lat structure over various time averages
                for t = {'ann', 'djf', 'jja', 'mam', 'son'}; time = t{1};
                    rcae_t.(land).(time) = def_rcae(type, flux_t.(land).(time), vh.(land).(time), par);
                    rcae_rc_t.(land).(time) = def_rcae_recomp_r1(type, flux_t.(land).(time), vh.(land).(time), par);
                end % time
            end % land

            printname = [foldername 'rcae_t.mat'];
            printname_rc = [foldername 'rcae_rc_t.mat'];
        elseif strcmp(ftype, 'flux_z') % show lat x mon structure of RCAE
            for l = {'lo', 'l', 'o'}; land = l{1};
                rcae_z.(land) = def_rcae(type, flux_z.(land), vh_mon.(land), par);
                rcae_rc_z.(land) = def_rcae_recomp_r1(type, flux_z.(land), vh_mon.(land), par);
                printname = [foldername 'rcae_z.mat'];
                printname_rc = [foldername 'rcae_rc_z.mat'];
            end
        end

        % save rcae data
        if ~exist(foldername, 'dir')
            mkdir(foldername)
        end
        if strcmp(ftype, 'flux'); save(printname, 'rcae', 'lat', '-v7.3'); save(printname_rc, 'rcae_rc', 'lat', '-v7.3');
        elseif strcmp(ftype, 'flux_t'); save(printname, 'rcae_t', 'lat', '-v7.3'); save(printname_rc, 'rcae_rc_t', 'lat', '-v7.3');
        elseif strcmp(ftype, 'flux_z'); save(printname, 'rcae_z', 'lat', '-v7.3'); save(printname_rc, 'rcae_rc_z', 'lat', '-v7.3'); end
    end

end % process RCE/RAE regimes
function proc_rcae_alt(type, par)
    for f = {'flux', 'flux_t', 'flux_z'}; ftype = f{1};
        if strcmp(type, 'era5') | strcmp(type, 'erai')
            filename = sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s.mat', type, par.lat_interp, ftype); % read ERA5 zonally averaged flux
            foldername = sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/eps_%g_ga_%g/', type, par.lat_interp, par.ep, par.ga);
            prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s', type);
            prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s', type, par.lat_interp);
        elseif strcmp(type, 'gcm')
            filename = sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s/%s/%s.mat', type, par.model, par.gcm.clim, par.lat_interp, ftype); % read gcm zonally averaged flux
            foldername = sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s/%s/eps_%g_ga_%g/', type, par.model, par.gcm.clim, par.lat_interp, par.ep, par.ga);
            prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s/%s', type, par.model, par.gcm.clim);
            prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s/%s', type, par.model, par.gcm.clim, par.lat_interp);
        elseif strcmp(type, 'echam')
            filename = sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s/%s.mat', type, par.echam.clim, par.lat_interp, ftype); % read echam zonally averaged flux
            foldername = sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s/eps_%g_ga_%g/', type, par.echam.clim, par.lat_interp, par.ep, par.ga);
            prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s', type, par.echam.clim);
            prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s', type, par.echam.clim, par.lat_interp);
        end
        if ~exist(filename); error(sprintf('Data does not exist. Please run proc_%s.m first.', ftype)); else
            load(filename);
        end

        load(sprintf('%s/masks.mat', prefix_proc)); % load land and ocean masks
        load(sprintf('%s/vh.mat', prefix_proc)); % load atmospheric heat transport
        load(sprintf('%s/vh_mon.mat', prefix_proc)); % load atmospheric heat transport

        % identify locations of RCE and RAE
        if strcmp(ftype, 'flux')
            rcae_alt = def_rcae_alt(type, flux, vh_mon, par); % lon x lat x time structure
            rcae_alt_rc = def_rcae_alt_recomp_r1(type, flux, vh_mon, par); % lon x lat x time structure
            printname = [foldername 'rcae_alt.mat'];
            printname_rc = [foldername 'rcae_alt_rc.mat'];
        elseif strcmp(ftype, 'flux_t')
            for l = {'lo', 'l', 'o'}; land = l{1};
                % lon x lat structure over various time averages
                for t = {'ann', 'djf', 'jja', 'mam', 'son'}; time = t{1};
                    rcae_alt_t.(land).(time) = def_rcae_alt(type, flux_t.(land).(time), vh.(land).(time), par);
                    rcae_alt_rc_t.(land).(time) = def_rcae_alt_recomp_r1(type, flux_t.(land).(time), vh.(land).(time), par);
                end % time
            end % land

            printname = [foldername 'rcae_alt_t.mat'];
            printname_rc = [foldername 'rcae_alt_rc_t.mat'];
        elseif strcmp(ftype, 'flux_z') % show lat x mon structure of RCAE_ALT
            for l = {'lo', 'l', 'o'}; land = l{1};
                rcae_alt_z.(land) = def_rcae_alt(type, flux_z.(land), vh_mon.(land), par);
                rcae_alt_rc_z.(land) = def_rcae_alt_recomp_r1(type, flux_z.(land), vh_mon.(land), par);
                printname = [foldername 'rcae_alt_z.mat'];
                printname_rc = [foldername 'rcae_alt_rc_z.mat'];
            end
        end

        % save rcae_alt data
        if ~exist(foldername, 'dir')
            mkdir(foldername)
        end
        if strcmp(ftype, 'flux'); save(printname, 'rcae_alt', 'lat', '-v7.3'); save(printname_rc, 'rcae_alt_rc', 'lat', '-v7.3');
        elseif strcmp(ftype, 'flux_t'); save(printname, 'rcae_alt_t', 'lat', '-v7.3'); save(printname_rc, 'rcae_alt_rc_t', 'lat', '-v7.3');
        elseif strcmp(ftype, 'flux_z'); save(printname, 'rcae_alt_z', 'lat', '-v7.3'); save(printname_rc, 'rcae_alt_rc_z', 'lat', '-v7.3'); end
    end

end % process RCE/RAE regimes (all divergence allowed for RCE)
function proc_temp(type, par)
    if strcmp(type, 'era5') | strcmp(type, 'erai')
        foldername = sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/eps_%g_ga_%g/', type, par.lat_interp, par.ep, par.ga);
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s', type);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s', type);
        file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/%s/temp/%s_temp_%s.ymonmean.nc', type, type, par.(type).yr_span));
        fullpath=sprintf('%s/%s', file.folder, file.name);
        temp = ncread(fullpath, 't');
    elseif strcmp(type, 'gcm')
        foldername = sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s/%s/eps_%g_ga_%g/', type, par.model, par.gcm.clim, par.lat_interp, par.ep, par.ga);
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s/%s', type, par.model, par.gcm.clim);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s', type, par.model, par.gcm.clim);
        var = 'ta';
        file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/gcm/%s/%s_Amon_%s_%s_r1i1p1_*.ymonmean.nc', par.model, var, par.model, par.gcm.clim));
        fullpath=sprintf('%s/%s', file.folder, file.name);
        temp = ncread(fullpath, var);
    elseif strcmp(type, 'echam')
        foldername = sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s/', type, par.echam.clim, par.lat_interp);
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s', type, par.echam.clim);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s', type, par.echam.clim);
        var = 't';
        file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/echam/ATM_rjg_%s_*.ymonmean.nc', par.echam.clim));
        fullpath=sprintf('%s/%s', file.folder, file.name);
        temp = double(ncread(fullpath, var));
    end
    load(sprintf('%s/grid.mat', prefix)); % read ERA5 grid data
    load(sprintf('%s/srfc.mat', prefix)); % load surface data
    load(sprintf('%s/%s/masks.mat', prefix_proc, par.lat_interp)); % load land and ocean masks
    load(sprintf('%s/%s/eps_%g_ga_%g/rcae_t.mat', prefix_proc, par.lat_interp, par.ep, par.ga)); % load rcae data

    lat = par.lat_std;

    % interpolate temp to standard lat grid
    temp = permute(temp, [2 1 3 4]);
    temp = interp1(grid.dim3.lat, temp, lat);
    temp = permute(temp, [2 1 3 4]);

    % create surface mask
    if strcmp(type, 'era5') | strcmp(type, 'erai')
        ps = permute(srfc.sp, [2 1 3]); % bring lat to front to interpolate
        ps = interp1(grid.dim2.lat, ps, lat, 'linear', 'extrap'); % interpolate to standard grid
        ps = permute(ps, [2 1 3]); % reorder to original
        ps_vert = repmat(ps, [1 1 1 size(temp, 3)]); % dims (lon x lat x time x plev)
        ps_vert = permute(ps_vert, [1 2 4 3]); % dims (lon x lat x plev x time)
        pa = double(permute(repmat(grid.dim3.plev, [1 size(ps)]), [2 3 1 4]));
        ts = permute(srfc.t2m, [2 1 3]); % bring lat to front to interpolate
        ts = interp1(grid.dim2.lat, ts, lat); % interpolate to standard grid
        ts = permute(ts, [2 1 3]); % reorder to original
        ts_vert = repmat(ts, [1 1 1 size(temp, 3)]); % dims (lon x lat x time x plev)
        ts_vert = permute(ts_vert, [1 2 4 3]); % dims (lon x lat x plev x time)
    elseif strcmp(type, 'gcm')
        ps = permute(srfc.ps, [2 1 3]); % bring lat to front to interpolate
        ps = interp1(grid.dim2.lat, ps, lat, 'linear', 'extrap'); % interpolate to standard grid
        ps = permute(ps, [2 1 3]); % reorder to original
        ps_vert = repmat(ps, [1 1 1 size(temp, 3)]); % dims (lon x lat x time x plev)
        ps_vert = permute(ps_vert, [1 2 4 3]); % dims (lon x lat x plev x time)
        pa = permute(repmat(grid.dim3.plev, [1 size(ps)]), [2 3 1 4]);
        ts = permute(srfc.tas, [2 1 3]); % bring lat to front to interpolate
        ts = interp1(grid.dim2.lat, ts, lat); % interpolate to standard grid
        ts = permute(ts, [2 1 3]); % reorder to original
        ts_vert = repmat(ts, [1 1 1 size(temp, 3)]); % dims (lon x lat x time x plev)
        ts_vert = permute(ts_vert, [1 2 4 3]); % dims (lon x lat x plev x time)
    elseif strcmp(type, 'echam')
        ps = permute(srfc.aps, [2 1 3]); % bring lat to front to interpolate
        ps = interp1(grid.dim2.lat, ps, lat, 'linear', 'extrap'); % interpolate to standard grid
        ps = permute(ps, [2 1 3]); % reorder to original
        ps_vert = repmat(ps, [1 1 1 size(temp, 3)]); % dims (lon x lat x time x plev)
        ps_vert = permute(ps_vert, [1 2 4 3]); % dims (lon x lat x plev x time)
        pa = permute(repmat(grid.dim3.plev, [1 size(ps)]), [2 3 1 4]);
        ts = permute(srfc.temp2, [2 1 3]); % bring lat to front to interpolate
        ts = interp1(grid.dim2.lat, ts, lat); % interpolate to standard grid
        ts = permute(ts, [2 1 3]); % reorder to original
        ts_vert = repmat(ts, [1 1 1 size(temp, 3)]); % dims (lon x lat x time x plev)
        ts_vert = permute(ts_vert, [1 2 4 3]); % dims (lon x lat x plev x time)
    end
    surface_mask = nan(size(temp));
    surface_mask(pa < ps_vert) = 1;

    temp_sm.lo = temp.*surface_mask; % filter temp with surface mask

    % add tsurf data and interpolate to higher resolution vertical grid
    [pa_plus ta_plus] = deal(nan([size(pa,1), size(pa,2) size(pa,3)+1 size(pa,4)])); % create empty grid with one extra vertical level
    pa_plus(:,:,1:end-1,:) = pa; % populate with standard pressure grid
    ta_plus(:,:,1:end-1,:) = temp_sm.lo; % populate with standard atmospheric temperature
    pa_plus(:,:,end,:) = ps_vert(:,:,1,:); % add surface pressure data into standard pressure grid
    ta_plus(:,:,end,:) = ts_vert(:,:,1,:); % add surface temperature data
    pa_plus = permute(pa_plus, [3 1 2 4]); % bring plev dimension to front
    ta_plus = permute(ta_plus, [3 1 2 4]); % bring plev dimension to front
    [pa_plus sort_index] = sort(pa_plus, 1, 'descend'); % sort added surface pressure such that pressure decreases monotonically
    tai_sm.lo = nan(length(par.pa), size(pa, 1), size(pa, 2), size(pa, 4));
    pb = CmdLineProgressBar("Sorting and interpolating temperature to new standard grid...");
    for ilon=1:size(pa_plus,2)
        pb.print(ilon, size(pa_plus,2));
        for ilat=1:size(pa_plus,3)
            for time=1:size(pa_plus,4)
                ta_plus(:,ilon,ilat,time) = ta_plus(sort_index(:,ilon,ilat,time),ilon,ilat,time); % sort temperature (has to be in loop because sort_index works for vector calls only)
                tai_sm.lo(:,ilon,ilat,time) = interp1(pa_plus(:,ilon,ilat,time), ta_plus(:,ilon,ilat,time), par.pa, 'linear'); % interpolate to higher resolution vertical grid
            end
        end
    end
    clear pa_plus ta_plus; % clear unneeded variables
    tai_sm.lo = permute(tai_sm.lo, [2 3 1 4]); % bring plev back to third dimension

    % apply masks to surface pressure and RCAE regimes
    ps_n.lo = ps;
    ps_n.l = ps .* mask.ocean;
    ps_n.o = ps .* mask.land;

    mask.land_vert = repmat(mask.land, [1 1 1 size(temp, 3)]); % expand land mask to vertical dim
    mask.land_vert = permute(mask.land, [1 2 4 3]); % place vertical dim where it belongs
    mask.ocean_vert = repmat(mask.ocean, [1 1 1 size(temp, 3)]); % expand ocean mask to vertical dim
    mask.ocean_vert = permute(mask.ocean, [1 2 4 3]); % place vertical dim where it belongs

    temp_sm.l = temp.*surface_mask.*mask.ocean_vert; % filter temp with surface mask
    temp_sm.o = temp.*surface_mask.*mask.land_vert; % filter temp with surface mask
    temp_sm.lo = permute(temp_sm.lo, [1 2 4 3]); % bring plev to last dimension
    temp_sm.l = permute(temp_sm.l, [1 2 4 3]); % bring plev to last dimension
    temp_sm.o = permute(temp_sm.o, [1 2 4 3]); % bring plev to last dimension
    pa_sm.lo = pa.*surface_mask; % filter pressure grid with surface mask
    pa_sm.l = pa.*surface_mask.*mask.ocean_vert; % filter pa with surface mask
    pa_sm.o = pa.*surface_mask.*mask.land_vert; % filter pa with surface mask

    % for l = {'lo', 'l', 'o'}; land = l{1};
    %     temp_smz.(land) = squeeze(nanmean(temp_sm.(land), 1)); % take zonal average
    %     temp_smz.(land) = permute(temp_smz.(land), [1 3 2]); % put plev at the last dimension
    %     ps_z.(land) = squeeze(nanmean(ps_n.(land), 1));
    % end

    mask_t.land = nanmean(mask.land, 3);
    mask_t.ocean = nanmean(mask.ocean, 3);
    if any(strcmp(type, {'era5', 'erai'})); f_vec = par.era.fw;
    elseif strcmp(type, 'gcm'); f_vec = par.gcm.fw;
    elseif strcmp(type, 'echam'); f_vec = par.echam.fw; end
    for f = f_vec; fw = f{1};
        for c = fieldnames(rcae_t.lo.ann.(fw))'; crit = c{1};
            for l = {'lo', 'l', 'o'}; land = l{1}; % over land, over ocean, or both
                for t = {'ann', 'djf', 'jja', 'mam', 'son'}; time = t{1};
                    for re = {'rce', 'rae'}; regime = re{1};
                        if strcmp(time, 'ann')
                            temp_n.(land).(time) = squeeze(nanmean(temp_sm.(land), 3));
                            tai_n.(land).(time) = squeeze(nanmean(tai_sm.(land), 3));
                        elseif strcmp(time, 'djf')
                            temp_shift = circshift(temp_sm.(land), 1, 3);
                            tai_shift = circshift(tai_sm.(land), 1, 3);
                            temp_n.(land).(time) = squeeze(nanmean(temp_shift(:,:,1:3,:), 3));
                            tai_n.(land).(time) = squeeze(nanmean(tai_shift(:,:,1:3,:), 3));
                        elseif strcmp(time, 'jja')
                            temp_n.(land).(time) = squeeze(nanmean(temp_sm.(land)(:,:,6:8,:), 3));
                            tai_n.(land).(time) = squeeze(nanmean(tai_sm.(land)(:,:,6:8,:), 3));
                        elseif strcmp(time, 'mam')
                            temp_n.(land).(time) = squeeze(nanmean(temp_sm.(land)(:,:,3:5,:), 3));
                            tai_n.(land).(time) = squeeze(nanmean(tai_sm.(land)(:,:,3:5,:), 3));
                        elseif strcmp(time, 'son')
                            temp_n.(land).(time) = squeeze(nanmean(temp_sm.(land)(:,:,9:11,:), 3));
                            tai_n.(land).(time) = squeeze(nanmean(tai_sm.(land)(:,:,9:11,:), 3));
                        end

                        filt.(land).(time).(fw).(crit).(regime) = nan(size(rcae_t.(land).(time).(fw).(crit))); % create empty arrays to store filtering array
                        if strcmp(regime, 'rce'); filt.(land).(time).(fw).(crit).(regime)(rcae_t.(land).(time).(fw).(crit)==1)=1; % set RCE=1, elsewhere nan
                        elseif strcmp(regime, 'rae'); filt.(land).(time).(fw).(crit).(regime)(rcae_t.(land).(time).(fw).(crit)==-1)=1; % set RAE=1, elsewhere nan
                        end

                        temp_t.(land).(time).(fw).(crit).(regime) = temp_n.(land).(time) .* filt.(land).(time).(fw).(crit).(regime);
                        tai_t.(land).(time).(fw).(crit).(regime) = tai_n.(land).(time) .* filt.(land).(time).(fw).(crit).(regime);

                        nanfilt.(regime) = nan(size(temp_t.(land).(time).(fw).(crit).(regime)));
                        nanfilt.(regime)(~isnan(temp_t.(land).(time).(fw).(crit).(regime))) = 1;

                        % take cosine-weighted meridional average
                        for d = {'all', 'nh', 'sh', 'tp'}; domain = d{1};
                            if strcmp(regime, 'rce')
                                if strcmp(domain, 'all')
                                    cosw = repmat(cosd(lat)', [size(nanfilt.(regime), 1) 1]);
                                    denm = temp_t.(land).(time).(fw).(crit).(regime);
                                    denmi = tai_t.(land).(time).(fw).(crit).(regime);
                                    nume = nanfilt.(regime);
                                elseif strcmp(domain, 'nh')
                                    cosw = repmat(cosd(lat(lat>30))', [size(nanfilt.(regime), 1) 1]);
                                    denm = temp_t.(land).(time).(fw).(crit).(regime)(:,lat>30,:);
                                    denmi = tai_t.(land).(time).(fw).(crit).(regime)(:,lat>30,:);
                                    nume = nanfilt.(regime)(:,lat>30,:);
                                elseif strcmp(domain, 'sh')
                                    cosw = repmat(cosd(lat(lat<-30))', [size(nanfilt.(regime), 1) 1]);
                                    denm = temp_t.(land).(time).(fw).(crit).(regime)(:,lat<-30,:);
                                    denmi = tai_t.(land).(time).(fw).(crit).(regime)(:,lat<-30,:);
                                    nume = nanfilt.(regime)(:,lat<-30,:);
                                elseif strcmp(domain, 'tp')
                                    cosw = repmat(cosd(lat(abs(lat)<10))', [size(nanfilt.(regime), 1) 1]);
                                    denm = temp_t.(land).(time).(fw).(crit).(regime)(:,abs(lat)<10,:);
                                    denmi = tai_t.(land).(time).(fw).(crit).(regime)(:,abs(lat)<10,:);
                                    nume = nanfilt.(regime)(:,abs(lat)<10,:);
                                end
                            elseif strcmp(regime, 'rae')
                                if strcmp(domain, 'all')
                                    cosw = repmat(cosd(lat)', [size(nanfilt.(regime), 1) 1]);
                                    denm = temp_t.(land).(time).(fw).(crit).(regime);
                                    denmi = tai_t.(land).(time).(fw).(crit).(regime);
                                    nume = nanfilt.(regime);
                                elseif strcmp(domain, 'nh')
                                    cosw = repmat(cosd(lat(lat>0))', [size(nanfilt.(regime), 1) 1]);
                                    denm = temp_t.(land).(time).(fw).(crit).(regime)(:,lat>0,:);
                                    denmi = tai_t.(land).(time).(fw).(crit).(regime)(:,lat>0,:);
                                    nume = nanfilt.(regime)(:,lat>0,:);
                                elseif strcmp(domain, 'sh')
                                    cosw = repmat(cosd(lat(lat<0))', [size(nanfilt.(regime), 1) 1]);
                                    denm = temp_t.(land).(time).(fw).(crit).(regime)(:,lat<0,:);
                                    denmi = tai_t.(land).(time).(fw).(crit).(regime)(:,lat<0,:);
                                    nume = nanfilt.(regime)(:,lat<0,:);
                                end
                            end

                            ta.(regime).(domain).(fw).(crit).(land).(time) = nansum(cosw.*denm, 2) ./ nansum(cosw.*nume, 2); % area-weighted meridional average
                            ta.(regime).(domain).(fw).(crit).(land).(time) = squeeze(nanmean(ta.(regime).(domain).(fw).(crit).(land).(time), 1)); % zonal average
                            size(denmi)
                            size(nume)
                            tai.(regime).(domain).(fw).(crit).(land).(time) = nansum(cosw.*denmi, 2) ./ nansum(cosw.*nume, 2); % area-weighted meridional average
                            tai.(regime).(domain).(fw).(crit).(land).(time) = squeeze(nanmean(tai.(regime).(domain).(fw).(crit).(land).(time), 1)); % zonal average
                        end
                    end
                end
            end
        end
    end

    % save filtered data
    printname = [foldername 'ta'];
    if ~exist(foldername, 'dir')
        mkdir(foldername)
    end
    save(printname, 'ta', 'tai');
end % process temperature in RCE/RAE regimes
function proc_ma(type, par)
% calculate moist adiabats
    if strcmp(type, 'era5') | strcmp(type, 'erai')
        foldername = sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/eps_%g_ga_%g/', type, par.lat_interp, par.ep, par.ga);
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s', type);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s', type);
    elseif strcmp(type, 'gcm')
        foldername = sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s/%s/eps_%g_ga_%g/', type, par.model, par.gcm.clim, par.lat_interp, par.ep, par.ga);
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s/%s', type, par.model, par.gcm.clim);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s', type, par.model, par.gcm.clim);
    elseif strcmp(type, 'echam')
        foldername = sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s/eps_%g_ga_%g/', type, par.echam.clim, par.lat_interp, par.ep, par.ga);
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s', type, par.echam.clim);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s', type, par.echam.clim, par.lat_interp);
    end
    load(sprintf('%s/grid.mat', prefix)); % read grid data
    load(sprintf('%s/srfc.mat', prefix)); % read surface variable data
    load(sprintf('%s/%s/masks.mat', prefix_proc, par.lat_interp)); % load land and ocean masks
    load(sprintf('%s/%s/eps_%g_ga_%g/rcae_t.mat', prefix_proc, par.lat_interp, par.ep, par.ga)); % load rcae data

    lat = par.lat_std;
    for fn = fieldnames(srfc)'
        srfc.(fn{1}) = permute(srfc.(fn{1}), [2 1 3]); % bring lat to front
        srfc.(fn{1}) = interp1(grid.dim2.lat, srfc.(fn{1}), lat); % interpolate to standard grid
        srfc.(fn{1}) = permute(srfc.(fn{1}), [2 1 3]); % reorder to original dims
        for l = {'lo', 'l', 'o'}; land = l{1};
            if strcmp(land, 'lo'); srfc_n.(fn{1}).(land) = srfc.(fn{1});
            elseif strcmp(land, 'l'); srfc_n.(fn{1}).(land) = srfc.(fn{1}).*mask.ocean; % filter out ocean
            elseif strcmp(land, 'o'); srfc_n.(fn{1}).(land) = srfc.(fn{1}).*mask.land; % filter out land
            end
        end
    end

    for f = {'mse', 'dse'}; fw = f{1};
        for c = fieldnames(rcae_t.lo.ann.(fw))'; crit = c{1};
            for l = {'lo', 'l', 'o'}; land = l{1};
                for t = {'ann', 'djf', 'jja', 'mam', 'son'}; time = t{1};
                    for re = {'rce', 'rae'}; regime = re{1};
                        for v = fieldnames(srfc)'; vname = v{1};
                            if strcmp(time, 'ann')
                                srfc_t.(land).(time).(vname) = squeeze(nanmean(srfc_n.(vname).(land), 3));
                            elseif strcmp(time, 'djf')
                                srfc_shift = circshift(srfc_n.(vname).(land), 1, 3);
                                srfc_t.(land).(time).(vname) = squeeze(nanmean(srfc_shift(:,:,1:3), 3));
                            elseif strcmp(time, 'jja')
                                srfc_t.(land).(time).(vname) = squeeze(nanmean(srfc_n.(vname).(land)(:,:,6:8), 3));
                            elseif strcmp(time, 'mam')
                                srfc_t.(land).(time).(vname) = squeeze(nanmean(srfc_n.(vname).(land)(:,:,3:5), 3));
                            elseif strcmp(time, 'son')
                                srfc_t.(land).(time).(vname) = squeeze(nanmean(srfc_n.(vname).(land)(:,:,9:11), 3));
                            end

                            filt.(land).(time).(fw).(crit).(regime) = nan(size(rcae_t.(land).(time).(fw).(crit)));
                            if strcmp(regime, 'rce'); filt.(land).(time).(fw).(crit).(regime)(rcae_t.(land).(time).(fw).(crit)==1)=1; % set RCE=1, elsewhere nan
                            elseif strcmp(regime, 'rae'); filt.(land).(time).(fw).(crit).(regime)(rcae_t.(land).(time).(fw).(crit)==-1)=1; % set RAE=1, elsewhere nan
                            end

                            srfc_tf.(land).(time).(fw).(crit).(regime).(vname) = srfc_t.(land).(time).(vname) .* filt.(land).(time).(fw).(crit).(regime);

                            nanfilt.(regime) = nan(size(srfc_tf.(land).(time).(fw).(crit).(regime).(vname)));
                            nanfilt.(regime)(~isnan(srfc_tf.(land).(time).(fw).(crit).(regime).(vname))) = 1;

                            % take cosine-weighted average
                            for d = {'all', 'nh', 'sh', 'tp'}; domain = d{1};
                                if strcmp(regime, 'rce')
                                    if strcmp(domain, 'all')
                                        cosw = repmat(cosd(lat)', [size(nanfilt.(regime), 1) 1]);
                                        denm = srfc_tf.(land).(time).(fw).(crit).(regime).(vname);
                                        nume = nanfilt.(regime);
                                    elseif strcmp(domain, 'nh')
                                        cosw = repmat(cosd(lat(lat>30))', [size(nanfilt.(regime), 1) 1]);
                                        denm = srfc_tf.(land).(time).(fw).(crit).(regime).(vname)(:, lat>30);
                                        nume = nanfilt.(regime)(:, lat>30);
                                    elseif strcmp(domain, 'sh')
                                        cosw = repmat(cosd(lat(lat<-30))', [size(nanfilt.(regime), 1) 1]);
                                        denm = srfc_tf.(land).(time).(fw).(crit).(regime).(vname)(:, lat<-30);
                                        nume = nanfilt.(regime)(:, lat<-30);
                                    elseif strcmp(domain, 'tp')
                                        cosw = repmat(cosd(lat(abs(lat)<10))', [size(nanfilt.(regime), 1) 1]);
                                        denm = srfc_tf.(land).(time).(fw).(crit).(regime).(vname)(:, abs(lat)<10);
                                        nume = nanfilt.(regime)(:, abs(lat)<10);
                                    end
                                elseif strcmp(regime, 'rae')
                                    if strcmp(domain, 'all')
                                        cosw = repmat(cosd(lat)', [size(nanfilt.(regime), 1) 1]);
                                        denm = srfc_tf.(land).(time).(fw).(crit).(regime).(vname);
                                        nume = nanfilt.(regime);
                                    elseif strcmp(domain, 'nh')
                                        cosw = repmat(cosd(lat(lat>0))', [size(nanfilt.(regime), 1) 1]);
                                        denm = srfc_tf.(land).(time).(fw).(crit).(regime).(vname)(:, lat>0);
                                        nume = nanfilt.(regime)(:, lat>0);
                                    elseif strcmp(domain, 'sh')
                                        cosw = repmat(cosd(lat(lat<0))', [size(nanfilt.(regime), 1) 1]);
                                        denm = srfc_tf.(land).(time).(fw).(crit).(regime).(vname)(:, lat<0);
                                        nume = nanfilt.(regime)(:, lat<0);
                                    end
                                end

                                ma.(regime).(domain).(fw).(crit).(land).(time).(vname) = nansum(cosw.*denm, 2) ./ nansum(cosw.*nume, 2); % weighted meridional average
                                ma.(regime).(domain).(fw).(crit).(land).(time).(vname) = squeeze(nanmean(ma.(regime).(domain).(fw).(crit).(land).(time).(vname), 1)); % zonal average

                            end % end domain loop
                        end % end srfc variables loop

                        if strcmp(type, 'era5') | strcmp(type, 'erai')
                            ma.rce.all.(fw).(crit).(land).(time).ta = calc_ma_dew(ma.rce.all.(fw).(crit).(land).(time), grid.dim3.plev, par, type); % compute moist adiabat with dew point temperature
                            ma.rce.tp.(fw).(crit).(land).(time).ta = calc_ma_dew(ma.rce.tp.(fw).(crit).(land).(time), grid.dim3.plev, par, type); % compute moist adiabat with dew point temperature
                            ma.rce.nh.(fw).(crit).(land).(time).ta = calc_ma_dew(ma.rce.nh.(fw).(crit).(land).(time), grid.dim3.plev, par, type); % compute moist adiabat with dew point temperature
                            ma.rce.sh.(fw).(crit).(land).(time).ta = calc_ma_dew(ma.rce.sh.(fw).(crit).(land).(time), grid.dim3.plev, par, type); % compute moist adiabat with dew point temperature
                        elseif strcmp(type, 'gcm')
                            ma.rce.all.(fw).(crit).(land).(time).ta = calc_ma_hurs(ma.rce.all.(fw).(crit).(land).(time), grid.dim3.plev, par); % compute moist adiabat with RH
                            ma.rce.tp.(fw).(crit).(land).(time).ta = calc_ma_hurs(ma.rce.tp.(fw).(crit).(land).(time), grid.dim3.plev, par); % compute moist adiabat with RH
                            ma.rce.nh.(fw).(crit).(land).(time).ta = calc_ma_hurs(ma.rce.nh.(fw).(crit).(land).(time), grid.dim3.plev, par); % compute moist adiabat with RH
                            ma.rce.sh.(fw).(crit).(land).(time).ta = calc_ma_hurs(ma.rce.sh.(fw).(crit).(land).(time), grid.dim3.plev, par); % compute moist adiabat with RH
                        end

                    end % end RCE/RAE regime loop
                end % end time average loop
            end % end land option loop
        end % end RCAE definition loop
    end % end MSE/DSE framework loop

    % save data into mat file
    printname = [foldername 'ma.mat'];
    if ~exist(foldername, 'dir')
        mkdir(foldername)
    end
    save(printname, 'ma');

end % process moist adiabat in RCE/RAE regimes

function comp_rad(par)
    root = '/project2/tas1/miyawaki/projects/002/data';
    don = load(sprintf('%s/raw/don/radiation_dynamics_climatology.mat', root));
    grid.don.lat = don.lat; grid.don.lon = don.lon;
    don.ttr = don.OLR; don = rmfield(don, 'OLR');
    don.tsr = don.RSDT - don.RSUT; don.net = don.tsr - don.ttr;
    don.ssr = -(don.total_surface_flux - don.surface_turbulent_plus_LW);
    don.swabs = don.tsr - don.ssr;
    tmp = load(sprintf('%s/read/ceres/rad.mat', root)); ceres.rad = tmp.rad; clear tmp;
    ceres.rad.tsr = ceres.rad.tsdr - ceres.rad.tsur; ceres.rad.net = ceres.rad.tsr - ceres.rad.ttr;
    ceres.rad.swabs = ceres.rad.tsr - ceres.rad.ssr;
    ceres.rad.str = -ceres.rad.str;
    ceres.rad.ra = ceres.rad.tsr - ceres.rad.ssr + ceres.rad.str - ceres.rad.ttr;
    tmp = load(sprintf('%s/read/ceres/grid.mat', root)); grid.ceres = tmp.grid; clear tmp;
    tmp = load(sprintf('%s/read/erai/rad.mat', root)); erai.rad = tmp.rad; clear tmp;
    tmp = load(sprintf('%s/read/erai/stf.mat', root)); erai.stf = tmp.stf; clear tmp;
    tmp = load(sprintf('%s/read/erai/grid.mat', root)); grid.erai = tmp.grid; clear tmp;
    erai.rad.ttr = -erai.rad.ttr; erai.rad.str = -erai.rad.str;
    erai.rad.net = erai.rad.tsr - erai.rad.ttr;
    erai.rad.swabs = erai.rad.tsr - erai.rad.ssr;
    erai.rad.ra = erai.rad.tsr - erai.rad.ssr + erai.rad.str - erai.rad.ttr;
    erai.rad.surface_turbulent_plus_LW = -erai.stf.sshf - erai.stf.slhf + erai.rad.str;
    tmp = load(sprintf('%s/read/era5/rad_2000_2012.mat', root)); era5.rad = tmp.rad; clear tmp;
    tmp = load(sprintf('%s/read/era5/stf_2000_2012.mat', root)); era5.stf = tmp.stf; clear tmp;
    tmp = load(sprintf('%s/read/era5/grid.mat', root)); grid.era5 = tmp.grid; clear tmp;
    era5.rad.ttr = -era5.rad.ttr; era5.rad.str = -era5.rad.str;
    era5.rad.net = era5.rad.tsr - era5.rad.ttr;
    era5.rad.swabs = era5.rad.tsr - era5.rad.ssr;
    era5.rad.ra = era5.rad.tsr - era5.rad.ssr + era5.rad.str - era5.rad.ttr;
    era5.rad.surface_turbulent_plus_LW = - era5.stf.sshf - era5.stf.slhf + era5.rad.str;

    lat = par.lat_std;

    % take zonal and time averages
    for fn = fieldnames(don)'; fname = fn{1};
        if ~any(strcmp(fname, {'lat', 'lon'}))
            don_z.(fname) = squeeze(nanmean(don.(fname), 3));
            don_z.(fname) = permute(don_z.(fname), [2 1]);
            don_z.(fname) = interp1(don.lat, don_z.(fname), lat);
            don_z.(fname) = permute(don_z.(fname), [2 1]);
            don_t.(fname) = interp1(don.lat, squeeze(nanmean(don.(fname), 1)), lat);
            don_zt.(fname) = squeeze(nanmean(don_z.(fname), 1));
            don_zt.(fname) = permute(don_zt.(fname), [2 1]);
        end
    end

    for d = {'rad'}; dtype = d{1};
        for fn = fieldnames(ceres.(dtype))'; fname = fn{1};
            ceres_z.(dtype).(fname) = interp1(grid.ceres.dim2.lat, squeeze(nanmean(ceres.(dtype).(fname), 1)), lat);
            ceres_t.(dtype).(fname) = squeeze(nanmean(ceres.(dtype).(fname), 3));
            ceres_t.(dtype).(fname) = permute(ceres_t.(dtype).(fname), [2 1]);
            ceres_t.(dtype).(fname) = interp1(grid.ceres.dim2.lat, ceres_t.(dtype).(fname), lat);
            ceres_t.(dtype).(fname) = permute(ceres_t.(dtype).(fname), [2 1]);
            ceres_zt.(dtype).(fname) = squeeze(nanmean(ceres_z.(dtype).(fname), 2));
        end
        for fn = fieldnames(erai.(dtype))'; fname = fn{1};
            erai_z.(dtype).(fname) = interp1(grid.erai.dim2.lat, squeeze(nanmean(erai.(dtype).(fname), 1)), lat);
            erai_t.(dtype).(fname) = squeeze(nanmean(erai.(dtype).(fname), 3));
            erai_t.(dtype).(fname) = permute(erai_t.(dtype).(fname), [2 1]);
            erai_t.(dtype).(fname) = interp1(grid.erai.dim2.lat, erai_t.(dtype).(fname), lat);
            erai_t.(dtype).(fname) = permute(erai_t.(dtype).(fname), [2 1]);
            erai_zt.(dtype).(fname) = squeeze(nanmean(erai_z.(dtype).(fname), 2));
        end
        for fn = fieldnames(era5.(dtype))'; fname = fn{1};
            era5_z.(dtype).(fname) = interp1(grid.era5.dim2.lat, squeeze(nanmean(era5.(dtype).(fname), 1)), lat);
            era5_t.(dtype).(fname) = squeeze(nanmean(era5.(dtype).(fname), 3));
            era5_t.(dtype).(fname) = permute(era5_t.(dtype).(fname), [2 1]);
            era5_t.(dtype).(fname) = interp1(grid.era5.dim2.lat, era5_t.(dtype).(fname), lat);
            era5_t.(dtype).(fname) = permute(era5_t.(dtype).(fname), [2 1]);
            era5_zt.(dtype).(fname) = squeeze(nanmean(era5_z.(dtype).(fname), 2));
        end
    end

    save(sprintf('%s/proc/comp/comp_zt', root), 'don_zt', 'ceres_zt', 'erai_zt', 'era5_zt', 'grid', 'lat')
    save(sprintf('%s/proc/comp/comp_z', root), 'don_z', 'ceres_z', 'erai_z', 'era5_z', 'grid', 'lat')
    save(sprintf('%s/proc/comp/comp_t', root), 'don_t', 'ceres_t', 'erai_t', 'era5_t', 'grid', 'lat')
end % radiative flux comparisons (ERA, CERES, DB13)
function ceres_rad(par)
    root = '/project2/tas1/miyawaki/projects/002/data';
    tmp = load(sprintf('%s/read/ceres/rad_2001_2009.mat', root)); ceres.rad = tmp.rad; clear tmp;
    ceres.rad.tsr = ceres.rad.tsdr - ceres.rad.tsur; ceres.rad.net = ceres.rad.tsr - ceres.rad.ttr;
    ceres.rad.swabs = ceres.rad.tsr - ceres.rad.ssr;
    ceres.rad.str = -ceres.rad.str;
    ceres.rad.ra = ceres.rad.tsr - ceres.rad.ssr + ceres.rad.str - ceres.rad.ttr;
    tmp = load(sprintf('%s/read/ceres/grid.mat', root)); grid.ceres = tmp.grid; clear tmp;

    lat = par.lat_std;

    % take zonal and time averages
    for d = {'rad'}; dtype = d{1};
        for fn = fieldnames(ceres.(dtype))'; fname = fn{1};
            ceres_z.(dtype).(fname) = interp1(grid.ceres.dim2.lat, squeeze(nanmean(ceres.(dtype).(fname), 1)), lat);
            ceres_t.(dtype).(fname) = squeeze(nanmean(ceres.(dtype).(fname), 3));
            ceres_t.(dtype).(fname) = permute(ceres_t.(dtype).(fname), [2 1]);
            ceres_t.(dtype).(fname) = interp1(grid.ceres.dim2.lat, ceres_t.(dtype).(fname), lat);
            ceres_t.(dtype).(fname) = permute(ceres_t.(dtype).(fname), [2 1]);
            ceres_zt.(dtype).(fname) = squeeze(nanmean(ceres_z.(dtype).(fname), 2));
        end
    end

    save(sprintf('%s/proc/comp/ceres_2001_2009.mat', root), 'ceres_z', 'ceres_t', 'ceres_zt', 'grid', 'lat')
end

function choose_disp(par)
    % disp_global_rad(par)
    disp_global_stf(par)
end % display globally-averaged radiation fluxes
function disp_global_rad(par)
    root = '/project2/tas1/miyawaki/projects/002/data';
    load(sprintf('%s/proc/comp/comp_zt', root));

    % output global averages
    for fn = {'tsr', 'ssr', 'ttr', 'str'}; fname = fn{1};
        disp( sprintf('CERES %s is %g Wm^-2', upper(fname), nansum(cosd(lat).*ceres_zt.rad.(fname))/nansum(cosd(lat))) );
        disp( sprintf('ERA-I %s is %g Wm^-2', upper(fname), nansum(cosd(lat).*erai_zt.rad.(fname))/nansum(cosd(lat))) );
        disp( sprintf('ERA5 %s is %g Wm^-2', upper(fname), nansum(cosd(lat).*era5_zt.rad.(fname))/nansum(cosd(lat))) );
    end
end
function disp_global_stf(par)
    root = '/project2/tas1/miyawaki/projects/002/data';

    for t = {'erai', 'era5'}; type = t{1};
        load(sprintf('%s/proc/%s/%s/flux_zt.mat', root, type, par.lat_interp));

        for l = {'lo', 'l', 'o'}; land = l{1}; disp(sprintf('---------- %s ----------',land));
            % output global averages
            for fn = {'slhf', 'sshf'}; fname = fn{1};
                disp( sprintf('%s %s is %g Wm^-2', upper(type), upper(fname), nansum(cosd(lat).*flux_zt.(land).ann.(fname)')/nansum(cosd(lat))) );
            end
        end
    end
end

% helper functions
function rcae = def_rcae(type, flux, vh, par)
    if any(strcmp(type, {'era5', 'erai'})); f_vec = par.era.fw;
    elseif strcmp(type, 'gcm'); f_vec = par.gcm.fw;
    elseif strcmp(type, 'echam'); f_vec = par.echam.fw; end
    for f = f_vec; fw = f{1};
        % identify locations of RCE using threshold epsilon (ep)
        rcae.(fw).def = zeros(size(flux.r1.(fw)));
        rcae.(fw).def(abs(flux.r1.(fw)) < par.ep) = 1;
        % identify locations of RAE using threshold gamma (ga)
        if strcmp(fw, 'dse'); rcae.(fw).def(flux.r1.mse > par.ga) = -1; % always use MSE for RAE
        else rcae.(fw).def(flux.r1.(fw) > par.ga) = -1; end;

        % identify locations of RCE following Jakob et al. (2019) (abs(TEDIV) < 50 W/m^2)
        rcae.(fw).jak = zeros(size(flux.r1.(fw)));
        rcae.(fw).jak(abs(flux.res.(fw)) < 50) = 1;
        % identify locations of RAE using threshold gamma (ga)
        if strcmp(fw, 'dse'); rcae.(fw).jak(flux.r1.mse > par.ga) = -1; % always use MSE for RAE
        else rcae.(fw).jak(flux.r1.(fw) > par.ga) = -1; end;

        % identify locations of rce following jakob et al. (2019) (abs(tediv) < 50 w/m^2)
        rcae.(fw).jak30 = zeros(size(flux.r1.(fw)));
        rcae.(fw).jak30(abs(flux.res.(fw)) < 30) = 1;
        % identify locations of rae using threshold gamma (ga)
        if strcmp(fw, 'dse'); rcae.(fw).jak30(flux.r1.mse > par.ga) = -1; % always use mse for rae
        else rcae.(fw).jak30(flux.r1.(fw) > par.ga) = -1; end;

        % identify locations of rce following jakob et al. (2019) (abs(tediv) < 50 w/m^2)
        rcae.(fw).jak10 = zeros(size(flux.r1.(fw)));
        rcae.(fw).jak10(abs(flux.res.(fw)) < 10) = 1;
        % identify locations of rae using threshold gamma (ga)
        if strcmp(fw, 'dse'); rcae.(fw).jak10(flux.r1.mse > par.ga) = -1; % always use mse for rae
        else rcae.(fw).jak10(flux.r1.(fw) > par.ga) = -1; end;

        % add additional criteria for RCE that P-E>0
        rcae.(fw).pe = zeros(size(flux.r1.(fw)));
        if strcmp(type, 'era5') | strcmp(type, 'erai')
            rcae.(fw).pe(abs(flux.r1.(fw)) < par.ep & (flux.cp+flux.lsp+flux.e>0)) = 1; % note that evap is defined negative into atmosphere
        elseif strcmp(type, 'gcm')
            rcae.(fw).pe(abs(flux.r1.(fw)) < par.ep & (flux.pr-flux.evspsbl>0)) = 1;
        elseif strcmp(type, 'echam')
            rcae.(fw).pe(abs(flux.r1.(fw)) < par.ep & (flux.aprc+flux.aprl-flux.evap>0)) = 1;
        end
        if strcmp(fw, 'dse'); rcae.(fw).pe(flux.r1.mse > par.ga) = -1;
        else rcae.(fw).pe(flux.r1.(fw) > par.ga) = -1; end;

        % add additional criteria for RCE that (P_ls - E)<<1
        rcae.(fw).cp = zeros(size(flux.r1.(fw)));
        if strcmp(type, 'era5') | strcmp(type, 'erai')
            rcae.(fw).cp(abs(flux.r1.(fw)) < par.ep & abs(flux.lsp./flux.cp<par.ep_cp)) = 1; % note that evap is defined negative into atmosphere
        elseif strcmp(type, 'gcm')
            rcae.(fw).cp(abs(flux.r1.(fw)) < par.ep & abs((flux.pr-flux.prc)./flux.prc<par.ep_cp)) = 1; % note that evap is defined negative into atmosphere
        elseif strcmp(type, 'echam')
            rcae.(fw).cp(abs(flux.r1.(fw)) < par.ep & abs(flux.aprl./flux.aprc<par.ep_cp)) = 1; % note that evap is defined negative into atmosphere
        end
        if strcmp(fw, 'dse') rcae.(fw).cp(flux.r1.mse > par.ga) = -1;
        else rcae.(fw).cp(flux.r1.(fw) > par.ga) = -1; end;

        % % add additional criteria for RCE that w500<0 (ascent)
        % rcae.(fw).w500 = zeros(size(flux.r1.(fw)));
        % if strcmp(type, 'era5') | strcmp(type, 'erai')
        %     rcae.(fw).w500(abs(flux.r1.(fw)) < par.ep & flux.w500 < 0) = 1; % note that evap is defined negative into atmosphere
        % elseif strcmp(type, 'gcm')
        %     rcae.(fw).w500(abs(flux.r1.(fw)) < par.ep & flux.w500 < 0) = 1; % note that evap is defined negative into atmosphere
        % end
        % if strcmp(fw, 'dse'); rcae.(fw).w500(flux.r1.mse > par.ga) = -1;
        % else rcae.(fw).w500(flux.r1.(fw) > par.ga) = -1; end;

        % % add additional criteria for RCE that vh<max(vh)/2 (weak meridional velocity)
        % rcae.(fw).vas2 = zeros(size(flux.r1.(fw)));
        % if strcmp(type, 'era5') | strcmp(type, 'erai')
        %     rcae.(fw).vas2(abs(flux.r1.(fw)) < par.ep & abs(flux.vas) < nanmax(nanmax(abs(flux.vas)))/2) = 1; % note that evap is defined negative into atmosphere
        % elseif strcmp(type, 'gcm')
        %     rcae.(fw).vas2(abs(flux.r1.(fw)) < par.ep & flux.vas < nanmax(nanmax(abs(flux.vas)))/2) = 1; % note that evap is defined negative into atmosphere
        % end
        % if strcmp(fw, 'dse'); rcae.(fw).vas2(flux.r1.mse > par.ga) = -1;
        % else rcae.(fw).vas2(flux.r1.(fw) > par.ga) = -1; end;

        % % add additional criteria for RCE that vh<max(vh)/4 (weak meridional velocity)
        % rcae.(fw).vas4 = zeros(size(flux.r1.(fw)));
        % if strcmp(type, 'era5') | strcmp(type, 'erai')
        %     rcae.(fw).vas4(abs(flux.r1.(fw)) < par.ep & abs(flux.vas) < nanmax(nanmax(abs(flux.vas)))/4) = 1; % note that evap is defined negative into atmosphere
        % elseif strcmp(type, 'gcm')
        %     rcae.(fw).vas4(abs(flux.r1.(fw)) < par.ep & flux.vas < nanmax(nanmax(abs(flux.vas)))/4) = 1; % note that evap is defined negative into atmosphere
        % end
        % if strcmp(fw, 'dse'); rcae.(fw).vas4(flux.r1.mse > par.ga) = -1;
        % else rcae.(fw).vas4(flux.r1.(fw) > par.ga) = -1; end;

        if size(flux.r1.(fw), 1)==length(par.lat_std) & size(flux.r1.(fw), 2)==12 % this criteria only works for zonally and time averaged data because vh is only a function of lat
            % add additional criteria for RCE that horizontal transport is weak
            rcae.(fw).vh2 = zeros(size(flux.r1.(fw)));
            % if strcmp(type, 'era5') | strcmp(type, 'erai')
                rcae.(fw).vh2(abs(flux.r1.(fw)) < par.ep & abs(vh.(fw)) < nanmax(nanmax(abs(vh.(fw))))/2 ) = 1;
            % elseif strcmp(type, 'gcm')
            %     rcae.(fw).vh2(abs(flux.r1.(fw)) < par.ep & abs(vh.(fw)) < nanmax(nanmax(abs(vh.(fw))))/2 ) = 1;
            % end
            if strcmp(fw, 'dse'); rcae.(fw).vh2(flux.r1.mse > par.ga) = -1;
            else; rcae.(fw).vh2(flux.r1.(fw) > par.ga) = -1; end;

            rcae.(fw).vh3 = zeros(size(flux.r1.(fw)));
            % if strcmp(type, 'era5') | strcmp(type, 'erai')
                rcae.(fw).vh3(abs(flux.r1.(fw)) < par.ep & abs(vh.(fw)) < nanmax(nanmax(abs(vh.(fw))))/3 ) = 1;
            % elseif strcmp(type, 'gcm')
            %     rcae.(fw).vh3(abs(flux.r1.(fw)) < par.ep & abs(vh.(fw)) < nanmax(nanmax(abs(vh.(fw))))/3 ) = 1;
            % end
            if strcmp(fw, 'dse') rcae.(fw).vh3(flux.r1.mse > par.ga) = -1;
            else rcae.(fw).vh3(flux.r1.(fw) > par.ga) = -1; end;

            rcae.(fw).vh4 = zeros(size(flux.r1.(fw)));
            % if strcmp(type, 'era5') | strcmp(type, 'erai')
                rcae.(fw).vh4(abs(flux.r1.(fw)) < par.ep & abs(vh.(fw)) < nanmax(nanmax(abs(vh.(fw))))/4 ) = 1;
            % elseif strcmp(type, 'gcm')
            %     rcae.(fw).vh4(abs(flux.r1.(fw)) < par.ep & abs(vh.(fw)) < nanmax(nanmax(abs(vh.(fw))))/4 ) = 1;
            % end
            if strcmp(fw, 'dse'); rcae.(fw).vh4(flux.r1.mse > par.ga) = -1;
            else; rcae.(fw).vh4(flux.r1.(fw) > par.ga) = -1; end;
        end

    end
end % define RCE and RAE
function rcae = def_rcae_recomp_r1(type, flux, vh, par) % recalculate R1 at this step
    if any(strcmp(type, {'era5', 'erai'})); f_vec = par.era.fw;
    elseif strcmp(type, 'gcm'); f_vec = par.gcm.fw;
    elseif strcmp(type, 'echam'); f_vec = par.echam.fw; end
    for f = f_vec; fw = f{1};
        % calculate R1 again to test importance of order of operations
        if strcmp(fw, 'mse2'); r1 = flux.res.(fw)./flux.lw;
        else r1 = flux.res.(fw)./flux.ra.(fw); end

        % identify locations of RCE using threshold epsilon (ep)
        rcae.(fw).def = zeros(size(r1));
        rcae.(fw).def(abs(r1) < par.ep) = 1;
        % identify locations of RAE using threshold gamma (ga)
        if strcmp(fw, 'dse'); rcae.(fw).def(flux.r1.mse > par.ga) = -1; % always use MSE for RAE
        else rcae.(fw).def(r1 > par.ga) = -1; end;

        % identify locations of RCE following Jakob et al. (2019) (abs(TEDIV) < 50 W/m^2)
        rcae.(fw).jak = zeros(size(r1));
        rcae.(fw).jak(abs(flux.res.(fw)) < 50) = 1;
        % identify locations of RAE using threshold gamma (ga)
        if strcmp(fw, 'dse'); rcae.(fw).jak(flux.r1.mse > par.ga) = -1; % always use MSE for RAE
        else rcae.(fw).jak(r1 > par.ga) = -1; end;

        % identify locations of rce following jakob et al. (2019) (abs(tediv) < 50 w/m^2)
        rcae.(fw).jak30 = zeros(size(r1));
        rcae.(fw).jak30(abs(flux.res.(fw)) < 30) = 1;
        % identify locations of rae using threshold gamma (ga)
        if strcmp(fw, 'dse'); rcae.(fw).jak30(flux.r1.mse > par.ga) = -1; % always use mse for rae
        else rcae.(fw).jak30(r1 > par.ga) = -1; end;

        % identify locations of rce following jakob et al. (2019) (abs(tediv) < 50 w/m^2)
        rcae.(fw).jak10 = zeros(size(r1));
        rcae.(fw).jak10(abs(flux.res.(fw)) < 10) = 1;
        % identify locations of rae using threshold gamma (ga)
        if strcmp(fw, 'dse'); rcae.(fw).jak10(flux.r1.mse > par.ga) = -1; % always use mse for rae
        else rcae.(fw).jak10(r1 > par.ga) = -1; end;

        % add additional criteria for RCE that P-E>0
        rcae.(fw).pe = zeros(size(r1));
        if strcmp(type, 'era5') | strcmp(type, 'erai')
            rcae.(fw).pe(abs(r1) < par.ep & (flux.cp+flux.lsp+flux.e>0)) = 1; % note that evap is defined negative into atmosphere
        elseif strcmp(type, 'gcm')
            rcae.(fw).pe(abs(r1) < par.ep & (flux.pr-flux.evspsbl>0)) = 1;
        elseif strcmp(type, 'echam')
            rcae.(fw).pe(abs(flux.r1.(fw)) < par.ep & (flux.aprc+flux.aprl-flux.evap>0)) = 1;
        end
        if strcmp(fw, 'dse'); rcae.(fw).pe(flux.r1.mse > par.ga) = -1;
        else rcae.(fw).pe(r1 > par.ga) = -1; end;

        % add additional criteria for RCE that (P_ls - E)<<1
        rcae.(fw).cp = zeros(size(r1));
        if strcmp(type, 'era5') | strcmp(type, 'erai')
            rcae.(fw).cp(abs(r1) < par.ep & abs(flux.lsp./flux.cp<par.ep_cp)) = 1; % note that evap is defined negative into atmosphere
        elseif strcmp(type, 'gcm')
            rcae.(fw).cp(abs(r1) < par.ep & abs((flux.pr-flux.prc)./flux.prc<par.ep_cp)) = 1; % note that evap is defined negative into atmosphere
        elseif strcmp(type, 'echam')
            rcae.(fw).cp(abs(flux.r1.(fw)) < par.ep & abs(flux.aprl./flux.aprc<par.ep_cp)) = 1; % note that evap is defined negative into atmosphere
        end
        if strcmp(fw, 'dse') rcae.(fw).cp(flux.r1.mse > par.ga) = -1;
        else rcae.(fw).cp(r1 > par.ga) = -1; end;

        % % add additional criteria for RCE that w500<0 (ascent)
        % rcae.(fw).w500 = zeros(size(r1));
        % if strcmp(type, 'era5') | strcmp(type, 'erai')
        %     rcae.(fw).w500(abs(r1) < par.ep & flux.w500 < 0) = 1; % note that evap is defined negative into atmosphere
        % elseif strcmp(type, 'gcm')
        %     rcae.(fw).w500(abs(r1) < par.ep & flux.w500 < 0) = 1; % note that evap is defined negative into atmosphere
        % end
        % if strcmp(fw, 'dse'); rcae.(fw).w500(flux.r1.mse > par.ga) = -1;
        % else rcae.(fw).w500(r1 > par.ga) = -1; end;

        % % add additional criteria for RCE that vh<max(vh)/2 (weak meridional velocity)
        % rcae.(fw).vas2 = zeros(size(r1));
        % if strcmp(type, 'era5') | strcmp(type, 'erai')
        %     rcae.(fw).vas2(abs(r1) < par.ep & abs(flux.vas) < nanmax(nanmax(abs(flux.vas)))/2) = 1; % note that evap is defined negative into atmosphere
        % elseif strcmp(type, 'gcm')
        %     rcae.(fw).vas2(abs(r1) < par.ep & flux.vas < nanmax(nanmax(abs(flux.vas)))/2) = 1; % note that evap is defined negative into atmosphere
        % end
        % if strcmp(fw, 'dse'); rcae.(fw).vas2(flux.r1.mse > par.ga) = -1;
        % else rcae.(fw).vas2(r1 > par.ga) = -1; end;

        % % add additional criteria for RCE that vh<max(vh)/4 (weak meridional velocity)
        % rcae.(fw).vas4 = zeros(size(r1));
        % if strcmp(type, 'era5') | strcmp(type, 'erai')
        %     rcae.(fw).vas4(abs(r1) < par.ep & abs(flux.vas) < nanmax(nanmax(abs(flux.vas)))/4) = 1; % note that evap is defined negative into atmosphere
        % elseif strcmp(type, 'gcm')
        %     rcae.(fw).vas4(abs(r1) < par.ep & flux.vas < nanmax(nanmax(abs(flux.vas)))/4) = 1; % note that evap is defined negative into atmosphere
        % end
        % if strcmp(fw, 'dse'); rcae.(fw).vas4(flux.r1.mse > par.ga) = -1;
        % else rcae.(fw).vas4(r1 > par.ga) = -1; end;

        if size(r1, 1)==length(par.lat_std) & size(r1, 2)==12 % this criteria only works for zonally and time averaged data because vh is only a function of lat
            % add additional criteria for RCE that horizontal transport is weak
            rcae.(fw).vh2 = zeros(size(r1));
            % if strcmp(type, 'era5') | strcmp(type, 'erai')
                rcae.(fw).vh2(abs(r1) < par.ep & abs(vh.(fw)) < nanmax(nanmax(abs(vh.(fw))))/2 ) = 1;
            % elseif strcmp(type, 'gcm')
            %     rcae.(fw).vh2(abs(r1) < par.ep & abs(vh.(fw)) < nanmax(nanmax(abs(vh.(fw))))/2 ) = 1;
            % end
            if strcmp(fw, 'dse'); rcae.(fw).vh2(flux.r1.mse > par.ga) = -1;
            else; rcae.(fw).vh2(r1 > par.ga) = -1; end;

            rcae.(fw).vh3 = zeros(size(r1));
            % if strcmp(type, 'era5') | strcmp(type, 'erai')
                rcae.(fw).vh3(abs(r1) < par.ep & abs(vh.(fw)) < nanmax(nanmax(abs(vh.(fw))))/3 ) = 1;
            % elseif strcmp(type, 'gcm')
            %     rcae.(fw).vh3(abs(r1) < par.ep & abs(vh.(fw)) < nanmax(nanmax(abs(vh.(fw))))/3 ) = 1;
            % end
            if strcmp(fw, 'dse') rcae.(fw).vh3(flux.r1.mse > par.ga) = -1;
            else rcae.(fw).vh3(r1 > par.ga) = -1; end;

            rcae.(fw).vh4 = zeros(size(r1));
            % if strcmp(type, 'era5') | strcmp(type, 'erai')
                rcae.(fw).vh4(abs(r1) < par.ep & abs(vh.(fw)) < nanmax(nanmax(abs(vh.(fw))))/4 ) = 1;
            % elseif strcmp(type, 'gcm')
            %     rcae.(fw).vh4(abs(r1) < par.ep & abs(vh.(fw)) < nanmax(nanmax(abs(vh.(fw))))/4 ) = 1;
            % end
            if strcmp(fw, 'dse'); rcae.(fw).vh4(flux.r1.mse > par.ga) = -1;
            else; rcae.(fw).vh4(r1 > par.ga) = -1; end;
        end

    end
end % define RCE and RAE
function rcae = def_rcae_alt(type, flux, vh, par)
    if any(strcmp(type, {'era5', 'erai'})); f_vec = par.era.fw;
    elseif strcmp(type, 'gcm'); f_vec = par.gcm.fw;
    elseif strcmp(type, 'echam'); f_vec = par.echam.fw; end
    for f = f_vec; fw = f{1};
        % identify locations of RCE using threshold epsilon (ep)
        rcae.(fw).def = zeros(size(flux.r1.(fw)));
        rcae.(fw).def(flux.r1.(fw) < par.ep) = 1;
        % identify locations of RAE using threshold gamma (ga)
        if strcmp(fw, 'dse'); rcae.(fw).def(flux.r1.mse > par.ga) = -1; % always use MSE for RAE
        else rcae.(fw).def(flux.r1.(fw) > par.ga) = -1; end;

        % identify locations of RCE following Jakob et al. (2019) (abs(TEDIV) < 50 W/m^2)
        rcae.(fw).jak = zeros(size(flux.r1.(fw)));
        rcae.(fw).jak(flux.res.(fw) < 50) = 1;
        % identify locations of RAE using threshold gamma (ga)
        if strcmp(fw, 'dse'); rcae.(fw).jak(flux.r1.mse > par.ga) = -1; % always use MSE for RAE
        else rcae.(fw).jak(flux.r1.(fw) > par.ga) = -1; end;

        % identify locations of rce following jakob et al. (2019) (abs(tediv) < 50 w/m^2)
        rcae.(fw).jak30 = zeros(size(flux.r1.(fw)));
        rcae.(fw).jak30(flux.res.(fw) < 30) = 1;
        % identify locations of rae using threshold gamma (ga)
        if strcmp(fw, 'dse'); rcae.(fw).jak30(flux.r1.mse > par.ga) = -1; % always use mse for rae
        else rcae.(fw).jak30(flux.r1.(fw) > par.ga) = -1; end;

        % identify locations of rce following jakob et al. (2019) (abs(tediv) < 50 w/m^2)
        rcae.(fw).jak10 = zeros(size(flux.r1.(fw)));
        rcae.(fw).jak10(flux.res.(fw) < 10) = 1;
        % identify locations of rae using threshold gamma (ga)
        if strcmp(fw, 'dse'); rcae.(fw).jak10(flux.r1.mse > par.ga) = -1; % always use mse for rae
        else rcae.(fw).jak10(flux.r1.(fw) > par.ga) = -1; end;

        % add additional criteria for RCE that P-E>0
        rcae.(fw).pe = zeros(size(flux.r1.(fw)));
        if strcmp(type, 'era5') | strcmp(type, 'erai')
            rcae.(fw).pe(flux.r1.(fw) < par.ep & (flux.cp+flux.lsp+flux.e>0)) = 1; % note that evap is defined negative into atmosphere
        elseif strcmp(type, 'gcm')
            rcae.(fw).pe(flux.r1.(fw) < par.ep & (flux.pr-flux.evspsbl>0)) = 1;
        elseif strcmp(type, 'echam')
            rcae.(fw).pe(flux.r1.(fw) < par.ep & (flux.aprc+flux.aprl-flux.evap>0)) = 1;
        end
        if strcmp(fw, 'dse'); rcae.(fw).pe(flux.r1.mse > par.ga) = -1;
        else rcae.(fw).pe(flux.r1.(fw) > par.ga) = -1; end;

        % add additional criteria for RCE that (P_ls - E)<<1
        rcae.(fw).cp = zeros(size(flux.r1.(fw)));
        if strcmp(type, 'era5') | strcmp(type, 'erai')
            rcae.(fw).cp(flux.r1.(fw) < par.ep & abs(flux.lsp./flux.cp<par.ep_cp)) = 1; % note that evap is defined negative into atmosphere
        elseif strcmp(type, 'gcm')
            rcae.(fw).cp(flux.r1.(fw) < par.ep & abs((flux.pr-flux.prc)./flux.prc<par.ep_cp)) = 1; % note that evap is defined negative into atmosphere
        elseif strcmp(type, 'echam')
            rcae.(fw).cp(flux.r1.(fw) < par.ep & abs(flux.aprl./flux.aprc<par.ep_cp)) = 1; % note that evap is defined negative into atmosphere
        end
        if strcmp(fw, 'dse') rcae.(fw).cp(flux.r1.mse > par.ga) = -1;
        else rcae.(fw).cp(flux.r1.(fw) > par.ga) = -1; end;

        % % add additional criteria for RCE that w500<0 (ascent)
        % rcae.(fw).w500 = zeros(size(flux.r1.(fw)));
        % if strcmp(type, 'era5') | strcmp(type, 'erai')
        %     rcae.(fw).w500(flux.r1.(fw) < par.ep & flux.w500 < 0) = 1; % note that evap is defined negative into atmosphere
        % elseif strcmp(type, 'gcm')
        %     rcae.(fw).w500(flux.r1.(fw) < par.ep & flux.w500 < 0) = 1; % note that evap is defined negative into atmosphere
        % end
        % if strcmp(fw, 'dse'); rcae.(fw).w500(flux.r1.mse > par.ga) = -1;
        % else rcae.(fw).w500(flux.r1.(fw) > par.ga) = -1; end;

        % % add additional criteria for RCE that vh<max(vh)/2 (weak meridional velocity)
        % rcae.(fw).vas2 = zeros(size(flux.r1.(fw)));
        % if strcmp(type, 'era5') | strcmp(type, 'erai')
        %     rcae.(fw).vas2(flux.r1.(fw) < par.ep & abs(flux.vas) < nanmax(nanmax(abs(flux.vas)))/2) = 1; % note that evap is defined negative into atmosphere
        % elseif strcmp(type, 'gcm')
        %     rcae.(fw).vas2(flux.r1.(fw) < par.ep & flux.vas < nanmax(nanmax(abs(flux.vas)))/2) = 1; % note that evap is defined negative into atmosphere
        % end
        % if strcmp(fw, 'dse'); rcae.(fw).vas2(flux.r1.mse > par.ga) = -1;
        % else rcae.(fw).vas2(flux.r1.(fw) > par.ga) = -1; end;

        % % add additional criteria for RCE that vh<max(vh)/4 (weak meridional velocity)
        % rcae.(fw).vas4 = zeros(size(flux.r1.(fw)));
        % if strcmp(type, 'era5') | strcmp(type, 'erai')
        %     rcae.(fw).vas4(flux.r1.(fw) < par.ep & abs(flux.vas) < nanmax(nanmax(abs(flux.vas)))/4) = 1; % note that evap is defined negative into atmosphere
        % elseif strcmp(type, 'gcm')
        %     rcae.(fw).vas4(flux.r1.(fw) < par.ep & flux.vas < nanmax(nanmax(abs(flux.vas)))/4) = 1; % note that evap is defined negative into atmosphere
        % end
        % if strcmp(fw, 'dse'); rcae.(fw).vas4(flux.r1.mse > par.ga) = -1;
        % else rcae.(fw).vas4(flux.r1.(fw) > par.ga) = -1; end;

        if size(flux.r1.(fw), 1)==length(par.lat_std) & size(flux.r1.(fw), 2)==12 % this criteria only works for zonally and time averaged data because vh is only a function of lat
            % add additional criteria for RCE that horizontal transport is weak
            rcae.(fw).vh2 = zeros(size(flux.r1.(fw)));
            % if strcmp(type, 'era5') | strcmp(type, 'erai')
                rcae.(fw).vh2(flux.r1.(fw) < par.ep & abs(vh.(fw)) < nanmax(nanmax(abs(vh.(fw))))/2 ) = 1;
            % elseif strcmp(type, 'gcm')
            %     rcae.(fw).vh2(flux.r1.(fw) < par.ep & abs(vh.(fw)) < nanmax(nanmax(abs(vh.(fw))))/2 ) = 1;
            % end
            if strcmp(fw, 'dse'); rcae.(fw).vh2(flux.r1.mse > par.ga) = -1;
            else; rcae.(fw).vh2(flux.r1.(fw) > par.ga) = -1; end;

            rcae.(fw).vh3 = zeros(size(flux.r1.(fw)));
            % if strcmp(type, 'era5') | strcmp(type, 'erai')
                rcae.(fw).vh3(flux.r1.(fw) < par.ep & abs(vh.(fw)) < nanmax(nanmax(abs(vh.(fw))))/3 ) = 1;
            % elseif strcmp(type, 'gcm')
            %     rcae.(fw).vh3(flux.r1.(fw) < par.ep & abs(vh.(fw)) < nanmax(nanmax(abs(vh.(fw))))/3 ) = 1;
            % end
            if strcmp(fw, 'dse') rcae.(fw).vh3(flux.r1.mse > par.ga) = -1;
            else rcae.(fw).vh3(flux.r1.(fw) > par.ga) = -1; end;

            rcae.(fw).vh4 = zeros(size(flux.r1.(fw)));
            % if strcmp(type, 'era5') | strcmp(type, 'erai')
                rcae.(fw).vh4(flux.r1.(fw) < par.ep & abs(vh.(fw)) < nanmax(nanmax(abs(vh.(fw))))/4 ) = 1;
            % elseif strcmp(type, 'gcm')
            %     rcae.(fw).vh4(flux.r1.(fw) < par.ep & abs(vh.(fw)) < nanmax(nanmax(abs(vh.(fw))))/4 ) = 1;
            % end
            if strcmp(fw, 'dse'); rcae.(fw).vh4(flux.r1.mse > par.ga) = -1;
            else; rcae.(fw).vh4(flux.r1.(fw) > par.ga) = -1; end;
        end

    end
end % define RCE and RAE
function rcae = def_rcae_alt_recomp_r1(type, flux, vh, par) % recalculate R1 at this step
    if any(strcmp(type, {'era5', 'erai'})); f_vec = par.era.fw;
    elseif strcmp(type, 'gcm'); f_vec = par.gcm.fw;
    elseif strcmp(type, 'echam'); f_vec = par.echam.fw; end
    for f = f_vec; fw = f{1};
        % calculate R1 again to test importance of order of operations
        if strcmp(fw, 'mse2'); r1 = flux.res.(fw)./flux.lw;
        else r1 = flux.res.(fw)./flux.ra.(fw); end

        % identify locations of RCE using threshold epsilon (ep)
        rcae.(fw).def = zeros(size(r1));
        rcae.(fw).def(r1 < par.ep) = 1;
        % identify locations of RAE using threshold gamma (ga)
        if strcmp(fw, 'dse'); rcae.(fw).def(flux.r1.mse > par.ga) = -1; % always use MSE for RAE
        else rcae.(fw).def(r1 > par.ga) = -1; end;

        % identify locations of RCE following Jakob et al. (2019) (abs(TEDIV) < 50 W/m^2)
        rcae.(fw).jak = zeros(size(r1));
        rcae.(fw).jak(flux.res.(fw) < 50) = 1;
        % identify locations of RAE using threshold gamma (ga)
        if strcmp(fw, 'dse'); rcae.(fw).jak(flux.r1.mse > par.ga) = -1; % always use MSE for RAE
        else rcae.(fw).jak(r1 > par.ga) = -1; end;

        % identify locations of rce following jakob et al. (2019) (abs(tediv) < 50 w/m^2)
        rcae.(fw).jak30 = zeros(size(r1));
        rcae.(fw).jak30(flux.res.(fw) < 30) = 1;
        % identify locations of rae using threshold gamma (ga)
        if strcmp(fw, 'dse'); rcae.(fw).jak30(flux.r1.mse > par.ga) = -1; % always use mse for rae
        else rcae.(fw).jak30(r1 > par.ga) = -1; end;

        % identify locations of rce following jakob et al. (2019) (abs(tediv) < 50 w/m^2)
        rcae.(fw).jak10 = zeros(size(r1));
        rcae.(fw).jak10(flux.res.(fw) < 10) = 1;
        % identify locations of rae using threshold gamma (ga)
        if strcmp(fw, 'dse'); rcae.(fw).jak10(flux.r1.mse > par.ga) = -1; % always use mse for rae
        else rcae.(fw).jak10(r1 > par.ga) = -1; end;

        % add additional criteria for RCE that P-E>0
        rcae.(fw).pe = zeros(size(r1));
        if strcmp(type, 'era5') | strcmp(type, 'erai')
            rcae.(fw).pe(r1 < par.ep & (flux.cp+flux.lsp+flux.e>0)) = 1; % note that evap is defined negative into atmosphere
        elseif strcmp(type, 'gcm')
            rcae.(fw).pe(r1 < par.ep & (flux.pr-flux.evspsbl>0)) = 1;
        elseif strcmp(type, 'echam')
            rcae.(fw).pe(r1 < par.ep & (flux.aprc + flux.aprl -flux.evap>0)) = 1;
        end
        if strcmp(fw, 'dse'); rcae.(fw).pe(flux.r1.mse > par.ga) = -1;
        else rcae.(fw).pe(r1 > par.ga) = -1; end;

        % add additional criteria for RCE that (P_ls - E)<<1
        rcae.(fw).cp = zeros(size(r1));
        if strcmp(type, 'era5') | strcmp(type, 'erai')
            rcae.(fw).cp(r1 < par.ep & abs(flux.lsp./flux.cp<par.ep_cp)) = 1; % note that evap is defined negative into atmosphere
        elseif strcmp(type, 'gcm')
            rcae.(fw).cp(r1 < par.ep & abs((flux.pr-flux.prc)./flux.prc<par.ep_cp)) = 1; % note that evap is defined negative into atmosphere
        elseif strcmp(type, 'echam')
            rcae.(fw).cp(r1 < par.ep & abs(flux.aprl./flux.aprc<par.ep_cp)) = 1; % note that evap is defined negative into atmosphere
        end
        if strcmp(fw, 'dse') rcae.(fw).cp(flux.r1.mse > par.ga) = -1;
        else rcae.(fw).cp(r1 > par.ga) = -1; end;

        % % add additional criteria for RCE that w500<0 (ascent)
        % rcae.(fw).w500 = zeros(size(r1));
        % if strcmp(type, 'era5') | strcmp(type, 'erai')
        %     rcae.(fw).w500(r1 < par.ep & flux.w500 < 0) = 1; % note that evap is defined negative into atmosphere
        % elseif strcmp(type, 'gcm')
        %     rcae.(fw).w500(r1 < par.ep & flux.w500 < 0) = 1; % note that evap is defined negative into atmosphere
        % end
        % if strcmp(fw, 'dse'); rcae.(fw).w500(flux.r1.mse > par.ga) = -1;
        % else rcae.(fw).w500(r1 > par.ga) = -1; end;

        % % add additional criteria for RCE that vh<max(vh)/2 (weak meridional velocity)
        % rcae.(fw).vas2 = zeros(size(r1));
        % if strcmp(type, 'era5') | strcmp(type, 'erai')
        %     rcae.(fw).vas2(r1 < par.ep & abs(flux.vas) < nanmax(nanmax(abs(flux.vas)))/2) = 1; % note that evap is defined negative into atmosphere
        % elseif strcmp(type, 'gcm')
        %     rcae.(fw).vas2(r1 < par.ep & flux.vas < nanmax(nanmax(abs(flux.vas)))/2) = 1; % note that evap is defined negative into atmosphere
        % end
        % if strcmp(fw, 'dse'); rcae.(fw).vas2(flux.r1.mse > par.ga) = -1;
        % else rcae.(fw).vas2(r1 > par.ga) = -1; end;

        % % add additional criteria for RCE that vh<max(vh)/4 (weak meridional velocity)
        % rcae.(fw).vas4 = zeros(size(r1));
        % if strcmp(type, 'era5') | strcmp(type, 'erai')
        %     rcae.(fw).vas4(r1 < par.ep & abs(flux.vas) < nanmax(nanmax(abs(flux.vas)))/4) = 1; % note that evap is defined negative into atmosphere
        % elseif strcmp(type, 'gcm')
        %     rcae.(fw).vas4(r1 < par.ep & flux.vas < nanmax(nanmax(abs(flux.vas)))/4) = 1; % note that evap is defined negative into atmosphere
        % end
        % if strcmp(fw, 'dse'); rcae.(fw).vas4(flux.r1.mse > par.ga) = -1;
        % else rcae.(fw).vas4(r1 > par.ga) = -1; end;

        if size(r1, 1)==length(par.lat_std) & size(r1, 2)==12 % this criteria only works for zonally and time averaged data because vh is only a function of lat
            % add additional criteria for RCE that horizontal transport is weak
            rcae.(fw).vh2 = zeros(size(r1));
            % if strcmp(type, 'era5') | strcmp(type, 'erai')
                rcae.(fw).vh2(r1 < par.ep & abs(vh.(fw)) < nanmax(nanmax(abs(vh.(fw))))/2 ) = 1;
            % elseif strcmp(type, 'gcm')
            %     rcae.(fw).vh2(r1 < par.ep & abs(vh.(fw)) < nanmax(nanmax(abs(vh.(fw))))/2 ) = 1;
            % end
            if strcmp(fw, 'dse'); rcae.(fw).vh2(flux.r1.mse > par.ga) = -1;
            else; rcae.(fw).vh2(r1 > par.ga) = -1; end;

            rcae.(fw).vh3 = zeros(size(r1));
            % if strcmp(type, 'era5') | strcmp(type, 'erai')
                rcae.(fw).vh3(r1 < par.ep & abs(vh.(fw)) < nanmax(nanmax(abs(vh.(fw))))/3 ) = 1;
            % elseif strcmp(type, 'gcm')
            %     rcae.(fw).vh3(r1 < par.ep & abs(vh.(fw)) < nanmax(nanmax(abs(vh.(fw))))/3 ) = 1;
            % end
            if strcmp(fw, 'dse') rcae.(fw).vh3(flux.r1.mse > par.ga) = -1;
            else rcae.(fw).vh3(r1 > par.ga) = -1; end;

            rcae.(fw).vh4 = zeros(size(r1));
            % if strcmp(type, 'era5') | strcmp(type, 'erai')
                rcae.(fw).vh4(r1 < par.ep & abs(vh.(fw)) < nanmax(nanmax(abs(vh.(fw))))/4 ) = 1;
            % elseif strcmp(type, 'gcm')
            %     rcae.(fw).vh4(r1 < par.ep & abs(vh.(fw)) < nanmax(nanmax(abs(vh.(fw))))/4 ) = 1;
            % end
            if strcmp(fw, 'dse'); rcae.(fw).vh4(flux.r1.mse > par.ga) = -1;
            else; rcae.(fw).vh4(r1 > par.ga) = -1; end;
        end

    end
end % define RCE and RAE
function land_mask = remove_land(lat, lon, nt)

    [lat2dgrid, lon2dgrid] = meshgrid(lat, lon);

    % land_mask =~ circshift(landmask(lat2dgrid, lon2dgrid, 100), length(lon)/2, 1);
    land_vec = land_or_ocean(lat2dgrid(:), lon2dgrid(:), 5, 0);
    land_mask = reshape(land_vec, [size(lat2dgrid)]);
    land_mask = double(land_mask);
    land_mask(land_mask==0) = nan;
    land_mask = repmat(land_mask, [1 1 nt]); % land mask in dims (lon x lat x time)

end % compute land mask (make land nan)
function ocean_mask = remove_ocean(lat, lon, nt)

	[lat2dgrid, lon2dgrid] = meshgrid(lat, lon);

	% ocean_mask = circshift(landmask(lat2dgrid, lon2dgrid, 100), length(lon)/2, 1);
    ocean_vec = ~land_or_ocean(lat2dgrid(:), lon2dgrid(:), 5, 0);
    ocean_mask = reshape(ocean_vec, [size(lat2dgrid)]);
	ocean_mask = double(ocean_mask);
	ocean_mask(ocean_mask==0) = nan;
	ocean_mask = repmat(ocean_mask, [1 1 nt]); % ocean mask in dims (lon x lat x time)

end % compute ocean mask (make ocean nan)

function ma_out = calc_ma_dew(ma_in, plev, par, type) % compute moist adiabat based on dew point temperature
    if any(strcmp(type, {'erai', 'era5'}))
        ps = ma_in.sp; tas = ma_in.t2m; tdas = ma_in.d2m;
    elseif strcmp(type, 'echam')
        ps = ma_in.aps; tas = ma_in.temp2; tdas = ma_in.dew2;
    end

    if isnan(ps) | isnan(tas) | isnan(tdas)
        ma_out = nan(size(plev)); % if any of the initial values are nan, ouput nan profile
    else
        pa_cd(1) = ps; % set initial pressure
        ta_cd(1) = tas; % set initial temperature
        esat_cd(1) = calc_esat(ta_cd(1), par.frz); % calculate initial saturation vapor pressure
        qsat_cd(1) = calc_q(pa_cd(1), esat_cd(1), par); % calculate initial specific humidity

        % convert dew point temperature to relative humidity
        e_cd(1) = calc_esat(tdas, par.frz); % calculate initial vapor pressure
        rh(1) = e_cd(1)/esat_cd(1); % set initial relative humidity

        r(1) = par.eps * rh(1) * esat_cd(1) / (pa_cd(1) - rh(1) * esat_cd(1));
        deriv_dry(1) = dry(pa_cd(1), ta_cd(1), par);
        [pa_d, ta_pa_d]=ode45(@(pa, T) dry(pa, T, par), par.pa_span, ta_cd(1));
        i = 2;
        while rh(i-1) < 1
            idx_lcl = i;
            pa_cd(i) = pa_cd(i-1) + par.dpa;
            ta_cd(i, 1) = ta_cd(i-1) + par.dpa * deriv_dry(i-1);
            esat_cd(i, 1) = calc_esat(ta_cd(i), par.frz);
            qsat_cd(i, 1) = calc_q(pa_cd(i), esat_cd(i), par);
            r(i, 1) = r(i-1);
            rh(i, 1) = r(i) * pa_cd(i) / (esat_cd(i) * (r(i) + par.eps));
            deriv_dry(i, 1) = dry(pa_cd(i), ta_cd(i), par);
            i = i + 1;
        end
        pa_lcl = pa_cd(end);
        if strcmp(par.ma_type, 'reversible')
            par.r_t = r(end); % moisture at LCL is the total moisture of rising parcel
            [pa_cs, ta_cs] = ode45(@(pa, T) sat_rev(pa, T, par.frz, par), [pa_lcl, par.pa_span(end)], ta_cd(end));
        elseif strcmp(par.ma_type, 'pseudo')
            [pa_cs, ta_cs] = ode45(@(pa, T) sat_pseudo(pa, T, par.frz, par), [pa_lcl, par.pa_span(end)], ta_cd(end));
        elseif strcmp(par.ma_type, 'std')
            [pa_cs, ta_cs] = ode45(@(pa, T) sat_std(pa, T, par.frz, par), [pa_lcl, par.pa_span(end)], ta_cd(end));
        else
            error(sprintf('The chosen moist adiabat type, %s, is not available. Choose between reversible, pseudo, or std.', par.ma_type));
        end
        for i = 1:length(pa_cs)
            if strcmp(par.ma_type, 'reversible')
                [deriv_cs(i, 1), qsat_cs(i, 1)] = sat_rev(pa_cs(i), ta_cs(i), par.frz, par);
            elseif strcmp(par.ma_type, 'pseudo')
                [deriv_cs(i, 1), qsat_cs(i, 1)] = sat_pseudo(pa_cs(i), ta_cs(i), par.frz, par);
            elseif strcmp(par.ma_type, 'std')
                [deriv_cs(i, 1), qsat_cs(i, 1)] = sat_std(pa_cs(i), ta_cs(i), par.frz, par);
            else
                error(sprintf('The chosen moist adiabat type, %s, is not available. Choose between reversible, pseudo, or std.', par.ma_type));
            end
            deriv_csd(i, 1) = dry(pa_cs(i), ta_cs(i), par);
        end
        ta_s(:, 1) = real([ta_cd; ta_cs(2:end)]);
        pa_s(:, 1) = [pa_cd(:); pa_cs(2:end)];
        qsat_s(:, 1) = [qsat_cd(:, 1); qsat_cs(2:end)];
        dtadpa_s(:, 1) = [deriv_dry; deriv_cs(2:end)];

        ta_s = interp1(pa_s, ta_s, plev, 'linear', nan);
        qsat_s = interp1(pa_s, qsat_s, plev, 'linear', nan);
        dtadpa_s = interp1(pa_s, dtadpa_s, plev, 'linear', nan);

        % OUTPUT
        ma_out = ta_s;
    end

end
function ma_out = calc_maz_dew(ma_in, z_int, par, type)
    if any(strcmp(type, {'erai', 'era5'}))
        ps = ma_in.sp; tas = ma_in.t2m; tdas = ma_in.d2m;
    elseif strcmp(type, 'echam')
        ps = ma_in.aps; tas = ma_in.temp2; tdas = ma_in.dew2;
    end

    if isnan(ps) | isnan(tas) | isnan(tdas)
        ma_out = nan(size(z_int)); % if any of the initial values are nan, ouput nan profile
    else
        % integrate temperature profiles over height
        z_cd(1) = ma_in.zs;

        pa_cd(1) = ps;
        ta_cd(1) = tas;
        esat_cd(1) = calc_esat(ta_cd(1), par.frz);
        qsat_cd(1) = calc_q(pa_cd(1), esat_cd(1), par);

        % convert dew point temperature to relative humidity
        e_cd(1) = calc_esat(tdas, par.frz); % calculate initial vapor pressure
        rh(1) = e_cd(1)/esat_cd(1); % set initial relative humidity

        r(1) = par.eps * rh(1) * esat_cd(1) / (pa_cd(1) - rh(1) * esat_cd(1));
        dtdz_dpdz(1,:) = dry_z([ta_cd(1), pa_cd(1)], par);
        T_p0 = [ta_cd(1), pa_cd(1)];
        [z_d, ta_pa_d]=ode45(@(z, T_p) dry_z(T_p, par), par.z_span, T_p0);
        i = 2;
        while rh(i-1) < 1
            idx_lcl = i;
            z_cd(i, 1) = z_cd(i-1) + par.dz;
            pa_cd(i, 1) = pa_cd(i-1) + par.dz * dtdz_dpdz(i-1, 2);
            ta_cd(i, 1) = ta_cd(i-1) + par.dz * dtdz_dpdz(i-1, 1);
            esat_cd(i, 1) = calc_esat(ta_cd(i), par.frz);
            qsat_cd(i, 1) = calc_q(pa_cd(i), esat_cd(i), par);
            r(i, 1) = r(i-1);
            rh(i, 1) = r(i) * pa_cd(i) / (esat_cd(i) * (r(i) + par.eps));
            dtdz_dpdz(i,:) = dry_z([ta_cd(i), pa_cd(i)], par);
            i = i + 1;
        end
        z_lcl = z_cd(end);
        pa_lcl = pa_cd(end);
        ta_pa_cd = [ta_cd, pa_cd];
        deriv_dry(:,1) = dtdz_dpdz(:,1);
        if strcmp(par.ma_type, 'reversible')
            par.r_t = r(end); % moisture at LCL is the total moisture of rising parcel
            [z_cs, ta_pa_cs] = ode45(@(z, T_p) sat_rev_z(T_p, par.frz, par), [z_lcl, par.z_span(end)], ta_pa_cd(end,:));
            % h1 = interp1(ta_pa_cs(:,1), z_cs, 240); % middle of troposphere
        elseif strcmp(par.ma_type, 'pseudo')
            [z_cs, ta_pa_cs] = ode45(@(z, T_p) sat_pseudo_z(T_p, par.frz, par), [z_lcl, par.z_span(end)], ta_pa_cd(end,:));
            % h1 = interp1(ta_pa_cs(:,1), z_cs, 240); % middle of troposphere
        elseif strcmp(par.ma_type, 'std')
            [z_cs, ta_pa_cs] = ode45(@(z, T_p) sat_z(T_p, par.frz, par), [z_lcl, par.z_span(end)], ta_pa_cd(end,:));
            % h1 = interp1(ta_pa_cs(:,1), z_cs, 240); % middle of troposphere
        else
            error(sprintf('The chosen moist adiabat type, %s, is not available. Choose between reversible, pseudo, or std.', zr.ma_type));
        end
        for i = 1:length(z_cs)
            if strcmp(par.ma_type, 'reversible')
                [dt_dp_cs(i, :), qsat_cs(i, 1)] = sat_rev_z(ta_pa_cs(i,:), par.frz, par);
            elseif strcmp(par.ma_type, 'pseudo')
                [dt_dp_cs(i, :), qsat_cs(i, 1)] = sat_pseudo_z(ta_pa_cs(i,:), par.frz, par);
            elseif strcmp(par.ma_type, 'std')
                [dt_dp_cs(i, :), qsat_cs(i, 1)] = sat_z(ta_pa_cs(i,:), par.frz, par);
            else
                error(sprintf('The chosen moist adiabat type, %s, is not available. Choose between reversible, pseudo, or std.', par.ma_type));
            end
            dt_dp_csd(i, :) = dry_z([ta_pa_cs(i,:)], par);
        end
        deriv_cs(:,1) = dt_dp_cs(:,1);
        deriv_csd(:,1) = dt_dp_csd(:,1);
        z_s(:, 1) = [z_cd; z_cs(2:end)];
        ta_s(:, 1) = real([ta_cd; ta_pa_cs(2:end,1)]);
        pa_s(:, 1) = [pa_cd(:); ta_pa_cs(2:end,2)];
        qsat_s(:, 1) = [qsat_cd(:, 1); qsat_cs(2:end)];
        dtadz_s(:, 1) = [deriv_dry; deriv_cs(2:end)];

        pa_s = interp1(z_s, pa_s, z_int, 'linear');
        ta_s = interp1(z_s, ta_s, z_int, 'linear');
        qsat_s = interp1(z_s, qsat_s, z_int, 'linear');
        dtadz_s = interp1(z_s, dtadz_s, z_int, 'linear');

        % output temperature
        ma_out = ta_s;
    end

end
function ma_out = calc_maz_dew_pa(ma_in, plev, par, type)
    if any(strcmp(type, {'erai', 'era5'}))
        ps = ma_in.sp; tas = ma_in.t2m; tdas = ma_in.d2m;
    elseif strcmp(type, 'echam')
        ps = ma_in.aps; tas = ma_in.temp2; tdas = ma_in.dew2;
    end

    if isnan(ps) | isnan(tas) | isnan(tdas)
        ma_out = nan(size(z_int)); % if any of the initial values are nan, ouput nan profile
    else
        % integrate temperature profiles over height
        z_cd(1) = ma_in.zs;

        pa_cd(1) = ps;
        ta_cd(1) = tas;
        esat_cd(1) = calc_esat(ta_cd(1), par.frz);
        qsat_cd(1) = calc_q(pa_cd(1), esat_cd(1), par);

        % convert dew point temperature to relative humidity
        e_cd(1) = calc_esat(tdas, par.frz); % calculate initial vapor pressure
        rh(1) = e_cd(1)/esat_cd(1); % set initial relative humidity

        r(1) = par.eps * rh(1) * esat_cd(1) / (pa_cd(1) - rh(1) * esat_cd(1));
        dtdz_dpdz(1,:) = dry_z([ta_cd(1), pa_cd(1)], par);
        T_p0 = [ta_cd(1), pa_cd(1)];
        [z_d, ta_pa_d]=ode45(@(z, T_p) dry_z(T_p, par), par.z_span, T_p0);
        i = 2;
        while rh(i-1) < 1
            idx_lcl = i;
            z_cd(i, 1) = z_cd(i-1) + par.dz;
            pa_cd(i, 1) = pa_cd(i-1) + par.dz * dtdz_dpdz(i-1, 2);
            ta_cd(i, 1) = ta_cd(i-1) + par.dz * dtdz_dpdz(i-1, 1);
            esat_cd(i, 1) = calc_esat(ta_cd(i), par.frz);
            qsat_cd(i, 1) = calc_q(pa_cd(i), esat_cd(i), par);
            r(i, 1) = r(i-1);
            rh(i, 1) = r(i) * pa_cd(i) / (esat_cd(i) * (r(i) + par.eps));
            dtdz_dpdz(i,:) = dry_z([ta_cd(i), pa_cd(i)], par);
            i = i + 1;
        end
        z_lcl = z_cd(end);
        pa_lcl = pa_cd(end);
        ta_pa_cd = [ta_cd, pa_cd];
        deriv_dry(:,1) = dtdz_dpdz(:,1);
        if strcmp(par.ma_type, 'reversible')
            par.r_t = r(end); % moisture at LCL is the total moisture of rising parcel
            [z_cs, ta_pa_cs] = ode45(@(z, T_p) sat_rev_z(T_p, par.frz, par), [z_lcl, par.z_span(end)], ta_pa_cd(end,:));
            % h1 = interp1(ta_pa_cs(:,1), z_cs, 240); % middle of troposphere
        elseif strcmp(par.ma_type, 'pseudo')
            [z_cs, ta_pa_cs] = ode45(@(z, T_p) sat_pseudo_z(T_p, par.frz, par), [z_lcl, par.z_span(end)], ta_pa_cd(end,:));
            % h1 = interp1(ta_pa_cs(:,1), z_cs, 240); % middle of troposphere
        elseif strcmp(par.ma_type, 'std')
            [z_cs, ta_pa_cs] = ode45(@(z, T_p) sat_z(T_p, par.frz, par), [z_lcl, par.z_span(end)], ta_pa_cd(end,:));
            % h1 = interp1(ta_pa_cs(:,1), z_cs, 240); % middle of troposphere
        else
            error(sprintf('The chosen moist adiabat type, %s, is not available. Choose between reversible, pseudo, or std.', zr.ma_type));
        end
        for i = 1:length(z_cs)
            if strcmp(par.ma_type, 'reversible')
                [dt_dp_cs(i, :), qsat_cs(i, 1)] = sat_rev_z(ta_pa_cs(i,:), par.frz, par);
            elseif strcmp(par.ma_type, 'pseudo')
                [dt_dp_cs(i, :), qsat_cs(i, 1)] = sat_pseudo_z(ta_pa_cs(i,:), par.frz, par);
            elseif strcmp(par.ma_type, 'std')
                [dt_dp_cs(i, :), qsat_cs(i, 1)] = sat_z(ta_pa_cs(i,:), par.frz, par);
            else
                error(sprintf('The chosen moist adiabat type, %s, is not available. Choose between reversible, pseudo, or std.', par.ma_type));
            end
            dt_dp_csd(i, :) = dry_z([ta_pa_cs(i,:)], par);
        end
        deriv_cs(:,1) = dt_dp_cs(:,1);
        deriv_csd(:,1) = dt_dp_csd(:,1);
        z_s(:, 1) = [z_cd; z_cs(2:end)];
        ta_s(:, 1) = real([ta_cd; ta_pa_cs(2:end,1)]);
        pa_s(:, 1) = [pa_cd(:); ta_pa_cs(2:end,2)];
        qsat_s(:, 1) = [qsat_cd(:, 1); qsat_cs(2:end)];
        dtadz_s(:, 1) = [deriv_dry; deriv_cs(2:end)];

        nanfilter = isnan(real(pa_s));
        pa_s(nanfilter) = [];
        dtadz_s(nanfilter) = [];

        dtadz_s = -1e3*interp1(real(pa_s), real(dtadz_s), plev, 'linear'); % output in K/km at standard pressure grid

        % output temperature
        ma_out = dtadz_s;
    end

end
function ma_out = calc_maz_dew_si(ma_in, plev, par, type, grid)
    if any(strcmp(type, {'erai', 'era5'}))
        ps = ma_in.sp; tas = ma_in.t2m; tdas = ma_in.d2m;
    elseif strcmp(type, 'echam')
        ps = ma_in.aps; tas = ma_in.temp2; tdas = ma_in.dew2;
    end

    if isnan(ps) | isnan(tas) | isnan(tdas)
        ma_out = nan(size(z_int)); % if any of the initial values are nan, ouput nan profile
    else
        % integrate temperature profiles over height
        z_cd(1) = ma_in.zs;

        pa_cd(1) = ps;
        ta_cd(1) = tas;
        esat_cd(1) = calc_esat(ta_cd(1), par.frz);
        qsat_cd(1) = calc_q(pa_cd(1), esat_cd(1), par);

        % convert dew point temperature to relative humidity
        e_cd(1) = calc_esat(tdas, par.frz); % calculate initial vapor pressure
        rh(1) = e_cd(1)/esat_cd(1); % set initial relative humidity

        r(1) = par.eps * rh(1) * esat_cd(1) / (pa_cd(1) - rh(1) * esat_cd(1));
        dtdz_dpdz(1,:) = dry_z([ta_cd(1), pa_cd(1)], par);
        T_p0 = [ta_cd(1), pa_cd(1)];
        [z_d, ta_pa_d]=ode45(@(z, T_p) dry_z(T_p, par), par.z_span, T_p0);
        i = 2;
        while rh(i-1) < 1
            idx_lcl = i;
            z_cd(i, 1) = z_cd(i-1) + par.dz;
            pa_cd(i, 1) = pa_cd(i-1) + par.dz * dtdz_dpdz(i-1, 2);
            ta_cd(i, 1) = ta_cd(i-1) + par.dz * dtdz_dpdz(i-1, 1);
            esat_cd(i, 1) = calc_esat(ta_cd(i), par.frz);
            qsat_cd(i, 1) = calc_q(pa_cd(i), esat_cd(i), par);
            r(i, 1) = r(i-1);
            rh(i, 1) = r(i) * pa_cd(i) / (esat_cd(i) * (r(i) + par.eps));
            dtdz_dpdz(i,:) = dry_z([ta_cd(i), pa_cd(i)], par);
            i = i + 1;
        end
        z_lcl = z_cd(end);
        pa_lcl = pa_cd(end);
        ta_pa_cd = [ta_cd, pa_cd];
        deriv_dry(:,1) = dtdz_dpdz(:,1);
        if strcmp(par.ma_type, 'reversible')
            par.r_t = r(end); % moisture at LCL is the total moisture of rising parcel
            [z_cs, ta_pa_cs] = ode45(@(z, T_p) sat_rev_z(T_p, par.frz, par), [z_lcl, par.z_span(end)], ta_pa_cd(end,:));
            % h1 = interp1(ta_pa_cs(:,1), z_cs, 240); % middle of troposphere
        elseif strcmp(par.ma_type, 'pseudo')
            [z_cs, ta_pa_cs] = ode45(@(z, T_p) sat_pseudo_z(T_p, par.frz, par), [z_lcl, par.z_span(end)], ta_pa_cd(end,:));
            % h1 = interp1(ta_pa_cs(:,1), z_cs, 240); % middle of troposphere
        elseif strcmp(par.ma_type, 'std')
            [z_cs, ta_pa_cs] = ode45(@(z, T_p) sat_z(T_p, par.frz, par), [z_lcl, par.z_span(end)], ta_pa_cd(end,:));
            % h1 = interp1(ta_pa_cs(:,1), z_cs, 240); % middle of troposphere
        else
            error(sprintf('The chosen moist adiabat type, %s, is not available. Choose between reversible, pseudo, or std.', zr.ma_type));
        end
        for i = 1:length(z_cs)
            if strcmp(par.ma_type, 'reversible')
                [dt_dp_cs(i, :), qsat_cs(i, 1)] = sat_rev_z(ta_pa_cs(i,:), par.frz, par);
            elseif strcmp(par.ma_type, 'pseudo')
                [dt_dp_cs(i, :), qsat_cs(i, 1)] = sat_pseudo_z(ta_pa_cs(i,:), par.frz, par);
            elseif strcmp(par.ma_type, 'std')
                [dt_dp_cs(i, :), qsat_cs(i, 1)] = sat_z(ta_pa_cs(i,:), par.frz, par);
            else
                error(sprintf('The chosen moist adiabat type, %s, is not available. Choose between reversible, pseudo, or std.', par.ma_type));
            end
            dt_dp_csd(i, :) = dry_z([ta_pa_cs(i,:)], par);
        end
        deriv_cs(:,1) = dt_dp_cs(:,1);
        deriv_csd(:,1) = dt_dp_csd(:,1);
        z_s(:, 1) = [z_cd; z_cs(2:end)];
        ta_s(:, 1) = real([ta_cd; ta_pa_cs(2:end,1)]);
        pa_s(:, 1) = [pa_cd(:); ta_pa_cs(2:end,2)];
        qsat_s(:, 1) = [qsat_cd(:, 1); qsat_cs(2:end)];
        dtadz_s(:, 1) = [deriv_dry; deriv_cs(2:end)];

        nanfilter = isnan(real(pa_s));
        pa_s(nanfilter) = [];
        dtadz_s(nanfilter) = [];

        dtadz_s = -1e3*interp1(real(pa_s)/ps, real(dtadz_s), grid.dim3.si, 'linear'); % output in K/km at standard pressure grid

        % output temperature
        ma_out = dtadz_s;
    end

end
function ma_out = calc_ma_hurs(ma_in, plev, par)
    if isnan(ma_in.ps) | isnan(ma_in.tas) | isnan(ma_in.hurs)
        ma_out = nan(size(plev)); % if any of the initial values are nan, ouput nan profile
    else
        pa_cd(1) = ma_in.ps; % set initial pressure
        ta_cd(1) = ma_in.tas; % set initial temperature
        esat_cd(1) = calc_esat(ta_cd(1), par.frz); % calculate initial saturation vapor pressure
        qsat_cd(1) = calc_q(pa_cd(1), esat_cd(1), par); % calculate initial specific humidity
        rh(1) = ma_in.hurs/100; % set initial relative humidity

        r(1) = par.eps * rh(1) * esat_cd(1) / (pa_cd(1) - rh(1) * esat_cd(1));
        deriv_dry(1) = dry(pa_cd(1), ta_cd(1), par);
        [pa_d, ta_pa_d]=ode45(@(pa, T) dry(pa, T, par), par.pa_span, ta_cd(1));
        i = 2;
        while rh(i-1) < 1
            idx_lcl = i;
            pa_cd(i) = pa_cd(i-1) + par.dpa;
            ta_cd(i, 1) = ta_cd(i-1) + par.dpa * deriv_dry(i-1);
            esat_cd(i, 1) = calc_esat(ta_cd(i), par.frz);
            qsat_cd(i, 1) = calc_q(pa_cd(i), esat_cd(i), par);
            r(i, 1) = r(i-1);
            rh(i, 1) = r(i) * pa_cd(i) / (esat_cd(i) * (r(i) + par.eps));
            deriv_dry(i, 1) = dry(pa_cd(i), ta_cd(i), par);
            i = i + 1;
        end
        pa_lcl = pa_cd(end);
        if strcmp(par.ma_type, 'reversible')
            par.r_t = r(end); % moisture at LCL is the total moisture of rising parcel
            [pa_cs, ta_cs] = ode45(@(pa, T) sat_rev(pa, T, par.frz, par), [pa_lcl, par.pa_span(end)], ta_cd(end));
        elseif strcmp(par.ma_type, 'pseudo')
            [pa_cs, ta_cs] = ode45(@(pa, T) sat_pseudo(pa, T, par.frz, par), [pa_lcl, par.pa_span(end)], ta_cd(end));
        elseif strcmp(par.ma_type, 'std')
            [pa_cs, ta_cs] = ode45(@(pa, T) sat_std(pa, T, par.frz, par), [pa_lcl, par.pa_span(end)], ta_cd(end));
        else
            error(sprintf('The chosen moist adiabat type, %s, is not available. Choose between reversible, pseudo, or std.', par.ma_type));
        end
        for i = 1:length(pa_cs)
            if strcmp(par.ma_type, 'reversible')
                [deriv_cs(i, 1), qsat_cs(i, 1)] = sat_rev(pa_cs(i), ta_cs(i), par.frz, par);
            elseif strcmp(par.ma_type, 'pseudo')
                [deriv_cs(i, 1), qsat_cs(i, 1)] = sat_pseudo(pa_cs(i), ta_cs(i), par.frz, par);
            elseif strcmp(par.ma_type, 'std')
                [deriv_cs(i, 1), qsat_cs(i, 1)] = sat_std(pa_cs(i), ta_cs(i), par.frz, par);
            else
                error(sprintf('The chosen moist adiabat type, %s, is not available. Choose between reversible, pseudo, or std.', par.ma_type));
            end
            deriv_csd(i, 1) = dry(pa_cs(i), ta_cs(i), par);
        end
        ta_s(:, 1) = [ta_cd; ta_cs(2:end)];
        pa_s(:, 1) = [pa_cd(:); pa_cs(2:end)];
        qsat_s(:, 1) = [qsat_cd(:, 1); qsat_cs(2:end)];
        dtadpa_s(:, 1) = [deriv_dry; deriv_cs(2:end)];

        ta_s = interp1(pa_s, ta_s, plev, 'linear', nan);
        qsat_s = interp1(pa_s, qsat_s, plev, 'linear', nan);
        dtadpa_s = interp1(pa_s, dtadpa_s, plev, 'linear', nan);

        % OUTPUT
        ma_out = ta_s;
    end

end % compute moist adiabat based on RH
function ma_out = calc_maz_hurs(ma_in, z_int, par)
    if isnan(ma_in.ps) | isnan(ma_in.tas) | isnan(ma_in.hurs) | isnan(ma_in.zs)
        ma_out = nan(size(z_int)); % if any of the initial values are nan, ouput nan profile
    else
        % integrate temperature profiles over height
        z_cd(1) = ma_in.zs;
        pa_cd(1) = ma_in.ps;
        ta_cd(1) = ma_in.tas;
        esat_cd(1) = calc_esat(ta_cd(1), par.frz);
        qsat_cd(1) = calc_q(pa_cd(1), esat_cd(1), par);
        rh(1) = ma_in.hurs / 100;
        r(1) = par.eps * rh(1) * esat_cd(1) / (pa_cd(1) - rh(1) * esat_cd(1));
        dtdz_dpdz(1,:) = dry_z([ta_cd(1), pa_cd(1)], par);
        T_p0 = [ta_cd(1), pa_cd(1)];
        [z_d, ta_pa_d]=ode45(@(z, T_p) dry_z(T_p, par), par.z_span, T_p0);
        i = 2;
        while rh(i-1) < 1
            idx_lcl = i;
            z_cd(i, 1) = z_cd(i-1) + par.dz;
            pa_cd(i, 1) = pa_cd(i-1) + par.dz * dtdz_dpdz(i-1, 2);
            ta_cd(i, 1) = ta_cd(i-1) + par.dz * dtdz_dpdz(i-1, 1);
            esat_cd(i, 1) = calc_esat(ta_cd(i), par.frz);
            qsat_cd(i, 1) = calc_q(pa_cd(i), esat_cd(i), par);
            r(i, 1) = r(i-1);
            rh(i, 1) = r(i) * pa_cd(i) / (esat_cd(i) * (r(i) + par.eps));
            dtdz_dpdz(i,:) = dry_z([ta_cd(i), pa_cd(i)], par);
            i = i + 1;
        end
        z_lcl = z_cd(end);
        pa_lcl = pa_cd(end);
        ta_pa_cd = [ta_cd, pa_cd];
        deriv_dry(:,1) = dtdz_dpdz(:,1);
        if strcmp(par.ma_type, 'reversible')
            par.r_t = r(end); % moisture at LCL is the total moisture of rising parcel
            [z_cs, ta_pa_cs] = ode45(@(z, T_p) sat_rev_z(T_p, par.frz, par), [z_lcl, par.z_span(end)], ta_pa_cd(end,:));
            % h1 = interp1(ta_pa_cs(:,1), z_cs, 240); % middle of troposphere
        elseif strcmp(par.ma_type, 'pseudo')
            [z_cs, ta_pa_cs] = ode45(@(z, T_p) sat_pseudo_z(T_p, par.frz, par), [z_lcl, par.z_span(end)], ta_pa_cd(end,:));
            % h1 = interp1(ta_pa_cs(:,1), z_cs, 240); % middle of troposphere
        elseif strcmp(par.ma_type, 'std')
            [z_cs, ta_pa_cs] = ode45(@(z, T_p) sat_z(T_p, par.frz, par), [z_lcl, par.z_span(end)], ta_pa_cd(end,:));
            % h1 = interp1(ta_pa_cs(:,1), z_cs, 240); % middle of troposphere
        else
            error(sprintf('The chosen moist adiabat type, %s, is not available. Choose between reversible, pseudo, or std.', zr.ma_type));
        end
        for i = 1:length(z_cs)
            if strcmp(par.ma_type, 'reversible')
                [dt_dp_cs(i, :), qsat_cs(i, 1)] = sat_rev_z(ta_pa_cs(i,:), par.frz, par);
            elseif strcmp(par.ma_type, 'pseudo')
                [dt_dp_cs(i, :), qsat_cs(i, 1)] = sat_pseudo_z(ta_pa_cs(i,:), par.frz, par);
            elseif strcmp(par.ma_type, 'std')
                [dt_dp_cs(i, :), qsat_cs(i, 1)] = sat_z(ta_pa_cs(i,:), par.frz, par);
            else
                error(sprintf('The chosen moist adiabat type, %s, is not available. Choose between reversible, pseudo, or std.', par.ma_type));
            end
            dt_dp_csd(i, :) = dry_z([ta_pa_cs(i,:)], par);
        end
        deriv_cs(:,1) = dt_dp_cs(:,1);
        deriv_csd(:,1) = dt_dp_csd(:,1);
        z_s(:, 1) = [z_cd; z_cs(2:end)];
        ta_s(:, 1) = real([ta_cd; ta_pa_cs(2:end,1)]);
        pa_s(:, 1) = [pa_cd(:); ta_pa_cs(2:end,2)];
        qsat_s(:, 1) = [qsat_cd(:, 1); qsat_cs(2:end)];
        dtadz_s(:, 1) = [deriv_dry; deriv_cs(2:end)];

        pa_s = interp1(z_s, pa_s, z_int, 'linear');
        ta_s = interp1(z_s, ta_s, z_int, 'linear');
        qsat_s = interp1(z_s, qsat_s, z_int, 'linear');
        dtadz_s = interp1(z_s, dtadz_s, z_int, 'linear');

        % output temperature
        ma_out = ta_s;
    end

end
function ma_out = calc_maz_hurs_pa(ma_in, plev, par)
    if isnan(ma_in.ps) | isnan(ma_in.tas) | isnan(ma_in.hurs) | isnan(ma_in.zs)
        ma_out = nan(size(z_int)); % if any of the initial values are nan, ouput nan profile
    else
        % integrate temperature profiles over height
        z_cd(1) = ma_in.zs;
        pa_cd(1) = ma_in.ps;
        ta_cd(1) = ma_in.tas;
        esat_cd(1) = calc_esat(ta_cd(1), par.frz);
        qsat_cd(1) = calc_q(pa_cd(1), esat_cd(1), par);
        rh(1) = ma_in.hurs / 100;
        r(1) = par.eps * rh(1) * esat_cd(1) / (pa_cd(1) - rh(1) * esat_cd(1));
        dtdz_dpdz(1,:) = dry_z([ta_cd(1), pa_cd(1)], par);
        T_p0 = [ta_cd(1), pa_cd(1)];
        [z_d, ta_pa_d]=ode45(@(z, T_p) dry_z(T_p, par), par.z_span, T_p0);
        i = 2;
        while rh(i-1) < 1
            idx_lcl = i;
            z_cd(i, 1) = z_cd(i-1) + par.dz;
            pa_cd(i, 1) = pa_cd(i-1) + par.dz * dtdz_dpdz(i-1, 2);
            ta_cd(i, 1) = ta_cd(i-1) + par.dz * dtdz_dpdz(i-1, 1);
            esat_cd(i, 1) = calc_esat(ta_cd(i), par.frz);
            qsat_cd(i, 1) = calc_q(pa_cd(i), esat_cd(i), par);
            r(i, 1) = r(i-1);
            rh(i, 1) = r(i) * pa_cd(i) / (esat_cd(i) * (r(i) + par.eps));
            dtdz_dpdz(i,:) = dry_z([ta_cd(i), pa_cd(i)], par);
            i = i + 1;
        end
        z_lcl = z_cd(end);
        pa_lcl = pa_cd(end);
        ta_pa_cd = [ta_cd, pa_cd];
        deriv_dry(:,1) = dtdz_dpdz(:,1);
        if strcmp(par.ma_type, 'reversible')
            par.r_t = r(end); % moisture at LCL is the total moisture of rising parcel
            [z_cs, ta_pa_cs] = ode45(@(z, T_p) sat_rev_z(T_p, par.frz, par), [z_lcl, par.z_span(end)], ta_pa_cd(end,:));
            % h1 = interp1(ta_pa_cs(:,1), z_cs, 240); % middle of troposphere
        elseif strcmp(par.ma_type, 'pseudo')
            [z_cs, ta_pa_cs] = ode45(@(z, T_p) sat_pseudo_z(T_p, par.frz, par), [z_lcl, par.z_span(end)], ta_pa_cd(end,:));
            % h1 = interp1(ta_pa_cs(:,1), z_cs, 240); % middle of troposphere
        elseif strcmp(par.ma_type, 'std')
            [z_cs, ta_pa_cs] = ode45(@(z, T_p) sat_z(T_p, par.frz, par), [z_lcl, par.z_span(end)], ta_pa_cd(end,:));
            % h1 = interp1(ta_pa_cs(:,1), z_cs, 240); % middle of troposphere
        else
            error(sprintf('The chosen moist adiabat type, %s, is not available. Choose between reversible, pseudo, or std.', zr.ma_type));
        end
        for i = 1:length(z_cs)
            if strcmp(par.ma_type, 'reversible')
                [dt_dp_cs(i, :), qsat_cs(i, 1)] = sat_rev_z(ta_pa_cs(i,:), par.frz, par);
            elseif strcmp(par.ma_type, 'pseudo')
                [dt_dp_cs(i, :), qsat_cs(i, 1)] = sat_pseudo_z(ta_pa_cs(i,:), par.frz, par);
            elseif strcmp(par.ma_type, 'std')
                [dt_dp_cs(i, :), qsat_cs(i, 1)] = sat_z(ta_pa_cs(i,:), par.frz, par);
            else
                error(sprintf('The chosen moist adiabat type, %s, is not available. Choose between reversible, pseudo, or std.', par.ma_type));
            end
            dt_dp_csd(i, :) = dry_z([ta_pa_cs(i,:)], par);
        end
        deriv_cs(:,1) = dt_dp_cs(:,1);
        deriv_csd(:,1) = dt_dp_csd(:,1);
        z_s(:, 1) = [z_cd; z_cs(2:end)];
        ta_s(:, 1) = real([ta_cd; ta_pa_cs(2:end,1)]);
        pa_s(:, 1) = [pa_cd(:); ta_pa_cs(2:end,2)];
        qsat_s(:, 1) = [qsat_cd(:, 1); qsat_cs(2:end)];
        dtadz_s(:, 1) = [deriv_dry; deriv_cs(2:end)];

        nanfilter = isnan(real(pa_s));
        pa_s(nanfilter) = [];
        dtadz_s(nanfilter) = [];

        dtadz_s = -1e3*interp1(real(pa_s), real(dtadz_s), plev, 'linear'); % output in K/km at standard pressure grid

        % output temperature
        ma_out = dtadz_s;
    end

end
function ma_out = calc_maz_hurs_si(ma_in, plev, par, grid)
    if isnan(ma_in.ps) | isnan(ma_in.tas) | isnan(ma_in.hurs) | isnan(ma_in.zs)
        ma_out = nan(size(z_int)); % if any of the initial values are nan, ouput nan profile
    else
        % integrate temperature profiles over height
        z_cd(1) = ma_in.zs;
        pa_cd(1) = ma_in.ps;
        ta_cd(1) = ma_in.tas;
        esat_cd(1) = calc_esat(ta_cd(1), par.frz);
        qsat_cd(1) = calc_q(pa_cd(1), esat_cd(1), par);
        rh(1) = ma_in.hurs / 100;
        r(1) = par.eps * rh(1) * esat_cd(1) / (pa_cd(1) - rh(1) * esat_cd(1));
        dtdz_dpdz(1,:) = dry_z([ta_cd(1), pa_cd(1)], par);
        T_p0 = [ta_cd(1), pa_cd(1)];
        [z_d, ta_pa_d]=ode45(@(z, T_p) dry_z(T_p, par), par.z_span, T_p0);
        i = 2;
        while rh(i-1) < 1
            idx_lcl = i;
            z_cd(i, 1) = z_cd(i-1) + par.dz;
            pa_cd(i, 1) = pa_cd(i-1) + par.dz * dtdz_dpdz(i-1, 2);
            ta_cd(i, 1) = ta_cd(i-1) + par.dz * dtdz_dpdz(i-1, 1);
            esat_cd(i, 1) = calc_esat(ta_cd(i), par.frz);
            qsat_cd(i, 1) = calc_q(pa_cd(i), esat_cd(i), par);
            r(i, 1) = r(i-1);
            rh(i, 1) = r(i) * pa_cd(i) / (esat_cd(i) * (r(i) + par.eps));
            dtdz_dpdz(i,:) = dry_z([ta_cd(i), pa_cd(i)], par);
            i = i + 1;
        end
        z_lcl = z_cd(end);
        pa_lcl = pa_cd(end);
        ta_pa_cd = [ta_cd, pa_cd];
        deriv_dry(:,1) = dtdz_dpdz(:,1);
        if strcmp(par.ma_type, 'reversible')
            par.r_t = r(end); % moisture at LCL is the total moisture of rising parcel
            [z_cs, ta_pa_cs] = ode45(@(z, T_p) sat_rev_z(T_p, par.frz, par), [z_lcl, par.z_span(end)], ta_pa_cd(end,:));
            % h1 = interp1(ta_pa_cs(:,1), z_cs, 240); % middle of troposphere
        elseif strcmp(par.ma_type, 'pseudo')
            [z_cs, ta_pa_cs] = ode45(@(z, T_p) sat_pseudo_z(T_p, par.frz, par), [z_lcl, par.z_span(end)], ta_pa_cd(end,:));
            % h1 = interp1(ta_pa_cs(:,1), z_cs, 240); % middle of troposphere
        elseif strcmp(par.ma_type, 'std')
            [z_cs, ta_pa_cs] = ode45(@(z, T_p) sat_z(T_p, par.frz, par), [z_lcl, par.z_span(end)], ta_pa_cd(end,:));
            % h1 = interp1(ta_pa_cs(:,1), z_cs, 240); % middle of troposphere
        else
            error(sprintf('The chosen moist adiabat type, %s, is not available. Choose between reversible, pseudo, or std.', zr.ma_type));
        end
        for i = 1:length(z_cs)
            if strcmp(par.ma_type, 'reversible')
                [dt_dp_cs(i, :), qsat_cs(i, 1)] = sat_rev_z(ta_pa_cs(i,:), par.frz, par);
            elseif strcmp(par.ma_type, 'pseudo')
                [dt_dp_cs(i, :), qsat_cs(i, 1)] = sat_pseudo_z(ta_pa_cs(i,:), par.frz, par);
            elseif strcmp(par.ma_type, 'std')
                [dt_dp_cs(i, :), qsat_cs(i, 1)] = sat_z(ta_pa_cs(i,:), par.frz, par);
            else
                error(sprintf('The chosen moist adiabat type, %s, is not available. Choose between reversible, pseudo, or std.', par.ma_type));
            end
            dt_dp_csd(i, :) = dry_z([ta_pa_cs(i,:)], par);
        end
        deriv_cs(:,1) = dt_dp_cs(:,1);
        deriv_csd(:,1) = dt_dp_csd(:,1);
        z_s(:, 1) = [z_cd; z_cs(2:end)];
        ta_s(:, 1) = real([ta_cd; ta_pa_cs(2:end,1)]);
        pa_s(:, 1) = [pa_cd(:); ta_pa_cs(2:end,2)];
        qsat_s(:, 1) = [qsat_cd(:, 1); qsat_cs(2:end)];
        dtadz_s(:, 1) = [deriv_dry; deriv_cs(2:end)];

        nanfilter = isnan(real(pa_s));
        pa_s(nanfilter) = [];
        dtadz_s(nanfilter) = [];

        dtadz_s = -1e3*interp1(real(pa_s)/ma_in.ps, real(dtadz_s), grid.dim3.si, 'linear'); % output in K/km at standard pressure grid

        % output temperature
        ma_out = dtadz_s;
    end

end

% old functions
function make_dtdz(type, par) % compute model lapse rate in lon x lat x plev x mon
% calculate moist adiabats
    if strcmp(type, 'era5') | strcmp(type, 'erai')
        foldername = sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/', type, par.lat_interp);
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s', type);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s', type);
    elseif strcmp(type, 'gcm')
        foldername = sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s/%s/', type, par.model, par.gcm.clim, par.lat_interp);
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s/%s', type, par.model, par.gcm.clim);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s', type, par.model, par.gcm.clim);
    elseif strcmp(type, 'echam')
        foldername = sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s/', type, par.echam.clim, par.lat_interp);
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s', type, par.echam.clim);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s', type, par.echam.clim);
    end
    load(sprintf('%s/grid.mat', prefix)); % read grid data
    load(sprintf('%s/tai.mat', prefix)); % read moist adiabat

    dtdz = -1e3*(tai(:,:,2:end,:)-tai(:,:,1:end-1,:))./(zgi(:,:,2:end,:)-zgi(:,:,1:end-1,:)); % lapse rate in K/km

    dtdz = permute(dtdz, [3 1 2 4]); % bring height forward
    dtdz = interp1(1/2*(par.pa(2:end)+par.pa(1:end-1)), dtdz, par.pa);
    dtdz = permute(dtdz, [2 3 1 4]); % bring height back to 3rd

    if any(strcmp(type, {'era5', 'erai'})); newdir=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s', type);
    elseif strcmp(type, 'gcm'); newdir=sprintf('/project2/tas1/miyawaki/projects/002/data/read/gcm/%s/%s', par.model, par.gcm.clim);
    elseif strcmp(type, 'echam'); newdir=sprintf('/project2/tas1/miyawaki/projects/002/data/read/echam/%s', par.echam.clim); end;
    if ~exist(newdir, 'dir'); mkdir(newdir); end
    filename='dtdz.mat';
    save(sprintf('%s/%s', newdir, filename), 'dtdz', '-v7.3');

end
function make_dtdzi(type, par) % compute model lapse rate in lon x lat x plev x mon
% calculate moist adiabats
    if strcmp(type, 'era5') | strcmp(type, 'erai')
        foldername = sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/', type, par.lat_interp);
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s', type);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s', type);
        file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/%s/temp/%s_temp_%s.ymonmean.nc', type, type, par.(type).yr_span));
        fullpath=sprintf('%s/%s', file.folder, file.name);
        temp = ncread(fullpath, 't');
        file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/%s/zg/%s_zg_%s.ymonmean.nc', type, type, par.(type).yr_span));
        fullpath=sprintf('%s/%s', file.folder, file.name);
        zg = ncread(fullpath, 'z'); zg = zg/par.g; % convert from geopotential to height in m
    elseif strcmp(type, 'gcm')
        foldername = sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s/%s/', type, par.model, par.gcm.clim, par.lat_interp);
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s/%s', type, par.model, par.gcm.clim);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s', type, par.model, par.gcm.clim);
        var = 'ta';
        file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/gcm/%s/%s_Amon_%s_%s_r1i1p1_*.ymonmean.nc', par.model, var, par.model, par.gcm.clim));
        fullpath=sprintf('%s/%s', file.folder, file.name);
        temp = ncread(fullpath, var);
        var = 'zg';
        file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/gcm/%s/%s_Amon_%s_%s_r1i1p1_*.ymonmean.nc', par.model, var, par.model, par.gcm.clim));
        fullpath=sprintf('%s/%s', file.folder, file.name);
        zg = ncread(fullpath, var);
    elseif strcmp(type, 'echam')
        foldername = sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s/', type, par.echam.clim, par.lat_interp);
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s', type, par.echam.clim);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s', type, par.echam.clim);
        var = 't';
        file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/echam/ATM_rjg_%s_*.ymonmean.nc', par.echam.clim));
        fullpath=sprintf('%s/%s', file.folder, file.name);
        temp = double(ncread(fullpath, var));
        var = 'geopoth';
        file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/echam/ATM_rjg_%s_*.ymonmean.nc', par.echam.clim));
        fullpath=sprintf('%s/%s', file.folder, file.name);
        zg = double(ncread(fullpath, var));
    end
    load(sprintf('%s/grid.mat', prefix)); % read grid data
    load(sprintf('%s/srfc.mat', prefix)); % read surface variable data

    dtdzi = -1e3*(temp(:,:,2:end,:)-temp(:,:,1:end-1,:))./(zg(:,:,2:end,:)-zg(:,:,1:end-1,:)); % lapse rate in K/km

    dtdzi = permute(dtdzi, [3 1 2 4]); % bring height forward
    dtdzi = interp1(1/2*(grid.dim3.plev(2:end)+grid.dim3.plev(1:end-1)), dtdzi, grid.dim3.plev);
    dtdzi = permute(dtdzi, [2 3 1 4]); % bring height back to 3rd

    % create surface mask
    if strcmp(type, 'era5') | strcmp(type, 'erai')
        ps_vert = repmat(srfc.sp, [1 1 1 size(dtdzi, 3)]); % dims (lon x lat x time x plev)
        ps_vert = permute(ps_vert, [1 2 4 3]); % dims (lon x lat x plev x time)
        pa = double(permute(repmat(grid.dim3.plev, [1 size(srfc.ps)]), [2 3 1 4]));
    elseif strcmp(type, 'gcm')
        ps_vert = repmat(srfc.ps, [1 1 1 size(dtdzi, 3)]); % dims (lon x lat x time x plev)
        ps_vert = permute(ps_vert, [1 2 4 3]); % dims (lon x lat x plev x time)
        pa = permute(repmat(grid.dim3.plev, [1 size(srfc.ps)]), [2 3 1 4]);
    elseif strcmp(type, 'echam')
        ps_vert = repmat(srfc.aps, [1 1 1 size(dtdzi, 3)]); % dims (lon x lat x time x plev)
        ps_vert = permute(ps_vert, [1 2 4 3]); % dims (lon x lat x plev x time)
        pa = permute(repmat(grid.dim3.plev, [1 size(srfc.aps)]), [2 3 1 4]);
    end
    surface_mask = nan(size(dtdzi));
    surface_mask(pa < ps_vert) = 1;

    dtdzi = dtdzi.*surface_mask; % filter dtdzi with surface mask

    dtdzi = permute(dtdzi, [3 1 2 4]); % bring height forward
    dtdzi = interp1(grid.dim3.plev, dtdzi, par.pa, 'spline', nan);
    dtdzi = permute(dtdzi, [2 3 1 4]); % bring height back to 3rd

    if any(strcmp(type, {'era5', 'erai'})); newdir=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s', type);
    elseif strcmp(type, 'gcm'); newdir=sprintf('/project2/tas1/miyawaki/projects/002/data/read/gcm/%s/%s', par.model, par.gcm.clim);
    elseif strcmp(type, 'echam'); newdir=sprintf('/project2/tas1/miyawaki/projects/002/data/read/echam/%s', par.echam.clim); end;
    if ~exist(newdir, 'dir'); mkdir(newdir); end
    filename='dtdzi.mat';
    save(sprintf('%s/%s', newdir, filename), 'dtdzi', '-v7.3');

end
function make_dtdzz(type, par) % compute model lapse rate in lat x plev x mon
% calculate moist adiabats
    if strcmp(type, 'era5') | strcmp(type, 'erai')
        foldername = sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/', type, par.lat_interp);
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s', type);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s', type);
        file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/%s/temp/%s_temp_%s.ymonmean.nc', type, type, par.(type).yr_span));
        fullpath=sprintf('%s/%s', file.folder, file.name);
        temp = ncread(fullpath, 't');
        file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/%s/zg/%s_zg_%s.ymonmean.nc', type, type, par.(type).yr_span));
        fullpath=sprintf('%s/%s', file.folder, file.name);
        zg = ncread(fullpath, 'z'); zg = zg/par.g; % convert from geopotential to height in m
    elseif strcmp(type, 'gcm')
        foldername = sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s/%s/', type, par.model, par.gcm.clim, par.lat_interp);
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s/%s', type, par.model, par.gcm.clim);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s', type, par.model, par.gcm.clim);
        var = 'ta';
        file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/gcm/%s/%s_Amon_%s_%s_r1i1p1_*.ymonmean.nc', par.model, var, par.model, par.gcm.clim));
        fullpath=sprintf('%s/%s', file.folder, file.name);
        temp = ncread(fullpath, var);
        var = 'zg';
        file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/gcm/%s/%s_Amon_%s_%s_r1i1p1_*.ymonmean.nc', par.model, var, par.model, par.gcm.clim));
        fullpath=sprintf('%s/%s', file.folder, file.name);
        zg = ncread(fullpath, var);
    elseif strcmp(type, 'echam')
        foldername = sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s/', type, par.echam.clim, par.lat_interp);
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s', type, par.echam.clim);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s', type, par.echam.clim);
        var = 't';
        file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/echam/ATM_rjg_%s_*.ymonmean.nc', par.echam.clim));
        fullpath=sprintf('%s/%s', file.folder, file.name);
        temp = double(ncread(fullpath, var));
        var = 'geopoth';
        file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/echam/ATM_rjg_%s_*.ymonmean.nc', par.echam.clim));
        fullpath=sprintf('%s/%s', file.folder, file.name);
        zg = double(ncread(fullpath, var));
    end
    load(sprintf('%s/grid.mat', prefix)); % read grid data
    load(sprintf('%s/srfc.mat', prefix)); % read surface variable data

    temp = squeeze(nanmean(temp, 1)); % zonal average
    zg = squeeze(nanmean(zg, 1)); % zonal average

    dtdzz = -1e3*(temp(:,2:end,:)-temp(:,1:end-1,:))./(zg(:,2:end,:)-zg(:,1:end-1,:)); % lapse rate in K/km

    dtdzz = permute(dtdzz, [2 1 3]); % bring height forward
    dtdzz = interp1(1/2*(grid.dim3.plev(2:end)+grid.dim3.plev(1:end-1)), dtdzz, grid.dim3.plev);
    dtdzz = permute(dtdzz, [2 1 3]); % bring height back to 2nd

    % create surface mask
    if strcmp(type, 'era5') | strcmp(type, 'erai')
        ps_vert = repmat(srfc.sp, [1 1 1 size(dtdzz, 2)]); % dims (lon x lat x time x plev)
        ps_vert = permute(ps_vert, [1 2 4 3]); % dims (lon x lat x plev x time)
        pa = double(permute(repmat(grid.dim3.plev, [1 size(srfc.sp)]), [2 3 1 4]));
    elseif strcmp(type, 'gcm')
        ps_vert = repmat(srfc.ps, [1 1 1 size(dtdzz, 2)]); % dims (lon x lat x time x plev)
        ps_vert = permute(ps_vert, [1 2 4 3]); % dims (lon x lat x plev x time)
        pa = permute(repmat(grid.dim3.plev, [1 size(srfc.ps)]), [2 3 1 4]);
    elseif strcmp(type, 'echam')
        ps_vert = repmat(srfc.aps, [1 1 1 size(dtdzz, 2)]); % dims (lon x lat x time x plev)
        ps_vert = permute(ps_vert, [1 2 4 3]); % dims (lon x lat x plev x time)
        pa = permute(repmat(grid.dim3.plev, [1 size(srfc.aps)]), [2 3 1 4]);
    end
    surface_mask = nan(size(pa));
    surface_mask(pa < ps_vert) = 1;

    dtdzz = repmat(dtdzz, [1 1 1 size(pa,1)]); % repeat in longitude
    dtdzz = permute(dtdzz, [4 1 2 3]); % bring lon to 1st

    dtdzz = dtdzz.*surface_mask; % filter dtdzz with surface mask

    dtdzz = permute(dtdzz, [3 1 2 4]); % bring height forward
    dtdzz = interp1(grid.dim3.plev, dtdzz, par.pa, 'spline', nan);
    dtdzz = permute(dtdzz, [2 3 1 4]); % bring height back to 3rd

    if any(strcmp(type, {'era5', 'erai'})); newdir=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s', type);
    elseif strcmp(type, 'gcm'); newdir=sprintf('/project2/tas1/miyawaki/projects/002/data/read/gcm/%s/%s', par.model, par.gcm.clim);
    elseif strcmp(type, 'echam'); newdir=sprintf('/project2/tas1/miyawaki/projects/002/data/read/echam/%s', par.echam.clim); end;
    if ~exist(newdir, 'dir'); mkdir(newdir); end
    filename='dtdzz.mat';
    save(sprintf('%s/%s', newdir, filename), 'dtdzz', '-v7.3');

end
function make_dtdzsi_old(type, par) % compute model lapse rate in lon x lat x plev x mon
% calculate moist adiabats
    if strcmp(type, 'era5') | strcmp(type, 'erai')
        foldername = sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/', type, par.lat_interp);
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s', type);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s', type);
        file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/%s/temp/%s_temp_%s.ymonmean.nc', type, type, par.(type).yr_span));
        fullpath=sprintf('%s/%s', file.folder, file.name);
        temp = ncread(fullpath, 't');
        file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/%s/zg/%s_zg_%s.ymonmean.nc', type, type, par.(type).yr_span));
        fullpath=sprintf('%s/%s', file.folder, file.name);
        zg = ncread(fullpath, 'z'); zg = zg/par.g; % convert from geopotential to height in m
    elseif strcmp(type, 'gcm')
        foldername = sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s/%s/', type, par.model, par.gcm.clim, par.lat_interp);
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s/%s', type, par.model, par.gcm.clim);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s', type, par.model, par.gcm.clim);
        var = 'ta';
        file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/gcm/%s/%s_Amon_%s_%s_r1i1p1_*.ymonmean.nc', par.model, var, par.model, par.gcm.clim));
        fullpath=sprintf('%s/%s', file.folder, file.name);
        temp = ncread(fullpath, var);
        var = 'zg';
        file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/gcm/%s/%s_Amon_%s_%s_r1i1p1_*.ymonmean.nc', par.model, var, par.model, par.gcm.clim));
        fullpath=sprintf('%s/%s', file.folder, file.name);
        zg = ncread(fullpath, var);
    elseif strcmp(type, 'echam')
        foldername = sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s/', type, par.echam.clim, par.lat_interp);
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s', type, par.echam.clim);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s', type, par.echam.clim);
        var = 't';
        file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/echam/ATM_rjg_%s_*.ymonmean.nc', par.echam.clim));
        fullpath=sprintf('%s/%s', file.folder, file.name);
        temp = double(ncread(fullpath, var));
        var = 'geopoth';
        file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/echam/ATM_rjg_%s_*.ymonmean.nc', par.echam.clim));
        fullpath=sprintf('%s/%s', file.folder, file.name);
        zg = double(ncread(fullpath, var));
    end
    load(sprintf('%s/grid.mat', prefix)); % read grid data
    load(sprintf('%s/srfc.mat', prefix)); % read surface variable data
    % load(sprintf('%s/orog.mat', prefix)); % read surface orography

    % create surface mask
    if strcmp(type, 'era5') | strcmp(type, 'erai')
        ps_vert = repmat(srfc.sp, [1 1 1 size(temp, 3)]); % dims (lon x lat x time x plev)
        ps_vert = permute(ps_vert, [1 2 4 3]); % dims (lon x lat x plev x time)
        pa = double(permute(repmat(grid.dim3.plev, [1 size(srfc.sp)]), [2 3 1 4]));
    elseif strcmp(type, 'gcm')
        ps_vert = repmat(srfc.ps, [1 1 1 size(temp, 3)]); % dims (lon x lat x time x plev)
        ps_vert = permute(ps_vert, [1 2 4 3]); % dims (lon x lat x plev x time)
        pa = permute(repmat(grid.dim3.plev, [1 size(srfc.ps)]), [2 3 1 4]);
    elseif strcmp(type, 'echam')
        ps_vert = repmat(srfc.aps, [1 1 1 size(temp, 3)]); % dims (lon x lat x time x plev)
        ps_vert = permute(ps_vert, [1 2 4 3]); % dims (lon x lat x plev x time)
        pa = permute(repmat(grid.dim3.plev, [1 size(srfc.aps)]), [2 3 1 4]);
    end
    surface_mask = nan(size(pa));
    surface_mask(pa < ps_vert) = 1;

    ta_sm = temp .* surface_mask;
    zg_sm = zg .* surface_mask;

    % add tsurf data and interpolate to higher resolution vertical grid
    [pa_plus ta_plus zg_plus] = deal(nan([size(pa,1), size(pa,2) size(pa,3)+1 size(pa,4)])); % create empty grid with one extra vertical level
    pa_plus(:,:,1:end-1,:) = pa; % populate with standard pressure grid
    ta_plus(:,:,1:end-1,:) = ta_sm; % populate with standard atmospheric temperature
    zg_plus(:,:,1:end-1,:) = zg_sm; % populate with szgndard atmospheric temperature
    pa_plus(:,:,end,:) = ps_vert(:,:,1,:); % add surface pressure data into standard pressure grid
    if any(strcmp(type, {'era5', 'erai'})); ta_plus(:,:,end,:) = srfc.t2m(:,:,:); % add surface temperature data
    elseif strcmp(type, 'gcm'); ta_plus(:,:,end,:) = srfc.tas(:,:,:); % add surface temperature data
    elseif strcmp(type, 'echam'); ta_plus(:,:,end,:) = srfc.temp2(:,:,:); end % add surface temperature data
    % zg_plus(:,:,end,:) = srfc.zs(:,:,:); % add surface height data
    % zg_plus(:,:,end,:) = repmat(orog, [1 1 12]); % add surface height data
    zg_plus(:,:,end,:) = nan(size(srfc.zs)); % add surface height data
    ps_vert = permute(ps_vert, [3 1 2 4]); % bring plev dimension to front
    pa_plus = permute(pa_plus, [3 1 2 4]); % bring plev dimension to front
    ta_plus = permute(ta_plus, [3 1 2 4]); % bring plev dimension to front
    zg_plus = permute(zg_plus, [3 1 2 4]); % bring plev dimension to front
    [pa_plus sort_index] = sort(pa_plus, 1, 'descend'); % sort added surface pressure such that pressure decreases monotonically
    pb = CmdLineProgressBar("Sorting temperature with surface data added...");
    for lo=1:size(pa_plus,2)
        pb.print(lo, size(pa_plus,2));
        for la=1:size(pa_plus,3)
            for mo=1:size(pa_plus,4)
                ta_plus(:,lo,la,mo) = ta_plus(sort_index(:,lo,la,mo),lo,la,mo); % sort temperature (has to be in loop because sort_index works for vector calls only)
                zg_plus(:,lo,la,mo) = zg_plus(sort_index(:,lo,la,mo),lo,la,mo); % sort temperature (has to be in loop because sort_index works for vector calls only)
            end
        end
    end

    dtdz = -1e3*(ta_plus(2:end,:,:,:)-ta_plus(1:end-1,:,:,:))./(zg_plus(2:end,:,:,:)-zg_plus(1:end-1,:,:,:)); % lapse rate in K/km

    pb = CmdLineProgressBar("Sorting and interpolating dtdz to new standard grid...");
    for lo=1:size(pa_plus,2)
        pb.print(lo, size(pa_plus,2));
        for la=1:size(pa_plus,3)
            for mo=1:size(pa_plus,4)
                dtdzsi(:,lo,la,mo) = interp1(1/2*(pa_plus(2:end,lo,la,mo)+pa_plus(1:end-1,lo,la,mo))./ps_vert(:,lo,la,mo), dtdz(:,lo,la,mo), grid.dim3.si, 'linear');
            end
        end
    end
    dtdzsi = permute(dtdzsi, [2 3 1 4]); % bring height to 3rd
    dtdzsi = fillmissing(dtdzsi, 'nearest');

    if any(strcmp(type, {'era5', 'erai'})); newdir=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s', type);
    elseif strcmp(type, 'gcm'); newdir=sprintf('/project2/tas1/miyawaki/projects/002/data/read/gcm/%s/%s', par.model, par.gcm.clim);
    elseif strcmp(type, 'echam'); newdir=sprintf('/project2/tas1/miyawaki/projects/002/data/read/echam/%s', par.echam.clim); end;
    if ~exist(newdir, 'dir'); mkdir(newdir); end
    filename='dtdzsi.mat';
    save(sprintf('%s/%s', newdir, filename), 'dtdzsi', '-v7.3');

end
