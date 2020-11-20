clc; close all; clear variables;

addpath(genpath('/project2/tas1/miyawaki/matlab'));

gcm_info
echam_info

%% set parameters
if 1
% par.erai.yr_span = '2000_2012'; % spanning years for ERA-Interim
par.erai.yr_span = '2000_2018'; % spanning years for ERA-Interim
par.era5.yr_span = '2000_2018'; % spanning years for ERA5
par.merra2.yr_span = '2000_2018'; % spanning years for MERRA22
par.gcm.clim = 'piControl'; % choose either piControl or abrupt4xCO2
par.echam_clims = {'rp000141'}; % par.echam.all_mld; % choose from 20170908 (snowball), 20170915_2 (modern), or rp000*** (various mixed layer depth and with/without sea ice)
par.lat_interp = 'native'; % which latitudinal grid to interpolate to: native (no interpolation), don (donohoe, coarse), era (native ERA-Interim, fine), or std (custom, very fine)
par.lat_std = transpose(-90:0.25:90); % define standard latitude grid for 'std' interpolation
par.ep_swp = 0.1; %[0.25 0.3 0.35]; % threshold for RCE definition. RCE is defined as where abs(R1) < ep
par.ga_swp = 0.9; % optional threshold for RAE. If undefined, the default value is 1-par.ep
par.ep_cp = 0.5; % additional flag for RCE definition using convective precipitation. RCE is defined as where lsp/cp < ep_cp
par.ma_type = 'reversible'; % choose the type of moist adiabat: reversible, pseudo, or std
par.ma_init = 0.95; % initial starting level for moist adiabat ('surf' = start with dry adiabat until LCL or enter sigma level for starting saturated adiabat)
par.frz = 0; % consider latent heat of fusion in moist adiabat?
par.pa_span = [1000 10]*100; % pressure range for calculating moist adiabat
par.pa = linspace(1000,10,100)*1e2; % high resolution vertical grid to interpolate to
par.si = 1e-5*par.pa; % high resolution vertical grid to interpolate to
par.dpa = -10; % pressure increment for integrating dry adiabat section of moist adiabat (matters for how accurately the LCL is computed)
par.z_span = [0 25]*10^3; % height range for calculating moist adiabat
par.dz = 10; % pressure increment for integrating dry adiabat section of moist adiabat (matters for how accurately the LCL is computed)
par.do_surf = 0; % whether or not to calculate temperature profile in pressure grids including ts and ps data
par.si_eval = [0.8 0.85 0.9]; % sigma level for evaluating inversion strength (T(si_eval) - T(surface))
par.si_bl_swp = [0.85 0.9 0.95]; % sigma level to separate vertical average for close to moist adiabatic and stable surface stratifications
par.si_up = 0.4; % sigma level for upper boundary of vertical average for close to moist adiabatic
% par.era.fw = {'mse', 'dse', 'db13', 'db13s', 'db13t', 'div', 'divt', 'div79'};
% par.era.fw = {'div79', 'mse', 'dse', 'db13', 'db13s', 'db13t', 'div', 'divt'};
par.era.fw = {'mse', 'dse'};
par.merra2.fw = {'mse', 'dse'};
par.gcm.fw = {'mse', 'dse'};
par.echam.fw = {'mse', 'dse'};
par.cpd = 1005.7; par.cpv = 1870; par.cpl = 4186; par.cpi = 2108; par.Rd = 287; par.Rv = 461; par.g = 9.81; par.L = 2.501e6; par.a = 6357e3; par.eps = 0.622; % common constants, all in SI units for brevity
end

%% call functions
% comp_flux(par)
% ceres_flux(par)
% choose_disp(par)

type = 'era5'; % data type to run analysis on
choose_proc(type, par)
for k=1:length(par.echam_clims); par.echam.clim=par.echam_clims{k};
    type='echam';
    % disp(par.echam.clim)
    % choose_proc(type, par);
end
for k=1:length(par.gcm_models); par.model = par.gcm_models{k};
    type = 'gcm';
    % disp(par.model)
    % choose_proc(type, par)
end

% for i=1:length(par.si_bl_swp); par.si_bl = par.si_bl_swp(i);
%     type = 'merra2'; % data type to run analysis on
%     % choose_proc_si_bl(type, par)
%     for k=1:length(par.echam_clims); par.echam.clim=par.echam_clims{k};
%         type='echam';
%         disp(par.echam.clim)
%         choose_proc_si_bl(type, par);
%     end
%     for k=1:length(par.gcm_models); par.model = par.gcm_models{k};
%         type = 'gcm';
%         % disp(par.model)
%         % choose_proc_si_bl(type, par)
%     end
% end

% for i=1:length(par.ep_swp); par.ep = par.ep_swp(i); par.ga = par.ga_swp(i);
%     type = 'era5';
%     choose_proc_ep(type, par)
%     for k=1:length(par.echam_clims); par.echam.clim=par.echam_clims{k};
%         type='echam';
%         % disp(par.echam.clim)
%         % choose_proc_ep(type, par);
%     end
%     for k = 1:length(par.gcm_models); par.model = par.gcm_models{k};
%         type = 'gcm';
%         % disp(par.model)
%         % choose_proc_ep(type, par)
%     end
% end

function choose_proc(type, par)
    % save_mask(type, par) % save land and ocean masks once (faster than creating mask every time I need it)
    % proc_flux(type, par) % calculate energy fluxes in the vertically-integrated MSE budget using ERA-Interim data
    % proc_net_flux(type, par) % calculate net energy fluxes at TOA and surface
    proc_temp_mon_lat(type, par) % calculate mon x lat temperature profiles
    % make_masi(type, par) % calculate moist adiabats at every lon x lat x mon
    % proc_ma_mon_lat(type, par) % calculate mon x lat moist adiabats
    % proc_thetaeq_mon_lat(type, par) % calculate mon x lat potential equivalent temp profiles
    % proc_thetaeq_mon_lat(type, par) % calculate mon x lat equivalent potential temperature profiles
    % proc_dtdz_mon_lat(type, par) % calculate mon x lat lapse rate
    % proc_malr_mon_lat(type, par) % calculate mon x lat moist adiabats
    % make_tai(type, par) % calculate moist adiabat in lon x lat x mon
    % make_dtdzsi(type, par) % calculate model lapse rate and interpolate to sigma coordinates
    % make_malrsi(type, par) % calculate moist adiabatic lapse rate of model temperature sigma coordinates

    % proc_temp_mon_lat_interp(type, par) % calculate mon x lat temperature profiles
    % proc_temp_mon_lat_interp_mean(type, par) % calculate mon x lat temperature profiles
end % select which functions to run at a time
function choose_proc_si_bl(type, par)
    proc_ga_malr_diff_si_mon_lat(type, par) % calculate mon x lat gamma percentage difference
    proc_ga_dalr_bl_diff_si_mon_lat(type, par) % calculate mon x lat gamma percentage difference
    % proc_ga_trop_malr_diff_si_mon_lat(type, par) % calculate mon x lat gamma percentage difference
end
function choose_proc_ep(type, par)
    % proc_rcae(type, par) % calculate RCE and RAE regimes
    % proc_rcae_alt(type, par) % calculate RCE and RAE regimes (all divergence allowed for RCE)
    proc_ta_si(type, par) % calculate RCE and RAE temperature profiles
    % proc_ma_si(type, par) % calculate moist adiabats corresponding to RCE profiles
end % select which ep-functions to run at a time

%% define functions
function save_mask(type, par)
    if strcmp(type, 'era5') | strcmp(type, 'erai')
        foldername = sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/', type, par.lat_interp);
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s', type);
    elseif strcmp(type, 'merra2')
        foldername = sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/', type, par.lat_interp);
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s', type);
    elseif strcmp(type, 'gcm')
        foldername = sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s/%s/', type, par.model, par.gcm.clim, par.lat_interp);
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s/%s', type, par.model, par.gcm.clim);
    elseif strcmp(type, 'echam')
        foldername = sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s/', type, par.echam.clim, par.lat_interp);
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s', type, par.echam.clim);
    elseif strcmp(type, 'echam_ml') | strcmp(type, 'echam_pl')
        foldername = sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/', type, par.lat_interp);
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s', type);
    end

    load(sprintf('%s/grid.mat', prefix)); % read grid data
    if strcmp(par.lat_interp, 'std')
        lat = par.lat_std;
    else
        lat = grid.dim3.lat;
    end

    if strcmp(type, 'gcm') | strcmp(type, 'echam')
        load(sprintf('%s/sftlf.mat', prefix)); % load land fraction data

        mask.ocean = nan(size(sftlf)); mask.ocean(sftlf>0.5) = 1; mask.ocean=repmat(mask.ocean,[1 1 12]);
        mask.land = nan(size(mask.ocean)); mask.land(isnan(mask.ocean))=1;
    else
        mask.land = remove_land(lat, grid.dim3.lon, 12);
        mask.ocean = remove_ocean(lat, grid.dim3.lon, 12);
    end

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
        prefix_don=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s', 'erai');
        if any(strcmp(type, {'erai', 'era5'})) & strcmp(par.erai.yr_span, '1979_2018')
            don79 = load(sprintf('%s/dondiv79.mat', prefix_don)); % read Donohoe data 1979--2018
        elseif any(strcmp(type, {'erai', 'era5'})) & strcmp(par.erai.yr_span, '2000_2018')
            don79 = load(sprintf('%s/dondiv00.mat', prefix_don)); % read Donohoe data 1979--2018
        end
    elseif strcmp(type, 'merra2')
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s', type);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s', type, par.lat_interp);
    elseif strcmp(type, 'gcm')
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s/%s', type, par.model, par.gcm.clim);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s/%s', type, par.model, par.gcm.clim, par.lat_interp);
    elseif strcmp(type, 'echam')
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s', type, par.echam.clim);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s', type, par.echam.clim, par.lat_interp);
    end

    load(sprintf('%s/grid.mat', prefix)) % read grid data
    load(sprintf('%s/rad.mat', prefix)) % read radiation data
    load(sprintf('%s/hydro.mat', prefix)) % read hydrology data
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
    for fn = hydro_vars; fname = fn{1}; % interpolate to std lat
        flux.(fname) = permute(hydro.(fname), [2 1 3]);
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
    elseif strcmp(type, 'merra2')
        flux.stf.mse = flux.HFLUX + flux.EFLUX;
        flux.stf.dse = par.L*flux.PRECTOT + flux.HFLUX;
    elseif strcmp(type, 'gcm')
        flux.stf.mse = flux.hfls + flux.hfss; flux.stf.mse2 = flux.stf.mse;
        flux.stf.dse = par.L*flux.pr + flux.hfss;
    elseif strcmp(type, 'echam')
        flux.stf.mse = -(flux.ahfl + flux.ahfs); flux.stf.mse2 = flux.stf.mse;
        flux.stf.dse = par.L*(flux.aprc+flux.aprl) - flux.ahfs;
    end

    if any(strcmp(type, {'era5', 'erai'})); f_vec = par.era.fw;
    elseif strcmp(type, 'merra2'); f_vec = par.merra2.fw;
    elseif strcmp(type, 'gcm'); f_vec = par.gcm.fw;
    elseif strcmp(type, 'echam'); f_vec = par.echam.fw; end
    for f = f_vec; fw = f{1};
        if any(strcmp(type, {'era5', 'erai'}));
            flux.rtoa = flux.tsr + flux.ttr; % net flux at TOA
            flux.olr = flux.ttr;
            flux.swsfc = -flux.ssr;
            flux.lwsfc = -flux.str;
            flux.lw = flux.ttr - flux.str; flux.sw = flux.tsr-flux.ssr; % compute net shortwave and longwave flux through atmosphere
            if contains(fw, 'ceresrad'); flux.ra.(fw) = ceres.ra;  % compute net radiative cooling from radiative fluxes
            else; flux.ra.(fw) = flux.tsr - flux.ssr + flux.ttr - flux.str; end % use radiative cooling from CERES data
        elseif strcmp(type, 'merra2');
            flux.rtoa = flux.SWTNT - flux.LWTUP; % net flux at TOA
            flux.olr = -flux.LWTUP;
            flux.swsfc = -flux.SWGNT;
            flux.lwsfc = -flux.LWGNT;
            flux.lw = -flux.LWTUP - flux.LWGNT; flux.sw = flux.SWTNT - flux.SWGNT;
            flux.ra.(fw) = flux.SWTNT - flux.LWTUP - flux.SWGNT - flux.LWGNT; % radiative cooling
        elseif strcmp(type, 'gcm');
            flux.rtoa = flux.rsdt - flux.rsut - flux.rlut; % net flux at TOA
            flux.olr = -flux.rlut;
            flux.swsfc = flux.rsus - flux.rsds;
            flux.lwsfc = flux.rlus - flux.rlds;
            flux.lw = flux.rlus - flux.rlds - flux.rlut; flux.sw = flux.rsdt - flux.rsut + flux.rsus - flux.rsds;
            flux.ra.(fw) = flux.rsdt - flux.rsut + flux.rsus - flux.rsds + flux.rlus - flux.rlds - flux.rlut;
        elseif strcmp(type, 'echam');
            flux.rtoa = flux.trad0 + flux.srad0; % net flux at TOA
            flux.olr = flux.trad0;
            flux.swsfc = -flux.srads;
            flux.lwsfc = -flux.trads;
            flux.lw = flux.trad0 - flux.trads; flux.sw = flux.srad0 - flux.srads;
            flux.ra.(fw) = flux.lw + flux.sw;
        end % calculate atmospheric radiative cooling
        flux.rsfc = flux.swsfc + flux.lwsfc;

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
            flux.res.(fw) = flux.TEDIV + flux.tend; % use MSE flux divergence from DB13, use MSE tendency from ERA
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
        elseif strcmp(fw, 'div00erarad')
            flux.res.(fw) = flux.don79div + flux.tend; % use MSE tendency from ERA output and MSE flux divergence from Donohoe MSE flux transport data
            flux.stf.(fw) = flux.res.(fw) - flux.ra.(fw); % infer the surface turbulent fluxes
        end

        flux.shf.(fw) = flux.lwsfc + flux.stf.(fw); % surface LW and surface turbulent fluxes
        flux.sfc.(fw) = flux.rsfc + flux.stf.(fw); % net flux at surface

        if strcmp(fw, 'mse2')
            flux.r1.(fw) = (flux.res.(fw))./flux.lw; % calculate nondimensional number R1 disregarding MSE budget closure
            flux.r2.(fw) = flux.stf.(fw)./flux.lw; % calculate nondimensional number R2 disregarding MSE budget closure
        else
            flux.r1.(fw) = (flux.res.(fw))./flux.ra.(fw); % calculate nondimensional number R1 disregarding MSE budget closure
            flux.r2.(fw) = flux.stf.(fw)./flux.ra.(fw); % calculate nondimensional number R2 disregarding MSE budget closure
        end
        if any(strcmp(type, {'era5', 'erai'}));
            flux.ftoa.(fw) = flux.tsr + flux.ttr; flux.fsfc.(fw) = -flux.ssr - flux.str + flux.stf.(fw);
        elseif strcmp(type, 'merra2')
            flux.ftoa.(fw) = flux.SWTNT - flux.LWTUP;
            flux.fsfc.(fw) = -flux.SWGNT - flux.LWGNT + flux.stf.(fw);
        elseif strcmp(type, 'gcm');
            flux.ftoa.(fw) = flux.rsdt - flux.rsut - flux.rlut;
            flux.fsfc.(fw) = flux.rsus - flux.rsds + flux.rlus - flux.rlds + flux.stf.(fw);
        elseif strcmp(type, 'echam');
            flux.ftoa.(fw) = flux.trad0 + flux.srad0;
            flux.fsfc.(fw) = -flux.trads - flux.srads + flux.stf.(fw);
        end

        % linear decomposition of R1 seasonality
        flux.comp1.(fw) = (flux.res.(fw)-nanmean(flux.res.(fw),3))./nanmean(flux.ra.(fw),3);
        flux.comp2.(fw) = - nanmean(flux.res.(fw),3)./nanmean(flux.ra.(fw),3).^2 .* (flux.ra.(fw)-nanmean(flux.ra.(fw),3));
        % flux.comp2.(fw) = - nanmean(flux.res.(fw)./flux.ra.(fw).^2, 3) .* (flux.ra.(fw)-nanmean(flux.ra.(fw),3));

    end


    if strcmp(type, 'era5') | strcmp(type, 'erai');
        var_vec = {'sshf', 'slhf', 'cp', 'lsp', 'e', 'lw', 'sw', 'rtoa', 'olr', 'lwsfc', 'swsfc', 'tend', 'divt', 'divg', 'divq', 'TETEN', 'TEDIV', 'don79div'};
        foldername = sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/', type, par.lat_interp);
    elseif strcmp(type, 'merra2')
        var_vec = {'EFLUX', 'HFLUX', 'PRECCON', 'PRECTOT', 'EVAP', 'lw', 'sw', 'rtoa', 'olr', 'lwsfc', 'swsfc'};
        foldername = sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/', type, par.lat_interp);
    elseif strcmp(type, 'gcm')
        var_vec = {'hfls', 'hfss', 'prc', 'pr', 'evspsbl', 'lw', 'sw', 'rtoa', 'olr', 'lwsfc', 'swsfc'};
        foldername = sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s/%s/', type, par.model, par.gcm.clim, par.lat_interp);
    elseif strcmp(type, 'echam')
        var_vec = {'ahfl', 'ahfs', 'aprc', 'aprl', 'evap', 'lw', 'sw', 'rtoa', 'olr', 'lwsfc', 'swsfc'};
        foldername = sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s/', type, par.echam.clim, par.lat_interp);
    end

    for fn = var_vec; fname = fn{1};
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

    for fn = {'ra', 'stf', 'res', 'r1', 'r2', 'ftoa', 'fsfc', 'sfc', 'shf', 'comp1', 'comp2'}; fname = fn{1};
        if any(strcmp(type, {'era5', 'erai'})); f_vec = par.era.fw;
        elseif strcmp(type, 'merra2'); f_vec = par.merra2.fw;
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

                    if any(strcmp(type, {'erai', 'era5'}))
                        % remove global mean from flux divergence
                        resmean = permute(vh_mon.(land).(fw)(end,:), [2 1]);
                        resmean = repmat(resmean, [1 length(grid.dim3.lon) length(grid.dim3.lat)]);
                        resmean = permute(resmean, [2 3 1]);
                        flux_n.(land).res.(fw) = flux_n.(land).res.(fw) + resmean/(4*pi*par.a^2);

                        % re-compute northward MSE transport using the residual data
                        tediv_0 = fillmissing(flux_n.(land).res.(fw), 'constant', 0); % replace nans with 0s so missing data doesn't influence transport
                        % tediv_z = squeeze(trapz(deg2rad(grid.dim2.lon), tediv_0, 1)); % zonally integrate
                        tediv_z = squeeze(nanmean(tediv_0, 1)); % zonally average
                        vh_mon.(land).(fw) = cumtrapz(deg2rad(lat), 2*pi*par.a^2*cosd(lat).*tediv_z, 1); % cumulatively integrate
                    end

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
    elseif strcmp(type, 'gcm')
        foldername = sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s/%s/', type, par.model, par.gcm.clim, par.lat_interp);
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s/%s', type, par.model, par.gcm.clim);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s', type, par.model, par.gcm.clim);
    elseif strcmp(type, 'echam')
        foldername = sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s/', type, par.echam.clim, par.lat_interp);
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s', type, par.echam.clim);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s', type, par.echam.clim);
    elseif strcmp(type, 'echam_ml')
        foldername = sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/', type, par.lat_interp);
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s', type);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s', type);
    elseif strcmp(type, 'echam_pl')
        foldername = sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/', type, par.lat_interp);
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s', type);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s', type);
    end
    load(sprintf('%s/grid.mat', prefix)); % read grid data
    load(sprintf('%s/ta_si.mat', prefix)); tasi_orig = ta_si.spl; clear ta_si; % read temp in si coordinates
    load(sprintf('%s/pa_si.mat', prefix)); pasi_orig = pa_si; clear pa_si; % read temp in si coordinates
    load(sprintf('%s/srfc.mat', prefix)); % load surface data
    load(sprintf('%s/%s/masks.mat', prefix_proc, par.lat_interp)); % load land and ocean masks

    if strcmp(par.lat_interp, 'std')
        lat = par.lat_std;
    else
        lat = grid.dim3.lat;
    end

    % interpolate ta to standard lat grid
    tasi_orig = permute(tasi_orig, [2 1 3 4]);
    tasi_orig = interp1(grid.dim3.lat, tasi_orig, lat);
    tasi_orig = permute(tasi_orig, [2 1 3 4]);

    pasi_orig = permute(pasi_orig, [2 1 3 4]);
    pasi_orig = interp1(grid.dim3.lat, pasi_orig, lat);
    pasi_orig = permute(pasi_orig, [2 1 3 4]);

    tasi_sm.lo = tasi_orig; % surface is already masked in standard sigma coordinates
    pasi_sm.lo = pasi_orig; % surface is already masked in spndard sigma coordinates

    tasi_sm.lo = permute(tasi_sm.lo, [1 2 4 3]); % bring plev to last dimension

    pasi_sm.lo = permute(pasi_sm.lo, [1 2 4 3]); % bring plev to last dimension

    mask_t.land = nanmean(mask.land, 3);
    mask_t.ocean = nanmean(mask.ocean, 3);

    for l = {'lo'}; land = l{1}; % over land, over ocean, or both
        tasi.(land)= squeeze(nanmean(tasi_sm.(land), 1)); % zonal average
        pasi.(land)= squeeze(nanmean(pasi_sm.(land), 1)); % zonal average
    end

    for l = {'lo'}; land = l{1}; % over land, over ocean, or both
        % take time averages
        for t = {'ann', 'djf', 'jja', 'mam', 'son'}; time = t{1};
            if strcmp(time, 'ann')
                tasi_t.(land).(time) = squeeze(nanmean(tasi_sm.(land), 3));
                pasi_t.(land).(time) = squeeze(nanmean(pasi_sm.(land), 3));
            elseif strcmp(time, 'djf')
                tasi_shift.(land) = circshift(tasi_sm.(land), 1, 3);
                pasi_shift.(land) = circshift(pasi_sm.(land), 1, 3);
                tasi_t.(land).(time) = squeeze(nanmean(tasi_shift.(land)(:,:,1:3,:), 3));
                pasi_t.(land).(time) = squeeze(nanmean(pasi_shift.(land)(:,:,1:3,:), 3));
            elseif strcmp(time, 'jja')
                tasi_t.(land).(time) = squeeze(nanmean(tasi_sm.(land)(:,:,6:8,:), 3));
                pasi_t.(land).(time) = squeeze(nanmean(pasi_sm.(land)(:,:,6:8,:), 3));
            elseif strcmp(time, 'mam')
                tasi_t.(land).(time) = squeeze(nanmean(tasi_sm.(land)(:,:,3:5,:), 3));
                pasi_t.(land).(time) = squeeze(nanmean(pasi_sm.(land)(:,:,3:5,:), 3));
            elseif strcmp(time, 'son')
                tasi_t.(land).(time) = squeeze(nanmean(tasi_sm.(land)(:,:,9:11,:), 3));
                pasi_t.(land).(time) = squeeze(nanmean(pasi_sm.(land)(:,:,9:11,:), 3));
            end
        end
    end

    % save filtered data
    printname = [foldername 'ta_mon_lat'];
    if ~exist(foldername, 'dir')
        mkdir(foldername)
    end
    if par.do_surf; save(printname, 'tasi', 'lat');
    else save(printname, 'tasi', 'pasi', 'lat', '-v7.3'); end

    printname = [foldername 'ta_lon_lat'];
    if ~exist(foldername, 'dir')
        mkdir(foldername)
    end
    if par.do_surf; save(printname, 'tasi_t', 'lat');
    else save(printname, 'tasi_t', 'pasi_t', 'lat', '-v7.3'); end
end % compute mon x lat temperature field
function make_masi(type, par)
% compute moist adiabats at every lon, lat, mon
    if any(strcmp(type, {'era5', 'erai'}))
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s', type);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s', type);
        file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/%s/temp/%s_temp_%s.ymonmean.nc', type, type, par.(type).yr_span));
        fullpath=sprintf('%s/%s', file.folder, file.name);
        temp = ncread(fullpath, 't');
        file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/%s/zg/%s_zg_%s.ymonmean.nc', type, type, par.(type).yr_span));
        fullpath=sprintf('%s/%s', file.folder, file.name);
        zg = ncread(fullpath, 'z'); zg = zg/par.g; % convert from geopotential to height in m
    elseif strcmp(type, 'merra2')
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s', type);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s', type);
    elseif strcmp(type, 'gcm')
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
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s', type, par.echam.clim);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s', type, par.echam.clim);
        var = 't';
        if contains(par.echam.clim, 'rp000')
            file=dir(sprintf('/project2/tas1/ockham/data11/tas/echam-aiv_rcc_6.1.00p1/%s/ATM_%s_0020_39.nc', par.echam.clim, par.echam.clim));
        else
            file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/echam/ATM_rjg_%s_*.ymonmean.nc', par.echam.clim));
        end
        fullpath=sprintf('%s/%s', file.folder, file.name);
        temp = double(ncread(fullpath, var));
        var = 'geopoth';
        fullpath=sprintf('%s/%s', file.folder, file.name);
        zg = double(ncread(fullpath, var));
    elseif contains(type, 'echam_pl')
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s', type);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s', type);
    end

    load(sprintf('%s/grid.mat', prefix)); % read grid data
    load(sprintf('%s/srfc.mat', prefix)); % load surface data

    ma_si = nan([length(grid.dim3.lon), length(grid.dim3.lat), length(grid.dim3.si), 12]);

    pb = CmdLineProgressBar("Calculating moist adiabats...");
    for ilon = 1:length(grid.dim3.lon)
        pb.print(ilon, length(grid.dim3.lon)); % output progress of moist adiabat calculonion
        for ilat = 1:length(grid.dim3.lat);
            for imon = 1:12;
                for fn = fieldnames(srfc)'; srfc_var = fn{1};
                    if strcmp(par.ma_init, 'surf')
                        ma_in.(srfc_var) = srfc.(srfc_var)(ilon,ilat,imon);
                    else
                        if any(strcmp(type, {'era5', 'erai'}))
                            ma_in.sp = srfc.sp(ilon,ilat,imon);
                            ma_in.pinit = par.ma_init*ma_in.sp;
                            ma_in.t2m = interp1(grid.dim3.plev, squeeze(temp(ilon,ilat,:,imon)), ma_in.pinit);
                            ma_in.d2m = ma_in.t2m; % saturated
                        elseif strcmp(type, 'gcm')
                            ma_in.ps = srfc.ps(ilon,ilat,imon);
                            ma_in.pinit = par.ma_init*ma_in.ps;
                            ma_in.tas = interp1(grid.dim3.plev, squeeze(temp(ilon,ilat,:,imon)), ma_in.pinit);
                            ma_in.hurs = 100; % saturated
                        end
                        ma_in.zs = interp1(grid.dim3.plev, squeeze(zg(ilon,ilat,:,imon)), ma_in.pinit);
                    end
                end

                if any(strcmp(type, {'era5', 'erai', 'echam'}));
                    ma_si(ilon,ilat,:,imon) = calc_ma_dew_si(ma_in, grid.dim3.plev, par, type, grid);
                elseif strcmp(type, 'gcm')
                    ma_si(ilon,ilat,:,imon) = calc_ma_hurs_si(ma_in, grid.dim3.plev, par, grid);
                end

            end
        end
    end

    if any(strcmp(type, {'era5', 'erai'})); newdir=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s', type);
    elseif strcmp(type, 'merra2'); newdir=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s', type);
    elseif strcmp(type, 'gcm'); newdir=sprintf('/project2/tas1/miyawaki/projects/002/data/read/gcm/%s/%s', par.model, par.gcm.clim);
    elseif strcmp(type, 'echam'); newdir=sprintf('/project2/tas1/miyawaki/projects/002/data/read/echam/%s', par.echam.clim);
    elseif contains(type, 'echam_pl'); newdir=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s', type); end;
    if ~exist(newdir, 'dir'); mkdir(newdir); end
    if strcmp(par.ma_init, 'surf')
        filename=sprintf('ma_si_%s.mat', par.ma_init);
    else
        filename=sprintf('ma_si_%g.mat', par.ma_init);
    end
    save(sprintf('%s/%s', newdir, filename), 'ma_si', '-v7.3');
end
function proc_ma_mon_lat(type, par)
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
    elseif strcmp(type, 'echam_ml')
        foldername = sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/', type, par.lat_interp);
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s', type);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s', type);
    elseif strcmp(type, 'echam_pl')
        foldername = sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/', type, par.lat_interp);
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s', type);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s', type);
    end
    load(sprintf('%s/grid.mat', prefix)); % read grid data
    if strcmp(par.ma_init, 'surf')
        load(sprintf('%s/ma_si_%s.mat', prefix, par.ma_init)); masi_orig = ma_si; clear ma_si; % read temp in si coordinates
    else
        load(sprintf('%s/ma_si_%g.mat', prefix, par.ma_init)); masi_orig = ma_si; clear ma_si; % read temp in si coordinates
    end
    load(sprintf('%s/pa_si.mat', prefix)); pasi_orig = pa_si; clear pa_si; % read temp in si coordinates
    load(sprintf('%s/srfc.mat', prefix)); % load surface data
    load(sprintf('%s/%s/masks.mat', prefix_proc, par.lat_interp)); % load land and ocean masks

    if strcmp(par.lat_interp, 'std')
        lat = par.lat_std;
    else
        lat = grid.dim3.lat;
    end

    % interpolate ta to standard lat grid
    masi_orig = permute(masi_orig, [2 1 3 4]);
    masi_orig = interp1(grid.dim3.lat, masi_orig, lat);
    masi_orig = permute(masi_orig, [2 1 3 4]);

    pasi_orig = permute(pasi_orig, [2 1 3 4]);
    pasi_orig = interp1(grid.dim3.lat, pasi_orig, lat);
    pasi_orig = permute(pasi_orig, [2 1 3 4]);

    masi_sm.lo = masi_orig; % surface is already masked in standard sigma coordinates
    pasi_sm.lo = pasi_orig; % surface is already masked in spndard sigma coordinates

    masi_sm.lo = permute(masi_sm.lo, [1 2 4 3]); % bring plev to last dimension

    pasi_sm.lo = permute(pasi_sm.lo, [1 2 4 3]); % bring plev to last dimension

    mask_t.land = nanmean(mask.land, 3);
    mask_t.ocean = nanmean(mask.ocean, 3);

    for l = {'lo'}; land = l{1}; % over land, over ocean, or both
        masi.(land)= squeeze(nanmean(masi_sm.(land), 1)); % zonal average
        pasi.(land)= squeeze(nanmean(pasi_sm.(land), 1)); % zonal average
    end

    for l = {'lo'}; land = l{1}; % over land, over ocean, or both
        % take time averages
        for t = {'ann', 'djf', 'jja', 'mam', 'son'}; time = t{1};
            if strcmp(time, 'ann')
                masi_t.(land).(time) = squeeze(nanmean(masi_sm.(land), 3));
                pasi_t.(land).(time) = squeeze(nanmean(pasi_sm.(land), 3));
            elseif strcmp(time, 'djf')
                masi_shift.(land) = circshift(masi_sm.(land), 1, 3);
                pasi_shift.(land) = circshift(pasi_sm.(land), 1, 3);
                masi_t.(land).(time) = squeeze(nanmean(masi_shift.(land)(:,:,1:3,:), 3));
                pasi_t.(land).(time) = squeeze(nanmean(pasi_shift.(land)(:,:,1:3,:), 3));
            elseif strcmp(time, 'jja')
                masi_t.(land).(time) = squeeze(nanmean(masi_sm.(land)(:,:,6:8,:), 3));
                pasi_t.(land).(time) = squeeze(nanmean(pasi_sm.(land)(:,:,6:8,:), 3));
            elseif strcmp(time, 'mam')
                masi_t.(land).(time) = squeeze(nanmean(masi_sm.(land)(:,:,3:5,:), 3));
                pasi_t.(land).(time) = squeeze(nanmean(pasi_sm.(land)(:,:,3:5,:), 3));
            elseif strcmp(time, 'son')
                masi_t.(land).(time) = squeeze(nanmean(masi_sm.(land)(:,:,9:11,:), 3));
                pasi_t.(land).(time) = squeeze(nanmean(pasi_sm.(land)(:,:,9:11,:), 3));
            end
        end
    end

    % save filtered data
    printname = [foldername 'ma_mon_lat'];
    if ~exist(foldername, 'dir')
        mkdir(foldername)
    end
    if par.do_surf; save(printname, 'masi', 'lat');
    else save(printname, 'masi', 'pasi', 'lat', '-v7.3'); end

    printname = [foldername 'ma_lon_lat'];
    if ~exist(foldername, 'dir')
        mkdir(foldername)
    end
    if par.do_surf; save(printname, 'masi_t', 'lat');
    else save(printname, 'masi_t', 'pasi_t', 'lat', '-v7.3'); end
end % compute mon x lat temperature field
function proc_thetaeq_mon_lat(type, par)
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
    elseif strcmp(type, 'echam_ml')
        foldername = sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/', type, par.lat_interp);
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s', type);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s', type);
    elseif strcmp(type, 'echam_pl')
        foldername = sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/', type, par.lat_interp);
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s', type);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s', type);
    end
    load(sprintf('%s/grid.mat', prefix)); % read grid data
    load(sprintf('%s/thetaeq_si.mat', prefix)); thetaeqsi_orig = thetaeq_si.spl; clear thetaeq_si; % read temp in si coordinates
    load(sprintf('%s/pa_si.mat', prefix)); pasi_orig = pa_si; clear pa_si; % read temp in si coordinates
    load(sprintf('%s/srfc.mat', prefix)); % load surface data
    load(sprintf('%s/%s/masks.mat', prefix_proc, par.lat_interp)); % load land and ocean masks

    if strcmp(par.lat_interp, 'std')
        lat = par.lat_std;
    else
        lat = grid.dim3.lat;
    end

    % interpolate thetaeq to sthetaeqndard lat grid
    thetaeqsi_orig = permute(thetaeqsi_orig, [2 1 3 4]);
    thetaeqsi_orig = interp1(grid.dim3.lat, thetaeqsi_orig, lat);
    thetaeqsi_orig = permute(thetaeqsi_orig, [2 1 3 4]);

    pasi_orig = permute(pasi_orig, [2 1 3 4]);
    pasi_orig = interp1(grid.dim3.lat, pasi_orig, lat);
    pasi_orig = permute(pasi_orig, [2 1 3 4]);

    thetaeqsi_sm.lo = thetaeqsi_orig; % surface is already masked in sthetaeqndard sigma coordinates
    pasi_sm.lo = pasi_orig; % surface is already masked in spndard sigma coordinates

    thetaeqsi_sm.lo = permute(thetaeqsi_sm.lo, [1 2 4 3]); % bring plev to last dimension

    pasi_sm.lo = permute(pasi_sm.lo, [1 2 4 3]); % bring plev to last dimension

    mask_t.land = nanmean(mask.land, 3);
    mask_t.ocean = nanmean(mask.ocean, 3);

    for l = {'lo'}; land = l{1}; % over land, over ocean, or both
        thetaeqsi.(land)= squeeze(mean(thetaeqsi_sm.(land), 1)); % zonal average
        pasi.(land)= squeeze(mean(pasi_sm.(land), 1)); % zonal average
    end

    for l = {'lo'}; land = l{1}; % over land, over ocean, or both
        % thetaeqke time averages
        for t = {'ann', 'djf', 'jja', 'mam', 'son'}; time = t{1};
            if strcmp(time, 'ann')
                thetaeqsi_t.(land).(time) = squeeze(mean(thetaeqsi_sm.(land), 3));
                pasi_t.(land).(time) = squeeze(mean(pasi_sm.(land), 3));
            elseif strcmp(time, 'djf')
                thetaeqsi_shift.(land) = circshift(thetaeqsi_sm.(land), 1, 3);
                pasi_shift.(land) = circshift(pasi_sm.(land), 1, 3);
                thetaeqsi_t.(land).(time) = squeeze(mean(thetaeqsi_shift.(land)(:,:,1:3,:), 3));
                pasi_t.(land).(time) = squeeze(mean(pasi_shift.(land)(:,:,1:3,:), 3));
            elseif strcmp(time, 'jja')
                thetaeqsi_t.(land).(time) = squeeze(mean(thetaeqsi_sm.(land)(:,:,6:8,:), 3));
                pasi_t.(land).(time) = squeeze(mean(pasi_sm.(land)(:,:,6:8,:), 3));
            elseif strcmp(time, 'mam')
                thetaeqsi_t.(land).(time) = squeeze(mean(thetaeqsi_sm.(land)(:,:,3:5,:), 3));
                pasi_t.(land).(time) = squeeze(mean(pasi_sm.(land)(:,:,3:5,:), 3));
            elseif strcmp(time, 'son')
                thetaeqsi_t.(land).(time) = squeeze(mean(thetaeqsi_sm.(land)(:,:,9:11,:), 3));
                pasi_t.(land).(time) = squeeze(mean(pasi_sm.(land)(:,:,9:11,:), 3));
            end
        end
    end

    % save filtered dathetaeq
    printname = [foldername 'thetaeq_mon_lat'];
    if ~exist(foldername, 'dir')
        mkdir(foldername)
    end
    if par.do_surf; save(printname, 'thetaeqsi', 'lat');
    else save(printname, 'thetaeqsi', 'pasi', 'lat', '-v7.3'); end

    printname = [foldername 'thetaeq_lon_lat'];
    if ~exist(foldername, 'dir')
        mkdir(foldername)
    end
    if par.do_surf; save(printname, 'thetaeqsi_t', 'lat');
    else save(printname, 'thetaeqsi_t', 'pasi_t', 'lat', '-v7.3'); end
end % compute mon x lat temperature field
function proc_dtdz_mon_lat(type, par)
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
    load(sprintf('%s/dtdzsi.mat', prefix));
    load(sprintf('%s/%s/masks.mat', prefix_proc, par.lat_interp)); % load land and ocean masks

    % Land/ocean filter 3D variables
    mask.land_vert = repmat(mask.land, [1 1 1 size(dtdzsi, 3)]); % expand land mask to vertical dim
    mask.land_vert = permute(mask.land, [1 2 4 3]); % place vertical dim where it belongs
    mask.ocean_vert = repmat(mask.ocean, [1 1 1 size(dtdzsi, 3)]); % expand ocean mask to vertical dim
    mask.ocean_vert = permute(mask.ocean, [1 2 4 3]); % place vertical dim where it belongs

    dtdz_sm.lo = dtdzsi;
    dtdz_sm.l = dtdz_sm.lo.*mask.ocean_vert; % filter dtdz with surface mask
    dtdz_sm.o = dtdz_sm.lo.*mask.land_vert; % filter dtdz with surface mask
    dtdz_sm.lo = permute(dtdz_sm.lo, [1 2 4 3]); % bring plev to last dimension
    dtdz_sm.l = permute(dtdz_sm.l, [1 2 4 3]); % bring plev to last dimension
    dtdz_sm.o = permute(dtdz_sm.o, [1 2 4 3]); % bring plev to last dimension

    if strcmp(par.lat_interp, 'std')
        lat = par.lat_std;
    else
        lat = grid.dim3.lat;
    end

    pb = CmdLineProgressBar("Calculating lapse rates...");
    for l = {'lo', 'l', 'o'}; land = l{1};
        for ilat = 1:length(lat);
            pb.print(ilat, length(lat));
            for imon = 1:12;
                dtdz.(land) = squeeze(nanmean(dtdz_sm.(land))); % zonal mean
            end

            for ilon = 1:length(grid.dim3.lon);
                for t = {'ann', 'djf', 'jja', 'mam', 'son'}; time = t{1};
                    if strcmp(time, 'ann')
                        dtdz_t.(land).(time) = squeeze(nanmean(dtdz_sm.(land), 3));
                    elseif strcmp(time, 'djf')
                        dtdz_shift.(land) = circshift(dtdz_sm.(land), 1, 3);
                        dtdz_t.(land).(time) = squeeze(nanmean(dtdz_shift.(land)(:,:,1:3,:), 3));
                    elseif strcmp(time, 'jja')
                        dtdz_t.(land).(time) = squeeze(nanmean(dtdz_sm.(land)(:,:,6:8,:), 3));
                    elseif strcmp(time, 'mam')
                        dtdz_t.(land).(time) = squeeze(nanmean(dtdz_sm.(land)(:,:,3:5,:), 3));
                    elseif strcmp(time, 'son')
                        dtdz_t.(land).(time) = squeeze(nanmean(dtdz_sm.(land)(:,:,9:11,:), 3));
                    end
                end
            end
        end
    end

    % save data into mat file
    printname = [foldername 'dtdz_mon_lat.mat'];
    if ~exist(foldername, 'dir')
        mkdir(foldername)
    end
    save(printname, 'dtdz');

    printname = [foldername 'dtdz_lon_lat.mat'];
    if ~exist(foldername, 'dir')
        mkdir(foldername)
    end
    save(printname, 'dtdz_t');

end % compute mon x lat moist adiabat field
function proc_malr_mon_lat(type, par)
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
    load(sprintf('%s/%s/ta_mon_lat.mat', prefix_proc, par.lat_interp)); % load taerature data
    load(sprintf('%s/%s/ta_lon_lat.mat', prefix_proc, par.lat_interp)); % load taerature data
    load(sprintf('%s/%s/masks.mat', prefix_proc, par.lat_interp)); % load land and ocean masks

    % calculate MALR following equations 3-5 in Stone and Carlson (1979)
    dalr = 9.8; % K/km
    eps = 0.622;
    R = 0.287; % J/g/K
    cp = 1.005; % J/g/K
    es0 = 6.11; % mb
    T0 = 273; % K

    if strcmp(par.lat_interp, 'std')
        lat = par.lat_std;
    else
        lat = grid.dim3.lat;
    end

    pb = CmdLineProgressBar("Calculating moist adiabatic lapse rate...");
    for l = {'lo', 'l', 'o'}; land = l{1};
        for ilat = 1:length(lat);
            pb.print(ilat, length(lat)); % output progress of moist adiabat calculation
            for imon = 1:12;
                tai = ta.(land);
                pa = 1e-2*grid.dim3.plev; pa = repmat(pa, [1 size(tai, 1) size(tai, 2)]); pa = permute(pa, [2 3 1]);
                L = 2510 - 2.38*(tai-T0);
                es = es0*exp(eps*L/R.*(1/T0-1./tai));
                desdT = eps*L.*es./(R*tai.^2);
                dtmdz.(land) = dalr * (1+eps*L.*es./(pa.*R.*tai))./(1+(eps.*L./(cp.*pa)).*(desdT));
                clear tai pa L es desdT;

                tai = taz.(land);
                pa = 1e-2*paz.(land);
                L = 2510 - 2.38*(tai-T0);
                es = es0*exp(eps*L/R.*(1/T0-1./tai));
                desdT = eps*L.*es./(R*tai.^2);
                dtmdz_z.(land) = dalr * (1+eps*L.*es./(pa.*R.*tai))./(1+(eps.*L./(cp.*pa)).*(desdT));
                clear tai pa L es desdT;

                tai = tasi.(land);
                pa = 1e-2*pasi.(land);
                L = 2510 - 2.38*(tai-T0);
                es = es0*exp(eps*L/R.*(1/T0-1./tai));
                desdT = eps*L.*es./(R*tai.^2);
                dtmdz_si.(land) = dalr * (1+eps*L.*es./(pa.*R.*tai))./(1+(eps.*L./(cp.*pa)).*(desdT));
                clear tai pa L es desdT;

            end
            for ilon = 1:length(grid.dim3.lon);
                for t = {'ann', 'djf', 'jja', 'mam', 'son'}; time = t{1};

                    tai = ta_t.(land).(time);
                    pa = 1e-2*grid.dim3.plev; pa = repmat(pa, [1 size(tai, 1) size(tai, 2)]); pa = permute(pa, [2 3 1]);
                    L = 2510 - 2.38*(tai-T0);
                    es = es0*exp(eps*L/R.*(1/T0-1./tai));
                    desdT = eps*L.*es./(R*tai.^2);
                    dtmdz_t.(land).(time) = dalr * (1+eps*L.*es./(pa.*R.*tai))./(1+(eps.*L./(cp.*pa)).*(desdT));
                    clear tai pa L es desdT;

                    tai = taz_t.(land).(time);
                    pa = 1e-2*paz_t.(land).(time);
                    L = 2510 - 2.38*(tai-T0);
                    es = es0*exp(eps*L/R.*(1/T0-1./tai));
                    desdT = eps*L.*es./(R*tai.^2);
                    dtmdz_z_t.(land).(time) = dalr * (1+eps*L.*es./(pa.*R.*tai))./(1+(eps.*L./(cp.*pa)).*(desdT));
                    clear tai pa L es desdT;

                    tai = tasi_t.(land).(time);
                    pa = 1e-2*pasi_t.(land).(time);
                    L = 2510 - 2.38*(tai-T0);
                    es = es0*exp(eps*L/R.*(1/T0-1./tai));
                    desdT = eps*L.*es./(R*tai.^2);
                    dtmdz_si_t.(land).(time) = dalr * (1+eps*L.*es./(pa.*R.*tai))./(1+(eps.*L./(cp.*pa)).*(desdT));
                    clear tai pa L es desdT;

                end
            end
        end
    end

    % save data into mat file
    printname = [foldername 'malr_mon_lat.mat'];
    if ~exist(foldername, 'dir')
        mkdir(foldername)
    end
    save(printname, 'dtmdz', 'dtmdz_z', 'dtmdz_si');

    printname = [foldername 'malr_lon_lat.mat'];
    if ~exist(foldername, 'dir')
        mkdir(foldername)
    end
    save(printname, 'dtmdz_t', 'dtmdz_z_t', 'dtmdz_si_t');

end % compute mon x lat moist adiabat field
function make_tai(type, par)
    % add surface data to temperature and interpolate to hi-res grid
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
    elseif strcmp(type, 'merra2')
        foldername = sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/', type, par.lat_interp);
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s', type);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s', type);
        file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/%s/temp/%s_temp_%s.ymonmean.nc', type, type, par.(type).yr_span));
        fullpath=sprintf('%s/%s', file.folder, file.name);
        temp = ncread(fullpath, 'T');
        file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/%s/zg/%s_zg_%s.ymonmean.nc', type, type, par.(type).yr_span));
        fullpath=sprintf('%s/%s', file.folder, file.name);
        zg = ncread(fullpath, 'H');
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
        if contains(par.echam.clim, 'rp000')
            file=dir(sprintf('/project2/tas1/ockham/data11/tas/echam-aiv_rcc_6.1.00p1/%s/ATM_%s_0020_39.nc', par.echam.clim, par.echam.clim));
        else
            file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/echam/ATM_rjg_%s_*.ymonmean.nc', par.echam.clim));
        end
        fullpath=sprintf('%s/%s', file.folder, file.name);
        temp = double(ncread(fullpath, var));
        var = 'geopoth';
        fullpath=sprintf('%s/%s', file.folder, file.name);
        zg = double(ncread(fullpath, var));
    elseif any(strcmp(type, {'echam_ml', 'echam_pl'}))
        foldername = sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/', type, par.lat_interp);
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s', type);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s', type);
        var = 't';
        file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/%s/ATM_*.ymonmean.nc', type));
        fullpath=sprintf('%s/%s', file.folder, file.name);
        temp = double(ncread(fullpath, var));
        var = 'geopoth';
        file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/%s/ATM_*.ymonmean.nc', type));
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
    elseif strcmp(type, 'merra2')
        ps_vert = repmat(srfc.PS, [1 1 1 size(temp, 3)]); % dims (lon x lat x time x plev)
        ps_vert = permute(ps_vert, [1 2 4 3]); % dims (lon x lat x plev x time)
        pa = double(permute(repmat(grid.dim3.plev, [1 size(srfc.PS)]), [2 3 1 4]));
        ts_vert = repmat(srfc.T2M, [1 1 1 size(temp, 3)]); % dims (lon x lat x time x plev)
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
    elseif any(strcmp(type, {'echam', 'echam_pl'}))
        ps_vert = repmat(srfc.aps, [1 1 1 size(temp, 3)]); % dims (lon x lat x time x plev)
        ps_vert = permute(ps_vert, [1 2 4 3]); % dims (lon x lat x plev x time)
        pa = permute(repmat(grid.dim3.plev, [1 size(srfc.aps)]), [2 3 1 4]);
        ts_vert = repmat(srfc.temp2, [1 1 1 size(temp, 3)]); % dims (lon x lat x time x plev)
        ts_vert = permute(ts_vert, [1 2 4 3]); % dims (lon x lat x plev x time)
        zs_vert = repmat(srfc.zs, [1 1 1 size(temp, 3)]); % dims (lon x lat x time x plev)
        zs_vert = permute(zs_vert, [1 2 4 3]); % dims (lon x lat x plev x time)
    elseif strcmp(type, 'echam_ml')
        ps_vert = repmat(srfc.aps, [1 1 1 size(temp, 3)]); % dims (lon x lat x time x plev)
        ps_vert = permute(ps_vert, [1 2 4 3]); % dims (lon x lat x plev x time)
        a = permute(repmat(grid.dim3.a, [1 size(srfc.aps)]), [2 3 1 4]);
        b = permute(repmat(grid.dim3.b, [1 size(srfc.aps)]), [2 3 1 4]);
        pa = a+b.*ps_vert;
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

                tai_sm.lo(:,ilon,ilat,time) = interp1(tapanan, tanan, par.pa, 'spline', nan); % interpolate to higher resolution vertical grid
                zgi_sm.lo(:,ilon,ilat,time) = interp1(zgpanan, zgnan, par.pa, 'spline', nan); % interpolate to higher resolution vertical grid
            end
        end
    end
    clear pa_plus ta_plus zg_plus; % clear unneeded variables
    tai_sm.lo = permute(tai_sm.lo, [2 3 1 4]); % bring plev back to third dimension
    tai = tai_sm.lo;
    zgi_sm.lo = permute(zgi_sm.lo, [2 3 1 4]); % bring plev back to third dimension
    zgi = zgi_sm.lo;

    if any(strcmp(type, {'era5', 'erai'})); newdir=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s', type);
    elseif strcmp(type, 'merra2'); newdir=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s', type);
    elseif strcmp(type, 'gcm'); newdir=sprintf('/project2/tas1/miyawaki/projects/002/data/read/gcm/%s/%s', par.model, par.gcm.clim);
    elseif strcmp(type, 'echam'); newdir=sprintf('/project2/tas1/miyawaki/projects/002/data/read/echam/%s', par.echam.clim);
    elseif any(strcmp(type, {'echam_ml', 'echam_pl'})); newdir=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s', type); end;
    if ~exist(newdir, 'dir'); mkdir(newdir); end
    filename='tai.mat';
    save(sprintf('%s/%s', newdir, filename), 'tai', 'zgi', '-v7.3');

end
function make_dtdzsi(type, par)
% compute model lapse rate in lat x plev x mon (add surface data after converting to sigma)
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
    elseif strcmp(type, 'merra2')
        foldername = sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/', type, par.lat_interp);
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s', type);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s', type);
        file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/%s/temp/%s_temp_%s.ymonmean.nc', type, type, par.(type).yr_span));
        fullpath=sprintf('%s/%s', file.folder, file.name);
        temp = ncread(fullpath, 'T');
        file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/%s/zg/%s_zg_%s.ymonmean.nc', type, type, par.(type).yr_span));
        fullpath=sprintf('%s/%s', file.folder, file.name);
        zg = ncread(fullpath, 'H');
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
        if contains(par.echam.clim, 'rp000')
            file=dir(sprintf('/project2/tas1/ockham/data11/tas/echam-aiv_rcc_6.1.00p1/%s/ATM_%s_0020_39.nc', par.echam.clim, par.echam.clim));
        else
            file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/echam/ATM_rjg_%s_*.ymonmean.nc', par.echam.clim));
        end
        fullpath=sprintf('%s/%s', file.folder, file.name);
        temp = double(ncread(fullpath, var));
        var = 'geopoth';
        fullpath=sprintf('%s/%s', file.folder, file.name);
        zg = double(ncread(fullpath, var));
    elseif any(strcmp(type, {'echam_ml', 'echam_pl'}))
        foldername = sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/', type, par.lat_interp);
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s', type);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s', type);
        var = 't';
        file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/%s/ATM_*.ymonmean.nc', type));
        fullpath=sprintf('%s/%s', file.folder, file.name);
        temp = double(ncread(fullpath, var));
        var = 'geopoth';
        file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/%s/ATM_*.ymonmean.nc', type));
        fullpath=sprintf('%s/%s', file.folder, file.name);
        zg = double(ncread(fullpath, var));
    end
    load(sprintf('%s/grid.mat', prefix)); % read grid data
    load(sprintf('%s/srfc.mat', prefix)); % read surface variable data
    % load(sprintf('%s/orog.mat', prefix)); % read surface orography

    coarse_si = 1e-5*grid.dim3.plev;
    [~,si0] = max(coarse_si); % surface sigma index

    % create surface mask
    if strcmp(type, 'era5') | strcmp(type, 'erai')
        ps_vert = repmat(srfc.sp, [1 1 1 size(temp, 3)]); % dims (lon x lat x time x plev)
        ps_vert = permute(ps_vert, [1 2 4 3]); % dims (lon x lat x plev x time)
        pa = double(permute(repmat(grid.dim3.plev, [1 size(srfc.sp)]), [2 3 1 4]));
    elseif strcmp(type, 'merra2')
        ps_vert = repmat(srfc.PS, [1 1 1 size(temp, 3)]); % dims (lon x lat x time x plev)
        ps_vert = permute(ps_vert, [1 2 4 3]); % dims (lon x lat x plev x time)
        pa = double(permute(repmat(grid.dim3.plev, [1 size(srfc.PS)]), [2 3 1 4]));
    elseif strcmp(type, 'gcm')
        ps_vert = repmat(srfc.ps, [1 1 1 size(temp, 3)]); % dims (lon x lat x time x plev)
        ps_vert = permute(ps_vert, [1 2 4 3]); % dims (lon x lat x plev x time)
        pa = permute(repmat(grid.dim3.plev, [1 size(srfc.ps)]), [2 3 1 4]);
    elseif any(strcmp(type, {'echam', 'echam_pl'}))
        ps_vert = repmat(srfc.aps, [1 1 1 size(temp, 3)]); % dims (lon x lat x time x plev)
        ps_vert = permute(ps_vert, [1 2 4 3]); % dims (lon x lat x plev x time)
        pa = permute(repmat(grid.dim3.plev, [1 size(srfc.aps)]), [2 3 1 4]);
    elseif strcmp(type, 'echam_ml')
        ps_vert = repmat(srfc.aps, [1 1 1 size(temp, 3)]); % dims (lon x lat x time x plev)
        ps_vert = permute(ps_vert, [1 2 4 3]); % dims (lon x lat x plev x time)
        a = permute(repmat(grid.dim3.a, [1 size(srfc.aps)]), [2 3 1 4]);
        b = permute(repmat(grid.dim3.b, [1 size(srfc.aps)]), [2 3 1 4]);
        pa = a+b.*ps_vert;
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

                ta_tmp = interp1(pa(:,lo,la,mo)./ps_vert(1,lo,la,mo), ta_sm(:,lo,la,mo), coarse_si, 'linear', nan);
                zg_tmp = interp1(pa(:,lo,la,mo)./ps_vert(1,lo,la,mo), zg_sm(:,lo,la,mo), coarse_si, 'linear', nan);

                % add surface data
                if strcmp(type, 'era5') | strcmp(type, 'erai')
                    ta_tmp(si0) = srfc.t2m(lo,la,mo);
                elseif strcmp(type, 'merra2')
                    ta_tmp(si0) = srfc.T2M(lo,la,mo);
                elseif strcmp(type, 'gcm')
                    ta_tmp(si0) = srfc.tas(lo,la,mo);
                elseif contains(type, 'echam')
                    ta_tmp(si0) = srfc.temp2(lo,la,mo);
                end
                zg_tmp(si0) = srfc.zs(lo,la,mo);

                % only keep nonnan data and redo interpolation
                notnan = find(~isnan(squeeze(ta_tmp)) & ~isnan(squeeze(zg_tmp)));
                ta_si(:,lo,la,mo) = interp1(coarse_si(notnan), ta_tmp(notnan), grid.dim3.si, 'spline', nan);
                zg_si(:,lo,la,mo) = interp1(coarse_si(notnan), zg_tmp(notnan), grid.dim3.si, 'spline', nan);

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
    elseif strcmp(type, 'merra2'); newdir=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s', type);
    elseif strcmp(type, 'gcm'); newdir=sprintf('/project2/tas1/miyawaki/projects/002/data/read/gcm/%s/%s', par.model, par.gcm.clim);
    elseif strcmp(type, 'echam'); newdir=sprintf('/project2/tas1/miyawaki/projects/002/data/read/echam/%s', par.echam.clim);
    elseif any(strcmp(type, {'echam_ml', 'echam_pl'})); newdir=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s', type); end;
    if ~exist(newdir, 'dir'); mkdir(newdir); end
    filename='dtdzsi.mat';
    save(sprintf('%s/%s', newdir, filename), 'dtdzsi', '-v7.3');

end
function make_malrsi(type, par)
% compute moist adiabatic lapse rate at every lat, time
% calculate moist adiabats
    if strcmp(type, 'era5') | strcmp(type, 'erai')
        foldername = sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/', type, par.lat_interp);
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s', type);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s', type);
    elseif strcmp(type, 'merra2')
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
    elseif any(strcmp(type, {'echam_ml', 'echam_pl'}))
        foldername = sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/', type, par.lat_interp);
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s', type);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s', type);
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
    if any(strcmp(type, {'era5', 'erai'}))
        tas = srfc.t2m;
        ps = srfc.sp;
    elseif strcmp(type, 'merra2')
        tas = srfc.T2M;
        ps = srfc.PS;
    elseif strcmp(type, 'gcm')
        tas = srfc.tas;
        ps = srfc.ps;
    elseif contains(type, 'echam')
        tas = srfc.temp2;
        ps = srfc.aps;
    end
    L = 2510 - 2.38*(tas-T0);
    es = es0*exp(eps*L/R.*(1/T0-1./tas));
    desdT = eps*L.*es./(R*tas.^2);
    dtasmdz = dalr * (1+eps*L.*es./(1e-2*ps.*R.*tas))./(1+(eps.*L./(cp.*1e-2*ps)).*(desdT));

    if strcmp(type, 'era5') | strcmp(type, 'erai')
        ps_vert = repmat(srfc.sp, [1 1 1 size(dtmdz, 3)]); % dims (lon x lat x time x plev)
        ps_vert = permute(ps_vert, [1 2 4 3]); % dims (lon x lat x plev x time)
        pa = double(permute(repmat(par.pa', [1 size(srfc.sp)]), [2 3 1 4]));
    elseif strcmp(type, 'merra2')
        ps_vert = repmat(srfc.PS, [1 1 1 size(dtmdz, 3)]); % dims (lon x lat x time x plev)
        ps_vert = permute(ps_vert, [1 2 4 3]); % dims (lon x lat x plev x time)
        pa = double(permute(repmat(par.pa', [1 size(srfc.PS)]), [2 3 1 4]));
    elseif strcmp(type, 'gcm')
        ps_vert = repmat(srfc.ps, [1 1 1 size(dtmdz, 3)]); % dims (lon x lat x time x plev)
        ps_vert = permute(ps_vert, [1 2 4 3]); % dims (lon x lat x plev x time)
        pa = permute(repmat(par.pa', [1 size(srfc.ps)]), [2 3 1 4]);
    elseif contains(type, 'echam')
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
    elseif strcmp(type, 'merra2'); newdir=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s', type);
    elseif strcmp(type, 'gcm'); newdir=sprintf('/project2/tas1/miyawaki/projects/002/data/read/gcm/%s/%s', par.model, par.gcm.clim);
    elseif strcmp(type, 'echam'); newdir=sprintf('/project2/tas1/miyawaki/projects/002/data/read/echam/%s', par.echam.clim);
    elseif any(strcmp(type, {'echam_ml', 'echam_pl'})); newdir=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s', type); end;
    if ~exist(newdir, 'dir'); mkdir(newdir); end
    filename='malrsi.mat';
    save(sprintf('%s/%s', newdir, filename), 'dtmdzsi', '-v7.3');

end

% SI BL PROC
function proc_ga_malr_diff_si_mon_lat(type, par)
    if strcmp(type, 'era5') | strcmp(type, 'erai')
        foldername = sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/si_bl_%g/', type, par.lat_interp, par.si_bl);
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s', type);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s', type);
        file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/%s/temp/%s_temp_%s.ymonmean.nc', type, type, par.(type).yr_span));
        fullpath=sprintf('%s/%s', file.folder, file.name);
        ta_orig = ncread(fullpath, 't');
    elseif strcmp(type, 'merra2')
        foldername = sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/si_bl_%g/', type, par.lat_interp, par.si_bl);
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s', type);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s', type);
        file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/%s/temp/%s_temp_%s.ymonmean.nc', type, type, par.(type).yr_span));
        fullpath=sprintf('%s/%s', file.folder, file.name);
        ta_orig = ncread(fullpath, 'T');
    elseif strcmp(type, 'gcm')
        foldername = sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s/%s/si_bl_%g/', type, par.model, par.gcm.clim, par.lat_interp, par.si_bl);
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s/%s', type, par.model, par.gcm.clim);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s', type, par.model, par.gcm.clim);
        var = 'ta';
        file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/gcm/%s/%s_Amon_%s_%s_r1i1p1_*.ymonmean.nc', par.model, var, par.model, par.gcm.clim));
        fullpath=sprintf('%s/%s', file.folder, file.name);
        ta_orig = ncread(fullpath, var);
    elseif strcmp(type, 'echam')
        foldername = sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s/si_bl_%g/', type, par.echam.clim, par.lat_interp, par.si_bl);
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s', type, par.echam.clim);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s', type, par.echam.clim);
        var = 't';
        if contains(par.echam.clim, 'rp000')
            file=dir(sprintf('/project2/tas1/ockham/data11/tas/echam-aiv_rcc_6.1.00p1/%s/ATM_%s_0020_39.nc', par.echam.clim, par.echam.clim));
        else
            file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/echam/ATM_rjg_%s_*.ymonmean.nc', par.echam.clim));
        end
        fullpath=sprintf('%s/%s', file.folder, file.name);
        ta_orig = double(ncread(fullpath, var));
    elseif any(strcmp(type, {'echam_ml', 'echam_pl'}))
        foldername = sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/si_bl_%g/', type, par.lat_interp, par.si_bl);
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s', type);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s', type);
        var = 't';
        file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/%s/ATM_*.ymonmean.nc', type));
        fullpath=sprintf('%s/%s', file.folder, file.name);
        ta_orig = double(ncread(fullpath, var));
    end
    load(sprintf('%s/grid.mat', prefix)); % read grid data
    % load(sprintf('%s/dtdzzsi.mat', prefix)); % read temp in si coordinates
    load(sprintf('%s/dtdzsi.mat', prefix)); dtdzzsi = dtdzsi; clear dtdzsi; % read temp in si coordinates
    % load(sprintf('%s/malrzsi.mat', prefix)); % read temp in si coordinates
    load(sprintf('%s/malrsi.mat', prefix)); dtmdzzsi = dtmdzsi; clear dtmdzsi; % read temp in si coordinates
    load(sprintf('%s/%s/masks.mat', prefix_proc, par.lat_interp)); % load land and ocean masks

    if strcmp(par.lat_interp, 'std')
        lat = par.lat_std;
    else
        lat = grid.dim3.lat;
    end

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
function proc_ga_dalr_bl_diff_si_mon_lat(type, par)
    if strcmp(type, 'era5') | strcmp(type, 'erai')
        foldername = sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/si_bl_%g/', type, par.lat_interp, par.si_bl);
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s', type);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s', type);
        file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/%s/temp/%s_temp_%s.ymonmean.nc', type, type, par.(type).yr_span));
        fullpath=sprintf('%s/%s', file.folder, file.name);
        ta_orig = ncread(fullpath, 't');
    elseif strcmp(type, 'merra2')
        foldername = sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/si_bl_%g/', type, par.lat_interp, par.si_bl);
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s', type);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s', type);
        file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/%s/temp/%s_temp_%s.ymonmean.nc', type, type, par.(type).yr_span));
        fullpath=sprintf('%s/%s', file.folder, file.name);
        ta_orig = ncread(fullpath, 'T');
    elseif strcmp(type, 'gcm')
        foldername = sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s/%s/si_bl_%g/', type, par.model, par.gcm.clim, par.lat_interp, par.si_bl);
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s/%s', type, par.model, par.gcm.clim);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s', type, par.model, par.gcm.clim);
        var = 'ta';
        file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/gcm/%s/%s_Amon_%s_%s_r1i1p1_*.ymonmean.nc', par.model, var, par.model, par.gcm.clim));
        fullpath=sprintf('%s/%s', file.folder, file.name);
        ta_orig = ncread(fullpath, var);
    elseif strcmp(type, 'echam')
        foldername = sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s/si_bl_%g/', type, par.echam.clim, par.lat_interp, par.si_bl);
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s', type, par.echam.clim);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s', type, par.echam.clim);
        var = 't';
        if contains(par.echam.clim, 'rp000')
            file=dir(sprintf('/project2/tas1/ockham/data11/tas/echam-aiv_rcc_6.1.00p1/%s/ATM_%s_0020_39.nc', par.echam.clim, par.echam.clim));
        else
            file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/echam/ATM_rjg_%s_*.ymonmean.nc', par.echam.clim));
        end
        fullpath=sprintf('%s/%s', file.folder, file.name);
        ta_orig = double(ncread(fullpath, var));
    elseif any(strcmp(type, {'echam_ml', 'echam_pl'}))
        foldername = sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/si_bl_%g/', type, par.lat_interp, par.si_bl);
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s', type);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s', type);
        var = 't';
        file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/%s/ATM_*.ymonmean.nc', type));
        fullpath=sprintf('%s/%s', file.folder, file.name);
        ta_orig = double(ncread(fullpath, var));
    end
    load(sprintf('%s/grid.mat', prefix)); % read grid data
    load(sprintf('%s/dtdzsi.mat', prefix)); dtdzzsi = dtdzsi; clear dtdzsi; % read temp in si coordinates
    % load(sprintf('%s/dtdzzsi.mat', prefix)); % read temp in si coordinates
    load(sprintf('%s/%s/masks.mat', prefix_proc, par.lat_interp)); % load land and ocean masks

    dtmdzzsi = 1e3*par.g/par.cpd*ones(size(dtdzzsi));

    if strcmp(par.lat_interp, 'std')
        lat = par.lat_std;
    else
        lat = grid.dim3.lat;
    end

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
        if contains(par.echam.clim, 'rp000')
            file=dir(sprintf('/project2/tas1/ockham/data11/tas/echam-aiv_rcc_6.1.00p1/%s/ATM_%s_0020_39.nc', par.echam.clim, par.echam.clim));
        else
            file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/echam/ATM_rjg_%s_*.ymonmean.nc', par.echam.clim));
        end
        fullpath=sprintf('%s/%s', file.folder, file.name);
        ta_orig = double(ncread(fullpath, var));
    end
    load(sprintf('%s/grid.mat', prefix)); % read grid data
    load(sprintf('%s/dtdzzsi.mat', prefix)); % read temp in si coordinates
    load(sprintf('%s/malrzsi.mat', prefix)); % read temp in si coordinates
    load(sprintf('%s/%s/masks.mat', prefix_proc, par.lat_interp)); % load land and ocean masks

    if strcmp(par.lat_interp, 'std')
        lat = par.lat_std;
    else
        lat = grid.dim3.lat;
    end

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

% EP PROC
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
            % for l = {'lo', 'l', 'o'}; land = l{1};
            for l = {'lo'}; land = l{1};
                % lon x lat structure over various time averages
                % for t = {'ann', 'djf', 'jja', 'mam', 'son'}; time = t{1};
                for t = {'ann'}; time = t{1};
                    rcae_alt_t.(land).(time) = def_rcae_alt(type, flux_t.(land).(time), vh.(land).(time), par);
                    rcae_alt_rc_t.(land).(time) = def_rcae_alt_recomp_r1(type, flux_t.(land).(time), vh.(land).(time), par);
                end % time
            end % land

            printname = [foldername 'rcae_alt_t.mat'];
            printname_rc = [foldername 'rcae_alt_rc_t.mat'];
        elseif strcmp(ftype, 'flux_z') % show lat x mon structure of RCAE_ALT
            % for l = {'lo', 'l', 'o'}; land = l{1};
            for l = {'lo'}; land = l{1};
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
function proc_ta_si(type, par)
    if strcmp(type, 'era5') | strcmp(type, 'erai')
        foldername = sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/eps_%g_ga_%g/', type, par.lat_interp, par.ep, par.ga);
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s', type);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s', type);
    elseif strcmp(type, 'merra2')
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
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s', type, par.echam.clim);
    end
    load(sprintf('%s/grid.mat', prefix)); % read ERA5 grid data
    load(sprintf('%s/srfc.mat', prefix)); % load surface data
    ta_orig = load(sprintf('%s/ta_si.mat', prefix)); ta_si=ta_orig.ta_si.spl; clear ta_orig; % load temperature in sigma
    load(sprintf('%s/%s/masks.mat', prefix_proc, par.lat_interp)); % load land and ocean masks
    load(sprintf('%s/%s/eps_%g_ga_%g/rcae_alt_t.mat', prefix_proc, par.lat_interp, par.ep, par.ga)); % load rcae data

    % mask.land_vert = repmat(mask.land, [1 1 1 size(ta_si, 3)]); % expand land mask to vertical dim
    % mask.land_vert = permute(mask.land, [1 2 4 3]); % place vertical dim where it belongs
    % mask.ocean_vert = repmat(mask.ocean, [1 1 1 size(ta_si, 3)]); % expand ocean mask to vertical dim
    % mask.ocean_vert = permute(mask.ocean, [1 2 4 3]); % place vertical dim where it belongs

    % ta_si_sm.l = ta_si.*mask.ocean_vert; % filter ta_si with surface mask
    % ta_si_sm.o = ta_si.*mask.land_vert; % filter ta_si with surface mask
    ta_si_sm.lo = permute(ta_si, [1 2 4 3]); % bring plev to last dimension
    % ta_si_sm.l = permute(ta_si_sm.l, [1 2 4 3]); % bring plev to last dimension
    % ta_si_sm.o = permute(ta_si_sm.o, [1 2 4 3]); % bring plev to last dimension

    clear ta_si;

    % mask_t.land = nanmean(mask.land, 3);
    % mask_t.ocean = nanmean(mask.ocean, 3);
    if any(strcmp(type, {'era5', 'erai'})); f_vec = par.era.fw;
    elseif any(strcmp(type, 'merra2')); f_vec = par.merra2.fw;
    elseif strcmp(type, 'gcm'); f_vec = par.gcm.fw;
    elseif strcmp(type, 'echam'); f_vec = par.echam.fw; end
    for f = f_vec; fw = f{1};
        for c = fieldnames(rcae_alt_t.lo.ann.(fw))'; crit = c{1};
            % for l = {'lo', 'l', 'o'}; land = l{1}; % over land, over ocean, or both
            for l = {'lo'}; land = l{1}; % over land, over ocean, or both
                % for t = {'ann', 'djf', 'jja', 'mam', 'son'}; time = t{1};
                for t = {'ann'}; time = t{1};
                    for re = {'rce', 'rae', 'rcae'}; regime = re{1};
                        if strcmp(time, 'ann')
                            ta_si_n.(land).(time) = squeeze(nanmean(ta_si_sm.(land), 3));
                        elseif strcmp(time, 'djf')
                            ta_si_shift = circshift(ta_si_sm.(land), 1, 3);
                            ta_si_n.(land).(time) = squeeze(nanmean(ta_si_shift(:,:,1:3,:), 3));
                        elseif strcmp(time, 'jja')
                            ta_si_n.(land).(time) = squeeze(nanmean(ta_si_sm.(land)(:,:,6:8,:), 3));
                        elseif strcmp(time, 'mam')
                            ta_si_n.(land).(time) = squeeze(nanmean(ta_si_sm.(land)(:,:,3:5,:), 3));
                        elseif strcmp(time, 'son')
                            ta_si_n.(land).(time) = squeeze(nanmean(ta_si_sm.(land)(:,:,9:11,:), 3));
                        end

                        filt.(land).(time).(fw).(crit).(regime) = nan(size(rcae_alt_t.(land).(time).(fw).(crit))); % create empty arrays to store filtering array
                        if strcmp(regime, 'rce'); filt.(land).(time).(fw).(crit).(regime)(rcae_alt_t.(land).(time).(fw).(crit)==1)=1; % set RCE=1, elsewhere nan
                        elseif strcmp(regime, 'rae'); filt.(land).(time).(fw).(crit).(regime)(rcae_alt_t.(land).(time).(fw).(crit)==-1)=1; % set RAE=1, elsewhere nan
                        elseif strcmp(regime, 'rcae'); filt.(land).(time).(fw).(crit).(regime)(rcae_alt_t.(land).(time).(fw).(crit)==0)=1; % set RCAE=1, elsewhere nan
                        end

                        ta_si_t.(land).(time).(fw).(crit).(regime) = ta_si_n.(land).(time) .* filt.(land).(time).(fw).(crit).(regime);

                        nanfilt.(regime) = nan(size(ta_si_t.(land).(time).(fw).(crit).(regime)));
                        nanfilt.(regime)(~isnan(ta_si_t.(land).(time).(fw).(crit).(regime))) = 1;

                        % take cosine-weighted meridional average
                        for d = {'all', 'nh', 'sh', 'tp'}; domain = d{1};
                            if strcmp(regime, 'rce')
                                if strcmp(domain, 'all')
                                    cosw = repmat(cosd(lat)', [size(nanfilt.(regime), 1) 1]);
                                    denm = ta_si_t.(land).(time).(fw).(crit).(regime);
                                    nume = nanfilt.(regime);
                                elseif strcmp(domain, 'nh')
                                    cosw = repmat(cosd(lat(lat>0))', [size(nanfilt.(regime), 1) 1]);
                                    denm = ta_si_t.(land).(time).(fw).(crit).(regime)(:,lat>0,:);
                                    nume = nanfilt.(regime)(:,lat>0,:);
                                elseif strcmp(domain, 'sh')
                                    cosw = repmat(cosd(lat(lat<0))', [size(nanfilt.(regime), 1) 1]);
                                    denm = ta_si_t.(land).(time).(fw).(crit).(regime)(:,lat<0,:);
                                    nume = nanfilt.(regime)(:,lat<0,:);
                                elseif strcmp(domain, 'tp')
                                    cosw = repmat(cosd(lat(abs(lat)<10))', [size(nanfilt.(regime), 1) 1]);
                                    denm = ta_si_t.(land).(time).(fw).(crit).(regime)(:,abs(lat)<10,:);
                                    nume = nanfilt.(regime)(:,abs(lat)<10,:);
                                end
                            elseif strcmp(regime, 'rae')
                                if strcmp(domain, 'all')
                                    cosw = repmat(cosd(lat)', [size(nanfilt.(regime), 1) 1]);
                                    denm = ta_si_t.(land).(time).(fw).(crit).(regime);
                                    nume = nanfilt.(regime);
                                elseif strcmp(domain, 'nh')
                                    cosw = repmat(cosd(lat(lat>0))', [size(nanfilt.(regime), 1) 1]);
                                    denm = ta_si_t.(land).(time).(fw).(crit).(regime)(:,lat>0,:);
                                    nume = nanfilt.(regime)(:,lat>0,:);
                                elseif strcmp(domain, 'sh')
                                    cosw = repmat(cosd(lat(lat<0))', [size(nanfilt.(regime), 1) 1]);
                                    denm = ta_si_t.(land).(time).(fw).(crit).(regime)(:,lat<0,:);
                                    nume = nanfilt.(regime)(:,lat<0,:);
                                end
                            elseif strcmp(regime, 'rcae')
                                if strcmp(domain, 'all')
                                    cosw = repmat(cosd(lat)', [size(nanfilt.(regime), 1) 1]);
                                    denm = ta_si_t.(land).(time).(fw).(crit).(regime);
                                    nume = nanfilt.(regime);
                                elseif strcmp(domain, 'nh')
                                    cosw = repmat(cosd(lat(lat>0))', [size(nanfilt.(regime), 1) 1]);
                                    denm = ta_si_t.(land).(time).(fw).(crit).(regime)(:,lat>0,:);
                                    nume = nanfilt.(regime)(:,lat>0,:);
                                elseif strcmp(domain, 'sh')
                                    cosw = repmat(cosd(lat(lat<0))', [size(nanfilt.(regime), 1) 1]);
                                    denm = ta_si_t.(land).(time).(fw).(crit).(regime)(:,lat<0,:);
                                    nume = nanfilt.(regime)(:,lat<0,:);
                                end
                            end

                            ta_si.(regime).(domain).(fw).(crit).(land).(time) = nansum(cosw.*denm, 2) ./ nansum(cosw.*nume, 2); % area-weighted meridional average
                            ta_si.(regime).(domain).(fw).(crit).(land).(time) = squeeze(nanmean(ta_si.(regime).(domain).(fw).(crit).(land).(time), 1)); % zonal average
                        end
                    end
                end
            end
        end
    end

    % save filtered data
    printname = [foldername 'ta_si'];
    if ~exist(foldername, 'dir')
        mkdir(foldername)
    end
    save(printname, 'ta_si');
end % process temperature in RCE/RAE regimes
function proc_ma_si(type, par)
% calculate moist adiabats
    if strcmp(type, 'era5') | strcmp(type, 'erai')
        foldername = sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/eps_%g_ga_%g/', type, par.lat_interp, par.ep, par.ga);
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s', type);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s', type);
    elseif strcmp(type, 'merra2')
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
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s', type, par.echam.clim);
    end
    load(sprintf('%s/grid.mat', prefix)); % read grid data
    load(sprintf('%s/srfc.mat', prefix)); % read surface variable data
    % load(sprintf('%s/%s/masks.mat', prefix_proc, par.lat_interp)); % load land and ocean masks
    load(sprintf('%s/%s/eps_%g_ga_%g/rcae_alt_t.mat', prefix_proc, par.lat_interp, par.ep, par.ga)); % load rcae data

    if strcmp(par.lat_interp, 'std')
        lat = par.lat_std;
    else
        lat = grid.dim3.lat;
    end

    for fn = fieldnames(srfc)'
        srfc.(fn{1}) = permute(srfc.(fn{1}), [2 1 3]); % bring lat to front
        srfc.(fn{1}) = interp1(grid.dim2.lat, srfc.(fn{1}), lat); % interpolate to standard grid
        srfc.(fn{1}) = permute(srfc.(fn{1}), [2 1 3]); % reorder to original dims
        % for l = {'lo', 'l', 'o'}; land = l{1};
        for l = {'lo'}; land = l{1};
            if strcmp(land, 'lo'); srfc_n.(fn{1}).(land) = srfc.(fn{1});
            % elseif strcmp(land, 'l'); srfc_n.(fn{1}).(land) = srfc.(fn{1}).*mask.ocean; % filter out ocean
            % elseif strcmp(land, 'o'); srfc_n.(fn{1}).(land) = srfc.(fn{1}).*mask.land; % filter out land
            end
        end
    end

    if any(strcmp(type, {'era5', 'erai'})); f_vec = par.era.fw;
    elseif any(strcmp(type, 'merra2')); f_vec = par.merra2.fw;
    elseif strcmp(type, 'gcm'); f_vec = par.gcm.fw;
    elseif strcmp(type, 'echam'); f_vec = par.echam.fw; end
    for f = f_vec; fw = f{1};
        for c = fieldnames(rcae_alt_t.lo.ann.(fw))'; crit = c{1};
            % for l = {'lo', 'l', 'o'}; land = l{1};
            for l = {'lo'}; land = l{1};
                % for t = {'ann', 'djf', 'jja', 'mam', 'son'}; time = t{1};
                for t = {'ann'}; time = t{1};
                    for re = {'rce', 'rae', 'rcae'}; regime = re{1};
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

                            filt.(land).(time).(fw).(crit).(regime) = nan(size(rcae_alt_t.(land).(time).(fw).(crit)));
                            if strcmp(regime, 'rce'); filt.(land).(time).(fw).(crit).(regime)(rcae_alt_t.(land).(time).(fw).(crit)==1)=1; % set RCE=1, elsewhere nan
                            elseif strcmp(regime, 'rae'); filt.(land).(time).(fw).(crit).(regime)(rcae_alt_t.(land).(time).(fw).(crit)==-1)=1; % set RAE=1, elsewhere nan
                            elseif strcmp(regime, 'rcae'); filt.(land).(time).(fw).(crit).(regime)(rcae_alt_t.(land).(time).(fw).(crit)==0)=1; % set RCAE=1, elsewhere nan
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
                                        cosw = repmat(cosd(lat(lat>0))', [size(nanfilt.(regime), 1) 1]);
                                        denm = srfc_tf.(land).(time).(fw).(crit).(regime).(vname)(:, lat>0);
                                        nume = nanfilt.(regime)(:, lat>0);
                                    elseif strcmp(domain, 'sh')
                                        cosw = repmat(cosd(lat(lat<0))', [size(nanfilt.(regime), 1) 1]);
                                        denm = srfc_tf.(land).(time).(fw).(crit).(regime).(vname)(:, lat<0);
                                        nume = nanfilt.(regime)(:, lat<0);
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
                                elseif strcmp(regime, 'rcae')
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

                                ma_si.(regime).(domain).(fw).(crit).(land).(time).(vname) = nansum(cosw.*denm, 2) ./ nansum(cosw.*nume, 2); % weighted meridional average
                                ma_si.(regime).(domain).(fw).(crit).(land).(time).(vname) = squeeze(nanmean(ma_si.(regime).(domain).(fw).(crit).(land).(time).(vname), 1)); % zonal average

                            end % end domain loop
                        end % end srfc variables loop

                    end % end RCE/RAE regime loop

                    if strcmp(type, 'era5') | strcmp(type, 'erai')
                        ma_si.rce.all.(fw).(crit).(land).(time).ta = calc_ma_dew_si(ma_si.rce.all.(fw).(crit).(land).(time), grid.dim3.plev, par, type, grid); % compute moist adiabat with dew_si point temperature
                        ma_si.rce.tp.(fw).(crit).(land).(time).ta = calc_ma_dew_si(ma_si.rce.tp.(fw).(crit).(land).(time), grid.dim3.plev, par, type, grid); % compute moist adiabat with dew_si point temperature
                        ma_si.rce.nh.(fw).(crit).(land).(time).ta = calc_ma_dew_si(ma_si.rce.nh.(fw).(crit).(land).(time), grid.dim3.plev, par, type, grid); % compute moist adiabat with dew_si point temperature
                        ma_si.rce.sh.(fw).(crit).(land).(time).ta = calc_ma_dew_si(ma_si.rce.sh.(fw).(crit).(land).(time), grid.dim3.plev, par, type, grid); % compute moist adiabat with dew_si point temperature
                        ma_si.rcae.all.(fw).(crit).(land).(time).ta = calc_ma_dew_si(ma_si.rcae.all.(fw).(crit).(land).(time), grid.dim3.plev, par, type, grid); % compute moist adiabat with dew_si point temperature
                        ma_si.rcae.nh.(fw).(crit).(land).(time).ta = calc_ma_dew_si(ma_si.rcae.nh.(fw).(crit).(land).(time), grid.dim3.plev, par, type, grid); % compute moist adiabat with dew_si point temperature
                        ma_si.rcae.sh.(fw).(crit).(land).(time).ta = calc_ma_dew_si(ma_si.rcae.sh.(fw).(crit).(land).(time), grid.dim3.plev, par, type, grid); % compute moist adiabat with dew_si point temperature
                    elseif strcmp(type, 'merra2')
                        ma_si.rce.all.(fw).(crit).(land).(time).ta = calc_ma_dew_si(ma_si.rce.all.(fw).(crit).(land).(time), grid.dim3.plev, par, type, grid); % compute moist adiabat with dew_si point temperature
                        ma_si.rce.tp.(fw).(crit).(land).(time).ta = calc_ma_dew_si(ma_si.rce.tp.(fw).(crit).(land).(time), grid.dim3.plev, par, type, grid); % compute moist adiabat with dew_si point temperature
                        ma_si.rce.nh.(fw).(crit).(land).(time).ta = calc_ma_dew_si(ma_si.rce.nh.(fw).(crit).(land).(time), grid.dim3.plev, par, type, grid); % compute moist adiabat with dew_si point temperature
                        ma_si.rce.sh.(fw).(crit).(land).(time).ta = calc_ma_dew_si(ma_si.rce.sh.(fw).(crit).(land).(time), grid.dim3.plev, par, type, grid); % compute moist adiabat with dew_si point temperature
                        ma_si.rcae.all.(fw).(crit).(land).(time).ta = calc_ma_dew_si(ma_si.rcae.all.(fw).(crit).(land).(time), grid.dim3.plev, par, type, grid); % compute moist adiabat with dew_si point temperature
                        ma_si.rcae.nh.(fw).(crit).(land).(time).ta = calc_ma_dew_si(ma_si.rcae.nh.(fw).(crit).(land).(time), grid.dim3.plev, par, type, grid); % compute moist adiabat with dew_si point temperature
                        ma_si.rcae.sh.(fw).(crit).(land).(time).ta = calc_ma_dew_si(ma_si.rcae.sh.(fw).(crit).(land).(time), grid.dim3.plev, par, type, grid); % compute moist adiabat with dew_si point temperature
                    elseif strcmp(type, 'gcm')
                        ma_si.rce.all.(fw).(crit).(land).(time).ta = calc_ma_hurs_si(ma_si.rce.all.(fw).(crit).(land).(time), grid.dim3.plev, par, grid); % compute moist adiabat with RH
                        ma_si.rce.tp.(fw).(crit).(land).(time).ta = calc_ma_hurs_si(ma_si.rce.tp.(fw).(crit).(land).(time), grid.dim3.plev, par, grid); % compute moist adiabat with RH
                        ma_si.rce.nh.(fw).(crit).(land).(time).ta = calc_ma_hurs_si(ma_si.rce.nh.(fw).(crit).(land).(time), grid.dim3.plev, par, grid); % compute moist adiabat with RH
                        ma_si.rce.sh.(fw).(crit).(land).(time).ta = calc_ma_hurs_si(ma_si.rce.sh.(fw).(crit).(land).(time), grid.dim3.plev, par, grid); % compute moist adiabat with RH
                        ma_si.rcae.all.(fw).(crit).(land).(time).ta = calc_ma_hurs_si(ma_si.rcae.all.(fw).(crit).(land).(time), grid.dim3.plev, par, grid); % compute moist adiabat with RH
                        ma_si.rcae.nh.(fw).(crit).(land).(time).ta = calc_ma_hurs_si(ma_si.rcae.nh.(fw).(crit).(land).(time), grid.dim3.plev, par, grid); % compute moist adiabat with RH
                        ma_si.rcae.sh.(fw).(crit).(land).(time).ta = calc_ma_hurs_si(ma_si.rcae.sh.(fw).(crit).(land).(time), grid.dim3.plev, par, grid); % compute moist adiabat with RH
                    end

                end % end time average loop
            end % end land option loop
        end % end RCAE definition loop
    end % end MSE/DSE framework loop

    % save data into mat file
    printname = [foldername 'ma_si.mat'];
    if ~exist(foldername, 'dir')
        mkdir(foldername)
    end
    save(printname, 'ma_si');

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

    if strcmp(par.lat_interp, 'std')
        lat = par.lat_std;
    else
        lat = grid.dim3.lat;
    end

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

    if strcmp(par.lat_interp, 'std')
        lat = par.lat_std;
    else
        lat = grid.dim3.lat;
    end

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

        % identify locations of rce following jakob et al. (2019) (abs(tediv) < 30 w/m^2)
        rcae.(fw).jak30 = zeros(size(flux.r1.(fw)));
        rcae.(fw).jak30(abs(flux.res.(fw)) < 30) = 1;
        % identify locations of rae using threshold gamma (ga)
        if strcmp(fw, 'dse'); rcae.(fw).jak30(flux.r1.mse > par.ga) = -1; % always use mse for rae
        else rcae.(fw).jak30(flux.r1.(fw) > par.ga) = -1; end;

        % identify locations of rce following jakob et al. (2019) (abs(tediv) < 10 w/m^2)
        rcae.(fw).jak10 = zeros(size(flux.r1.(fw)));
        rcae.(fw).jak10(abs(flux.res.(fw)) < 10) = 1;
        % identify locations of rae using threshold gamma (ga)
        if strcmp(fw, 'dse'); rcae.(fw).jak10(flux.r1.mse > par.ga) = -1; % always use mse for rae
        else rcae.(fw).jak10(flux.r1.(fw) > par.ga) = -1; end;

        % % add additional criteria for RCE that P-E>0
        % rcae.(fw).pe = zeros(size(flux.r1.(fw)));
        % if strcmp(type, 'era5') | strcmp(type, 'erai')
        %     rcae.(fw).pe(abs(flux.r1.(fw)) < par.ep & (flux.cp+flux.lsp+flux.e>0)) = 1; % note that evap is defined negative into atmosphere
        % elseif strcmp(type, 'gcm')
        %     rcae.(fw).pe(abs(flux.r1.(fw)) < par.ep & (flux.pr-flux.evspsbl>0)) = 1;
        % elseif strcmp(type, 'echam')
        %     rcae.(fw).pe(abs(flux.r1.(fw)) < par.ep & (flux.aprc+flux.aprl-flux.evap>0)) = 1;
        % end
        % if strcmp(fw, 'dse'); rcae.(fw).pe(flux.r1.mse > par.ga) = -1;
        % else rcae.(fw).pe(flux.r1.(fw) > par.ga) = -1; end;

        % % add additional criteria for RCE that (P_ls - E)<<1
        % rcae.(fw).cp = zeros(size(flux.r1.(fw)));
        % if strcmp(type, 'era5') | strcmp(type, 'erai')
        %     rcae.(fw).cp(abs(flux.r1.(fw)) < par.ep & abs(flux.lsp./flux.cp<par.ep_cp)) = 1; % note that evap is defined negative into atmosphere
        % elseif strcmp(type, 'gcm')
        %     rcae.(fw).cp(abs(flux.r1.(fw)) < par.ep & abs((flux.pr-flux.prc)./flux.prc<par.ep_cp)) = 1; % note that evap is defined negative into atmosphere
        % elseif strcmp(type, 'echam')
        %     rcae.(fw).cp(abs(flux.r1.(fw)) < par.ep & abs(flux.aprl./flux.aprc<par.ep_cp)) = 1; % note that evap is defined negative into atmosphere
        % end
        % if strcmp(fw, 'dse') rcae.(fw).cp(flux.r1.mse > par.ga) = -1;
        % else rcae.(fw).cp(flux.r1.(fw) > par.ga) = -1; end;

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

        % if size(flux.r1.(fw), 1)==length(par.lat_std) & size(flux.r1.(fw), 2)==12 % this criteria only works for zonally and time averaged data because vh is only a function of lat
        %     % add additional criteria for RCE that horizontal transport is weak
        %     rcae.(fw).vh2 = zeros(size(flux.r1.(fw)));
        %     % if strcmp(type, 'era5') | strcmp(type, 'erai')
        %         rcae.(fw).vh2(abs(flux.r1.(fw)) < par.ep & abs(vh.(fw)) < nanmax(nanmax(abs(vh.(fw))))/2 ) = 1;
        %     % elseif strcmp(type, 'gcm')
        %     %     rcae.(fw).vh2(abs(flux.r1.(fw)) < par.ep & abs(vh.(fw)) < nanmax(nanmax(abs(vh.(fw))))/2 ) = 1;
        %     % end
        %     if strcmp(fw, 'dse'); rcae.(fw).vh2(flux.r1.mse > par.ga) = -1;
        %     else; rcae.(fw).vh2(flux.r1.(fw) > par.ga) = -1; end;

        %     rcae.(fw).vh3 = zeros(size(flux.r1.(fw)));
        %     % if strcmp(type, 'era5') | strcmp(type, 'erai')
        %         rcae.(fw).vh3(abs(flux.r1.(fw)) < par.ep & abs(vh.(fw)) < nanmax(nanmax(abs(vh.(fw))))/3 ) = 1;
        %     % elseif strcmp(type, 'gcm')
        %     %     rcae.(fw).vh3(abs(flux.r1.(fw)) < par.ep & abs(vh.(fw)) < nanmax(nanmax(abs(vh.(fw))))/3 ) = 1;
        %     % end
        %     if strcmp(fw, 'dse') rcae.(fw).vh3(flux.r1.mse > par.ga) = -1;
        %     else rcae.(fw).vh3(flux.r1.(fw) > par.ga) = -1; end;

        %     rcae.(fw).vh4 = zeros(size(flux.r1.(fw)));
        %     % if strcmp(type, 'era5') | strcmp(type, 'erai')
        %         rcae.(fw).vh4(abs(flux.r1.(fw)) < par.ep & abs(vh.(fw)) < nanmax(nanmax(abs(vh.(fw))))/4 ) = 1;
        %     % elseif strcmp(type, 'gcm')
        %     %     rcae.(fw).vh4(abs(flux.r1.(fw)) < par.ep & abs(vh.(fw)) < nanmax(nanmax(abs(vh.(fw))))/4 ) = 1;
        %     % end
        %     if strcmp(fw, 'dse'); rcae.(fw).vh4(flux.r1.mse > par.ga) = -1;
        %     else; rcae.(fw).vh4(flux.r1.(fw) > par.ga) = -1; end;
        % end

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

        % identify locations of rce following jakob et al. (2019) (abs(tediv) < 30 w/m^2)
        rcae.(fw).jak30 = zeros(size(r1));
        rcae.(fw).jak30(abs(flux.res.(fw)) < 30) = 1;
        % identify locations of rae using threshold gamma (ga)
        if strcmp(fw, 'dse'); rcae.(fw).jak30(flux.r1.mse > par.ga) = -1; % always use mse for rae
        else rcae.(fw).jak30(r1 > par.ga) = -1; end;

        % identify locations of rce following jakob et al. (2019) (abs(tediv) < 10 w/m^2)
        rcae.(fw).jak10 = zeros(size(r1));
        rcae.(fw).jak10(abs(flux.res.(fw)) < 10) = 1;
        % identify locations of rae using threshold gamma (ga)
        if strcmp(fw, 'dse'); rcae.(fw).jak10(flux.r1.mse > par.ga) = -1; % always use mse for rae
        else rcae.(fw).jak10(r1 > par.ga) = -1; end;

        % % add additional criteria for RCE that P-E>0
        % rcae.(fw).pe = zeros(size(r1));
        % if strcmp(type, 'era5') | strcmp(type, 'erai')
        %     rcae.(fw).pe(abs(r1) < par.ep & (flux.cp+flux.lsp+flux.e>0)) = 1; % note that evap is defined negative into atmosphere
        % elseif strcmp(type, 'gcm')
        %     rcae.(fw).pe(abs(r1) < par.ep & (flux.pr-flux.evspsbl>0)) = 1;
        % elseif strcmp(type, 'echam')
        %     rcae.(fw).pe(abs(flux.r1.(fw)) < par.ep & (flux.aprc+flux.aprl-flux.evap>0)) = 1;
        % end
        % if strcmp(fw, 'dse'); rcae.(fw).pe(flux.r1.mse > par.ga) = -1;
        % else rcae.(fw).pe(r1 > par.ga) = -1; end;

        % % add additional criteria for RCE that (P_ls - E)<<1
        % rcae.(fw).cp = zeros(size(r1));
        % if strcmp(type, 'era5') | strcmp(type, 'erai')
        %     rcae.(fw).cp(abs(r1) < par.ep & abs(flux.lsp./flux.cp<par.ep_cp)) = 1; % note that evap is defined negative into atmosphere
        % elseif strcmp(type, 'gcm')
        %     rcae.(fw).cp(abs(r1) < par.ep & abs((flux.pr-flux.prc)./flux.prc<par.ep_cp)) = 1; % note that evap is defined negative into atmosphere
        % elseif strcmp(type, 'echam')
        %     rcae.(fw).cp(abs(flux.r1.(fw)) < par.ep & abs(flux.aprl./flux.aprc<par.ep_cp)) = 1; % note that evap is defined negative into atmosphere
        % end
        % if strcmp(fw, 'dse') rcae.(fw).cp(flux.r1.mse > par.ga) = -1;
        % else rcae.(fw).cp(r1 > par.ga) = -1; end;

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

        % if size(r1, 1)==length(par.lat_std) & size(r1, 2)==12 % this criteria only works for zonally and time averaged data because vh is only a function of lat
        %     % add additional criteria for RCE that horizontal transport is weak
        %     rcae.(fw).vh2 = zeros(size(r1));
        %     % if strcmp(type, 'era5') | strcmp(type, 'erai')
        %         rcae.(fw).vh2(abs(r1) < par.ep & abs(vh.(fw)) < nanmax(nanmax(abs(vh.(fw))))/2 ) = 1;
        %     % elseif strcmp(type, 'gcm')
        %     %     rcae.(fw).vh2(abs(r1) < par.ep & abs(vh.(fw)) < nanmax(nanmax(abs(vh.(fw))))/2 ) = 1;
        %     % end
        %     if strcmp(fw, 'dse'); rcae.(fw).vh2(flux.r1.mse > par.ga) = -1;
        %     else; rcae.(fw).vh2(r1 > par.ga) = -1; end;

        %     rcae.(fw).vh3 = zeros(size(r1));
        %     % if strcmp(type, 'era5') | strcmp(type, 'erai')
        %         rcae.(fw).vh3(abs(r1) < par.ep & abs(vh.(fw)) < nanmax(nanmax(abs(vh.(fw))))/3 ) = 1;
        %     % elseif strcmp(type, 'gcm')
        %     %     rcae.(fw).vh3(abs(r1) < par.ep & abs(vh.(fw)) < nanmax(nanmax(abs(vh.(fw))))/3 ) = 1;
        %     % end
        %     if strcmp(fw, 'dse') rcae.(fw).vh3(flux.r1.mse > par.ga) = -1;
        %     else rcae.(fw).vh3(r1 > par.ga) = -1; end;

        %     rcae.(fw).vh4 = zeros(size(r1));
        %     % if strcmp(type, 'era5') | strcmp(type, 'erai')
        %         rcae.(fw).vh4(abs(r1) < par.ep & abs(vh.(fw)) < nanmax(nanmax(abs(vh.(fw))))/4 ) = 1;
        %     % elseif strcmp(type, 'gcm')
        %     %     rcae.(fw).vh4(abs(r1) < par.ep & abs(vh.(fw)) < nanmax(nanmax(abs(vh.(fw))))/4 ) = 1;
        %     % end
        %     if strcmp(fw, 'dse'); rcae.(fw).vh4(flux.r1.mse > par.ga) = -1;
        %     else; rcae.(fw).vh4(r1 > par.ga) = -1; end;
        % end

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

        % identify locations of rce following jakob et al. (2019) (abs(tediv) < 30 w/m^2)
        rcae.(fw).jak30 = zeros(size(flux.r1.(fw)));
        rcae.(fw).jak30(flux.res.(fw) < 30) = 1;
        % identify locations of rae using threshold gamma (ga)
        if strcmp(fw, 'dse'); rcae.(fw).jak30(flux.r1.mse > par.ga) = -1; % always use mse for rae
        else rcae.(fw).jak30(flux.r1.(fw) > par.ga) = -1; end;

        % identify locations of rce following jakob et al. (2019) (abs(tediv) < 10 w/m^2)
        rcae.(fw).jak10 = zeros(size(flux.r1.(fw)));
        rcae.(fw).jak10(flux.res.(fw) < 10) = 1;
        % identify locations of rae using threshold gamma (ga)
        if strcmp(fw, 'dse'); rcae.(fw).jak10(flux.r1.mse > par.ga) = -1; % always use mse for rae
        else rcae.(fw).jak10(flux.r1.(fw) > par.ga) = -1; end;

        % % add additional criteria for RCE that P-E>0
        % rcae.(fw).pe = zeros(size(flux.r1.(fw)));
        % if strcmp(type, 'era5') | strcmp(type, 'erai')
        %     rcae.(fw).pe(flux.r1.(fw) < par.ep & (flux.cp+flux.lsp+flux.e>0)) = 1; % note that evap is defined negative into atmosphere
        % elseif strcmp(type, 'gcm')
        %     rcae.(fw).pe(flux.r1.(fw) < par.ep & (flux.pr-flux.evspsbl>0)) = 1;
        % elseif strcmp(type, 'echam')
        %     rcae.(fw).pe(flux.r1.(fw) < par.ep & (flux.aprc+flux.aprl-flux.evap>0)) = 1;
        % end
        % if strcmp(fw, 'dse'); rcae.(fw).pe(flux.r1.mse > par.ga) = -1;
        % else rcae.(fw).pe(flux.r1.(fw) > par.ga) = -1; end;

        % % add additional criteria for RCE that (P_ls - E)<<1
        % rcae.(fw).cp = zeros(size(flux.r1.(fw)));
        % if strcmp(type, 'era5') | strcmp(type, 'erai')
        %     rcae.(fw).cp(flux.r1.(fw) < par.ep & abs(flux.lsp./flux.cp<par.ep_cp)) = 1; % note that evap is defined negative into atmosphere
        % elseif strcmp(type, 'gcm')
        %     rcae.(fw).cp(flux.r1.(fw) < par.ep & abs((flux.pr-flux.prc)./flux.prc<par.ep_cp)) = 1; % note that evap is defined negative into atmosphere
        % elseif strcmp(type, 'echam')
        %     rcae.(fw).cp(flux.r1.(fw) < par.ep & abs(flux.aprl./flux.aprc<par.ep_cp)) = 1; % note that evap is defined negative into atmosphere
        % end
        % if strcmp(fw, 'dse') rcae.(fw).cp(flux.r1.mse > par.ga) = -1;
        % else rcae.(fw).cp(flux.r1.(fw) > par.ga) = -1; end;

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

        % identify locations of rce following jakob et al. (2019) (abs(tediv) < 30 w/m^2)
        rcae.(fw).jak30 = zeros(size(r1));
        rcae.(fw).jak30(flux.res.(fw) < 30) = 1;
        % identify locations of rae using threshold gamma (ga)
        if strcmp(fw, 'dse'); rcae.(fw).jak30(flux.r1.mse > par.ga) = -1; % always use mse for rae
        else rcae.(fw).jak30(r1 > par.ga) = -1; end;

        % identify locations of rce following jakob et al. (2019) (abs(tediv) < 10 w/m^2)
        rcae.(fw).jak10 = zeros(size(r1));
        rcae.(fw).jak10(flux.res.(fw) < 10) = 1;
        % identify locations of rae using threshold gamma (ga)
        if strcmp(fw, 'dse'); rcae.(fw).jak10(flux.r1.mse > par.ga) = -1; % always use mse for rae
        else rcae.(fw).jak10(r1 > par.ga) = -1; end;

        % % add additional criteria for RCE that P-E>0
        % rcae.(fw).pe = zeros(size(r1));
        % if strcmp(type, 'era5') | strcmp(type, 'erai')
        %     rcae.(fw).pe(r1 < par.ep & (flux.cp+flux.lsp+flux.e>0)) = 1; % note that evap is defined negative into atmosphere
        % elseif strcmp(type, 'gcm')
        %     rcae.(fw).pe(r1 < par.ep & (flux.pr-flux.evspsbl>0)) = 1;
        % elseif strcmp(type, 'echam')
        %     rcae.(fw).pe(r1 < par.ep & (flux.aprc + flux.aprl -flux.evap>0)) = 1;
        % end
        % if strcmp(fw, 'dse'); rcae.(fw).pe(flux.r1.mse > par.ga) = -1;
        % else rcae.(fw).pe(r1 > par.ga) = -1; end;

        % % add additional criteria for RCE that (P_ls - E)<<1
        % rcae.(fw).cp = zeros(size(r1));
        % if strcmp(type, 'era5') | strcmp(type, 'erai')
        %     rcae.(fw).cp(r1 < par.ep & abs(flux.lsp./flux.cp<par.ep_cp)) = 1; % note that evap is defined negative into atmosphere
        % elseif strcmp(type, 'gcm')
        %     rcae.(fw).cp(r1 < par.ep & abs((flux.pr-flux.prc)./flux.prc<par.ep_cp)) = 1; % note that evap is defined negative into atmosphere
        % elseif strcmp(type, 'echam')
        %     rcae.(fw).cp(r1 < par.ep & abs(flux.aprl./flux.aprc<par.ep_cp)) = 1; % note that evap is defined negative into atmosphere
        % end
        % if strcmp(fw, 'dse') rcae.(fw).cp(flux.r1.mse > par.ga) = -1;
        % else rcae.(fw).cp(r1 > par.ga) = -1; end;

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

function ma_out = calc_ma_dew(ma_in, plev, par, type)
% compute moist adiabat based on dew point temperature
    if any(strcmp(type, {'erai', 'era5'}))
        ps = ma_in.sp; tas = ma_in.t2m; tdas = ma_in.d2m;
    elseif strcmp(type, 'echam')
        ps = ma_in.aps; tas = ma_in.temp2; tdas = ma_in.dew2;
    end

    if isnan(ps) | isnan(tas) | isnan(tdas)
        ma_out = nan(size(plev)); % if any of the initial values are nan, ouput nan profile
    else
    if isfield(ma_in, 'pinit')
        pa_cd(1) = ma_in.pinit;
    else
        pa_cd(1) = ps;
    end
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

        ta_s = interp1(pa_s, ta_s, plev, 'spline', nan);
        qsat_s = interp1(pa_s, qsat_s, plev, 'spline', nan);
        dtadpa_s = interp1(pa_s, dtadpa_s, plev, 'spline', nan);

        % OUTPUT
        ma_out = ta_s;
    end

end
function ma_out = calc_ma_dew_si(ma_in, plev, par, type, grid)
    if any(strcmp(type, {'erai', 'era5'}))
        ps = ma_in.sp; tas = ma_in.t2m; tdas = ma_in.d2m;
    elseif strcmp(type, 'echam')
        ps = ma_in.aps; tas = ma_in.temp2; tdas = ma_in.dew2;
    end

    if isnan(ps) | isnan(tas) | isnan(tdas)
        ma_out = nan(size(grid.dim3.si)); % if any of the initial values are nan, ouput nan profile
    else
        % integrate temperature profiles over height
        z_cd(1) = ma_in.zs;

    if isfield(ma_in, 'pinit')
        pa_cd(1) = ma_in.pinit;
    else
        pa_cd(1) = ps;
    end
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
        ta_s(nanfilter) = [];

        ma_si = interp1(real(pa_s)/ps, real(ta_s), grid.dim3.si, 'spline', nan);

        % output temperature
        ma_out = ma_si;
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

    if isfield(ma_in, 'pinit')
        pa_cd(1) = ma_in.pinit;
    else
        pa_cd(1) = ps;
    end
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

        pa_s = interp1(z_s, pa_s, z_int, 'spline', nan);
        ta_s = interp1(z_s, ta_s, z_int, 'spline', nan);
        qsat_s = interp1(z_s, qsat_s, z_int, 'spline', nan);
        dtadz_s = interp1(z_s, dtadz_s, z_int, 'spline', nan);

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

    if isfield(ma_in, 'pinit')
        pa_cd(1) = ma_in.pinit;
    else
        pa_cd(1) = ps;
    end
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

        dtadz_s = -1e3*interp1(real(pa_s), real(dtadz_s), plev, 'spline', nan); % output in K/km at standard pressure grid

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

    if isfield(ma_in, 'pinit')
        pa_cd(1) = ma_in.pinit;
    else
        pa_cd(1) = ps;
    end
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

        dtadz_s = interp1(real(pa_s)/ps, real(dtadz_s), grid.dim3.si, 'spline', nan);

        % output temperature
        ma_out = dtadz_s;
    end

end
function ma_out = calc_ma_hurs(ma_in, plev, par)
    if isnan(ma_in.ps) | isnan(ma_in.tas) | isnan(ma_in.hurs)
        ma_out = nan(size(plev)); % if any of the initial values are nan, ouput nan profile
    else
        if isfield(ma_in, 'pinit')
            pa_cd(1) = ma_in.pinit;
        else
            pa_cd(1) = ma_in.ps;
        end
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

        ta_s = interp1(pa_s, ta_s, plev, 'spline', nan);
        qsat_s = interp1(pa_s, qsat_s, plev, 'spline', nan);
        dtadpa_s = interp1(pa_s, dtadpa_s, plev, 'spline', nan);

        % OUTPUT
        ma_out = ta_s;
    end

end % compute moist adiabat based on RH
function ma_out = calc_ma_hurs_si(ma_in, plev, par, grid)
    if isnan(ma_in.ps) | isnan(ma_in.tas) | isnan(ma_in.hurs) | isnan(ma_in.zs)
        ma_out = nan(size(grid.dim3.si)); % if any of the initial values are nan, ouput nan profile
    else
        % integrate temperature profiles over height
        z_cd(1) = ma_in.zs;
        if isfield(ma_in, 'pinit')
            pa_cd(1) = ma_in.pinit;
        else
            pa_cd(1) = ma_in.ps;
        end
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
        % ta_s(:, 1) = real([ta_cd; ta_pa_cs(2:end,1)]);
        ta_s(:, 1) = [ta_cd; ta_pa_cs(2:end,1)];
        pa_s(:, 1) = [pa_cd(:); ta_pa_cs(2:end,2)];
        qsat_s(:, 1) = [qsat_cd(:, 1); qsat_cs(2:end)];
        dtadz_s(:, 1) = [deriv_dry; deriv_cs(2:end)];

        nanfilter = isnan(real(pa_s));
        pa_s(nanfilter) = [];
        ta_s(nanfilter) = [];

        ma_si = interp1(real(pa_s)/ma_in.ps, real(ta_s), grid.dim3.si, 'spline', nan);

        % output temperature
        ma_out = ma_si;
    end

end
function ma_out = calc_maz_hurs(ma_in, z_int, par)
    if isnan(ma_in.ps) | isnan(ma_in.tas) | isnan(ma_in.hurs) | isnan(ma_in.zs)
        ma_out = nan(size(z_int)); % if any of the initial values are nan, ouput nan profile
    else
        % integrate temperature profiles over height
        z_cd(1) = ma_in.zs;
        if isfield(ma_in, 'pinit')
            pa_cd(1) = ma_in.pinit;
        else
            pa_cd(1) = ma_in.ps;
        end
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

        pa_s = interp1(z_s, pa_s, z_int, 'spline', nan);
        ta_s = interp1(z_s, ta_s, z_int, 'spline', nan);
        qsat_s = interp1(z_s, qsat_s, z_int, 'spline', nan);
        dtadz_s = interp1(z_s, dtadz_s, z_int, 'spline', nan);

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
        if isfield(ma_in, 'pinit')
            pa_cd(1) = ma_in.pinit;
        else
            pa_cd(1) = ma_in.ps;
        end
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

        dtadz_s = -1e3*interp1(real(pa_s), real(dtadz_s), plev, 'spline', nan); % output in K/km at standard pressure grid

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
        if isfield(ma_in, 'pinit')
            pa_cd(1) = ma_in.pinit;
        else
            pa_cd(1) = ma_in.ps;
        end
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

        dtadz_s = -1e3*interp1(real(pa_s)/ma_in.ps, real(dtadz_s), grid.dim3.si, 'spline', nan); % output in K/km at standard pressure grid

        % output temperature
        ma_out = dtadz_s;
    end

end

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

    if strcmp(par.lat_interp, 'std')
        lat = par.lat_std;
    else
        lat = grid.dim3.lat;
    end

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

    tai_sm.l = tai_sm.lo.*surface_mask.*mask.ocean_vert; % filter tai with surface mask
    tai_sm.o = tai_sm.lo.*surface_mask.*mask.land_vert; % filter tai with surface mask
    tai_sm.lo = permute(tai_sm.lo, [1 2 4 3]); % bring plev to last dimension
    tai_sm.l = permute(tai_sm.l, [1 2 4 3]); % bring plev to last dimension
    tai_sm.o = permute(tai_sm.o, [1 2 4 3]); % bring plev to last dimension

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

    if strcmp(par.lat_interp, 'std')
        lat = par.lat_std;
    else
        lat = grid.dim3.lat;
    end

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
