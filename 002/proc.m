clc; close all; clear variables;

addpath(genpath('/project2/tas1/miyawaki/matlab'));

%% set parameters
par.lat_interp = 'std'; % which latitudinal grid to interpolate to: don (donohoe, coarse), era (native ERA-Interim, fine), or std (custom, very fine)
par.lat_std = transpose(-90:0.25:90); % define standard latitude grid for 'std' interpolation
par.ep_swp = 0.3; % threshold for RCE definition. RCE is defined as where abs(R1) < ep
par.ep_cp = 0.5; % threshold for RCE definition. RCE is defined as where abs(R1) < ep
par.ma_type = 'std'; % choose the type of moist adiabat: reversible, pseudo, or std
par.frz = 0; % consider latent heat of fusion in moist adiabat?
par.pa_span = [1000 100]*100; % pressure range for calculating moist adiabat
par.dpa = -10; % pressure increment for integrating dry adiabat section of moist adiabat (matters for how accurately the LCL is computed)
par.cpd = 1005.7; par.Rd = 287; par.Rv = 461; par.g = 9.81; par.L = 2.501e6; par.a = 6357e3; par.eps = 0.622; % common constants, all in SI units for brevity
gcm_info

%% call functions
type = 'era5'; % data type to run analysis on
% choose_proc(type, par)
for k=1:length(par.gcm_models); par.model = par.gcm_models{k};
    type = 'gcm';
    % choose_proc(type, par)
end

for i=1:length(par.ep_swp); par.ep = par.ep_swp(i); par.ga = 1-par.ep;
    type = 'era5';
    % choose_proc_ep(type, par)
    for k = 1:length(par.gcm_models); par.model = par.gcm_models{k};
        type = 'gcm';
        choose_proc_ep(type, par)
    end
end

%% define functions
function choose_proc(type, par)
    % proc_flux(type, par) % calculate energy fluxes in the vertically-integrated MSE budget using ERA-Interim data
    % proc_net_flux(type, par) % calculate net energy fluxes at TOA and surface
    save_mask(type, par) % save land and ocean masks once (faster than creating mask every time I need it)
end
function proc_flux(type, par)
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
        flux.(fn{1}) = permute(rad.(fn{1}), [2 1 3]);
        flux.(fn{1}) = interp1(grid.dim2.lat, flux.(fn{1}), par.lat_std, 'spline');
        flux.(fn{1}) = permute(flux.(fn{1}), [2 1 3]);
    end
    for fn = pe_vars % interpolate to std lat
        flux.(fn{1}) = permute(pe.(fn{1}), [2 1 3]);
        flux.(fn{1}) = interp1(grid.dim2.lat, flux.(fn{1}), par.lat_std, 'spline');
        flux.(fn{1}) = permute(flux.(fn{1}), [2 1 3]);
    end
    for fn = stf_vars % interpolate to std lat
        flux.(fn{1}) = permute(stf.(fn{1}), [2 1 3]);
        flux.(fn{1}) = interp1(grid.dim2.lat, flux.(fn{1}), par.lat_std, 'spline');
        flux.(fn{1}) = permute(flux.(fn{1}), [2 1 3]);
    end

    if strcmp(type, 'era5') | strcmp(type, 'erai')
        flux.ra = flux.tsr - flux.ssr + flux.ttr - flux.str; % compute net radiative cooling from radiative fluxes
        % compute surface turbulent fluxes directly from INTP data
        % multiply by negative to define flux from surface to atmosphere as positive
        flux.stf = -( flux.sshf + flux.slhf );
    elseif strcmp(type, 'gcm')
        flux.ra = flux.rsdt - flux.rsut + flux.rsus - flux.rsds + flux.rlus - flux.rlds - flux.rlut; % calculate atmospheric radiative cooling
        flux.stf = flux.hfls + flux.hfss;
    end

    flux.res = flux.ra + flux.stf; % infer MSE tendency and flux divergence as residuals
    % compute northward MSE transport using the residual data
    tediv_t = squeeze(nanmean(flux.res, 3)); % take time average
    tediv_tz = trapz(deg2rad(grid.dim2.lon), tediv_t, 1); % zonally integrate
    vh = cumtrapz(deg2rad(lat), par.a^2*cosd(lat).*tediv_tz'); % cumulatively integrate in latitude
    flux.r1 = (flux.res)./flux.ra; % calculate nondimensional number R1 disregarding MSE budget closure
    flux.r2 = flux.stf./flux.ra; % calculate nondimensional number R2 disregarding MSE budget closure

    if strcmp(type, 'era5') | strcmp(type, 'erai')
        % take zonal averages
        for fn = {'ra', 'stf', 'sshf', 'slhf', 'res', 'r1', 'r2', 'cp', 'lsp', 'e'}
            flux_z.(fn{1}) = squeeze(nanmean(flux.(fn{1}), 1));
        end
        foldername = sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/', type, par.lat_interp);
    elseif strcmp(type, 'gcm')
        % take zonal averages
        for fn = {'ra', 'stf', 'hfls', 'hfss', 'res', 'r1', 'r2', 'prc', 'pr', 'evspsbl'}
            flux_z.(fn{1}) = squeeze(nanmean(flux.(fn{1}), 1));
        end
        foldername = sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s/', type, par.model, par.lat_interp);
    end

    % save energy flux data into mat file
    printname = [foldername 'flux.mat'];
    if ~exist(foldername, 'dir')
        mkdir(foldername)
    end
    % save only what is necessary
    save(printname, 'flux', 'lat');

    % save energy flux data into mat file
    printname = [foldername 'flux_z.mat'];
    if ~exist(foldername, 'dir')
        mkdir(foldername)
    end
    save(printname, 'flux_z', 'vh', 'lat');

end
function proc_net_flux(type, par)
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
function save_mask(type, par)
    if strcmp(type, 'era5') | strcmp(type, 'erai')
        foldername = sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/', type, par.lat_interp);
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s', type);
    elseif strcmp(type, 'gcm')
        foldername = sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s/', type, par.model, par.lat_interp);
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s', type, par.model);
    end

    load(sprintf('%s/grid.mat', prefix)); % read grid data

    lat = par.lat_std;

    mask.land = remove_land(lat, grid.dim3.lon, 12);
    mask.ocean = remove_ocean(lat, grid.dim3.lon, 12);

    % save masks
    printname = [foldername 'masks'];
    if ~exist(foldername, 'dir')
        mkdir(foldername)
    end
    save(printname, 'mask');

end

function choose_proc_ep(type, par)
    % proc_rcae(type, par) % calculate RCE and RAE regimes
    % proc_temp(type, par) % calculate RCE and RAE temperature profiles
    % proc_ma(type, par) % calculate moist adiabats corresponding to RCE profiles
    proc_temp_mon_lat(type, par) % calculate mon x lat temperature profiles
    % proc_ma_mon_lat(type, par) % calculate mon x lat moist adiabats
end
function proc_rcae(type, par)
    for f = {'flux', 'flux_z'}; ftype = f{1};
        if strcmp(type, 'era5') | strcmp(type, 'erai')
            filename = sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s.mat', type, par.lat_interp, ftype); % read ERA5 zonally averaged flux
            foldername = sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/eps_%g/', type, par.lat_interp, par.ep);
            prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s', type);
            prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s', type, par.lat_interp);
        elseif strcmp(type, 'gcm')
            filename = sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s/%s.mat', type, par.model, par.lat_interp, ftype); % read gcm zonally averaged flux
            foldername = sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s/eps_%g/', type, par.model, par.lat_interp, par.ep);
            prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s', type, par.model);
            prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s', type, par.model, par.lat_interp);
        end
        if ~exist(filename); error(sprintf('Data does not exist. Please run proc_%s.m first.', ftype)); else
            load(filename);
        end

        load(sprintf('%s/masks.mat', prefix_proc)); % load land and ocean masks

        % identify locations of RCE and RAE
        if strcmp(ftype, 'flux')
            rcae = def_rcae(type, flux, par); % lon x lat x time structure
            printname = [foldername 'rcae.mat'];

            for l = {'lo', 'l', 'o'}; land = l{1};
                for fn = fieldnames(flux)'; fname = fn{1};
                    if strcmp(land, 'lo'); flux_n.(land).(fname) = flux.(fname);
                    elseif strcmp(land, 'l'); flux_n.(land).(fname) = flux.(fname) .*mask.ocean;
                    elseif strcmp(land, 'o'); flux_n.(land).(fname) = flux.(fname) .*mask.land;
                    end
                end

                % lon x lat structure over various time averages
                for t = {'ann', 'djf', 'jja', 'mam', 'son'}; time = t{1};
                    for fn = fieldnames(flux)'; fname = fn{1};
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
                    end
                    rcae_t.(land).(time) = def_rcae(type, flux_t.(land).(time), par);
                end
            end

            printname_t = [foldername 'rcae_t.mat'];
        elseif strcmp(ftype, 'flux_z') % show lat x mon structure of RCAE
            rcae_z = def_rcae(type, flux_z, par);
            printname = [foldername 'rcae_z.mat'];
        end

        % save rcae data
        if ~exist(foldername, 'dir')
            mkdir(foldername)
        end
        if strcmp(ftype, 'flux')
            save(printname, 'rcae', 'lat');
            save(printname_t, 'rcae_t', 'lat');
        elseif strcmp(ftype, 'flux_z')
            save(printname, 'rcae_z', 'lat');
        end
    end

end
function proc_temp(type, par)
    if strcmp(type, 'era5') | strcmp(type, 'erai')
        foldername = sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/eps_%g/', type, par.lat_interp, par.ep);
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s', type);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s', type);
        file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/%s/temp/%s_temp_*.ymonmean.nc', type, type));
        fullpath=sprintf('%s/%s', file.folder, file.name);
        temp = ncread(fullpath, 't');
    elseif strcmp(type, 'gcm')
        foldername = sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s/eps_%g/', type, par.model, par.lat_interp, par.ep);
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s', type, par.model);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s', type, par.model);
        var = 'ta';
        file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/gcm/%s/%s_Amon_%s_piControl_r1i1p1_*.ymonmean.nc', par.model, var, par.model));
        fullpath=sprintf('%s/%s', file.folder, file.name);
        temp = ncread(fullpath, var);
    end
    load(sprintf('%s/grid.mat', prefix)); % read ERA5 grid data
    load(sprintf('%s/srfc.mat', prefix)); % load surface data
    load(sprintf('%s/%s/masks.mat', prefix_proc, par.lat_interp)); % load land and ocean masks
    load(sprintf('%s/%s/eps_%g/rcae_t.mat', prefix_proc, par.lat_interp, par.ep)); % load rcae data

    lat = par.lat_std;

    % interpolate temp to standard lat grid
    temp = permute(temp, [2 1 3 4]);
    temp = interp1(grid.dim3.lat, temp, lat);
    temp = permute(temp, [2 1 3 4]);

    % create surface mask
    if strcmp(type, 'era5') | strcmp(type, 'erai')
        ps = permute(srfc.sp, [2 1 3]); % bring lat to front to interpolate
        ps = interp1(grid.dim2.lat, ps, lat); % interpolate to standard grid
        ps = permute(ps, [2 1 3]); % reorder to original
        ps_vert = repmat(ps, [1 1 1 size(temp, 3)]); % dims (lon x lat x time x plev)
        ps_vert = permute(ps_vert, [1 2 4 3]); % dims (lon x lat x plev x time)
        pa = double(permute(repmat(grid.dim3.plev, [1 size(sp)]), [2 3 1 4]))*10^2; % convert from hPa to Pa
    elseif strcmp(type, 'gcm')
        ps = permute(srfc.ps, [2 1 3]); % bring lat to front to interpolate
        ps = interp1(grid.dim2.lat, ps, lat); % interpolate to standard grid
        ps = permute(ps, [2 1 3]); % reorder to original
        ps_vert = repmat(ps, [1 1 1 size(temp, 3)]); % dims (lon x lat x time x plev)
        ps_vert = permute(ps_vert, [1 2 4 3]); % dims (lon x lat x plev x time)
        pa = permute(repmat(grid.dim3.plev, [1 size(ps)]), [2 3 1 4]); % convert from hPa to Pa
    end
    surface_mask = nan(size(temp));
    surface_mask(pa < ps_vert) = 1;

    % apply masks to surface pressure and RCAE regimes
    ps_n.lo = ps;
    ps_n.l = ps .* mask.ocean;
    ps_n.o = ps .* mask.land;

    mask.land_vert = repmat(mask.land, [1 1 1 size(temp, 3)]); % expand land mask to vertical dim
    mask.land_vert = permute(mask.land, [1 2 4 3]); % place vertical dim where it belongs
    mask.ocean_vert = repmat(mask.ocean, [1 1 1 size(temp, 3)]); % expand ocean mask to vertical dim
    mask.ocean_vert = permute(mask.ocean, [1 2 4 3]); % place vertical dim where it belongs

    temp_sm.lo = temp.*surface_mask; % filter temp with surface mask
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
    for fn = fieldnames(rcae_t.lo.ann)'
        for l = {'lo', 'l', 'o'}; land = l{1}; % over land, over ocean, or both
            for t = {'ann', 'djf', 'jja', 'mam', 'son'}; time = t{1};
                for re = {'rce', 'rae'}; regime = re{1};
                    if strcmp(time, 'ann')
                        temp_n.(land).(time) = squeeze(nanmean(temp_sm.(land), 3));
                    elseif strcmp(time, 'djf')
                        temp_shift = circshift(temp_sm.(land), 1, 3);
                        temp_n.(land).(time) = squeeze(nanmean(temp_shift(:,:,1:3,:), 3));
                    elseif strcmp(time, 'jja')
                        temp_n.(land).(time) = squeeze(nanmean(temp_sm.(land)(:,:,6:8,:), 3));
                    elseif strcmp(time, 'mam')
                        temp_n.(land).(time) = squeeze(nanmean(temp_sm.(land)(:,:,3:5,:), 3));
                    elseif strcmp(time, 'son')
                        temp_n.(land).(time) = squeeze(nanmean(temp_sm.(land)(:,:,9:11,:), 3));
                    end

                    filt.(land).(time).(fn{1}).(regime) = nan(size(rcae_t.(land).(time).(fn{1}))); % create empty arrays to store filtering array
                    if strcmp(regime, 'rce'); filt.(land).(time).(fn{1}).(regime)(rcae_t.(land).(time).(fn{1})==1)=1; % set RCE=1, elsewhere nan
                    elseif strcmp(regime, 'rae'); filt.(land).(time).(fn{1}).(regime)(rcae_t.(land).(time).(fn{1})==-1)=1; % set RAE=1, elsewhere nan
                    end

                    temp_t.(land).(time).(fn{1}).(regime) = temp_n.(land).(time) .* filt.(land).(time).(fn{1}).(regime);

                    nanfilt.(regime) = nan(size(temp_t.(land).(time).(fn{1}).(regime)));
                    nanfilt.(regime)(~isnan(temp_t.(land).(time).(fn{1}).(regime))) = 1;

                    % take cosine-weighted meridional average
                    for d = {'all', 'nh', 'sh', 'tp'}; domain = d{1};
                        if strcmp(regime, 'rce')
                            if strcmp(domain, 'all')
                                cosw = repmat(cosd(lat)', [size(nanfilt.(regime), 1) 1]);
                                denm = temp_t.(land).(time).(fn{1}).(regime);
                                nume = nanfilt.(regime);
                            elseif strcmp(domain, 'nh')
                                cosw = repmat(cosd(lat(lat>30))', [size(nanfilt.(regime), 1) 1]);
                                denm = temp_t.(land).(time).(fn{1}).(regime)(:,lat>30,:);
                                nume = nanfilt.(regime)(:,lat>30,:);
                            elseif strcmp(domain, 'sh')
                                cosw = repmat(cosd(lat(lat<-30))', [size(nanfilt.(regime), 1) 1]);
                                denm = temp_t.(land).(time).(fn{1}).(regime)(:,lat<-30,:);
                                nume = nanfilt.(regime)(:,lat<-30,:);
                            elseif strcmp(domain, 'tp')
                                cosw = repmat(cosd(lat(abs(lat)<30))', [size(nanfilt.(regime), 1) 1]);
                                denm = temp_t.(land).(time).(fn{1}).(regime)(:,abs(lat)<30,:);
                                nume = nanfilt.(regime)(:,abs(lat)<30,:);
                            end
                        elseif strcmp(regime, 'rae')
                            if strcmp(domain, 'all')
                                cosw = repmat(cosd(lat)', [size(nanfilt.(regime), 1) 1]);
                                denm = temp_t.(land).(time).(fn{1}).(regime);
                                nume = nanfilt.(regime);
                            elseif strcmp(domain, 'nh')
                                cosw = repmat(cosd(lat(lat>0))', [size(nanfilt.(regime), 1) 1]);
                                denm = temp_t.(land).(time).(fn{1}).(regime)(:,lat>0,:);
                                nume = nanfilt.(regime)(:,lat>0,:);
                            elseif strcmp(domain, 'sh')
                                cosw = repmat(cosd(lat(lat<-0))', [size(nanfilt.(regime), 1) 1]);
                                denm = temp_t.(land).(time).(fn{1}).(regime)(:,lat<-0,:);
                                nume = nanfilt.(regime)(:,lat<-0,:);
                            end
                        end

                        ta.(regime).(domain).(fn{1}).(land).(time) = nansum(cosw.*denm, 2) ./ nansum(cosw.*nume, 2); % area-weighted meridional average
                        ta.(regime).(domain).(fn{1}).(land).(time) = squeeze(nanmean(ta.(regime).(domain).(fn{1}).(land).(time), 1));
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
    save(printname, 'ta');
end
function proc_temp_mon_lat(type, par)
    if strcmp(type, 'era5') | strcmp(type, 'erai')
        foldername = sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/eps_%g/', type, par.lat_interp, par.ep);
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s', type);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s', type);
        file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/%s/temp/%s_temp_*.ymonmean.nc', type, type));
        fullpath=sprintf('%s/%s', file.folder, file.name);
        temp = ncread(fullpath, 't');
    elseif strcmp(type, 'gcm')
        foldername = sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s/eps_%g/', type, par.model, par.lat_interp, par.ep);
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s', type, par.model);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s', type, par.model);
        var = 'ta';
        file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/gcm/%s/%s_Amon_%s_piControl_r1i1p1_*.ymonmean.nc', par.model, var, par.model));
        fullpath=sprintf('%s/%s', file.folder, file.name);
        temp = ncread(fullpath, var);
    end
    load(sprintf('%s/grid.mat', prefix)); % read ERA5 grid data
    load(sprintf('%s/srfc.mat', prefix)); % load surface data
    load(sprintf('%s/%s/masks.mat', prefix_proc, par.lat_interp)); % load land and ocean masks

    lat = par.lat_std;

    % interpolate temp to standard lat grid
    temp = permute(temp, [2 1 3 4]);
    temp = interp1(grid.dim3.lat, temp, lat);
    temp = permute(temp, [2 1 3 4]);

    % create surface mask
    if strcmp(type, 'era5') | strcmp(type, 'erai')
        ps = permute(srfc.sp, [2 1 3]); % bring lat to front to interpolate
        ps = interp1(grid.dim2.lat, ps, lat); % interpolate to standard grid
        ps = permute(ps, [2 1 3]); % reorder to original
        ps_vert = repmat(ps, [1 1 1 size(temp, 3)]); % dims (lon x lat x time x plev)
        ps_vert = permute(ps_vert, [1 2 4 3]); % dims (lon x lat x plev x time)
        pa = double(permute(repmat(grid.dim3.plev, [1 size(sp)]), [2 3 1 4]))*10^2; % convert from hPa to Pa
    elseif strcmp(type, 'gcm')
        ps = permute(srfc.ps, [2 1 3]); % bring lat to front to interpolate
        ps = interp1(grid.dim2.lat, ps, lat); % interpolate to standard grid
        ps = permute(ps, [2 1 3]); % reorder to original
        ps_vert = repmat(ps, [1 1 1 size(temp, 3)]); % dims (lon x lat x time x plev)
        ps_vert = permute(ps_vert, [1 2 4 3]); % dims (lon x lat x plev x time)
        pa = permute(repmat(grid.dim3.plev, [1 size(ps)]), [2 3 1 4]); % convert from hPa to Pa
    end
    surface_mask = nan(size(temp));
    surface_mask(pa < ps_vert) = 1;

    mask.land_vert = repmat(mask.land, [1 1 1 size(temp, 3)]); % expand land mask to vertical dim
    mask.land_vert = permute(mask.land, [1 2 4 3]); % place vertical dim where it belongs
    mask.ocean_vert = repmat(mask.ocean, [1 1 1 size(temp, 3)]); % expand ocean mask to vertical dim
    mask.ocean_vert = permute(mask.ocean, [1 2 4 3]); % place vertical dim where it belongs

    temp_sm.lo = temp.*surface_mask; % filter temp with surface mask
    temp_sm.l = temp.*surface_mask.*mask.ocean_vert; % filter temp with surface mask
    temp_sm.o = temp.*surface_mask.*mask.land_vert; % filter temp with surface mask
    temp_sm.lo = permute(temp_sm.lo, [1 2 4 3]); % bring plev to last dimension
    temp_sm.l = permute(temp_sm.l, [1 2 4 3]); % bring plev to last dimension
    temp_sm.o = permute(temp_sm.o, [1 2 4 3]); % bring plev to last dimension

    mask_t.land = nanmean(mask.land, 3);
    mask_t.ocean = nanmean(mask.ocean, 3);

    for l = {'lo', 'l', 'o'}; land = l{1}; % over land, over ocean, or both
        ta.(land)= squeeze(nanmean(temp_sm.(land), 1)); % zonal average
    end

    % save filtered data
    printname = [foldername 'ta_mon_lat'];
    if ~exist(foldername, 'dir')
        mkdir(foldername)
    end
    save(printname, 'ta');
end
function proc_ma(type, par)
% calculate moist adiabats
    if strcmp(type, 'era5') | strcmp(type, 'erai')
        foldername = sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/eps_%g/', type, par.lat_interp, par.ep);
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s', type);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s', type);
    elseif strcmp(type, 'gcm')
        foldername = sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s/eps_%g/', type, par.model, par.lat_interp, par.ep);
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s', type, par.model);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s', type, par.model);
    end
    load(sprintf('%s/grid.mat', prefix)); % read grid data
    load(sprintf('%s/srfc.mat', prefix)); % read surface variable data
    load(sprintf('%s/%s/masks.mat', prefix_proc, par.lat_interp)); % load land and ocean masks
    load(sprintf('%s/%s/eps_%g/rcae_t.mat', prefix_proc, par.lat_interp, par.ep)); % load rcae data

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

            % srfc_z.(fn{1}).(land) = squeeze(nanmean(srfc_n.(fn{1}).(land), 1)); % zonal mean
            % srfc_zi.(fn{1}).(land) = interp1(grid.dim2.lat, srfc_z.(fn{1}).(land), lat); % interpolate to standard lat grid
        end
    end

    for fn = fieldnames(rcae_t.lo.ann)'
        for l = {'lo', 'l', 'o'}; land = l{1};
            for t = {'ann', 'djf', 'jja', 'mam', 'son'}; time = t{1};
                for re = {'rce', 'rae'}; regime = re{1};
                    for fn_var = fieldnames(srfc)'
                        if strcmp(time, 'ann')
                            srfc_t.(land).(time).(fn_var{1}) = squeeze(nanmean(srfc_n.(fn_var{1}).(land), 3));
                        elseif strcmp(time, 'djf')
                            srfc_shift = circshift(srfc_n.(fn_var{1}).(land), 1, 3);
                            srfc_t.(land).(time).(fn_var{1}) = squeeze(nanmean(srfc_shift(:,:,1:3), 3));
                        elseif strcmp(time, 'jja')
                            srfc_t.(land).(time).(fn_var{1}) = squeeze(nanmean(srfc_n.(fn_var{1}).(land)(:,:,6:8), 3));
                        elseif strcmp(time, 'mam')
                            srfc_t.(land).(time).(fn_var{1}) = squeeze(nanmean(srfc_n.(fn_var{1}).(land)(:,:,3:5), 3));
                        elseif strcmp(time, 'son')
                            srfc_t.(land).(time).(fn_var{1}) = squeeze(nanmean(srfc_n.(fn_var{1}).(land)(:,:,9:11), 3));
                        end

                        filt.(land).(time).(fn{1}).(regime) = nan(size(rcae_t.(land).(time).(fn{1})));
                        if strcmp(regime, 'rce'); filt.(land).(time).(fn{1}).(regime)(rcae_t.(land).(time).(fn{1})==1)=1; % set RCE=1, elsewhere nan
                        elseif strcmp(regime, 'rae'); filt.(land).(time).(fn{1}).(regime)(rcae_t.(land).(time).(fn{1})==-1)=1; % set RAE=1, elsewhere nan
                        end

                        srfc_tf.(land).(time).(fn{1}).(regime).(fn_var{1}) = srfc_t.(land).(time).(fn_var{1}) .* filt.(land).(time).(fn{1}).(regime);

                        nanfilt.(regime) = nan(size(srfc_tf.(land).(time).(fn{1}).(regime).(fn_var{1})));
                        nanfilt.(regime)(~isnan(srfc_tf.(land).(time).(fn{1}).(regime).(fn_var{1}))) = 1;
                      
                        % take cosine-weighted average
                        for d = {'all', 'nh', 'sh', 'tp'}; domain = d{1};
                            if strcmp(regime, 'rce')
                                if strcmp(domain, 'all')
                                    cosw = repmat(cosd(lat)', [size(nanfilt.(regime), 1) 1]);
                                    denm = srfc_tf.(land).(time).(fn{1}).(regime).(fn_var{1});
                                    nume = nanfilt.(regime);
                                elseif strcmp(domain, 'nh')
                                    cosw = repmat(cosd(lat(lat>30))', [size(nanfilt.(regime), 1) 1]);
                                    denm = srfc_tf.(land).(time).(fn{1}).(regime).(fn_var{1})(:, lat>30);
                                    nume = nanfilt.(regime)(:, lat>30);
                                elseif strcmp(domain, 'sh')
                                    cosw = repmat(cosd(lat(lat<-30))', [size(nanfilt.(regime), 1) 1]);
                                    denm = srfc_tf.(land).(time).(fn{1}).(regime).(fn_var{1})(:, lat<-30);
                                    nume = nanfilt.(regime)(:, lat<-30);
                                elseif strcmp(domain, 'tp')
                                    cosw = repmat(cosd(lat(abs(lat)<30))', [size(nanfilt.(regime), 1) 1]);
                                    denm = srfc_tf.(land).(time).(fn{1}).(regime).(fn_var{1})(:, abs(lat)<30);
                                    nume = nanfilt.(regime)(:, abs(lat)<30);
                                end
                            elseif strcmp(regime, 'rae')
                                if strcmp(domain, 'all')
                                    cosw = repmat(cosd(lat)', [size(nanfilt.(regime), 1) 1]);
                                    denm = srfc_tf.(land).(time).(fn{1}).(regime).(fn_var{1});
                                    nume = nanfilt.(regime);
                                elseif strcmp(domain, 'nh')
                                    cosw = repmat(cosd(lat(lat>0))', [size(nanfilt.(regime), 1) 1]);
                                    denm = srfc_tf.(land).(time).(fn{1}).(regime).(fn_var{1})(:, lat>0);
                                    nume = nanfilt.(regime)(:, lat>0);
                                elseif strcmp(domain, 'sh')
                                    cosw = repmat(cosd(lat(lat<0))', [size(nanfilt.(regime), 1) 1]);
                                    denm = srfc_tf.(land).(time).(fn{1}).(regime).(fn_var{1})(:, lat<0);
                                    nume = nanfilt.(regime)(:, lat<0);
                                end
                            end

                            ma.(regime).(domain).(fn{1}).(land).(time).(fn_var{1}) = nansum(cosw.*denm, 2) ./ nansum(cosw.*nume, 2); % weighted meridional average
                            ma.(regime).(domain).(fn{1}).(land).(time).(fn_var{1}) = squeeze(nanmean(ma.(regime).(domain).(fn{1}).(land).(time).(fn_var{1}), 1)); % zonal average

                        end % end domain loop
                    end % end srfc variables loop

                    if strcmp(type, 'era5') | strcmp(type, 'erai')
                        % ma = calc_ma_dew(ma, grid.dim3.plev, par); % compute moist adiabat with dew point temperature
                    elseif strcmp(type, 'gcm')
                        ma.rce.all.(fn{1}).(land).(time).ta = calc_ma_hurs(ma.rce.all.(fn{1}).(land).(time), grid.dim3.plev, par); % compute moist adiabat with RH
                        ma.rce.tp.(fn{1}).(land).(time).ta = calc_ma_hurs(ma.rce.tp.(fn{1}).(land).(time), grid.dim3.plev, par); % compute moist adiabat with RH
                        ma.rce.nh.(fn{1}).(land).(time).ta = calc_ma_hurs(ma.rce.nh.(fn{1}).(land).(time), grid.dim3.plev, par); % compute moist adiabat with RH
                        ma.rce.sh.(fn{1}).(land).(time).ta = calc_ma_hurs(ma.rce.sh.(fn{1}).(land).(time), grid.dim3.plev, par); % compute moist adiabat with RH
                    end

                end % end RCE/RAE regime loop
            end % end time average loop
        end % end land option loop
    end % end RCAE definition loop


    % for fn = fieldnames(rcae)'
    %     rce_filt.(fn{1}) = nan(size(rcae.(fn{1}))); % create empty arrays to store filtering array
    %     rae_filt.(fn{1}) = nan(size(rcae.(fn{1}))); % dims (lon x lat x time)
    %     rce_filt.(fn{1})(rcae.(fn{1})==1)=1; % set RCE=1, elsewhere nan
    %     rae_filt.(fn{1})(rcae.(fn{1})==-1)=1; % set RAE=1, elsewhere nan

    %     for l = {'lo', 'l', 'o'}; land = l{1}; % over land, ocean, or both?

    %         for i = {'ann', 'djf', 'jja', 'mam', 'son'}; time = i{1};

    %             for fn_var = fieldnames(srfc)'

    %                 srfc_zi.rce.(fn{1}).(land).(fn_var{1}) = rce_filt.(fn{1}) .* srfc_zi.(fn_var{1}).(land); % RCE-filtered values
    %                 srfc_zi.rae.(fn{1}).(land).(fn_var{1}) = rae_filt.(fn{1}) .* srfc_zi.(fn_var{1}).(land); % RAE-filtered values
    %                 if strcmp(time, 'ann')
    %                     rce_filt_t.(fn{1}).(land) = squeeze(nanmean(rce_filt.(fn{1}), 2)); % take time mean of RCE filter
    %                     rae_filt_t.(fn{1}).(land) = squeeze(nanmean(rae_filt.(fn{1}), 2)); % take time mean of RAE filter
    %                     srfc_zit.rce.(fn{1}).(land).(time).(fn_var{1}) = squeeze(nanmean(srfc_zi.rce.(fn{1}).(land).(fn_var{1}), 2)); % take time mean of temperature
    %                     srfc_zit.rae.(fn{1}).(land).(time).(fn_var{1}) = squeeze(nanmean(srfc_zi.rae.(fn{1}).(land).(fn_var{1}), 2)); % take time mean of temperature
    %                 elseif strcmp(time, 'djf')
    %                     rce_filt_shift.(fn{1}).(land) = circshift(rce_filt.(fn{1}), 1, 2); % shift months by one so I can take DJF average in one go (1:3)
    %                     rae_filt_shift.(fn{1}).(land) = circshift(rae_filt.(fn{1}), 1, 2);
    %                     srfc_zi_shift.rce.(fn{1}).(land).(fn_var{1}) = circshift(srfc_zi.rce.(fn{1}).(land).(fn_var{1}), 1, 2);
    %                     srfc_zi_shift.rae.(fn{1}).(land).(fn_var{1}) = circshift(srfc_zi.rae.(fn{1}).(land).(fn_var{1}), 1, 2);
    %                     rce_filt_t.(fn{1}).(land) = squeeze(nanmean(rce_filt_shift.(fn{1}).(land)(:,1:3,:), 2)); % take time mean of RCE filter
    %                     rae_filt_t.(fn{1}).(land) = squeeze(nanmean(rae_filt_shift.(fn{1}).(land)(:,1:3,:), 2)); % take time mean of RAE filter
    %                     srfc_zit.rce.(fn{1}).(land).(time).(fn_var{1}) = squeeze(nanmean(srfc_zi_shift.rce.(fn{1}).(land).(fn_var{1})(:,1:3,:), 2)); % take time mean of temperature
    %                     srfc_zit.rae.(fn{1}).(land).(time).(fn_var{1}) = squeeze(nanmean(srfc_zi_shift.rae.(fn{1}).(land).(fn_var{1})(:,1:3,:), 2)); % take time mean of temperature
    %                 elseif strcmp(time, 'jja')
    %                     rce_filt_t.(fn{1}).(land) = squeeze(nanmean(rce_filt.(fn{1})(:,6:8,:), 2)); % take time mean of RCE filter
    %                     rae_filt_t.(fn{1}).(land) = squeeze(nanmean(rae_filt.(fn{1})(:,6:8,:), 2)); % take time mean of RAE filter
    %                     srfc_zit.rce.(fn{1}).(land).(time).(fn_var{1}) = squeeze(nanmean(srfc_zi.rce.(fn{1}).(land).(fn_var{1})(:,6:8,:), 2)); % take time mean of temperature
    %                     srfc_zit.rae.(fn{1}).(land).(time).(fn_var{1}) = squeeze(nanmean(srfc_zi.rae.(fn{1}).(land).(fn_var{1})(:,6:8,:), 2)); % take time mean of temperature
    %                 elseif strcmp(time, 'mam')
    %                     rce_filt_t.(fn{1}).(land) = squeeze(nanmean(rce_filt.(fn{1})(:,3:5,:), 2)); % take time mean of RCE filter
    %                     rae_filt_t.(fn{1}).(land) = squeeze(nanmean(rae_filt.(fn{1})(:,3:5,:), 2)); % take time mean of RAE filter
    %                     srfc_zit.rce.(fn{1}).(land).(time).(fn_var{1}) = squeeze(nanmean(srfc_zi.rce.(fn{1}).(land).(fn_var{1})(:,3:5,:), 2)); % take time mean of temperature
    %                     srfc_zit.rae.(fn{1}).(land).(time).(fn_var{1}) = squeeze(nanmean(srfc_zi.rae.(fn{1}).(land).(fn_var{1})(:,3:5,:), 2)); % take time mean of temperature
    %                 elseif strcmp(time, 'son')
    %                     rce_filt_t.(fn{1}).(land) = squeeze(nanmean(rce_filt.(fn{1})(:,9:11,:), 2)); % take time mean of RCE filter
    %                     rae_filt_t.(fn{1}).(land) = squeeze(nanmean(rae_filt.(fn{1})(:,9:11,:), 2)); % take time mean of RAE filter
    %                     srfc_zit.rce.(fn{1}).(land).(time).(fn_var{1}) = squeeze(nanmean(srfc_zi.rce.(fn{1}).(land).(fn_var{1})(:,9:11,:), 2)); % take time mean of temperature
    %                     srfc_zit.rae.(fn{1}).(land).(time).(fn_var{1}) = squeeze(nanmean(srfc_zi.rae.(fn{1}).(land).(fn_var{1})(:,9:11,:), 2)); % take time mean of temperature
    %                 end

    %                 ma.rce.(fn{1}).(land).(time).(fn_var{1}) = nansum(cosd(lat).*srfc_zit.rce.(fn{1}).(land).(time).(fn_var{1})) / nansum(cosd(lat).*rce_filt_t.(fn{1}).(land)); % area-weighted meridional average
    %                 ma.rce.tp.(fn{1}).(land).(time).(fn_var{1}) = nansum(cosd(lat(abs(lat)<30)).*srfc_zit.rce.(fn{1}).(land).(time).(fn_var{1})(abs(lat)<30,:)) / nansum(cosd(lat(abs(lat)<30)).*rce_filt_t.(fn{1}).(land)(abs(lat)<30)); % area-weighted meridional average
    %                 ma.rce.nh.(fn{1}).(land).(time).(fn_var{1}) = nansum(cosd(lat(lat>30)).*srfc_zit.rce.(fn{1}).(land).(time).(fn_var{1})(lat>30,:)) / nansum(cosd(lat(lat>30)).*rce_filt_t.(fn{1}).(land)(lat>30)); % area-weighted meridional average
    %                 ma.rce.sh.(fn{1}).(land).(time).(fn_var{1}) = nansum(cosd(lat(lat<-30)).*srfc_zit.rce.(fn{1}).(land).(time).(fn_var{1})(lat<-30,:)) / nansum(cosd(lat(lat<-30)).*rce_filt_t.(fn{1}).(land)(lat<-30)); % area-weighted meridional average
    %                 ma.rae.(fn{1}).(land).(time).(fn_var{1}) = nansum(cosd(lat).*srfc_zit.rae.(fn{1}).(land).(time).(fn_var{1})) / nansum(cosd(lat).*rae_filt_t.(fn{1}).(land)); % area-weighted meridional average
    %                 ma.rae.nh.(fn{1}).(land).(time).(fn_var{1}) = nansum(cosd(lat(lat>0)).*srfc_zit.rae.(fn{1}).(land).(time).(fn_var{1})(lat>0,:)) / nansum(cosd(lat(lat>0)).*rae_filt_t.(fn{1}).(land)(lat>0)); % area-weighted meridional average
    %                 ma.rae.sh.(fn{1}).(land).(time).(fn_var{1}) = nansum(cosd(lat(lat<0)).*srfc_zit.rae.(fn{1}).(land).(time).(fn_var{1})(lat<0,:)) / nansum(cosd(lat(lat<0)).*rae_filt_t.(fn{1}).(land)(lat<0)); % area-weighted meridional average

    %             end

    %             if strcmp(type, 'era5') | strcmp(type, 'erai')
    %                 ma = calc_ma_dew(ma, grid.dim3.plev, par); % compute moist adiabat with dew point temperature
    %             elseif strcmp(type, 'gcm')
    %                 ma.rce.(fn{1}).(land).(time).ta = calc_ma_hurs(ma.rce.(fn{1}).(land).(time), grid.dim3.plev, par); % compute moist adiabat with RH
    %                 ma.rce.tp.(fn{1}).(land).(time).ta = calc_ma_hurs(ma.rce.tp.(fn{1}).(land).(time), grid.dim3.plev, par); % compute moist adiabat with RH
    %                 ma.rce.nh.(fn{1}).(land).(time).ta = calc_ma_hurs(ma.rce.nh.(fn{1}).(land).(time), grid.dim3.plev, par); % compute moist adiabat with RH
    %                 ma.rce.sh.(fn{1}).(land).(time).ta = calc_ma_hurs(ma.rce.sh.(fn{1}).(land).(time), grid.dim3.plev, par); % compute moist adiabat with RH
    %             end

    %         end
    %     end
    % end

    % save data into mat file
    printname = [foldername 'ma.mat'];
    if ~exist(foldername, 'dir')
        mkdir(foldername)
    end
    save(printname, 'ma');

end
function proc_ma_mon_lat(type, par)
% calculate moist adiabats
    if strcmp(type, 'era5') | strcmp(type, 'erai')
        foldername = sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/eps_%g/', type, par.lat_interp, par.ep);
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s', type);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s', type);
    elseif strcmp(type, 'gcm')
        foldername = sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s/eps_%g/', type, par.model, par.lat_interp, par.ep);
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s', type, par.model);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s', type, par.model);
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

        end % end srfc variables loop

        if strcmp(type, 'era5') | strcmp(type, 'erai')
            % ma = calc_ma_dew(ma, grid.dim3.plev, par); % compute moist adiabat with dew point temperature
        elseif strcmp(type, 'gcm')
            pb = CmdLineProgressBar("Calculating moist adiabats...");
            for ilat = 1:length(lat);
                pb.print(ilat, length(lat)); % output progress of moist adiabat calculation
                for imon = 1:12;
                    for fn_var = fieldnames(srfc)'
                        ima.(fn_var{1}) = ma.(land).(fn_var{1})(ilat, imon);
                    end
                    ma.(land).ta(ilat, imon, :) = calc_ma_hurs(ima, grid.dim3.plev, par); % compute moist adiabat with RH
                end
            end
        end
    end % end land option loop

    % save data into mat file
    printname = [foldername 'ma_mon_lat.mat'];
    if ~exist(foldername, 'dir')
        mkdir(foldername)
    end
    save(printname, 'ma');

end

% helper functions
function rcae = def_rcae(type, flux, par)
    % identify locations of RCE using threshold epsilon (ep)
    rcae.def = zeros(size(flux.r1));
    rcae.def(abs(flux.r1) < par.ep) = 1;
    % identify locations of RAE using threshold gamma (ga)
    rcae.def(flux.r1 > par.ga) = -1;

    % add additional criteria for RCE that P-E>0
    rcae.pe = zeros(size(flux.r1));
    if strcmp(type, 'era5') | strcmp(type, 'erai')
        rcae.pe(abs(flux.r1) < par.ep & (flux.cp+flux.lsp+flux.e>0)) = 1; % note that evap is defined negative into atmosphere
    elseif strcmp(type, 'gcm')
        rcae.pe(abs(flux.r1) < par.ep & (flux.pr-flux.evspsbl>0)) = 1;
    end
    rcae.pe(flux.r1 > par.ga) = -1;

    % add additional criteria for RCE that (P_ls - E)<<1
    rcae.cp = zeros(size(flux.r1));
    if strcmp(type, 'era5') | strcmp(type, 'erai')
        rcae.cp(abs(flux.r1) < par.ep & abs(flux.lsp./flux.cp<par.ep_cp)) = 1; % note that evap is defined negative into atmosphere
    elseif strcmp(type, 'gcm')
        rcae.cp(abs(flux.r1) < par.ep & abs((flux.pr-flux.prc)./flux.prc<par.ep_cp)) = 1; % note that evap is defined negative into atmosphere
    end
    rcae.cp(flux.r1 > par.ga) = -1;
end
function land_mask = remove_land(lat, lon, nt)

	[lat2dgrid, lon2dgrid] = meshgrid(lat, lon - 180);

	land_mask =~ circshift(landmask(lat2dgrid, lon2dgrid), length(lon)/2, 1);
	land_mask = double(land_mask);
	land_mask(land_mask==0) = nan;
	land_mask = repmat(land_mask, [1 1 nt]); % land mask in dims (lon x lat x time)

end
function ocean_mask = remove_ocean(lat, lon, nt)

	[lat2dgrid, lon2dgrid] = meshgrid(lat, lon - 180);

	ocean_mask = circshift(landmask(lat2dgrid, lon2dgrid), length(lon)/2, 1);
	ocean_mask = double(ocean_mask);
	ocean_mask(ocean_mask==0) = nan;
	ocean_mask = repmat(ocean_mask, [1 1 nt]); % ocean mask in dims (lon x lat x time)

end

function ma = calc_ma_dew(ma_in, plev, par)
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

        ta_s = interp1(pa_s, ta_s, plev, 'spline', nan);
        qsat_s = interp1(pa_s, qsat_s, plev, 'spline', nan);
        dtadpa_s = interp1(pa_s, dtadpa_s, plev, 'spline', nan);

        % OUTPUT
        ma_out = ta_s;
    end

end

