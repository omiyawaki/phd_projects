clc; close all; clear variables;

addpath(genpath('/project2/tas1/miyawaki/matlab'));

%% set parameters
par.lat_interp = 'std'; % which latitudinal grid to interpolate to: don (donohoe, coarse), era (native ERA-Interim, fine), or std (custom, very fine)
par.lat_std = transpose(-90:0.25:90); % define standard latitude grid for 'std' interpolation
par.ep_swp = [0.25 0.3 0.35]; % threshold for RCE definition. RCE is defined as where abs(R1) < ep
par.ep_cp = 0.5; % threshold for RCE definition. RCE is defined as where abs(R1) < ep
par.ma_type = 'std'; % choose the type of moist adiabat: reversible, pseudo, or std
par.frz = 0; % consider latent heat of fusion in moist adiabat?
par.pa_span = [1000 100]*100; % pressure range for calculating moist adiabat
par.dpa = -10; % pressure increment for integrating dry adiabat section of moist adiabat (matters for how accurately the LCL is computed)
par.cpd = 1005.7; par.Rd = 287; par.Rv = 461; par.g = 9.81; par.L = 2.501e6; par.a = 6357e3; par.eps = 0.622; % common constants, all in SI units for brevity
gcm_info

%% call functions
% comp_flux(par)
disp_global_flux(par)

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
        % choose_proc_ep(type, par)
    end
end

%% define functions
function choose_proc(type, par)
    proc_flux(type, par) % calculate energy fluxes in the vertically-integrated MSE budget using ERA-Interim data
    % proc_net_flux(type, par) % calculate net energy fluxes at TOA and surface
    % save_mask(type, par) % save land and ocean masks once (faster than creating mask every time I need it)
end
function proc_flux(type, par)
    if strcmp(type, 'era5')
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s', type);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s', type, par.lat_interp);
    elseif strcmp(type, 'gcm')
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s', type, par.model);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s', type, par.model, par.lat_interp);
    end

    load(sprintf('%s/grid.mat', prefix)) % read grid data
    load(sprintf('%s/rad.mat', prefix)) % read radiation data
    load(sprintf('%s/pe.mat', prefix)) % read hydrology data
    load(sprintf('%s/stf.mat', prefix)) % read surface turbulent flux data
    load(sprintf('%s/masks.mat', prefix_proc)); % load land and ocean masks

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
        flux.stf.mse = -( flux.sshf + flux.slhf );
        flux.stf.dse = par.L*(flux.cp+flux.lsp) - flux.slhf;
    elseif strcmp(type, 'gcm')
        flux.ra = flux.rsdt - flux.rsut + flux.rsus - flux.rsds + flux.rlus - flux.rlds - flux.rlut; % calculate atmospheric radiative cooling
        flux.stf.mse = flux.hfls + flux.hfss;
        flux.stf.dse = par.L*flux.pr + flux.hfss;
    end

    for f = {'mse', 'dse'}; fw = f{1};
        flux.res.(fw) = flux.ra + flux.stf.(fw); % infer MSE tendency and flux divergence as residuals
        % compute northward MSE transport using the residual data
        for t = {'ann', 'djf', 'jja', 'mam', 'son'}; time = t{1};
            if strcmp(time, 'ann')
                tediv_t = squeeze(nanmean(flux.res.(fw), 3)); % take time average
            elseif strcmp(time, 'djf')
                flux_shift.res = circshift(flux.res.(fw), 1, 3); % shift month by 1 to allow for djf average
                tediv_t = squeeze(nanmean(flux_shift.res(:,:,1:3), 3)); % take time average
            elseif strcmp(time, 'jja')
                tediv_t = squeeze(nanmean(flux.res.(fw)(:,:,6:8), 3)); % take time average
            elseif strcmp(time, 'mam')
                tediv_t = squeeze(nanmean(flux.res.(fw)(:,:,3:5), 3)); % take time average
            elseif strcmp(time, 'son')
                tediv_t = squeeze(nanmean(flux.res.(fw)(:,:,9:11), 3)); % take time average
            end
            tediv_tz = trapz(deg2rad(grid.dim2.lon), tediv_t, 1); % zonally integrate
            vh.(time).(fw) = cumtrapz(deg2rad(lat), par.a^2*cosd(lat).*tediv_tz'); % cumulatively integrate in latitude
        end
        flux.r1.(fw) = (flux.res.(fw))./flux.ra; % calculate nondimensional number R1 disregarding MSE budget closure
        flux.r2.(fw) = flux.stf.(fw)./flux.ra; % calculate nondimensional number R2 disregarding MSE budget closure
    end

    if strcmp(type, 'era5') | strcmp(type, 'erai')
        for fn = {'ra', 'sshf', 'slhf', 'cp', 'lsp', 'e'}; fname = fn{1};
            % take zonal averages
            flux_z.(fname) = squeeze(nanmean(flux.(fname), 1));
            % take time averages
            for l = {'lo', 'l', 'o'}; land = l{1};
                if strcmp(land, 'lo'); flux_n.(land).(fname) = flux.(fname);
                elseif strcmp(land, 'l'); flux_n.(land).(fname) = flux.(fname) .*mask.ocean;
                elseif strcmp(land, 'o'); flux_n.(land).(fname) = flux.(fname) .*mask.land;
                end

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
        for fn = {'ra', 'hfls', 'hfss', 'prc', 'pr', 'evspsbl'}; fname = fn{1};
            % take zonal averages
            flux_z.(fname) = squeeze(nanmean(flux.(fname), 1));

            % take time averages
            for l = {'lo', 'l', 'o'}; land = l{1};
                if strcmp(land, 'lo'); flux_n.(land).(fname) = flux.(fname);
                elseif strcmp(land, 'l'); flux_n.(land).(fname) = flux.(fname) .*mask.ocean;
                elseif strcmp(land, 'o'); flux_n.(land).(fname) = flux.(fname) .*mask.land;
                end

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
        foldername = sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s/', type, par.model, par.lat_interp);
    end
    for fn = {'stf', 'res', 'r1', 'r2'}; fname = fn{1};
        for f = {'mse', 'dse'}; fw = f{1};
            % take zonal average
            flux_z.(fname).(fw) = squeeze(nanmean(flux.(fname).(fw), 1));
            % take time averages
            for l = {'lo', 'l', 'o'}; land = l{1};
                if strcmp(land, 'lo'); flux_n.(land).(fname).(fw) = flux.(fname).(fw);
                elseif strcmp(land, 'l'); flux_n.(land).(fname).(fw) = flux.(fname).(fw) .*mask.ocean;
                elseif strcmp(land, 'o'); flux_n.(land).(fname).(fw) = flux.(fname).(fw) .*mask.land;
                end

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
    for v = {'flux', 'flux_t', 'flux_z', 'flux_zt', 'vh'}; varname = v{1};
        save(sprintf('%s%s', foldername, varname), varname, 'lat');
    end

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
    proc_temp(type, par) % calculate RCE and RAE temperature profiles
    proc_ma(type, par) % calculate moist adiabats corresponding to RCE profiles
    % proc_temp_mon_lat(type, par) % calculate mon x lat temperature profiles
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
                    if any(strcmp(fname, {'stf', 'res', 'r1', 'r2'}))
                        for f = {'mse', 'dse'}; fw = f{1};
                            if strcmp(land, 'lo'); flux_n.(land).(fname).(fw) = flux.(fname).(fw);
                            elseif strcmp(land, 'l'); flux_n.(land).(fname).(fw) = flux.(fname).(fw) .*mask.ocean;
                            elseif strcmp(land, 'o'); flux_n.(land).(fname).(fw) = flux.(fname).(fw) .*mask.land;
                            end
                        end
                    else
                        if strcmp(land, 'lo'); flux_n.(land).(fname) = flux.(fname);
                        elseif strcmp(land, 'l'); flux_n.(land).(fname) = flux.(fname) .*mask.ocean;
                        elseif strcmp(land, 'o'); flux_n.(land).(fname) = flux.(fname) .*mask.land;
                        end
                    end
                end

                % lon x lat structure over various time averages
                for t = {'ann', 'djf', 'jja', 'mam', 'son'}; time = t{1};
                    for fn = fieldnames(flux)'; fname = fn{1};
                        if any(strcmp(fname, {'stf', 'res', 'r1', 'r2'}))
                            for f = {'mse', 'dse'}; fw = f{1};
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
                            end
                        else
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
    for f = {'mse', 'dse'}; fw = f{1};
        for c = fieldnames(rcae_t.lo.ann.(fw))'; crit = c{1};
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

                        filt.(land).(time).(fw).(crit).(regime) = nan(size(rcae_t.(land).(time).(fw).(crit))); % create empty arrays to store filtering array
                        if strcmp(regime, 'rce'); filt.(land).(time).(fw).(crit).(regime)(rcae_t.(land).(time).(fw).(crit)==1)=1; % set RCE=1, elsewhere nan
                        elseif strcmp(regime, 'rae'); filt.(land).(time).(fw).(crit).(regime)(rcae_t.(land).(time).(fw).(crit)==-1)=1; % set RAE=1, elsewhere nan
                        end

                        temp_t.(land).(time).(fw).(crit).(regime) = temp_n.(land).(time) .* filt.(land).(time).(fw).(crit).(regime);

                        nanfilt.(regime) = nan(size(temp_t.(land).(time).(fw).(crit).(regime)));
                        nanfilt.(regime)(~isnan(temp_t.(land).(time).(fw).(crit).(regime))) = 1;

                        % take cosine-weighted meridional average
                        for d = {'all', 'nh', 'sh', 'tp'}; domain = d{1};
                            if strcmp(regime, 'rce')
                                if strcmp(domain, 'all')
                                    cosw = repmat(cosd(lat)', [size(nanfilt.(regime), 1) 1]);
                                    denm = temp_t.(land).(time).(fw).(crit).(regime);
                                    nume = nanfilt.(regime);
                                elseif strcmp(domain, 'nh')
                                    cosw = repmat(cosd(lat(lat>30))', [size(nanfilt.(regime), 1) 1]);
                                    denm = temp_t.(land).(time).(fw).(crit).(regime)(:,lat>30,:);
                                    nume = nanfilt.(regime)(:,lat>30,:);
                                elseif strcmp(domain, 'sh')
                                    cosw = repmat(cosd(lat(lat<-30))', [size(nanfilt.(regime), 1) 1]);
                                    denm = temp_t.(land).(time).(fw).(crit).(regime)(:,lat<-30,:);
                                    nume = nanfilt.(regime)(:,lat<-30,:);
                                elseif strcmp(domain, 'tp')
                                    cosw = repmat(cosd(lat(abs(lat)<30))', [size(nanfilt.(regime), 1) 1]);
                                    denm = temp_t.(land).(time).(fw).(crit).(regime)(:,abs(lat)<30,:);
                                    nume = nanfilt.(regime)(:,abs(lat)<30,:);
                                end
                            elseif strcmp(regime, 'rae')
                                if strcmp(domain, 'all')
                                    cosw = repmat(cosd(lat)', [size(nanfilt.(regime), 1) 1]);
                                    denm = temp_t.(land).(time).(fw).(crit).(regime);
                                    nume = nanfilt.(regime);
                                elseif strcmp(domain, 'nh')
                                    cosw = repmat(cosd(lat(lat>0))', [size(nanfilt.(regime), 1) 1]);
                                    denm = temp_t.(land).(time).(fw).(crit).(regime)(:,lat>0,:);
                                    nume = nanfilt.(regime)(:,lat>0,:);
                                elseif strcmp(domain, 'sh')
                                    cosw = repmat(cosd(lat(lat<-0))', [size(nanfilt.(regime), 1) 1]);
                                    denm = temp_t.(land).(time).(fw).(crit).(regime)(:,lat<-0,:);
                                    nume = nanfilt.(regime)(:,lat<-0,:);
                                end
                            end

                            ta.(regime).(domain).(fw).(crit).(land).(time) = nansum(cosw.*denm, 2) ./ nansum(cosw.*nume, 2); % area-weighted meridional average
                            ta.(regime).(domain).(fw).(crit).(land).(time) = squeeze(nanmean(ta.(regime).(domain).(fw).(crit).(land).(time), 1));
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
                                        cosw = repmat(cosd(lat(abs(lat)<30))', [size(nanfilt.(regime), 1) 1]);
                                        denm = srfc_tf.(land).(time).(fw).(crit).(regime).(vname)(:, abs(lat)<30);
                                        nume = nanfilt.(regime)(:, abs(lat)<30);
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
                            % ma = calc_ma_dew(ma, grid.dim3.plev, par); % compute moist adiabat with dew point temperature
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

function comp_flux(par)
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
end
function disp_global_flux(par)
    root = '/project2/tas1/miyawaki/projects/002/data';
    load(sprintf('%s/proc/comp/comp_zt', root));
    load(sprintf('%s/proc/comp/comp_z', root));
    load(sprintf('%s/proc/comp/comp_t', root));

    % output global averages
    for fn = {'tsr', 'ssr', 'ttr', 'str'}; fname = fn{1};
        disp( sprintf('CERES %s is %g Wm^-2', upper(fname), nansum(cosd(lat).*ceres_zt.rad.(fname))/nansum(cosd(lat))) );
        disp( sprintf('ERA-I %s is %g Wm^-2', upper(fname), nansum(cosd(lat).*erai_zt.rad.(fname))/nansum(cosd(lat))) );
        disp( sprintf('ERA5 %s is %g Wm^-2', upper(fname), nansum(cosd(lat).*era5_zt.rad.(fname))/nansum(cosd(lat))) );
    end
end

% helper functions
function rcae = def_rcae(type, flux, par)
    for f = {'mse', 'dse'}; fw = f{1};
        % identify locations of RCE using threshold epsilon (ep)
        rcae.(fw).def = zeros(size(flux.r1.(fw)));
        rcae.(fw).def(abs(flux.r1.(fw)) < par.ep) = 1;
        % identify locations of RAE using threshold gamma (ga)
        rcae.(fw).def(flux.r1.(fw) > par.ga) = -1;

        % add additional criteria for RCE that P-E>0
        rcae.(fw).pe = zeros(size(flux.r1.(fw)));
        if strcmp(type, 'era5') | strcmp(type, 'erai')
            rcae.(fw).pe(abs(flux.r1.(fw)) < par.ep & (flux.cp+flux.lsp+flux.e>0)) = 1; % note that evap is defined negative into atmosphere
        elseif strcmp(type, 'gcm')
            rcae.(fw).pe(abs(flux.r1.(fw)) < par.ep & (flux.pr-flux.evspsbl>0)) = 1;
        end
        rcae.(fw).pe(flux.r1.(fw) > par.ga) = -1;

        % add additional criteria for RCE that (P_ls - E)<<1
        rcae.(fw).cp = zeros(size(flux.r1.(fw)));
        if strcmp(type, 'era5') | strcmp(type, 'erai')
            rcae.(fw).cp(abs(flux.r1.(fw)) < par.ep & abs(flux.lsp./flux.cp<par.ep_cp)) = 1; % note that evap is defined negative into atmosphere
        elseif strcmp(type, 'gcm')
            rcae.(fw).cp(abs(flux.r1.(fw)) < par.ep & abs((flux.pr-flux.prc)./flux.prc<par.ep_cp)) = 1; % note that evap is defined negative into atmosphere
        end
        rcae.(fw).cp(flux.r1.(fw) > par.ga) = -1;
    end
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

