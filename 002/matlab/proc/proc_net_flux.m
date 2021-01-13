function proc_net_flux(type, par)
% calculates the global TOA energy imbalance using ERA-Interim data

    prefix = make_prefix(type, par);
    prefix_proc = make_prefix_proc(type, par);

    load(sprintf('%s/grid.mat', prefix)); % read grid data
    load(sprintf('%s/rad.mat', prefix)); % read radiation data
    load(sprintf('%s/stf.mat', prefix)); % read surface turbulent flux data

    lat = grid.dim2.lat;

    if strcmp(type, 'era5') | strcmp(type, 'erai') | strcmp(type, 'era5c')
        net_toa_raw = rad.tsr + rad.ttr; % compute net radiative fluxes at TOA, positive down
    elseif strcmp(type, 'merra2');
        net_toa_raw = rad.SWTNT - rad.LWTUP; % net flux at TOA
    elseif any(strcmp(type, {'gcm', 'jra55'}))
        net_toa_raw = - rad.rsut + rad.rsdt - rad.rlut;
    elseif strcmp(type, 'echam')
        net_toa_raw = rad.trad0 + rad.srad0;
    end

    net_toa_tz = squeeze(nanmean(nanmean( net_toa_raw, 1 ), 3))'; % take zonal and time avera5ges and transpose to have same dimensions as latitude grid
    if grid.dim2.lat(1) > 0
        vh_mon = cumtrapz(flip(deg2rad(lat),1), flip(2*pi*par.a^2*cosd(lat).*net_toa_tz,1), 1); % cumulatively integrate
    else
        vh_mon = cumtrapz(deg2rad(lat), 2*pi*par.a^2*cosd(lat).*net_toa_tz, 1); % cumulatively integrate
    end
    net_toa = vh_mon(end)/(4*pi*par.a^2);

        % net_toa = nansum(cosd(grid.dim2.lat).*net_toa_tz) / nansum(cosd(grid.dim2.lat));
    if any(strcmp(type, {'era5', 'era5c', 'erai', 'merra2', 'jra55'}))
        disp( sprintf('The net radiative imbalance in %s at TOA is %g Wm^-2.', upper(type), net_toa) );
    elseif strcmp(type, 'gcm')
        disp( sprintf('The net radiative imbalance in %s at TOA is %g Wm^-2.', par.model, net_toa) );
    elseif strcmp(type, 'echam')
        disp( sprintf('The net radiative imbalance in ECHAM at TOA is %g Wm^-2.', net_toa) );
    end

    if strcmp(type, 'era5') | strcmp(type, 'erai') | strcmp(type, 'era5c')
        net_sfc_raw = rad.ssr + rad.str + stf.sshf + stf.slhf; % compute net radiative fluxes at surface, positive down
    elseif strcmp(type, 'merra2')
        net_sfc_raw = -rad.LWGNT - rad.SWGNT + stf.HFLUX + stf.EFLUX; % compute net radiative fluxes at surface, positive down
    elseif any(strcmp(type, {'gcm', 'jra55'}))
        net_sfc_raw = - rad.rsus + rad.rsds - rad.rlus + rad.rlds - stf.hfss - stf.hfls; % compute net radiative fluxes at surface, positive down
    elseif strcmp(type, 'echam')
        net_sfc_raw = rad.srads + rad.trads + stf.ahfl + stf.ahfs; % compute net radiative fluxes at surface, positive down
    end
    net_sfc_tz = squeeze(nanmean(nanmean( net_sfc_raw, 1 ), 3))'; % take zonal and time avera5ges and transpose to have same dimensions as latitude grid
    net_sfc = nansum(cosd(grid.dim2.lat).*net_sfc_tz) / nansum(cosd(grid.dim2.lat));
    if any(strcmp(type, {'era5', 'era5c', 'erai', 'merra2', 'jra55'}))
        disp( sprintf('The net radiative imbalance in %s at the surface is %g Wm^-2.', upper(type), net_sfc) );
    elseif strcmp(type, 'gcm')
        disp( sprintf('The net radiative imbalance in %s at the surface is %g Wm^-2.', par.model, net_sfc) );
    elseif strcmp(type, 'echam')
        disp( sprintf('The net radiative imbalance in ECHAM at the surface is %g Wm^-2.', net_sfc) );
    end
end
