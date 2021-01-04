function proc_net_flux(type, par)
% calculates the global TOA energy imbalance using ERA-Interim data
    if strcmp(type, 'era5') | strcmp(type, 'erai') | strcmp(type, 'era5c')
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s', type, par.(type).yr_span);
    elseif strcmp(type, 'gcm')
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s/%s', type, par.model, par.gcm.clim);
    elseif strcmp(type, 'echam')
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s', type, par.echam.clim);
    end
    load(sprintf('%s/grid.mat', prefix)); % read grid data
    load(sprintf('%s/rad.mat', prefix)); % read radiation data
    load(sprintf('%s/stf.mat', prefix)); % read surface turbulent flux data

    if strcmp(type, 'era5') | strcmp(type, 'erai') | strcmp(type, 'era5c')
        net_toa_raw = rad.tsr + rad.ttr; % compute net radiative fluxes at TOA, positive down
    elseif strcmp(type, 'gcm')
        net_toa_raw = - rad.rsut + rad.rsdt - rad.rlut;
    elseif strcmp(type, 'echam')
        net_toa_raw = rad.trad0 + rad.srad0;
    end

    net_toa_tz = squeeze(nanmean(nanmean( net_toa_raw, 1 ), 3))'; % take zonal and time avera5ges and transpose to have same dimensions as latitude grid
    net_toa = nansum(cosd(grid.dim2.lat).*net_toa_tz) / nansum(cosd(grid.dim2.lat));
    if strcmp(type, 'era5') | strcmp(type, 'erai') | strcmp(type, 'era5c')
        disp( sprintf('The net radiative imbalance in %s at TOA is %g Wm^-2.', type, net_toa) );
    elseif strcmp(type, 'gcm')
        disp( sprintf('The net radiative imbalance in %s at TOA is %g Wm^-2.', par.model, net_toa) );
    elseif strcmp(type, 'echam')
        disp( sprintf('The net radiative imbalance in ECHAM at TOA is %g Wm^-2.', net_toa) );
    end

    if strcmp(type, 'era5') | strcmp(type, 'erai') | strcmp(type, 'era5c')
        net_sfc_raw = rad.ssr + rad.str + stf.sshf + stf.slhf; % compute net radiative fluxes at surface, positive down
    elseif strcmp(type, 'gcm')
        net_sfc_raw = - rad.rsus + rad.rsds - rad.rlus + rad.rlds - stf.hfss - stf.hfls; % compute net radiative fluxes at surface, positive down
    elseif strcmp(type, 'echam')
        net_sfc_raw = rad.srads + rad.trads + stf.ahfl + stf.ahfs; % compute net radiative fluxes at surface, positive down
    end
    net_sfc_tz = squeeze(nanmean(nanmean( net_sfc_raw, 1 ), 3))'; % take zonal and time avera5ges and transpose to have same dimensions as latitude grid
    net_sfc = nansum(cosd(grid.dim2.lat).*net_sfc_tz) / nansum(cosd(grid.dim2.lat));
    if strcmp(type, 'era5') | strcmp(type, 'erai') | strcmp(type, 'era5c')
        disp( sprintf('The net radiative imbalance in %s at the surface is %g Wm^-2.', type, net_sfc) );
    elseif strcmp(type, 'gcm')
        disp( sprintf('The net radiative imbalance in %s at the surface is %g Wm^-2.', par.model, net_sfc) );
    elseif strcmp(type, 'echam')
        disp( sprintf('The net radiative imbalance in ECHAM at the surface is %g Wm^-2.', net_sfc) );
    end
end
