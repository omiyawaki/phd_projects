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
