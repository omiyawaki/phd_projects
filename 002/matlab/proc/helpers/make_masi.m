function make_masi(type, par)
% compute moist adiabats at every lon, lat, mon

    prefix = make_prefix(type, par);
    foldername = make_savedir(type, par);
    temp = load_temp(type, par);
    zg = load_zg(type, par);

    load(sprintf('%s/grid.mat', prefix)); % read grid data
    load(sprintf('%s/srfc.mat', prefix)); % load surface data

    if strcmp(type, 'gcm') & any(contains(par.model, {'GISS-E2-H', 'GISS-E2-R'})) % zg data in GISS-E2-H has an anomalous lat grid
        var = 'zg';
        file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/gcm/%s/%s_Amon_%s_%s_r1i1p1_*.ymonmean.nc', par.model, var, par.model, par.gcm.clim));
        fullpath=sprintf('%s/%s', file.folder, file.name);
        zg_lat = double(ncread(fullpath, 'lat'));
        zg = permute(zg, [2 1 3 4]);
        zg = interp1(zg_lat, zg, grid.dim3.lat);
        zg = permute(zg, [2 1 3 4]);
    end

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
                        if any(strcmp(type, {'era5', 'erai', 'era5c'}))
                            ma_in.sp = srfc.sp(ilon,ilat,imon);
                            ma_in.pinit = par.ma_init*ma_in.sp;
                            ma_in.t2m = interp1(grid.dim3.plev, squeeze(temp(ilon,ilat,:,imon)), ma_in.pinit);
                            ma_in.d2m = ma_in.t2m; % saturated
                        elseif strcmp(type, 'merra2')
                            ma_in.PS = srfc.PS(ilon,ilat,imon);
                            ma_in.pinit = par.ma_init*ma_in.PS;
                            ma_in.T2M = interp1(grid.dim3.plev, squeeze(temp(ilon,ilat,:,imon)), ma_in.pinit);
                            ma_in.D2M = ma_in.T2M; % note that merra2 actually uses specific humidity; this is a quick workaround for now
                        elseif strcmp(type, 'echam')
                            ma_in.aps = srfc.aps(ilon,ilat,imon);
                            ma_in.pinit = par.ma_init*ma_in.aps;
                            ma_in.temp2 = interp1(grid.dim3.plev, squeeze(temp(ilon,ilat,:,imon)), ma_in.pinit);
                            ma_in.dew2 = ma_in.temp2; % saturated
                        elseif any(strcmp(type, {'gcm', 'jra55'}))
                            ma_in.ps = srfc.ps(ilon,ilat,imon);
                            ma_in.pinit = par.ma_init*ma_in.ps;
                            ma_in.tas = interp1(grid.dim3.plev, squeeze(temp(ilon,ilat,:,imon)), ma_in.pinit);
                            ma_in.hurs = 100; % saturated
                        end
                        ma_in.zs = interp1(grid.dim3.plev, squeeze(zg(ilon,ilat,:,imon)), ma_in.pinit);
                    end
                end

                if any(strcmp(type, {'era5', 'era5c', 'erai', 'echam', 'merra2'}));
                    if strcmp(type, 'merra2') & strcmp(par.ma_init, 'surf')
                        ma_si(ilon,ilat,:,imon) = calc_ma_q_si(ma_in, grid.dim3.plev, par, type, grid);
                    else
                        ma_si(ilon,ilat,:,imon) = calc_ma_dew_si(ma_in, grid.dim3.plev, par, type, grid);
                    end
                elseif any(strcmp(type, {'gcm', 'jra55'}))
                    ma_si(ilon,ilat,:,imon) = calc_ma_hurs_si(ma_in, grid.dim3.plev, par, grid);
                end

            end
        end
    end

    if strcmp(par.ma_init, 'surf')
        filename=sprintf('ma_si_%s.mat', par.ma_init);
    else
        filename=sprintf('ma_si_%g.mat', par.ma_init);
    end
    save(sprintf('%s/%s', foldername, filename), 'ma_si', '-v7.3');
end
