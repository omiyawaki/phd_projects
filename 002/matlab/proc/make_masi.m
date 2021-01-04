function make_masi(type, par)
% compute moist adiabats at every lon, lat, mon
    if any(strcmp(type, {'era5', 'erai', 'era5c'}))
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s', type, par.(type).yr_span);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s', type, par.(type).yr_span);
        file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/%s/temp/%s_temp_%s.ymonmean.nc', type, type, par.(type).yr_span));
        fullpath=sprintf('%s/%s', file.folder, file.name);
        temp = double(ncread(fullpath, 't'));
        file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/%s/zg/%s_zg_%s.ymonmean.nc', type, type, par.(type).yr_span));
        fullpath=sprintf('%s/%s', file.folder, file.name);
        zg = double(ncread(fullpath, 'z')); zg = zg/par.g; % convert from geopotential to height in m
    elseif strcmp(type, 'merra2')
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s', type, par.(type).yr_span);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s', type, par.(type).yr_span);
        file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/%s/temp/%s_temp_%s.ymonmean.nc', type, type, par.(type).yr_span));
        fullpath=sprintf('%s/%s', file.folder, file.name);
        temp = double(ncread(fullpath, 'T'));
        file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/%s/zg/%s_zg_%s.ymonmean.nc', type, type, par.(type).yr_span));
        fullpath=sprintf('%s/%s', file.folder, file.name);
        zg = double(ncread(fullpath, 'H')); % convert from geopotential to height in m
    elseif strcmp(type, 'gcm')
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s/%s', type, par.model, par.gcm.clim);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s', type, par.model, par.gcm.clim);
        var = 'ta';
        file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/gcm/%s/%s_Amon_%s_%s_r1i1p1_*.ymonmean.nc', par.model, var, par.model, par.gcm.clim));
        fullpath=sprintf('%s/%s', file.folder, file.name);
        temp = double(ncread(fullpath, var));
        var = 'zg';
        file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/gcm/%s/%s_Amon_%s_%s_r1i1p1_*.ymonmean.nc', par.model, var, par.model, par.gcm.clim));
        fullpath=sprintf('%s/%s', file.folder, file.name);
        zg = double(ncread(fullpath, var));
    elseif strcmp(type, 'echam')
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s', type, par.echam.clim);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s', type, par.echam.clim);
        var = 't';
        if contains(par.echam.clim, 'rp000')
            file=dir(sprintf('/project2/tas1/ockham/data11/tas/echam-aiv_rcc_6.1.00p1/%s/ATM_%s_0020_39.nc', par.echam.clim, par.echam.clim));
        else
            file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/echam/ATM*_%s_*.ymonmean.nc', par.echam.clim));
        end
        fullpath=sprintf('%s/%s', file.folder, file.name);
        temp = double(ncread(fullpath, var));
        var = 'geopoth';
        fullpath=sprintf('%s/%s', file.folder, file.name);
        zg = double(ncread(fullpath, var));
    elseif contains(type, 'echam_pl')
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s', type, par.(type).yr_span);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s', type, par.(type).yr_span);
    end

    load(sprintf('%s/grid.mat', prefix)); % read grid data
    load(sprintf('%s/srfc.mat', prefix)); % load surface data

    if strcmp(type, 'gcm') & any(contains(par.model, {'GISS-E2-H', 'GISS-E2-R'})) % zg data in GISS-E2-H has an anomalous lat grid
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
                        elseif strcmp(type, 'gcm')
                            ma_in.ps = srfc.ps(ilon,ilat,imon);
                            ma_in.pinit = par.ma_init*ma_in.ps;
                            ma_in.tas = interp1(grid.dim3.plev, squeeze(temp(ilon,ilat,:,imon)), ma_in.pinit);
                            ma_in.hurs = 100; % saturated
                        end
                        ma_in.zs = interp1(grid.dim3.plev, squeeze(zg(ilon,ilat,:,imon)), ma_in.pinit);
                    end
                end

                if any(strcmp(type, {'era5', 'era5c', 'erai', 'echam', 'merra2'}));
                    ma_si(ilon,ilat,:,imon) = calc_ma_dew_si(ma_in, grid.dim3.plev, par, type, grid);
                elseif strcmp(type, 'gcm')
                    ma_si(ilon,ilat,:,imon) = calc_ma_hurs_si(ma_in, grid.dim3.plev, par, grid);
                end

            end
        end
    end

    if any(strcmp(type, {'era5', 'era5c', 'erai'})); newdir=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s', type, par.(type).yr_span);
    elseif strcmp(type, 'merra2'); newdir=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s', type, par.(type).yr_span);
    elseif strcmp(type, 'gcm'); newdir=sprintf('/project2/tas1/miyawaki/projects/002/data/read/gcm/%s/%s', par.model, par.gcm.clim);
    elseif strcmp(type, 'echam'); newdir=sprintf('/project2/tas1/miyawaki/projects/002/data/read/echam/%s', par.echam.clim);
    elseif contains(type, 'echam_pl'); newdir=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s', type, par.(type).yr_span); end;
    if ~exist(newdir, 'dir'); mkdir(newdir); end
    if strcmp(par.ma_init, 'surf')
        filename=sprintf('ma_si_%s.mat', par.ma_init);
    else
        filename=sprintf('ma_si_%g.mat', par.ma_init);
    end
    save(sprintf('%s/%s', newdir, filename), 'ma_si', '-v7.3');
end
