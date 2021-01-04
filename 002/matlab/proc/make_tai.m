function make_tai(type, par)
    % add surface data to temperature and interpolate to hi-res grid
    if strcmp(type, 'era5') | strcmp(type, 'erai') | strcmp(type, 'era5c')
        foldername = sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s', type, par.(type).yr_span, par.lat_interp);
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s', type, par.(type).yr_span);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s', type, par.(type).yr_span);
        file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/%s/temp/%s_temp_%s.ymonmean.nc', type, type, par.(type).yr_span));
        fullpath=sprintf('%s/%s', file.folder, file.name);
        temp = double(ncread(fullpath, 't'));
        file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/%s/zg/%s_zg_%s.ymonmean.nc', type, type, par.(type).yr_span));
        fullpath=sprintf('%s/%s', file.folder, file.name);
        zg = double(ncread(fullpath, 'z')); zg = zg/par.g; % convert from geopotential to height in m
    elseif strcmp(type, 'merra2')
        foldername = sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s', type, par.(type).yr_span, par.lat_interp);
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s', type, par.(type).yr_span);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s', type, par.(type).yr_span);
        file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/%s/temp/%s_temp_%s.ymonmean.nc', type, type, par.(type).yr_span));
        fullpath=sprintf('%s/%s', file.folder, file.name);
        temp = double(ncread(fullpath, 'T'));
        file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/%s/zg/%s_zg_%s.ymonmean.nc', type, type, par.(type).yr_span));
        fullpath=sprintf('%s/%s', file.folder, file.name);
        zg = double(ncread(fullpath, 'H'));
    elseif strcmp(type, 'gcm')
        foldername = sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s/%s/', type, par.model, par.gcm.clim, par.lat_interp);
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
        foldername = sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s/', type, par.echam.clim, par.lat_interp);
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
    elseif any(strcmp(type, {'echam_ml', 'echam_pl'}))
        foldername = sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s', type, par.(type).yr_span, par.lat_interp);
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s', type, par.(type).yr_span);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s', type, par.(type).yr_span);
        var = 't';
        file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/%s/ATM*.ymonmean.nc', type));
        fullpath=sprintf('%s/%s', file.folder, file.name);
        temp = double(ncread(fullpath, var));
        var = 'geopoth';
        file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/%s/ATM_*.ymonmean.nc', type));
        fullpath=sprintf('%s/%s', file.folder, file.name);
        zg = double(ncread(fullpath, var));
    end
    load(sprintf('%s/grid.mat', prefix)); % read grid data
    load(sprintf('%s/srfc.mat', prefix)); % read surface variable data

    if strcmp(type, 'gcm') & any(contains(par.model, {'GISS-E2-H', 'GISS-E2-R'})) % zg data in GISS-E2-H has an anomalous lat grid
        zg_lat = double(ncread(fullpath, 'lat'));
        zg = permute(zg, [2 1 3 4]);
        zg = interp1(zg_lat, zg, grid.dim3.lat);
        zg = permute(zg, [2 1 3 4]);
    end

    % create surface mask
    if strcmp(type, 'era5') | strcmp(type, 'erai') | strcmp(type, 'era5c')
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
        if contains(type, 'era5')
            pb2 = CmdLineProgressBar("Lat...");
        end
        for ilat=1:size(pa_plus,3)
            if contains(type, 'era5')
                pb2.print(ilat, size(pa_plus,3));
            end
            for time=1:size(pa_plus,4)
                ta_plus(:,ilon,ilat,time) = ta_plus(sort_index(:,ilon,ilat,time),ilon,ilat,time); % sort temperature (has to be in loop because sort_index works for vector calls only)
                zg_plus(:,ilon,ilat,time) = zg_plus(sort_index(:,ilon,ilat,time),ilon,ilat,time); % sort temperature (has to be in loop because sort_index works for vector calls only)

                if all(isnan(ta_plus)) | sum(isnan(ta_plus)) == 1
                    tai_sm.lo(:,ilon,ilat,time) = nan(size(par.pa));
                end

                if all(isnan(zg_plus)) | sum(isnan(zg_plus)) == 1
                    zg_sm.lo(:,ilon,ilat,time) = nan(size(par.pa));
                end

                if sum(isnan(ta_plus))>1 & sum(isnan(zg_plus))>1

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
    end
    clear pa_plus ta_plus zg_plus; % clear unneeded variables
    tai_sm.lo = permute(tai_sm.lo, [2 3 1 4]); % bring plev back to third dimension
    tai = tai_sm.lo;
    zgi_sm.lo = permute(zgi_sm.lo, [2 3 1 4]); % bring plev back to third dimension
    zgi = zgi_sm.lo;

    if any(strcmp(type, {'era5', 'era5c', 'erai'})); newdir=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s', type, par.(type).yr_span);
    elseif strcmp(type, 'merra2'); newdir=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s', type, par.(type).yr_span);
    elseif strcmp(type, 'gcm'); newdir=sprintf('/project2/tas1/miyawaki/projects/002/data/read/gcm/%s/%s', par.model, par.gcm.clim);
    elseif strcmp(type, 'echam'); newdir=sprintf('/project2/tas1/miyawaki/projects/002/data/read/echam/%s', par.echam.clim);
    elseif any(strcmp(type, {'echam_ml', 'echam_pl'})); newdir=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s', type, par.(type).yr_span); end;
    if ~exist(newdir, 'dir'); mkdir(newdir); end
    filename='tai.mat';
    save(sprintf('%s/%s', newdir, filename), 'tai', 'zgi', '-v7.3');

end
