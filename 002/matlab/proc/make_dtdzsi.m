function make_dtdzsi(type, par)
% compute model lapse rate in lat x plev x mon (add surface data after converting to sigma)
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

    if strcmp(type, 'gcm') & any(contains(par.model, {'GISS-E2-H', 'GISS-E2-R'})) % zg data in GISS-E2-H has an anomalous lat grid
        zg_lat = double(ncread(fullpath, 'lat'));
        zg = permute(zg, [2 1 3 4]);
        zg = interp1(zg_lat, zg, grid.dim3.lat);
        zg = permute(zg, [2 1 3 4]);
    end

    coarse_si = 1e-5*grid.dim3.plev;
    [~,si0] = max(coarse_si); % surface sigma index

    % create surface mask
    if strcmp(type, 'era5') | strcmp(type, 'erai') | strcmp(type, 'era5c')
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

                if all(isnan(ta_tmp)) | all(isnan(zg_tmp))
                    ta_si(:,lo,la,mo) = nan(size(grid.dim3.si));
                    zg_si(:,lo,la,mo) = nan(size(grid.dim3.si));
                else

                    % add surface data
                    if strcmp(type, 'era5') | strcmp(type, 'erai') | strcmp(type, 'era5c')
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
    end
    clear ta_sm zg_sm

    % % interpolate to higher resolution grid
    % ta_si = interp1(grid.dim3.si, ta_si, )

    % calculate lapse rate before taking zonal average
    dtdzsi = -1e3*(ta_si(2:end,:,:,:)-ta_si(1:end-1,:,:,:))./(zg_si(2:end,:,:,:)-zg_si(1:end-1,:,:,:)); % lapse rate in K/km
    dtdzsi = interp1(1/2*(grid.dim3.si(2:end)+grid.dim3.si(1:end-1)), dtdzsi, grid.dim3.si);
    dtdzsi(1,:,:,:) = -1e3*(ta_si(2,:,:,:)-ta_si(1,:,:,:))./(zg_si(2,:,:,:)-zg_si(1,:,:,:));
    dtdzsi = permute(dtdzsi, [2 3 1 4]);

    if any(strcmp(type, {'era5', 'era5c', 'erai'})); newdir=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s', type, par.(type).yr_span);
    elseif strcmp(type, 'merra2'); newdir=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s', type, par.(type).yr_span);
    elseif strcmp(type, 'gcm'); newdir=sprintf('/project2/tas1/miyawaki/projects/002/data/read/gcm/%s/%s', par.model, par.gcm.clim);
    elseif strcmp(type, 'echam'); newdir=sprintf('/project2/tas1/miyawaki/projects/002/data/read/echam/%s', par.echam.clim);
    elseif any(strcmp(type, {'echam_ml', 'echam_pl'})); newdir=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s', type, par.(type).yr_span); end;
    if ~exist(newdir, 'dir'); mkdir(newdir); end
    filename='dtdzsi.mat';
    save(sprintf('%s/%s', newdir, filename), 'dtdzsi', '-v7.3');

end
