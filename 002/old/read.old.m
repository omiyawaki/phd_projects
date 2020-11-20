% old functions
function read_albedo(type, par) % read surface albedo data from Tiffany's ECHAM6 file
    if any(strcmp(type, {'era5', 'erai'})) % only for MPI-ESM-LR
        % dimensions are (lon x lat x time)
        albedo = double(ncread(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/%s/albedo/%s_albedo_%s.ymonmean.nc', type, type, par.(type).yr_span), 'fal'));
        save(sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/albedo.mat', type), 'albedo');
    elseif strcmp(type, 'gcm') & strcmp(par.model, 'MPI-ESM-LR') % only for MPI-ESM-LR
        for i=1:length(par.gcm.vars.rad); var = par.gcm.vars.rad{i};
            prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s/%s', type, par.model, par.gcm.clim);
            file=dir(sprintf('/project2/tas1/echam6/piControl/BOT_piControl_r1i1p1-P_1900_1949.nc'));
            fullpath=sprintf('%s/%s', file.folder, file.name);
            albedo=ncread(fullpath, 'albedo');
            lat=ncread(fullpath, 'lat');
            lon=ncread(fullpath, 'lon');
        end
        load(sprintf('%s/grid.mat', prefix)); % read grid data
        albedo = interp1(lon, albedo, grid.dim2.lon); % interpolate to MPI-ESM-LR longitude grid
        albedo = permute(albedo, [2 1 3]); % bring lat to 1st
        albedo = interp1(lat, albedo, grid.dim2.lat); % interpolate to MPI-ESM-LR latitude grid
        albedo = permute(albedo, [2 1 3]); % bring lat back to 2nd
        newdir=sprintf('/project2/tas1/miyawaki/projects/002/data/read/gcm/%s', par.model);
        if ~exist(newdir, 'dir'); mkdir(newdir); end
        filename='albedo.mat';
        save(sprintf('%s/%s', newdir, filename), 'albedo');
    end
end
function read_snow(type, par) % read snow depth data from Tiffany's ECHAM6 file
    if any(strcmp(type, {'era5', 'erai'})) % only for MPI-ESM-LR
        % dimensions are (lon x lat x time)
        % sn = double(ncread(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/%s/sn/%s_sn_%s.ymonmean.nc', type, type, par.(type).yr_span), 'fal'));
        % save(sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/sn.mat', type), 'sn');
    elseif strcmp(type, 'gcm') & strcmp(par.model, 'MPI-ESM-LR') % only for MPI-ESM-LR
        for i=1:length(par.gcm.vars.rad); var = par.gcm.vars.rad{i};
            prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s', type, par.model);
            file=dir(sprintf('/project2/tas1/echam6/piControl/BOT_piControl_r1i1p1-P_1900_1949.nc'));
            fullpath=sprintf('%s/%s', file.folder, file.name);
            sn=ncread(fullpath, 'siced');
            lat=ncread(fullpath, 'lat');
            lon=ncread(fullpath, 'lon');
        end
        load(sprintf('%s/grid.mat', prefix)); % read grid data
        sn = interp1(lon, sn, grid.dim2.lon); % interpolate to MPI-ESM-LR longitude grid
        sn = permute(sn, [2 1 3]); % bring lat to 1st
        sn = interp1(lat, sn, grid.dim2.lat); % interpolate to MPI-ESM-LR latitude grid
        sn = permute(sn, [2 1 3]); % bring lat back to 2nd
        newdir=sprintf('/project2/tas1/miyawaki/projects/002/data/read/gcm/%s', par.model);
        if ~exist(newdir, 'dir'); mkdir(newdir); end
        filename='sn.mat';
        save(sprintf('%s/%s', newdir, filename), 'sn');
    end
end

function make_mlev(type, par)
    if any(strcmp(type, {'echam_ml', 'echam_pl'}))
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s', type);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s', type);
        var = 'aps';
        file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/echam_ml/ATM_*.ymonmean.nc'));
        fullpath=sprintf('%s/%s', file.folder, file.name);
        ps_orig = double(ncread(fullpath, var));
    else
        error('This code only works for data output in the model vertical grid.')
    end

    load(sprintf('%s/grid.mat', prefix)); % read grid data
    load(sprintf('%s/srfc.mat', prefix)); % load surface data

    % compute sigma from a and b
    ps_vert = repmat(ps_orig, [1 1 1 length(grid.dim3.a)]); % dims (lon x lat x time x plev)
    ps_vert = permute(ps_vert, [1 2 4 3]); % dims (lon x lat x plev x time)
    a = permute(repmat(grid.dim3.a, [1 size(ps_orig)]), [2 3 1 4]);
    b = permute(repmat(grid.dim3.b, [1 size(ps_orig)]), [2 3 1 4]);
    pa = a + b.*ps_vert;

    pb = CmdLineProgressBar("Calculaing mlev...");
    for lo=1:size(pa,2)
        pb.print(lo, size(pa,2));
        for la=1:size(pa,3)
            for mo=1:size(pa,4)
                pa_si(:,lo,la,mo) = interp1(pa(:,lo,la,mo), grid.dim3.si, 'linear');


                % only keep nonnan data and redo interpolation
                notnan = find(~isnan(squeeze(pa_si(:,lo,la,mo))));
                pa_si(:,lo,la,mo) = interp1(grid.dim3.si(notnan), pa_si(notnan,lo,la,mo), grid.dim3.si);

            end
        end
    end

    pa_si = permute(pa_si, [2 3 1 4]); % reorder to lon x lat x si x mon

    if any(strcmp(type, {'era5', 'erai'})); newdir=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s', type);
    elseif strcmp(type, 'gcm'); newdir=sprintf('/project2/tas1/miyawaki/projects/002/data/read/gcm/%s/%s', par.model, par.gcm.clim);
    elseif strcmp(type, 'echam'); newdir=sprintf('/project2/tas1/miyawaki/projects/002/data/read/echam/%s', par.echam.clim);
    elseif any(strcmp(type, {'echam_ml','echam_pl'})); newdir=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s', type); end;
    if ~exist(newdir, 'dir'); mkdir(newdir); end
    filename='pa_si.mat';
    save(sprintf('%s/%s', newdir, filename), 'pa_si', '-v7.3');
end
function make_dtdz(type, par)
    if any(strcmp(type, {'era5', 'erai'}))
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s', type);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s', type);
        file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/%s/temp/%s_temp_%s.ymonmean.nc', type, type, par.(type).yr_span));
        fullpath=sprintf('%s/%s', file.folder, file.name);
        temp = ncread(fullpath, 't');
        file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/%s/zg/%s_zg_%s.ymonmean.nc', type, type, par.(type).yr_span));
        fullpath=sprintf('%s/%s', file.folder, file.name);
        zg = ncread(fullpath, 'z'); zg = zg/par.g; % convert from geopotential to height in m
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
        file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/echam/ATM_rjg_%s_*.ymonmean.nc', par.echam.clim));
        fullpath=sprintf('%s/%s', file.folder, file.name);
        temp = double(ncread(fullpath, var));
        var = 'geopoth';
        file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/echam/ATM_rjg_%s_*.ymonmean.nc', par.echam.clim));
        fullpath=sprintf('%s/%s', file.folder, file.name);
        zg = double(ncread(fullpath, var));
    end

    load(sprintf('%s/grid.mat', prefix)); % read grid data

    dtdz_half = -1e3*(temp(:,:,2:end,:)-temp(:,:,1:end-1,:))./(zg(:,:,2:end,:)-zg(:,:,1:end-1,:)); % calculate lapse rate in K/km
    p_half = 1/2*(grid.dim3.plev(2:end)+grid.dim3.plev(1:end-1));

    dtdz_half = permute(dtdz_half, [3 1 2 4]); % bring height front
    dtdz = interp1(p_half, dtdz_half, par.pa);
    dtdz = permute(dtdz, [2 3 1 4]); % bring height back to 3rd

    if any(strcmp(type, {'era5', 'erai'})); newdir=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s', type);
    elseif strcmp(type, 'gcm'); newdir=sprintf('/project2/tas1/miyawaki/projects/002/data/read/gcm/%s/%s', par.model, par.gcm.clim);
    elseif strcmp(type, 'echam'); newdir=sprintf('/project2/tas1/miyawaki/projects/002/data/read/echam/%s', par.echam.clim); end;
    if ~exist(newdir, 'dir'); mkdir(newdir); end
    filename='dtdz.mat';
    save(sprintf('%s/%s', newdir, filename), 'dtdz', '-v7.3');
end
function make_dtdz_z(type, par)
    if strcmp(type, 'era5') | strcmp(type, 'erai')
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s', type);
    elseif strcmp(type, 'gcm')
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s/%s', type, par.model, par.gcm.clim);
    elseif strcmp(type, 'echam')
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s', type, par.echam.clim);
    end
    load(sprintf('%s/grid.mat', prefix)); % read grid data
    load(sprintf('%s/tempz.mat', prefix)); % read temp in z coordinates

    dtdz_z = nan(size(tempz));

    tempz = permute(tempz, [3 1 2 4]); % bring height forward
    dtdz_z = -1e3*(tempz(2:end,:,:,:)-tempz(1:end-1,:,:,:))./repmat(grid.dim3.z(2:end)-grid.dim3.z(1:end-1), [1,size(tempz,2),size(tempz,3),size(tempz,4)]); % lapse rate in K/km
    dtdz_z = interp1(1/2*(grid.dim3.z(2:end)+grid.dim3.z(1:end-1)), dtdz_z, grid.dim3.z);

    dtdz_z = permute(dtdz_z, [2 3 1 4]); % bring height back to 3rd dimension

    if any(strcmp(type, {'era5', 'erai'})); newdir=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s', type);
    elseif strcmp(type, 'gcm'); newdir=sprintf('/project2/tas1/miyawaki/projects/002/data/read/gcm/%s/%s', par.model, par.gcm.clim);
    elseif strcmp(type, 'echam'); newdir=sprintf('/project2/tas1/miyawaki/projects/002/data/read/echam/%s', par.echam.clim); end;
    if ~exist(newdir, 'dir'); mkdir(newdir); end
    filename='dtdz_z.mat';
    save(sprintf('%s/%s', newdir, filename), 'dtdz_z', '-v7.3');
end
function make_ztrop(type, par) % calculate WMO tropopause
    if strcmp(type, 'era5') | strcmp(type, 'erai')
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s', type);
    elseif strcmp(type, 'gcm')
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s/%s', type, par.model, par.gcm.clim);
    elseif strcmp(type, 'echam')
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s', type, par.echam.clim);
    end
    load(sprintf('%s/grid.mat', prefix)); % read grid data
    load(sprintf('%s/tempz.mat', prefix)); % read temp in z coordinates

    ga = nan(size(tempz));
    ztrop = nan([size(tempz,1),size(tempz,2),size(tempz,4)]);

    tempz = permute(tempz, [3 1 2 4]);
    ga = -1e3*(tempz(2:end,:,:,:)-tempz(1:end-1,:,:,:))./repmat(grid.dim3.z(2:end)-grid.dim3.z(1:end-1), [1,size(tempz,2),size(tempz,3),size(tempz,4)]); % lapse rate in K/km
    ga = interp1(1/2*(grid.dim3.z(2:end)+grid.dim3.z(1:end-1)), ga, grid.dim3.z);

    clear tempz

    pb=CmdLineProgressBar("Calculating ztrop..."); % track progress of this loop
    for lo = 1:length(grid.dim3.lon)
        pb.print(lo, length(grid.dim3.lon));
        for la = 1:length(grid.dim3.lat)
            for mo = 1:12
                ga_hires = interp1(grid.dim3.z, ga(:,lo,la,mo), par.z_hires); % make high-resolution grid for computing tropopause
                cand = squeeze(ga_hires<=2); % all levels where lapse rate is less than 2 K/km
                if any(cand)
                    idx_cand = find(cand); % indices of candidates
                    i = 1;
                    while isnan(ztrop(lo,la,mo)) & i<=length(idx_cand)
                        idx = idx_cand(i);
                        ztrop_tmp = par.z_hires(idx);
                        ga_2km = nanmean(interp1(par.z_hires, ga_hires, ztrop_tmp:100:ztrop_tmp+2e3)); % compute average lapse rate from this point to 2 km above it
                        if ga_2km < 2
                            ztrop(lo,la,mo) = ztrop_tmp;
                        end
                        i=i+1;
                    end
                end
            end
        end
    end

    if any(strcmp(type, {'era5', 'erai'})); newdir=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s', type);
    elseif strcmp(type, 'gcm'); newdir=sprintf('/project2/tas1/miyawaki/projects/002/data/read/gcm/%s/%s', par.model, par.gcm.clim);
    elseif strcmp(type, 'echam'); newdir=sprintf('/project2/tas1/miyawaki/projects/002/data/read/echam/%s', par.echam.clim); end;
    if ~exist(newdir, 'dir'); mkdir(newdir); end
    filename='ztrop.mat';
    save(sprintf('%s/%s', newdir, filename), 'ztrop', '-v7.3');
end
function make_ztrop_z(type, par) % calculate WMO tropopause of latitudinally-averaged data
    if strcmp(type, 'era5') | strcmp(type, 'erai')
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s', type);
    elseif strcmp(type, 'gcm')
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s/%s', type, par.model, par.gcm.clim);
    elseif strcmp(type, 'echam')
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s', type, par.echam.clim);
    end
    load(sprintf('%s/grid.mat', prefix)); % read grid data
    load(sprintf('%s/tempz.mat', prefix)); % read temp in z coordinates

    % zonal average
    tempz = squeeze(nanmean(tempz, 1));

    ga = nan(size(tempz));
    ztrop_z = nan([size(tempz,1),size(tempz,3)]);

    tempz = permute(tempz, [2 1 3]); % bring levels to the front
    tempz = fillmissing(tempz, 'nearest');
    ga = -1e3*(tempz(2:end,:,:)-tempz(1:end-1,:,:))./repmat(grid.dim3.z(2:end)-grid.dim3.z(1:end-1), [1,size(tempz,2),size(tempz,3)]); % lapse rate in K/km
    ga = interp1(1/2*(grid.dim3.z(2:end)+grid.dim3.z(1:end-1)), ga, grid.dim3.z);

    clear tempz

    pb=CmdLineProgressBar("Calculating ztrop_z..."); % track progress of this loop
    for la = 1:length(grid.dim3.lat)
        pb.print(la, length(grid.dim3.lat));
        for mo = 1:12
            ga_hires = interp1(grid.dim3.z, ga(:,la,mo), par.z_hires); % make high-resolution grid for computing tropopause
            cand = squeeze(ga_hires<=2); % all levels where lapse rate is less than 2 K/km
            if any(cand)
                idx_cand = find(cand); % indices of candidates
                i = 1;
                while isnan(ztrop_z(la,mo)) & i<=length(idx_cand)
                    idx = idx_cand(i);
                    ztrop_tmp = par.z_hires(idx);
                    ga_2km = nanmean(interp1(par.z_hires, ga_hires, ztrop_tmp:100:ztrop_tmp+2e3)); % compute average lapse rate from this point to 2 km above it
                    if ga_2km < 2
                        ztrop_z(la,mo) = ztrop_tmp;
                    end
                    i=i+1;
                end
            end
        end
    end

    if any(strcmp(type, {'era5', 'erai'})); newdir=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s', type);
    elseif strcmp(type, 'gcm'); newdir=sprintf('/project2/tas1/miyawaki/projects/002/data/read/gcm/%s/%s', par.model, par.gcm.clim);
    elseif strcmp(type, 'echam'); newdir=sprintf('/project2/tas1/miyawaki/projects/002/data/read/echam/%s', par.echam.clim); end;
    if ~exist(newdir, 'dir'); mkdir(newdir); end
    filename='ztrop_z.mat';
    save(sprintf('%s/%s', newdir, filename), 'ztrop_z', '-v7.3');
end
function make_ptrop_convert(type, par) % calculate WMO tropopause
    if strcmp(type, 'era5') | strcmp(type, 'erai')
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s', type);
        file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/%s/zg/%s_zg_%s.ymonmean.nc', type, type, par.(type).yr_span));
        fullpath=sprintf('%s/%s', file.folder, file.name);
        zg = ncread(fullpath, 'z')/par.g;
    elseif strcmp(type, 'gcm')
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s/%s', type, par.model, par.gcm.clim);
        var = 'zg';
        file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/gcm/%s/%s_Amon_%s_%s_r1i1p1_*.ymonmean.nc', par.model, var, par.model, par.gcm.clim));
        fullpath=sprintf('%s/%s', file.folder, file.name);
        zg = ncread(fullpath, var);
    elseif strcmp(type, 'echam')
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s', type, par.echam.clim);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s', type, par.echam.clim);
        var = 'geopoth';
        file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/echam/ATM_rjg_%s_*.ymonmean.nc', par.echam.clim));
        fullpath=sprintf('%s/%s', file.folder, file.name);
        zg = double(ncread(fullpath, var));
    end
    load(sprintf('%s/grid.mat', prefix)); % read grid data
    load(sprintf('%s/ztrop.mat', prefix)); % read z tropopause data

    ptrop = nan(size(ztrop));

    pb=CmdLineProgressBar("Calculating ptrop..."); % track progress of this loop
    for lo = 1:length(grid.dim3.lon)
        pb.print(lo, length(grid.dim3.lon));
        for la = 1:length(grid.dim3.lat)
            for mo = 1:12
                ptrop(lo,la,mo) = interp1(squeeze(zg(lo,la,:,mo)), grid.dim3.plev, ztrop(lo,la,mo)); % make high-resolution grid for computing tropopause
            end
        end
    end

    if any(strcmp(type, {'era5', 'erai'})); newdir=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s', type);
    elseif strcmp(type, 'gcm'); newdir=sprintf('/project2/tas1/miyawaki/projects/002/data/read/gcm/%s/%s', par.model, par.gcm.clim);
    elseif strcmp(type, 'echam'); newdir=sprintf('/project2/tas1/miyawaki/projects/002/data/read/echam/%s', par.echam.clim); end;
    if ~exist(newdir, 'dir'); mkdir(newdir); end
    filename='ptrop.mat';
    save(sprintf('%s/%s', newdir, filename), 'ptrop', '-v7.3');
end
function make_ptrop_z_convert(type, par) % calculate WMO tropopause
    if strcmp(type, 'era5') | strcmp(type, 'erai')
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s', type);
        file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/%s/zg/%s_zg_%s.ymonmean.nc', type, type, par.(type).yr_span));
        fullpath=sprintf('%s/%s', file.folder, file.name);
        zg = ncread(fullpath, 'z')/par.g;
    elseif strcmp(type, 'gcm')
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s/%s', type, par.model, par.gcm.clim);
        var = 'zg';
        file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/gcm/%s/%s_Amon_%s_%s_r1i1p1_*.ymonmean.nc', par.model, var, par.model, par.gcm.clim));
        fullpath=sprintf('%s/%s', file.folder, file.name);
        zg = ncread(fullpath, var);
    elseif strcmp(type, 'echam')
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s', type, par.echam.clim);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s', type, par.echam.clim);
        var = 'geopoth';
        file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/echam/ATM_rjg_%s_*.ymonmean.nc', par.echam.clim));
        fullpath=sprintf('%s/%s', file.folder, file.name);
        zg = double(ncread(fullpath, var));
    end
    load(sprintf('%s/grid.mat', prefix)); % read grid data
    load(sprintf('%s/ztrop_z.mat', prefix)); % read z tropopause data

    zg_z = squeeze(nanmean(zg,1)); clear zg;
    ptrop_z = nan(size(ztrop_z));

    pb=CmdLineProgressBar("Calculating ptrop..."); % track progress of this loop
    for la = 1:length(grid.dim3.lat)
    pb.print(la, length(grid.dim3.lat));
        for mo = 1:12
            ptrop_z(la,mo) = interp1(squeeze(zg_z(la,:,mo)), grid.dim3.plev, ztrop_z(la,mo)); % make high-resolution grid for computing tropopause
        end
    end

    if any(strcmp(type, {'era5', 'erai'})); newdir=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s', type);
    elseif strcmp(type, 'gcm'); newdir=sprintf('/project2/tas1/miyawaki/projects/002/data/read/gcm/%s/%s', par.model, par.gcm.clim);
    elseif strcmp(type, 'echam'); newdir=sprintf('/project2/tas1/miyawaki/projects/002/data/read/echam/%s', par.echam.clim); end;
    if ~exist(newdir, 'dir'); mkdir(newdir); end
    filename='ptrop_z.mat';
    save(sprintf('%s/%s', newdir, filename), 'ptrop_z', '-v7.3');
end
function make_ptrop(type, par) % calculate WMO tropopause
    if strcmp(type, 'era5') | strcmp(type, 'erai')
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s', type);
        file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/%s/temp/%s_temp_%s.ymonmean.nc', type, type, par.(type).yr_span));
        fullpath=sprintf('%s/%s', file.folder, file.name);
        temp = ncread(fullpath, 't');
        file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/%s/zg/%s_zg_%s.ymonmean.nc', type, type, par.(type).yr_span));
        fullpath=sprintf('%s/%s', file.folder, file.name);
        zg = ncread(fullpath, 'z')/par.g;
    elseif strcmp(type, 'gcm')
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s/%s', type, par.model, par.gcm.clim);
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
        file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/echam/ATM_rjg_%s_*.ymonmean.nc', par.echam.clim));
        fullpath=sprintf('%s/%s', file.folder, file.name);
        temp = double(ncread(fullpath, var));
        var = 'geopoth';
        file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/echam/ATM_rjg_%s_*.ymonmean.nc', par.echam.clim));
        fullpath=sprintf('%s/%s', file.folder, file.name);
        zg = double(ncread(fullpath, var));
    end
    load(sprintf('%s/grid.mat', prefix)); % read grid data
    load(sprintf('%s/ztrop.mat', prefix)); % read z tropopause data

    ga = nan(size(temp));
    ptrop = nan([size(temp,1),size(temp,2),size(temp,4)]);

    temp = permute(temp, [3 1 2 4]); % bring levels to the front
    zg = permute(zg, [3 1 2 4]); % bring levels to the front
    temp = fillmissing(temp, 'nearest');
    ga = -1e3*(temp(2:end,:,:,:)-temp(1:end-1,:,:,:))./(zg(2:end,:,:,:)-zg(1:end-1,:,:,:)); % lapse rate in K/km

    clear temp

    pb=CmdLineProgressBar("Calculating ptrop..."); % track progress of this loop
    for lo = 1:length(grid.dim3.lon)
        pb.print(lo, length(grid.dim3.lon));
        for la = 1:length(grid.dim3.lat)
            for mo = 1:12
                ga_hires = interp1(1/2*(zg(2:end,lo,la,mo)+zg(1:end-1,lo,la,mo)), ga(:,lo,la,mo), zg(:,lo,la,mo));

                cand = squeeze(ga_hires<=2) & grid.dim3.plev<=500e2; % all levels where lapse rate is less than 2 K/km and above 500 hPa
                if any(cand)
                    idx_cand = find(cand); % indices of candidates
                    i = 1;
                    while isnan(ptrop(lo,la,mo)) & i<=length(idx_cand)
                        idx = idx_cand(i);
                        ztrop_tmp = zg(idx,lo,la,mo);
                        ptrop_tmp = grid.dim3.plev(idx);
                        ga_2km = nanmean(interp1(zg(:,lo,la,mo), ga_hires, ztrop_tmp:100:ztrop_tmp+2e3)); % compute average lapse rate from this point to 2 km above it
                        if ga_2km < 2
                            ptrop(lo,la,mo) = ptrop_tmp;
                        end
                        i=i+1;
                    end
                end
            end
        end
    end

    if any(strcmp(type, {'era5', 'erai'})); newdir=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s', type);
    elseif strcmp(type, 'gcm'); newdir=sprintf('/project2/tas1/miyawaki/projects/002/data/read/gcm/%s/%s', par.model, par.gcm.clim);
    elseif strcmp(type, 'echam'); newdir=sprintf('/project2/tas1/miyawaki/projects/002/data/read/echam/%s', par.echam.clim); end;
    if ~exist(newdir, 'dir'); mkdir(newdir); end
    filename='ptrop.mat';
    save(sprintf('%s/%s', newdir, filename), 'ptrop', '-v7.3');
end
function make_ptrop_z(type, par) % calculate WMO tropopause
    if strcmp(type, 'era5') | strcmp(type, 'erai')
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s', type);
        file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/%s/temp/%s_temp_%s.ymonmean.nc', type, type, par.(type).yr_span));
        fullpath=sprintf('%s/%s', file.folder, file.name);
        temp = ncread(fullpath, 't');
        file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/%s/zg/%s_zg_%s.ymonmean.nc', type, type, par.(type).yr_span));
        fullpath=sprintf('%s/%s', file.folder, file.name);
        zg = ncread(fullpath, 'z')/par.g;
    elseif strcmp(type, 'gcm')
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s/%s', type, par.model, par.gcm.clim);
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
        file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/echam/ATM_rjg_%s_*.ymonmean.nc', par.echam.clim));
        fullpath=sprintf('%s/%s', file.folder, file.name);
        temp = double(ncread(fullpath, var));
        var = 'geopoth';
        file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/echam/ATM_rjg_%s_*.ymonmean.nc', par.echam.clim));
        fullpath=sprintf('%s/%s', file.folder, file.name);
        zg = double(ncread(fullpath, var));
    end
    load(sprintf('%s/grid.mat', prefix)); % read grid data

    temp_z = squeeze(nanmean(temp,1)); clear temp;
    zg_z = squeeze(nanmean(zg,1)); clear zg;

    ga = nan(size(temp_z));
    ptrop_z = nan([size(temp_z,1),size(temp_z,3)]);

    temp_z = permute(temp_z, [2 1 3]); % bring levels to the front
    zg_z = permute(zg_z, [2 1 3]); % bring levels to the front
    temp_z = fillmissing(temp_z, 'nearest');
    ga = -1e3*(temp_z(2:end,:,:)-temp_z(1:end-1,:,:))./(zg_z(2:end,:,:)-zg_z(1:end-1,:,:)); % lapse rate in K/km

    clear temp_z

    pb=CmdLineProgressBar("Calculating ptrop..."); % track progress of this loop
    for la = 1:length(grid.dim3.lat)
    pb.print(la, length(grid.dim3.lat));
        for mo = 1:12
            ga_hires = interp1(1/2*(zg_z(2:end,la,mo)+zg_z(1:end-1,la,mo)), ga(:,la,mo), zg_z(:,la,mo));

            cand = squeeze(ga_hires<=2) & grid.dim3.plev<=500e2; % all levels where lapse rate is less than 2 K/km and above 500 hPa
            if any(cand)
                idx_cand = find(cand); % indices of candidates
                i = 1;
                while isnan(ptrop_z(la,mo)) & i<=length(idx_cand)
                    idx = idx_cand(i);
                    ztrop_tmp = zg_z(idx,la,mo);
                    ptrop_tmp = grid.dim3.plev(idx);
                    ga_2km = nanmean(interp1(zg_z(:,la,mo), ga_hires, ztrop_tmp:100:ztrop_tmp+2e3)); % compute average lapse rate from this point to 2 km above it
                    if ga_2km < 2
                        ptrop_z(la,mo) = ptrop_tmp;
                    end
                    i=i+1;
                end
            end
        end
    end

    if any(strcmp(type, {'era5', 'erai'})); newdir=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s', type);
    elseif strcmp(type, 'gcm'); newdir=sprintf('/project2/tas1/miyawaki/projects/002/data/read/gcm/%s/%s', par.model, par.gcm.clim);
    elseif strcmp(type, 'echam'); newdir=sprintf('/project2/tas1/miyawaki/projects/002/data/read/echam/%s', par.echam.clim); end;
    if ~exist(newdir, 'dir'); mkdir(newdir); end
    filename='ptrop_z.mat';
    save(sprintf('%s/%s', newdir, filename), 'ptrop_z', '-v7.3');
end
function make_sitrop(type, par) % calculate WMO tropopause
    if strcmp(type, 'era5') | strcmp(type, 'erai')
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s', type);
    elseif strcmp(type, 'gcm')
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s/%s', type, par.model, par.gcm.clim);
    elseif strcmp(type, 'echam')
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s', type, par.echam.clim);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s', type, par.echam.clim);
    end
    load(sprintf('%s/grid.mat', prefix)); % read grid data
    load(sprintf('%s/ptrop.mat', prefix)); % read z tropopause data

    ptrop = nan(size(ztrop));

    pb=CmdLineProgressBar("Calculating ptrop..."); % track progress of this loop
    for lo = 1:length(grid.dim3.lon)
        pb.print(lo, length(grid.dim3.lon));
        for la = 1:length(grid.dim3.lat)
            for mo = 1:12
                ptrop(lo,la,mo) = interp1(squeeze(zg(lo,la,:,mo)), grid.dim3.plev, ztrop(lo,la,mo)); % make high-resolution grid for computing tropopause
            end
        end
    end

    if any(strcmp(type, {'era5', 'erai'})); newdir=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s', type);
    elseif strcmp(type, 'gcm'); newdir=sprintf('/project2/tas1/miyawaki/projects/002/data/read/gcm/%s/%s', par.model, par.gcm.clim);
    elseif strcmp(type, 'echam'); newdir=sprintf('/project2/tas1/miyawaki/projects/002/data/read/echam/%s', par.echam.clim); end;
    if ~exist(newdir, 'dir'); mkdir(newdir); end
    filename='ptrop.mat';
    save(sprintf('%s/%s', newdir, filename), 'ptrop', '-v7.3');
end
function make_sitrop_z(type, par) % calculate WMO tropopause
    if strcmp(type, 'era5') | strcmp(type, 'erai')
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s', type);
    elseif strcmp(type, 'gcm')
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s/%s', type, par.model, par.gcm.clim);
    elseif strcmp(type, 'echam')
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s', type, par.echam.clim);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s', type, par.echam.clim);
    end
    load(sprintf('%s/grid.mat', prefix)); % read grid data
    load(sprintf('%s/ptrop_z.mat', prefix)); % read z tropopause data

    zg_z = squeeze(nanmean(zg,1)); clear zg;
    ptrop_z = nan(size(ztrop_z));

    pb=CmdLineProgressBar("Calculating ptrop..."); % track progress of this loop
    for la = 1:length(grid.dim3.lat)
    pb.print(la, length(grid.dim3.lat));
        for mo = 1:12
            ptrop_z(la,mo) = interp1(squeeze(zg_z(la,:,mo)), grid.dim3.plev, ztrop_z(la,mo)); % make high-resolution grid for computing tropopause
        end
    end

    if any(strcmp(type, {'era5', 'erai'})); newdir=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s', type);
    elseif strcmp(type, 'gcm'); newdir=sprintf('/project2/tas1/miyawaki/projects/002/data/read/gcm/%s/%s', par.model, par.gcm.clim);
    elseif strcmp(type, 'echam'); newdir=sprintf('/project2/tas1/miyawaki/projects/002/data/read/echam/%s', par.echam.clim); end;
    if ~exist(newdir, 'dir'); mkdir(newdir); end
    filename='ptrop_z.mat';
    save(sprintf('%s/%s', newdir, filename), 'ptrop_z', '-v7.3');
end
function make_alb(type, par) % calculate surface albedo
    if any(strcmp(type, {'era5', 'erai'}))
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s', type);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s', type);
    elseif strcmp(type, 'gcm')
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s/%s', type, par.model, par.gcm.clim);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s', type, par.model, par.gcm.clim);
    end
    load(sprintf('%s/grid.mat', prefix)); % read grid data
    load(sprintf('%s/rad.mat', prefix)); % read radiation data

    if strcmp(type, 'erai') | strcmp(type, 'era5')
        alb = rad.rsus./rad.rsds; % TODO download shortwave radiation up and down separately for ERA
    elseif strcmp(type, 'gcm')
        alb = rad.rsus./rad.rsds;
    end

    if any(strcmp(type, {'era5', 'erai'})); newdir=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s', type);
    elseif strcmp(type, 'gcm'); newdir=sprintf('/project2/tas1/miyawaki/projects/002/data/read/gcm/%s/%s', par.model, par.gcm.clim); end;
    if ~exist(newdir, 'dir'); mkdir(newdir); end
    filename='alb.mat';
    save(sprintf('%s/%s', newdir, filename), 'alb', '-v7.3');

end
function make_albcs(type, par) % calculate clear sky surface albedo
    if any(strcmp(type, {'era5', 'erai'}))
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s', type);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s', type);
    elseif strcmp(type, 'gcm')
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s/%s', type, par.model, par.gcm.clim);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s', type, par.model, par.gcm.clim);
    end
    load(sprintf('%s/grid.mat', prefix)); % read grid data
    load(sprintf('%s/radcs.mat', prefix)); % read radcsiation data

    if strcmp(type, 'erai') | strcmp(type, 'era5')
        albcs = radcs.rsus./radcs.rsds; % TODO download shortwave radcsiation up and down separately for ERA
    elseif strcmp(type, 'gcm')
        albcs = radcs.rsuscs./radcs.rsdscs;
    end

    if any(strcmp(type, {'era5', 'erai'})); newdir=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s', type);
    elseif strcmp(type, 'gcm'); newdir=sprintf('/project2/tas1/miyawaki/projects/002/data/read/gcm/%s/%s', par.model, par.gcm.clim); end;
    if ~exist(newdir, 'dir'); mkdir(newdir); end
    filename='albcs.mat';
    save(sprintf('%s/%s', newdir, filename), 'albcs', '-v7.3');

end
function make_palb(type, par) % calculate planetary albedo
    if any(strcmp(type, {'era5', 'erai'}))
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s', type);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s', type);
    elseif strcmp(type, 'gcm')
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s/%s', type, par.model, par.gcm.clim);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s', type, par.model, par.gcm.clim);
    end
    load(sprintf('%s/grid.mat', prefix)); % read grid data
    load(sprintf('%s/rad.mat', prefix)); % read radiation data

    if strcmp(type, 'erai') | strcmp(type, 'era5')
        palb = rad.rsus./rad.rsds; % TODO download shortwave radiation up and down separately for ERA
    elseif strcmp(type, 'gcm')
        palb = rad.rsut./rad.rsdt;
    end

    if any(strcmp(type, {'era5', 'erai'})); newdir=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s', type);
    elseif strcmp(type, 'gcm'); newdir=sprintf('/project2/tas1/miyawaki/projects/002/data/read/gcm/%s/%s', par.model, par.gcm.clim); end;
    if ~exist(newdir, 'dir'); mkdir(newdir); end
    filename='palb.mat';
    save(sprintf('%s/%s', newdir, filename), 'palb', '-v7.3');

end
