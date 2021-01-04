function make_tempz(type, par)
    if any(strcmp(type, {'era5', 'era5c', 'erai'}))
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s', type);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s', type);
        file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/%s/temp/%s_temp_%s.ymonmean.nc', type, type, par.(type).yr_span));
        fullpath=sprintf('%s/%s', file.folder, file.name);
        temp = double(ncread(fullpath, 't'));
        file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/%s/zg/%s_zg_%s.ymonmean.nc', type, type, par.(type).yr_span));
        fullpath=sprintf('%s/%s', file.folder, file.name);
        zg = double(ncread(fullpath, 'z')); zg = zg/par.g; % convert from geopotential to height in m
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
    elseif strcmp(type, 'echam_ml')
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s', type);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s', type);
        var = 't';
        file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/echam_ml/ATM*.ymonmean.nc'));
        fullpath=sprintf('%s/%s', file.folder, file.name);
        temp = double(ncread(fullpath, var));
        var = 'geopoth';
        file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/echam_ml/ATM*.ymonmean.nc'));
        fullpath=sprintf('%s/%s', file.folder, file.name);
        zg = double(ncread(fullpath, var));
    elseif strcmp(type, 'echam_pl')
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s', type);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s', type);
        var = 't';
        file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/echam_pl/ATM*.ymonmean.nc'));
        fullpath=sprintf('%s/%s', file.folder, file.name);
        temp = double(ncread(fullpath, var));
        var = 'geopoth';
        file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/echam_pl/ATM*.ymonmean.nc'));
        fullpath=sprintf('%s/%s', file.folder, file.name);
        zg = double(ncread(fullpath, var));
    end

    load(sprintf('%s/grid.mat', prefix)); % read grid data

    temp = permute(temp, [3 1 2 4]);
    zg = permute(zg, [3 1 2 4]);

    pb=CmdLineProgressBar("Calculating tempz..."); % track progress of this loop
    for lo = 1:length(grid.dim3.lon)
        pb.print(lo, length(grid.dim3.lon));
        for la = 1:length(grid.dim3.lat)
            for mo = 1:12
                % only keep nonnan data and interpolate
                notnan = find(~isnan(zg(:,lo,la,mo)));

                tempz(:,lo,la,mo) = interp1(zg(notnan,lo,la,mo), temp(notnan,lo,la,mo), grid.dim3.z);
            end
        end
    end

    tempz = permute(tempz, [2 3 1 4]); % reorder to lon x lat x z x mon

    if any(strcmp(type, {'era5', 'era5c', 'erai'})); newdir=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s', type);
    elseif strcmp(type, 'gcm'); newdir=sprintf('/project2/tas1/miyawaki/projects/002/data/read/gcm/%s/%s', par.model, par.gcm.clim);
    elseif strcmp(type, 'echam'); newdir=sprintf('/project2/tas1/miyawaki/projects/002/data/read/echam/%s', par.echam.clim);
    elseif strcmp(type, 'echam_ml'); newdir=sprintf('/project2/tas1/miyawaki/projects/002/data/read/echam_ml');
    elseif strcmp(type, 'echam_pl'); newdir=sprintf('/project2/tas1/miyawaki/projects/002/data/read/echam_pl'); end;
    if ~exist(newdir, 'dir'); mkdir(newdir); end
    filename='tempz.mat';
    save(sprintf('%s/%s', newdir, filename), 'tempz', '-v7.3');
end
