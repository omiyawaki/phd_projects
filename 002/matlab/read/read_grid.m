function read_grid(type, par)
    filename = 'grid.mat';
    % read data net SW and LW radiation data downloaded from Era5
    % first read lon and lat vectors since this is different from the Donohoe grid
    if any(strcmp(type, {'era5', 'erai'}))
        grid.dim2.lon = double(ncread(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/%s/rad/%s_rad_%s.ymonmean.nc', type, type, par.(type).yr_span), 'longitude'));
        grid.dim2.lat = double(ncread(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/%s/rad/%s_rad_%s.ymonmean.nc', type, type, par.(type).yr_span), 'latitude'));
        grid.dim3 = grid.dim2;
        grid.dim3.plev = 10^2*double(ncread(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/%s/temp/%s_temp_%s.ymonmean.nc', type, type, par.(type).yr_span), 'level')); % multiply by 100 to convert hPa to Pa
        grid.dim3.z = par.z;
        grid.dim3.si = 1e-5*par.pa;
        newdir = sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s', type, par.(type).yr_span);
        if ~exist(newdir, 'dir'); mkdir(newdir); end
        save(sprintf('%s/%s', newdir, filename), 'grid')
    elseif any(strcmp(type, {'era5c'}))
        grid.dim2.lon = double(ncread(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/%s/rad/%s_rad_%s.ymonmean.nc', type, type, par.(type).yr_span), 'lon'));
        grid.dim2.lat = double(ncread(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/%s/rad/%s_rad_%s.ymonmean.nc', type, type, par.(type).yr_span), 'lat'));
        grid.dim3 = grid.dim2;
        grid.dim3.plev = 10^2*double(ncread(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/%s/temp/%s_temp_%s.ymonmean.nc', type, type, par.(type).yr_span), 'level')); % multiply by 100 to convert hPa to Pa
        grid.dim3.z = par.z;
        grid.dim3.si = 1e-5*par.pa;
        newdir = sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s', type, par.(type).yr_span);
        if ~exist(newdir, 'dir'); mkdir(newdir); end
        save(sprintf('%s/%s', newdir, filename), 'grid')
    elseif strcmp(type, 'merra2')
        grid.dim2.lon = double(ncread(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/%s/rad/%s_rad_%s.ymonmean.nc', type, type, par.(type).yr_span), 'lon'));
        grid.dim2.lat = double(ncread(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/%s/rad/%s_rad_%s.ymonmean.nc', type, type, par.(type).yr_span), 'lat'));
        grid.dim3 = grid.dim2;
        grid.dim3.plev = 1e2*double(ncread(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/%s/temp/%s_temp_%s.ymonmean.nc', type, type, par.(type).yr_span), 'lev')); % multiply by 100 to convert hPa to Pa
        grid.dim3.z = par.z;
        grid.dim3.si = 1e-5*par.pa;
        newdir = sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s', type, par.(type).yr_span);
        if ~exist(newdir, 'dir'); mkdir(newdir); end
        save(sprintf('%s/%s', newdir, filename), 'grid')
    elseif strcmp(type, 'jra55')
        grid.dim2.lon = double(ncread(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/%s/rad/%s_dswrf_%s.ymonmean.nc', type, type, par.(type).yr_span), 'g0_lon_2'));
        grid.dim2.lat = double(ncread(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/%s/rad/%s_dswrf_%s.ymonmean.nc', type, type, par.(type).yr_span), 'g0_lat_1'));
        grid.dim3 = grid.dim2;
        grid.dim3.plev = 1e2*double(ncread(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/%s/temp/%s_tmp_%s.ymonmean.nc', type, type, par.(type).yr_span), 'lv_ISBL1')); % multiply by 100 to convert hPa to Pa
        grid.dim3.z = par.z;
        grid.dim3.si = 1e-5*par.pa;
        newdir = sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s', type, par.(type).yr_span);
        if ~exist(newdir, 'dir'); mkdir(newdir); end
        save(sprintf('%s/%s', newdir, filename), 'grid')
    elseif strcmp(type, 'gcm')
        file.dim2=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/gcm/%s/%s_Amon_%s_%s_r1i1p1_*.nc', par.model, 'tas', par.model, par.gcm.clim));
        file.dim3=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/gcm/%s/%s_Amon_%s_%s_r1i1p1_*.nc', par.model, 'ta', par.model, par.gcm.clim));
        fullpath.dim2=sprintf('%s/%s', file.dim2.folder, file.dim2.name);
        fullpath.dim3=sprintf('%s/%s', file.dim3.folder, file.dim3.name);
        grid.dim2.lon=double(ncread(fullpath.dim2, 'lon'));
        grid.dim3.lon=double(ncread(fullpath.dim3, 'lon'));
        grid.dim2.lat=double(ncread(fullpath.dim2, 'lat'));
        grid.dim3.lat=double(ncread(fullpath.dim3, 'lat'));
        grid.dim3.plev=double(ncread(fullpath.dim3, 'plev'));
        grid.dim3.z = par.z;
        grid.dim3.si = 1e-5*par.pa;
        newdir=sprintf('/project2/tas1/miyawaki/projects/002/data/read/gcm/%s/%s', par.model, par.gcm.clim);
        if ~exist(newdir, 'dir'); mkdir(newdir); end
        save(sprintf('%s/%s', newdir, filename), 'grid');
    elseif strcmp(type, 'echam')
        if contains(par.echam.clim, 'echr000') | any(strcmp(par.echam.clim, par.echam.exceptions))
            file.dim2=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/echam/BOT*_%s_*.ymonmean.nc', par.echam.clim));
            file.dim3=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/echam/ATM*_%s_*.ymonmean.nc', par.echam.clim));
        else
            file.dim2=dir(sprintf('/project2/tas1/ockham/data11/tas/echam-aiv_rcc_6.1.00p1/%s/BOT_%s_0020_39.nc', par.echam.clim, par.echam.clim));
            file.dim3=dir(sprintf('/project2/tas1/ockham/data11/tas/echam-aiv_rcc_6.1.00p1/%s/ATM_%s_0020_39.nc', par.echam.clim, par.echam.clim));
        end
        fullpath.dim2=sprintf('%s/%s', file.dim2.folder, file.dim2.name);
        fullpath.dim3=sprintf('%s/%s', file.dim3.folder, file.dim3.name);
        grid.dim2.lon=double(ncread(fullpath.dim2, 'lon'));
        grid.dim3.lon=double(ncread(fullpath.dim3, 'lon'));
        grid.dim2.lat=double(ncread(fullpath.dim2, 'lat'));
        grid.dim3.lat=double(ncread(fullpath.dim3, 'lat'));
        if strcmp(par.echam.clim, 'echr0025')
                grid.dim3.plev=double(ncread(fullpath.dim3, 'plev'));
        else
                grid.dim3.plev=double(ncread(fullpath.dim3, 'lev'));
        end
        grid.dim3.z = par.z;
        grid.dim3.si = 1e-5*par.pa;
        newdir=sprintf('/project2/tas1/miyawaki/projects/002/data/read/echam/%s', par.echam.clim);
        if ~exist(newdir, 'dir'); mkdir(newdir); end
        filename='grid.mat';
        save(sprintf('%s/%s', newdir, filename), 'grid');
    elseif strcmp(type, 'hahn')
        fprefix = make_hahn_fprefix(par);
        file.dim2=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/hahn/lapserateclima/%s.TS.nc', fprefix));
        file.dim3=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/hahn/lapserateclima/%s.T.nc', fprefix));
        fullpath.dim2=sprintf('%s/%s', file.dim2.folder, file.dim2.name);
        fullpath.dim3=sprintf('%s/%s', file.dim3.folder, file.dim3.name);
        grid.dim2.lon=double(ncread(fullpath.dim2, 'lon'));
        grid.dim3.lon=double(ncread(fullpath.dim3, 'lon'));
        grid.dim2.lat=double(ncread(fullpath.dim2, 'lat'));
        grid.dim3.lat=double(ncread(fullpath.dim3, 'lat'));
        grid.dim3.plev=1e2*double(ncread(fullpath.dim3, 'lev')); % convert hPa to Pa
        grid.dim3.z = par.z;
        grid.dim3.si = 1e-5*par.pa;
        newdir=sprintf('/project2/tas1/miyawaki/projects/002/data/read/hahn/%s', par.hahn.clim);
        if ~exist(newdir, 'dir'); mkdir(newdir); end
        filename='grid.mat';
        save(sprintf('%s/%s', newdir, filename), 'grid');
    elseif strcmp(type, 'echam_ml')
        file.dim2=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/echam_ml/BOT_*.ymonmean.nc'));
        file.dim3=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/echam_ml/ATM_*.ymonmean.nc'));
        fullpath.dim2=sprintf('%s/%s', file.dim2.folder, file.dim2.name);
        fullpath.dim3=sprintf('%s/%s', file.dim3.folder, file.dim3.name);
        grid.dim2.lon=double(ncread(fullpath.dim2, 'lon'));
        grid.dim3.lon=double(ncread(fullpath.dim3, 'lon'));
        grid.dim2.lat=double(ncread(fullpath.dim2, 'lat'));
        grid.dim3.lat=double(ncread(fullpath.dim3, 'lat'));
        grid.dim3.a=double(ncread(fullpath.dim3, 'hyam'));
        grid.dim3.b=double(ncread(fullpath.dim3, 'hybm'));
        grid.dim3.z = par.z;
        grid.dim3.plev = par.pa;
        % grid.dim3.si = 1e-5*([grid.dim3.a+grid.dim3.b*1e5; 1e5]);
        grid.dim3.si = 1e-5*par.pa;
        newdir=sprintf('/project2/tas1/miyawaki/projects/002/data/read/echam_ml');
        if ~exist(newdir, 'dir'); mkdir(newdir); end
        filename='grid.mat';
        save(sprintf('%s/%s', newdir, filename), 'grid');
    elseif strcmp(type, 'echam_pl')
        file.dim2=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/echam_pl/BOT_*.ymonmean.nc'));
        file.dim3=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/echam_pl/ATM_*.ymonmean.nc'));
        fullpath.dim2=sprintf('%s/%s', file.dim2.folder, file.dim2.name);
        fullpath.dim3=sprintf('%s/%s', file.dim3.folder, file.dim3.name);
        grid.dim2.lon=double(ncread(fullpath.dim2, 'lon'));
        grid.dim3.lon=double(ncread(fullpath.dim3, 'lon'));
        grid.dim2.lat=double(ncread(fullpath.dim2, 'lat'));
        grid.dim3.lat=double(ncread(fullpath.dim3, 'lat'));
        grid.dim3.plev=double(ncread(fullpath.dim3, 'lev'));
        grid.dim3.z = par.z;
        % grid.dim3.si = 1e-5*grid.dim3.plev;
        grid.dim3.si = 1e-5*par.pa;
        newdir=sprintf('/project2/tas1/miyawaki/projects/002/data/read/echam_pl');
        if ~exist(newdir, 'dir'); mkdir(newdir); end
        filename='grid.mat';
        save(sprintf('%s/%s', newdir, filename), 'grid');
    elseif strcmp(type, 'ceres')
        grid.dim2.lon = double(ncread(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/%s/CERES_EBAF_Ed4.1_Subset_%s.ymonmean.nc', type, par.(type).yr_span), 'lon'));
        grid.dim2.lat = double(ncread(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/%s/CERES_EBAF_Ed4.1_Subset_%s.ymonmean.nc', type, par.(type).yr_span), 'lat'));
        save(sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/grid.mat', type), 'grid')
    end

    % save grid
end
