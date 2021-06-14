function read_hus_ml0(type, ymonmean, par)

if strcmp(ymonmean, 'ymonmean')
    ymm_in = '.ymonmean';
    ymm_out = '';
elseif strcmp(ymonmean, 'mon')
    ymm_in = '';
    ymm_out = '_mon';
end

% land fraction
    if any(strcmp(type, {'era5', 'era5c', 'erai'}))
        hus_ml0 = double(ncread(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/%s/hus_ml0/%s_hus_ml0_%s%s.nc', type, type, par.(type).yr_span, ymm_in), 'si10'));
        save(sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s/hus_ml0%s.mat', type, par.(type).yr_span, ymm_out), 'hus_ml0');
    elseif strcmp(type, 'gcm')
        file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/gcm/%s/hus_ml0_*%s*.nc', par.model, par.gcm.clim));
        fullpath=sprintf('%s/%s', file.folder, file.name);
        hus_ml0=double(ncread(fullpath, 'hus_ml0'));
        newdir=sprintf('/project2/tas1/miyawaki/projects/002/data/read/gcm/%s/%s/%s', par.model, par.gcm.clim, par.gcm.yr_span);
        if ~exist(newdir, 'dir'); mkdir(newdir); end
        filename='hus_ml0.mat';
        save(sprintf('%s/%s', newdir, filename), 'hus_ml0');
    elseif any(strcmp(type, {'merra2', 'merra2c'}))
        hus_ml0 = double(ncread(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/%s/hus_ml0/%s_hus_ml0_%s%s.nc', type, type, par.(type).yr_span, ymm_in), 'hus_ml0'));
        % uWind = double(ncread(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/%s/hus_ml0/%s_hus_ml0_%s%s.nc', type, type, par.(type).yr_span, ymm_in), 'U10M'));
        % vWind = double(ncread(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/%s/hus_ml0/%s_hus_ml0_%s%s.nc', type, type, par.(type).yr_span, ymm_in), 'V10M'));
        % hus_ml0 = sqrt(uWind.^2 + vWind.^2);
        save(sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s/hus_ml0%s.mat', type, par.(type).yr_span, ymm_out), 'hus_ml0');
    elseif strcmp(type, 'jra55')
        hus_ml0 = double(ncread(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/%s/hus_ml0/%s_hus_ml0_1980_2005%s.nc', type, type, ymm_in), 'SPFH_GDS4_HYBL_S123'));
        hus_ml0_lat = double(ncread(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/%s/hus_ml0/%s_hus_ml0_1980_2005%s.nc', type, type, ymm_in), 'lat'));
        hus_ml0_lon = double(ncread(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/%s/hus_ml0/%s_hus_ml0_1980_2005%s.nc', type, type, ymm_in), 'lon'));
        save(sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s/hus_ml0%s.mat', type, par.(type).yr_span, ymm_out), 'hus_ml0', 'hus_ml0_lat', 'hus_ml0_lon');
    elseif strcmp(type, 'echam')
        % file=dir(sprintf('/project2/tas1/CMIP5_piControl/%s/sice_*.nc', 'MPI-ESM-LR'));
        % fullpath=sprintf('%s/%s', file.folder, file.name);
        % sice=double(ncread(fullpath, 'sice'));
        % newdir=sprintf('/project2/tas1/miyawaki/projects/002/data/read/echam/%s', par.echam.clim);
        % if ~exist(newdir, 'dir'); mkdir(newdir); end
        % filename='sice.mat';
        % save(sprintf('%s/%s', newdir, filename), 'sice');
        if contains(par.echam.clim, 'rp000')
            file=dir(sprintf('/project2/tas1/ockham/data11/tas/echam-aiv_rcc_6.1.00p1/%s/BOT_%s_0020_39.nc', par.echam.clim, par.echam.clim));
        else
            file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/echam/BOT*_%s_*.ymonmean.nc', par.echam.clim));
        end
        fullpath=sprintf('%s/%s', file.folder, file.name);
        sice=double(ncread(fullpath, 'friac'));
        newdir=sprintf('/project2/tas1/miyawaki/projects/002/data/read/echam/%s', par.echam.clim);
        if ~exist(newdir, 'dir'); mkdir(newdir); end
        filename='sice.mat';
        save(sprintf('%s/%s', newdir, filename), 'sice');
    end
end
