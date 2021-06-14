function read_sfcWind(type, ymonmean, par)

if strcmp(ymonmean, 'ymonmean')
    ymm_in = '.ymonmean';
    ymm_out = '';
elseif strcmp(ymonmean, 'mon')
    ymm_in = '';
    ymm_out = '_mon';
end

% land fraction
    if any(strcmp(type, {'era5', 'era5c', 'erai'}))
        sfcWind = double(ncread(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/%s/sfcWind/%s_sfcWind_%s%s.nc', type, type, par.(type).yr_span, ymm_in), 'si10'));
        save(sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s/sfcWind%s.mat', type, par.(type).yr_span, ymm_out), 'sfcWind');
    elseif strcmp(type, 'gcm')
        file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/gcm/%s/sfcWind_*%s*.nc', par.model, par.gcm.clim));
        fullpath=sprintf('%s/%s', file.folder, file.name);
        sfcWind=double(ncread(fullpath, 'sfcWind'));
        newdir=sprintf('/project2/tas1/miyawaki/projects/002/data/read/gcm/%s/%s/%s', par.model, par.gcm.clim, par.gcm.yr_span);
        if ~exist(newdir, 'dir'); mkdir(newdir); end
        filename='sfcWind.mat';
        save(sprintf('%s/%s', newdir, filename), 'sfcWind');
    elseif any(strcmp(type, {'merra2', 'merra2c'}))
        sfcWind = double(ncread(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/%s/sfcWind/%s_sfcWind_%s%s.nc', type, type, par.(type).yr_span, ymm_in), 'sfcWind'));
        % uWind = double(ncread(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/%s/sfcWind/%s_sfcWind_%s%s.nc', type, type, par.(type).yr_span, ymm_in), 'U10M'));
        % vWind = double(ncread(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/%s/sfcWind/%s_sfcWind_%s%s.nc', type, type, par.(type).yr_span, ymm_in), 'V10M'));
        % sfcWind = sqrt(uWind.^2 + vWind.^2);
        save(sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s/sfcWind%s.mat', type, par.(type).yr_span, ymm_out), 'sfcWind');
    elseif strcmp(type, 'jra55')
        % sfcWind = double(ncread(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/%s/sfcWind/%s_sfcWind_%s%s.nc', type, type, par.(type).yr_span, ymm_in), 'sfcWind'));
        sfcWind = double(ncread(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/%s/sfcWind/%s_sfcWind_1979_2005%s.nc', type, type, ymm_in), 'sfcWind'));
        % uWind = double(ncread(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/%s/sfcWind/%s_ugrd_%s%s.nc', type, type, par.(type).yr_span, ymm_in), 'UGRD_GDS0_HTGL_S123'));
        % vWind = double(ncread(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/%s/sfcWind/%s_vgrd_%s%s.nc', type, type, par.(type).yr_span, ymm_in), 'VGRD_GDS0_HTGL_S123'));
        % sfcWind = sqrt(uWind.^2 + vWind.^2);
        save(sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s/sfcWind%s.mat', type, par.(type).yr_span, ymm_out), 'sfcWind');
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
