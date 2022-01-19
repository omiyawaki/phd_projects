function read_sice(type, par)
% land fraction
    if any(strcmp(type, {'era5', 'era5c', 'erai'}))
        sice = double(ncread(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/%s/sice/%s_sice_%s.ymonmean.nc', type, type, par.(type).yr_span), 'siconc'));
        save(sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s/sice.mat', type, par.(type).yr_span), 'sice');
    elseif strcmp(type, 'gcm')
        if any(strcmp(par.gcm.clim, {'piControl', 'abrupt4xCO2'}))
            file=dir(sprintf('/project2/tas1/CMIP5_piControl/%s/sic_*.nc', par.model));
        else
            file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/gcm/%s/sic_*%s*.nc', par.model, par.gcm.clim));
        end
        fullpath=sprintf('%s/%s', file.folder, file.name);
        sice=double(ncread(fullpath, 'sic'));
        newdir=sprintf('/project2/tas1/miyawaki/projects/002/data/read/gcm/%s/%s', par.model, par.gcm.clim);
        if ~exist(newdir, 'dir'); mkdir(newdir); end
        filename='sice.mat';
        save(sprintf('%s/%s', newdir, filename), 'sice');
    elseif any(strcmp(type, {'merra2', 'merra2c'}))
        sice = double(ncread(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/%s/sice/%s_sice_1980_2005.ymonmean.nc', type, type), 'FRSEAICE'));
        save(sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s/sice.mat', type, par.(type).yr_span), 'sice');
    elseif strcmp(type, 'jra55')
        sice = double(ncread(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/%s/sice/jra55_icec_1979_2005.ymonmean.nc', type), 'ICEC_GDS0_SFC_S113'));
        save(sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s/sice.mat', type, par.(type).yr_span), 'sice');
    elseif strcmp(type, 'echam')
        % file=dir(sprintf('/project2/tas1/CMIP5_piControl/%s/sice_*.nc', 'MPI-ESM-LR'));
        % fullpath=sprintf('%s/%s', file.folder, file.name);
        % sice=double(ncread(fullpath, 'sice'));
        % newdir=sprintf('/project2/tas1/miyawaki/projects/002/data/read/echam/%s', par.echam.clim);
        % if ~exist(newdir, 'dir'); mkdir(newdir); end
        % filename='sice.mat';
        % save(sprintf('%s/%s', newdir, filename), 'sice');
        if contains(par.echam.clim, 'echr000') | any(strcmp(par.echam.clim, par.echam.exceptions))
            file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/echam/BOT*_%s_*.ymonmean.nc', par.echam.clim));
        else
            file=dir(sprintf('/project2/tas1/ockham/data11/tas/echam-aiv_rcc_6.1.00p1/%s/BOT_%s_0020_39.nc', par.echam.clim, par.echam.clim));
        end
        fullpath=sprintf('%s/%s', file.folder, file.name);
        sice=double(ncread(fullpath, 'friac'));
        newdir=sprintf('/project2/tas1/miyawaki/projects/002/data/read/echam/%s', par.echam.clim);
        if ~exist(newdir, 'dir'); mkdir(newdir); end
        filename='sice.mat';
        save(sprintf('%s/%s', newdir, filename), 'sice');
    end
end
