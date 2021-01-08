function read_lfrac(type, par)
% land fraction
    if any(strcmp(type, {'era5', 'era5c', 'erai'}))
        sftlf = double(ncread(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/%s/lfrac/%s_lfrac_%s.ymonmean.nc', type, type, par.(type).yr_span), 'lsm'));
        save(sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s/sftlf.mat', type, par.(type).yr_span), 'sftlf');
    elseif strcmp(type, 'gcm')
        if any(strcmp(par.gcm.clim, {'piControl', 'abrupt4xCO2'}))
            file=dir(sprintf('/project2/tas1/CMIP5_piControl/%s/sftlf_*.nc', par.model));
        else
            file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/%s_raw/%s/sftlf_*.nc', par.gcm.clim, par.model));
        end
        fullpath=sprintf('%s/%s', file.folder, file.name);
        sftlf=double(ncread(fullpath, 'sftlf'));
        newdir=sprintf('/project2/tas1/miyawaki/projects/002/data/read/gcm/%s/%s', par.model, par.gcm.clim);
        if ~exist(newdir, 'dir'); mkdir(newdir); end
        filename='sftlf.mat';
        save(sprintf('%s/%s', newdir, filename), 'sftlf');
    elseif strcmp(type, 'echam')
        file=dir(sprintf('/project2/tas1/CMIP5_piControl/%s/sftlf_*.nc', 'MPI-ESM-LR'));
        fullpath=sprintf('%s/%s', file.folder, file.name);
        sftlf=double(ncread(fullpath, 'sftlf'));
        newdir=sprintf('/project2/tas1/miyawaki/projects/002/data/read/echam/%s', par.echam.clim);
        if ~exist(newdir, 'dir'); mkdir(newdir); end
        filename='sftlf.mat';
        save(sprintf('%s/%s', newdir, filename), 'sftlf');
    end
end
