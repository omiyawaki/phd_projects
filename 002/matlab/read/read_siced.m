function read_siced(type, par)
    if strcmp(type, 'echam')
        if contains(par.echam.clim, 'rp000')
            file=dir(sprintf('/project2/tas1/ockham/data11/tas/echam-aiv_rcc_6.1.00p1/%s/BOT_%s_0020_39.nc', par.echam.clim, par.echam.clim));
        else
            file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/echam/BOT*_%s_*.ymonmean.nc', par.echam.clim));
        end
        fullpath=sprintf('%s/%s', file.folder, file.name);
        siced=double(ncread(fullpath, 'siced'));
        newdir=sprintf('/project2/tas1/miyawaki/projects/002/data/read/echam/%s', par.echam.clim);
        if ~exist(newdir, 'dir'); mkdir(newdir); end
        filename='siced.mat';
        save(sprintf('%s/%s', newdir, filename), 'siced');
    elseif strcmp(type, 'gcm')
        if any(strcmp(par.gcm.clim, {'piControl', 'abrupt4xCO2'}))
            file=dir(sprintf('/project2/tas1/CMIP5_piControl/%s/sit_*.nc', par.model));
        else
            file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/gcm/%s/sit_*%s*.nc', par.model, par.gcm.clim));
        end
        fullpath=sprintf('%s/%s', file.folder, file.name);
        siced=double(ncread(fullpath, 'sit'));
        newdir=sprintf('/project2/tas1/miyawaki/projects/002/data/read/gcm/%s/%s', par.model, par.gcm.clim);
        if ~exist(newdir, 'dir'); mkdir(newdir); end
        filename='siced.mat';
        save(sprintf('%s/%s', newdir, filename), 'siced');
    else
        error('This function only works for ECHAM and CMIP data.');
    end
end
