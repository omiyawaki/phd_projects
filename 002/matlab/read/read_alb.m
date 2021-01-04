function read_alb(type, par)
% read sea ice depth
    if strcmp(type, 'echam')
        file=dir(sprintf('/project2/tas1/ockham/data11/tas/echam-aiv_rcc_6.1.00p1/%s/BOT_%s_0020_39.nc', par.echam.clim, par.echam.clim));
        fullpath=sprintf('%s/%s', file.folder, file.name);
        alb=double(ncread(fullpath, 'albedo'));
        newdir=sprintf('/project2/tas1/miyawaki/projects/002/data/read/echam/%s', par.echam.clim);
        if ~exist(newdir, 'dir'); mkdir(newdir); end
        filename='alb.mat';
        save(sprintf('%s/%s', newdir, filename), 'alb');
    else
        if any(strcmp(type, {'era5', 'era5c', 'erai'}))
            prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s', type);
            prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s', type);
        elseif strcmp(type, 'gcm')
            prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s/%s', type, par.model, par.gcm.clim);
            prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s', type, par.model, par.gcm.clim);
        end
        load(sprintf('%s/grid.mat', prefix)); % read grid data
        load(sprintf('%s/rad.mat', prefix)); % read radiation data

        if strcmp(type, 'erai') | strcmp(type, 'era5') | strcmp(type, 'era5c')
            alb = rad.rsus./rad.rsds; % TODO download shortwave radiation up and down separately for ERA
        elseif strcmp(type, 'gcm')
            alb = rad.rsus./rad.rsds;
        end

        if any(strcmp(type, {'era5', 'era5c', 'erai'})); newdir=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s', type);
        elseif strcmp(type, 'gcm'); newdir=sprintf('/project2/tas1/miyawaki/projects/002/data/read/gcm/%s/%s', par.model, par.gcm.clim); end;
        if ~exist(newdir, 'dir'); mkdir(newdir); end
        filename='alb.mat';
        save(sprintf('%s/%s', newdir, filename), 'alb', '-v7.3');

    end
end
