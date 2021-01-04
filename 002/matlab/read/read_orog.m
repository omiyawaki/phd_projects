function read_orog(type, par)
% orography
    if any(strcmp(type, {'era5', 'erai', 'era5c'}))
        % dimensions are (lon x lat x time)
        % albedo = double(ncread(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/%s/albedo/%s_albedo_%s.ymonmean.nc', type, type, par.(type).yr_span), 'fal'));
        % save(sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/albedo.mat', type), 'albedo');
    elseif strcmp(type, 'gcm')
        for i=1:length(par.gcm.vars.rad); var = par.gcm.vars.rad{i};
            prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s/%s', type, par.model, par.gcm.clim);
            file=dir(sprintf('/project2/tas1/CMIP5_%s/%s/orog_*.nc', par.gcm.clim, par.model));
            fullpath=sprintf('%s/%s', file.folder, file.name);
            orog=double(ncread(fullpath, 'orog'));
        end
        newdir=sprintf('/project2/tas1/miyawaki/projects/002/data/read/gcm/%s/%s', par.model, par.gcm.clim);
        if ~exist(newdir, 'dir'); mkdir(newdir); end
        filename='orog.mat';
        save(sprintf('%s/%s', newdir, filename), 'orog');
    elseif strcmp(type, 'echam')
        file=dir(sprintf('/project2/tas1/CMIP5_piControl/%s/orog_*.nc', 'MPI-ESM-LR'));
        fullpath=sprintf('%s/%s', file.folder, file.name);
        orog=double(ncread(fullpath, 'orog'));
        newdir=sprintf('/project2/tas1/miyawaki/projects/002/data/read/echam/%s', par.echam.clim);
        if ~exist(newdir, 'dir'); mkdir(newdir); end
        filename='orog.mat';
        save(sprintf('%s/%s', newdir, filename), 'orog');
    end
end
