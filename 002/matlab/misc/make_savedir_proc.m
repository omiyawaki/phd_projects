function savedir = make_savedir_proc(type, par)
    if any(strcmp(type, {'era5', 'era5c', 'erai', 'merra2', 'jra55'}))
        savedir = sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s/', type, par.(type).yr_span, par.lat_interp);
    elseif strcmp(type, 'gcm')
        savedir = sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s/%s/', type, par.model, par.gcm.clim, par.lat_interp);
    elseif any(strcmp(type, {'echam', 'hahn'}))
        savedir = sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s/', type, par.(type).clim, par.lat_interp);
    end

    if ~exist(savedir, 'dir'); mkdir(savedir); end
end
