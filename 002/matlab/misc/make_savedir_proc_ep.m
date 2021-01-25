function savedir = make_savedir_proc_ep(type, par)
    if any(strcmp(type, {'era5', 'era5c', 'erai', 'merra2', 'jra55'}))
        savedir = sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s/eps_%g_ga_%g/', type, par.(type).yr_span, par.lat_interp, par.ep, par.ga);
    elseif strcmp(type, 'gcm')
        foldername = sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s/%s/eps_%g_ga_%g/', type, par.model, par.gcm.clim, par.lat_interp, par.ep, par.ga);
    elseif any(strcmp(type, {'echam', 'hahn'}))
        foldername = sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s/eps_%g_ga_%g/', type, par.(type).clim, par.lat_interp, par.ep, par.ga);
    end

    if ~exist(savedir, 'dir'); mkdir(savedir); end
end
