function savedir = make_savedir_si_bl(type, par)
    if any(strcmp(type, {'rea', 'era5', 'era5c', 'erai', 'merra2', 'merra2c', 'jra55'}))
        savedir = sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s/si_bl_%g/', type, par.(type).yr_span, par.lat_interp, par.si_bl);
    elseif strcmp(type, 'gcm')
        savedir = sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s/%s/%s/si_bl_%g/', type, par.model, par.(type).clim, par.(type).yr_span, par.lat_interp, par.si_bl);
    elseif any(strcmp(type, {'echam', 'hahn'}))
        savedir = sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s/si_bl_%g/', type, par.(type).clim, par.lat_interp, par.si_bl);
    end

    if ~exist(savedir, 'dir'); mkdir(savedir); end
end
