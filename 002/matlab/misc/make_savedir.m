function savedir = make_savedir(type, par)
    if any(strcmp(type, {'era5', 'era5c', 'erai', 'merra2', 'jra55'}))
        savedir = sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s', type, par.(type).yr_span);
    elseif strcmp(type, 'gcm')
        savedir = sprintf('/project2/tas1/miyawaki/projects/002/data/read/gcm/%s/%s', par.model, par.gcm.clim);
    elseif strcmp(type, 'echam')
        savedir = sprintf('/project2/tas1/miyawaki/projects/002/data/read/echam/%s', par.echam.clim);
    end
end
