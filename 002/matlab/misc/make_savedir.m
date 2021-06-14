function savedir = make_savedir(type, par)
    if any(strcmp(type, {'era5', 'era5c', 'erai', 'merra2', 'merra2c', 'jra55', 'rea'}))
        savedir = sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s', type, par.(type).yr_span);
    elseif strcmp(type, 'gcm')
        if ~isfield(par, 'model')
            par.model = 'mmm';
        end
        savedir = sprintf('/project2/tas1/miyawaki/projects/002/data/read/gcm/%s/%s/%s', par.model, par.(type).clim, par.(type).yr_span);
    elseif any(strcmp(type, {'echam', 'hahn'}))
        savedir = sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s', type, par.(type).clim);
    end

    if ~exist(savedir, 'dir'); mkdir(savedir); end
end
