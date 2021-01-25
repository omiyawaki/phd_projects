function plotdir = make_plotdir(type, par)

    if any(strcmp(type, {'era5', 'erai', 'era5c', 'jra55', 'merra2'}))
        plotdir = sprintf('./figures/%s/%s/%s', type, par.(type).yr_span, par.lat_interp);
    elseif strcmp(type, 'gcm')
        plotdir = sprintf('./figures/%s/%s/%s/%s', type, par.model, par.gcm.clim, par.lat_interp);
    elseif any(strcmp(type, {'echam', 'hahn'}))
        plotdir = sprintf('./figures/%s/%s/%s', type, par.(type).clim, par.lat_interp);
    end

end
