function make_title_type_lat(type, lat1, lat2, par)

    if any(strcmp(type, {'era5', 'erai', 'merra2', 'jra55'}));
        title(sprintf('%s, $\\phi=%g^\\circ$ to $%g^\\circ$', upper(type), lat1, lat2));
    elseif any(strcmp(type, {'era5c'}));
        title(sprintf('%s, $\\phi=%g^\\circ$ to $%g^\\circ$', upper('era5'), lat1, lat2));
    elseif strcmp(type, 'gcm');
        if contains(par.model, 'mmm')
            title(sprintf('CMIP5 %s, $\\phi=%g^\\circ$ to $%g^\\circ$', par.gcm.clim, lat1, lat2));
        else
            title(sprintf('%s, $\\phi=%g^\\circ$ to $%g^\\circ$', par.model, lat1, lat2));
        end
    elseif any(strcmp(type, {'echam'}))
        title(sprintf('%s, %s, $\\phi=%g^\\circ$ to $%g^\\circ$', upper(type), par.(type).(par.(type).clim), lat1, lat2));
    elseif any(strcmp(type, {'hahn'}))
        title(sprintf('%s, %s, $\\phi=%g^\\circ$ to $%g^\\circ$', upper('CESM'), par.(type).(par.(type).clim), lat1, lat2));
    end;

end
