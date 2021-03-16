function make_title_type_lat(type, lat, par)

    if any(strcmp(type, {'era5', 'erai', 'merra2', 'jra55'}));
        title(sprintf('%s, $\\phi=%g^\\circ$', upper(type), lat));
    elseif any(strcmp(type, {'era5c'}));
        title(sprintf('%s, $\\phi=%g^\\circ$', upper('era5'), lat));
    elseif strcmp(type, 'rea');
        title(sprintf('Reanalysis mean, $\\phi=%g^\\circ$', lat));
    elseif strcmp(type, 'gcm');
        if contains(par.model, 'mmm')
            title(sprintf('CMIP5 %s, $\\phi=%g^\\circ$', par.gcm.clim, lat));
        else
            title(sprintf('%s, $\\phi=%g^\\circ$', par.model, lat));
        end
    elseif any(strcmp(type, {'echam'}))
        title(sprintf('%s, %s, $\\phi=%g^\\circ$', upper(type), par.(type).(par.(type).clim), lat));
    elseif any(strcmp(type, {'hahn'}))
        title(sprintf('%s, %s, $\\phi=%g^\\circ$', upper('CESM'), par.(type).(par.(type).clim), lat));
    end;

end
