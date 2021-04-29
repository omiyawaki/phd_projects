function make_title_type_mon_lat(type, mon_str, lat, par)

    if any(strcmp(type, {'era5', 'erai', 'merra2', 'jra55'}));
        title(sprintf('%s, %s, $\\phi=%g^\\circ$', upper(type), mon_str, lat));
    elseif strcmp(type, 'era5c')
        title(sprintf('%s, %s, $\\phi=%g^\\circ$', upper('era5'), mon_str, lat));
    elseif strcmp(type, 'merra2c')
        title(sprintf('%s, %s, $\\phi=%g^\\circ$', upper('merra2'), mon_str, lat));
    elseif strcmp(type, 'rea');
        title(sprintf('Reanalysis mean, %s, $\\phi=%g^\\circ$', mon_str, lat));
    elseif strcmp(type, 'gcm');
        if contains(par.model, 'mmm')
            title(sprintf('CMIP5 %s, %s, $\\phi=%g^\\circ$', par.gcm.clim, mon_str, lat));
        else
            title(sprintf('%s, %s, $\\phi=%g^\\circ$', par.model, mon_str, lat));
        end
    elseif any(strcmp(type, {'echam'}));
        title(sprintf('%s, %s, %s, $\\phi=%g^\\circ$', upper(type), par.(type).(par.(type).clim), mon_str, lat));
    elseif any(strcmp(type, {'hahn'}));
        title(sprintf('%s, %s, %s, $\\phi=%g^\\circ$', upper('CESM'), par.(type).(par.(type).clim), mon_str, lat));
    end

end
