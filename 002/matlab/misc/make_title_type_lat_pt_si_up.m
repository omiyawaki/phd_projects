function make_title_type_lat_pt_si_up(type, lat, par)

    if any(strcmp(type, {'era5', 'erai', 'merra2', 'jra55'}));
        title(sprintf('%s, $\\phi=%g^\\circ$, $\\sigma_u=%g$', upper(type), lat, par.si_up));
    elseif any(strcmp(type, {'era5c'}));
        title(sprintf('%s, $\\phi=%g^\\circ$', upper('era5'), lat));
    elseif any(strcmp(type, {'merra2c'}));
        title(sprintf('%s, $\\phi=%g^\\circ$', upper('merra2'), lat));
    elseif strcmp(type, 'rea');
        title(sprintf('Reanalysis mean, $\\phi=%g^\\circ$, $\\sigma_u=%g$', lat, par.si_up));
    elseif strcmp(type, 'gcm');
        if contains(par.model, 'mmm')
            title(sprintf('CMIP5 %s, $\\phi=%g^\\circ$, $\\sigma_u=%g$', par.gcm.clim, lat, par.si_up));
        else
            title(sprintf('%s, $\\phi=%g^\\circ$, $\\sigma_u=%g$', par.model, lat, par.si_up));
        end
    elseif any(strcmp(type, {'echam'}))
        title(sprintf('AQUA, %s, $\\phi=%g^\\circ$, $\\sigma_u=%g$', par.(type).(par.(type).clim), lat, par.si_up));
        % title(sprintf('%s, %s, $\\phi=%g^\\circ$', upper(type), par.(type).(par.(type).clim), lat));
    elseif any(strcmp(type, {'hahn'}))
        title(sprintf('%s, %s, $\\phi=%g^\\circ$, $\\sigma_u=%g$', upper('CESM'), par.(type).(par.(type).clim), lat, par.si_up));
    end;

end
