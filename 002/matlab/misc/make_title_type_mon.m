function make_title_type_mon(type, mon_str, par)

    if any(strcmp(type, {'era5', 'erai', 'merra2', 'jra55'}));
        title(sprintf('%s, %s', upper(type), mon_str));
    elseif strcmp(type, 'era5c')
        title(sprintf('%s, %s', upper('era5'), mon_str));
    elseif strcmp(type, 'merra2c')
        title(sprintf('%s, %s', upper('merra2'), mon_str));
    elseif strcmp(type, 'rea');
        title(sprintf('Reanalysis mean, %s', mon_str));
    elseif strcmp(type, 'gcm');
        if contains(par.model, 'mmm')
            title(sprintf('CMIP5 %s, %s', par.gcm.clim, mon_str));
        else
            title(sprintf('%s, %s', par.model, mon_str));
        end
    elseif any(strcmp(type, {'echam'}));
        title(sprintf('%s, %s, %s', upper(type), par.(type).(par.(type).clim), mon_str));
    elseif any(strcmp(type, {'hahn'}));
        title(sprintf('%s, %s, %s', upper('CESM'), par.(type).(par.(type).clim), mon_str));
    end

end
