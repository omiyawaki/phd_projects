function make_title_type_time(type, time, par)

    if any(strcmp(type, {'era5', 'erai', 'merra2', 'jra55'}));
        title(sprintf('%s, %s', upper(type), upper(time)));
    elseif strcmp(type, 'era5c')
        title(sprintf('%s, %s', upper('era5'), upper(time)));
    elseif strcmp(type, 'gcm');
        if contains(par.model, 'mmm')
            title(sprintf('CMIP5 %s, %s', par.gcm.clim, upper(time)));
        else
            title(sprintf('%s, %s', par.model, upper(time)));
        end
    elseif any(strcmp(type, {'echam', 'hahn'}));
        title(sprintf('%s, %s, %s', upper(type), par.(type).(par.(type).clim), upper(time)));
    end

end
