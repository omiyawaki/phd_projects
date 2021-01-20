function make_title_type(type, par)

    if any(strcmp(type, {'era5', 'erai', 'merra2', 'jra55'}));
        title(sprintf('%s', upper(type)));
    elseif strcmp(type, 'era5c')
        title(sprintf('%s', upper('era5')));
    elseif strcmp(type, 'gcm');
        if contains(par.model, 'mmm')
            title(sprintf('CMIP5 %s', par.gcm.clim));
        else
            title(sprintf('%s', par.model));
        end
    elseif strcmp(type, 'echam'); title(sprintf('%s, %s', upper(type), par.echam.(par.echam.clim))); end;

end
