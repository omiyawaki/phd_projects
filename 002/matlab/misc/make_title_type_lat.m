function make_title_type_lat(type, lat1, lat2, par)

    if any(strcmp(type, {'era5', 'era5c', 'erai'})); title(sprintf('%s, $\\phi=%g^\\circ$ to $%g^\\circ$', upper(type), lat1, lat2));
    elseif any(strcmp(type, 'merra2')); title(sprintf('%s, $\\phi=%g^\\circ$ to $%g^\\circ$', upper(type), lat1, lat2));
    elseif any(strcmp(type, 'echam')); title(sprintf('%s, $\\phi=%g^\\circ$ to $%g^\\circ$', upper(type), lat1, lat2));
    elseif strcmp(type, 'gcm');
        if contains(par.model, 'mmm')
            title(sprintf('CMIP5 %s, $\\phi=%g^\\circ$ to $%g^\\circ$', par.gcm.clim, lat1, lat2));
        else
            title(sprintf('%s, $\\phi=%g^\\circ$ to $%g^\\circ$', par.model, lat1, lat2));
        end
    end;

end
