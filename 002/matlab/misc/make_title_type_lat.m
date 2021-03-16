function make_title_type_lat(type, lat1, lat2, par)

    if any(strcmp(type, {'era5', 'erai', 'merra2', 'jra55'}));
        title(sprintf('%s, $\\phi=%g^\\circ$ to $%g^\\circ$', upper(type), lat1, lat2));
    elseif any(strcmp(type, {'era5c'}));
        title(sprintf('%s, $\\phi=%g^\\circ$ to $%g^\\circ$', upper('era5'), lat1, lat2));
    elseif strcmp(type, 'rea');
        % title(sprintf('Reanalysis mean, $\\phi=%g^\\circ$ to $%g^\\circ$', lat1, lat2));
        if lat1>0
            if lat2==90
                title(sprintf('Reanalysis mean, NH High Latitudes'));
            else
                title(sprintf('Reanalysis mean, NH Midlatitudes'));
            end
        else
            if lat2==-90
                title(sprintf('Reanalysis mean, SH High Latitudes'));
            else
                title(sprintf('Reanalysis mean, SH Midlatitudes'));
            end
        end
    elseif strcmp(type, 'gcm');
        if contains(par.model, 'mmm')
            % title(sprintf('CMIP5 %s, $\\phi=%g^\\circ$ to $%g^\\circ$', par.gcm.clim, lat1, lat2));
            if lat1>0
                if lat2==90
                    title(sprintf('CMIP5 %s, NH High Latitudes'));
                else
                    title(sprintf('CMIP5 %s, NH Midlatitudes'));
                end
            else
                if lat2==-90
                    title(sprintf('CMIP5 %s, SH High Latitudes'));
                else
                    title(sprintf('CMIP5 %s, SH Midlatitudes'));
                end
            end
        else
            title(sprintf('%s, $\\phi=%g^\\circ$ to $%g^\\circ$', par.model, lat1, lat2));
        end
    elseif any(strcmp(type, {'echam'}))
        % title(sprintf('%s, %s, $\\phi=%g^\\circ$ to $%g^\\circ$', upper(type), par.(type).(par.(type).clim), lat1, lat2));
        if lat1>0
            if lat2==90
                title(sprintf('%s, %s, NH High Latitudes', upper(type), par.(type).(par.(type).clim)));
            else
                title(sprintf('%s, %s, NH Midlatitudes', upper(type), par.(type).(par.(type).clim)));
            end
        else
            if lat2==-90
                title(sprintf('%s, %s, SH High Latitudes', upper(type), par.(type).(par.(type).clim)));
            else
                title(sprintf('%s, %s, SH Midlatitudes', upper(type), par.(type).(par.(type).clim)));
            end
        end
    elseif any(strcmp(type, {'hahn'}))
        title(sprintf('%s, %s, $\\phi=%g^\\circ$ to $%g^\\circ$', upper('CESM'), par.(type).(par.(type).clim), lat1, lat2));
    end;

end
