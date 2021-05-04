function make_title_type_lat(type, lat1, lat2, par)

    if any(strcmp(type, {'era5', 'erai', 'merra2', 'jra55'}));
        title(sprintf('%s, $\\phi=%g^\\circ$ to $%g^\\circ$', upper(type), lat1, lat2));
    elseif any(strcmp(type, {'era5c'}));
        title(sprintf('%s, $\\phi=%g^\\circ$ to $%g^\\circ$', upper('era5'), lat1, lat2));
    elseif any(strcmp(type, {'merra2c'}));
        title(sprintf('%s, $\\phi=%g^\\circ$ to $%g^\\circ$', upper('merra2'), lat1, lat2));
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
                    title(sprintf('CMIP5 %s, NH High Latitudes', par.(type).clim));
                else
                    title(sprintf('CMIP5 %s, NH Midlatitudes', par.(type).clim));
                end
            else
                if lat2==-90
                    title(sprintf('CMIP5 %s, SH High Latitudes', par.(type).clim));
                else
                    title(sprintf('CMIP5 %s, SH Midlatitudes', par.(type).clim));
                end
            end
        else
            title(sprintf('%s, $\\phi=%g^\\circ$ to $%g^\\circ$', par.model, lat1, lat2));
        end
    elseif any(strcmp(type, {'echam'}))
        % title(sprintf('%s, %s, $\\phi=%g^\\circ$ to $%g^\\circ$', upper(type), par.(type).(par.(type).clim), lat1, lat2));
        if any(strcmp(par.echam.clim, {'rp000133', 'rp000141', 'rp000034', 'rp000086', 'rp000172', 'rp000092'}))
            if lat2==90 | lat2==-90
                title(sprintf('AQUA, %s, No RAE', par.(type).(par.(type).clim)));
                % title(sprintf('%s, %s, No RAE', upper(type), par.(type).(par.(type).clim)));
            else
                title(sprintf('AQUA, %s, NH-ML-like', par.(type).(par.(type).clim)));
                % title(sprintf('%s, %s, NH-like', upper(type), par.(type).(par.(type).clim)));
            end
        elseif any(strcmp(par.echam.clim, {'rp000145', 'rp000131', 'rp000147', 'rp000135', 'rp000149', 'rp000046'}))
            if lat2==90 | lat2==-90
                title(sprintf('AQUA, %s, No RAE', par.(type).(par.(type).clim)));
                % title(sprintf('%s, %s, No RAE', upper(type), par.(type).(par.(type).clim)));
            else
                title(sprintf('AQUA, %s, SH-ML-like', par.(type).(par.(type).clim)));
                % title(sprintf('%s, %s, SH-like', upper(type), par.(type).(par.(type).clim)));
            end
        elseif any(strcmp(par.echam.clim, {'rp000126', 'rp000148', 'rp000134', 'rp000146', 'rp000130', 'rp000144'}))
            title(sprintf('AQUA, %s, NH-HL-like', par.(type).(par.(type).clim)));
            % title(sprintf('%s, %s, NH-like', upper(type), par.(type).(par.(type).clim)));
        end
    elseif any(strcmp(type, {'hahn'}))
        title(sprintf('%s, %s, $\\phi=%g^\\circ$ to $%g^\\circ$', upper('CESM'), par.(type).(par.(type).clim), lat1, lat2));
    end;

end
