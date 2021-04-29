function prefix_proc = make_prefix_proc(type, par, clim)
    if nargin < 3
        if any(strcmp(type, {'gcm', 'echam', 'hahn'}))
            clim = par.(type).clim;
        end
    end

    if any(strcmp(type, {'rea', 'era5', 'era5c', 'erai', 'merra2', 'merra2c', 'jra55', 'echam_ml', 'echam_pl'}))
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s', type, par.(type).yr_span, par.lat_interp);
    elseif strcmp(type, 'gcm')
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s/%s/%s', type, par.model, clim, par.(type).yr_span, par.lat_interp);
    elseif any(strcmp(type, {'echam', 'hahn'}))
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s', type, clim, par.lat_interp);
    end
end
