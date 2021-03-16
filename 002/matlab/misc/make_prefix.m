function prefix = make_prefix(type, par, clim)
    if nargin < 3
        if any(strcmp(type, {'gcm', 'echam', 'hahn'}))
            clim = par.(type).clim;
        end
    end

    if any(strcmp(type, {'rea', 'era5', 'era5c', 'erai', 'merra2', 'jra55', 'echam_ml', 'echam_pl'}))
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s', type, par.(type).yr_span);
    elseif strcmp(type, 'gcm')
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s/%s', type, par.model, clim);
    elseif any(strcmp(type, {'echam', 'hahn'}))
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s', type, clim);
    end
end
