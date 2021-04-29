function save_mask(type, par)
    if any(strcmp(type, {'era5', 'erai', 'era5c', 'jra55', 'merra2'}))
        foldername = sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s/', type, par.(type).yr_span, par.lat_interp);
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s', type, par.(type).yr_span);
    elseif strcmp(type, 'gcm')
        foldername = sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s/%s/', type, par.model, par.gcm.clim, par.lat_interp);
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s/%s', type, par.model, par.gcm.clim);
    elseif strcmp(type, 'echam')
        foldername = sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s/', type, par.echam.clim, par.lat_interp);
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s', type, par.echam.clim);
    elseif strcmp(type, 'echam_ml') | strcmp(type, 'echam_pl')
        foldername = sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s/', type, par.(type).yr_span, par.lat_interp);
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s', type, par.(type).yr_span);
    end

    load(sprintf('%s/grid.mat', prefix)); % read grid data
    if strcmp(par.lat_interp, 'std')
        lat = par.lat_std;
    else
        lat = grid.dim3.lat;
    end

    if strcmp(type, 'echam')
        if any(contains(par.echam.clim, {'rp000'}))
            ;
        else
            load(sprintf('%s/sftlf.mat', prefix)); % load land fraction data
        end
    else
        load(sprintf('%s/sftlf.mat', prefix)); % load land fraction data
    end

    if any(strcmp(type, {'era5c', 'erai', 'era5'}))
        mask.ocean = nan(size(sftlf)); mask.ocean(sftlf>0.5) = 1;
        mask.land = nan(size(mask.ocean)); mask.land(isnan(mask.ocean))=1;
    elseif strcmp(type, 'echam') & any(contains(par.echam.clim, {'rp000'}))
        mask.ocean = zeros([length(grid.dim2.lon) length(grid.dim2.lat)]);
        mask.land = ones([length(grid.dim2.lon) length(grid.dim2.lat)]);
    elseif any(strcmp(type, {'gcm', 'merra2', 'jra55', 'echam'})) % repeat to monthly dimension
        mask.ocean = nan(size(sftlf)); mask.ocean(sftlf>0.5) = 1; mask.ocean=repmat(mask.ocean,[1 1 12]);
        mask.land = nan(size(mask.ocean)); mask.land(isnan(mask.ocean))=1;
    else
        mask.land = remove_land(lat, grid.dim3.lon, 12);
        mask.ocean = remove_ocean(lat, grid.dim3.lon, 12);
    end

    % save masks
    printname = [foldername 'masks'];
    if ~exist(foldername, 'dir')
        mkdir(foldername)
    end
    save(printname, 'mask');

end % compute masks once and save it for later use
