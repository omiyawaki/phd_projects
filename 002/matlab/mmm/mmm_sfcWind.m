function mmm_sfcWind(type, par)
    lat = par.lat;

    % output info
    foldername = make_savedir(type, par);
    filename = 'sfcWind';
    load(sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s/%s/grid.mat', type, par.outname, par.clim));

    par.lat_interp = 'native'; % input files will be in native grid

    sfcWind_list = nan(length(par.model_list), length(grid.dim2.lon), length(grid.dim2.lat), 12);
    sfcWind_mmm = nan(length(grid.dim2.lon), length(grid.dim2.lat), 12);
    sfcWind_std = nan(length(grid.dim2.lon), length(grid.dim2.lat), 12);
    sfcWind_min = nan(length(grid.dim2.lon), length(grid.dim2.lat), 12);
    sfcWind_max = nan(length(grid.dim2.lon), length(grid.dim2.lat), 12);
    sfcWind_25 = nan(length(grid.dim2.lon), length(grid.dim2.lat), 12);
    sfcWind_75 = nan(length(grid.dim2.lon), length(grid.dim2.lat), 12);

    pb = CmdLineProgressBar("Creating the multi-model mean...");
    for k=1:length(par.model_list); par.model = par.model_list{k};
        pb.print(k, length(par.model_list)); % output progress of moist adiabat calculonion

        if strcmp(type, 'gcm')
            type_in = type;
        else
            type_in = par.model;
        end

        % input info
        prefix = make_prefix(type_in, par);
        grid0 = load(sprintf('%s/grid.mat', prefix));
        sfcWind0 = load(sprintf('%s/sfcWind.mat', prefix));

        tmp = sfcWind0.sfcWind;
        tmp = interp1(grid0.grid.dim2.lon, tmp, grid.dim2.lon);
        tmp = permute(tmp, [2 1 3]);
        if (size(grid0.grid.dim2.lat,1) ~= size(tmp,1)) & (size(grid0.grid.dim3.lat,1) ~= size(tmp,1))
            tmp = interp1(grid0.grid.dim2.lat_sfcWind, tmp, grid.dim2.lat);
        elseif size(grid0.grid.dim2.lat,1) ~= size(tmp,1)
            tmp = interp1(grid0.grid.dim3.lat, tmp, grid.dim2.lat);
        else
            tmp = interp1(grid0.grid.dim2.lat, tmp, grid.dim2.lat);
        end
        tmp = permute(tmp, [2 1 3]);
        sfcWind_list(k,:,:,:) = tmp;

    end % models

    sfcWind_mmm = squeeze(nanmean(sfcWind_list, 1));
    sfcWind_std = squeeze(nanstd(sfcWind_list, 1));
    sfcWind_min = squeeze(min(sfcWind_list, [], 1));
    sfcWind_max = squeeze(max(sfcWind_list, [], 1));
    sfcWind_25 = squeeze(prctile(sfcWind_list, [25], 1));
    sfcWind_75 = squeeze(prctile(sfcWind_list, [75], 1));

    sfcWind = sfcWind_mmm;

    % save data into mat file
    if ~exist(foldername, 'dir')
        mkdir(foldername)
    end
    save(sprintf('%s/%s', foldername, filename), 'sfcWind', 'sfcWind_std', 'sfcWind_min', 'sfcWind_max', 'sfcWind_25', 'sfcWind_75', 'lat', '-v7.3');

end
 
