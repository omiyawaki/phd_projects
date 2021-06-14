function mmm_ga_frac_midlatitude(type, par)
    lat = par.lat;

    % output info
    foldername = make_savedir_proc(type, par);
    load(sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s/%s/grid.mat', type, par.outname, par.clim));

    par.lat_interp = 'native'; % input files will be in native grid

    par.lat_bound_list = [-80 80];

    for lb = 1:length(par.lat_bound_list); par.lat_bound = par.lat_bound_list(lb);

        [lat, clat, clat_mon, par] = make_polar_lat(par);

        filename = sprintf('ga_frac_poleward_of_lat_%g.mat', par.lat_bound);

        for l = par.land_list; land=l{1};
            ga_frac_lat_list.(land) = nan(length(par.model_list), 12, length(grid.dim3.si));
            ga_frac_lat_mmm.(land) = nan(1, 12, length(grid.dim3.si));
            ga_frac_lat_std.(land) = nan(1, 12, length(grid.dim3.si));
            ga_frac_lat_min.(land) = nan(1, 12, length(grid.dim3.si));
            ga_frac_lat_max.(land) = nan(1, 12, length(grid.dim3.si));
            ga_frac_lat_25.(land) = nan(1, 12, length(grid.dim3.si));
            ga_frac_lat_75.(land) = nan(1, 12, length(grid.dim3.si));
        end % land


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
            prefix_proc = make_prefix_proc(type_in, par);
            grid0 = load(sprintf('%s/grid.mat', prefix));
            ga_frac_lat_0 = load(sprintf('%s/%s', prefix_proc, filename));

            for l = par.land_list; land=l{1};
                ga_frac_lat_list.(land)(k,:,:) = ga_frac_lat_0.ga_frac_lat.(land);
            end % land
        end % models

    for l = par.land_list; land=l{1};
        ga_frac_lat_mmm.(land) = squeeze(nanmean(ga_frac_lat_list.(land), 1));
        ga_frac_lat_std.(land) = squeeze(nanstd(ga_frac_lat_list.(land), 1));
        ga_frac_lat_min.(land) = squeeze(min(ga_frac_lat_list.(land), [], 1));
        ga_frac_lat_max.(land) = squeeze(max(ga_frac_lat_list.(land), [], 1));
        ga_frac_lat_25.(land) = squeeze(prctile(ga_frac_lat_list.(land), [25], 1));
        ga_frac_lat_75.(land) = squeeze(prctile(ga_frac_lat_list.(land), [75], 1));
    end % land

    ga_frac_lat = ga_frac_lat_mmm;

    % save energy flux data into mat file
    if ~exist(foldername, 'dir')
        mkdir(foldername)
    end
    save(sprintf('%s%s', foldername, filename), 'ga_frac_lat', 'ga_frac_lat_std', 'ga_frac_lat_min', 'ga_frac_lat_max', 'ga_frac_lat_25', 'ga_frac_lat_75', 'lat', '-v7.3');

    end % lat bound

end
 
