function mmm_ga_diff_mon_lat(type, par)

    % output info
    foldername = make_savedir_proc(type, par);
    load(sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s/%s/grid.mat', type, par.outname, par.clim));

    par.lat_interp = 'native';

    lat = par.lat;
    for l = par.land_list; land=l{1};
        ga_diff_list.(land) = nan(length(par.model_list), length(par.lat), 12, length(par.si));
        ga_diff_mmm.(land) = nan(length(par.lat), 12, length(par.si));
        ga_diff_std.(land) = nan(length(par.lat), 12, length(par.si));
        ga_diff_min.(land) = nan(length(par.lat), 12, length(par.si));
        ga_diff_max.(land) = nan(length(par.lat), 12, length(par.si));
    end

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
        ga_diff0 = load(sprintf('%s/ga_diff_mon_lat.mat', prefix_proc));

        for l = par.land_list; land=l{1};
        % for l = {'lo'}; land=l{1};
            % interpolate native grid data to standard grid
            ga_diff0i = interp1(grid0.grid.dim3.lat, ga_diff0.ga_diff.(land), grid.dim3.lat);
            ga_diff_list.(land)(k,:,:,:) = ga_diff0i;
        end % land

    end % models

    for l = par.land_list; land=l{1};
    % for l = {'lo'}; land=l{1};
        ga_diff_mmm.(land) = squeeze(nanmean(ga_diff_list.(land),1));
        ga_diff_std.(land) = squeeze(nanstd(ga_diff_list.(land),1));
        ga_diff_min.(land) = squeeze(min(ga_diff_list.(land), [], 1));
        ga_diff_max.(land) = squeeze(max(ga_diff_list.(land), [], 1));
    end

    ga_diff = ga_diff_mmm;

    printname = [foldername 'ga_diff_mon_lat'];
    if ~exist(foldername, 'dir')
        mkdir(foldername)
    end
    save(printname, 'ga_diff', 'ga_diff_std', 'ga_diff_min', 'ga_diff_max', 'lat', '-v7.3');

end
