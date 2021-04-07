function mmm_ga_malr_diff_si_mon_lat(type, par)
    lat = par.lat;

    % output info
    par.model = 'mmm';
    foldername = make_savedir_si_bl(type, par);
    load(sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s/%s/grid.mat', type, par.outname, par.clim));

    par.lat_interp = 'native';

    for l = par.land_list; land=l{1};
        ga_malr_diff_list.(land) = nan(length(par.model_list), length(par.lat), 12);
        ga_malr_diff_mmm.(land) = nan(length(par.lat), 12);
        ga_malr_diff_std.(land) = nan(length(par.lat), 12);
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
        ga_malr_diff0 = load(sprintf('%s/si_bl_%g/ga_malr_diff_si_mon_lat_%g.mat', prefix_proc, par.si_bl, par.si_up));


        for l = par.land_list; land=l{1};
            ga_malr_diff0i.ga_malr_diff.(land) = interp1(grid0.grid.dim3.lat, ga_malr_diff0.ga_malr_diff.(land), grid.dim3.lat);
            ga_malr_diff_list.(land)(k,:,:) = ga_malr_diff0i.ga_malr_diff.(land);
        end % land
    end % models

    for l = par.land_list; land=l{1};
        ga_malr_diff_mmm.(land) = squeeze(nanmean(ga_malr_diff_list.(land),1));
        ga_malr_diff_std.(land) = squeeze(nanstd(ga_malr_diff_list.(land),1));
    end

    ga_malr_diff = ga_malr_diff_mmm;

    % save energy flux data into mat file
    if ~exist(foldername, 'dir')
        mkdir(foldername)
    end
    save(sprintf('%sga_malr_diff_si_mon_lat_%g.mat', foldername, par.si_up), 'ga_malr_diff', 'lat', '-v7.3');

end

 
