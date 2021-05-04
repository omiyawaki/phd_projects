function mmm_ga_malr_diff_midlatitude_line(type, par)
    lat = par.lat;

    % output info
    par.model = 'mmm';
    foldername = make_savedir_si_bl(type, par);
    load(sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s/%s/grid.mat', type, par.outname, par.clim));

    par.lat_interp = 'native';

    lat_bound_list = [-10 10];
    lat_center = 50;

    for lb = 1:length(lat_bound_list); par.lat_bound = lat_bound_list(lb);

        [lat, clat, clat_mon, par] = make_midlatitude_lat(lat_center, par);

        for l = par.land_list; land=l{1};
            ga_frac_lat_list.(land) = nan(length(par.model_list), 12);
            ga_frac_lat_mmm.(land) = nan(12);
            ga_frac_lat_std.(land) = nan(12);
            ga_frac_lat_min.(land) = nan(12);
            ga_frac_lat_max.(land) = nan(12);
            ga_frac_lat_25.(land) = nan(12);
            ga_frac_lat_75.(land) = nan(12);
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
            ga_frac_lat0 = load(sprintf('%s/si_bl_%g/ga_malr_diff_midlatitude_lat_%g_to_%g_%g.mat', prefix_proc, par.si_bl, par.lat_center-par.lat_bound, par.lat_center+par.lat_bound, par.si_up));

            for l = par.land_list; land=l{1};
                ga_frac_lat_list.(land)(k,:) = ga_frac_lat0.ga_frac_lat.(land);
            end % land

        end % models

        for l = par.land_list; land=l{1};
            ga_frac_lat_mmm.(land) = squeeze(nanmean(ga_frac_lat_list.(land),1));
            ga_frac_lat_std.(land) = squeeze(nanstd(ga_frac_lat_list.(land),1));
            ga_frac_lat_min.(land) = squeeze(min(ga_frac_lat_list.(land), [], 1));
            ga_frac_lat_max.(land) = squeeze(max(ga_frac_lat_list.(land), [], 1));
            ga_frac_lat_25.(land) = squeeze(prctile(ga_frac_lat_list.(land), [25], 1));
            ga_frac_lat_75.(land) = squeeze(prctile(ga_frac_lat_list.(land), [75], 1));
        end

        ga_frac_lat = ga_frac_lat_mmm;

        % save energy flux data into mat file
        if ~exist(foldername, 'dir')
            mkdir(foldername)
        end
        save(sprintf('%sga_malr_diff_midlatitude_lat_%g_to_%g_%g.mat', foldername, par.lat_center-par.lat_bound, par.lat_center+par.lat_bound, par.si_up), 'ga_frac_lat_mmm', 'ga_frac_lat_std', 'ga_frac_lat_min', 'ga_frac_lat_max', 'ga_frac_lat_25', 'ga_frac_lat_75', 'lat', '-v7.3');

    end

end

 
