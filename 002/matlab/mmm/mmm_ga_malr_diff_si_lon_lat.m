function mmm_ga_malr_diff_t(type, par)
    lat = par.lat;

    % output info
    foldername = make_savedir_proc(type, par);
    load(sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s/%s/grid.mat', type, par.outname, par.clim));

    par.lat_interp = 'native'; % input files will be in native grid

    for l = par.land_list; land=l{1};
        for t = {'ann', 'djf', 'mam', 'jja', 'son'}; time = t{1};
            ga_malr_diff_t_list.(land).(time) = nan(length(par.model_list), length(par.lon), length(par.lat));
            ga_malr_diff_t_mmm.(land).(time) = nan(length(par.lon), length(par.lat));
            ga_malr_diff_t_std.(land).(time) = nan(length(par.lon), length(par.lat));
        end
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
        ga_malr_diff_t0 = load(sprintf('%s/si_bl_%g/ga_malr_diff_si_lon_lat_%g.mat', prefix_proc, par.si_bl, par.si_up));

        for l = par.land_list; land=l{1};
            for t = {'ann', 'djf', 'mam', 'jja', 'son'}; time = t{1};
                ga_malr_diff_t0i.ga_malr_diff_t.(land).(time) = interp1(grid0.grid.dim3.lon, ga_malr_diff_t0.ga_malr_diff_t.(land).(time), grid.dim3.lon);
                ga_malr_diff_t0i.ga_malr_diff_t.(land).(time) = permute(ga_malr_diff_t0i.ga_malr_diff_t.(land).(time), [2 1]);
                ga_malr_diff_t0i.ga_malr_diff_t.(land).(time) = interp1(grid0.grid.dim3.lat, ga_malr_diff_t0i.ga_malr_diff_t.(land).(time), grid.dim3.lat);
                ga_malr_diff_t0i.ga_malr_diff_t.(land).(time) = permute(ga_malr_diff_t0i.ga_malr_diff_t.(land).(time), [2 1]);
                
                ga_malr_diff_t_list.(land).(time)(k,:,:) = ga_malr_diff_t0i.ga_malr_diff_t.(land).(time);
            end % time
        end % land
    end % models

    for l = par.land_list; land=l{1};
        for t = {'ann', 'djf', 'mam', 'jja', 'son'}; time = t{1};
                ga_malr_diff_t_mmm.(land).(time) = squeeze(nanmean(ga_malr_diff_t_list.(land).(time), 1));
                if strcmp(type, 'rea')
                    ga_malr_diff_t_std.(land).(time) = squeeze(range(ga_malr_diff_t_list.(land).(time), 1));
                else
                    ga_malr_diff_t_std.(land).(time) = squeeze(nanstd(ga_malr_diff_t_list.(land).(time), 1));
                end

        end
    end

    ga_malr_diff_t = ga_malr_diff_t_mmm;

    % save energy ga_malr_diff_t data into mat file
    if ~exist(foldername, 'dir')
        mkdir(foldername)
    end
    save(sprintf('%ssi_bl_%g/%s_%g.mat', foldername, par.si_bl, 'ga_malr_diff_t', par.si_up), 'ga_malr_diff_t', 'ga_malr_diff_t_std', 'lat', '-v7.3');

end
 
