function mmm_ga_dalr_bl_diff_si_lat(type, par)
    lat = par.lat;

    % output info
    par.model = 'mmm';
    foldername = make_savedir_si_bl(type, par);
    load(sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s/%s/grid.mat', type, par.outname, par.clim));

    par.lat_interp = 'native'; % input files will be in native grid

    for l = par.land_list; land=l{1};
        for t = {'ann', 'djf', 'mam', 'jja', 'son'}; time = t{1};
            ga_dalr_bl_diff_zt_list.(land).(time) = nan(length(par.model_list), length(par.lat));
            ga_dalr_bl_diff_zt_mmm.(land).(time) = nan(1, length(par.lat));
            ga_dalr_bl_diff_zt_std.(land).(time) = nan(1, length(par.lat));
            ga_dalr_bl_diff_zt_min.(land).(time) = nan(1, length(par.lat));
            ga_dalr_bl_diff_zt_max.(land).(time) = nan(1, length(par.lat));
            ga_dalr_bl_diff_zt_25.(land).(time) = nan(1, length(par.lat));
            ga_dalr_bl_diff_zt_75.(land).(time) = nan(1, length(par.lat));
        end
    end

    pb = CmdLineProgressBar("Creating the multi-model mean...");
    for k=1:length(par.model_list); par.model = par.model_list{k};
        pb.print(k, length(par.model_list)); % output progress of model cycle

        if strcmp(type, 'gcm')
            type_in = type;
        else
            type_in = par.model;
        end

        % input info
        prefix = make_prefix(type_in, par);
        prefix_proc = make_prefix_proc(type_in, par);
        grid0 = load(sprintf('%s/grid.mat', prefix));
        ga_dalr_bl_diff_zt0 = load(sprintf('%s/si_bl_%g/ga_dalr_bl_diff_si_lat.mat', prefix_proc, par.si_bl));


        for l = par.land_list; land=l{1};
            for t = {'ann', 'djf', 'mam', 'jja', 'son'}; time = t{1};
                ga_dalr_bl_diff_zt0i.ga_dalr_bl_diff_zt.(land).(time) = interp1(grid0.grid.dim3.lat, ga_dalr_bl_diff_zt0.ga_dalr_bl_diff_zt.(land).(time), grid.dim3.lat);
                ga_dalr_bl_diff_zt_list.(land).(time)(k,:) = ga_dalr_bl_diff_zt0i.ga_dalr_bl_diff_zt.(land).(time);
            end % time
        end % land
    end % models

    for l = par.land_list; land=l{1};
        for t = {'ann', 'djf', 'mam', 'jja', 'son'}; time = t{1};
            ga_dalr_bl_diff_zt_mmm.(land).(time) = squeeze(nanmean(ga_dalr_bl_diff_zt_list.(land).(time), 1));
            ga_dalr_bl_diff_zt_std.(land).(time) = squeeze(nanstd(ga_dalr_bl_diff_zt_list.(land).(time), 1));
            ga_dalr_bl_diff_zt_min.(land).(time) = squeeze(min(ga_dalr_bl_diff_zt_list.(land).(time), [], 1));
            ga_dalr_bl_diff_zt_max.(land).(time) = squeeze(max(ga_dalr_bl_diff_zt_list.(land).(time), [], 1));
            ga_dalr_bl_diff_zt_25.(land).(time) = squeeze(prctile(ga_dalr_bl_diff_zt_list.(land).(time), [25], 1));
            ga_dalr_bl_diff_zt_75.(land).(time) = squeeze(prctile(ga_dalr_bl_diff_zt_list.(land).(time), [75], 1));
        end % time
    end % land

    ga_dalr_bl_diff_zt = ga_dalr_bl_diff_zt_mmm;

    % save energy ga_dalr_bl_diff data into mat file
    if ~exist(foldername, 'dir')
        mkdir(foldername)
    end
    save(sprintf('%sga_dalr_bl_diff_si_lat.mat', foldername), 'ga_dalr_bl_diff_zt', 'ga_dalr_bl_diff_zt_std', 'ga_dalr_bl_diff_zt_min', 'ga_dalr_bl_diff_zt_max', 'ga_dalr_bl_diff_zt_25' ,'ga_dalr_bl_diff_zt_75', 'lat', '-v7.3');

end

 
