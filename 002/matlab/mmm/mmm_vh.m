function mmm_vh(type, par)
    lat = par.lat;

    % output info
    foldername = make_savedir_proc(type, par);
    load(sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s/%s/grid.mat', type, par.outname, par.clim));

    par.lat_interp = 'native';

    for l = par.land_list; land=l{1};
        for t = {'ann', 'djf', 'mam', 'jja', 'son'}; time = t{1};
            f_vec = par.(type).fw;
            for f = f_vec; fw = f{1};
                vh_list.(land).(time).(fw) = nan(length(par.model_list),length(par.lat));
                vh_mmm.(land).(time).(fw) = nan(1,length(par.lat));
                vh_std.(land).(time).(fw) = nan(1,length(par.lat));
                vh_min.(land).(time).(fw) = nan(1,length(par.lat));
                vh_max.(land).(time).(fw) = nan(1,length(par.lat));
                vh_25.(land).(time).(fw) = nan(1,length(par.lat));
                vh_75.(land).(time).(fw) = nan(1,length(par.lat));
            end
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
        vh0 = load(sprintf('%s/vh.mat', prefix_proc));

        for l = par.land_list; land=l{1};
            for t = {'ann', 'djf', 'mam', 'jja', 'son'}; time = t{1};
                f_vec = par.(type).fw;
                for f = f_vec; fw = f{1};
                    vh0i.vh.(land).(time).(fw) = interp1(grid0.grid.dim3.lat, vh0.vh.(land).(time).(fw), grid.dim3.lat);
                    vh_list.(land).(time).(fw)(k,:) = vh0i.vh.(land).(time).(fw);
                end
            end % time
        end % land
    end % models

    for l = par.land_list; land=l{1};
        for t = {'ann', 'djf', 'mam', 'jja', 'son'}; time = t{1};
            f_vec = par.(type).fw;
            for f = f_vec; fw = f{1};
                vh_mmm.(land).(time).(fw) = squeeze(nanmean(vh_list.(land).(time).(fw), 1));
                vh_std.(land).(time).(fw) = squeeze(nanstd(vh_list.(land).(time).(fw), 1));
                vh_min.(land).(time).(fw) = squeeze(min(vh_list.(land).(time).(fw), [], 1));
                vh_max.(land).(time).(fw) = squeeze(max(vh_list.(land).(time).(fw), [], 1));
                vh_25.(land).(time).(fw) = squeeze(prctile(vh_list.(land).(time).(fw), [25], 1));
                vh_75.(land).(time).(fw) = squeeze(prctile(vh_list.(land).(time).(fw), [75], 1));
            end
        end
    end

    vh = vh_mmm;

    % save energy flux data into mat file
    if ~exist(foldername, 'dir')
        mkdir(foldername)
    end
    save(sprintf('%s%s', foldername, 'vh'), 'vh', 'vh_std', 'vh_min', 'vh_max', 'vh_25', 'vh_75', 'lat', '-v7.3');

end
 
