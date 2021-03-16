function mmm_vh_mon(type, par)
    lat = par.lat;

    % output info
    foldername = make_savedir_proc(type, par);
    load(sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s/%s/grid.mat', type, par.outname, par.clim));

    par.lat_interp = 'native'; % input files will be in native grid

    for l = {'lo'}; land=l{1};
        f_vec = par.gcm.fw;
        for f = f_vec; fw = f{1};
            vh_mon_list.(land).(fw) = nan(length(par.model_list), length(par.lat), 12);
            vh_mon_mmm.(land).(fw) = nan(length(par.lat), 12);
            vh_mon_std.(land).(fw) = nan(length(par.lat), 12);
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
        vh_mon0 = load(sprintf('%s/vh_mon.mat', prefix_proc));

        for l = {'lo'}; land=l{1};
            f_vec = par.(type).fw;
            for f = f_vec; fw = f{1};
                vh_mon0i.vh_mon.(land).(fw) = interp1(grid0.grid.dim3.lat, vh_mon0.vh_mon.(land).(fw), grid.dim3.lat);
                vh_mon_list.(land).(fw)(k,:,:) = vh_mon0i.vh_mon.(land).(fw);
            end
        end % land
    end % models

    for l = {'lo'}; land=l{1};
        f_vec = par.(type).fw;
        for f = f_vec; fw = f{1};
            vh_mon_mmm.(land).(fw) = squeeze(nanmean(vh_mon_list.(land).(fw), 1));
            if strcmp(type, 'rea')
                vh_mon_std.(land).(fw) = squeeze(range(vh_mon_list.(land).(fw), 1));
            else
                vh_mon_std.(land).(fw) = squeeze(nanstd(vh_mon_list.(land).(fw), 1));
            end
        end
    end

    vh_mon = vh_mon_mmm;

    % save energy flux data into mat file
    if ~exist(foldername, 'dir')
        mkdir(foldername)
    end
    save(sprintf('%s%s', foldername, 'vh_mon'), 'vh_mon', 'lat', '-v7.3');

end

 
