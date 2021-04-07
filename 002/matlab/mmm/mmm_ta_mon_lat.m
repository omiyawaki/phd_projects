function mmm_ta_mon_lat(type, par)

    % output info
    foldername = make_savedir_proc(type, par);
    load(sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s/%s/grid.mat', type, par.outname, par.clim));

    par.lat_interp = 'native';

    lat = par.lat;
    for l = par.land_list; land=l{1};
        tasi_list.(land) = nan(length(par.model_list), length(par.lat), 12, length(par.si));
        tasi_mmm.(land) = nan(length(par.lat), 12, length(par.si));
        tasi_std.(land) = nan(length(par.lat), 12, length(par.si));
        tasi_min.(land) = nan(length(par.lat), 12, length(par.si));
        tasi_max.(land) = nan(length(par.lat), 12, length(par.si));
        tasi_25.(land) = nan(length(par.lat), 12, length(par.si));
        tasi_75.(land) = nan(length(par.lat), 12, length(par.si));
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
        tasi0 = load(sprintf('%s/ta_mon_lat.mat', prefix_proc));

        for l = par.land_list; land=l{1};
            % interpolate native grid data to standard grid
            tasi0i = interp1(grid0.grid.dim3.lat, tasi0.tasi.(land), grid.dim3.lat);
            tasi_list.(land)(k,:,:,:) = tasi0i;
        end % land

    end % models

    for l = par.land_list; land=l{1};
        tasi_mmm.(land) = squeeze(nanmean(tasi_list.(land),1));
        tasi_std.(land) = squeeze(nanstd(tasi_list.(land),1));
        tasi_min.(land) = squeeze(min(tasi_list.(land), [], 1));
        tasi_max.(land) = squeeze(max(tasi_list.(land), [], 1));
        tasi_25.(land) = squeeze(prctile(tasi_list.(land), [25], 1));
        tasi_75.(land) = squeeze(prctile(tasi_list.(land), [75], 1));
    end

    tasi = tasi_mmm;

    printname = [foldername 'ta_mon_lat'];
    if ~exist(foldername, 'dir')
        mkdir(foldername)
    end
    save(printname, 'tasi', 'tasi_std', 'tasi_min', 'tasi_max', 'tasi_25', 'tasi_75', 'lat', '-v7.3');

end
