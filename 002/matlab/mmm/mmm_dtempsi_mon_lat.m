function mmm_dtempsi_mon_lat(type, par)

    % output info
    foldername = make_savedir_proc(type, par);
    load(sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s/%s/%s/grid.mat', type, par.outname, par.clim, par.(type).yr_span));

    par.lat_interp = 'native';

    lat = par.lat;
    for l = par.land_list; land=l{1};
        dtempsi_list.(land) = nan(length(par.model_list), length(par.lat), 12, length(par.si));
        dtempsi_mmm.(land) = nan(length(par.lat), 12, length(par.si));
        dtempsi_std.(land) = nan(length(par.lat), 12, length(par.si));
        dtempsi_min.(land) = nan(length(par.lat), 12, length(par.si));
        dtempsi_max.(land) = nan(length(par.lat), 12, length(par.si));
        dtempsi_25.(land) = nan(length(par.lat), 12, length(par.si));
        dtempsi_75.(land) = nan(length(par.lat), 12, length(par.si));
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
        dtempsi0 = load(sprintf('%s/dtempsi_mon_lat.mat', prefix_proc));

        for l = par.land_list; land=l{1};
            % interpolate native grid data to standard grid
            dtempsi0i = interp1(grid0.grid.dim3.lat, dtempsi0.dtempsi.(land), grid.dim3.lat);
            % dtempsi0i = permute(dtempsi0i, [3,2,1]); % bring plev to front
            % dtempsi0ii = interp1(1e-5*grid0.grid.dim3.plev, dtempsi0i, par.si);
            % dtempsi0ii = permute(dtempsi0ii, [3,2,1]); % bring lat back to front
            dtempsi_list.(land)(k,:,:,:) = dtempsi0i;
        end % land

    end % models

    for l = par.land_list; land=l{1};
        dtempsi_mmm.(land) = squeeze(nanmean(dtempsi_list.(land),1));
        dtempsi_std.(land) = squeeze(nanstd(dtempsi_list.(land),1));
        dtempsi_min.(land) = squeeze(min(dtempsi_list.(land), [], 1));
        dtempsi_max.(land) = squeeze(max(dtempsi_list.(land), [], 1));
        dtempsi_25.(land) = squeeze(prctile(dtempsi_list.(land), [25], 1));
        dtempsi_75.(land) = squeeze(prctile(dtempsi_list.(land), [75], 1));
    end

    dtempsi = dtempsi_mmm;

    printname = [foldername 'dtempsi_mon_lat'];
    if ~exist(foldername, 'dir')
        mkdir(foldername)
    end
    save(printname, 'dtempsi', 'dtempsi_std', 'dtempsi_min', 'dtempsi_max', 'dtempsi_25', 'dtempsi_75', 'lat', '-v7.3');

end
