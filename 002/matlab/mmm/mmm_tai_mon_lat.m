function mmm_tai_mon_lat(type, par)

    % output info
    foldername = make_savedir_proc(type, par);
    load(sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s/%s/grid.mat', type, par.outname, par.clim));

    par.lat_interp = 'native';

    lat = par.lat;
    for l = par.land_list; land=l{1};
        tai_list.(land) = nan(length(par.model_list), length(par.lat), 12, length(par.si));
        tai_mmm.(land) = nan(length(par.lat), 12, length(par.si));
        tai_std.(land) = nan(length(par.lat), 12, length(par.si));
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
        tai0 = load(sprintf('%s/tai_mon_lat.mat', prefix_proc));

        for l = par.land_list; land=l{1};
            % interpolate native grid data to standard grid
            tai0i = interp1(grid0.grid.dim3.lat, tai0.tai.(land), grid.dim3.lat);
            tai_list.(land)(k,:,:,:) = tai0i;
        end % land

    end % models

    for l = par.land_list; land=l{1};
        tai_mmm.(land) = squeeze(nanmean(tai_list.(land),1));
        if strcmp(type, 'rea')
            tai_std.(land) = squeeze(range(tai_list.(land),1));
        else
            tai_std.(land) = squeeze(nanstd(tai_list.(land),1));
        end
    end

    tai = tai_mmm;

    printname = [foldername 'tai_mon_lat'];
    if ~exist(foldername, 'dir')
        mkdir(foldername)
    end
    save(printname, 'tai', 'tai_std', 'lat', '-v7.3');

end
