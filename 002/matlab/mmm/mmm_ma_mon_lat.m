function mmm_ma_mon_lat(type, par)

    % output info
    foldername = make_savedir_proc(type, par);
    load(sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s/%s/grid.mat', type, par.outname, par.clim));

    par.lat_interp = 'native';

    lat = par.lat;
    for l = {'lo'}; land=l{1};
        masi_list.(land) = nan(length(par.model_list), length(par.lat), 12, length(par.si));
        masi_mmm.(land) = nan(length(par.lat), 12, length(par.si));
        masi_std.(land) = nan(length(par.lat), 12, length(par.si));
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
        masi0 = load(sprintf('%s/ma_mon_lat_%s.mat', prefix_proc, num2str(par.ma_init)));

        for l = {'lo'}; land=l{1};
            % interpolate native grid data to smandard grid
            masi0i = interp1(grid0.grid.dim3.lat, masi0.masi.(land), grid.dim3.lat);
            masi_list.(land)(k,:,:,:) = masi0i;
        end % land

    end % models

    for l = {'lo'}; land=l{1};
        masi_mmm.(land) = squeeze(nanmean(masi_list.(land),1));
        if strcmp(type, 'rea')
            masi_std.(land) = squeeze(range(masi_list.(land),1));
        else
            masi_std.(land) = squeeze(nanstd(masi_list.(land),1));
        end
    end

    masi = masi_mmm;

    printname = [foldername 'ma_mon_lat_' num2str(par.ma_init) '.mat'];
    if ~exist(foldername, 'dir')
        mkdir(foldername)
    end
    save(printname, 'masi', 'masi_std', 'lat', '-v7.3');

end

