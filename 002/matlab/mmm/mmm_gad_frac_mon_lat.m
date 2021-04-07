function mmm_gad_frac_mon_lat(type, par)

    % output info
    foldername = make_savedir_proc(type, par);
    load(sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s/%s/grid.mat', type, par.outname, par.clim));

    par.lat_interp = 'native';

    lat = par.lat;
    for l = par.land_list; land=l{1};
        gad_frac_list.(land) = nan(length(par.model_list), length(par.lat), 12, length(par.si));
        gad_frac_mmm.(land) = nan(length(par.lat), 12, length(par.si));
        gad_frac_std.(land) = nan(length(par.lat), 12, length(par.si));
        gad_frac_min.(land) = nan(length(par.lat), 12, length(par.si));
        gad_frac_max.(land) = nan(length(par.lat), 12, length(par.si));
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
        gad_frac0 = load(sprintf('%s/gad_frac_mon_lat.mat', prefix_proc));

        for l = par.land_list; land=l{1};
        % for l = {'lo'}; land=l{1};
            % interpolate native grid data to standard grid
            gad_frac0i = interp1(grid0.grid.dim3.lat, gad_frac0.gad_frac.(land), grid.dim3.lat);
            gad_frac_list.(land)(k,:,:,:) = gad_frac0i;
        end % land

    end % models

    for l = par.land_list; land=l{1};
    % for l = {'lo'}; land=l{1};
        gad_frac_mmm.(land) = squeeze(nanmean(gad_frac_list.(land),1));
        gad_frac_std.(land) = squeeze(nanstd(gad_frac_list.(land),1));
        gad_frac_min.(land) = squeeze(min(gad_frac_list.(land), [], 1));
        gad_frac_max.(land) = squeeze(max(gad_frac_list.(land), [], 1));
    end

    gad_frac = gad_frac_mmm;

    printname = [foldername 'gad_frac_mon_lat'];
    if ~exist(foldername, 'dir')
        mkdir(foldername)
    end
    save(printname, 'gad_frac', 'gad_frac_std', 'gad_frac_min', 'gad_frac_max', 'lat', '-v7.3');

end
