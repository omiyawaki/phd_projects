function mmm_dr2_midlatitude_line(type, par)
    lat = par.lat;

    % output info
    foldername = make_savedir_proc(type, par);
    load(sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s/%s/grid.mat', type, par.outname, par.clim));

    par.lat_interp = 'native'; % input files will be in native grid

    par.lat_bound_list = [-10 10];
    center = 50;

    for lb = 1:length(par.lat_bound_list); par.lat_bound = par.lat_bound_list(lb);

        dlat = 0.25; % step size for standard lat grid
        if par.lat_bound>0; par.lat_center=center; lat = [-par.lat_bound:dlat:par.lat_bound]+par.lat_center; par.shiftby=0;
        else; par.lat_center=-center; lat = [-par.lat_bound:-dlat:par.lat_bound]+par.lat_center; par.shiftby=6; end;
        clat = cosd(lat); % cosine of latitude for cosine weighting
        clat_mon = repmat(clat', [1 12]);

        filename = sprintf('dr1_midlatitude_lat_%g_to_%g.mat', par.lat_center-par.lat_bound, par.lat_center+par.lat_bound);
        savename = sprintf('dr2_midlatitude_lat_%g_to_%g.mat', par.lat_center-par.lat_bound, par.lat_center+par.lat_bound);

        for l = {'lo'}; land=l{1};
            for fn = {'r2z_lat', 'r2z_ann_lat', 'dr2z_lat', 'comp1_lat', 'comp2_lat'}; fname = fn{1};
                f_vec = par.gcm.fw;
                for f = f_vec; fw = f{1};
                    dr2_list.(fname).(land).(fw) = nan(length(par.model_list), 12);
                    dr2_mmm.(fname).(land).(fw) = nan(1, 12);
                    dr2_std.(fname).(land).(fw) = nan(1, 12);
                    dr2_min.(fname).(land).(fw) = nan(1, 12);
                    dr2_max.(fname).(land).(fw) = nan(1, 12);
                end
            end % fnames
        end % land


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
            dr2_0 = load(sprintf('%s/%s', prefix_proc, filename));

            for l = {'lo'}; land=l{1};
                for fn = {'r2z_lat', 'r2z_ann_lat', 'dr2z_lat', 'comp1_lat', 'comp2_lat'}; fname = fn{1};
                    f_vec = par.gcm.fw;
                    for f = f_vec; fw = f{1};
                        dr2_list.(fname).(land).(fw)(k,:) = dr2_0.dr2.(fname).(land).(fw);
                    end
                end % fnames
            end % land
        end % models

    for l = {'lo'}; land=l{1};
        for fn = {'r2z_lat', 'r2z_ann_lat', 'dr2z_lat', 'comp1_lat', 'comp2_lat'}; fname = fn{1};
            f_vec = par.gcm.fw;
            for f = f_vec; fw = f{1};
                dr2_mmm.(fname).(land).(fw) = squeeze(nanmean(dr2_list.(fname).(land).(fw), 1));
                dr2_std.(fname).(land).(fw) = squeeze(nanstd(dr2_list.(fname).(land).(fw), 1));
                dr2_min.(fname).(land).(fw) = squeeze(min(dr2_list.(fname).(land).(fw), [], 1));
                dr2_max.(fname).(land).(fw) = squeeze(max(dr2_list.(fname).(land).(fw), [], 1));
            end
        end % fnames
    end % land

    dr2 = dr2_mmm;

    % save energy flux data into mat file
    if ~exist(foldername, 'dir')
        mkdir(foldername)
    end
    save(sprintf('%s%s', foldername, savename), 'dr2', 'dr2_std', 'dr2_min', 'dr2_max', 'lat', '-v7.3');

    end % lat bound

end
 
