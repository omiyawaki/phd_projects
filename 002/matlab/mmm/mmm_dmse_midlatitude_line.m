function mmm_dmse_midlatitude_line(type, par)
    lat = par.lat;

    % output info
    foldername = make_savedir_proc(type, par);
    load(sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s/%s/grid.mat', type, par.outname, par.clim));

    par.lat_interp = 'native'; % input files will be in native grid

    par.lat_bound_list = [-10 10];
    lat_center = 50;

    for lb = 1:length(par.lat_bound_list); par.lat_bound = par.lat_bound_list(lb);

        [lat, clat, clat_mon, par] = make_midlatitude_lat(lat_center, par);

        filename = sprintf('dmse_midlatitude_lat_%g_to_%g.mat', par.lat_center-par.lat_bound, par.lat_center+par.lat_bound);

        for l = {'lo'}; land=l{1};
            for fn = {'ra_lat', 'res_lat', 'lh_lat', 'sh_lat'}; fname = fn{1};
                f_vec = par.gcm.fw;
                for f = f_vec; fw = f{1};
                    dmse_list.(fname).(land).(fw) = nan(length(par.model_list), 12);
                    dmse_mmm.(fname).(land).(fw) = nan(1, 12);
                    dmse_std.(fname).(land).(fw) = nan(1, 12);
                    dmse_min.(fname).(land).(fw) = nan(1, 12);
                    dmse_max.(fname).(land).(fw) = nan(1, 12);
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
            dmse_0 = load(sprintf('%s/%s', prefix_proc, filename));

            for l = {'lo'}; land=l{1};
                for fn = {'ra_lat', 'res_lat', 'lh_lat', 'sh_lat'}; fname = fn{1};
                    f_vec = par.gcm.fw;
                    for f = f_vec; fw = f{1};
                        dmse_list.(fname).(land).(fw)(k,:) = dmse_0.dmse.(fname).(land).(fw);
                    end
                end % fnames
            end % land
        end % models

    for l = {'lo'}; land=l{1};
        for fn = {'ra_lat', 'res_lat', 'lh_lat', 'sh_lat'}; fname = fn{1};
            f_vec = par.gcm.fw;
            for f = f_vec; fw = f{1};
                dmse_mmm.(fname).(land).(fw) = squeeze(nanmean(dmse_list.(fname).(land).(fw), 1));
                dmse_std.(fname).(land).(fw) = squeeze(nanstd(dmse_list.(fname).(land).(fw), 1));
                dmse_min.(fname).(land).(fw) = squeeze(min(dmse_list.(fname).(land).(fw), [], 1));
                dmse_max.(fname).(land).(fw) = squeeze(max(dmse_list.(fname).(land).(fw), [], 1));
            end
        end % fnames
    end % land

    dmse = dmse_mmm;

    % save energy flux data into mat file
    if ~exist(foldername, 'dir')
        mkdir(foldername)
    end
    save(sprintf('%s%s', foldername, filename), 'dmse', 'dmse_std', 'dmse_min', 'dmse_max', 'lat', '-v7.3');

    end % lat bound

end
 
