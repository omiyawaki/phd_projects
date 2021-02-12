function mmm_dmse_midlatitude_line(type, par)
    lat = par.lat;

    par.lat_bound_list = [-10 10];
    lat_center = 50;

    for lb = 1:length(par.lat_bound_list); par.lat_bound = par.lat_bound_list(lb);

        [lat, clat, clat_mon, par] = make_midlatitude_lat(lat_center, par);

        filename = sprintf('dmse_midlatitude_lat_%g_to_%g.mat', par.lat_center-par.lat_bound, par.lat_center+par.lat_bound);

        for l = {'lo'}; land=l{1};
            for fn = {'ra_lat', 'res_lat', 'lh_lat', 'sh_lat'}; fname = fn{1};
                f_vec = par.gcm.fw;
                for f = f_vec; fw = f{1};
                    dmse_list.(fname).(land).(fw) = nan(length(par.gcm_models), 12);
                end
            end % fnames
        end % land


        pb = CmdLineProgressBar("Creating the multi-model mean...");
        for k=1:length(par.gcm_models); par.model = par.gcm_models{k};
            pb.print(k, length(par.gcm_models)); % output progress of moist adiabat calculonion

            prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s/%s', type, par.model, par.gcm.clim);
            prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s', type, par.model, par.gcm.clim);
            grid0 = load(sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s/%s/grid.mat', type, par.model, par.gcm.clim));
            dmse_0 = load(sprintf('%s/%s/%s', prefix_proc, 'native', filename));

            foldername = sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s/%s/', type, par.outname, par.gcm.clim, par.lat_interp);
            load(sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s/%s/grid.mat', type, par.outname, par.gcm.clim));

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
                dmse_mmm.(fname).(land).(fw) = nanmean(dmse_list.(fname).(land).(fw), 1);
                dmse_std.(fname).(land).(fw) = nanstd(dmse_list.(fname).(land).(fw), 1);
            end
        end % fnames
    end % land

    dmse = dmse_mmm;

    % save energy flux data into mat file
    if ~exist(foldername, 'dir')
        mkdir(foldername)
    end
    save(sprintf('%s%s', foldername, filename), 'dmse', 'dmse_std', 'lat', '-v7.3');

    end % lat bound

end
 
