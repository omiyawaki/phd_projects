function mmm_dr1_midlatitude_line(type, par)
    lat = par.lat;

    par.lat_bound_list = [-10 10];
    center = 50;

    for lb = 1:length(par.lat_bound_list); par.lat_bound = par.lat_bound_list(lb);

        dlat = 0.25; % step size for standard lat grid
        if par.lat_bound>0; par.lat_center=center; lat = [-par.lat_bound:dlat:par.lat_bound]+par.lat_center; par.shiftby=0;
        else; par.lat_center=-center; lat = [-par.lat_bound:-dlat:par.lat_bound]+par.lat_center; par.shiftby=6; end;
        clat = cosd(lat); % cosine of latitude for cosine weighting
        clat_mon = repmat(clat', [1 12]);

        filename = sprintf('dr1_midlatitude_lat_%g_to_%g.mat', par.lat_center-par.lat_bound, par.lat_center+par.lat_bound);

        for l = {'lo'}; land=l{1};
            for fn = {'r1z_lat', 'r1z_ann_lat', 'dr1z_lat', 'comp1s_lat', 'comp2s_lat'}; fname = fn{1};
                f_vec = par.gcm.fw;
                for f = f_vec; fw = f{1};
                    dr1_list.(fname).(land).(fw) = nan(length(par.gcm_models), 12);
                end
            end % fnames
        end % land


        pb = CmdLineProgressBar("Creating the multi-model mean...");
        for k=1:length(par.gcm_models); par.model = par.gcm_models{k};
            pb.print(k, length(par.gcm_models)); % output progress of moist adiabat calculonion

            prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s/%s', type, par.model, par.gcm.clim);
            prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s', type, par.model, par.gcm.clim);
            grid0 = load(sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s/%s/grid.mat', type, par.model, par.gcm.clim));
            dr1_0 = load(sprintf('%s/%s/%s', prefix_proc, 'native', filename));

            foldername = sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s/%s/', type, par.outname, par.gcm.clim, par.lat_interp);
            load(sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s/%s/grid.mat', type, par.outname, par.gcm.clim));

            for l = {'lo'}; land=l{1};
                for fn = {'r1z_lat', 'r1z_ann_lat', 'dr1z_lat', 'comp1s_lat', 'comp2s_lat'}; fname = fn{1};
                    f_vec = par.gcm.fw;
                    for f = f_vec; fw = f{1};
                        dr1_list.(fname).(land).(fw)(k,:) = dr1_0.dr1.(fname).(land).(fw);
                    end
                end % fnames
            end % land
        end % models

    for l = {'lo'}; land=l{1};
        for fn = {'r1z_lat', 'r1z_ann_lat', 'dr1z_lat', 'comp1s_lat', 'comp2s_lat'}; fname = fn{1};
            f_vec = par.gcm.fw;
            for f = f_vec; fw = f{1};
                dr1_mmm.(fname).(land).(fw) = nanmean(dr1_list.(fname).(land).(fw), 1);
                dr1_std.(fname).(land).(fw) = nanstd(dr1_list.(fname).(land).(fw), 1);
            end
        end % fnames
    end % land

    dr1 = dr1_mmm;

    % save energy flux data into mat file
    if ~exist(foldername, 'dir')
        mkdir(foldername)
    end
    save(sprintf('%s%s', foldername, filename), 'dr1', 'dr1_std', 'lat', '-v7.3');

    end % lat bound

end
 
