function mmm_vh_mon(type, par)
    lat = par.lat;

    for l = {'lo'}; land=l{1};
        f_vec = par.gcm.fw;
        for f = f_vec; fw = f{1};
            vh_mon_mmm.(land).(fw) = nan(length(par.lat), 12);
        end
    end

    pb = CmdLineProgressBar("Creating the multi-model mean...");
    for k=1:length(par.gcm_models); par.model = par.gcm_models{k};
        pb.print(k, length(par.gcm_models)); % output progress of moist adiabat calculonion

        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s/%s', type, par.model, par.gcm.clim);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s', type, par.model, par.gcm.clim);
        grid0 = load(sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s/%s/grid.mat', type, par.model, par.gcm.clim));
        vh_mon0 = load(sprintf('%s/%s/vh_mon.mat', prefix_proc, 'native')); % load lat x mon RCAE data

        foldername = sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s/%s/', type, par.outname, par.gcm.clim, par.lat_interp);
        load(sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s/%s/grid.mat', type, par.outname, par.gcm.clim));

        for l = {'lo'}; land=l{1};
            f_vec = par.gcm.fw;
            for f = f_vec; fw = f{1};
                vh_mon0i.vh_mon.(land).(fw) = interp1(grid0.grid.dim3.lat, vh_mon0.vh_mon.(land).(fw), grid.dim3.lat);
                vh_mon_mmm.(land).(fw) = nanmean(cat(3, vh_mon0i.vh_mon.(land).(fw), vh_mon_mmm.(land).(fw)), 3);
            end
        end % land
    end % models

    vh_mon = vh_mon_mmm;

    % save energy flux data into mat file
    if ~exist(foldername, 'dir')
        mkdir(foldername)
    end
    save(sprintf('%s%s', foldername, 'vh_mon'), 'vh_mon', 'lat', '-v7.3');

end

 
