function mmm_vh(type, par)
    lat = par.lat;

    for l = {'lo'}; land=l{1};
        for t = {'ann', 'djf', 'mam', 'jja', 'son'}; time = t{1};
            f_vec = par.gcm.fw;
            for f = f_vec; fw = f{1};
                vh_mmm.(land).(time).(fw) = nan(1,length(par.lat));
            end
        end
    end

    pb = CmdLineProgressBar("Creating the multi-model mean...");
    for k=1:length(par.gcm_models); par.model = par.gcm_models{k};
        pb.print(k, length(par.gcm_models)); % output progress of moist adiabat calculonion

        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s/%s', type, par.model, par.gcm.clim);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s', type, par.model, par.gcm.clim);
        grid0 = load(sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s/%s/grid.mat', type, par.model, par.gcm.clim));
        vh0 = load(sprintf('%s/%s/vh.mat', prefix_proc, 'native')); % load lat x mon RCAE data

        foldername = sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s/%s/', type, par.outname, par.gcm.clim, par.lat_interp);
        load(sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s/%s/grid.mat', type, par.outname, par.gcm.clim));

        for l = {'lo'}; land=l{1};
            for t = {'ann', 'djf', 'mam', 'jja', 'son'}; time = t{1};
                f_vec = par.gcm.fw;
                for f = f_vec; fw = f{1};
                    vh0i.vh.(land).(time).(fw) = interp1(grid0.grid.dim3.lat, vh0.vh.(land).(time).(fw), grid.dim3.lat);
                    vh_mmm.(land).(time).(fw) = nanmean(cat(1, vh0i.vh.(land).(time).(fw), vh_mmm.(land).(time).(fw)), 1);
                end
            end % time
        end % land
    end % models

    vh = vh_mmm;

    % save energy flux data into mat file
    if ~exist(foldername, 'dir')
        mkdir(foldername)
    end
    save(sprintf('%s%s', foldername, 'vh'), 'vh', 'lat', '-v7.3');

end
 
