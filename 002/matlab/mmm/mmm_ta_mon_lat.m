function mmm_ta_mon_lat(type, par)

    lat = par.lat;
    for l = {'lo'}; land=l{1};
        tasi_mmm.(land) = nan(length(par.lat), 12, length(par.si));
    end

    pb = CmdLineProgressBar("Creating the multi-model mean...");
    for k=1:length(par.gcm_models); par.model = par.gcm_models{k};
        pb.print(k, length(par.gcm_models)); % output progress of moist adiabat calculonion

        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s/%s', type, par.model, par.gcm.clim);
        grid0=load(sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s/%s/grid.mat', type, par.model, par.gcm.clim));
        tasi0=load(sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s/%s/ta_mon_lat.mat', type, par.model, par.gcm.clim, 'native'));

        foldername = sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s/%s/', type, par.outname, par.gcm.clim, par.lat_interp);
        load(sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s/%s/grid.mat', type, par.outname, par.gcm.clim));

        for l = {'lo'}; land=l{1};
            % interpolate native grid data to standard grid
            tasi0i = interp1(grid0.grid.dim3.lat, tasi0.tasi.(land), grid.dim3.lat);
            tasi_mmm.(land) = nanmean(cat(4,tasi0i,tasi_mmm.(land)),4);
        end % land

    end % models

    tasi = tasi_mmm;

    printname = [foldername 'ta_mon_lat'];
    if ~exist(foldername, 'dir')
        mkdir(foldername)
    end
    save(printname, 'tasi', 'lat', '-v7.3');

end
