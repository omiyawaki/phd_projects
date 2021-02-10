function mmm_ga_dalr_bl_diff_si_mon_lat(type, par)
    lat = par.lat;

    for l = {'lo'}; land=l{1};
        ga_dalr_bl_diff_mmm.(land) = nan(length(par.lat), 12);
    end

    pb = CmdLineProgressBar("Creating the multi-model mean...");
    for k=1:length(par.gcm_models); par.model = par.gcm_models{k};
        pb.print(k, length(par.gcm_models)); % output progress of moist adiabat calculonion

        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s/%s', type, par.model, par.gcm.clim);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s', type, par.model, par.gcm.clim);
        grid0 = load(sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s/%s/grid.mat', type, par.model, par.gcm.clim));
        ga_dalr_bl_diff0 = load(sprintf('%s/%s/si_bl_%g/ga_dalr_bl_diff_si_mon_lat.mat', prefix_proc, 'native', par.si_bl)); % load lat x mon RCAE data

        foldername = sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s/%s/si_bl_%g/', type, par.outname, par.gcm.clim, par.lat_interp, par.si_bl);
        load(sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s/%s/grid.mat', type, par.outname, par.gcm.clim));

        for l = {'lo'}; land=l{1};
            ga_dalr_bl_diff0i.ga_dalr_bl_diff.(land) = interp1(grid0.grid.dim3.lat, ga_dalr_bl_diff0.ga_dalr_bl_diff.(land), grid.dim3.lat);
            ga_dalr_bl_diff_mmm.(land) = nanmean(cat(3, ga_dalr_bl_diff0i.ga_dalr_bl_diff.(land), ga_dalr_bl_diff_mmm.(land)), 3);
        end % land
    end % models

    ga_dalr_bl_diff = ga_dalr_bl_diff_mmm;

    % save energy flux data into mat file
    if ~exist(foldername, 'dir')
        mkdir(foldername)
    end
    save(sprintf('%s%s', foldername, 'ga_dalr_bl_diff_si_mon_lat'), 'ga_dalr_bl_diff', 'lat', '-v7.3');

end

 
