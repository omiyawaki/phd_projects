function mmm_flux_z(type, par)
    lat = par.lat;

    for l = {'lo'}; land=l{1};
        var_vec = {'hfls', 'hfss', 'pr', 'evspsbl', 'lw', 'sw', 'rtoa', 'olr', 'lwsfc', 'swsfc'};
        for fn = var_vec; fname = fn{1};
            flux_z_mmm.(land).(fname) = nan(length(par.lat), 12);
        end
        for fn = {'ra', 'stf', 'res', 'r1', 'r2', 'ftoa', 'fsfc', 'sfc', 'shf', 'comp1', 'comp2'}; fname = fn{1};
            f_vec = par.gcm.fw;
            for f = f_vec; fw = f{1};
                flux_z_mmm.(land).(fname).(fw) = nan(length(par.lat), 12);
            end
        end
    end

    pb = CmdLineProgressBar("Creating the multi-model mean...");
    for k=1:length(par.gcm_models); par.model = par.gcm_models{k};
        pb.print(k, length(par.gcm_models)); % output progress of moist adiabat calculonion

        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s/%s', type, par.model, par.gcm.clim);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s', type, par.model, par.gcm.clim);
        grid0 = load(sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s/%s/grid.mat', type, par.model, par.gcm.clim));
        flux_z0 = load(sprintf('%s/%s/flux_z.mat', prefix_proc, 'native')); % load lat x mon RCAE data

        foldername = sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s/%s/', type, par.outname, par.gcm.clim, par.lat_interp);
        load(sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s/%s/grid.mat', type, par.outname, par.gcm.clim));

        for l = {'lo'}; land=l{1};
            var_vec = {'hfls', 'hfss', 'pr', 'evspsbl', 'lw', 'sw', 'rtoa', 'olr', 'lwsfc', 'swsfc'};
            for fn = var_vec; fname = fn{1};
                flux_z0i.flux_z.(land).(fname) = interp1(grid0.grid.dim3.lat, flux_z0.flux_z.(land).(fname), grid.dim3.lat);
                flux_z_mmm.(land).(fname) = nanmean(cat(3, flux_z0i.flux_z.(land).(fname), flux_z_mmm.(land).(fname)), 3);
            end
            for fn = {'ra', 'stf', 'res', 'r1', 'r2', 'ftoa', 'fsfc', 'sfc', 'shf', 'comp1', 'comp2'}; fname = fn{1};
                f_vec = par.gcm.fw;
                for f = f_vec; fw = f{1};
                    flux_z0i.flux_z.(land).(fname).(fw) = interp1(grid0.grid.dim3.lat, flux_z0.flux_z.(land).(fname).(fw), grid.dim3.lat);
                    flux_z_mmm.(land).(fname).(fw) = nanmean(cat(3, flux_z0i.flux_z.(land).(fname).(fw), flux_z_mmm.(land).(fname).(fw)), 3);
                end
            end
        end
    end

    flux_z = flux_z_mmm;

    % save energy flux data into mat file
    if ~exist(foldername, 'dir')
        mkdir(foldername)
    end
    save(sprintf('%s%s', foldername, 'flux_z'), 'flux_z', 'lat', '-v7.3');

end
 
