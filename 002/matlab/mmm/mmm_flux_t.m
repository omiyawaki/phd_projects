function mmm_flux_t(type, par)
    lat = par.lat;

    for l = {'lo'}; land=l{1};
        for t = {'ann', 'djf', 'mam', 'jja', 'son'}; time = t{1};
            var_vec = {'hfls', 'hfss', 'pr', 'evspsbl', 'lw', 'sw', 'rtoa', 'olr', 'lwsfc', 'swsfc'};
            for fn = var_vec; fname = fn{1};
                flux_t_mmm.(land).(time).(fname) = nan(length(par.lat), length(par.lon));
                flux_t_std.(land).(time).(fname) = nan(length(par.lat), length(par.lon));
            end
            for fn = {'ra', 'stf', 'res', 'r1', 'r2', 'ftoa', 'fsfc', 'sfc', 'shf', 'comp1', 'comp2'}; fname = fn{1};
                f_vec = par.gcm.fw;
                for f = f_vec; fw = f{1};
                    flux_t_mmm.(land).(time).(fname).(fw) = nan(length(par.lat), length(par.lon));
                    flux_t_std.(land).(time).(fname).(fw) = nan(length(par.lat), length(par.lon));
                end
            end
        end
    end

    pb = CmdLineProgressBar("Creating the multi-model mean...");
    for k=1:length(par.gcm_models); par.model = par.gcm_models{k};
        pb.print(k, length(par.gcm_models)); % output progress of moist adiabat calculonion

        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s/%s', type, par.model, par.gcm.clim);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s', type, par.model, par.gcm.clim);
        grid0 = load(sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s/%s/grid.mat', type, par.model, par.gcm.clim));
        flux_t0 = load(sprintf('%s/%s/flux_t.mat', prefix_proc, 'native')); % load lat x mon RCAE data

        foldername = sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s/%s/', type, par.outname, par.gcm.clim, par.lat_interp);
        load(sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s/%s/grid.mat', type, par.outname, par.gcm.clim));

        for l = {'lo'}; land=l{1};
            for t = {'ann', 'djf', 'mam', 'jja', 'son'}; time = t{1};
                var_vec = {'hfls', 'hfss', 'pr', 'evspsbl', 'lw', 'sw', 'rtoa', 'olr', 'lwsfc', 'swsfc'};
                for fn = var_vec; fname = fn{1};
                    flux_t0i.flux_t.(land).(time).(fname) = interp1(grid0.grid.dim3.lon, flux_t0.flux_t.(land).(time).(fname), grid.dim3.lon);
                    flux_t0i.flux_t.(land).(time).(fname) = permute(flux_t0i.flux_t.(land).(time).(fname), [2 1]);
                    flux_t0i.flux_t.(land).(time).(fname) = interp1(grid0.grid.dim3.lat, flux_t0i.flux_t.(land).(time).(fname), grid.dim3.lat);
                    flux_t0i.flux_t.(land).(time).(fname) = permute(flux_t0i.flux_t.(land).(time).(fname), [2 1]);
                    flux_t_mmm.(land).(time).(fname) = nanmean(cat(3, flux_t0i.flux_t.(land).(time).(fname), flux_t_mmm.(land).(time).(fname)), 3);
                    flux_t_std.(land).(time).(fname) = nanstd(cat(3, flux_t0i.flux_t.(land).(time).(fname), flux_t_std.(land).(time).(fname)), 3);
                end
                for fn = {'ra', 'stf', 'res', 'r1', 'r2', 'ftoa', 'fsfc', 'sfc', 'shf', 'comp1', 'comp2'}; fname = fn{1};
                    f_vec = par.gcm.fw;
                    for f = f_vec; fw = f{1};
                        flux_t0i.flux_t.(land).(time).(fname).(fw) = interp1(grid0.grid.dim3.lon, flux_t0.flux_t.(land).(time).(fname).(fw), grid.dim3.lon);
                        flux_t0i.flux_t.(land).(time).(fname).(fw) = permute(flux_t0i.flux_t.(land).(time).(fname).(fw), [2 1]);
                        flux_t0i.flux_t.(land).(time).(fname).(fw) = interp1(grid0.grid.dim3.lat, flux_t0i.flux_t.(land).(time).(fname).(fw), grid.dim3.lat);
                        flux_t0i.flux_t.(land).(time).(fname).(fw) = permute(flux_t0i.flux_t.(land).(time).(fname).(fw), [2 1]);
                        flux_t_mmm.(land).(time).(fname).(fw) = nanmean(cat(3, flux_t0i.flux_t.(land).(time).(fname).(fw), flux_t_mmm.(land).(time).(fname).(fw)), 3);
                        flux_t_std.(land).(time).(fname).(fw) = nanstd(cat(3, flux_t0i.flux_t.(land).(time).(fname).(fw), flux_t_std.(land).(time).(fname).(fw)), 3);
                    end
                end
            end % time
        end % land
    end % models

    flux_t = flux_t_mmm;

    % save energy flux data into mat file
    if ~exist(foldername, 'dir')
        mkdir(foldername)
    end
    save(sprintf('%s%s', foldername, 'flux_t'), 'flux_t', 'flux_t_std', 'lat', '-v7.3');

end
 
