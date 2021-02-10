function mmm_flux_zt(type, par)
    lat = par.lat;

    for l = {'lo'}; land=l{1};
        for t = {'ann', 'djf', 'mam', 'jja', 'son'}; time = t{1};
            var_vec = {'hfls', 'hfss', 'pr', 'evspsbl', 'lw', 'sw', 'rtoa', 'olr', 'lwsfc', 'swsfc'};
            for fn = var_vec; fname = fn{1};
                flux_zt_list.(land).(time).(fname) = nan(length(par.gcm_models), length(par.lat));
                flux_zt_mmm.(land).(time).(fname) = nan(1, length(par.lat));
                flux_zt_std.(land).(time).(fname) = nan(1, length(par.lat));
            end
            for fn = {'ra', 'stf', 'res', 'r1', 'r1z', 'r2', 'ftoa', 'fsfc', 'sfc', 'shf', 'comp1', 'comp2'}; fname = fn{1};
                f_vec = par.gcm.fw;
                for f = f_vec; fw = f{1};
                    flux_zt_list.(land).(time).(fname).(fw) = nan(length(par.gcm_models), length(par.lat));
                    flux_zt_mmm.(land).(time).(fname).(fw) = nan(1, length(par.lat));
                    flux_zt_std.(land).(time).(fname).(fw) = nan(1, length(par.lat));
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
        flux_zt0 = load(sprintf('%s/%s/flux_zt.mat', prefix_proc, 'native')); % load lat x mon RCAE data

        foldername = sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s/%s/', type, par.outname, par.gcm.clim, par.lat_interp);
        load(sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s/%s/grid.mat', type, par.outname, par.gcm.clim));

        for l = {'lo'}; land=l{1};
            for t = {'ann', 'djf', 'mam', 'jja', 'son'}; time = t{1};
                var_vec = {'hfls', 'hfss', 'pr', 'evspsbl', 'lw', 'sw', 'rtoa', 'olr', 'lwsfc', 'swsfc'};
                for fn = var_vec; fname = fn{1};
                    flux_zt0i.flux_zt.(land).(time).(fname) = interp1(grid0.grid.dim3.lat, flux_zt0.flux_zt.(land).(time).(fname), grid.dim3.lat);
                    %flux_zt_mmm.(land).(time).(fname) = nanmean(cat(1, flux_zt0i.flux_zt.(land).(time).(fname), flux_zt_mmm.(land).(time).(fname)), 1);
                    %flux_zt_std.(land).(time).(fname) = nanstd(cat(1, flux_zt0i.flux_zt.(land).(time).(fname), flux_zt_std.(land).(time).(fname)), 1);
                    flux_zt_list.(land).(time).(fname)(k,:) = flux_zt0i.flux_zt.(land).(time).(fname);
                end
                for fn = {'ra', 'stf', 'res', 'r1', 'r2', 'ftoa', 'fsfc', 'sfc', 'shf', 'comp1', 'comp2'}; fname = fn{1};
                    f_vec = par.gcm.fw;
                    for f = f_vec; fw = f{1};
                        flux_zt0i.flux_zt.(land).(time).(fname).(fw) = interp1(grid0.grid.dim3.lat, flux_zt0.flux_zt.(land).(time).(fname).(fw), grid.dim3.lat);
                        %flux_zt_mmm.(land).(time).(fname).(fw) = nanmean(cat(1, flux_zt0i.flux_zt.(land).(time).(fname).(fw), flux_zt_mmm.(land).(time).(fname).(fw)), 1);
                        %flux_zt_std.(land).(time).(fname).(fw) = nanstd(cat(1, flux_zt0i.flux_zt.(land).(time).(fname).(fw), flux_zt_std.(land).(time).(fname).(fw)), 1);
                        flux_zt_list.(land).(time).(fname).(fw)(k,:) = flux_zt0i.flux_zt.(land).(time).(fname).(fw);
                        
                        if strcmp(fname, 'r1')
                            %flux_zt_mmm.(land).(time).r1z.(fw) = nanmean(cat(1, flux_zt0i.flux_zt.(land).(time).res.(fw)./flux_zt0i.flux_zt.(land).(time).ra.(fw), flux_zt_mmm.(land).(time).r1z.(fw)), 1);
                            %flux_zt_std.(land).(time).r1z.(fw) = nanstd(cat(1, flux_zt0i.flux_zt.(land).(time).res.(fw)./flux_zt0i.flux_zt.(land).(time).ra.(fw), flux_zt_std.(land).(time).r1z.(fw)), 1);
                            flux_zt_list.(land).(time).r1z.(fw)(k,:) = flux_zt0i.flux_zt.(land).(time).res.(fw)./flux_zt0i.flux_zt.(land).(time).ra.(fw);
                        end
                    end
                end
            end % time
        end % land
    end % models
    
    for l = {'lo'}; land=l{1};
        for t = {'ann', 'djf', 'mam', 'jja', 'son'}; time = t{1};
            var_vec = {'hfls', 'hfss', 'pr', 'evspsbl', 'lw', 'sw', 'rtoa', 'olr', 'lwsfc', 'swsfc'};
            for fn = var_vec; fname = fn{1};
                flux_zt_mmm.(land).(time).(fname) = nanmean(flux_zt_list.(land).(time).(fname), 1);
                flux_zt_std.(land).(time).(fname) = nanstd(flux_zt_list.(land).(time).(fname), 1);
            end
            for fn = {'ra', 'stf', 'res', 'r1', 'r1z', 'r2', 'ftoa', 'fsfc', 'sfc', 'shf', 'comp1', 'comp2'}; fname = fn{1};
                f_vec = par.gcm.fw;
                for f = f_vec; fw = f{1};
                    flux_zt_mmm.(land).(time).(fname).(fw) = nanmean(flux_zt_list.(land).(time).(fname).(fw), 1);
                    flux_zt_std.(land).(time).(fname).(fw) = nanstd(flux_zt_list.(land).(time).(fname).(fw), 1);
                end
            end
        end % time
    end % land

    flux_zt = flux_zt_mmm;

    % save energy flux data into mat file
    if ~exist(foldername, 'dir')
        mkdir(foldername)
    end
    save(sprintf('%s%s', foldername, 'flux_zt'), 'flux_zt', 'flux_zt_std', 'lat', '-v7.3');

end

 
