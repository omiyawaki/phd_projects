function mmm_flux_t(type, par)
    lat = par.lat;

    % output info
    foldername = make_savedir_proc(type, par);
    load(sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s/%s/grid.mat', type, par.outname, par.clim));

    par.lat_interp = 'native'; % input files will be in native grid

    for l = par.land_list; land=l{1};
        for t = {'ann', 'djf', 'mam', 'jja', 'son'}; time = t{1};
            % var_vec = {'hfls', 'hfss', 'lw', 'sw', 'rtoa', 'olr', 'lwsfc', 'swsfc', 'tend'};
            var_vec = {'hfls', 'hfss', 'lw', 'sw', 'rtoa', 'olr', 'lwsfc', 'swsfc'};
            for fn = var_vec; fname = fn{1};
                flux_t_list.(land).(time).(fname) = nan(length(par.model_list), length(par.lon), length(par.lat));
                flux_t_mmm.(land).(time).(fname) = nan(length(par.lon), length(par.lat));
                flux_t_std.(land).(time).(fname) = nan(length(par.lon), length(par.lat));
            end
            for fn = {'ra', 'stf', 'res', 'r1', 'r2', 'ftoa', 'fsfc', 'sfc', 'shf', 'comp1', 'comp2'}; fname = fn{1};
                f_vec = par.gcm.fw;
                for f = f_vec; fw = f{1};
                    flux_t_list.(land).(time).(fname).(fw) = nan(length(par.model_list), length(par.lon), length(par.lat));
                    flux_t_mmm.(land).(time).(fname).(fw) = nan(length(par.lon), length(par.lat));
                    flux_t_std.(land).(time).(fname).(fw) = nan(length(par.lon), length(par.lat));
                end
            end
        end
    end

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
        flux_t0 = load(sprintf('%s/flux_t.mat', prefix_proc));

        for l = par.land_list; land=l{1};
            for t = {'ann', 'djf', 'mam', 'jja', 'son'}; time = t{1};
                f_vec = par.(type).fw;
                for f = f_vec; fw = f{1};
                    % var_vec = {'hfls', 'hfss', 'lw', 'sw', 'rtoa', 'olr', 'lwsfc', 'swsfc', 'tend'};
                    var_vec = {'hfls', 'hfss', 'lw', 'sw', 'rtoa', 'olr', 'lwsfc', 'swsfc'};
                    var_vec_in = make_varvec(type_in, fw);
                    for fn = 1:length(var_vec); fname = var_vec{fn}; fname_in = var_vec_in{fn};
                        if any(strcmp(fname_in, {'sshf', 'slhf'}))
                            flux_t0i.flux_t.(land).(time).(fname_in) = -interp1(grid0.grid.dim3.lon, flux_t0.flux_t.(land).(time).(fname_in), grid.dim3.lon);
                        else
                            flux_t0i.flux_t.(land).(time).(fname_in) = interp1(grid0.grid.dim3.lon, flux_t0.flux_t.(land).(time).(fname_in), grid.dim3.lon);
                        end
                        flux_t0i.flux_t.(land).(time).(fname_in) = permute(flux_t0i.flux_t.(land).(time).(fname_in), [2 1]);
                        flux_t0i.flux_t.(land).(time).(fname_in) = interp1(grid0.grid.dim3.lat, flux_t0i.flux_t.(land).(time).(fname_in), grid.dim3.lat);
                        flux_t0i.flux_t.(land).(time).(fname) = permute(flux_t0i.flux_t.(land).(time).(fname_in), [2 1]);
                        
                        flux_t_list.(land).(time).(fname)(k,:,:) = flux_t0i.flux_t.(land).(time).(fname);
                    end
                for fn = {'ra', 'stf', 'res', 'r1', 'r2', 'ftoa', 'fsfc', 'sfc', 'shf', 'comp1', 'comp2'}; fname = fn{1};
                    flux_t0i.flux_t.(land).(time).(fname).(fw) = interp1(grid0.grid.dim3.lon, flux_t0.flux_t.(land).(time).(fname).(fw), grid.dim3.lon);
                    flux_t0i.flux_t.(land).(time).(fname).(fw) = permute(flux_t0i.flux_t.(land).(time).(fname).(fw), [2 1]);
                    flux_t0i.flux_t.(land).(time).(fname).(fw) = interp1(grid0.grid.dim3.lat, flux_t0i.flux_t.(land).(time).(fname).(fw), grid.dim3.lat);
                    flux_t0i.flux_t.(land).(time).(fname).(fw) = permute(flux_t0i.flux_t.(land).(time).(fname).(fw), [2 1]);

                    flux_t_list.(land).(time).(fname).(fw)(k,:,:) = flux_t0i.flux_t.(land).(time).(fname).(fw);
                end
            end % time
        end % land
    end % models

    for l = par.land_list; land=l{1};
        for t = {'ann', 'djf', 'mam', 'jja', 'son'}; time = t{1};
            % var_vec = {'hfls', 'hfss', 'lw', 'sw', 'rtoa', 'olr', 'lwsfc', 'swsfc', 'tend'};
            var_vec = {'hfls', 'hfss', 'lw', 'sw', 'rtoa', 'olr', 'lwsfc', 'swsfc'};
            for fn = 1:length(var_vec); fname = var_vec{fn};
                flux_t_mmm.(land).(time).(fname) = squeeze(nanmean(flux_t_list.(land).(time).(fname), 1));
                if strcmp(type, 'rea')
                    flux_t_std.(land).(time).(fname) = squeeze(range(flux_t_list.(land).(time).(fname), 1));
                else
                    flux_t_std.(land).(time).(fname) = squeeze(nanstd(flux_t_list.(land).(time).(fname), 1));
                end
            end

            for fn = {'ra', 'stf', 'res', 'r1', 'r2', 'ftoa', 'fsfc', 'sfc', 'shf', 'comp1', 'comp2'}; fname = fn{1};
                f_vec = par.(type).fw;
                for f = f_vec; fw = f{1};
                    flux_t_mmm.(land).(time).(fname).(fw) = squeeze(nanmean(flux_t_list.(land).(time).(fname).(fw), 1));
                    if strcmp(type, 'rea')
                        flux_t_std.(land).(time).(fname).(fw) = squeeze(range(flux_t_list.(land).(time).(fname).(fw), 1));
                    else
                        flux_t_std.(land).(time).(fname).(fw) = squeeze(nanstd(flux_t_list.(land).(time).(fname).(fw), 1));
                    end
                end
            end
        end
    end

    flux_t = flux_t_mmm;

    % save energy flux data into mat file
    if ~exist(foldername, 'dir')
        mkdir(foldername)
    end
    save(sprintf('%s%s', foldername, 'flux_t'), 'flux_t', 'flux_t_std', 'lat', '-v7.3');

end
 
