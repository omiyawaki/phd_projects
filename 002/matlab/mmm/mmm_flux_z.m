function mmm_flux_z(type, par)
    lat = par.lat;

    % output info
    foldername = make_savedir_proc(type, par);
    load(sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s/%s/grid.mat', type, par.outname, par.clim));

    par.lat_interp = 'native'; % input files will be in native grid

    for l = {'lo'}; land=l{1};
        var_vec = {'hfls', 'hfss', 'pr', 'evspsbl', 'lw', 'sw', 'rtoa', 'olr', 'lwsfc', 'swsfc'};
        for fn = var_vec; fname = fn{1};
            flux_z_list.(land).(fname) = nan(length(par.model_list), length(par.lat), 12);
            flux_z_mmm.(land).(fname) = nan(length(par.lat), 12);
            flux_z_std.(land).(fname) = nan(length(par.lat), 12);
        end
        for fn = {'ra', 'stf', 'res', 'r1', 'r2', 'ftoa', 'fsfc', 'sfc', 'shf', 'comp1', 'comp2'}; fname = fn{1};
            f_vec = par.(type).fw;
            for f = f_vec; fw = f{1};
                flux_z_list.(land).(fname).(fw) = nan(length(par.model_list), length(par.lat), 12);
                flux_z_mmm.(land).(fname).(fw) = nan(length(par.lat), 12);
                flux_z_std.(land).(fname).(fw) = nan(length(par.lat), 12);
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
        flux_z0 = load(sprintf('%s/flux_z.mat', prefix_proc));

        for l = {'lo'}; land=l{1};
            f_vec = par.(type).fw;
            for f = f_vec; fw = f{1};

                var_vec = {'hfls', 'hfss', 'pr', 'evspsbl', 'lw', 'sw', 'rtoa', 'olr', 'lwsfc', 'swsfc'};
                var_vec_in = make_varvec(type_in, fw);
                for fn = 1:length(var_vec); fname = var_vec{fn}; fname_in = var_vec_in{fn};
                    if any(strcmp(fname_in, {'sshf', 'slhf'}))
                        flux_z0i.flux_z.(land).(fname_in) = -interp1(grid0.grid.dim3.lat, flux_z0.flux_z.(land).(fname_in), grid.dim3.lat);
                    else
                        flux_z0i.flux_z.(land).(fname_in) = interp1(grid0.grid.dim3.lat, flux_z0.flux_z.(land).(fname_in), grid.dim3.lat);
                    end
                    flux_z_list.(land).(fname)(k,:,:) = flux_z0i.flux_z.(land).(fname_in);
                end
                for fn = {'ra', 'stf', 'res', 'r1', 'r2', 'ftoa', 'fsfc', 'sfc', 'shf', 'comp1', 'comp2'}; fname = fn{1};
                    flux_z0i.flux_z.(land).(fname).(fw) = interp1(grid0.grid.dim3.lat, flux_z0.flux_z.(land).(fname).(fw), grid.dim3.lat);
                    flux_z_list.(land).(fname).(fw)(k,:,:) = flux_z0i.flux_z.(land).(fname).(fw);
                end

            end
        end

    end

    for l = {'lo'}; land=l{1};
        var_vec = {'hfls', 'hfss', 'pr', 'evspsbl', 'lw', 'sw', 'rtoa', 'olr', 'lwsfc', 'swsfc'};
        for fn = 1:length(var_vec); fname = var_vec{fn};
            flux_z_mmm.(land).(fname) = squeeze(nanmean(flux_z_list.(land).(fname), 1));
            if strcmp(type, 'rea')
                flux_z_std.(land).(fname) = squeeze(range(flux_z_list.(land).(fname), 1));
            else
                flux_z_std.(land).(fname) = squeeze(nanstd(flux_z_list.(land).(fname), 1));
            end
        end
        for fn = {'ra', 'stf', 'res', 'r1', 'r2', 'ftoa', 'fsfc', 'sfc', 'shf', 'comp1', 'comp2'}; fname = fn{1};
            f_vec = par.(type).fw;
            for f = f_vec; fw = f{1};
                flux_z_mmm.(land).(fname).(fw) = squeeze(nanmean(flux_z_list.(land).(fname).(fw), 1));
                if strcmp(type, 'rea')
                    flux_z_std.(land).(fname).(fw) = squeeze(range(flux_z_list.(land).(fname).(fw), 1));
                else
                    flux_z_std.(land).(fname).(fw) = squeeze(nanstd(flux_z_list.(land).(fname).(fw), 1));
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