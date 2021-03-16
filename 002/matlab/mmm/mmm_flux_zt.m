function mmm_flux_zt(type, par)
    lat = par.lat;

    % output info
    foldername = make_savedir_proc(type, par);
    load(sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s/%s/grid.mat', type, par.outname, par.clim));

    par.lat_interp = 'native'; % input files will be in native grid

    for l = {'lo'}; land=l{1};
        for t = {'ann', 'djf', 'mam', 'jja', 'son'}; time = t{1};
            var_vec = {'hfls', 'hfss', 'pr', 'evspsbl', 'lw', 'sw', 'rtoa', 'olr', 'lwsfc', 'swsfc'};
            for fn = var_vec; fname = fn{1};
                flux_zt_list.(land).(time).(fname) = nan(length(par.model_list), length(par.lat));
                flux_zt_mmm.(land).(time).(fname) = nan(1, length(par.lat));
                flux_zt_std.(land).(time).(fname) = nan(1, length(par.lat));
            end
            for fn = {'ra', 'stf', 'res', 'r1', 'r1z', 'r2', 'ftoa', 'fsfc', 'sfc', 'shf', 'comp1', 'comp2'}; fname = fn{1};
                f_vec = par.gcm.fw;
                for f = f_vec; fw = f{1};
                    flux_zt_list.(land).(time).(fname).(fw) = nan(length(par.model_list), length(par.lat));
                    flux_zt_mmm.(land).(time).(fname).(fw) = nan(1, length(par.lat));
                    flux_zt_std.(land).(time).(fname).(fw) = nan(1, length(par.lat));
                end
            end
        end
    end

    pb = CmdLineProgressBar("Creating the multi-model mean...");
    for k=1:length(par.model_list); par.model = par.model_list{k};
        pb.print(k, length(par.model_list)); % output progress of model cycle

        if strcmp(type, 'gcm')
            type_in = type;
        else
            type_in = par.model;
        end

        % input info
        prefix = make_prefix(type_in, par);
        prefix_proc = make_prefix_proc(type_in, par);
        grid0 = load(sprintf('%s/grid.mat', prefix));
        flux_zt0 = load(sprintf('%s/flux_zt.mat', prefix_proc)); % load lat x mon RCAE data


        for l = {'lo'}; land=l{1};
            for t = {'ann', 'djf', 'mam', 'jja', 'son'}; time = t{1};
                f_vec = par.(type).fw;
                for f = f_vec; fw = f{1};

                    var_vec = {'hfls', 'hfss', 'pr', 'evspsbl', 'lw', 'sw', 'rtoa', 'olr', 'lwsfc', 'swsfc'};
                    var_vec_in = make_varvec(type_in, fw);
                    for fn = 1:length(var_vec); fname = var_vec{fn}; fname_in = var_vec_in{fn};
                        if any(strcmp(fname_in, {'sshf', 'slhf'}))
                            flux_zt0i.flux_zt.(land).(time).(fname_in) = -interp1(grid0.grid.dim3.lat, flux_zt0.flux_zt.(land).(time).(fname_in), grid.dim3.lat);
                        else
                            flux_zt0i.flux_zt.(land).(time).(fname_in) = interp1(grid0.grid.dim3.lat, flux_zt0.flux_zt.(land).(time).(fname_in), grid.dim3.lat);
                        end
                        flux_zt_list.(land).(time).(fname)(k,:) = flux_zt0i.flux_zt.(land).(time).(fname_in);
                    end
                    for fn = {'ra', 'stf', 'res', 'r1', 'r2', 'ftoa', 'fsfc', 'sfc', 'shf', 'comp1', 'comp2'}; fname = fn{1};
                        flux_zt0i.flux_zt.(land).(time).(fname).(fw) = interp1(grid0.grid.dim3.lat, flux_zt0.flux_zt.(land).(time).(fname).(fw), grid.dim3.lat);
                        flux_zt_list.(land).(time).(fname).(fw)(k,:) = flux_zt0i.flux_zt.(land).(time).(fname).(fw);
                        
                        if strcmp(fname, 'r1')
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
            for fn = 1:length(var_vec); fname = var_vec{fn};
                flux_zt_mmm.(land).(time).(fname) = squeeze(nanmean(flux_zt_list.(land).(time).(fname), 1));
                if strcmp(type, 'rea')
                    flux_zt_std.(land).(time).(fname) = squeeze(range(flux_zt_list.(land).(time).(fname), 1));
                else
                    flux_zt_std.(land).(time).(fname) = squeeze(nanstd(flux_zt_list.(land).(time).(fname), 1));
                end
            end
            for fn = {'ra', 'stf', 'res', 'r1', 'r1z', 'r2', 'ftoa', 'fsfc', 'sfc', 'shf', 'comp1', 'comp2'}; fname = fn{1};
                f_vec = par.(type).fw;
                for f = f_vec; fw = f{1};
                    flux_zt_mmm.(land).(time).(fname).(fw) = squeeze(nanmean(flux_zt_list.(land).(time).(fname).(fw), 1));
                    if strcmp(type, 'rea')
                        flux_zt_std.(land).(time).(fname).(fw) = squeeze(range(flux_zt_list.(land).(time).(fname).(fw), 1));
                    else
                        flux_zt_std.(land).(time).(fname).(fw) = squeeze(nanstd(flux_zt_list.(land).(time).(fname).(fw), 1));
                    end
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

 
