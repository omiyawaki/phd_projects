function proc_ma(type, par)
% calculate moist adiabats
    if strcmp(type, 'era5') | strcmp(type, 'erai') | strcmp(type, 'era5c')
        foldername = sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/eps_%g_ga_%g/', type, par.lat_interp, par.ep, par.ga);
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s', type, par.(type).yr_span);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s', type, par.(type).yr_span);
    elseif strcmp(type, 'gcm')
        foldername = sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s/%s/eps_%g_ga_%g/', type, par.model, par.gcm.clim, par.lat_interp, par.ep, par.ga);
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s/%s', type, par.model, par.gcm.clim);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s', type, par.model, par.gcm.clim);
    elseif strcmp(type, 'echam')
        foldername = sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s/eps_%g_ga_%g/', type, par.echam.clim, par.lat_interp, par.ep, par.ga);
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s', type, par.echam.clim);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s', type, par.echam.clim, par.lat_interp);
    end
    load(sprintf('%s/grid.mat', prefix)); % read grid data
    load(sprintf('%s/srfc.mat', prefix)); % read surface variable data
    load(sprintf('%s/%s/masks.mat', prefix_proc, par.lat_interp)); % load land and ocean masks
    load(sprintf('%s/%s/eps_%g_ga_%g/rcae_t.mat', prefix_proc, par.lat_interp, par.ep, par.ga)); % load rcae data

    if strcmp(par.lat_interp, 'std')
        lat = par.lat_std;
    else
        lat = grid.dim3.lat;
    end

    for fn = fieldnames(srfc)'
        srfc.(fn{1}) = permute(srfc.(fn{1}), [2 1 3]); % bring lat to front
        srfc.(fn{1}) = interp1(grid.dim2.lat, srfc.(fn{1}), lat); % interpolate to standard grid
        srfc.(fn{1}) = permute(srfc.(fn{1}), [2 1 3]); % reorder to original dims
        for l = {'lo', 'l', 'o'}; land = l{1};
            if strcmp(land, 'lo'); srfc_n.(fn{1}).(land) = srfc.(fn{1});
            elseif strcmp(land, 'l'); srfc_n.(fn{1}).(land) = srfc.(fn{1}).*mask.ocean; % filter out ocean
            elseif strcmp(land, 'o'); srfc_n.(fn{1}).(land) = srfc.(fn{1}).*mask.land; % filter out land
            end
        end
    end

    for f = {'mse', 'dse'}; fw = f{1};
        for c = fieldnames(rcae_t.lo.ann.(fw))'; crit = c{1};
            for l = {'lo', 'l', 'o'}; land = l{1};
                for t = {'ann', 'djf', 'jja', 'mam', 'son'}; time = t{1};
                    for re = {'rce', 'rae'}; regime = re{1};
                        for v = fieldnames(srfc)'; vname = v{1};
                            if strcmp(time, 'ann')
                                srfc_t.(land).(time).(vname) = squeeze(nanmean(srfc_n.(vname).(land), 3));
                            elseif strcmp(time, 'djf')
                                srfc_shift = circshift(srfc_n.(vname).(land), 1, 3);
                                srfc_t.(land).(time).(vname) = squeeze(nanmean(srfc_shift(:,:,1:3), 3));
                            elseif strcmp(time, 'jja')
                                srfc_t.(land).(time).(vname) = squeeze(nanmean(srfc_n.(vname).(land)(:,:,6:8), 3));
                            elseif strcmp(time, 'mam')
                                srfc_t.(land).(time).(vname) = squeeze(nanmean(srfc_n.(vname).(land)(:,:,3:5), 3));
                            elseif strcmp(time, 'son')
                                srfc_t.(land).(time).(vname) = squeeze(nanmean(srfc_n.(vname).(land)(:,:,9:11), 3));
                            end

                            filt.(land).(time).(fw).(crit).(regime) = nan(size(rcae_t.(land).(time).(fw).(crit)));
                            if strcmp(regime, 'rce'); filt.(land).(time).(fw).(crit).(regime)(rcae_t.(land).(time).(fw).(crit)==1)=1; % set RCE=1, elsewhere nan
                            elseif strcmp(regime, 'rae'); filt.(land).(time).(fw).(crit).(regime)(rcae_t.(land).(time).(fw).(crit)==-1)=1; % set RAE=1, elsewhere nan
                            end

                            srfc_tf.(land).(time).(fw).(crit).(regime).(vname) = srfc_t.(land).(time).(vname) .* filt.(land).(time).(fw).(crit).(regime);

                            nanfilt.(regime) = nan(size(srfc_tf.(land).(time).(fw).(crit).(regime).(vname)));
                            nanfilt.(regime)(~isnan(srfc_tf.(land).(time).(fw).(crit).(regime).(vname))) = 1;

                            % take cosine-weighted average
                            for d = {'all', 'nh', 'sh', 'tp'}; domain = d{1};
                                if strcmp(regime, 'rce')
                                    if strcmp(domain, 'all')
                                        cosw = repmat(cosd(lat)', [size(nanfilt.(regime), 1) 1]);
                                        denm = srfc_tf.(land).(time).(fw).(crit).(regime).(vname);
                                        nume = nanfilt.(regime);
                                    elseif strcmp(domain, 'nh')
                                        cosw = repmat(cosd(lat(lat>30))', [size(nanfilt.(regime), 1) 1]);
                                        denm = srfc_tf.(land).(time).(fw).(crit).(regime).(vname)(:, lat>30);
                                        nume = nanfilt.(regime)(:, lat>30);
                                    elseif strcmp(domain, 'sh')
                                        cosw = repmat(cosd(lat(lat<-30))', [size(nanfilt.(regime), 1) 1]);
                                        denm = srfc_tf.(land).(time).(fw).(crit).(regime).(vname)(:, lat<-30);
                                        nume = nanfilt.(regime)(:, lat<-30);
                                    elseif strcmp(domain, 'tp')
                                        cosw = repmat(cosd(lat(abs(lat)<10))', [size(nanfilt.(regime), 1) 1]);
                                        denm = srfc_tf.(land).(time).(fw).(crit).(regime).(vname)(:, abs(lat)<10);
                                        nume = nanfilt.(regime)(:, abs(lat)<10);
                                    end
                                elseif strcmp(regime, 'rae')
                                    if strcmp(domain, 'all')
                                        cosw = repmat(cosd(lat)', [size(nanfilt.(regime), 1) 1]);
                                        denm = srfc_tf.(land).(time).(fw).(crit).(regime).(vname);
                                        nume = nanfilt.(regime);
                                    elseif strcmp(domain, 'nh')
                                        cosw = repmat(cosd(lat(lat>0))', [size(nanfilt.(regime), 1) 1]);
                                        denm = srfc_tf.(land).(time).(fw).(crit).(regime).(vname)(:, lat>0);
                                        nume = nanfilt.(regime)(:, lat>0);
                                    elseif strcmp(domain, 'sh')
                                        cosw = repmat(cosd(lat(lat<0))', [size(nanfilt.(regime), 1) 1]);
                                        denm = srfc_tf.(land).(time).(fw).(crit).(regime).(vname)(:, lat<0);
                                        nume = nanfilt.(regime)(:, lat<0);
                                    end
                                end

                                ma.(regime).(domain).(fw).(crit).(land).(time).(vname) = nansum(cosw.*denm, 2) ./ nansum(cosw.*nume, 2); % weighted meridional average
                                ma.(regime).(domain).(fw).(crit).(land).(time).(vname) = squeeze(nanmean(ma.(regime).(domain).(fw).(crit).(land).(time).(vname), 1)); % zonal average

                            end % end domain loop
                        end % end srfc variables loop

                        if strcmp(type, 'era5') | strcmp(type, 'erai') | strcmp(type, 'era5c')
                            ma.rce.all.(fw).(crit).(land).(time).ta = calc_ma_dew(ma.rce.all.(fw).(crit).(land).(time), grid.dim3.plev, par, type); % compute moist adiabat with dew point temperature
                            ma.rce.tp.(fw).(crit).(land).(time).ta = calc_ma_dew(ma.rce.tp.(fw).(crit).(land).(time), grid.dim3.plev, par, type); % compute moist adiabat with dew point temperature
                            ma.rce.nh.(fw).(crit).(land).(time).ta = calc_ma_dew(ma.rce.nh.(fw).(crit).(land).(time), grid.dim3.plev, par, type); % compute moist adiabat with dew point temperature
                            ma.rce.sh.(fw).(crit).(land).(time).ta = calc_ma_dew(ma.rce.sh.(fw).(crit).(land).(time), grid.dim3.plev, par, type); % compute moist adiabat with dew point temperature
                        elseif strcmp(type, 'gcm')
                            ma.rce.all.(fw).(crit).(land).(time).ta = calc_ma_hurs(ma.rce.all.(fw).(crit).(land).(time), grid.dim3.plev, par); % compute moist adiabat with RH
                            ma.rce.tp.(fw).(crit).(land).(time).ta = calc_ma_hurs(ma.rce.tp.(fw).(crit).(land).(time), grid.dim3.plev, par); % compute moist adiabat with RH
                            ma.rce.nh.(fw).(crit).(land).(time).ta = calc_ma_hurs(ma.rce.nh.(fw).(crit).(land).(time), grid.dim3.plev, par); % compute moist adiabat with RH
                            ma.rce.sh.(fw).(crit).(land).(time).ta = calc_ma_hurs(ma.rce.sh.(fw).(crit).(land).(time), grid.dim3.plev, par); % compute moist adiabat with RH
                        end

                    end % end RCE/RAE regime loop
                end % end time average loop
            end % end land option loop
        end % end RCAE definition loop
    end % end MSE/DSE framework loop

    % save data into mat file
    printname = [foldername 'ma.mat'];
    if ~exist(foldername, 'dir')
        mkdir(foldername)
    end
    save(printname, 'ma');

end % process moist adiabat in RCE/RAE regimes
