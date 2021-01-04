function proc_ta_si(type, par)
    if strcmp(type, 'era5') | strcmp(type, 'erai') | strcmp(type, 'era5c')
        foldername = sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/eps_%g_ga_%g/', type, par.lat_interp, par.ep, par.ga);
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s', type, par.(type).yr_span);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s', type, par.(type).yr_span);
    elseif strcmp(type, 'merra2')
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
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s', type, par.echam.clim);
    end
    load(sprintf('%s/grid.mat', prefix)); % read ERA5 grid data
    load(sprintf('%s/srfc.mat', prefix)); % load surface data
    ta_orig = load(sprintf('%s/ta_si.mat', prefix)); ta_si=ta_orig.ta_si.spl; clear ta_orig; % load temperature in sigma
    % load(sprintf('%s/%s/masks.mat', prefix_proc, par.lat_interp)); % load land and ocean masks
    load(sprintf('%s/%s/eps_%g_ga_%g/rcae_alt_t.mat', prefix_proc, par.lat_interp, par.ep, par.ga)); % load rcae data

    % mask.land_vert = repmat(mask.land, [1 1 1 size(ta_si, 3)]); % expand land mask to vertical dim
    % mask.land_vert = permute(mask.land, [1 2 4 3]); % place vertical dim where it belongs
    % mask.ocean_vert = repmat(mask.ocean, [1 1 1 size(ta_si, 3)]); % expand ocean mask to vertical dim
    % mask.ocean_vert = permute(mask.ocean, [1 2 4 3]); % place vertical dim where it belongs

    % ta_si_sm.l = ta_si.*mask.ocean_vert; % filter ta_si with surface mask
    % ta_si_sm.o = ta_si.*mask.land_vert; % filter ta_si with surface mask
    ta_si_sm.lo = permute(ta_si, [1 2 4 3]); % bring plev to last dimension
    % ta_si_sm.l = permute(ta_si_sm.l, [1 2 4 3]); % bring plev to last dimension
    % ta_si_sm.o = permute(ta_si_sm.o, [1 2 4 3]); % bring plev to last dimension

    clear ta_si;

    % mask_t.land = nanmean(mask.land, 3);
    % mask_t.ocean = nanmean(mask.ocean, 3);
    f_vec = assign_fw(type, par);
    for f = f_vec; fw = f{1};
        for c = fieldnames(rcae_alt_t.lo.ann.(fw))'; crit = c{1};
            % for l = {'lo', 'l', 'o'}; land = l{1}; % over land, over ocean, or both
            for l = {'lo'}; land = l{1}; % over land, over ocean, or both
                % for t = {'ann', 'djf', 'jja', 'mam', 'son'}; time = t{1};
                for t = {'ann'}; time = t{1};
                    for re = {'rce', 'rae', 'rcae'}; regime = re{1};
                        if strcmp(time, 'ann')
                            ta_si_n.(land).(time) = squeeze(nanmean(ta_si_sm.(land), 3));
                        elseif strcmp(time, 'djf')
                            ta_si_shift = circshift(ta_si_sm.(land), 1, 3);
                            ta_si_n.(land).(time) = squeeze(nanmean(ta_si_shift(:,:,1:3,:), 3));
                        elseif strcmp(time, 'jja')
                            ta_si_n.(land).(time) = squeeze(nanmean(ta_si_sm.(land)(:,:,6:8,:), 3));
                        elseif strcmp(time, 'mam')
                            ta_si_n.(land).(time) = squeeze(nanmean(ta_si_sm.(land)(:,:,3:5,:), 3));
                        elseif strcmp(time, 'son')
                            ta_si_n.(land).(time) = squeeze(nanmean(ta_si_sm.(land)(:,:,9:11,:), 3));
                        end

                        filt.(land).(time).(fw).(crit).(regime) = nan(size(rcae_alt_t.(land).(time).(fw).(crit))); % create empty arrays to store filtering array
                        if strcmp(regime, 'rce'); filt.(land).(time).(fw).(crit).(regime)(rcae_alt_t.(land).(time).(fw).(crit)==1)=1; % set RCE=1, elsewhere nan
                        elseif strcmp(regime, 'rae'); filt.(land).(time).(fw).(crit).(regime)(rcae_alt_t.(land).(time).(fw).(crit)==-1)=1; % set RAE=1, elsewhere nan
                        elseif strcmp(regime, 'rcae'); filt.(land).(time).(fw).(crit).(regime)(rcae_alt_t.(land).(time).(fw).(crit)==0)=1; % set RCAE=1, elsewhere nan
                        end

                        ta_si_t.(land).(time).(fw).(crit).(regime) = ta_si_n.(land).(time) .* filt.(land).(time).(fw).(crit).(regime);

                        nanfilt.(regime) = nan(size(ta_si_t.(land).(time).(fw).(crit).(regime)));
                        nanfilt.(regime)(~isnan(ta_si_t.(land).(time).(fw).(crit).(regime))) = 1;

                        % take cosine-weighted meridional average
                        for d = {'all', 'nh', 'sh', 'tp'}; domain = d{1};
                            if strcmp(regime, 'rce')
                                if strcmp(domain, 'all')
                                    cosw = repmat(cosd(lat)', [size(nanfilt.(regime), 1) 1]);
                                    denm = ta_si_t.(land).(time).(fw).(crit).(regime);
                                    nume = nanfilt.(regime);
                                elseif strcmp(domain, 'nh')
                                    cosw = repmat(cosd(lat(lat>0))', [size(nanfilt.(regime), 1) 1]);
                                    denm = ta_si_t.(land).(time).(fw).(crit).(regime)(:,lat>0,:);
                                    nume = nanfilt.(regime)(:,lat>0,:);
                                elseif strcmp(domain, 'sh')
                                    cosw = repmat(cosd(lat(lat<0))', [size(nanfilt.(regime), 1) 1]);
                                    denm = ta_si_t.(land).(time).(fw).(crit).(regime)(:,lat<0,:);
                                    nume = nanfilt.(regime)(:,lat<0,:);
                                elseif strcmp(domain, 'tp')
                                    cosw = repmat(cosd(lat(abs(lat)<10))', [size(nanfilt.(regime), 1) 1]);
                                    denm = ta_si_t.(land).(time).(fw).(crit).(regime)(:,abs(lat)<10,:);
                                    nume = nanfilt.(regime)(:,abs(lat)<10,:);
                                end
                            elseif strcmp(regime, 'rae')
                                if strcmp(domain, 'all')
                                    cosw = repmat(cosd(lat)', [size(nanfilt.(regime), 1) 1]);
                                    denm = ta_si_t.(land).(time).(fw).(crit).(regime);
                                    nume = nanfilt.(regime);
                                elseif strcmp(domain, 'nh')
                                    cosw = repmat(cosd(lat(lat>0))', [size(nanfilt.(regime), 1) 1]);
                                    denm = ta_si_t.(land).(time).(fw).(crit).(regime)(:,lat>0,:);
                                    nume = nanfilt.(regime)(:,lat>0,:);
                                elseif strcmp(domain, 'sh')
                                    cosw = repmat(cosd(lat(lat<0))', [size(nanfilt.(regime), 1) 1]);
                                    denm = ta_si_t.(land).(time).(fw).(crit).(regime)(:,lat<0,:);
                                    nume = nanfilt.(regime)(:,lat<0,:);
                                end
                            elseif strcmp(regime, 'rcae')
                                if strcmp(domain, 'all')
                                    cosw = repmat(cosd(lat)', [size(nanfilt.(regime), 1) 1]);
                                    denm = ta_si_t.(land).(time).(fw).(crit).(regime);
                                    nume = nanfilt.(regime);
                                elseif strcmp(domain, 'nh')
                                    cosw = repmat(cosd(lat(lat>0))', [size(nanfilt.(regime), 1) 1]);
                                    denm = ta_si_t.(land).(time).(fw).(crit).(regime)(:,lat>0,:);
                                    nume = nanfilt.(regime)(:,lat>0,:);
                                elseif strcmp(domain, 'sh')
                                    cosw = repmat(cosd(lat(lat<0))', [size(nanfilt.(regime), 1) 1]);
                                    denm = ta_si_t.(land).(time).(fw).(crit).(regime)(:,lat<0,:);
                                    nume = nanfilt.(regime)(:,lat<0,:);
                                end
                            end

                            ta_si.(regime).(domain).(fw).(crit).(land).(time) = nansum(cosw.*denm, 2) ./ nansum(cosw.*nume, 2); % area-weighted meridional average
                            ta_si.(regime).(domain).(fw).(crit).(land).(time) = squeeze(nanmean(ta_si.(regime).(domain).(fw).(crit).(land).(time), 1)); % zonal average
                        end
                    end
                end
            end
        end
    end

    % save filtered data
    printname = [foldername 'ta_si'];
    if ~exist(foldername, 'dir')
        mkdir(foldername)
    end
    save(printname, 'ta_si');
end % process temperature in RCE/RAE regimes
