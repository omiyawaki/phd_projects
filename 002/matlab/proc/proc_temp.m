function proc_temp(type, par)
    if strcmp(type, 'era5') | strcmp(type, 'erai') | strcmp(type, 'era5c')
        foldername = sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/eps_%g_ga_%g/', type, par.lat_interp, par.ep, par.ga);
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s', type, par.(type).yr_span);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s', type, par.(type).yr_span);
        file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/%s/temp/%s_temp_%s.ymonmean.nc', type, type, par.(type).yr_span));
        fullpath=sprintf('%s/%s', file.folder, file.name);
        temp = double(ncread(fullpath, 't'));
    elseif strcmp(type, 'gcm')
        foldername = sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s/%s/eps_%g_ga_%g/', type, par.model, par.gcm.clim, par.lat_interp, par.ep, par.ga);
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s/%s', type, par.model, par.gcm.clim);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s', type, par.model, par.gcm.clim);
        var = 'ta';
        file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/gcm/%s/%s_Amon_%s_%s_r1i1p1_*.ymonmean.nc', par.model, var, par.model, par.gcm.clim));
        fullpath=sprintf('%s/%s', file.folder, file.name);
        temp = double(ncread(fullpath, var));
    elseif strcmp(type, 'echam')
        foldername = sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s/', type, par.echam.clim, par.lat_interp);
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s', type, par.echam.clim);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s', type, par.echam.clim);
        var = 't';
        file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/echam/ATM*_%s_*.ymonmean.nc', par.echam.clim));
        fullpath=sprintf('%s/%s', file.folder, file.name);
        temp = double(ncread(fullpath, var));
    end
    load(sprintf('%s/grid.mat', prefix)); % read ERA5 grid data
    load(sprintf('%s/srfc.mat', prefix)); % load surface data
    load(sprintf('%s/%s/masks.mat', prefix_proc, par.lat_interp)); % load land and ocean masks
    load(sprintf('%s/%s/eps_%g_ga_%g/rcae_t.mat', prefix_proc, par.lat_interp, par.ep, par.ga)); % load rcae data

    if strcmp(par.lat_interp, 'std')
        lat = par.lat_std;
    else
        lat = grid.dim3.lat;
    end

    % interpolate temp to standard lat grid
    temp = permute(temp, [2 1 3 4]);
    temp = interp1(grid.dim3.lat, temp, lat);
    temp = permute(temp, [2 1 3 4]);

    % create surface mask
    if strcmp(type, 'era5') | strcmp(type, 'erai') | strcmp(type, 'era5c')
        ps = permute(srfc.sp, [2 1 3]); % bring lat to front to interpolate
        ps = interp1(grid.dim2.lat, ps, lat, 'linear', 'extrap'); % interpolate to standard grid
        ps = permute(ps, [2 1 3]); % reorder to original
        ps_vert = repmat(ps, [1 1 1 size(temp, 3)]); % dims (lon x lat x time x plev)
        ps_vert = permute(ps_vert, [1 2 4 3]); % dims (lon x lat x plev x time)
        pa = double(permute(repmat(grid.dim3.plev, [1 size(ps)]), [2 3 1 4]));
        ts = permute(srfc.t2m, [2 1 3]); % bring lat to front to interpolate
        ts = interp1(grid.dim2.lat, ts, lat); % interpolate to standard grid
        ts = permute(ts, [2 1 3]); % reorder to original
        ts_vert = repmat(ts, [1 1 1 size(temp, 3)]); % dims (lon x lat x time x plev)
        ts_vert = permute(ts_vert, [1 2 4 3]); % dims (lon x lat x plev x time)
    elseif strcmp(type, 'gcm')
        ps = permute(srfc.ps, [2 1 3]); % bring lat to front to interpolate
        ps = interp1(grid.dim2.lat, ps, lat, 'linear', 'extrap'); % interpolate to standard grid
        ps = permute(ps, [2 1 3]); % reorder to original
        ps_vert = repmat(ps, [1 1 1 size(temp, 3)]); % dims (lon x lat x time x plev)
        ps_vert = permute(ps_vert, [1 2 4 3]); % dims (lon x lat x plev x time)
        pa = permute(repmat(grid.dim3.plev, [1 size(ps)]), [2 3 1 4]);
        ts = permute(srfc.tas, [2 1 3]); % bring lat to front to interpolate
        ts = interp1(grid.dim2.lat, ts, lat); % interpolate to standard grid
        ts = permute(ts, [2 1 3]); % reorder to original
        ts_vert = repmat(ts, [1 1 1 size(temp, 3)]); % dims (lon x lat x time x plev)
        ts_vert = permute(ts_vert, [1 2 4 3]); % dims (lon x lat x plev x time)
    elseif strcmp(type, 'echam')
        ps = permute(srfc.aps, [2 1 3]); % bring lat to front to interpolate
        ps = interp1(grid.dim2.lat, ps, lat, 'linear', 'extrap'); % interpolate to standard grid
        ps = permute(ps, [2 1 3]); % reorder to original
        ps_vert = repmat(ps, [1 1 1 size(temp, 3)]); % dims (lon x lat x time x plev)
        ps_vert = permute(ps_vert, [1 2 4 3]); % dims (lon x lat x plev x time)
        pa = permute(repmat(grid.dim3.plev, [1 size(ps)]), [2 3 1 4]);
        ts = permute(srfc.temp2, [2 1 3]); % bring lat to front to interpolate
        ts = interp1(grid.dim2.lat, ts, lat); % interpolate to standard grid
        ts = permute(ts, [2 1 3]); % reorder to original
        ts_vert = repmat(ts, [1 1 1 size(temp, 3)]); % dims (lon x lat x time x plev)
        ts_vert = permute(ts_vert, [1 2 4 3]); % dims (lon x lat x plev x time)
    end
    surface_mask = nan(size(temp));
    surface_mask(pa < ps_vert) = 1;

    temp_sm.lo = temp.*surface_mask; % filter temp with surface mask

    % add tsurf data and interpolate to higher resolution vertical grid
    [pa_plus ta_plus] = deal(nan([size(pa,1), size(pa,2) size(pa,3)+1 size(pa,4)])); % create empty grid with one extra vertical level
    pa_plus(:,:,1:end-1,:) = pa; % populate with standard pressure grid
    ta_plus(:,:,1:end-1,:) = temp_sm.lo; % populate with standard atmospheric temperature
    pa_plus(:,:,end,:) = ps_vert(:,:,1,:); % add surface pressure data into standard pressure grid
    ta_plus(:,:,end,:) = ts_vert(:,:,1,:); % add surface temperature data
    pa_plus = permute(pa_plus, [3 1 2 4]); % bring plev dimension to front
    ta_plus = permute(ta_plus, [3 1 2 4]); % bring plev dimension to front
    [pa_plus sort_index] = sort(pa_plus, 1, 'descend'); % sort added surface pressure such that pressure decreases monotonically
    tai_sm.lo = nan(length(par.pa), size(pa, 1), size(pa, 2), size(pa, 4));
    pb = CmdLineProgressBar("Sorting and interpolating temperature to new standard grid...");
    for ilon=1:size(pa_plus,2)
        pb.print(ilon, size(pa_plus,2));
        for ilat=1:size(pa_plus,3)
            for time=1:size(pa_plus,4)
                ta_plus(:,ilon,ilat,time) = ta_plus(sort_index(:,ilon,ilat,time),ilon,ilat,time); % sort temperature (has to be in loop because sort_index works for vector calls only)
                tai_sm.lo(:,ilon,ilat,time) = interp1(pa_plus(:,ilon,ilat,time), ta_plus(:,ilon,ilat,time), par.pa, 'linear'); % interpolate to higher resolution vertical grid
            end
        end
    end
    clear pa_plus ta_plus; % clear unneeded variables
    tai_sm.lo = permute(tai_sm.lo, [2 3 1 4]); % bring plev back to third dimension

    % apply masks to surface pressure and RCAE regimes
    ps_n.lo = ps;
    ps_n.l = ps .* mask.ocean;
    ps_n.o = ps .* mask.land;

    mask.land_vert = repmat(mask.land, [1 1 1 size(temp, 3)]); % expand land mask to vertical dim
    mask.land_vert = permute(mask.land, [1 2 4 3]); % place vertical dim where it belongs
    mask.ocean_vert = repmat(mask.ocean, [1 1 1 size(temp, 3)]); % expand ocean mask to vertical dim
    mask.ocean_vert = permute(mask.ocean, [1 2 4 3]); % place vertical dim where it belongs

    temp_sm.l = temp.*surface_mask.*mask.ocean_vert; % filter temp with surface mask
    temp_sm.o = temp.*surface_mask.*mask.land_vert; % filter temp with surface mask
    temp_sm.lo = permute(temp_sm.lo, [1 2 4 3]); % bring plev to last dimension
    temp_sm.l = permute(temp_sm.l, [1 2 4 3]); % bring plev to last dimension
    temp_sm.o = permute(temp_sm.o, [1 2 4 3]); % bring plev to last dimension

    tai_sm.l = tai_sm.lo.*surface_mask.*mask.ocean_vert; % filter tai with surface mask
    tai_sm.o = tai_sm.lo.*surface_mask.*mask.land_vert; % filter tai with surface mask
    tai_sm.lo = permute(tai_sm.lo, [1 2 4 3]); % bring plev to last dimension
    tai_sm.l = permute(tai_sm.l, [1 2 4 3]); % bring plev to last dimension
    tai_sm.o = permute(tai_sm.o, [1 2 4 3]); % bring plev to last dimension

    pa_sm.lo = pa.*surface_mask; % filter pressure grid with surface mask
    pa_sm.l = pa.*surface_mask.*mask.ocean_vert; % filter pa with surface mask
    pa_sm.o = pa.*surface_mask.*mask.land_vert; % filter pa with surface mask

    % for l = {'lo', 'l', 'o'}; land = l{1};
    %     temp_smz.(land) = squeeze(nanmean(temp_sm.(land), 1)); % take zonal average
    %     temp_smz.(land) = permute(temp_smz.(land), [1 3 2]); % put plev at the last dimension
    %     ps_z.(land) = squeeze(nanmean(ps_n.(land), 1));
    % end

    mask_t.land = nanmean(mask.land, 3);
    mask_t.ocean = nanmean(mask.ocean, 3);
    f_vec = assign_fw(type, par);
    for f = f_vec; fw = f{1};
        for c = fieldnames(rcae_t.lo.ann.(fw))'; crit = c{1};
            for l = {'lo', 'l', 'o'}; land = l{1}; % over land, over ocean, or both
                for t = {'ann', 'djf', 'jja', 'mam', 'son'}; time = t{1};
                    for re = {'rce', 'rae'}; regime = re{1};
                        if strcmp(time, 'ann')
                            temp_n.(land).(time) = squeeze(nanmean(temp_sm.(land), 3));
                            tai_n.(land).(time) = squeeze(nanmean(tai_sm.(land), 3));
                        elseif strcmp(time, 'djf')
                            temp_shift = circshift(temp_sm.(land), 1, 3);
                            tai_shift = circshift(tai_sm.(land), 1, 3);
                            temp_n.(land).(time) = squeeze(nanmean(temp_shift(:,:,1:3,:), 3));
                            tai_n.(land).(time) = squeeze(nanmean(tai_shift(:,:,1:3,:), 3));
                        elseif strcmp(time, 'jja')
                            temp_n.(land).(time) = squeeze(nanmean(temp_sm.(land)(:,:,6:8,:), 3));
                            tai_n.(land).(time) = squeeze(nanmean(tai_sm.(land)(:,:,6:8,:), 3));
                        elseif strcmp(time, 'mam')
                            temp_n.(land).(time) = squeeze(nanmean(temp_sm.(land)(:,:,3:5,:), 3));
                            tai_n.(land).(time) = squeeze(nanmean(tai_sm.(land)(:,:,3:5,:), 3));
                        elseif strcmp(time, 'son')
                            temp_n.(land).(time) = squeeze(nanmean(temp_sm.(land)(:,:,9:11,:), 3));
                            tai_n.(land).(time) = squeeze(nanmean(tai_sm.(land)(:,:,9:11,:), 3));
                        end

                        filt.(land).(time).(fw).(crit).(regime) = nan(size(rcae_t.(land).(time).(fw).(crit))); % create empty arrays to store filtering array
                        if strcmp(regime, 'rce'); filt.(land).(time).(fw).(crit).(regime)(rcae_t.(land).(time).(fw).(crit)==1)=1; % set RCE=1, elsewhere nan
                        elseif strcmp(regime, 'rae'); filt.(land).(time).(fw).(crit).(regime)(rcae_t.(land).(time).(fw).(crit)==-1)=1; % set RAE=1, elsewhere nan
                        end

                        temp_t.(land).(time).(fw).(crit).(regime) = temp_n.(land).(time) .* filt.(land).(time).(fw).(crit).(regime);
                        tai_t.(land).(time).(fw).(crit).(regime) = tai_n.(land).(time) .* filt.(land).(time).(fw).(crit).(regime);

                        nanfilt.(regime) = nan(size(temp_t.(land).(time).(fw).(crit).(regime)));
                        nanfilt.(regime)(~isnan(temp_t.(land).(time).(fw).(crit).(regime))) = 1;

                        % take cosine-weighted meridional average
                        for d = {'all', 'nh', 'sh', 'tp'}; domain = d{1};
                            if strcmp(regime, 'rce')
                                if strcmp(domain, 'all')
                                    cosw = repmat(cosd(lat)', [size(nanfilt.(regime), 1) 1]);
                                    denm = temp_t.(land).(time).(fw).(crit).(regime);
                                    denmi = tai_t.(land).(time).(fw).(crit).(regime);
                                    nume = nanfilt.(regime);
                                elseif strcmp(domain, 'nh')
                                    cosw = repmat(cosd(lat(lat>30))', [size(nanfilt.(regime), 1) 1]);
                                    denm = temp_t.(land).(time).(fw).(crit).(regime)(:,lat>30,:);
                                    denmi = tai_t.(land).(time).(fw).(crit).(regime)(:,lat>30,:);
                                    nume = nanfilt.(regime)(:,lat>30,:);
                                elseif strcmp(domain, 'sh')
                                    cosw = repmat(cosd(lat(lat<-30))', [size(nanfilt.(regime), 1) 1]);
                                    denm = temp_t.(land).(time).(fw).(crit).(regime)(:,lat<-30,:);
                                    denmi = tai_t.(land).(time).(fw).(crit).(regime)(:,lat<-30,:);
                                    nume = nanfilt.(regime)(:,lat<-30,:);
                                elseif strcmp(domain, 'tp')
                                    cosw = repmat(cosd(lat(abs(lat)<10))', [size(nanfilt.(regime), 1) 1]);
                                    denm = temp_t.(land).(time).(fw).(crit).(regime)(:,abs(lat)<10,:);
                                    denmi = tai_t.(land).(time).(fw).(crit).(regime)(:,abs(lat)<10,:);
                                    nume = nanfilt.(regime)(:,abs(lat)<10,:);
                                end
                            elseif strcmp(regime, 'rae')
                                if strcmp(domain, 'all')
                                    cosw = repmat(cosd(lat)', [size(nanfilt.(regime), 1) 1]);
                                    denm = temp_t.(land).(time).(fw).(crit).(regime);
                                    denmi = tai_t.(land).(time).(fw).(crit).(regime);
                                    nume = nanfilt.(regime);
                                elseif strcmp(domain, 'nh')
                                    cosw = repmat(cosd(lat(lat>0))', [size(nanfilt.(regime), 1) 1]);
                                    denm = temp_t.(land).(time).(fw).(crit).(regime)(:,lat>0,:);
                                    denmi = tai_t.(land).(time).(fw).(crit).(regime)(:,lat>0,:);
                                    nume = nanfilt.(regime)(:,lat>0,:);
                                elseif strcmp(domain, 'sh')
                                    cosw = repmat(cosd(lat(lat<0))', [size(nanfilt.(regime), 1) 1]);
                                    denm = temp_t.(land).(time).(fw).(crit).(regime)(:,lat<0,:);
                                    denmi = tai_t.(land).(time).(fw).(crit).(regime)(:,lat<0,:);
                                    nume = nanfilt.(regime)(:,lat<0,:);
                                end
                            end

                            ta.(regime).(domain).(fw).(crit).(land).(time) = nansum(cosw.*denm, 2) ./ nansum(cosw.*nume, 2); % area-weighted meridional average
                            ta.(regime).(domain).(fw).(crit).(land).(time) = squeeze(nanmean(ta.(regime).(domain).(fw).(crit).(land).(time), 1)); % zonal average
                            size(denmi)
                            size(nume)
                            tai.(regime).(domain).(fw).(crit).(land).(time) = nansum(cosw.*denmi, 2) ./ nansum(cosw.*nume, 2); % area-weighted meridional average
                            tai.(regime).(domain).(fw).(crit).(land).(time) = squeeze(nanmean(tai.(regime).(domain).(fw).(crit).(land).(time), 1)); % zonal average
                        end
                    end
                end
            end
        end
    end

    % save filtered data
    printname = [foldername 'ta'];
    if ~exist(foldername, 'dir')
        mkdir(foldername)
    end
    save(printname, 'ta', 'tai');
end % process temperature in RCE/RAE regimes
