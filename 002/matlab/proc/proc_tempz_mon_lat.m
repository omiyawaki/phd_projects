function proc_tempz_mon_lat(type, par)
    if strcmp(type, 'era5') | strcmp(type, 'erai') | strcmp(type, 'era5c')
        foldername = sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s', type, par.(type).yr_span, par.lat_interp);
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s', type, par.(type).yr_span);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s', type, par.(type).yr_span);
    elseif strcmp(type, 'gcm')
        foldername = sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s/%s/', type, par.model, par.gcm.clim, par.lat_interp);
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s/%s', type, par.model, par.gcm.clim);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s', type, par.model, par.gcm.clim);
    elseif strcmp(type, 'echam')
        foldername = sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s/', type, par.echam.clim, par.lat_interp);
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s', type, par.echam.clim);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s', type, par.echam.clim);
    elseif strcmp(type, 'echam_ml')
        foldername = sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s', type, par.(type).yr_span, par.lat_interp);
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s', type, par.(type).yr_span);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s', type, par.(type).yr_span);
    elseif strcmp(type, 'echam_pl')
        foldername = sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s', type, par.(type).yr_span, par.lat_interp);
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s', type, par.(type).yr_span);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s', type, par.(type).yr_span);
    end
    load(sprintf('%s/grid.mat', prefix)); % read grid data
    load(sprintf('%s/tempz.mat', prefix)); taz_orig = tempz; clear tempz; % read temp in z coordinates
    load(sprintf('%s/pz.mat', prefix)); paz_orig = pz; clear pz; % read pa in z coordinates
    load(sprintf('%s/srfc.mat', prefix)); % load surface data
    load(sprintf('%s/%s/masks.mat', prefix_proc, par.lat_interp)); % load land and ocean masks

    if strcmp(par.lat_interp, 'std')
        lat = par.lat_std;
    else
        lat = grid.dim3.lat;
    end

    % interpolate ta to standard lat grid
    taz_orig = permute(taz_orig, [2 1 3 4]);
    taz_orig = interp1(grid.dim3.lat, taz_orig, lat);
    taz_orig = permute(taz_orig, [2 1 3 4]);

    paz_orig = permute(paz_orig, [2 1 3 4]);
    paz_orig = interp1(grid.dim3.lat, paz_orig, lat);
    paz_orig = permute(paz_orig, [2 1 3 4]);

    % create surface mask
    if strcmp(type, 'era5') | strcmp(type, 'erai') | strcmp(type, 'era5c')
        zs = permute(srfc.zs, [2 1 3]); % bring lat to front to interpolate
        zs = interp1(grid.dim2.lat, zs, lat); % interpolate to standard grid
        zs = permute(zs, [2 1 3]); % reorder to original
        zs_vert = repmat(zs, [1 1 1 size(taz_orig, 3)]); % dims (lon x lat x time x plev)
        zs_vert = permute(zs_vert, [1 2 4 3]); % dims (lon x lat x plev x time)
        za = permute(repmat(grid.dim3.z, [1 size(zs)]), [2 3 1 4]);
        ts = permute(srfc.t2m, [2 1 3]); % bring lat to front to interpolate
        ts = interp1(grid.dim2.lat, ts, lat); % interpolate to standard grid
        ts = permute(ts, [2 1 3]); % reorder to original
        ts_vert = repmat(ts, [1 1 1 size(taz_orig, 3)]); % dims (lon x lat x time x plev)
        ts_vert = permute(ts_vert, [1 2 4 3]); % dims (lon x lat x plev x time)
    elseif strcmp(type, 'gcm')
        zs = permute(srfc.zs, [2 1 3]); % bring lat to front to interpolate
        zs = interp1(grid.dim2.lat, zs, lat); % interpolate to standard grid
        zs = permute(zs, [2 1 3]); % reorder to original
        zs_vert = repmat(zs, [1 1 1 size(taz_orig, 3)]); % dims (lon x lat x time x plev)
        zs_vert = permute(zs_vert, [1 2 4 3]); % dims (lon x lat x plev x time)
        za = permute(repmat(grid.dim3.z, [1 size(zs)]), [2 3 1 4]);
        ts = permute(srfc.tas, [2 1 3]); % bring lat to front to interpolate
        ts = interp1(grid.dim2.lat, ts, lat); % interpolate to standard grid
        ts = permute(ts, [2 1 3]); % reorder to original
        ts_vert = repmat(ts, [1 1 1 size(taz_orig, 3)]); % dims (lon x lat x time x plev)
        ts_vert = permute(ts_vert, [1 2 4 3]); % dims (lon x lat x plev x time)
    elseif strcmp(type, 'echam') | strcmp(type, 'echam_ml') | strcmp(type, 'echam_pl')
        zs = permute(srfc.zs, [2 1 3]); % bring lat to front to interpolate
        zs = interp1(grid.dim2.lat, zs, lat); % interpolate to standard grid
        zs = permute(zs, [2 1 3]); % reorder to original
        zs_vert = repmat(zs, [1 1 1 size(taz_orig, 3)]); % dims (lon x lat x time x plev)
        zs_vert = permute(zs_vert, [1 2 4 3]); % dims (lon x lat x plev x time)
        za = permute(repmat(grid.dim3.z, [1 size(zs)]), [2 3 1 4]);
        ts = permute(srfc.temp2, [2 1 3]); % bring lat to front to interpolate
        ts = interp1(grid.dim2.lat, ts, lat); % interpolate to standard grid
        ts = permute(ts, [2 1 3]); % reorder to original
        ts_vert = repmat(ts, [1 1 1 size(taz_orig, 3)]); % dims (lon x lat x time x plev)
        ts_vert = permute(ts_vert, [1 2 4 3]); % dims (lon x lat x plev x time)
    end

    if strcmp(type, 'echam_ml')
        smz = ones(size(taz_orig));
    else
        smz = nan(size(taz_orig));
        smz(za > zs_vert) = 1;
    end

    taz_sm.lo = taz_orig.*smz; % filter taz with surface mask
    paz_sm.lo = paz_orig.*smz; % filter paz with surface mask

    if par.do_surf
        % add tsurf data and interpolate to higher resolution vertical grid
        [pa_plus ta_plus] = deal(nan([size(pa,1), size(pa,2) size(pa,3)+1 size(pa,4)])); % create empty grid with one extra vertical level
        pa_plus(:,:,1:end-1,:) = pa; % populate with standard pressure grid
        ta_plus(:,:,1:end-1,:) = ta_sm.lo; % populate with standard atmospheric temperature
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
    end

    % Land/ocean filter 2D variables
    tsurf_sm.lo = ts;
    psurf_sm.lo = ps;
    zsurf_sm.lo = zs;
    tsurf_sm.l = ts.*mask.ocean; %filter out ocean
    psurf_sm.l = ps.*mask.ocean; %filter out ocean
    zsurf_sm.l = zs.*mask.ocean; %filter out ocean
    tsurf_sm.o = ts.*mask.land; %filter out land
    psurf_sm.o = ps.*mask.land; %filter out land
    zsurf_sm.o = zs.*mask.land; %filter out land

    % Land/ocean filter 3D variables
    mask.land_vert = repmat(mask.land, [1 1 1 size(ta_orig, 3)]); % expand land mask to vertical dim
    mask.land_vert = permute(mask.land, [1 2 4 3]); % place vertical dim where it belongs
    mask.ocean_vert = repmat(mask.ocean, [1 1 1 size(ta_orig, 3)]); % expand ocean mask to vertical dim
    mask.ocean_vert = permute(mask.ocean, [1 2 4 3]); % place vertical dim where it belongs
    mask.land_verti = repmat(mask.land, [1 1 1 length(par.pa)]); % expand land mask to vertiical dim
    mask.land_verti = permute(mask.land, [1 2 4 3]); % place vertiical dim where it belongs
    mask.ocean_verti = repmat(mask.ocean, [1 1 1 length(par.pa)]); % expand ocean mask to vertiical dim
    mask.ocean_verti = permute(mask.ocean, [1 2 4 3]); % place vertiical dim where it belongs

    ta_sm.l = ta_sm.lo.*mask.ocean_vert; % filter ta with surface mask
    ta_sm.o = ta_sm.lo.*mask.land_vert; % filter ta with surface mask
    ta_sm.lo = permute(ta_sm.lo, [1 2 4 3]); % bring plev to last dimension
    ta_sm.l = permute(ta_sm.l, [1 2 4 3]); % bring plev to last dimension
    ta_sm.o = permute(ta_sm.o, [1 2 4 3]); % bring plev to last dimension
    taz_sm.l = taz_sm.lo.*mask.ocean_vert; % filter taz with surface mask
    taz_sm.o = taz_sm.lo.*mask.land_vert; % filter taz with surface mask
    taz_sm.lo = permute(taz_sm.lo, [1 2 4 3]); % bring plev to last dimension
    taz_sm.l = permute(taz_sm.l, [1 2 4 3]); % bring plev to last dimension
    taz_sm.o = permute(taz_sm.o, [1 2 4 3]); % bring plev to last dimension
    tasi_sm.l = tasi_sm.lo.*mask.ocean_vert; % filter tasi with surface mask
    tasi_sm.o = tasi_sm.lo.*mask.land_vert; % filter tasi with surface mask
    tasi_sm.lo = permute(tasi_sm.lo, [1 2 4 3]); % bring plev to last dimension
    tasi_sm.l = permute(tasi_sm.l, [1 2 4 3]); % bring plev to last dimension
    tasi_sm.o = permute(tasi_sm.o, [1 2 4 3]); % bring plev to last dimension

    paz_sm.l = paz_sm.lo.*mask.ocean_vert; % filter paz with surface mask
    paz_sm.o = paz_sm.lo.*mask.land_vert; % filter paz with surface mask
    paz_sm.lo = permute(paz_sm.lo, [1 2 4 3]); % bring plev to last dimension
    paz_sm.l = permute(paz_sm.l, [1 2 4 3]); % bring plev to last dimension
    paz_sm.o = permute(paz_sm.o, [1 2 4 3]); % bring plev to last dimension
    pasi_sm.l = pasi_sm.lo.*mask.ocean_vert; % filter pasi with surface mask
    pasi_sm.o = pasi_sm.lo.*mask.land_vert; % filter pasi with surface mask
    pasi_sm.lo = permute(pasi_sm.lo, [1 2 4 3]); % bring plev to last dimension
    pasi_sm.l = permute(pasi_sm.l, [1 2 4 3]); % bring plev to last dimension
    pasi_sm.o = permute(pasi_sm.o, [1 2 4 3]); % bring plev to last dimension
    if par.do_surf
        tai_sm.l = tai_sm.lo.*mask.ocean_verti; % filter tai with surface mask
        tai_sm.o = tai_sm.lo.*mask.land_verti; % filter tai with surface mask
        tai_sm.lo = permute(tai_sm.lo, [1 2 4 3]); % bring plev to last dimension
        tai_sm.l = permute(tai_sm.l, [1 2 4 3]); % bring plev to last dimension
        tai_sm.o = permute(tai_sm.o, [1 2 4 3]); % bring plev to last dimension
    end

    mask_t.land = nanmean(mask.land, 3);
    mask_t.ocean = nanmean(mask.ocean, 3);

    for l = {'lo', 'l', 'o'}; land = l{1}; % over land, over ocean, or both
        ta.(land)= squeeze(mean(ta_sm.(land), 1)); % zonal average
        taz.(land)= squeeze(mean(taz_sm.(land), 1)); % zonal average
        tasi.(land)= squeeze(mean(tasi_sm.(land), 1)); % zonal average

        paz.(land)= squeeze(mean(paz_sm.(land), 1)); % zonal average
        pasi.(land)= squeeze(mean(pasi_sm.(land), 1)); % zonal average

        if par.do_surf; tai.(land)= squeeze(mean(tai_sm.(land), 1)); end % zonal average
        tsurf.(land)= squeeze(mean(tsurf_sm.(land), 1)); % zonal average
        psurf.(land)= squeeze(mean(psurf_sm.(land), 1)); % zonal average
        zsurf.(land)= squeeze(mean(zsurf_sm.(land), 1)); % zonal average
    end

    for l = {'lo', 'l', 'o'}; land = l{1}; % over land, over ocean, or both
        % take time averages
        for t = {'ann', 'djf', 'jja', 'mam', 'son'}; time = t{1};
            if strcmp(time, 'ann')
                ta_t.(land).(time) = squeeze(mean(ta_sm.(land), 3));
                taz_t.(land).(time) = squeeze(mean(taz_sm.(land), 3));
                tasi_t.(land).(time) = squeeze(mean(tasi_sm.(land), 3));

                paz_t.(land).(time) = squeeze(mean(paz_sm.(land), 3));
                pasi_t.(land).(time) = squeeze(mean(pasi_sm.(land), 3));

                if par.do_surf; tai_t.(land).(time) = squeeze(mean(tai_sm.(land), 3)); end
                tsurf_t.(land).(time) = squeeze(mean(tsurf_sm.(land), 3));
                psurf_t.(land).(time) = squeeze(mean(psurf_sm.(land), 3));
                zsurf_t.(land).(time) = squeeze(mean(zsurf_sm.(land), 3));
            elseif strcmp(time, 'djf')
                ta_shift.(land) = circshift(ta_sm.(land), 1, 3);
                taz_shift.(land) = circshift(taz_sm.(land), 1, 3);
                tasi_shift.(land) = circshift(tasi_sm.(land), 1, 3);

                paz_shift.(land) = circshift(paz_sm.(land), 1, 3);
                pasi_shift.(land) = circshift(pasi_sm.(land), 1, 3);

                if par.do_surf; tai_shift.(land) = circshift(tai_sm.(land), 1, 3); end
                tsurf_shift.(land) = circshift(tsurf_sm.(land), 1, 3);
                psurf_shift.(land) = circshift(psurf_sm.(land), 1, 3);
                zsurf_shift.(land) = circshift(zsurf_sm.(land), 1, 3);
                ta_t.(land).(time) = squeeze(mean(ta_shift.(land)(:,:,1:3,:), 3));
                taz_t.(land).(time) = squeeze(mean(taz_shift.(land)(:,:,1:3,:), 3));
                tasi_t.(land).(time) = squeeze(mean(tasi_shift.(land)(:,:,1:3,:), 3));

                paz_t.(land).(time) = squeeze(mean(paz_shift.(land)(:,:,1:3,:), 3));
                pasi_t.(land).(time) = squeeze(mean(pasi_shift.(land)(:,:,1:3,:), 3));

                if par.do_surf; tai_t.(land).(time) = squeeze(mean(tai_shift.(land)(:,:,1:3,:), 3)); end
                tsurf_t.(land).(time) = squeeze(mean(tsurf_shift.(land)(:,:,1:3,:), 3));
                psurf_t.(land).(time) = squeeze(mean(psurf_shift.(land)(:,:,1:3,:), 3));
                zsurf_t.(land).(time) = squeeze(mean(zsurf_shift.(land)(:,:,1:3,:), 3));
            elseif strcmp(time, 'jja')
                ta_t.(land).(time) = squeeze(mean(ta_sm.(land)(:,:,6:8,:), 3));
                taz_t.(land).(time) = squeeze(mean(taz_sm.(land)(:,:,6:8,:), 3));
                tasi_t.(land).(time) = squeeze(mean(tasi_sm.(land)(:,:,6:8,:), 3));

                paz_t.(land).(time) = squeeze(mean(paz_sm.(land)(:,:,6:8,:), 3));
                pasi_t.(land).(time) = squeeze(mean(pasi_sm.(land)(:,:,6:8,:), 3));

                if par.do_surf; tai_t.(land).(time) = squeeze(mean(tai_sm.(land)(:,:,6:8,:), 3)); end
                tsurf_t.(land).(time) = squeeze(mean(tsurf_sm.(land)(:,:,6:8,:), 3));
                psurf_t.(land).(time) = squeeze(mean(psurf_sm.(land)(:,:,6:8,:), 3));
                zsurf_t.(land).(time) = squeeze(mean(zsurf_sm.(land)(:,:,6:8,:), 3));
            elseif strcmp(time, 'mam')
                ta_t.(land).(time) = squeeze(mean(ta_sm.(land)(:,:,3:5,:), 3));
                taz_t.(land).(time) = squeeze(mean(taz_sm.(land)(:,:,3:5,:), 3));
                tasi_t.(land).(time) = squeeze(mean(tasi_sm.(land)(:,:,3:5,:), 3));

                paz_t.(land).(time) = squeeze(mean(paz_sm.(land)(:,:,3:5,:), 3));
                pasi_t.(land).(time) = squeeze(mean(pasi_sm.(land)(:,:,3:5,:), 3));

                if par.do_surf; tai_t.(land).(time) = squeeze(mean(tai_sm.(land)(:,:,3:5,:), 3)); end
                tsurf_t.(land).(time) = squeeze(mean(tsurf_sm.(land)(:,:,3:5,:), 3));
                psurf_t.(land).(time) = squeeze(mean(psurf_sm.(land)(:,:,3:5,:), 3));
                zsurf_t.(land).(time) = squeeze(mean(zsurf_sm.(land)(:,:,3:5,:), 3));
            elseif strcmp(time, 'son')
                ta_t.(land).(time) = squeeze(mean(ta_sm.(land)(:,:,9:11,:), 3));
                taz_t.(land).(time) = squeeze(mean(taz_sm.(land)(:,:,9:11,:), 3));
                tasi_t.(land).(time) = squeeze(mean(tasi_sm.(land)(:,:,9:11,:), 3));

                paz_t.(land).(time) = squeeze(mean(paz_sm.(land)(:,:,9:11,:), 3));
                pasi_t.(land).(time) = squeeze(mean(pasi_sm.(land)(:,:,9:11,:), 3));

                if par.do_surf; tai_t.(land).(time) = squeeze(mean(tai_sm.(land)(:,:,9:11,:), 3)); end
                tsurf_t.(land).(time) = squeeze(mean(tsurf_sm.(land)(:,:,9:11,:), 3));
                psurf_t.(land).(time) = squeeze(mean(psurf_sm.(land)(:,:,9:11,:), 3));
                zsurf_t.(land).(time) = squeeze(mean(zsurf_sm.(land)(:,:,9:11,:), 3));
            end
        end
    end

    % save filtered data
    printname = [foldername 'ta_mon_lat'];
    if ~exist(foldername, 'dir')
        mkdir(foldername)
    end
    if par.do_surf; save(printname, 'ta', 'taz', 'tasi', 'tai', 'lat', 'tsurf', 'psurf', 'zsurf');
    else save(printname, 'ta', 'taz', 'tasi', 'paz', 'pasi', 'lat', 'tsurf', 'psurf', 'zsurf', '-v7.3'); end

    printname = [foldername 'ta_lon_lat'];
    if ~exist(foldername, 'dir')
        mkdir(foldername)
    end
    if par.do_surf; save(printname, 'ta_t', 'taz_t', 'tasi_t', 'tai_t', 'lat', 'tsurf_t', 'psurf_t', 'zsurf_t');
    else save(printname, 'ta_t', 'taz_t', 'tasi_t', 'paz_t', 'pasi_t', 'lat', 'tsurf_t', 'psurf_t', 'zsurf_t', '-v7.3'); end
end % compute mon x lat temperature field
