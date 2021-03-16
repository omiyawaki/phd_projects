function make_tai(type, par)
    % add surface data to temperature and interpolate to hi-res grid

    prefix = make_prefix(type, par);
    prefix_proc = make_prefix(type, par);
    foldername = make_savedir(type, par);
    temp = load_temp(type, par);
    zg = load_zg(type, par);

    load(sprintf('%s/grid.mat', prefix)); % read grid data
    load(sprintf('%s/srfc.mat', prefix)); % read surface variable data

    if strcmp(type, 'gcm') & any(contains(par.model, {'GISS-E2-H', 'GISS-E2-R'})) % zg data in GISS-E2-H has an anomalous lat grid
        var = 'zg';
        file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/gcm/%s/%s_Amon_%s_%s_r1i1p1_*.ymonmean.nc', par.model, var, par.model, par.gcm.clim));
        fullpath=sprintf('%s/%s', file.folder, file.name);
        zg_lat = double(ncread(fullpath, 'lat'));
        zg = permute(zg, [2 1 3 4]);
        zg = interp1(zg_lat, zg, grid.dim3.lat);
        zg = permute(zg, [2 1 3 4]);
    end

    % create surface mask
    if strcmp(type, 'era5') | strcmp(type, 'erai') | strcmp(type, 'era5c')
        ps_vert = repmat(srfc.sp, [1 1 1 size(temp, 3)]); % dims (lon x lat x time x plev)
        ps_vert = permute(ps_vert, [1 2 4 3]); % dims (lon x lat x plev x time)
        pa = double(permute(repmat(grid.dim3.plev, [1 size(srfc.sp)]), [2 3 1 4]));
        ts_vert = repmat(srfc.t2m, [1 1 1 size(temp, 3)]); % dims (lon x lat x time x plev)
        ts_vert = permute(ts_vert, [1 2 4 3]); % dims (lon x lat x plev x time)
        zs_vert = repmat(srfc.zs, [1 1 1 size(temp, 3)]); % dims (lon x lat x time x plev)
        zs_vert = permute(zs_vert, [1 2 4 3]); % dims (lon x lat x plev x time)
    elseif strcmp(type, 'merra2')
        ps_vert = repmat(srfc.PS, [1 1 1 size(temp, 3)]); % dims (lon x lat x time x plev)
        ps_vert = permute(ps_vert, [1 2 4 3]); % dims (lon x lat x plev x time)
        pa = double(permute(repmat(grid.dim3.plev, [1 size(srfc.PS)]), [2 3 1 4]));
        ts_vert = repmat(srfc.T2M, [1 1 1 size(temp, 3)]); % dims (lon x lat x time x plev)
        ts_vert = permute(ts_vert, [1 2 4 3]); % dims (lon x lat x plev x time)
        zs_vert = repmat(srfc.zs, [1 1 1 size(temp, 3)]); % dims (lon x lat x time x plev)
        zs_vert = permute(zs_vert, [1 2 4 3]); % dims (lon x lat x plev x time)
    elseif any(strcmp(type, {'gcm', 'jra55'}))
        ps_vert = repmat(srfc.ps, [1 1 1 size(temp, 3)]); % dims (lon x lat x time x plev)
        ps_vert = permute(ps_vert, [1 2 4 3]); % dims (lon x lat x plev x time)
        pa = permute(repmat(grid.dim3.plev, [1 size(srfc.ps)]), [2 3 1 4]);
        ts_vert = repmat(srfc.tas, [1 1 1 size(temp, 3)]); % dims (lon x lat x time x plev)
        ts_vert = permute(ts_vert, [1 2 4 3]); % dims (lon x lat x plev x time)
        zs_vert = repmat(srfc.zs, [1 1 1 size(temp, 3)]); % dims (lon x lat x time x plev)
        zs_vert = permute(zs_vert, [1 2 4 3]); % dims (lon x lat x plev x time)
    elseif any(strcmp(type, {'echam', 'echam_pl'}))
        ps_vert = repmat(srfc.aps, [1 1 1 size(temp, 3)]); % dims (lon x lat x time x plev)
        ps_vert = permute(ps_vert, [1 2 4 3]); % dims (lon x lat x plev x time)
        pa = permute(repmat(grid.dim3.plev, [1 size(srfc.aps)]), [2 3 1 4]);
        ts_vert = repmat(srfc.temp2, [1 1 1 size(temp, 3)]); % dims (lon x lat x time x plev)
        ts_vert = permute(ts_vert, [1 2 4 3]); % dims (lon x lat x plev x time)
        zs_vert = repmat(srfc.zs, [1 1 1 size(temp, 3)]); % dims (lon x lat x time x plev)
        zs_vert = permute(zs_vert, [1 2 4 3]); % dims (lon x lat x plev x time)
    elseif strcmp(type, 'echam_ml')
        ps_vert = repmat(srfc.aps, [1 1 1 size(temp, 3)]); % dims (lon x lat x time x plev)
        ps_vert = permute(ps_vert, [1 2 4 3]); % dims (lon x lat x plev x time)
        a = permute(repmat(grid.dim3.a, [1 size(srfc.aps)]), [2 3 1 4]);
        b = permute(repmat(grid.dim3.b, [1 size(srfc.aps)]), [2 3 1 4]);
        pa = a+b.*ps_vert;
        ts_vert = repmat(srfc.temp2, [1 1 1 size(temp, 3)]); % dims (lon x lat x time x plev)
        ts_vert = permute(ts_vert, [1 2 4 3]); % dims (lon x lat x plev x time)
        zs_vert = repmat(srfc.zs, [1 1 1 size(temp, 3)]); % dims (lon x lat x time x plev)
        zs_vert = permute(zs_vert, [1 2 4 3]); % dims (lon x lat x plev x time)
    end
    surface_mask = nan(size(temp));
    surface_mask(pa < ps_vert) = 1;

    temp_sm.lo = temp.*surface_mask; % filter temp with surface mask
    zg_sm.lo = zg.*surface_mask; % filter zg with surface mask

    % add tsurf data and interpolate to higher resolution vertical grid
    [pa_plus ta_plus zg_plus] = deal(nan([size(pa,1), size(pa,2) size(pa,3)+1 size(pa,4)])); % create empty grid with one extra vertical level
    pa_plus(:,:,1:end-1,:) = pa; % populate with standard pressure grid
    ta_plus(:,:,1:end-1,:) = temp_sm.lo; % populate with standard atmospheric temperature
    zg_plus(:,:,1:end-1,:) = zg_sm.lo; % populate with szgndard atmospheric temperature
    pa_plus(:,:,end,:) = ps_vert(:,:,1,:); % add surface pressure data into standard pressure grid
    ta_plus(:,:,end,:) = ts_vert(:,:,1,:); % add surface temperature data
    zg_plus(:,:,end,:) = zs_vert(:,:,1,:); % add surface zg data
    pa_plus = permute(pa_plus, [3 1 2 4]); % bring plev dimension to front
    ta_plus = permute(ta_plus, [3 1 2 4]); % bring plev dimension to front
    zg_plus = permute(zg_plus, [3 1 2 4]); % bring plev dimension to front
    [pa_plus sort_index] = sort(pa_plus, 1, 'descend'); % sort added surface pressure such that pressure decreases monotonically
    % [pa_plus_sorted sort_index] = sort(pa_plus, 1, 'descend'); % sort added surface pressure such that pressure decreases monotonically
    tai_sm.lo = nan(length(par.pa), size(pa, 1), size(pa, 2), size(pa, 4));
    zgi_sm.lo = nan(length(par.pa), size(pa, 1), size(pa, 2), size(pa, 4));
    pb = CmdLineProgressBar("Sorting and interpolating temperature to new standard grid...");
    for ilon=1:size(pa_plus,2)
        pb.print(ilon, size(pa_plus,2));
        if contains(type, 'era5')
            pb2 = CmdLineProgressBar("Lat...");
        end
        for ilat=1:size(pa_plus,3)
            if contains(type, 'era5')
                pb2.print(ilat, size(pa_plus,3));
            end
            for time=1:size(pa_plus,4)
                warning off; 
                % ta_plus(:,ilon,ilat,time) = ta_plus(sort_index(:,ilon,ilat,time),ilon,ilat,time); % sort temperature (has to be in loop because sort_index works for vector calls only)
                % zg_plus(:,ilon,ilat,time) = zg_plus(sort_index(:,ilon,ilat,time),ilon,ilat,time); % sort temperature (has to be in loop because sort_index works for vector calls only)
                pa_plus_col = pa_plus(:,ilon,ilat,time);
                ta_plus_col = ta_plus(sort_index(:,ilon,ilat,time),ilon,ilat,time); % sort temperature (has to be in loop because sort_index works for vector calls only)
                zg_plus_col = zg_plus(sort_index(:,ilon,ilat,time),ilon,ilat,time); % sort zg (has to be in loop because sort_index works for vector calls only)
                [~,iuq] = unique(pa_plus_col);

                if all(isnan(ta_plus_col)) | sum(isnan(ta_plus_col)) == length(ta_plus_col)-1
                    tai_sm.lo(:,ilon,ilat,time) = nan(size(par.pa'));
                else
                    % tai_sm.lo(:,ilon,ilat,time) = interp1(pa_plus_col(iuq), ta_plus_col(iuq), par.pa, 'spline', nan); % interpolate to higher resolution vertical grid
                    tai_sm.lo(:,ilon,ilat,time) = interp1(pa_plus_col(iuq), ta_plus_col(iuq), par.pa, 'pchip'); % interpolate to higher resolution vertical grid
                    tai_sm.lo(par.pa>ps_vert(ilon,ilat,1,time),ilon,ilat,time) = nan; % mask out data below surface pressure
                end

                if all(isnan(zg_plus_col)) | sum(isnan(zg_plus_col)) == length(ta_plus_col)-1
                    zgi_sm.lo(:,ilon,ilat,time) = nan(size(par.pa'));
                else
                    % zgi_sm.lo(:,ilon,ilat,time) = interp1(pa_plus_col(iuq), zg_plus_col(iuq), par.pa, 'spline', nan); % interpolate to higher resolution vertical grid
                    zgi_sm.lo(:,ilon,ilat,time) = interp1(pa_plus_col(iuq), zg_plus_col(iuq), par.pa, 'pchip'); % interpolate to higher resolution vertical grid
                    zgi_sm.lo(par.pa>ps_vert(ilon,ilat,1,time),ilon,ilat,time) = nan; % mask out data below surface pressure
                end

                % disp('Checkpoint 1')

                % if all(isnan(ta_plus)) | sum(isnan(ta_plus)) == 1
                %     tai_sm.lo(:,ilon,ilat,time) = nan(size(par.pa));
                % end

                % if all(isnan(zg_plus)) | sum(isnan(zg_plus)) == 1
                %     zg_sm.lo(:,ilon,ilat,time) = nan(size(par.pa));
                % end

                % disp('Checkpoint 2')

                % if sum(isnan(ta_plus))>1 & sum(isnan(zg_plus))>1

                %     tapanan = squeeze(pa_plus(:,ilon,ilat,time));
                %     tananfi = isnan(squeeze(pa_plus(:,ilon,ilat,time))) | isnan(squeeze(ta_plus(:,ilon,ilat,time)));
                %     tanan = squeeze(ta_plus(:,ilon,ilat,time));
                %     zgpanan = squeeze(pa_plus(:,ilon,ilat,time));
                %     zgnanfi = isnan(squeeze(pa_plus(:,ilon,ilat,time))) | isnan(squeeze(zg_plus(:,ilon,ilat,time)));
                %     zgnan = squeeze(zg_plus(:,ilon,ilat,time));

                %     tapanan(tananfi) = [];
                %     tanan(tananfi) = [];
                %     zgpanan(zgnanfi) = [];
                %     zgnan(zgnanfi) = [];

                %     disp('Checkpoint 3')

                %     % tai_sm.lo(:,ilon,ilat,time) = interp1(tapanan, tanan, par.pa, 'spline', nan); % interpolate to higher resolution vertical grid
                %     % zgi_sm.lo(:,ilon,ilat,time) = interp1(zgpanan, zgnan, par.pa, 'spline', nan); % interpolate to higher resolution vertical grid
                %     tai_sm.lo(:,ilon,ilat,time) = interp1(tapanan, tanan, par.pa, 'linear'); % interpolate to higher resolution vertical grid
                %     zgi_sm.lo(:,ilon,ilat,time) = interp1(zgpanan, zgnan, par.pa, 'linear'); % interpolate to higher resolution vertical grid
                %     disp('Checkpoint 4')
                % end

            end
        end
    end
    % tai_sm.lo = interp1(pa_plus, tanan, par.pa, 'spline', nan); % interpolate to higher resolution vertical grid
    % zgi_sm.lo = interp1(zgpanan, zgnan, par.pa, 'spline', nan); % interpolate to higher resolution vertical grid

    clear pa_plus ta_plus zg_plus; % clear unneeded variables
    tai_sm.lo = permute(tai_sm.lo, [2 3 1 4]); % bring plev back to third dimension
    tai = tai_sm.lo;
    zgi_sm.lo = permute(zgi_sm.lo, [2 3 1 4]); % bring plev back to third dimension
    zgi = zgi_sm.lo;

    filename='tai_simp.mat';
    save(sprintf('%s/%s', foldername, filename), 'tai', 'zgi', '-v7.3');

end
