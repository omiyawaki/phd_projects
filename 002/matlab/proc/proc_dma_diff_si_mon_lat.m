function proc_dma_diff_si_mon_lat(type, par)

    prefix = make_prefix(type, par);
    prefix_proc = make_prefix_proc(type, par);
    foldername = make_savedir_si_bl(type, par);
    ta_orig = load_temp(type, par);

    load(sprintf('%s/grid.mat', prefix)); % read grid data
    % load(sprintf('%s/dtdzzsi.mat', prefix)); % read temp in si coordinates
    load(sprintf('%s/dta_si.mat', prefix)); % read temp in si coordinates
    % load(sprintf('%s/malrzsi.mat', prefix)); % read temp in si coordinates
    load(sprintf('%s/dma_si.mat', prefix)); % read temp in si coordinates
    % load(sprintf('%s/masks.mat', prefix_proc)); % load land and ocean masks

    if contains(par.model, 'GISS-E2')
        dtdzzsi = permute(dtdzzsi, [2 1 3 4]);
        dtdzzsi = interp1(grid.dim3.lat_zg, dtdzzsi, grid.dim3.lat);
        dtdzzsi = permute(dtdzzsi, [2 1 3 4]);
    end

    if strcmp(par.lat_interp, 'std')
        lat = par.lat_std;
    else
        lat = grid.dim3.lat;
    end

    dt = (dma_si - dta_si.spl);
    dt = permute(dt,[2 1 3 4]); % bring lat to 1st
    dt = interp1(grid.dim3.lat, dt, lat); % interpolate to standard lat
    dt = permute(dt, [3 2 1 4]); % bring height front
    dt = interp1(grid.dim3.si, dt, linspace(par.si_bl,par.si_up,100)); % prepare to average between 1000-200 hPa
    dt = squeeze(nanmean(dt,1)); % take vertical average

    dma_si = permute(dma_si,[2 1 3 4]); % bring lat to 1st
    dma_si = interp1(grid.dim3.lat, dma_si, lat); % interpolate to standard lat
    dma_si = permute(dma_si, [3 2 1 4]); % bring height front
    dma_si = interp1(grid.dim3.si, dma_si, linspace(par.si_bl,par.si_up,100)); % prepare to average between 1000-200 hPa
    dma_si = squeeze(nanmean(dma_si,1)); % take vertical average

    dt_norm_orig = dt./dma_si;

    dt_norm0.lo = dt_norm_orig;
    % dt_norm0.l = dt_norm0.lo.*mask.ocean; % filter dt_norm0 with surface mask
    % dt_norm0.o = dt_norm0.lo.*mask.land; % filter dt_norm0 with surface mask

    for l = par.land_list; land = l{1}; % over land, over ocean, or both
        dt_norm.(land)= squeeze(nanmean(dt_norm0.(land), 1)); % zonal average
        % take time averages
        for t = {'ann', 'djf', 'jja', 'mam', 'son'}; time = t{1};
            if strcmp(time, 'ann')
                dt_norm_t.(land).(time) = squeeze(nanmean(dt_norm0.(land), 3));
            elseif strcmp(time, 'djf')
                dt_norm_shift.(land) = circshift(dt_norm0.(land), 1, 3);
                dt_norm_t.(land).(time) = squeeze(nanmean(dt_norm_shift.(land)(:,:,1:3,:), 3));
            elseif strcmp(time, 'jja')
                dt_norm_t.(land).(time) = squeeze(nanmean(dt_norm0.(land)(:,:,6:8,:), 3));
            elseif strcmp(time, 'mam')
                dt_norm_t.(land).(time) = squeeze(nanmean(dt_norm0.(land)(:,:,3:5,:), 3));
            elseif strcmp(time, 'son')
                dt_norm_t.(land).(time) = squeeze(nanmean(dt_norm0.(land)(:,:,9:11,:), 3));
            end
            dt_norm_zt.(land).(time) = squeeze(nanmean(dt_norm_t.(land).(time), 1));
        end
    end

    % save filtered data
    printname = sprintf('%sdt_norm_si_mon_lat_%g.mat', foldername, par.si_up);
    if ~exist(foldername, 'dir')
        mkdir(foldername)
    end
    save(printname, 'dt_norm', 'lat');

    printname = sprintf('%sdt_norm_si_lon_lat_%g.mat', foldername, par.si_up);
    save(printname, 'dt_norm_t', 'lat');

    printname = sprintf('%sdt_norm_si_lat_%g.mat', foldername, par.si_up);
    save(printname, 'dt_norm_zt', 'lat');

end % compute mon x lat gamma percentage difference field with land/ocean masking
