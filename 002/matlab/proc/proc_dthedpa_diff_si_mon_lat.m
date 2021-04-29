function proc_dthedpa_si_mon_lat(type, par)

    prefix = make_prefix(type, par);
    prefix_proc = make_prefix_proc(type, par);
    foldername = make_savedir_si_bl(type, par);

    load(sprintf('%s/grid.mat', prefix)); % read grid data
    load(sprintf('%s/dthedpasi.mat', prefix)); % read temp in si coordinates
    % load(sprintf('%s/masks.mat', prefix_proc)); % load land and ocean masks

    if strcmp(par.lat_interp, 'std')
        lat = par.lat_std;
    else
        lat = grid.dim3.lat;
    end

    dthedpa_orig = permute(dthedpasi,[2 1 3 4]); % bring lat to 1st
    dthedpa_orig = interp1(grid.dim3.lat, dthedpa_orig, lat); % interpolate to standard lat
    dthedpa_orig = permute(dthedpa_orig, [3 2 1 4]); % bring height front
    dthedpa_orig = interp1(grid.dim3.si, dthedpa_orig, linspace(par.si_bl,par.si_up,100)); % prepare to average between 1000-200 hPa
    dthedpa_orig = squeeze(nanmean(dthedpa_orig,1)); % take vertical average

    dthedpa0.lo = dthedpa_orig;
    % dthedpa0.l = dthedpa0.lo.*mask.ocean; % filter dthedpa0 with surface mask
    % dthedpa0.o = dthedpa0.lo.*mask.land; % filter dthedpa0 with surface mask

    for l = par.land_list; land = l{1}; % over land, over ocean, or both
        dthedpa.(land)= squeeze(nanmean(dthedpa0.(land), 1)); % zonal average
        % take time averages
        for t = {'ann', 'djf', 'jja', 'mam', 'son'}; time = t{1};
            if strcmp(time, 'ann')
                dthedpa_t.(land).(time) = squeeze(nanmean(dthedpa0.(land), 3));
            elseif strcmp(time, 'djf')
                dthedpa_shift.(land) = circshift(dthedpa0.(land), 1, 3);
                dthedpa_t.(land).(time) = squeeze(nanmean(dthedpa_shift.(land)(:,:,1:3,:), 3));
            elseif strcmp(time, 'jja')
                dthedpa_t.(land).(time) = squeeze(nanmean(dthedpa0.(land)(:,:,6:8,:), 3));
            elseif strcmp(time, 'mam')
                dthedpa_t.(land).(time) = squeeze(nanmean(dthedpa0.(land)(:,:,3:5,:), 3));
            elseif strcmp(time, 'son')
                dthedpa_t.(land).(time) = squeeze(nanmean(dthedpa0.(land)(:,:,9:11,:), 3));
            end
        end
    end

    % save filtered data
    printname = sprintf('%sdthedpa_si_mon_lat_%g.mat', foldername, par.si_up);
    if ~exist(foldername, 'dir')
        mkdir(foldername)
    end
    save(printname, 'dthedpa', 'lat');

    printname = sprintf('%sdthedpa_si_lon_lat_%g.mat', foldername, par.si_up);
    if ~exist(foldername, 'dir')
        mkdir(foldername)
    end
    save(printname, 'dthedpa_t', 'lat');
end % compute mon x lat gamma percentage difference field with land/ocean masking
