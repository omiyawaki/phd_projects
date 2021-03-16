function proc_ma_mon_lat(type, par)

    prefix = make_prefix(type, par);
    foldername = make_savedir_proc(type, par);

    load(sprintf('%s/grid.mat', prefix)); % read grid data
    if strcmp(par.ma_init, 'surf')
        printname = [foldername 'ma_mon_lat_' par.ma_init '.mat'];
        printname2 = [foldername 'ma_lon_lat_' par.ma_init '.mat'];
        load(sprintf('%s/ma_si_%s.mat', prefix, par.ma_init)); masi_orig = ma_si; clear ma_si; % read temp in si coordinates
    else
        printname = [foldername 'ma_mon_lat_' num2str(par.ma_init) '.mat'];
        printname2 = [foldername 'ma_lon_lat_' num2str(par.ma_init) '.mat'];
        load(sprintf('%s/ma_si_%g.mat', prefix, par.ma_init)); masi_orig = ma_si; clear ma_si; % read temp in si coordinates
    end
    load(sprintf('%s/pa_si.mat', prefix)); pasi_orig = pa_si; clear pa_si; % read temp in si coordinates
    load(sprintf('%s/srfc.mat', prefix)); % load surface data
    % load(sprintf('%s/%s/masks.mat', prefix_proc, par.lat_interp)); % load land and ocean masks

    if strcmp(par.lat_interp, 'std')
        lat = par.lat_std;
    else
        lat = grid.dim3.lat;
    end

    % interpolate ta to standard lat grid
    masi_orig = permute(masi_orig, [2 1 3 4]);
    masi_orig = interp1(grid.dim3.lat, masi_orig, lat);
    masi_orig = permute(masi_orig, [2 1 3 4]);

    pasi_orig = permute(pasi_orig, [2 1 3 4]);
    pasi_orig = interp1(grid.dim3.lat, pasi_orig, lat);
    pasi_orig = permute(pasi_orig, [2 1 3 4]);

    masi_sm.lo = masi_orig; % surface is already masked in standard sigma coordinates
    pasi_sm.lo = pasi_orig; % surface is already masked in spndard sigma coordinates

    masi_sm.lo = permute(masi_sm.lo, [1 2 4 3]); % bring plev to last dimension

    pasi_sm.lo = permute(pasi_sm.lo, [1 2 4 3]); % bring plev to last dimension

    % mask_t.land = nanmean(mask.land, 3);
    % mask_t.ocean = nanmean(mask.ocean, 3);

    for l = {'lo'}; land = l{1}; % over land, over ocean, or both
        masi.(land)= squeeze(nanmean(masi_sm.(land), 1)); % zonal average
        pasi.(land)= squeeze(nanmean(pasi_sm.(land), 1)); % zonal average
    end

    for l = {'lo'}; land = l{1}; % over land, over ocean, or both
        % take time averages
        for t = {'ann', 'djf', 'jja', 'mam', 'son'}; time = t{1};
            if strcmp(time, 'ann')
                masi_t.(land).(time) = squeeze(nanmean(masi_sm.(land), 3));
                pasi_t.(land).(time) = squeeze(nanmean(pasi_sm.(land), 3));
            elseif strcmp(time, 'djf')
                masi_shift.(land) = circshift(masi_sm.(land), 1, 3);
                pasi_shift.(land) = circshift(pasi_sm.(land), 1, 3);
                masi_t.(land).(time) = squeeze(nanmean(masi_shift.(land)(:,:,1:3,:), 3));
                pasi_t.(land).(time) = squeeze(nanmean(pasi_shift.(land)(:,:,1:3,:), 3));
            elseif strcmp(time, 'jja')
                masi_t.(land).(time) = squeeze(nanmean(masi_sm.(land)(:,:,6:8,:), 3));
                pasi_t.(land).(time) = squeeze(nanmean(pasi_sm.(land)(:,:,6:8,:), 3));
            elseif strcmp(time, 'mam')
                masi_t.(land).(time) = squeeze(nanmean(masi_sm.(land)(:,:,3:5,:), 3));
                pasi_t.(land).(time) = squeeze(nanmean(pasi_sm.(land)(:,:,3:5,:), 3));
            elseif strcmp(time, 'son')
                masi_t.(land).(time) = squeeze(nanmean(masi_sm.(land)(:,:,9:11,:), 3));
                pasi_t.(land).(time) = squeeze(nanmean(pasi_sm.(land)(:,:,9:11,:), 3));
            end
        end
    end

    % save filtered data

    if par.do_surf; save(printname, 'masi', 'lat');
    else save(printname, 'masi', 'pasi', 'lat', '-v7.3'); end

    if par.do_surf; save(printname2, 'masi_t', 'lat');
    else save(printname2, 'masi_t', 'pasi_t', 'lat', '-v7.3'); end
end % compute mon x lat temperature field
