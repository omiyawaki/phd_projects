function proc_temp_pl_mon_lat(type, par)

    prefix = make_prefix(type, par);
    foldername = make_savedir_proc(type, par);

    load(sprintf('%s/grid.mat', prefix)); % read grid data
    load(sprintf('%s/srfc.mat', prefix)); % read grid data
    ta_orig = load_temp(type, par);
    
    % load(sprintf('%s/%s/masks.mat', prefix_proc, par.lat_interp)); % load land and ocean masks

    if strcmp(par.lat_interp, 'std')
        lat = par.lat_std;
    else
        lat = grid.dim3.lat;
    end

    % interpolate ta to standard lat grid
    ta_orig = permute(ta_orig, [2 1 3 4]);
    ta_orig = interp1(grid.dim3.lat, ta_orig, lat);
    ta_orig = permute(ta_orig, [2 1 3 4]);

    [ps_vert, pa] = make_surfmask_vars(grid, type, srfc);
    sm = nan(size(ta_orig));
    sm(pa < ps_vert) = 1;

    ta_sm.lo = ta_orig.*sm; % mask out data below surface pressure

    ta_sm.lo = permute(ta_sm.lo, [1 2 4 3]); % bring plev to last dimension


    % mask_t.land = nanmean(mask.land, 3);
    % mask_t.ocean = nanmean(mask.ocean, 3);

    for l = {'lo'}; land = l{1}; % over land, over ocean, or both
        ta.(land)= squeeze(nanmean(ta_sm.(land), 1)); % zonal average
    end

    for l = {'lo'}; land = l{1}; % over land, over ocean, or both
        % take time averages
        for t = {'ann', 'djf', 'jja', 'mam', 'son'}; time = t{1};
            if strcmp(time, 'ann')
                ta_t.(land).(time) = squeeze(nanmean(ta_sm.(land), 3));
            elseif strcmp(time, 'djf')
                ta_shift.(land) = circshift(ta_sm.(land), 1, 3);
                ta_t.(land).(time) = squeeze(nanmean(ta_shift.(land)(:,:,1:3,:), 3));
            elseif strcmp(time, 'jja')
                ta_t.(land).(time) = squeeze(nanmean(ta_sm.(land)(:,:,6:8,:), 3));
            elseif strcmp(time, 'mam')
                ta_t.(land).(time) = squeeze(nanmean(ta_sm.(land)(:,:,3:5,:), 3));
            elseif strcmp(time, 'son')
                ta_t.(land).(time) = squeeze(nanmean(ta_sm.(land)(:,:,9:11,:), 3));
            end
        end
    end

    % save filtered data
    printname = [foldername 'ta_pl_mon_lat'];
    save(printname, 'ta', 'lat');

    printname = [foldername 'ta_pl_lon_lat'];
    save(printname, 'ta_t', 'lat');
end % compute mon x lat temperature field
