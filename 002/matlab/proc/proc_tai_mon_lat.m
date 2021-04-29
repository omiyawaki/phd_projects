function proc_tai_mon_lat(type, par)

    prefix = make_prefix(type, par);
    foldername = make_savedir_proc(type, par);

    load(sprintf('%s/grid.mat', prefix)); % read grid data
    load(sprintf('%s/srfc.mat', prefix)); % read srfc data
    tmp = load(sprintf('%s/tai_simp.mat', prefix)); tai_orig = tmp.tai; clear tmp; % read tai data
    
    % load(sprintf('%s/%s/masks.mat', prefix_proc, par.lat_interp)); % load land and ocean masks

    if strcmp(par.lat_interp, 'std')
        lat = par.lat_std;
    else
        lat = grid.dim3.lat;
    end

    % interpolate tai to staindard lat grid
    tai_orig = permute(tai_orig, [2 1 3 4]);
    tai_orig = interp1(grid.dim3.lat, tai_orig, lat);
    tai_orig = permute(tai_orig, [2 1 3 4]);

    tai_sm.lo = tai_orig; % tai data is already surface masked

    tai_sm.lo = permute(tai_sm.lo, [1 2 4 3]); % bring plev to last dimension


    % mask_t.land = nanmean(mask.land, 3);
    % mask_t.ocean = nanmean(mask.ocean, 3);

    for l = par.land_list; land = l{1}; % over land, over ocean, or both
        tai.(land)= squeeze(nanmean(tai_sm.(land), 1)); % zonal average
    end

    for l = par.land_list; land = l{1}; % over land, over ocean, or both
        % taike time averages
        for t = {'ann', 'djf', 'jja', 'mam', 'son'}; time = t{1};
            if strcmp(time, 'ann')
                tai_t.(land).(time) = squeeze(nanmean(tai_sm.(land), 3));
            elseif strcmp(time, 'djf')
                tai_shift.(land) = circshift(tai_sm.(land), 1, 3);
                tai_t.(land).(time) = squeeze(nanmean(tai_shift.(land)(:,:,1:3,:), 3));
            elseif strcmp(time, 'jja')
                tai_t.(land).(time) = squeeze(nanmean(tai_sm.(land)(:,:,6:8,:), 3));
            elseif strcmp(time, 'mam')
                tai_t.(land).(time) = squeeze(nanmean(tai_sm.(land)(:,:,3:5,:), 3));
            elseif strcmp(time, 'son')
                tai_t.(land).(time) = squeeze(nanmean(tai_sm.(land)(:,:,9:11,:), 3));
            end
        end
    end

    % save filtered datai
    printname = [foldername 'tai_mon_lat'];
    save(printname, 'tai', 'lat');

    printname = [foldername 'tai_lon_lat'];
    save(printname, 'tai_t', 'lat');
end % compute mon x lat temperature field
