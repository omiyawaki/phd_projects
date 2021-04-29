function proc_mse_mon_lat(type, par)
    
    prefix = make_prefix(type, par);
    foldername = make_savedir_proc(type, par);

    load(sprintf('%s/grid.mat', prefix)); % read grid data
    load(sprintf('%s/mse_si.mat', prefix)); msesi_orig = mse_si.spl; clear mse_si; % read temp in si coordinates
    load(sprintf('%s/pa_si.mat', prefix)); pasi_orig = pa_si; clear pa_si; % read temp in si coordinates
    load(sprintf('%s/srfc.mat', prefix)); % load surface data
    % load(sprintf('%s/%s/masks.mat', prefix_proc, par.lat_interp)); % load land and ocean masks

    if strcmp(par.lat_interp, 'std')
        lat = par.lat_std;
    else
        lat = grid.dim3.lat;
    end

    % interpolate mse to smsendard lat grid
    msesi_orig = permute(msesi_orig, [2 1 3 4]);
    msesi_orig = interp1(grid.dim3.lat, msesi_orig, lat);
    msesi_orig = permute(msesi_orig, [2 1 3 4]);

    pasi_orig = permute(pasi_orig, [2 1 3 4]);
    pasi_orig = interp1(grid.dim3.lat, pasi_orig, lat);
    pasi_orig = permute(pasi_orig, [2 1 3 4]);

    msesi_sm.lo = msesi_orig; % surface is already masked in smsendard sigma coordinates
    pasi_sm.lo = pasi_orig; % surface is already masked in spndard sigma coordinates

    msesi_sm.lo = permute(msesi_sm.lo, [1 2 4 3]); % bring plev to last dimension

    pasi_sm.lo = permute(pasi_sm.lo, [1 2 4 3]); % bring plev to last dimension

    % mask_t.land = nanmean(mask.land, 3);
    % mask_t.ocean = nanmean(mask.ocean, 3);

    for l = {'lo'}; land = l{1}; % over land, over ocean, or both
        msesi.(land)= squeeze(mean(msesi_sm.(land), 1)); % zonal average
        pasi.(land)= squeeze(mean(pasi_sm.(land), 1)); % zonal average
    end

    for l = {'lo'}; land = l{1}; % over land, over ocean, or both
        % mseke time averages
        for t = {'ann', 'djf', 'jja', 'mam', 'son'}; time = t{1};
            if strcmp(time, 'ann')
                msesi_t.(land).(time) = squeeze(mean(msesi_sm.(land), 3));
                pasi_t.(land).(time) = squeeze(mean(pasi_sm.(land), 3));
            elseif strcmp(time, 'djf')
                msesi_shift.(land) = circshift(msesi_sm.(land), 1, 3);
                pasi_shift.(land) = circshift(pasi_sm.(land), 1, 3);
                msesi_t.(land).(time) = squeeze(mean(msesi_shift.(land)(:,:,1:3,:), 3));
                pasi_t.(land).(time) = squeeze(mean(pasi_shift.(land)(:,:,1:3,:), 3));
            elseif strcmp(time, 'jja')
                msesi_t.(land).(time) = squeeze(mean(msesi_sm.(land)(:,:,6:8,:), 3));
                pasi_t.(land).(time) = squeeze(mean(pasi_sm.(land)(:,:,6:8,:), 3));
            elseif strcmp(time, 'mam')
                msesi_t.(land).(time) = squeeze(mean(msesi_sm.(land)(:,:,3:5,:), 3));
                pasi_t.(land).(time) = squeeze(mean(pasi_sm.(land)(:,:,3:5,:), 3));
            elseif strcmp(time, 'son')
                msesi_t.(land).(time) = squeeze(mean(msesi_sm.(land)(:,:,9:11,:), 3));
                pasi_t.(land).(time) = squeeze(mean(pasi_sm.(land)(:,:,9:11,:), 3));
            end
        end
    end

    % save filtered damse
    printname = [foldername 'mse_mon_lat'];
    if ~exist(foldername, 'dir')
        mkdir(foldername)
    end
    if par.do_surf; save(printname, 'msesi', 'lat');
    else save(printname, 'msesi', 'pasi', 'lat', '-v7.3'); end

    printname = [foldername 'mse_lon_lat'];
    if ~exist(foldername, 'dir')
        mkdir(foldername)
    end
    if par.do_surf; save(printname, 'msesi_t', 'lat');
    else save(printname, 'msesi_t', 'pasi_t', 'lat', '-v7.3'); end
end % compute mon x lat temperature field
