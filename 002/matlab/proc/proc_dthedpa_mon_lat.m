function proc_dthedpa_mon_lat(type, par)

    prefix = make_prefix(type, par);
    prefix_proc = make_prefix_proc(type, par);
    foldername = make_savedir_proc(type, par);

    load(sprintf('%s/grid.mat', prefix)); % read grid data
    tmp = load(sprintf('%s/dthedpasi.mat', prefix)); dthedpa_orig = tmp.dthedpasi; clear tmp;
    
    %load(sprintf('%s/pa_si.mat', prefix)); pasi_orig = pa_si; clear pa_si; % read temp in si coordinates
    % load(sprintf('%s/masks.mat', prefix_proc)); % load land and ocean masks

    if strcmp(par.lat_interp, 'std')
        lat = par.lat_std;
    else
        lat = grid.dim3.lat;
    end

    % interpolate ta to standard lat grid
    dthedpa_orig = permute(dthedpa_orig, [2 1 3 4]);
    dthedpa_orig = interp1(grid.dim3.lat, dthedpa_orig, lat);
    dthedpa_orig = permute(dthedpa_orig, [2 1 3 4]);

    %pasi_orig = permute(pasi_orig, [2 1 3 4]);
    %pasi_orig = interp1(grid.dim3.lat, pasi_orig, lat);
    %pasi_orig = permute(pasi_orig, [2 1 3 4]);

    dthedpa_sm.lo = dthedpa_orig; % surface is already masked in standard sigma coordinates
    %pasi_sm.lo = pasi_orig; % surface is already masked in spndard sigma coordinates

    dthedpa_sm.lo = permute(dthedpa_sm.lo, [1 2 4 3]); % bring plev to last dimension

    %pasi_sm.lo = permute(pasi_sm.lo, [1 2 4 3]); % bring plev to last dimension

    % mask_vert.land = repmat(mask.land, [1 1 1 size(dthedpa_sm.lo, 4)]);
    % mask_vert.ocean = repmat(mask.ocean, [1 1 1 size(dthedpa_sm.lo, 4)]);

    % dthedpa_sm.l = dthedpa_sm.lo .* mask_vert.ocean;
    % dthedpa_sm.o = dthedpa_sm.lo .* mask_vert.land;

    for l = par.land_list; land = l{1}; % over land, over ocean, or both
    % for l = {'lo'}; land = l{1}; % over land, over ocean, or both
        dthedpa.(land)= squeeze(nanmean(dthedpa_sm.(land), 1)); % zonal average
        %pasi.(land)= squeeze(nanmean(pasi_sm.(land), 1)); % zonal average
    end

    for l = par.land_list; land = l{1}; % over land, over ocean, or both
    % for l = {'lo'}; land = l{1}; % over land, over ocean, or both
        % take time averages
        for t = {'ann', 'djf', 'jja', 'mam', 'son'}; time = t{1};
            if strcmp(time, 'ann')
                dthedpa_t.(land).(time) = squeeze(nanmean(dthedpa_sm.(land), 3));
                %pasi_t.(land).(time) = squeeze(nanmean(pasi_sm.(land), 3));
            elseif strcmp(time, 'djf')
                dthedpa_shift.(land) = circshift(dthedpa_sm.(land), 1, 3);
                %pasi_shift.(land) = circshift(pasi_sm.(land), 1, 3);
                dthedpa_t.(land).(time) = squeeze(nanmean(dthedpa_shift.(land)(:,:,1:3,:), 3));
                %pasi_t.(land).(time) = squeeze(nanmean(pasi_shift.(land)(:,:,1:3,:), 3));
            elseif strcmp(time, 'jja')
                dthedpa_t.(land).(time) = squeeze(nanmean(dthedpa_sm.(land)(:,:,6:8,:), 3));
                %pasi_t.(land).(time) = squeeze(nanmean(pasi_sm.(land)(:,:,6:8,:), 3));
            elseif strcmp(time, 'mam')
                dthedpa_t.(land).(time) = squeeze(nanmean(dthedpa_sm.(land)(:,:,3:5,:), 3));
                %pasi_t.(land).(time) = squeeze(nanmean(pasi_sm.(land)(:,:,3:5,:), 3));
            elseif strcmp(time, 'son')
                dthedpa_t.(land).(time) = squeeze(nanmean(dthedpa_sm.(land)(:,:,9:11,:), 3));
                %pasi_t.(land).(time) = squeeze(nanmean(pasi_sm.(land)(:,:,9:11,:), 3));
            end
        end
    end

    % save filtered data
    printname = [foldername 'dthedpa_mon_lat'];
    save(printname, 'dthedpa', 'lat', '-v7.3');
    %if par.do_surf; save(printname, 'dthedpa', 'lat');
    %else save(printname, 'dthedpa', 'pasi', 'lat', '-v7.3'); end

    printname = [foldername 'dthedpa_lon_lat'];
    save(printname, 'dthedpa_t', 'lat', '-v7.3');
    %else save(printname, 'dthedpa_t', 'pasi_t', 'lat', '-v7.3'); end
end % compute mon x lat temperature field
