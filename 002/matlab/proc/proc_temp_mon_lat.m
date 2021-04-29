function proc_temp_mon_lat(type, par)

    prefix = make_prefix(type, par);
    prefix_proc = make_prefix_proc(type, par);
    foldername = make_savedir_proc(type, par);

    load(sprintf('%s/grid.mat', prefix)); % read grid data
    if par.levtype == 'ml'
        load(sprintf('%s/taml_si.mat', prefix)); tasi_orig = ta_si.spl; clear ta_si; % read temp in si coordinates
    elseif par.levtype == 'pl'
        load(sprintf('%s/ta_si.mat', prefix)); tasi_orig = ta_si.spl; clear ta_si; % read temp in si coordinates
    end
    
    %load(sprintf('%s/pa_si.mat', prefix)); pasi_orig = pa_si; clear pa_si; % read temp in si coordinates
    load(sprintf('%s/srfc.mat', prefix)); % load surface data
    % load(sprintf('%s/masks.mat', prefix_proc)); % load land and ocean masks

    if strcmp(par.lat_interp, 'std')
        lat = par.lat_std;
    else
        lat = grid.dim3.lat;
    end

    % interpolate ta to standard lat grid
    tasi_orig = permute(tasi_orig, [2 1 3 4]);
    tasi_orig = interp1(grid.dim3.lat, tasi_orig, lat);
    tasi_orig = permute(tasi_orig, [2 1 3 4]);

    %pasi_orig = permute(pasi_orig, [2 1 3 4]);
    %pasi_orig = interp1(grid.dim3.lat, pasi_orig, lat);
    %pasi_orig = permute(pasi_orig, [2 1 3 4]);

    tasi_sm.lo = tasi_orig; % surface is already masked in standard sigma coordinates
    %pasi_sm.lo = pasi_orig; % surface is already masked in spndard sigma coordinates

    tasi_sm.lo = permute(tasi_sm.lo, [1 2 4 3]); % bring plev to last dimension

    %pasi_sm.lo = permute(pasi_sm.lo, [1 2 4 3]); % bring plev to last dimension

    % mask_vert.land = repmat(mask.land, [1 1 1 size(tasi_sm.lo, 4)]);
    % mask_vert.ocean = repmat(mask.ocean, [1 1 1 size(tasi_sm.lo, 4)]);

    % tasi_sm.l = tasi_sm.lo .* mask_vert.ocean;
    % tasi_sm.o = tasi_sm.lo .* mask_vert.land;

    for l = par.land_list; land = l{1}; % over land, over ocean, or both
    % for l = {'lo'}; land = l{1}; % over land, over ocean, or both
        tasi.(land)= squeeze(nanmean(tasi_sm.(land), 1)); % zonal average
        %pasi.(land)= squeeze(nanmean(pasi_sm.(land), 1)); % zonal average
    end

    for l = par.land_list; land = l{1}; % over land, over ocean, or both
    % for l = {'lo'}; land = l{1}; % over land, over ocean, or both
        % take time averages
        for t = {'ann', 'djf', 'jja', 'mam', 'son'}; time = t{1};
            if strcmp(time, 'ann')
                tasi_t.(land).(time) = squeeze(nanmean(tasi_sm.(land), 3));
                %pasi_t.(land).(time) = squeeze(nanmean(pasi_sm.(land), 3));
            elseif strcmp(time, 'djf')
                tasi_shift.(land) = circshift(tasi_sm.(land), 1, 3);
                %pasi_shift.(land) = circshift(pasi_sm.(land), 1, 3);
                tasi_t.(land).(time) = squeeze(nanmean(tasi_shift.(land)(:,:,1:3,:), 3));
                %pasi_t.(land).(time) = squeeze(nanmean(pasi_shift.(land)(:,:,1:3,:), 3));
            elseif strcmp(time, 'jja')
                tasi_t.(land).(time) = squeeze(nanmean(tasi_sm.(land)(:,:,6:8,:), 3));
                %pasi_t.(land).(time) = squeeze(nanmean(pasi_sm.(land)(:,:,6:8,:), 3));
            elseif strcmp(time, 'mam')
                tasi_t.(land).(time) = squeeze(nanmean(tasi_sm.(land)(:,:,3:5,:), 3));
                %pasi_t.(land).(time) = squeeze(nanmean(pasi_sm.(land)(:,:,3:5,:), 3));
            elseif strcmp(time, 'son')
                tasi_t.(land).(time) = squeeze(nanmean(tasi_sm.(land)(:,:,9:11,:), 3));
                %pasi_t.(land).(time) = squeeze(nanmean(pasi_sm.(land)(:,:,9:11,:), 3));
            end
        end
    end

    % save filtered data
    if par.levtype == 'ml'
        printname = [foldername 'taml_mon_lat'];
    elseif par.levtype == 'pl'
        printname = [foldername 'ta_mon_lat'];
    end
    save(printname, 'tasi', 'lat');
    %if par.do_surf; save(printname, 'tasi', 'lat');
    %else save(printname, 'tasi', 'pasi', 'lat', '-v7.3'); end

    if par.levtype == 'ml'
        printname = [foldername 'taml_lon_lat'];
    elseif par.levtype == 'pl'
        printname = [foldername 'ta_lon_lat'];
    end
    save(printname, 'tasi_t', 'lat');
    %else save(printname, 'tasi_t', 'pasi_t', 'lat', '-v7.3'); end
end % compute mon x lat temperature field
