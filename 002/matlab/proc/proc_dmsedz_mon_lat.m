function proc_dmsedz_mon_lat(type, par)

    prefix = make_prefix(type, par);
    prefix_proc = make_prefix_proc(type, par);
    foldername = make_savedir_proc(type, par);

    load(sprintf('%s/grid.mat', prefix)); % read grid data
    tmp = load(sprintf('%s/dmsedzsi.mat', prefix)); dmsedz_orig = tmp.dmsedzsi; clear tmp;
    
    %load(sprintf('%s/pa_si.mat', prefix)); pasi_orig = pa_si; clear pa_si; % read temp in si coordinates
    % load(sprintf('%s/masks.mat', prefix_proc)); % load land and ocean masks

    if strcmp(par.lat_interp, 'std')
        lat = par.lat_std;
    else
        lat = grid.dim3.lat;
    end

    % interpolate ta to standard lat grid
    dmsedz_orig = permute(dmsedz_orig, [2 1 3 4]);
    dmsedz_orig = interp1(grid.dim3.lat, dmsedz_orig, lat);
    dmsedz_orig = permute(dmsedz_orig, [2 1 3 4]);

    %pasi_orig = permute(pasi_orig, [2 1 3 4]);
    %pasi_orig = interp1(grid.dim3.lat, pasi_orig, lat);
    %pasi_orig = permute(pasi_orig, [2 1 3 4]);

    dmsedz_sm.lo = dmsedz_orig; % surface is already masked in standard sigma coordinates
    %pasi_sm.lo = pasi_orig; % surface is already masked in spndard sigma coordinates

    dmsedz_sm.lo = permute(dmsedz_sm.lo, [1 2 4 3]); % bring plev to last dimension

    %pasi_sm.lo = permute(pasi_sm.lo, [1 2 4 3]); % bring plev to last dimension

    % mask_vert.land = repmat(mask.land, [1 1 1 size(dmsedz_sm.lo, 4)]);
    % mask_vert.ocean = repmat(mask.ocean, [1 1 1 size(dmsedz_sm.lo, 4)]);

    % dmsedz_sm.l = dmsedz_sm.lo .* mask_vert.ocean;
    % dmsedz_sm.o = dmsedz_sm.lo .* mask_vert.land;

    for l = par.land_list; land = l{1}; % over land, over ocean, or both
    % for l = {'lo'}; land = l{1}; % over land, over ocean, or both
        dmsedz.(land)= squeeze(nanmean(dmsedz_sm.(land), 1)); % zonal average
        %pasi.(land)= squeeze(nanmean(pasi_sm.(land), 1)); % zonal average
    end

    for l = par.land_list; land = l{1}; % over land, over ocean, or both
    % for l = {'lo'}; land = l{1}; % over land, over ocean, or both
        % take time averages
        for t = {'ann', 'djf', 'jja', 'mam', 'son'}; time = t{1};
            if strcmp(time, 'ann')
                dmsedz_t.(land).(time) = squeeze(nanmean(dmsedz_sm.(land), 3));
                %pasi_t.(land).(time) = squeeze(nanmean(pasi_sm.(land), 3));
            elseif strcmp(time, 'djf')
                dmsedz_shift.(land) = circshift(dmsedz_sm.(land), 1, 3);
                %pasi_shift.(land) = circshift(pasi_sm.(land), 1, 3);
                dmsedz_t.(land).(time) = squeeze(nanmean(dmsedz_shift.(land)(:,:,1:3,:), 3));
                %pasi_t.(land).(time) = squeeze(nanmean(pasi_shift.(land)(:,:,1:3,:), 3));
            elseif strcmp(time, 'jja')
                dmsedz_t.(land).(time) = squeeze(nanmean(dmsedz_sm.(land)(:,:,6:8,:), 3));
                %pasi_t.(land).(time) = squeeze(nanmean(pasi_sm.(land)(:,:,6:8,:), 3));
            elseif strcmp(time, 'mam')
                dmsedz_t.(land).(time) = squeeze(nanmean(dmsedz_sm.(land)(:,:,3:5,:), 3));
                %pasi_t.(land).(time) = squeeze(nanmean(pasi_sm.(land)(:,:,3:5,:), 3));
            elseif strcmp(time, 'son')
                dmsedz_t.(land).(time) = squeeze(nanmean(dmsedz_sm.(land)(:,:,9:11,:), 3));
                %pasi_t.(land).(time) = squeeze(nanmean(pasi_sm.(land)(:,:,9:11,:), 3));
            end
        end
    end

    % save filtered data
    printname = [foldername 'dmsedz_mon_lat'];
    save(printname, 'dmsedz', 'lat', '-v7.3');
    %if par.do_surf; save(printname, 'dmsedz', 'lat');
    %else save(printname, 'dmsedz', 'pasi', 'lat', '-v7.3'); end

    printname = [foldername 'dmsedz_lon_lat'];
    save(printname, 'dmsedz_t', 'lat', '-v7.3');
    %else save(printname, 'dmsedz_t', 'pasi_t', 'lat', '-v7.3'); end
end % compute mon x lat temperature field
