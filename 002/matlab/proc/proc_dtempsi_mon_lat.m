function proc_dtempsi_mon_lat(type, par)

    prefix = make_prefix(type, par);
    prefix_proc = make_prefix_proc(type, par);
    foldername = make_savedir_proc(type, par);

    load(sprintf('%s/grid.mat', prefix)); % read grid data
    load(sprintf('%s/dta_si.mat', prefix));

    dta_si = dta_si.spl;

    if strcmp(type, 'gcm')
        if contains(par.model, 'GISS-E2')
            dta_si = permute(dta_si, [2 1 3 4]);
            dta_si = interp1(grid.dim3.lat, dta_si, grid.dim3.lat);
            dta_si = permute(dta_si, [2 1 3 4]);
        end
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
    dta_si = permute(dta_si, [2 1 3 4]);
    dta_si = interp1(grid.dim3.lat, dta_si, lat);
    dta_si = permute(dta_si, [2 1 3 4]);

    %pasi_orig = permute(pasi_orig, [2 1 3 4]);
    %pasi_orig = interp1(grid.dim3.lat, pasi_orig, lat);
    %pasi_orig = permute(pasi_orig, [2 1 3 4]);

    dtempsi_sm.lo = dta_si; % surface is already masked in standard sigma coordinates
    %pasi_sm.lo = pasi_orig; % surface is already masked in spndard sigma coordinates

    dtempsi_sm.lo = permute(dtempsi_sm.lo, [1 2 4 3]); % bring plev to last dimension

    %pasi_sm.lo = permute(pasi_sm.lo, [1 2 4 3]); % bring plev to last dimension

    % mask_vert.land = repmat(mask.land, [1 1 1 size(dtempsi_sm.lo, 4)]);
    % mask_vert.ocean = repmat(mask.ocean, [1 1 1 size(dtempsi_sm.lo, 4)]);

    % dtempsi_sm.l = dtempsi_sm.lo .* mask_vert.ocean;
    % dtempsi_sm.o = dtempsi_sm.lo .* mask_vert.land;

    for l = par.land_list; land = l{1}; % over land, over ocean, or both
    % for l = {'lo'}; land = l{1}; % over land, over ocean, or both
        dtempsi.(land)= squeeze(nanmean(dtempsi_sm.(land), 1)); % zonal average
        %pasi.(land)= squeeze(nanmean(pasi_sm.(land), 1)); % zonal average
    end

    for l = par.land_list; land = l{1}; % over land, over ocean, or both
    % for l = {'lo'}; land = l{1}; % over land, over ocean, or both
        % take time averages
        for t = {'ann', 'djf', 'jja', 'mam', 'son'}; time = t{1};
            if strcmp(time, 'ann')
                dtempsi_t.(land).(time) = squeeze(nanmean(dtempsi_sm.(land), 3));
                %pasi_t.(land).(time) = squeeze(nanmean(pasi_sm.(land), 3));
            elseif strcmp(time, 'djf')
                dtempsi_shift.(land) = circshift(dtempsi_sm.(land), 1, 3);
                %pasi_shift.(land) = circshift(pasi_sm.(land), 1, 3);
                dtempsi_t.(land).(time) = squeeze(nanmean(dtempsi_shift.(land)(:,:,1:3,:), 3));
                %pasi_t.(land).(time) = squeeze(nanmean(pasi_shift.(land)(:,:,1:3,:), 3));
            elseif strcmp(time, 'jja')
                dtempsi_t.(land).(time) = squeeze(nanmean(dtempsi_sm.(land)(:,:,6:8,:), 3));
                %pasi_t.(land).(time) = squeeze(nanmean(pasi_sm.(land)(:,:,6:8,:), 3));
            elseif strcmp(time, 'mam')
                dtempsi_t.(land).(time) = squeeze(nanmean(dtempsi_sm.(land)(:,:,3:5,:), 3));
                %pasi_t.(land).(time) = squeeze(nanmean(pasi_sm.(land)(:,:,3:5,:), 3));
            elseif strcmp(time, 'son')
                dtempsi_t.(land).(time) = squeeze(nanmean(dtempsi_sm.(land)(:,:,9:11,:), 3));
                %pasi_t.(land).(time) = squeeze(nanmean(pasi_sm.(land)(:,:,9:11,:), 3));
            end
        end
    end

    % save filtered data
    printname = [foldername 'dtempsi_mon_lat'];
    save(printname, 'dtempsi', 'lat', '-v7.3');
    %if par.do_surf; save(printname, 'dtempsi', 'lat');
    %else save(printname, 'dtempsi', 'pasi', 'lat', '-v7.3'); end

    printname = [foldername 'dtempsi_lon_lat'];
    save(printname, 'dtempsi_t', 'lat', '-v7.3');
    %else save(printname, 'dtempsi_t', 'pasi_t', 'lat', '-v7.3'); end
end % compute mon x lat temperature field
