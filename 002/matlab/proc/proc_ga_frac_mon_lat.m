function proc_ga_frac_mon_lat(type, par)

    prefix = make_prefix(type, par);
    prefix_proc = make_prefix_proc(type, par);
    foldername = make_savedir_proc(type, par);

    load(sprintf('%s/grid.mat', prefix)); % read grid data
    load(sprintf('%s/dtdzsi.mat', prefix));
    load(sprintf('%s/malrsi.mat', prefix));
    ga_frac_orig = (dtmdzsi-dtdzsi)./dtmdzsi; % moist adiabatic lapse rate minus actual lapse rate
    
    %load(sprintf('%s/pa_si.mat', prefix)); pasi_orig = pa_si; clear pa_si; % read temp in si coordinates
    load(sprintf('%s/srfc.mat', prefix)); % load surface data
    load(sprintf('%s/masks.mat', prefix_proc)); % load land and ocean masks

    if strcmp(par.lat_interp, 'std')
        lat = par.lat_std;
    else
        lat = grid.dim3.lat;
    end

    % interpolate ta to standard lat grid
    ga_frac_orig = permute(ga_frac_orig, [2 1 3 4]);
    ga_frac_orig = interp1(grid.dim3.lat, ga_frac_orig, lat);
    ga_frac_orig = permute(ga_frac_orig, [2 1 3 4]);

    %pasi_orig = permute(pasi_orig, [2 1 3 4]);
    %pasi_orig = interp1(grid.dim3.lat, pasi_orig, lat);
    %pasi_orig = permute(pasi_orig, [2 1 3 4]);

    ga_frac_sm.lo = ga_frac_orig; % surface is already masked in standard sigma coordinates
    %pasi_sm.lo = pasi_orig; % surface is already masked in spndard sigma coordinates

    ga_frac_sm.lo = permute(ga_frac_sm.lo, [1 2 4 3]); % bring plev to last dimension

    %pasi_sm.lo = permute(pasi_sm.lo, [1 2 4 3]); % bring plev to last dimension

    mask_vert.land = repmat(mask.land, [1 1 1 size(ga_frac_sm.lo, 4)]);
    mask_vert.ocean = repmat(mask.ocean, [1 1 1 size(ga_frac_sm.lo, 4)]);

    ga_frac_sm.l = ga_frac_sm.lo .* mask_vert.ocean;
    ga_frac_sm.o = ga_frac_sm.lo .* mask_vert.land;

    for l = {'lo', 'l', 'o'}; land = l{1}; % over land, over ocean, or both
    % for l = {'lo'}; land = l{1}; % over land, over ocean, or both
        ga_frac.(land)= squeeze(nanmean(ga_frac_sm.(land), 1)); % zonal average
        %pasi.(land)= squeeze(nanmean(pasi_sm.(land), 1)); % zonal average
    end

    for l = {'lo', 'l', 'o'}; land = l{1}; % over land, over ocean, or both
    % for l = {'lo'}; land = l{1}; % over land, over ocean, or both
        % take time averages
        for t = {'ann', 'djf', 'jja', 'mam', 'son'}; time = t{1};
            if strcmp(time, 'ann')
                ga_frac_t.(land).(time) = squeeze(nanmean(ga_frac_sm.(land), 3));
                %pasi_t.(land).(time) = squeeze(nanmean(pasi_sm.(land), 3));
            elseif strcmp(time, 'djf')
                ga_frac_shift.(land) = circshift(ga_frac_sm.(land), 1, 3);
                %pasi_shift.(land) = circshift(pasi_sm.(land), 1, 3);
                ga_frac_t.(land).(time) = squeeze(nanmean(ga_frac_shift.(land)(:,:,1:3,:), 3));
                %pasi_t.(land).(time) = squeeze(nanmean(pasi_shift.(land)(:,:,1:3,:), 3));
            elseif strcmp(time, 'jja')
                ga_frac_t.(land).(time) = squeeze(nanmean(ga_frac_sm.(land)(:,:,6:8,:), 3));
                %pasi_t.(land).(time) = squeeze(nanmean(pasi_sm.(land)(:,:,6:8,:), 3));
            elseif strcmp(time, 'mam')
                ga_frac_t.(land).(time) = squeeze(nanmean(ga_frac_sm.(land)(:,:,3:5,:), 3));
                %pasi_t.(land).(time) = squeeze(nanmean(pasi_sm.(land)(:,:,3:5,:), 3));
            elseif strcmp(time, 'son')
                ga_frac_t.(land).(time) = squeeze(nanmean(ga_frac_sm.(land)(:,:,9:11,:), 3));
                %pasi_t.(land).(time) = squeeze(nanmean(pasi_sm.(land)(:,:,9:11,:), 3));
            end
        end
    end

    % save filtered data
    printname = [foldername 'ga_frac_mon_lat'];
    save(printname, 'ga_frac', 'lat', '-v7.3');
    %if par.do_surf; save(printname, 'ga_frac', 'lat');
    %else save(printname, 'ga_frac', 'pasi', 'lat', '-v7.3'); end

    printname = [foldername 'ga_frac_lon_lat'];
    save(printname, 'ga_frac_t', 'lat', '-v7.3');
    %else save(printname, 'ga_frac_t', 'pasi_t', 'lat', '-v7.3'); end
end % compute mon x lat temperature field
