function proc_ga_malr_diff_si_mon_lat(type, par)

    prefix = make_prefix(type, par);
    prefix_proc = make_prefix_proc(type, par);
    foldername = make_savedir_si_bl(type, par);
    ta_orig = load_temp(type, par);

    load(sprintf('%s/grid.mat', prefix)); % read grid data
    % load(sprintf('%s/dtdzzsi.mat', prefix)); % read temp in si coordinates
    load(sprintf('%s/dtdzsi.mat', prefix)); dtdzzsi = dtdzsi; clear dtdzsi; % read temp in si coordinates
    % load(sprintf('%s/malrzsi.mat', prefix)); % read temp in si coordinates
    load(sprintf('%s/malrsi.mat', prefix)); dtmdzzsi = dtmdzsi; clear dtmdzsi; % read temp in si coordinates
    % load(sprintf('%s/%s/masks.mat', prefix_proc, par.lat_interp)); % load land and ocean masks
    
    dtmdzzsi(1,72,:,3)
    return

    if strcmp(par.lat_interp, 'std')
        lat = par.lat_std;
    else
        lat = grid.dim3.lat;
    end

    ga_malr_diff_orig = (dtmdzzsi - dtdzzsi)./dtmdzzsi * 1e2; % percentage difference of moist adiabatic vs GCM lapse rates
    clear dtmdz dtdz;
    ga_malr_diff_orig = permute(ga_malr_diff_orig,[2 1 3 4]); % bring lat to 1st
    ga_malr_diff_orig = interp1(grid.dim3.lat, ga_malr_diff_orig, lat); % interpolate to standard lat
    ga_malr_diff_orig = permute(ga_malr_diff_orig, [3 2 1 4]); % bring height front
    ga_malr_diff_orig = interp1(grid.dim3.si, ga_malr_diff_orig, linspace(par.si_bl,par.si_up,100)); % prepare to average between 1000-200 hPa
    ga_malr_diff_orig = squeeze(nanmean(ga_malr_diff_orig,1)); % take vertical average

    ga_malr_diff0.lo = ga_malr_diff_orig;
    % ga_malr_diff0.l = ga_malr_diff0.lo.*mask.ocean; % filter ga_malr_diff0 with surface mask
    % ga_malr_diff0.o = ga_malr_diff0.lo.*mask.land; % filter ga_malr_diff0 with surface mask

    % for l = {'lo', 'l', 'o'}; land = l{1}; % over land, over ocean, or both
    for l = {'lo'}; land = l{1}; % over land, over ocean, or both
        ga_malr_diff.(land)= squeeze(nanmean(ga_malr_diff0.(land), 1)); % zonal average
        % take time averages
        for t = {'ann', 'djf', 'jja', 'mam', 'son'}; time = t{1};
            if strcmp(time, 'ann')
                ga_malr_diff_t.(land).(time) = squeeze(nanmean(ga_malr_diff0.(land), 3));
            elseif strcmp(time, 'djf')
                ga_malr_diff_shift.(land) = circshift(ga_malr_diff0.(land), 1, 3);
                ga_malr_diff_t.(land).(time) = squeeze(nanmean(ga_malr_diff_shift.(land)(:,:,1:3,:), 3));
            elseif strcmp(time, 'jja')
                ga_malr_diff_t.(land).(time) = squeeze(nanmean(ga_malr_diff0.(land)(:,:,6:8,:), 3));
            elseif strcmp(time, 'mam')
                ga_malr_diff_t.(land).(time) = squeeze(nanmean(ga_malr_diff0.(land)(:,:,3:5,:), 3));
            elseif strcmp(time, 'son')
                ga_malr_diff_t.(land).(time) = squeeze(nanmean(ga_malr_diff0.(land)(:,:,9:11,:), 3));
            end
        end
    end

    % save filtered data
    printname = [foldername 'ga_malr_diff_si_mon_lat'];
    if ~exist(foldername, 'dir')
        mkdir(foldername)
    end
    save(printname, 'ga_malr_diff', 'lat');

    printname = [foldername 'ga_malr_diff_si_lon_lat'];
    if ~exist(foldername, 'dir')
        mkdir(foldername)
    end
    save(printname, 'ga_malr_diff_t', 'lat');
end % compute mon x lat gamma percentage difference field with land/ocean masking
