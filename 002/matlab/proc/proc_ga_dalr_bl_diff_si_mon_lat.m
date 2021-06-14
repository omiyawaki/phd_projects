function proc_ga_dalr_bl_diff_si_mon_lat(type, par)

    prefix = make_prefix(type, par);
    prefix_proc = make_prefix_proc(type, par);
    foldername = make_savedir_si_bl(type, par);
    ta_orig = load_temp(type, par);
    
    load(sprintf('%s/grid.mat', prefix)); % read grid data
    load(sprintf('%s/dtdzsi.mat', prefix)); dtdzzsi = dtdzsi; clear dtdzsi; % read temp in si coordinates
    % load(sprintf('%s/dtdzzsi.mat', prefix)); % read temp in si coordinates
    % load(sprintf('%s/masks.mat', prefix_proc)); % load land and ocean masks

    if contains(par.model, 'GISS-E2')
        dtdzzsi = permute(dtdzzsi, [2 1 3 4]);
        dtdzzsi = interp1(grid.dim3.lat_zg, dtdzzsi, grid.dim3.lat);
        dtdzzsi = permute(dtdzzsi, [2 1 3 4]);
    end

    dtmdzzsi = 1e3*par.g/par.cpd*ones(size(dtdzzsi));

    if strcmp(par.lat_interp, 'std')
        lat = par.lat_std;
    else
        lat = grid.dim3.lat;
    end

    ga_dalr_bl_diff_orig = (dtmdzzsi - dtdzzsi)./dtmdzzsi * 1e2; % percentage difference of dry adiabatic vs GCM lapse rates
    clear dtmdz dtdz;
    ga_dalr_bl_diff_orig = permute(ga_dalr_bl_diff_orig,[2 1 3 4]); % bring lat to 1st
    ga_dalr_bl_diff_orig = interp1(grid.dim3.lat, ga_dalr_bl_diff_orig, lat); % interpolate to standard lat
    ga_dalr_bl_diff_orig = permute(ga_dalr_bl_diff_orig, [3 2 1 4]); % bring height front
    ga_dalr_bl_diff_orig = interp1(grid.dim3.si, ga_dalr_bl_diff_orig, linspace(1,par.si_bl,100)); % prepare to average between 1000 hPa to top of boundary layer as specified by par.si_bl
    ga_dalr_bl_diff_orig = squeeze(nanmean(ga_dalr_bl_diff_orig,1)); % take vertical average
    
    ga_dalr_bl_diff0.lo = ga_dalr_bl_diff_orig;
    % ga_dalr_bl_diff0.l = ga_dalr_bl_diff0.lo.*mask.ocean; % filter ga_dalr_bl_diff0 with surface mask
    % ga_dalr_bl_diff0.o = ga_dalr_bl_diff0.lo.*mask.land; % filter ga_dalr_bl_diff0 with surface mask

    for l = par.land_list; land = l{1}; % over land, over ocean, or both
    % for l = {'lo'}; land = l{1}; % over land, over ocean, or both
        ga_dalr_bl_diff.(land)= squeeze(nanmean(ga_dalr_bl_diff0.(land), 1)); % zonal average

        % take time averages
        for t = {'ann', 'djf', 'jja', 'mam', 'son'}; time = t{1};
            if strcmp(time, 'ann')
                ga_dalr_bl_diff_t.(land).(time) = squeeze(nanmean(ga_dalr_bl_diff0.(land), 3));
            elseif strcmp(time, 'djf')
                ga_dalr_bl_diff_shift.(land) = circshift(ga_dalr_bl_diff0.(land), 1, 3);
                ga_dalr_bl_diff_t.(land).(time) = squeeze(nanmean(ga_dalr_bl_diff_shift.(land)(:,:,1:3), 3));
            elseif strcmp(time, 'jja')
                ga_dalr_bl_diff_t.(land).(time) = squeeze(nanmean(ga_dalr_bl_diff0.(land)(:,:,6:8), 3));
            elseif strcmp(time, 'mam')
                ga_dalr_bl_diff_t.(land).(time) = squeeze(nanmean(ga_dalr_bl_diff0.(land)(:,:,3:5), 3));
            elseif strcmp(time, 'son')
                ga_dalr_bl_diff_t.(land).(time) = squeeze(nanmean(ga_dalr_bl_diff0.(land)(:,:,9:11), 3));
            end
            ga_dalr_bl_diff_zt.(land).(time)= squeeze(nanmean(ga_dalr_bl_diff_t.(land).(time), 1)); % zonal average
        end
    end

    % save filtered data
    printname = [foldername 'ga_dalr_bl_diff_si_mon_lat'];
    if ~exist(foldername, 'dir')
        mkdir(foldername)
    end
    save(printname, 'ga_dalr_bl_diff', 'lat');

    printname = [foldername 'ga_dalr_bl_diff_si_lon_lat'];
    save(printname, 'ga_dalr_bl_diff_t', 'lat');

    printname = [foldername 'ga_dalr_bl_diff_si_lat'];
    save(printname, 'ga_dalr_bl_diff_zt', 'lat');

end % compute mon x lat gamma percentage difference field with land/ocean masking
