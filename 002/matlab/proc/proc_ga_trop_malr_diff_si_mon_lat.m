function proc_ga_trop_malr_diff_si_mon_lat(type, par)
    if strcmp(type, 'era5') | strcmp(type, 'erai') | strcmp(type, 'era5c')
        foldername = sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s', type, par.(type).yr_span, par.lat_interp);
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s', type, par.(type).yr_span);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s', type, par.(type).yr_span);
        file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/%s/temp/%s_temp_%s.ymonmean.nc', type, type, par.(type).yr_span));
        fullpath=sprintf('%s/%s', file.folder, file.name);
        ta_orig = double(ncread(fullpath, 't'));
    elseif strcmp(type, 'gcm')
        foldername = sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s/%s/', type, par.model, par.gcm.clim, par.lat_interp);
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s/%s', type, par.model, par.gcm.clim);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s', type, par.model, par.gcm.clim);
        var = 'ta';
        file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/gcm/%s/%s_Amon_%s_%s_r1i1p1_*.ymonmean.nc', par.model, var, par.model, par.gcm.clim));
        fullpath=sprintf('%s/%s', file.folder, file.name);
        ta_orig = double(ncread(fullpath, var));
    elseif strcmp(type, 'echam')
        foldername = sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s/', type, par.echam.clim, par.lat_interp);
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s', type, par.echam.clim);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s', type, par.echam.clim);
        var = 't';
        if contains(par.echam.clim, 'rp000')
            file=dir(sprintf('/project2/tas1/ockham/data11/tas/echam-aiv_rcc_6.1.00p1/%s/ATM_%s_0020_39.nc', par.echam.clim, par.echam.clim));
        else
            file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/echam/ATM*_%s_*.ymonmean.nc', par.echam.clim));
        end
        fullpath=sprintf('%s/%s', file.folder, file.name);
        ta_orig = double(ncread(fullpath, var));
    end
    load(sprintf('%s/grid.mat', prefix)); % read grid data
    load(sprintf('%s/dtdzzsi.mat', prefix)); % read temp in si coordinates
    load(sprintf('%s/malrzsi.mat', prefix)); % read temp in si coordinates
    % load(sprintf('%s/%s/masks.mat', prefix_proc, par.lat_interp)); % load land and ocean masks

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
