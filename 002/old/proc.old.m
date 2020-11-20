% old functions
function proc_ma_mon_lat_old(type, par)
% calculate moist adiabats
    if strcmp(type, 'era5') | strcmp(type, 'erai')
        foldername = sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/', type, par.lat_interp);
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s', type);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s', type);
    elseif strcmp(type, 'gcm')
        foldername = sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s/%s/', type, par.model, par.gcm.clim, par.lat_interp);
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s/%s', type, par.model, par.gcm.clim);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s', type, par.model, par.gcm.clim);
    elseif strcmp(type, 'echam')
        foldername = sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s/', type, par.echam.clim, par.lat_interp);
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s', type, par.echam.clim);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s', type, par.echam.clim);
    end
    load(sprintf('%s/grid.mat', prefix)); % read grid data
    load(sprintf('%s/srfc.mat', prefix)); % read surface variable data
    load(sprintf('%s/%s/masks.mat', prefix_proc, par.lat_interp)); % load land and ocean masks

    if strcmp(par.lat_interp, 'std')
        lat = par.lat_std;
    else
        lat = grid.dim3.lat;
    end

    for fn = fieldnames(srfc)'
        srfc.(fn{1}) = permute(srfc.(fn{1}), [2 1 3]); % bring lat to front
        srfc.(fn{1}) = interp1(grid.dim2.lat, srfc.(fn{1}), lat); % interpolate to standard grid
        srfc.(fn{1}) = permute(srfc.(fn{1}), [2 1 3]); % reorder to original dims
        % for l = {'lo', 'l', 'o'}; land = l{1};
        for l = {'lo'}; land = l{1};
            if strcmp(land, 'lo'); srfc_n.(fn{1}).(land) = srfc.(fn{1});
            elseif strcmp(land, 'l'); srfc_n.(fn{1}).(land) = srfc.(fn{1}).*mask.ocean; % filter out ocean
            elseif strcmp(land, 'o'); srfc_n.(fn{1}).(land) = srfc.(fn{1}).*mask.land; % filter out land
            end
        end
    end

    % for l = {'lo', 'l', 'o'}; land = l{1};
    for l = {'lo'}; land = l{1};
        for fn_var = fieldnames(srfc)'
            ma.(land).(fn_var{1}) = squeeze(nanmean(srfc_n.(fn_var{1}).(land), 1)); % zonal average
            ma_t.(land).(fn_var{1}) = squeeze(nanmean(srfc_n.(fn_var{1}).(land), 3)); % time average

        end % end srfc variables loop

        if strcmp(type, 'era5') | strcmp(type, 'erai')
            pb = CmdLineProgressBar("Calculating moist adiabats...");
            for ilat = 1:length(lat);
                pb.print(ilat, length(lat)); % output progress of moist adiabat calculation
                for imon = 1:12;
                    for fn_var = fieldnames(srfc)'
                        ima.(fn_var{1}) = ma.(land).(fn_var{1})(ilat, imon);
                    end
                    ma.(land).ta(ilat, imon, :) = calc_ma_dew(ima, grid.dim3.plev, par, type); % compute moist adiabat with RH
                    maz.(land).ta(ilat, imon, :) = calc_maz_dew(ima, grid.dim3.z, par, type); % compute moist adiabat with RH
                    masi.(land).ta(ilat, imon, :) = calc_ma_dew_si(ima, grid.dim3.plev, par, type, grid); % compute moist adiabat with RH
                end
                for ilon = 1:length(grid.dim3.lon);
                    for fn_var = fieldnames(srfc)'
                        ima_t.(fn_var{1}) = ma_t.(land).(fn_var{1})(ilon, ilat);
                    end
                    ma_t.(land).ta(ilat, imon, :) = calc_ma_dew(ima_t, grid.dim3.plev, par, type); % compute moist adiabat with RH
                    maz_t.(land).ta(ilat, imon, :) = calc_maz_dew(ima_t, grid.dim3.z, par, type); % compute moist adiabat with RH
                    masi_t.(land).ta(ilat, imon, :) = calc_ma_dew_si(ima_t, grid.dim3.plev, par, type, grid); % compute moist adiabat with RH
                end
            end
        elseif strcmp(type, 'gcm')
            pb = CmdLineProgressBar("Calculating moist adiabats...");
            for ilat = 1:length(lat);
                pb.print(ilat, length(lat)); % output progress of moist adiabat calculation
                for imon = 1:12;
                    for fn_var = fieldnames(srfc)'
                        ima.(fn_var{1}) = ma.(land).(fn_var{1})(ilat, imon);
                    end
                    ma.(land).ta(ilat, imon, :) = calc_ma_hurs(ima, grid.dim3.plev, par); % compute moist adiabat with RH
                    maz.(land).ta(ilat, imon, :) = calc_maz_hurs(ima, grid.dim3.z, par); % compute moist adiabat with RH
                    masi.(land).ta(ilat, imon, :) = calc_ma_hurs_si(ima, grid.dim3.plev, par, grid); % compute moist adiabat with RH
                end
                for ilon = 1:length(grid.dim3.lon);
                    for fn_var = fieldnames(srfc)'
                        ima_t.(fn_var{1}) = ma_t.(land).(fn_var{1})(ilon, ilat);
                    end
                    ma_t.(land).ta(ilat, imon, :) = calc_ma_hurs(ima_t, grid.dim3.plev, par); % compute moist adiabat with RH
                    maz_t.(land).ta(ilat, imon, :) = calc_maz_hurs(ima_t, grid.dim3.z, par); % compute moist adiabat with RH
                    masi_t.(land).ta(ilat, imon, :) = calc_ma_hurs_si(ima_t, grid.dim3.plev, par, grid); % compute moist adiabat with RH
                end
            end
        elseif strcmp(type, 'echam')
            pb = CmdLineProgressBar("Calculating moist adiabats...");
            for ilat = 1:length(lat);
                pb.print(ilat, length(lat)); % output progress of moist adiabat calculation
                for imon = 1:12;
                    for fn_var = fieldnames(srfc)'
                        ima.(fn_var{1}) = ma.(land).(fn_var{1})(ilat, imon);
                    end
                    ma.(land).ta(ilat, imon, :) = calc_ma_dew(ima, grid.dim3.plev, par, type); % compute moist adiabat with RH
                    maz.(land).ta(ilat, imon, :) = calc_maz_dew(ima, grid.dim3.z, par, type); % compute moist adiabat with RH
                    masi.(land).ta(ilat, imon, :) = calc_ma_dew_si(ima, grid.dim3.plev, par, type, grid); % compute moist adiabat with RH
                end
                for ilon = 1:length(grid.dim3.lon);
                    for fn_var = fieldnames(srfc)'
                        ima_t.(fn_var{1}) = ma_t.(land).(fn_var{1})(ilon, ilat);
                    end
                    ma_t.(land).ta(ilat, imon, :) = calc_ma_dew(ima_t, grid.dim3.plev, par, type); % compute moist adiabat with RH
                    maz_t.(land).ta(ilat, imon, :) = calc_maz_dew(ima_t, grid.dim3.z, par, type); % compute moist adiabat with RH
                    masi_t.(land).ta(ilat, imon, :) = calc_ma_dew_si(ima_t, grid.dim3.plev, par, type, grid); % compute moist adiabat with RH
                end
            end
        end
    end % end land option loop

    % save data into mat file
    printname = [foldername 'ma_mon_lat.mat'];
    if ~exist(foldername, 'dir')
        mkdir(foldername)
    end
    save(printname, 'ma', 'maz', 'masi');

    printname = [foldername 'ma_lon_lat.mat'];
    if ~exist(foldername, 'dir')
        mkdir(foldername)
    end
    save(printname, 'ma_t', 'maz_t', 'masi_t');

end % compute mon x lat moist adiabat field
function proc_inv_str_mon_lat(type, par)
    if strcmp(type, 'era5') | strcmp(type, 'erai')
        foldername = sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/', type, par.lat_interp);
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s', type);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s', type);
        file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/%s/temp/%s_temp_%s.ymonmean.nc', type, type, par.(type).yr_span));
        fullpath=sprintf('%s/%s', file.folder, file.name);
        ta_orig = ncread(fullpath, 't');
    elseif strcmp(type, 'gcm')
        foldername = sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s/%s/', type, par.model, par.gcm.clim, par.lat_interp);
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s/%s', type, par.model, par.gcm.clim);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s', type, par.model, par.gcm.clim);
        var = 'ta';
        file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/gcm/%s/%s_Amon_%s_%s_r1i1p1_*.ymonmean.nc', par.model, var, par.model, par.gcm.clim));
        fullpath=sprintf('%s/%s', file.folder, file.name);
        ta_orig = ncread(fullpath, var);
    elseif strcmp(type, 'echam')
        foldername = sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s/', type, par.echam.clim, par.lat_interp);
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s', type, par.echam.clim);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s', type, par.echam.clim);
        var = 't';
        file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/echam/ATM_rjg_%s_*.ymonmean.nc', par.echam.clim));
        fullpath=sprintf('%s/%s', file.folder, file.name);
        ta_orig = double(ncread(fullpath, var));
    end
    load(sprintf('%s/grid.mat', prefix)); % read grid data
    load(sprintf('%s/ta_si.mat', prefix));tasi_orig = ta_si; clear ta_si; % read temp in si coordinates
    load(sprintf('%s/%s/masks.mat', prefix_proc, par.lat_interp)); % load land and ocean masks

    if strcmp(par.lat_interp, 'std')
        lat = par.lat_std;
    else
        lat = grid.dim3.lat;
    end

    % interpolate ta to standard lat grid
    tasi_orig = permute(tasi_orig, [2 1 3 4]);
    tasi_orig = interp1(grid.dim3.lat, tasi_orig, lat);
    tasi_orig = permute(tasi_orig, [2 1 3 4]);

    tasi_orig = permute(tasi_orig, [3 1 2 4]); % bring sigma forward

    for isig = 1:length(par.si_eval); sig = par.si_eval(isig);
        inv_orig(:,:,:,isig) = squeeze(interp1(grid.dim3.si, tasi_orig, sig) - tasi_orig(1,:,:,:)); % evaluate difference between surface and si_eval
    end

    % Land/ocean filter 2D variables
    mask.land_vert = repmat(mask.land, [1 1 1 length(par.si_eval)]); % expand land mask to si_eval dim
    mask.ocean_vert = repmat(mask.ocean, [1 1 1 length(par.si_eval)]); % expand ocean mask to si_eval dim

    inv0.lo = inv_orig;
    inv0.l = inv0.lo.*mask.ocean_vert; % filter inv0 with surface mask
    inv0.o = inv0.lo.*mask.land_vert; % filter inv0 with surface mask

    for l = {'lo', 'l', 'o'}; land = l{1}; % over land, over ocean, or both
        inv.(land)= squeeze(nanmean(inv0.(land), 1)); % zonal average
        % take time averages
        for t = {'ann', 'djf', 'jja', 'mam', 'son'}; time = t{1};
            if strcmp(time, 'ann')
                inv_t.(land).(time) = squeeze(nanmean(inv0.(land), 3));
            elseif strcmp(time, 'djf')
                inv_shift.(land) = circshift(inv0.(land), 1, 3);
                inv_t.(land).(time) = squeeze(nanmean(inv_shift.(land)(:,:,1:3,:), 3));
            elseif strcmp(time, 'jja')
                inv_t.(land).(time) = squeeze(nanmean(inv0.(land)(:,:,6:8,:), 3));
            elseif strcmp(time, 'mam')
                inv_t.(land).(time) = squeeze(nanmean(inv0.(land)(:,:,3:5,:), 3));
            elseif strcmp(time, 'son')
                inv_t.(land).(time) = squeeze(nanmean(inv0.(land)(:,:,9:11,:), 3));
            end
        end
    end

    % save filtered data
    printname = [foldername 'inv_mon_lat'];
    if ~exist(foldername, 'dir')
        mkdir(foldername)
    end
    save(printname, 'inv', 'lat');

    printname = [foldername 'inv_lon_lat'];
    if ~exist(foldername, 'dir')
        mkdir(foldername)
    end
    save(printname, 'inv_t', 'lat');
end % compute mon x lat inversion strength field with land/ocean masking
function make_dtdz(type, par) % compute model lapse rate in lon x lat x plev x mon
% calculate moist adiabats
    if strcmp(type, 'era5') | strcmp(type, 'erai')
        foldername = sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/', type, par.lat_interp);
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s', type);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s', type);
    elseif strcmp(type, 'gcm')
        foldername = sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s/%s/', type, par.model, par.gcm.clim, par.lat_interp);
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s/%s', type, par.model, par.gcm.clim);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s', type, par.model, par.gcm.clim);
    elseif strcmp(type, 'echam')
        foldername = sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s/', type, par.echam.clim, par.lat_interp);
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s', type, par.echam.clim);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s', type, par.echam.clim);
    end
    load(sprintf('%s/grid.mat', prefix)); % read grid data
    load(sprintf('%s/tai.mat', prefix)); % read moist adiabat

    dtdz = -1e3*(tai(:,:,2:end,:)-tai(:,:,1:end-1,:))./(zgi(:,:,2:end,:)-zgi(:,:,1:end-1,:)); % lapse rate in K/km

    dtdz = permute(dtdz, [3 1 2 4]); % bring height forward
    dtdz = interp1(1/2*(par.pa(2:end)+par.pa(1:end-1)), dtdz, par.pa);
    dtdz = permute(dtdz, [2 3 1 4]); % bring height back to 3rd

    if any(strcmp(type, {'era5', 'erai'})); newdir=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s', type);
    elseif strcmp(type, 'gcm'); newdir=sprintf('/project2/tas1/miyawaki/projects/002/data/read/gcm/%s/%s', par.model, par.gcm.clim);
    elseif strcmp(type, 'echam'); newdir=sprintf('/project2/tas1/miyawaki/projects/002/data/read/echam/%s', par.echam.clim); end;
    if ~exist(newdir, 'dir'); mkdir(newdir); end
    filename='dtdz.mat';
    save(sprintf('%s/%s', newdir, filename), 'dtdz', '-v7.3');

end
function make_dtdzi(type, par) % compute model lapse rate in lon x lat x plev x mon
% calculate moist adiabats
    if strcmp(type, 'era5') | strcmp(type, 'erai')
        foldername = sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/', type, par.lat_interp);
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s', type);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s', type);
        file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/%s/temp/%s_temp_%s.ymonmean.nc', type, type, par.(type).yr_span));
        fullpath=sprintf('%s/%s', file.folder, file.name);
        temp = ncread(fullpath, 't');
        file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/%s/zg/%s_zg_%s.ymonmean.nc', type, type, par.(type).yr_span));
        fullpath=sprintf('%s/%s', file.folder, file.name);
        zg = ncread(fullpath, 'z'); zg = zg/par.g; % convert from geopotential to height in m
    elseif strcmp(type, 'gcm')
        foldername = sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s/%s/', type, par.model, par.gcm.clim, par.lat_interp);
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s/%s', type, par.model, par.gcm.clim);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s', type, par.model, par.gcm.clim);
        var = 'ta';
        file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/gcm/%s/%s_Amon_%s_%s_r1i1p1_*.ymonmean.nc', par.model, var, par.model, par.gcm.clim));
        fullpath=sprintf('%s/%s', file.folder, file.name);
        temp = ncread(fullpath, var);
        var = 'zg';
        file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/gcm/%s/%s_Amon_%s_%s_r1i1p1_*.ymonmean.nc', par.model, var, par.model, par.gcm.clim));
        fullpath=sprintf('%s/%s', file.folder, file.name);
        zg = ncread(fullpath, var);
    elseif strcmp(type, 'echam')
        foldername = sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s/', type, par.echam.clim, par.lat_interp);
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s', type, par.echam.clim);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s', type, par.echam.clim);
        var = 't';
        file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/echam/ATM_rjg_%s_*.ymonmean.nc', par.echam.clim));
        fullpath=sprintf('%s/%s', file.folder, file.name);
        temp = double(ncread(fullpath, var));
        var = 'geopoth';
        file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/echam/ATM_rjg_%s_*.ymonmean.nc', par.echam.clim));
        fullpath=sprintf('%s/%s', file.folder, file.name);
        zg = double(ncread(fullpath, var));
    end
    load(sprintf('%s/grid.mat', prefix)); % read grid data
    load(sprintf('%s/srfc.mat', prefix)); % read surface variable data

    dtdzi = -1e3*(temp(:,:,2:end,:)-temp(:,:,1:end-1,:))./(zg(:,:,2:end,:)-zg(:,:,1:end-1,:)); % lapse rate in K/km

    dtdzi = permute(dtdzi, [3 1 2 4]); % bring height forward
    dtdzi = interp1(1/2*(grid.dim3.plev(2:end)+grid.dim3.plev(1:end-1)), dtdzi, grid.dim3.plev);
    dtdzi = permute(dtdzi, [2 3 1 4]); % bring height back to 3rd

    % create surface mask
    if strcmp(type, 'era5') | strcmp(type, 'erai')
        ps_vert = repmat(srfc.sp, [1 1 1 size(dtdzi, 3)]); % dims (lon x lat x time x plev)
        ps_vert = permute(ps_vert, [1 2 4 3]); % dims (lon x lat x plev x time)
        pa = double(permute(repmat(grid.dim3.plev, [1 size(srfc.ps)]), [2 3 1 4]));
    elseif strcmp(type, 'gcm')
        ps_vert = repmat(srfc.ps, [1 1 1 size(dtdzi, 3)]); % dims (lon x lat x time x plev)
        ps_vert = permute(ps_vert, [1 2 4 3]); % dims (lon x lat x plev x time)
        pa = permute(repmat(grid.dim3.plev, [1 size(srfc.ps)]), [2 3 1 4]);
    elseif strcmp(type, 'echam')
        ps_vert = repmat(srfc.aps, [1 1 1 size(dtdzi, 3)]); % dims (lon x lat x time x plev)
        ps_vert = permute(ps_vert, [1 2 4 3]); % dims (lon x lat x plev x time)
        pa = permute(repmat(grid.dim3.plev, [1 size(srfc.aps)]), [2 3 1 4]);
    end
    surface_mask = nan(size(dtdzi));
    surface_mask(pa < ps_vert) = 1;

    dtdzi = dtdzi.*surface_mask; % filter dtdzi with surface mask

    dtdzi = permute(dtdzi, [3 1 2 4]); % bring height forward
    dtdzi = interp1(grid.dim3.plev, dtdzi, par.pa, 'spline', nan);
    dtdzi = permute(dtdzi, [2 3 1 4]); % bring height back to 3rd

    if any(strcmp(type, {'era5', 'erai'})); newdir=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s', type);
    elseif strcmp(type, 'gcm'); newdir=sprintf('/project2/tas1/miyawaki/projects/002/data/read/gcm/%s/%s', par.model, par.gcm.clim);
    elseif strcmp(type, 'echam'); newdir=sprintf('/project2/tas1/miyawaki/projects/002/data/read/echam/%s', par.echam.clim); end;
    if ~exist(newdir, 'dir'); mkdir(newdir); end
    filename='dtdzi.mat';
    save(sprintf('%s/%s', newdir, filename), 'dtdzi', '-v7.3');

end
function make_dtdzz(type, par) % compute model lapse rate in lat x plev x mon
% calculate moist adiabats
    if strcmp(type, 'era5') | strcmp(type, 'erai')
        foldername = sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/', type, par.lat_interp);
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s', type);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s', type);
        file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/%s/temp/%s_temp_%s.ymonmean.nc', type, type, par.(type).yr_span));
        fullpath=sprintf('%s/%s', file.folder, file.name);
        temp = ncread(fullpath, 't');
        file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/%s/zg/%s_zg_%s.ymonmean.nc', type, type, par.(type).yr_span));
        fullpath=sprintf('%s/%s', file.folder, file.name);
        zg = ncread(fullpath, 'z'); zg = zg/par.g; % convert from geopotential to height in m
    elseif strcmp(type, 'gcm')
        foldername = sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s/%s/', type, par.model, par.gcm.clim, par.lat_interp);
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s/%s', type, par.model, par.gcm.clim);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s', type, par.model, par.gcm.clim);
        var = 'ta';
        file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/gcm/%s/%s_Amon_%s_%s_r1i1p1_*.ymonmean.nc', par.model, var, par.model, par.gcm.clim));
        fullpath=sprintf('%s/%s', file.folder, file.name);
        temp = ncread(fullpath, var);
        var = 'zg';
        file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/gcm/%s/%s_Amon_%s_%s_r1i1p1_*.ymonmean.nc', par.model, var, par.model, par.gcm.clim));
        fullpath=sprintf('%s/%s', file.folder, file.name);
        zg = ncread(fullpath, var);
    elseif strcmp(type, 'echam')
        foldername = sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s/', type, par.echam.clim, par.lat_interp);
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s', type, par.echam.clim);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s', type, par.echam.clim);
        var = 't';
        file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/echam/ATM_rjg_%s_*.ymonmean.nc', par.echam.clim));
        fullpath=sprintf('%s/%s', file.folder, file.name);
        temp = double(ncread(fullpath, var));
        var = 'geopoth';
        file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/echam/ATM_rjg_%s_*.ymonmean.nc', par.echam.clim));
        fullpath=sprintf('%s/%s', file.folder, file.name);
        zg = double(ncread(fullpath, var));
    end
    load(sprintf('%s/grid.mat', prefix)); % read grid data
    load(sprintf('%s/srfc.mat', prefix)); % read surface variable data

    temp = squeeze(nanmean(temp, 1)); % zonal average
    zg = squeeze(nanmean(zg, 1)); % zonal average

    dtdzz = -1e3*(temp(:,2:end,:)-temp(:,1:end-1,:))./(zg(:,2:end,:)-zg(:,1:end-1,:)); % lapse rate in K/km

    dtdzz = permute(dtdzz, [2 1 3]); % bring height forward
    dtdzz = interp1(1/2*(grid.dim3.plev(2:end)+grid.dim3.plev(1:end-1)), dtdzz, grid.dim3.plev);
    dtdzz = permute(dtdzz, [2 1 3]); % bring height back to 2nd

    % create surface mask
    if strcmp(type, 'era5') | strcmp(type, 'erai')
        ps_vert = repmat(srfc.sp, [1 1 1 size(dtdzz, 2)]); % dims (lon x lat x time x plev)
        ps_vert = permute(ps_vert, [1 2 4 3]); % dims (lon x lat x plev x time)
        pa = double(permute(repmat(grid.dim3.plev, [1 size(srfc.sp)]), [2 3 1 4]));
    elseif strcmp(type, 'gcm')
        ps_vert = repmat(srfc.ps, [1 1 1 size(dtdzz, 2)]); % dims (lon x lat x time x plev)
        ps_vert = permute(ps_vert, [1 2 4 3]); % dims (lon x lat x plev x time)
        pa = permute(repmat(grid.dim3.plev, [1 size(srfc.ps)]), [2 3 1 4]);
    elseif strcmp(type, 'echam')
        ps_vert = repmat(srfc.aps, [1 1 1 size(dtdzz, 2)]); % dims (lon x lat x time x plev)
        ps_vert = permute(ps_vert, [1 2 4 3]); % dims (lon x lat x plev x time)
        pa = permute(repmat(grid.dim3.plev, [1 size(srfc.aps)]), [2 3 1 4]);
    end
    surface_mask = nan(size(pa));
    surface_mask(pa < ps_vert) = 1;

    dtdzz = repmat(dtdzz, [1 1 1 size(pa,1)]); % repeat in longitude
    dtdzz = permute(dtdzz, [4 1 2 3]); % bring lon to 1st

    dtdzz = dtdzz.*surface_mask; % filter dtdzz with surface mask

    dtdzz = permute(dtdzz, [3 1 2 4]); % bring height forward
    dtdzz = interp1(grid.dim3.plev, dtdzz, par.pa, 'spline', nan);
    dtdzz = permute(dtdzz, [2 3 1 4]); % bring height back to 3rd

    if any(strcmp(type, {'era5', 'erai'})); newdir=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s', type);
    elseif strcmp(type, 'gcm'); newdir=sprintf('/project2/tas1/miyawaki/projects/002/data/read/gcm/%s/%s', par.model, par.gcm.clim);
    elseif strcmp(type, 'echam'); newdir=sprintf('/project2/tas1/miyawaki/projects/002/data/read/echam/%s', par.echam.clim); end;
    if ~exist(newdir, 'dir'); mkdir(newdir); end
    filename='dtdzz.mat';
    save(sprintf('%s/%s', newdir, filename), 'dtdzz', '-v7.3');

end
function make_dtdzsi_old(type, par) % compute model lapse rate in lon x lat x plev x mon
% calculate moist adiabats
    if strcmp(type, 'era5') | strcmp(type, 'erai')
        foldername = sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/', type, par.lat_interp);
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s', type);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s', type);
        file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/%s/temp/%s_temp_%s.ymonmean.nc', type, type, par.(type).yr_span));
        fullpath=sprintf('%s/%s', file.folder, file.name);
        temp = ncread(fullpath, 't');
        file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/%s/zg/%s_zg_%s.ymonmean.nc', type, type, par.(type).yr_span));
        fullpath=sprintf('%s/%s', file.folder, file.name);
        zg = ncread(fullpath, 'z'); zg = zg/par.g; % convert from geopotential to height in m
    elseif strcmp(type, 'gcm')
        foldername = sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s/%s/', type, par.model, par.gcm.clim, par.lat_interp);
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s/%s', type, par.model, par.gcm.clim);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s', type, par.model, par.gcm.clim);
        var = 'ta';
        file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/gcm/%s/%s_Amon_%s_%s_r1i1p1_*.ymonmean.nc', par.model, var, par.model, par.gcm.clim));
        fullpath=sprintf('%s/%s', file.folder, file.name);
        temp = ncread(fullpath, var);
        var = 'zg';
        file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/gcm/%s/%s_Amon_%s_%s_r1i1p1_*.ymonmean.nc', par.model, var, par.model, par.gcm.clim));
        fullpath=sprintf('%s/%s', file.folder, file.name);
        zg = ncread(fullpath, var);
    elseif strcmp(type, 'echam')
        foldername = sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s/', type, par.echam.clim, par.lat_interp);
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s', type, par.echam.clim);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s', type, par.echam.clim);
        var = 't';
        file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/echam/ATM_rjg_%s_*.ymonmean.nc', par.echam.clim));
        fullpath=sprintf('%s/%s', file.folder, file.name);
        temp = double(ncread(fullpath, var));
        var = 'geopoth';
        file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/echam/ATM_rjg_%s_*.ymonmean.nc', par.echam.clim));
        fullpath=sprintf('%s/%s', file.folder, file.name);
        zg = double(ncread(fullpath, var));
    end
    load(sprintf('%s/grid.mat', prefix)); % read grid data
    load(sprintf('%s/srfc.mat', prefix)); % read surface variable data
    % load(sprintf('%s/orog.mat', prefix)); % read surface orography

    % create surface mask
    if strcmp(type, 'era5') | strcmp(type, 'erai')
        ps_vert = repmat(srfc.sp, [1 1 1 size(temp, 3)]); % dims (lon x lat x time x plev)
        ps_vert = permute(ps_vert, [1 2 4 3]); % dims (lon x lat x plev x time)
        pa = double(permute(repmat(grid.dim3.plev, [1 size(srfc.sp)]), [2 3 1 4]));
    elseif strcmp(type, 'gcm')
        ps_vert = repmat(srfc.ps, [1 1 1 size(temp, 3)]); % dims (lon x lat x time x plev)
        ps_vert = permute(ps_vert, [1 2 4 3]); % dims (lon x lat x plev x time)
        pa = permute(repmat(grid.dim3.plev, [1 size(srfc.ps)]), [2 3 1 4]);
    elseif strcmp(type, 'echam')
        ps_vert = repmat(srfc.aps, [1 1 1 size(temp, 3)]); % dims (lon x lat x time x plev)
        ps_vert = permute(ps_vert, [1 2 4 3]); % dims (lon x lat x plev x time)
        pa = permute(repmat(grid.dim3.plev, [1 size(srfc.aps)]), [2 3 1 4]);
    end
    surface_mask = nan(size(pa));
    surface_mask(pa < ps_vert) = 1;

    ta_sm = temp .* surface_mask;
    zg_sm = zg .* surface_mask;

    % add tsurf data and interpolate to higher resolution vertical grid
    [pa_plus ta_plus zg_plus] = deal(nan([size(pa,1), size(pa,2) size(pa,3)+1 size(pa,4)])); % create empty grid with one extra vertical level
    pa_plus(:,:,1:end-1,:) = pa; % populate with standard pressure grid
    ta_plus(:,:,1:end-1,:) = ta_sm; % populate with standard atmospheric temperature
    zg_plus(:,:,1:end-1,:) = zg_sm; % populate with szgndard atmospheric temperature
    pa_plus(:,:,end,:) = ps_vert(:,:,1,:); % add surface pressure data into standard pressure grid
    if any(strcmp(type, {'era5', 'erai'})); ta_plus(:,:,end,:) = srfc.t2m(:,:,:); % add surface temperature data
    elseif strcmp(type, 'gcm'); ta_plus(:,:,end,:) = srfc.tas(:,:,:); % add surface temperature data
    elseif strcmp(type, 'echam'); ta_plus(:,:,end,:) = srfc.temp2(:,:,:); end % add surface temperature data
    % zg_plus(:,:,end,:) = srfc.zs(:,:,:); % add surface height data
    % zg_plus(:,:,end,:) = repmat(orog, [1 1 12]); % add surface height data
    zg_plus(:,:,end,:) = nan(size(srfc.zs)); % add surface height data
    ps_vert = permute(ps_vert, [3 1 2 4]); % bring plev dimension to front
    pa_plus = permute(pa_plus, [3 1 2 4]); % bring plev dimension to front
    ta_plus = permute(ta_plus, [3 1 2 4]); % bring plev dimension to front
    zg_plus = permute(zg_plus, [3 1 2 4]); % bring plev dimension to front
    [pa_plus sort_index] = sort(pa_plus, 1, 'descend'); % sort added surface pressure such that pressure decreases monotonically
    pb = CmdLineProgressBar("Sorting temperature with surface data added...");
    for lo=1:size(pa_plus,2)
        pb.print(lo, size(pa_plus,2));
        for la=1:size(pa_plus,3)
            for mo=1:size(pa_plus,4)
                ta_plus(:,lo,la,mo) = ta_plus(sort_index(:,lo,la,mo),lo,la,mo); % sort temperature (has to be in loop because sort_index works for vector calls only)
                zg_plus(:,lo,la,mo) = zg_plus(sort_index(:,lo,la,mo),lo,la,mo); % sort temperature (has to be in loop because sort_index works for vector calls only)
            end
        end
    end

    dtdz = -1e3*(ta_plus(2:end,:,:,:)-ta_plus(1:end-1,:,:,:))./(zg_plus(2:end,:,:,:)-zg_plus(1:end-1,:,:,:)); % lapse rate in K/km

    pb = CmdLineProgressBar("Sorting and interpolating dtdz to new standard grid...");
    for lo=1:size(pa_plus,2)
        pb.print(lo, size(pa_plus,2));
        for la=1:size(pa_plus,3)
            for mo=1:size(pa_plus,4)
                dtdzsi(:,lo,la,mo) = interp1(1/2*(pa_plus(2:end,lo,la,mo)+pa_plus(1:end-1,lo,la,mo))./ps_vert(:,lo,la,mo), dtdz(:,lo,la,mo), grid.dim3.si, 'linear');
            end
        end
    end
    dtdzsi = permute(dtdzsi, [2 3 1 4]); % bring height to 3rd
    dtdzsi = fillmissing(dtdzsi, 'nearest');

    if any(strcmp(type, {'era5', 'erai'})); newdir=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s', type);
    elseif strcmp(type, 'gcm'); newdir=sprintf('/project2/tas1/miyawaki/projects/002/data/read/gcm/%s/%s', par.model, par.gcm.clim);
    elseif strcmp(type, 'echam'); newdir=sprintf('/project2/tas1/miyawaki/projects/002/data/read/echam/%s', par.echam.clim); end;
    if ~exist(newdir, 'dir'); mkdir(newdir); end
    filename='dtdzsi.mat';
    save(sprintf('%s/%s', newdir, filename), 'dtdzsi', '-v7.3');

end
function proc_ga_diff_mon_lat(type, par)
    if strcmp(type, 'era5') | strcmp(type, 'erai')
        foldername = sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/', type, par.lat_interp);
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s', type);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s', type);
        file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/%s/temp/%s_temp_%s.ymonmean.nc', type, type, par.(type).yr_span));
        fullpath=sprintf('%s/%s', file.folder, file.name);
        ta_orig = ncread(fullpath, 't');
    elseif strcmp(type, 'gcm')
        foldername = sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s/%s/', type, par.model, par.gcm.clim, par.lat_interp);
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s/%s', type, par.model, par.gcm.clim);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s', type, par.model, par.gcm.clim);
        var = 'ta';
        file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/gcm/%s/%s_Amon_%s_%s_r1i1p1_*.ymonmean.nc', par.model, var, par.model, par.gcm.clim));
        fullpath=sprintf('%s/%s', file.folder, file.name);
        ta_orig = ncread(fullpath, var);
    elseif strcmp(type, 'echam')
        foldername = sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s/', type, par.echam.clim, par.lat_interp);
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s', type, par.echam.clim);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s', type, par.echam.clim);
        var = 't';
        file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/echam/ATM_rjg_%s_*.ymonmean.nc', par.echam.clim));
        fullpath=sprintf('%s/%s', file.folder, file.name);
        ta_orig = double(ncread(fullpath, var));
    end
    load(sprintf('%s/grid.mat', prefix)); % read grid data
    load(sprintf('%s/dtdz.mat', prefix)); % read temp in si coordinates
    load(sprintf('%s/dtmdz.mat', prefix)); % read temp in si coordinates
    load(sprintf('%s/%s/masks.mat', prefix_proc, par.lat_interp)); % load land and ocean masks

    if strcmp(par.lat_interp, 'std')
        lat = par.lat_std;
    else
        lat = grid.dim3.lat;
    end

    ga_diff_orig = (dtmdz - dtdz)./dtmdz * 1e2; % percentage difference of moist adiabatic vs GCM lapse rates
    clear dtmdz dtdz;
    ga_diff_orig = permute(ga_diff_orig,[2 1 3 4]); % bring lat to 1st
    ga_diff_orig = interp1(grid.dim3.lat, ga_diff_orig, lat); % interpolate to standard lat
    ga_diff_orig = permute(ga_diff_orig, [3 2 1 4]); % bring height front
    ga_diff_orig = interp1(par.pa, ga_diff_orig, 1e2*linspace(1000,200,100)); % prepare to average between 1000-200 hPa
    ga_diff_orig = squeeze(nanmean(ga_diff_orig,1)); % take vertical average

    ga_diff0.lo = ga_diff_orig;
    ga_diff0.l = ga_diff0.lo.*mask.ocean; % filter ga_diff0 with surface mask
    ga_diff0.o = ga_diff0.lo.*mask.land; % filter ga_diff0 with surface mask

    for l = {'lo', 'l', 'o'}; land = l{1}; % over land, over ocean, or both
        ga_diff.(land)= squeeze(nanmean(ga_diff0.(land), 1)); % zonal average
        % take time averages
        for t = {'ann', 'djf', 'jja', 'mam', 'son'}; time = t{1};
            if strcmp(time, 'ann')
                ga_diff_t.(land).(time) = squeeze(nanmean(ga_diff0.(land), 3));
            elseif strcmp(time, 'djf')
                ga_diff_shift.(land) = circshift(ga_diff0.(land), 1, 3);
                ga_diff_t.(land).(time) = squeeze(nanmean(ga_diff_shift.(land)(:,:,1:3,:), 3));
            elseif strcmp(time, 'jja')
                ga_diff_t.(land).(time) = squeeze(nanmean(ga_diff0.(land)(:,:,6:8,:), 3));
            elseif strcmp(time, 'mam')
                ga_diff_t.(land).(time) = squeeze(nanmean(ga_diff0.(land)(:,:,3:5,:), 3));
            elseif strcmp(time, 'son')
                ga_diff_t.(land).(time) = squeeze(nanmean(ga_diff0.(land)(:,:,9:11,:), 3));
            end
        end
    end

    % save filtered data
    printname = [foldername 'ga_diff_mon_lat'];
    if ~exist(foldername, 'dir')
        mkdir(foldername)
    end
    save(printname, 'ga_diff', 'lat');

    printname = [foldername 'ga_diff_lon_lat'];
    if ~exist(foldername, 'dir')
        mkdir(foldername)
    end
    save(printname, 'ga_diff_t', 'lat');
end % compute mon x lat gamma percentage difference field with land/ocean masking
function proc_ga_diff_si_mon_lat(type, par)
    if strcmp(type, 'era5') | strcmp(type, 'erai')
        foldername = sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/', type, par.lat_interp);
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s', type);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s', type);
        file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/%s/temp/%s_temp_%s.ymonmean.nc', type, type, par.(type).yr_span));
        fullpath=sprintf('%s/%s', file.folder, file.name);
        ta_orig = ncread(fullpath, 't');
    elseif strcmp(type, 'gcm')
        foldername = sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s/%s/', type, par.model, par.gcm.clim, par.lat_interp);
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s/%s', type, par.model, par.gcm.clim);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s', type, par.model, par.gcm.clim);
        var = 'ta';
        file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/gcm/%s/%s_Amon_%s_%s_r1i1p1_*.ymonmean.nc', par.model, var, par.model, par.gcm.clim));
        fullpath=sprintf('%s/%s', file.folder, file.name);
        ta_orig = ncread(fullpath, var);
    elseif strcmp(type, 'echam')
        foldername = sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s/', type, par.echam.clim, par.lat_interp);
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s', type, par.echam.clim);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s', type, par.echam.clim);
        var = 't';
        file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/echam/ATM_rjg_%s_*.ymonmean.nc', par.echam.clim));
        fullpath=sprintf('%s/%s', file.folder, file.name);
        ta_orig = double(ncread(fullpath, var));
    end
    load(sprintf('%s/grid.mat', prefix)); % read grid data
    load(sprintf('%s/dtdzzsi.mat', prefix)); % read temp in si coordinates
    load(sprintf('%s/dtmdzzsi.mat', prefix)); % read temp in si coordinates
    load(sprintf('%s/%s/masks.mat', prefix_proc, par.lat_interp)); % load land and ocean masks

    if strcmp(par.lat_interp, 'std')
        lat = par.lat_std;
    else
        lat = grid.dim3.lat;
    end

    ga_diff_orig = (dtmdzzsi - dtdzzsi)./dtmdzzsi * 1e2; % percentage difference of moist adiabatic vs GCM lapse rates
    clear dtmdz dtdz;
    ga_diff_orig = permute(ga_diff_orig,[2 1 3 4]); % bring lat to 1st
    ga_diff_orig = interp1(grid.dim3.lat, ga_diff_orig, lat); % interpolate to standard lat
    ga_diff_orig = permute(ga_diff_orig, [3 2 1 4]); % bring height front
    ga_diff_orig = interp1(grid.dim3.si, ga_diff_orig, linspace(par.si_bl,par.si_up,100)); % prepare to average between 1000-200 hPa
    ga_diff_orig = squeeze(nanmean(ga_diff_orig,1)); % take vertical average

    ga_diff0.lo = ga_diff_orig;
    ga_diff0.l = ga_diff0.lo.*mask.ocean; % filter ga_diff0 with surface mask
    ga_diff0.o = ga_diff0.lo.*mask.land; % filter ga_diff0 with surface mask

    for l = {'lo', 'l', 'o'}; land = l{1}; % over land, over ocean, or both
        ga_diff.(land)= squeeze(nanmean(ga_diff0.(land), 1)); % zonal average
        % take time averages
        for t = {'ann', 'djf', 'jja', 'mam', 'son'}; time = t{1};
            if strcmp(time, 'ann')
                ga_diff_t.(land).(time) = squeeze(nanmean(ga_diff0.(land), 3));
            elseif strcmp(time, 'djf')
                ga_diff_shift.(land) = circshift(ga_diff0.(land), 1, 3);
                ga_diff_t.(land).(time) = squeeze(nanmean(ga_diff_shift.(land)(:,:,1:3,:), 3));
            elseif strcmp(time, 'jja')
                ga_diff_t.(land).(time) = squeeze(nanmean(ga_diff0.(land)(:,:,6:8,:), 3));
            elseif strcmp(time, 'mam')
                ga_diff_t.(land).(time) = squeeze(nanmean(ga_diff0.(land)(:,:,3:5,:), 3));
            elseif strcmp(time, 'son')
                ga_diff_t.(land).(time) = squeeze(nanmean(ga_diff0.(land)(:,:,9:11,:), 3));
            end
        end
    end

    % save filtered data
    printname = [foldername 'ga_diff_si_mon_lat'];
    if ~exist(foldername, 'dir')
        mkdir(foldername)
    end
    save(printname, 'ga_diff', 'lat');

    printname = [foldername 'ga_diff_si_lon_lat'];
    if ~exist(foldername, 'dir')
        mkdir(foldername)
    end
    save(printname, 'ga_diff_t', 'lat');
end % compute mon x lat gamma percentage difference field with land/ocean masking
function proc_ga_bl_diff_si_mon_lat(type, par)
    if strcmp(type, 'era5') | strcmp(type, 'erai')
        foldername = sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/', type, par.lat_interp);
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s', type);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s', type);
        file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/%s/temp/%s_temp_%s.ymonmean.nc', type, type, par.(type).yr_span));
        fullpath=sprintf('%s/%s', file.folder, file.name);
        ta_orig = ncread(fullpath, 't');
    elseif strcmp(type, 'gcm')
        foldername = sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s/%s/', type, par.model, par.gcm.clim, par.lat_interp);
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s/%s', type, par.model, par.gcm.clim);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s', type, par.model, par.gcm.clim);
        var = 'ta';
        file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/gcm/%s/%s_Amon_%s_%s_r1i1p1_*.ymonmean.nc', par.model, var, par.model, par.gcm.clim));
        fullpath=sprintf('%s/%s', file.folder, file.name);
        ta_orig = ncread(fullpath, var);
    elseif strcmp(type, 'echam')
        foldername = sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s/', type, par.echam.clim, par.lat_interp);
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s', type, par.echam.clim);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s', type, par.echam.clim);
        var = 't';
        file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/echam/ATM_rjg_%s_*.ymonmean.nc', par.echam.clim));
        fullpath=sprintf('%s/%s', file.folder, file.name);
        ta_orig = double(ncread(fullpath, var));
    end
    load(sprintf('%s/grid.mat', prefix)); % read grid data
    load(sprintf('%s/dtdzzsi.mat', prefix)); % read temp in si coordinates
    load(sprintf('%s/dtmdzzsi.mat', prefix)); % read temp in si coordinates
    load(sprintf('%s/%s/masks.mat', prefix_proc, par.lat_interp)); % load land and ocean masks

    if strcmp(par.lat_interp, 'std')
        lat = par.lat_std;
    else
        lat = grid.dim3.lat;
    end

    ga_bl_diff_orig = (dtmdzzsi - dtdzzsi)./dtmdzzsi * 1e2; % percentage difference of moist adiabatic vs GCM lapse rates
    clear dtmdz dtdz;
    ga_bl_diff_orig = permute(ga_bl_diff_orig,[2 1 3 4]); % bring lat to 1st
    ga_bl_diff_orig = interp1(grid.dim3.lat, ga_bl_diff_orig, lat); % interpolate to standard lat
    ga_bl_diff_orig = permute(ga_bl_diff_orig, [3 2 1 4]); % bring height front
    ga_bl_diff_orig = interp1(grid.dim3.si, ga_bl_diff_orig, linspace(1,par.si_bl,100)); % prepare to average between 1000-200 hPa
    ga_bl_diff_orig = squeeze(nanmean(ga_bl_diff_orig,1)); % take vertical average

    ga_bl_diff0.lo = ga_bl_diff_orig;
    ga_bl_diff0.l = ga_bl_diff0.lo.*mask.ocean; % filter ga_bl_diff0 with surface mask
    ga_bl_diff0.o = ga_bl_diff0.lo.*mask.land; % filter ga_bl_diff0 with surface mask

    for l = {'lo', 'l', 'o'}; land = l{1}; % over land, over ocean, or both
        ga_bl_diff.(land)= squeeze(nanmean(ga_bl_diff0.(land), 1)); % zonal average
        % take time averages
        for t = {'ann', 'djf', 'jja', 'mam', 'son'}; time = t{1};
            if strcmp(time, 'ann')
                ga_bl_diff_t.(land).(time) = squeeze(nanmean(ga_bl_diff0.(land), 3));
            elseif strcmp(time, 'djf')
                ga_bl_diff_shift.(land) = circshift(ga_bl_diff0.(land), 1, 3);
                ga_bl_diff_t.(land).(time) = squeeze(nanmean(ga_bl_diff_shift.(land)(:,:,1:3,:), 3));
            elseif strcmp(time, 'jja')
                ga_bl_diff_t.(land).(time) = squeeze(nanmean(ga_bl_diff0.(land)(:,:,6:8,:), 3));
            elseif strcmp(time, 'mam')
                ga_bl_diff_t.(land).(time) = squeeze(nanmean(ga_bl_diff0.(land)(:,:,3:5,:), 3));
            elseif strcmp(time, 'son')
                ga_bl_diff_t.(land).(time) = squeeze(nanmean(ga_bl_diff0.(land)(:,:,9:11,:), 3));
            end
        end
    end

    % save filtered data
    printname = [foldername 'ga_bl_diff_si_mon_lat'];
    if ~exist(foldername, 'dir')
        mkdir(foldername)
    end
    save(printname, 'ga_bl_diff', 'lat');

    printname = [foldername 'ga_bl_diff_si_lon_lat'];
    if ~exist(foldername, 'dir')
        mkdir(foldername)
    end
    save(printname, 'ga_bl_diff_t', 'lat');
end % compute mon x lat gamma percentage difference field with land/ocean masking
function proc_ga_dalr_diff_si_mon_lat(type, par)
    if strcmp(type, 'era5') | strcmp(type, 'erai')
        foldername = sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/', type, par.lat_interp);
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s', type);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s', type);
        file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/%s/temp/%s_temp_%s.ymonmean.nc', type, type, par.(type).yr_span));
        fullpath=sprintf('%s/%s', file.folder, file.name);
        ta_orig = ncread(fullpath, 't');
    elseif strcmp(type, 'gcm')
        foldername = sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s/%s/', type, par.model, par.gcm.clim, par.lat_interp);
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s/%s', type, par.model, par.gcm.clim);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s', type, par.model, par.gcm.clim);
        var = 'ta';
        file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/gcm/%s/%s_Amon_%s_%s_r1i1p1_*.ymonmean.nc', par.model, var, par.model, par.gcm.clim));
        fullpath=sprintf('%s/%s', file.folder, file.name);
        ta_orig = ncread(fullpath, var);
    elseif strcmp(type, 'echam')
        foldername = sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s/', type, par.echam.clim, par.lat_interp);
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s', type, par.echam.clim);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s', type, par.echam.clim);
        var = 't';
        file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/echam/ATM_rjg_%s_*.ymonmean.nc', par.echam.clim));
        fullpath=sprintf('%s/%s', file.folder, file.name);
        ta_orig = double(ncread(fullpath, var));
    end
    load(sprintf('%s/grid.mat', prefix)); % read grid data
    load(sprintf('%s/dtdzzsi.mat', prefix)); % read temp in si coordinates
    load(sprintf('%s/%s/masks.mat', prefix_proc, par.lat_interp)); % load land and ocean masks

    if strcmp(par.lat_interp, 'std')
        lat = par.lat_std;
    else
        lat = grid.dim3.lat;
    end

    dtmdzzsi = 1e3*par.g/par.cpd*ones(size(dtdzzsi));

    ga_dalr_diff_orig = (dtmdzzsi - dtdzzsi)./dtmdzzsi * 1e2; % percentage difference of moist adiabatic vs GCM lapse rates
    clear dtmdz dtdz;
    ga_dalr_diff_orig = permute(ga_dalr_diff_orig,[2 1 3 4]); % bring lat to 1st
    ga_dalr_diff_orig = interp1(grid.dim3.lat, ga_dalr_diff_orig, lat); % interpolate to standard lat
    ga_dalr_diff_orig = permute(ga_dalr_diff_orig, [3 2 1 4]); % bring height front
    ga_dalr_diff_orig = interp1(grid.dim3.si, ga_dalr_diff_orig, linspace(par.si_bl,par.si_up,100)); % prepare to average between 1000-200 hPa
    ga_dalr_diff_orig = squeeze(nanmean(ga_dalr_diff_orig,1)); % take vertical average

    ga_dalr_diff0.lo = ga_dalr_diff_orig;
    ga_dalr_diff0.l = ga_dalr_diff0.lo.*mask.ocean; % filter ga_dalr_diff0 with surface mask
    ga_dalr_diff0.o = ga_dalr_diff0.lo.*mask.land; % filter ga_dalr_diff0 with surface mask

    for l = {'lo', 'l', 'o'}; land = l{1}; % over land, over ocean, or both
        ga_dalr_diff.(land)= squeeze(nanmean(ga_dalr_diff0.(land), 1)); % zonal average
        % take time averages
        for t = {'ann', 'djf', 'jja', 'mam', 'son'}; time = t{1};
            if strcmp(time, 'ann')
                ga_dalr_diff_t.(land).(time) = squeeze(nanmean(ga_dalr_diff0.(land), 3));
            elseif strcmp(time, 'djf')
                ga_dalr_diff_shift.(land) = circshift(ga_dalr_diff0.(land), 1, 3);
                ga_dalr_diff_t.(land).(time) = squeeze(nanmean(ga_dalr_diff_shift.(land)(:,:,1:3,:), 3));
            elseif strcmp(time, 'jja')
                ga_dalr_diff_t.(land).(time) = squeeze(nanmean(ga_dalr_diff0.(land)(:,:,6:8,:), 3));
            elseif strcmp(time, 'mam')
                ga_dalr_diff_t.(land).(time) = squeeze(nanmean(ga_dalr_diff0.(land)(:,:,3:5,:), 3));
            elseif strcmp(time, 'son')
                ga_dalr_diff_t.(land).(time) = squeeze(nanmean(ga_dalr_diff0.(land)(:,:,9:11,:), 3));
            end
        end
    end

    % save filtered data
    printname = [foldername 'ga_dalr_diff_si_mon_lat'];
    if ~exist(foldername, 'dir')
        mkdir(foldername)
    end
    save(printname, 'ga_dalr_diff', 'lat');

    printname = [foldername 'ga_dalr_diff_si_lon_lat'];
    if ~exist(foldername, 'dir')
        mkdir(foldername)
    end
    save(printname, 'ga_dalr_diff_t', 'lat');
end % compute mon x lat gamma percentage difference field with land/ocean masking
function proc_temp_mon_lat_interp(type, par)
    if strcmp(type, 'era5') | strcmp(type, 'erai')
        foldername = sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/', type, par.lat_interp);
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s', type);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s', type);
        file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/%s/temp/%s_temp_%s.ymonmean.nc', type, type, par.(type).yr_span));
        fullpath=sprintf('%s/%s', file.folder, file.name);
        ta_orig = ncread(fullpath, 't');
    elseif strcmp(type, 'gcm')
        foldername = sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s/%s/', type, par.model, par.gcm.clim, par.lat_interp);
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s/%s', type, par.model, par.gcm.clim);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s', type, par.model, par.gcm.clim);
        var = 'ta';
        file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/gcm/%s/%s_Amon_%s_%s_r1i1p1_*.ymonmean.nc', par.model, var, par.model, par.gcm.clim));
        fullpath=sprintf('%s/%s', file.folder, file.name);
        ta_orig = ncread(fullpath, var);
    elseif strcmp(type, 'echam')
        foldername = sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s/', type, par.echam.clim, par.lat_interp);
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s', type, par.echam.clim);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s', type, par.echam.clim);
        var = 't';
        if contains(par.echam.clim, 'rp000')
            file=dir(sprintf('/project2/tas1/ockham/data11/tas/echam-aiv_rcc_6.1.00p1/%s/ATM_%s_0020_39.nc', par.echam.clim, par.echam.clim));
        else
            file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/echam/ATM_rjg_%s_*.ymonmean.nc', par.echam.clim));
        end
        fullpath=sprintf('%s/%s', file.folder, file.name);
        ta_orig = double(ncread(fullpath, var));
    elseif strcmp(type, 'echam_ml')
        foldername = sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/', type, par.lat_interp);
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s', type);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s', type);
        var = 't';
        file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/echam_ml/ATM_*.ymonmean.nc'));
        fullpath=sprintf('%s/%s', file.folder, file.name);
        ta_orig = double(ncread(fullpath, var));
    elseif strcmp(type, 'echam_pl')
        foldername = sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/', type, par.lat_interp);
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s', type);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s', type);
        var = 't';
        file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/echam_pl/ATM_*.ymonmean.nc'));
        fullpath=sprintf('%s/%s', file.folder, file.name);
        ta_orig = double(ncread(fullpath, var));
    end
    load(sprintf('%s/grid.mat', prefix)); % read grid data
    load(sprintf('%s/srfc.mat', prefix)); % load surface data
    load(sprintf('%s/%s/masks.mat', prefix_proc, par.lat_interp)); % load land and ocean masks

    if strcmp(par.lat_interp, 'std')
        lat = par.lat_std;
    else
        lat = grid.dim3.lat;
    end


    for i = {'lin', 'cub', 'spl', 'mak'}; itp=i{1};
        clear tasi_orig pasi_orig

        load(sprintf('%s/ta_si.mat', prefix)); tasi_orig = ta_si.(itp); clear ta_si; % read temp in si coordinates
        load(sprintf('%s/pa_si.mat', prefix)); pasi_orig = pa_si; clear pa_si; % read temp in si coordinates

        % interpolate ta to standard lat grid
        tasi_orig = permute(tasi_orig, [2 1 3 4]);
        tasi_orig = interp1(grid.dim3.lat, tasi_orig, lat);
        tasi_orig = permute(tasi_orig, [2 1 3 4]);

        pasi_orig = permute(pasi_orig, [2 1 3 4]);
        pasi_orig = interp1(grid.dim3.lat, pasi_orig, lat);
        pasi_orig = permute(pasi_orig, [2 1 3 4]);

        % create surface mask
        if strcmp(type, 'era5') | strcmp(type, 'erai')
            ps = permute(srfc.sp, [2 1 3]); % bring lat to front to interpolate
            ps = interp1(grid.dim2.lat, ps, lat, 'linear', 'extrap'); % interpolate to standard grid
            ps = permute(ps, [2 1 3]); % reorder to original
            ps_vert = repmat(ps, [1 1 1 size(ta_orig, 3)]); % dims (lon x lat x time x plev)
            ps_vert = permute(ps_vert, [1 2 4 3]); % dims (lon x lat x plev x time)
            pa = double(permute(repmat(grid.dim3.plev, [1 size(ps)]), [2 3 1 4]));
            ts = permute(srfc.t2m, [2 1 3]); % bring lat to front to interpolate
            ts = interp1(grid.dim2.lat, ts, lat); % interpolate to standard grid
            ts = permute(ts, [2 1 3]); % reorder to original
        elseif strcmp(type, 'gcm')
            ps = permute(srfc.ps, [2 1 3]); % bring lat to front to interpolate
            ps = interp1(grid.dim2.lat, ps, lat, 'linear', 'extrap'); % interpolate to standard grid
            ps = permute(ps, [2 1 3]); % reorder to original
            ps_vert = repmat(ps, [1 1 1 size(ta_orig, 3)]); % dims (lon x lat x time x plev)
            ps_vert = permute(ps_vert, [1 2 4 3]); % dims (lon x lat x plev x time)
            pa = permute(repmat(grid.dim3.plev, [1 size(ps)]), [2 3 1 4]);
            ts = permute(srfc.tas, [2 1 3]); % bring lat to front to interpolate
            ts = interp1(grid.dim2.lat, ts, lat); % interpolate to standard grid
            ts = permute(ts, [2 1 3]); % reorder to original
        elseif strcmp(type, 'echam') | strcmp(type, 'echam_ml') | strcmp(type, 'echam_pl')
            ps = permute(srfc.aps, [2 1 3]); % bring lat to front to interpolate
            ps = interp1(grid.dim2.lat, ps, lat, 'linear', 'extrap'); % interpolate to standard grid
            ps = permute(ps, [2 1 3]); % reorder to original
            ps_vert = repmat(ps, [1 1 1 size(ta_orig, 3)]); % dims (lon x lat x time x plev)
            ps_vert = permute(ps_vert, [1 2 4 3]); % dims (lon x lat x plev x time)
            pa = permute(repmat(grid.dim3.plev, [1 size(ps)]), [2 3 1 4]);
            ts = permute(srfc.temp2, [2 1 3]); % bring lat to front to interpolate
            ts = interp1(grid.dim2.lat, ts, lat); % interpolate to standard grid
            ts = permute(ts, [2 1 3]); % reorder to original
        end

        if strcmp(type, 'echam_ml')
            sm = ones(size(ta_orig));
        else
            sm = nan(size(ta_orig));
            sm(pa < ps_vert) = 1;
        end

        tasi_sm.lo = tasi_orig; % surface is already masked in standard sigma coordinates

        pasi_sm.lo = pasi_orig; % surface is already masked in spndard sigma coordinates

        % Land/ocean filter 3D variables
        mask.land_vert = repmat(mask.land, [1 1 1 size(ta_orig, 3)]); % expand land mask to vertical dim
        mask.land_vert = permute(mask.land, [1 2 4 3]); % place vertical dim where it belongs
        mask.ocean_vert = repmat(mask.ocean, [1 1 1 size(ta_orig, 3)]); % expand ocean mask to vertical dim
        mask.ocean_vert = permute(mask.ocean, [1 2 4 3]); % place vertical dim where it belongs

        tasi_sm.l = tasi_sm.lo.*mask.ocean_vert; % filter tasi with surface mask
        tasi_sm.o = tasi_sm.lo.*mask.land_vert; % filter tasi with surface mask
        tasi_sm.lo = permute(tasi_sm.lo, [1 2 4 3]); % bring plev to last dimension
        tasi_sm.l = permute(tasi_sm.l, [1 2 4 3]); % bring plev to last dimension
        tasi_sm.o = permute(tasi_sm.o, [1 2 4 3]); % bring plev to last dimension

        pasi_sm.l = pasi_sm.lo.*mask.ocean_vert; % filter pasi with surface mask
        pasi_sm.o = pasi_sm.lo.*mask.land_vert; % filter pasi with surface mask
        pasi_sm.lo = permute(pasi_sm.lo, [1 2 4 3]); % bring plev to last dimension
        pasi_sm.l = permute(pasi_sm.l, [1 2 4 3]); % bring plev to last dimension
        pasi_sm.o = permute(pasi_sm.o, [1 2 4 3]); % bring plev to last dimension

        mask_t.land = nanmean(mask.land, 3);
        mask_t.ocean = nanmean(mask.ocean, 3);

        for l = {'lo', 'l', 'o'}; land = l{1}; % over land, over ocean, or both
            tasi.(land)= squeeze(nanmean(tasi_sm.(land), 1)); % zonal average
            pasi.(land)= squeeze(nanmean(pasi_sm.(land), 1)); % zonal average
        end

        for l = {'lo', 'l', 'o'}; land = l{1}; % over land, over ocean, or both
            % take time averages
            for t = {'ann', 'djf', 'jja', 'mam', 'son'}; time = t{1};
                if strcmp(time, 'ann')
                    tasi_t.(land).(time) = squeeze(nanmean(tasi_sm.(land), 3));
                    pasi_t.(land).(time) = squeeze(nanmean(pasi_sm.(land), 3));
                elseif strcmp(time, 'djf')
                    tasi_shift.(land) = circshift(tasi_sm.(land), 1, 3);
                    pasi_shift.(land) = circshift(pasi_sm.(land), 1, 3);

                    tasi_t.(land).(time) = squeeze(nanmean(tasi_shift.(land)(:,:,1:3,:), 3));
                    pasi_t.(land).(time) = squeeze(nanmean(pasi_shift.(land)(:,:,1:3,:), 3));

                elseif strcmp(time, 'jja')
                    tasi_t.(land).(time) = squeeze(nanmean(tasi_sm.(land)(:,:,6:8,:), 3));
                    pasi_t.(land).(time) = squeeze(nanmean(pasi_sm.(land)(:,:,6:8,:), 3));

                elseif strcmp(time, 'mam')
                    tasi_t.(land).(time) = squeeze(nanmean(tasi_sm.(land)(:,:,3:5,:), 3));
                    pasi_t.(land).(time) = squeeze(nanmean(pasi_sm.(land)(:,:,3:5,:), 3));

                elseif strcmp(time, 'son')
                    tasi_t.(land).(time) = squeeze(nanmean(tasi_sm.(land)(:,:,9:11,:), 3));
                    pasi_t.(land).(time) = squeeze(nanmean(pasi_sm.(land)(:,:,9:11,:), 3));

                end
            end
        end

        % save filtered data
        printname = [foldername 'ta_mon_lat_' itp];
        if ~exist(foldername, 'dir')
            mkdir(foldername)
        end
        if par.do_surf; save(printname, 'tasi', 'lat');
        else save(printname, 'tasi', 'pasi', 'lat', '-v7.3'); end

        printname = [foldername 'ta_lon_lat_' itp];
        if ~exist(foldername, 'dir')
            mkdir(foldername)
        end
        if par.do_surf; save(printname, 'tasi_t', 'lat');
        else save(printname, 'tasi_t', 'pasi_t', 'lat', '-v7.3'); end

    end
end % compute mon x lat temperature field
function proc_temp_mon_lat_interp_mean(type, par)
    if strcmp(type, 'era5') | strcmp(type, 'erai')
        foldername = sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/', type, par.lat_interp);
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s', type);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s', type);
        file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/%s/temp/%s_temp_%s.ymonmean.nc', type, type, par.(type).yr_span));
        fullpath=sprintf('%s/%s', file.folder, file.name);
        ta_orig = ncread(fullpath, 't');
    elseif strcmp(type, 'gcm')
        foldername = sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s/%s/', type, par.model, par.gcm.clim, par.lat_interp);
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s/%s', type, par.model, par.gcm.clim);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s', type, par.model, par.gcm.clim);
        var = 'ta';
        file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/gcm/%s/%s_Amon_%s_%s_r1i1p1_*.ymonmean.nc', par.model, var, par.model, par.gcm.clim));
        fullpath=sprintf('%s/%s', file.folder, file.name);
        ta_orig = ncread(fullpath, var);
    elseif strcmp(type, 'echam')
        foldername = sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s/', type, par.echam.clim, par.lat_interp);
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s', type, par.echam.clim);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s', type, par.echam.clim);
        var = 't';
        if contains(par.echam.clim, 'rp000')
            file=dir(sprintf('/project2/tas1/ockham/data11/tas/echam-aiv_rcc_6.1.00p1/%s/ATM_%s_0020_39.nc', par.echam.clim, par.echam.clim));
        else
            file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/echam/ATM_rjg_%s_*.ymonmean.nc', par.echam.clim));
        end
        fullpath=sprintf('%s/%s', file.folder, file.name);
        ta_orig = double(ncread(fullpath, var));
    elseif strcmp(type, 'echam_ml')
        foldername = sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/', type, par.lat_interp);
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s', type);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s', type);
        var = 't';
        file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/echam_ml/ATM_*.ymonmean.nc'));
        fullpath=sprintf('%s/%s', file.folder, file.name);
        ta_orig = double(ncread(fullpath, var));
    elseif strcmp(type, 'echam_pl')
        foldername = sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/', type, par.lat_interp);
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s', type);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s', type);
        var = 't';
        file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/echam_pl/ATM_*.ymonmean.nc'));
        fullpath=sprintf('%s/%s', file.folder, file.name);
        ta_orig = double(ncread(fullpath, var));
    end
    load(sprintf('%s/grid.mat', prefix)); % read grid data
    load(sprintf('%s/srfc.mat', prefix)); % load surface data
    load(sprintf('%s/%s/masks.mat', prefix_proc, par.lat_interp)); % load land and ocean masks

    if strcmp(par.lat_interp, 'std')
        lat = par.lat_std;
    else
        lat = grid.dim3.lat;
    end


    for i = {'lin', 'cub', 'spl', 'mak'}; itp=i{1};
        clear tasi_orig pasi_orig

        load(sprintf('%s/ta_si.mat', prefix)); tasi_orig = ta_si.(itp); clear ta_si; % read temp in si coordinates
        load(sprintf('%s/pa_si.mat', prefix)); pasi_orig = pa_si; clear pa_si; % read temp in si coordinates

        % interpolate ta to standard lat grid
        tasi_orig = permute(tasi_orig, [2 1 3 4]);
        tasi_orig = interp1(grid.dim3.lat, tasi_orig, lat);
        tasi_orig = permute(tasi_orig, [2 1 3 4]);

        pasi_orig = permute(pasi_orig, [2 1 3 4]);
        pasi_orig = interp1(grid.dim3.lat, pasi_orig, lat);
        pasi_orig = permute(pasi_orig, [2 1 3 4]);

        % create surface mask
        if strcmp(type, 'era5') | strcmp(type, 'erai')
            ps = permute(srfc.sp, [2 1 3]); % bring lat to front to interpolate
            ps = interp1(grid.dim2.lat, ps, lat, 'linear', 'extrap'); % interpolate to standard grid
            ps = permute(ps, [2 1 3]); % reorder to original
            ps_vert = repmat(ps, [1 1 1 size(ta_orig, 3)]); % dims (lon x lat x time x plev)
            ps_vert = permute(ps_vert, [1 2 4 3]); % dims (lon x lat x plev x time)
            pa = double(permute(repmat(grid.dim3.plev, [1 size(ps)]), [2 3 1 4]));
            ts = permute(srfc.t2m, [2 1 3]); % bring lat to front to interpolate
            ts = interp1(grid.dim2.lat, ts, lat); % interpolate to standard grid
            ts = permute(ts, [2 1 3]); % reorder to original
        elseif strcmp(type, 'gcm')
            ps = permute(srfc.ps, [2 1 3]); % bring lat to front to interpolate
            ps = interp1(grid.dim2.lat, ps, lat, 'linear', 'extrap'); % interpolate to standard grid
            ps = permute(ps, [2 1 3]); % reorder to original
            ps_vert = repmat(ps, [1 1 1 size(ta_orig, 3)]); % dims (lon x lat x time x plev)
            ps_vert = permute(ps_vert, [1 2 4 3]); % dims (lon x lat x plev x time)
            pa = permute(repmat(grid.dim3.plev, [1 size(ps)]), [2 3 1 4]);
            ts = permute(srfc.tas, [2 1 3]); % bring lat to front to interpolate
            ts = interp1(grid.dim2.lat, ts, lat); % interpolate to standard grid
            ts = permute(ts, [2 1 3]); % reorder to original
        elseif strcmp(type, 'echam') | strcmp(type, 'echam_ml') | strcmp(type, 'echam_pl')
            ps = permute(srfc.aps, [2 1 3]); % bring lat to front to interpolate
            ps = interp1(grid.dim2.lat, ps, lat, 'linear', 'extrap'); % interpolate to standard grid
            ps = permute(ps, [2 1 3]); % reorder to original
            ps_vert = repmat(ps, [1 1 1 size(ta_orig, 3)]); % dims (lon x lat x time x plev)
            ps_vert = permute(ps_vert, [1 2 4 3]); % dims (lon x lat x plev x time)
            pa = permute(repmat(grid.dim3.plev, [1 size(ps)]), [2 3 1 4]);
            ts = permute(srfc.temp2, [2 1 3]); % bring lat to front to interpolate
            ts = interp1(grid.dim2.lat, ts, lat); % interpolate to standard grid
            ts = permute(ts, [2 1 3]); % reorder to original
        end

        if strcmp(type, 'echam_ml')
            sm = ones(size(ta_orig));
        else
            sm = nan(size(ta_orig));
            sm(pa < 0.9961*ps_vert) = 1;
        end

        tasi_sm.lo = tasi_orig; % surface is already masked in standard sigma coordinates

        pasi_sm.lo = pasi_orig; % surface is already masked in spndard sigma coordinates

        % Land/ocean filter 3D variables
        mask.land_vert = repmat(mask.land, [1 1 1 size(ta_orig, 3)]); % expand land mask to vertical dim
        mask.land_vert = permute(mask.land, [1 2 4 3]); % place vertical dim where it belongs
        mask.ocean_vert = repmat(mask.ocean, [1 1 1 size(ta_orig, 3)]); % expand ocean mask to vertical dim
        mask.ocean_vert = permute(mask.ocean, [1 2 4 3]); % place vertical dim where it belongs

        tasi_sm.l = tasi_sm.lo.*mask.ocean_vert; % filter tasi with surface mask
        tasi_sm.o = tasi_sm.lo.*mask.land_vert; % filter tasi with surface mask
        tasi_sm.lo = permute(tasi_sm.lo, [1 2 4 3]); % bring plev to last dimension
        tasi_sm.l = permute(tasi_sm.l, [1 2 4 3]); % bring plev to last dimension
        tasi_sm.o = permute(tasi_sm.o, [1 2 4 3]); % bring plev to last dimension

        pasi_sm.l = pasi_sm.lo.*mask.ocean_vert; % filter pasi with surface mask
        pasi_sm.o = pasi_sm.lo.*mask.land_vert; % filter pasi with surface mask
        pasi_sm.lo = permute(pasi_sm.lo, [1 2 4 3]); % bring plev to last dimension
        pasi_sm.l = permute(pasi_sm.l, [1 2 4 3]); % bring plev to last dimension
        pasi_sm.o = permute(pasi_sm.o, [1 2 4 3]); % bring plev to last dimension

        mask_t.land = nanmean(mask.land, 3);
        mask_t.ocean = nanmean(mask.ocean, 3);

        for l = {'lo', 'l', 'o'}; land = l{1}; % over land, over ocean, or both
            tasi.(land)= squeeze(mean(tasi_sm.(land), 1)); % zonal average
            pasi.(land)= squeeze(mean(pasi_sm.(land), 1)); % zonal average
        end

        % save filtered data
        printname = [foldername 'ta_mean_mon_lat_' itp];
        if ~exist(foldername, 'dir')
            mkdir(foldername)
        end
        if par.do_surf; save(printname, 'tasi', 'lat');
        else save(printname, 'tasi', 'pasi', 'lat', '-v7.3'); end

    end
end % compute mon x lat temperature field
function proc_temp_mon_lat_old(type, par)
    if strcmp(type, 'era5') | strcmp(type, 'erai')
        foldername = sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/', type, par.lat_interp);
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s', type);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s', type);
        file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/%s/temp/%s_temp_%s.ymonmean.nc', type, type, par.(type).yr_span));
        fullpath=sprintf('%s/%s', file.folder, file.name);
        ta_orig = ncread(fullpath, 't');
    elseif strcmp(type, 'gcm')
        foldername = sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s/%s/', type, par.model, par.gcm.clim, par.lat_interp);
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s/%s', type, par.model, par.gcm.clim);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s', type, par.model, par.gcm.clim);
        var = 'ta';
        file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/gcm/%s/%s_Amon_%s_%s_r1i1p1_*.ymonmean.nc', par.model, var, par.model, par.gcm.clim));
        fullpath=sprintf('%s/%s', file.folder, file.name);
        ta_orig = ncread(fullpath, var);
    elseif strcmp(type, 'echam')
        foldername = sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s/', type, par.echam.clim, par.lat_interp);
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s', type, par.echam.clim);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s', type, par.echam.clim);
        var = 't';
        if contains(par.echam.clim, 'rp000')
            file=dir(sprintf('/project2/tas1/ockham/data11/tas/echam-aiv_rcc_6.1.00p1/%s/ATM_%s_0020_39.nc', par.echam.clim, par.echam.clim));
        else
            file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/echam/ATM_rjg_%s_*.ymonmean.nc', par.echam.clim));
        end
        fullpath=sprintf('%s/%s', file.folder, file.name);
        ta_orig = double(ncread(fullpath, var));
    elseif strcmp(type, 'echam_ml')
        foldername = sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/', type, par.lat_interp);
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s', type);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s', type);
        var = 't';
        file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/echam_ml/ATM_*.ymonmean.nc'));
        fullpath=sprintf('%s/%s', file.folder, file.name);
        ta_orig = double(ncread(fullpath, var));
    elseif strcmp(type, 'echam_pl')
        foldername = sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/', type, par.lat_interp);
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s', type);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s', type);
        var = 't';
        file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/echam_pl/ATM_*.ymonmean.nc'));
        fullpath=sprintf('%s/%s', file.folder, file.name);
        ta_orig = double(ncread(fullpath, var));
    end
    load(sprintf('%s/grid.mat', prefix)); % read grid data
    load(sprintf('%s/tempz.mat', prefix)); taz_orig = tempz; clear tempz; % read temp in z coordinates
    load(sprintf('%s/pz.mat', prefix)); paz_orig = pz; clear pz; % read pa in z coordinates
    load(sprintf('%s/ta_si.mat', prefix)); tasi_orig = ta_si.spl; clear ta_si; % read temp in si coordinates
    load(sprintf('%s/pa_si.mat', prefix)); pasi_orig = pa_si; clear pa_si; % read temp in si coordinates
    load(sprintf('%s/srfc.mat', prefix)); % load surface data
    load(sprintf('%s/%s/masks.mat', prefix_proc, par.lat_interp)); % load land and ocean masks

    if strcmp(par.lat_interp, 'std')
        lat = par.lat_std;
    else
        lat = grid.dim3.lat;
    end

    % interpolate ta to standard lat grid
    ta_orig = permute(ta_orig, [2 1 3 4]);
    ta_orig = interp1(grid.dim3.lat, ta_orig, lat);
    ta_orig = permute(ta_orig, [2 1 3 4]);
    taz_orig = permute(taz_orig, [2 1 3 4]);
    taz_orig = interp1(grid.dim3.lat, taz_orig, lat);
    taz_orig = permute(taz_orig, [2 1 3 4]);
    tasi_orig = permute(tasi_orig, [2 1 3 4]);
    tasi_orig = interp1(grid.dim3.lat, tasi_orig, lat);
    tasi_orig = permute(tasi_orig, [2 1 3 4]);

    paz_orig = permute(paz_orig, [2 1 3 4]);
    paz_orig = interp1(grid.dim3.lat, paz_orig, lat);
    paz_orig = permute(paz_orig, [2 1 3 4]);
    pasi_orig = permute(pasi_orig, [2 1 3 4]);
    pasi_orig = interp1(grid.dim3.lat, pasi_orig, lat);
    pasi_orig = permute(pasi_orig, [2 1 3 4]);

    % create surface mask
    if strcmp(type, 'era5') | strcmp(type, 'erai')
        ps = permute(srfc.sp, [2 1 3]); % bring lat to front to interpolate
        ps = interp1(grid.dim2.lat, ps, lat, 'linear', 'extrap'); % interpolate to standard grid
        ps = permute(ps, [2 1 3]); % reorder to original
        ps_vert = repmat(ps, [1 1 1 size(ta_orig, 3)]); % dims (lon x lat x time x plev)
        ps_vert = permute(ps_vert, [1 2 4 3]); % dims (lon x lat x plev x time)
        pa = double(permute(repmat(grid.dim3.plev, [1 size(ps)]), [2 3 1 4]));
        zs = permute(srfc.zs, [2 1 3]); % bring lat to front to interpolate
        zs = interp1(grid.dim2.lat, zs, lat); % interpolate to standard grid
        zs = permute(zs, [2 1 3]); % reorder to original
        zs_vert = repmat(zs, [1 1 1 size(taz_orig, 3)]); % dims (lon x lat x time x plev)
        zs_vert = permute(zs_vert, [1 2 4 3]); % dims (lon x lat x plev x time)
        za = permute(repmat(grid.dim3.z, [1 size(zs)]), [2 3 1 4]);
        ts = permute(srfc.t2m, [2 1 3]); % bring lat to front to interpolate
        ts = interp1(grid.dim2.lat, ts, lat); % interpolate to standard grid
        ts = permute(ts, [2 1 3]); % reorder to original
        ts_vert = repmat(ts, [1 1 1 size(taz_orig, 3)]); % dims (lon x lat x time x plev)
        ts_vert = permute(ts_vert, [1 2 4 3]); % dims (lon x lat x plev x time)
    elseif strcmp(type, 'gcm')
        ps = permute(srfc.ps, [2 1 3]); % bring lat to front to interpolate
        ps = interp1(grid.dim2.lat, ps, lat, 'linear', 'extrap'); % interpolate to standard grid
        ps = permute(ps, [2 1 3]); % reorder to original
        ps_vert = repmat(ps, [1 1 1 size(ta_orig, 3)]); % dims (lon x lat x time x plev)
        ps_vert = permute(ps_vert, [1 2 4 3]); % dims (lon x lat x plev x time)
        pa = permute(repmat(grid.dim3.plev, [1 size(ps)]), [2 3 1 4]);
        zs = permute(srfc.zs, [2 1 3]); % bring lat to front to interpolate
        zs = interp1(grid.dim2.lat, zs, lat); % interpolate to standard grid
        zs = permute(zs, [2 1 3]); % reorder to original
        zs_vert = repmat(zs, [1 1 1 size(taz_orig, 3)]); % dims (lon x lat x time x plev)
        zs_vert = permute(zs_vert, [1 2 4 3]); % dims (lon x lat x plev x time)
        za = permute(repmat(grid.dim3.z, [1 size(zs)]), [2 3 1 4]);
        ts = permute(srfc.tas, [2 1 3]); % bring lat to front to interpolate
        ts = interp1(grid.dim2.lat, ts, lat); % interpolate to standard grid
        ts = permute(ts, [2 1 3]); % reorder to original
        ts_vert = repmat(ts, [1 1 1 size(taz_orig, 3)]); % dims (lon x lat x time x plev)
        ts_vert = permute(ts_vert, [1 2 4 3]); % dims (lon x lat x plev x time)
    elseif strcmp(type, 'echam') | strcmp(type, 'echam_ml') | strcmp(type, 'echam_pl')
        ps = permute(srfc.aps, [2 1 3]); % bring lat to front to interpolate
        ps = interp1(grid.dim2.lat, ps, lat, 'linear', 'extrap'); % interpolate to standard grid
        ps = permute(ps, [2 1 3]); % reorder to original
        ps_vert = repmat(ps, [1 1 1 size(ta_orig, 3)]); % dims (lon x lat x time x plev)
        ps_vert = permute(ps_vert, [1 2 4 3]); % dims (lon x lat x plev x time)
        pa = permute(repmat(grid.dim3.plev, [1 size(ps)]), [2 3 1 4]);
        zs = permute(srfc.zs, [2 1 3]); % bring lat to front to interpolate
        zs = interp1(grid.dim2.lat, zs, lat); % interpolate to standard grid
        zs = permute(zs, [2 1 3]); % reorder to original
        zs_vert = repmat(zs, [1 1 1 size(taz_orig, 3)]); % dims (lon x lat x time x plev)
        zs_vert = permute(zs_vert, [1 2 4 3]); % dims (lon x lat x plev x time)
        za = permute(repmat(grid.dim3.z, [1 size(zs)]), [2 3 1 4]);
        ts = permute(srfc.temp2, [2 1 3]); % bring lat to front to interpolate
        ts = interp1(grid.dim2.lat, ts, lat); % interpolate to standard grid
        ts = permute(ts, [2 1 3]); % reorder to original
        ts_vert = repmat(ts, [1 1 1 size(taz_orig, 3)]); % dims (lon x lat x time x plev)
        ts_vert = permute(ts_vert, [1 2 4 3]); % dims (lon x lat x plev x time)
    end

    if strcmp(type, 'echam_ml')
        sm = ones(size(ta_orig));
        smz = ones(size(taz_orig));
    else
        sm = nan(size(ta_orig));
        sm(pa < ps_vert) = 1;
        smz = nan(size(taz_orig));
        smz(za > zs_vert) = 1;
    end

    ta_sm.lo = ta_orig.*sm; % filter ta with surface mask
    taz_sm.lo = taz_orig.*smz; % filter taz with surface mask
    tasi_sm.lo = tasi_orig; % surface is already masked in standard sigma coordinates

    paz_sm.lo = paz_orig.*smz; % filter paz with surface mask
    pasi_sm.lo = pasi_orig; % surface is already masked in spndard sigma coordinates

    if par.do_surf
        % add tsurf data and interpolate to higher resolution vertical grid
        [pa_plus ta_plus] = deal(nan([size(pa,1), size(pa,2) size(pa,3)+1 size(pa,4)])); % create empty grid with one extra vertical level
        pa_plus(:,:,1:end-1,:) = pa; % populate with standard pressure grid
        ta_plus(:,:,1:end-1,:) = ta_sm.lo; % populate with standard atmospheric temperature
        pa_plus(:,:,end,:) = ps_vert(:,:,1,:); % add surface pressure data into standard pressure grid
        ta_plus(:,:,end,:) = ts_vert(:,:,1,:); % add surface temperature data
        pa_plus = permute(pa_plus, [3 1 2 4]); % bring plev dimension to front
        ta_plus = permute(ta_plus, [3 1 2 4]); % bring plev dimension to front
        [pa_plus sort_index] = sort(pa_plus, 1, 'descend'); % sort added surface pressure such that pressure decreases monotonically
        tai_sm.lo = nan(length(par.pa), size(pa, 1), size(pa, 2), size(pa, 4));
        pb = CmdLineProgressBar("Sorting and interpolating temperature to new standard grid...");
        for ilon=1:size(pa_plus,2)
            pb.print(ilon, size(pa_plus,2));
            for ilat=1:size(pa_plus,3)
                for time=1:size(pa_plus,4)
                    ta_plus(:,ilon,ilat,time) = ta_plus(sort_index(:,ilon,ilat,time),ilon,ilat,time); % sort temperature (has to be in loop because sort_index works for vector calls only)
                    tai_sm.lo(:,ilon,ilat,time) = interp1(pa_plus(:,ilon,ilat,time), ta_plus(:,ilon,ilat,time), par.pa, 'linear'); % interpolate to higher resolution vertical grid
                end
            end
        end
        clear pa_plus ta_plus; % clear unneeded variables
        tai_sm.lo = permute(tai_sm.lo, [2 3 1 4]); % bring plev back to third dimension
    end

    % Land/ocean filter 2D variables
    tsurf_sm.lo = ts;
    psurf_sm.lo = ps;
    zsurf_sm.lo = zs;
    tsurf_sm.l = ts.*mask.ocean; %filter out ocean
    psurf_sm.l = ps.*mask.ocean; %filter out ocean
    zsurf_sm.l = zs.*mask.ocean; %filter out ocean
    tsurf_sm.o = ts.*mask.land; %filter out land
    psurf_sm.o = ps.*mask.land; %filter out land
    zsurf_sm.o = zs.*mask.land; %filter out land

    % Land/ocean filter 3D variables
    mask.land_vert = repmat(mask.land, [1 1 1 size(ta_orig, 3)]); % expand land mask to vertical dim
    mask.land_vert = permute(mask.land, [1 2 4 3]); % place vertical dim where it belongs
    mask.ocean_vert = repmat(mask.ocean, [1 1 1 size(ta_orig, 3)]); % expand ocean mask to vertical dim
    mask.ocean_vert = permute(mask.ocean, [1 2 4 3]); % place vertical dim where it belongs
    mask.land_verti = repmat(mask.land, [1 1 1 length(par.pa)]); % expand land mask to vertiical dim
    mask.land_verti = permute(mask.land, [1 2 4 3]); % place vertiical dim where it belongs
    mask.ocean_verti = repmat(mask.ocean, [1 1 1 length(par.pa)]); % expand ocean mask to vertiical dim
    mask.ocean_verti = permute(mask.ocean, [1 2 4 3]); % place vertiical dim where it belongs

    ta_sm.l = ta_sm.lo.*mask.ocean_vert; % filter ta with surface mask
    ta_sm.o = ta_sm.lo.*mask.land_vert; % filter ta with surface mask
    ta_sm.lo = permute(ta_sm.lo, [1 2 4 3]); % bring plev to last dimension
    ta_sm.l = permute(ta_sm.l, [1 2 4 3]); % bring plev to last dimension
    ta_sm.o = permute(ta_sm.o, [1 2 4 3]); % bring plev to last dimension
    taz_sm.l = taz_sm.lo.*mask.ocean_vert; % filter taz with surface mask
    taz_sm.o = taz_sm.lo.*mask.land_vert; % filter taz with surface mask
    taz_sm.lo = permute(taz_sm.lo, [1 2 4 3]); % bring plev to last dimension
    taz_sm.l = permute(taz_sm.l, [1 2 4 3]); % bring plev to last dimension
    taz_sm.o = permute(taz_sm.o, [1 2 4 3]); % bring plev to last dimension
    tasi_sm.l = tasi_sm.lo.*mask.ocean_vert; % filter tasi with surface mask
    tasi_sm.o = tasi_sm.lo.*mask.land_vert; % filter tasi with surface mask
    tasi_sm.lo = permute(tasi_sm.lo, [1 2 4 3]); % bring plev to last dimension
    tasi_sm.l = permute(tasi_sm.l, [1 2 4 3]); % bring plev to last dimension
    tasi_sm.o = permute(tasi_sm.o, [1 2 4 3]); % bring plev to last dimension

    paz_sm.l = paz_sm.lo.*mask.ocean_vert; % filter paz with surface mask
    paz_sm.o = paz_sm.lo.*mask.land_vert; % filter paz with surface mask
    paz_sm.lo = permute(paz_sm.lo, [1 2 4 3]); % bring plev to last dimension
    paz_sm.l = permute(paz_sm.l, [1 2 4 3]); % bring plev to last dimension
    paz_sm.o = permute(paz_sm.o, [1 2 4 3]); % bring plev to last dimension
    pasi_sm.l = pasi_sm.lo.*mask.ocean_vert; % filter pasi with surface mask
    pasi_sm.o = pasi_sm.lo.*mask.land_vert; % filter pasi with surface mask
    pasi_sm.lo = permute(pasi_sm.lo, [1 2 4 3]); % bring plev to last dimension
    pasi_sm.l = permute(pasi_sm.l, [1 2 4 3]); % bring plev to last dimension
    pasi_sm.o = permute(pasi_sm.o, [1 2 4 3]); % bring plev to last dimension
    if par.do_surf
        tai_sm.l = tai_sm.lo.*mask.ocean_verti; % filter tai with surface mask
        tai_sm.o = tai_sm.lo.*mask.land_verti; % filter tai with surface mask
        tai_sm.lo = permute(tai_sm.lo, [1 2 4 3]); % bring plev to last dimension
        tai_sm.l = permute(tai_sm.l, [1 2 4 3]); % bring plev to last dimension
        tai_sm.o = permute(tai_sm.o, [1 2 4 3]); % bring plev to last dimension
    end

    mask_t.land = nanmean(mask.land, 3);
    mask_t.ocean = nanmean(mask.ocean, 3);

    for l = {'lo', 'l', 'o'}; land = l{1}; % over land, over ocean, or both
        ta.(land)= squeeze(mean(ta_sm.(land), 1)); % zonal average
        taz.(land)= squeeze(mean(taz_sm.(land), 1)); % zonal average
        tasi.(land)= squeeze(mean(tasi_sm.(land), 1)); % zonal average

        paz.(land)= squeeze(mean(paz_sm.(land), 1)); % zonal average
        pasi.(land)= squeeze(mean(pasi_sm.(land), 1)); % zonal average

        if par.do_surf; tai.(land)= squeeze(mean(tai_sm.(land), 1)); end % zonal average
        tsurf.(land)= squeeze(mean(tsurf_sm.(land), 1)); % zonal average
        psurf.(land)= squeeze(mean(psurf_sm.(land), 1)); % zonal average
        zsurf.(land)= squeeze(mean(zsurf_sm.(land), 1)); % zonal average
    end

    for l = {'lo', 'l', 'o'}; land = l{1}; % over land, over ocean, or both
        % take time averages
        for t = {'ann', 'djf', 'jja', 'mam', 'son'}; time = t{1};
            if strcmp(time, 'ann')
                ta_t.(land).(time) = squeeze(mean(ta_sm.(land), 3));
                taz_t.(land).(time) = squeeze(mean(taz_sm.(land), 3));
                tasi_t.(land).(time) = squeeze(mean(tasi_sm.(land), 3));

                paz_t.(land).(time) = squeeze(mean(paz_sm.(land), 3));
                pasi_t.(land).(time) = squeeze(mean(pasi_sm.(land), 3));

                if par.do_surf; tai_t.(land).(time) = squeeze(mean(tai_sm.(land), 3)); end
                tsurf_t.(land).(time) = squeeze(mean(tsurf_sm.(land), 3));
                psurf_t.(land).(time) = squeeze(mean(psurf_sm.(land), 3));
                zsurf_t.(land).(time) = squeeze(mean(zsurf_sm.(land), 3));
            elseif strcmp(time, 'djf')
                ta_shift.(land) = circshift(ta_sm.(land), 1, 3);
                taz_shift.(land) = circshift(taz_sm.(land), 1, 3);
                tasi_shift.(land) = circshift(tasi_sm.(land), 1, 3);

                paz_shift.(land) = circshift(paz_sm.(land), 1, 3);
                pasi_shift.(land) = circshift(pasi_sm.(land), 1, 3);

                if par.do_surf; tai_shift.(land) = circshift(tai_sm.(land), 1, 3); end
                tsurf_shift.(land) = circshift(tsurf_sm.(land), 1, 3);
                psurf_shift.(land) = circshift(psurf_sm.(land), 1, 3);
                zsurf_shift.(land) = circshift(zsurf_sm.(land), 1, 3);
                ta_t.(land).(time) = squeeze(mean(ta_shift.(land)(:,:,1:3,:), 3));
                taz_t.(land).(time) = squeeze(mean(taz_shift.(land)(:,:,1:3,:), 3));
                tasi_t.(land).(time) = squeeze(mean(tasi_shift.(land)(:,:,1:3,:), 3));

                paz_t.(land).(time) = squeeze(mean(paz_shift.(land)(:,:,1:3,:), 3));
                pasi_t.(land).(time) = squeeze(mean(pasi_shift.(land)(:,:,1:3,:), 3));

                if par.do_surf; tai_t.(land).(time) = squeeze(mean(tai_shift.(land)(:,:,1:3,:), 3)); end
                tsurf_t.(land).(time) = squeeze(mean(tsurf_shift.(land)(:,:,1:3,:), 3));
                psurf_t.(land).(time) = squeeze(mean(psurf_shift.(land)(:,:,1:3,:), 3));
                zsurf_t.(land).(time) = squeeze(mean(zsurf_shift.(land)(:,:,1:3,:), 3));
            elseif strcmp(time, 'jja')
                ta_t.(land).(time) = squeeze(mean(ta_sm.(land)(:,:,6:8,:), 3));
                taz_t.(land).(time) = squeeze(mean(taz_sm.(land)(:,:,6:8,:), 3));
                tasi_t.(land).(time) = squeeze(mean(tasi_sm.(land)(:,:,6:8,:), 3));

                paz_t.(land).(time) = squeeze(mean(paz_sm.(land)(:,:,6:8,:), 3));
                pasi_t.(land).(time) = squeeze(mean(pasi_sm.(land)(:,:,6:8,:), 3));

                if par.do_surf; tai_t.(land).(time) = squeeze(mean(tai_sm.(land)(:,:,6:8,:), 3)); end
                tsurf_t.(land).(time) = squeeze(mean(tsurf_sm.(land)(:,:,6:8,:), 3));
                psurf_t.(land).(time) = squeeze(mean(psurf_sm.(land)(:,:,6:8,:), 3));
                zsurf_t.(land).(time) = squeeze(mean(zsurf_sm.(land)(:,:,6:8,:), 3));
            elseif strcmp(time, 'mam')
                ta_t.(land).(time) = squeeze(mean(ta_sm.(land)(:,:,3:5,:), 3));
                taz_t.(land).(time) = squeeze(mean(taz_sm.(land)(:,:,3:5,:), 3));
                tasi_t.(land).(time) = squeeze(mean(tasi_sm.(land)(:,:,3:5,:), 3));

                paz_t.(land).(time) = squeeze(mean(paz_sm.(land)(:,:,3:5,:), 3));
                pasi_t.(land).(time) = squeeze(mean(pasi_sm.(land)(:,:,3:5,:), 3));

                if par.do_surf; tai_t.(land).(time) = squeeze(mean(tai_sm.(land)(:,:,3:5,:), 3)); end
                tsurf_t.(land).(time) = squeeze(mean(tsurf_sm.(land)(:,:,3:5,:), 3));
                psurf_t.(land).(time) = squeeze(mean(psurf_sm.(land)(:,:,3:5,:), 3));
                zsurf_t.(land).(time) = squeeze(mean(zsurf_sm.(land)(:,:,3:5,:), 3));
            elseif strcmp(time, 'son')
                ta_t.(land).(time) = squeeze(mean(ta_sm.(land)(:,:,9:11,:), 3));
                taz_t.(land).(time) = squeeze(mean(taz_sm.(land)(:,:,9:11,:), 3));
                tasi_t.(land).(time) = squeeze(mean(tasi_sm.(land)(:,:,9:11,:), 3));

                paz_t.(land).(time) = squeeze(mean(paz_sm.(land)(:,:,9:11,:), 3));
                pasi_t.(land).(time) = squeeze(mean(pasi_sm.(land)(:,:,9:11,:), 3));

                if par.do_surf; tai_t.(land).(time) = squeeze(mean(tai_sm.(land)(:,:,9:11,:), 3)); end
                tsurf_t.(land).(time) = squeeze(mean(tsurf_sm.(land)(:,:,9:11,:), 3));
                psurf_t.(land).(time) = squeeze(mean(psurf_sm.(land)(:,:,9:11,:), 3));
                zsurf_t.(land).(time) = squeeze(mean(zsurf_sm.(land)(:,:,9:11,:), 3));
            end
        end
    end

    % save filtered data
    printname = [foldername 'ta_mon_lat'];
    if ~exist(foldername, 'dir')
        mkdir(foldername)
    end
    if par.do_surf; save(printname, 'ta', 'taz', 'tasi', 'tai', 'lat', 'tsurf', 'psurf', 'zsurf');
    else save(printname, 'ta', 'taz', 'tasi', 'paz', 'pasi', 'lat', 'tsurf', 'psurf', 'zsurf', '-v7.3'); end

    printname = [foldername 'ta_lon_lat'];
    if ~exist(foldername, 'dir')
        mkdir(foldername)
    end
    if par.do_surf; save(printname, 'ta_t', 'taz_t', 'tasi_t', 'tai_t', 'lat', 'tsurf_t', 'psurf_t', 'zsurf_t');
    else save(printname, 'ta_t', 'taz_t', 'tasi_t', 'paz_t', 'pasi_t', 'lat', 'tsurf_t', 'psurf_t', 'zsurf_t', '-v7.3'); end
end % compute mon x lat temperature field

function make_dtdzzsi_old(type, par) % compute model lapse rate in lat x plev x mon (add surface data before converting to sigma)
% calculate moist adiabats
    if strcmp(type, 'era5') | strcmp(type, 'erai')
        foldername = sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/', type, par.lat_interp);
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s', type);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s', type);
        file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/%s/temp/%s_temp_%s.ymonmean.nc', type, type, par.(type).yr_span));
        fullpath=sprintf('%s/%s', file.folder, file.name);
        temp = ncread(fullpath, 't');
        file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/%s/zg/%s_zg_%s.ymonmean.nc', type, type, par.(type).yr_span));
        fullpath=sprintf('%s/%s', file.folder, file.name);
        zg = ncread(fullpath, 'z'); zg = zg/par.g; % convert from geopotential to height in m
    elseif strcmp(type, 'gcm')
        foldername = sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s/%s/', type, par.model, par.gcm.clim, par.lat_interp);
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s/%s', type, par.model, par.gcm.clim);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s', type, par.model, par.gcm.clim);
        var = 'ta';
        file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/gcm/%s/%s_Amon_%s_%s_r1i1p1_*.ymonmean.nc', par.model, var, par.model, par.gcm.clim));
        fullpath=sprintf('%s/%s', file.folder, file.name);
        temp = ncread(fullpath, var);
        var = 'zg';
        file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/gcm/%s/%s_Amon_%s_%s_r1i1p1_*.ymonmean.nc', par.model, var, par.model, par.gcm.clim));
        fullpath=sprintf('%s/%s', file.folder, file.name);
        zg = ncread(fullpath, var);
    elseif strcmp(type, 'echam')
        foldername = sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s/', type, par.echam.clim, par.lat_interp);
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s', type, par.echam.clim);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s', type, par.echam.clim);
        var = 't';
        file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/echam/ATM_rjg_%s_*.ymonmean.nc', par.echam.clim));
        fullpath=sprintf('%s/%s', file.folder, file.name);
        temp = double(ncread(fullpath, var));
        var = 'geopoth';
        file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/echam/ATM_rjg_%s_*.ymonmean.nc', par.echam.clim));
        fullpath=sprintf('%s/%s', file.folder, file.name);
        zg = double(ncread(fullpath, var));
    end
    load(sprintf('%s/grid.mat', prefix)); % read grid data
    load(sprintf('%s/srfc.mat', prefix)); % read surface variable data
    % load(sprintf('%s/orog.mat', prefix)); % read surface orography

    % create surface mask
    if strcmp(type, 'era5') | strcmp(type, 'erai')
        ps_vert = repmat(srfc.sp, [1 1 1 size(temp, 3)]); % dims (lon x lat x time x plev)
        ps_vert = permute(ps_vert, [1 2 4 3]); % dims (lon x lat x plev x time)
        pa = double(permute(repmat(grid.dim3.plev, [1 size(srfc.sp)]), [2 3 1 4]));
    elseif strcmp(type, 'gcm')
        ps_vert = repmat(srfc.ps, [1 1 1 size(temp, 3)]); % dims (lon x lat x time x plev)
        ps_vert = permute(ps_vert, [1 2 4 3]); % dims (lon x lat x plev x time)
        pa = permute(repmat(grid.dim3.plev, [1 size(srfc.ps)]), [2 3 1 4]);
    elseif strcmp(type, 'echam')
        ps_vert = repmat(srfc.aps, [1 1 1 size(temp, 3)]); % dims (lon x lat x time x plev)
        ps_vert = permute(ps_vert, [1 2 4 3]); % dims (lon x lat x plev x time)
        pa = permute(repmat(grid.dim3.plev, [1 size(srfc.aps)]), [2 3 1 4]);
    end
    surface_mask = nan(size(pa));
    surface_mask(pa < ps_vert) = 1;

    ta_sm = temp .* surface_mask;
    zg_sm = zg .* surface_mask;

    % add tsurf data and interpolate to higher resolution vertical grid
    [pa_plus ta_plus zg_plus] = deal(nan([size(pa,1), size(pa,2) size(pa,3)+1 size(pa,4)])); % create empty grid with one extra vertical level
    pa_plus(:,:,1:end-1,:) = pa; % populate with standard pressure grid
    ta_plus(:,:,1:end-1,:) = ta_sm; % populate with standard atmospheric temperature
    zg_plus(:,:,1:end-1,:) = zg_sm; % populate with szgndard atmospheric temperature
    pa_plus(:,:,end,:) = ps_vert(:,:,1,:); % add surface pressure data into standard pressure grid
    if any(strcmp(type, {'era5', 'erai'})); ta_plus(:,:,end,:) = srfc.t2m(:,:,:); % add surface temperature data
    elseif strcmp(type, 'gcm'); ta_plus(:,:,end,:) = srfc.tas(:,:,:); % add surface temperature data
    elseif strcmp(type, 'echam'); ta_plus(:,:,end,:) = srfc.temp2(:,:,:); end % add surface temperature data
    % zg_plus(:,:,end,:) = srfc.zs(:,:,:); % add surface height data
    % zg_plus(:,:,end,:) = repmat(orog, [1 1 12]); % add surface height data
    zg_plus(:,:,end,:) = nan(size(srfc.zs)); % add surface height data
    ps_vert = permute(ps_vert, [3 1 2 4]); % bring plev dimension to front
    pa_plus = permute(pa_plus, [3 1 2 4]); % bring plev dimension to front
    ta_plus = permute(ta_plus, [3 1 2 4]); % bring plev dimension to front
    zg_plus = permute(zg_plus, [3 1 2 4]); % bring plev dimension to front
    [pa_plus sort_index] = sort(pa_plus, 1, 'descend'); % sort added surface pressure such that pressure decreases monotonically
    pb = CmdLineProgressBar("Sorting temperature with surface data added...");
    for lo=1:size(pa_plus,2)
        pb.print(lo, size(pa_plus,2));
        for la=1:size(pa_plus,3)
            for mo=1:size(pa_plus,4)
                ta_plus(:,lo,la,mo) = ta_plus(sort_index(:,lo,la,mo),lo,la,mo); % sort temperature (has to be in loop because sort_index works for vector calls only)
                zg_plus(:,lo,la,mo) = zg_plus(sort_index(:,lo,la,mo),lo,la,mo); % sort temperature (has to be in loop because sort_index works for vector calls only)
            end
        end
    end

    % convert to sigma coord.
    pb = CmdLineProgressBar("Sorting and interpolating dtdzz to new standard grid...");
    for lo=1:size(pa_plus,2)
        pb.print(lo, size(pa_plus,2));
        for la=1:size(pa_plus,3)
            for mo=1:size(pa_plus,4)
                ta_si(:,lo,la,mo) = interp1(pa_plus(:,lo,la,mo)./ps_vert(1,lo,la,mo), ta_plus(:,lo,la,mo), grid.dim3.si, 'linear');
                zg_si(:,lo,la,mo) = interp1(pa_plus(:,lo,la,mo)./ps_vert(1,lo,la,mo), zg_plus(:,lo,la,mo), grid.dim3.si, 'linear');
            end
        end
    end
    clear ta_plus zg_plus

    % take zonal average before computing lapse rate
    ta_si_z = squeeze(nanmean(ta_si,2)); clear ta_si;
    zg_si_z = squeeze(nanmean(zg_si,2)); clear zg_si;

    dtdzzsi = -1e3*(ta_si_z(2:end,:,:)-ta_si_z(1:end-1,:,:))./(zg_si_z(2:end,:,:)-zg_si_z(1:end-1,:,:)); % lapse rate in K/km
    dtdzzsi = interp1(1/2*(grid.dim3.si(2:end)+grid.dim3.si(1:end-1)), dtdzzsi, grid.dim3.si);
    dtdzzsi = repmat(dtdzzsi, [1 1 1 length(grid.dim3.lon)]); % repeat in longitude

    dtdzzsi = permute(dtdzzsi, [4 2 1 3]); % bring height to 3rd

    if any(strcmp(type, {'era5', 'erai'})); newdir=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s', type);
    elseif strcmp(type, 'gcm'); newdir=sprintf('/project2/tas1/miyawaki/projects/002/data/read/gcm/%s/%s', par.model, par.gcm.clim);
    elseif strcmp(type, 'echam'); newdir=sprintf('/project2/tas1/miyawaki/projects/002/data/read/echam/%s', par.echam.clim); end;
    if ~exist(newdir, 'dir'); mkdir(newdir); end
    filename='dtdzzsi.mat';
    save(sprintf('%s/%s', newdir, filename), 'dtdzzsi', '-v7.3');

end
function make_dtdzzsi(type, par) % compute model lapse rate in lat x plev x mon (add surface data after converting to sigma)
% calculate moist adiabats
    if strcmp(type, 'era5') | strcmp(type, 'erai')
        foldername = sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/', type, par.lat_interp);
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s', type);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s', type);
        file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/%s/temp/%s_temp_%s.ymonmean.nc', type, type, par.(type).yr_span));
        fullpath=sprintf('%s/%s', file.folder, file.name);
        temp = ncread(fullpath, 't');
        file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/%s/zg/%s_zg_%s.ymonmean.nc', type, type, par.(type).yr_span));
        fullpath=sprintf('%s/%s', file.folder, file.name);
        zg = ncread(fullpath, 'z'); zg = zg/par.g; % convert from geopotential to height in m
    elseif strcmp(type, 'gcm')
        foldername = sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s/%s/', type, par.model, par.gcm.clim, par.lat_interp);
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s/%s', type, par.model, par.gcm.clim);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s', type, par.model, par.gcm.clim);
        var = 'ta';
        file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/gcm/%s/%s_Amon_%s_%s_r1i1p1_*.ymonmean.nc', par.model, var, par.model, par.gcm.clim));
        fullpath=sprintf('%s/%s', file.folder, file.name);
        temp = ncread(fullpath, var);
        var = 'zg';
        file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/gcm/%s/%s_Amon_%s_%s_r1i1p1_*.ymonmean.nc', par.model, var, par.model, par.gcm.clim));
        fullpath=sprintf('%s/%s', file.folder, file.name);
        zg = ncread(fullpath, var);
    elseif strcmp(type, 'echam')
        foldername = sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s/', type, par.echam.clim, par.lat_interp);
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s', type, par.echam.clim);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s', type, par.echam.clim);
        var = 't';
        file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/echam/ATM_rjg_%s_*.ymonmean.nc', par.echam.clim));
        fullpath=sprintf('%s/%s', file.folder, file.name);
        temp = double(ncread(fullpath, var));
        var = 'geopoth';
        file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/echam/ATM_rjg_%s_*.ymonmean.nc', par.echam.clim));
        fullpath=sprintf('%s/%s', file.folder, file.name);
        zg = double(ncread(fullpath, var));
    end
    load(sprintf('%s/grid.mat', prefix)); % read grid data
    load(sprintf('%s/srfc.mat', prefix)); % read surface variable data
    % load(sprintf('%s/orog.mat', prefix)); % read surface orography

    % create surface mask
    if strcmp(type, 'era5') | strcmp(type, 'erai')
        ps_vert = repmat(srfc.sp, [1 1 1 size(temp, 3)]); % dims (lon x lat x time x plev)
        ps_vert = permute(ps_vert, [1 2 4 3]); % dims (lon x lat x plev x time)
        pa = double(permute(repmat(grid.dim3.plev, [1 size(srfc.sp)]), [2 3 1 4]));
    elseif strcmp(type, 'gcm')
        ps_vert = repmat(srfc.ps, [1 1 1 size(temp, 3)]); % dims (lon x lat x time x plev)
        ps_vert = permute(ps_vert, [1 2 4 3]); % dims (lon x lat x plev x time)
        pa = permute(repmat(grid.dim3.plev, [1 size(srfc.ps)]), [2 3 1 4]);
    elseif strcmp(type, 'echam')
        ps_vert = repmat(srfc.aps, [1 1 1 size(temp, 3)]); % dims (lon x lat x time x plev)
        ps_vert = permute(ps_vert, [1 2 4 3]); % dims (lon x lat x plev x time)
        pa = permute(repmat(grid.dim3.plev, [1 size(srfc.aps)]), [2 3 1 4]);
    end
    surface_mask = nan(size(pa));
    surface_mask(pa < ps_vert) = 1;

    ta_sm = temp .* surface_mask;
    zg_sm = zg .* surface_mask;

    ps_vert = permute(ps_vert, [3 1 2 4]);
    pa = permute(pa, [3 1 2 4]);
    ta_sm = permute(ta_sm, [3 1 2 4]);
    zg_sm = permute(zg_sm, [3 1 2 4]);
    % convert to sigma coord.
    pb = CmdLineProgressBar("Sorting and interpolating dtdzz to new standard grid...");
    for lo=1:size(pa,2)
        pb.print(lo, size(pa,2));
        for la=1:size(pa,3)
            for mo=1:size(pa,4)
                ta_si(:,lo,la,mo) = interp1(pa(:,lo,la,mo)./ps_vert(1,lo,la,mo), ta_sm(:,lo,la,mo), grid.dim3.si, 'linear');
                zg_si(:,lo,la,mo) = interp1(pa(:,lo,la,mo)./ps_vert(1,lo,la,mo), zg_sm(:,lo,la,mo), grid.dim3.si, 'linear');
            end
        end
    end
    clear ta_sm zg_sm

    % add surface data
    ta_si(1,:,:,:) = srfc.tas;
    zg_si(1,:,:,:) = srfc.zs;

    % take zonal average before computing lapse rate
    ta_si_z = squeeze(nanmean(ta_si,2)); clear ta_si;
    zg_si_z = squeeze(nanmean(zg_si,2)); clear zg_si;

    dtdzzsi = -1e3*(ta_si_z(2:end,:,:)-ta_si_z(1:end-1,:,:))./(zg_si_z(2:end,:,:)-zg_si_z(1:end-1,:,:)); % lapse rate in K/km
    dtdzzsi = interp1(1/2*(grid.dim3.si(2:end)+grid.dim3.si(1:end-1)), dtdzzsi, grid.dim3.si);
    dtdzzsi = repmat(dtdzzsi, [1 1 1 length(grid.dim3.lon)]); % repeat in longitude
    dtdzzsi = permute(dtdzzsi, [4 2 1 3]); % bring height to 3rd

    if any(strcmp(type, {'era5', 'erai'})); newdir=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s', type);
    elseif strcmp(type, 'gcm'); newdir=sprintf('/project2/tas1/miyawaki/projects/002/data/read/gcm/%s/%s', par.model, par.gcm.clim);
    elseif strcmp(type, 'echam'); newdir=sprintf('/project2/tas1/miyawaki/projects/002/data/read/echam/%s', par.echam.clim); end;
    if ~exist(newdir, 'dir'); mkdir(newdir); end
    filename='dtdzzsi.mat';
    save(sprintf('%s/%s', newdir, filename), 'dtdzzsi', '-v7.3');

end

function make_malr(type, par) % compute moist adiabatic lapse rate at every lat, lon, time
% calculate moist adiabats
    if strcmp(type, 'era5') | strcmp(type, 'erai')
        foldername = sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/', type, par.lat_interp);
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s', type);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s', type);
    elseif strcmp(type, 'gcm')
        foldername = sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s/%s/', type, par.model, par.gcm.clim, par.lat_interp);
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s/%s', type, par.model, par.gcm.clim);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s', type, par.model, par.gcm.clim);
    elseif strcmp(type, 'echam')
        foldername = sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s/', type, par.echam.clim, par.lat_interp);
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s', type, par.echam.clim);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s', type, par.echam.clim);
    end
    load(sprintf('%s/grid.mat', prefix)); % read grid data
    load(sprintf('%s/tai.mat', prefix)); clear zgi; % read interpolated temperature

    pa = repmat(par.pa', [1 length(grid.dim3.lon) length(grid.dim3.lat) 12]);
    pa = permute(pa, [2 3 1 4]);
    pa = pa/1e2; % convert to mb

    dtmdz = nan([length(grid.dim3.lon) length(grid.dim3.lat) length(par.pa) 12]);

    % calculate MALR following equations 3-5 in Stone and Carlson (1979)
    dalr = 9.8; % K/km
    eps = 0.622;
    R = 0.287; % J/g/K
    cp = 1.005; % J/g/K
    es0 = 6.11; % mb
    T0 = 273; % K
    L = 2510 - 2.38*(tai-T0);
    es = es0*exp(eps*L/R.*(1/T0-1./tai));
    desdT = eps*L.*es./(R*tai.^2);

    dtmdz = dalr * (1+eps*L.*es./(pa.*R.*tai))./(1+(eps.*L./(cp.*pa)).*(desdT));

    if any(strcmp(type, {'era5', 'erai'})); newdir=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s', type);
    elseif strcmp(type, 'gcm'); newdir=sprintf('/project2/tas1/miyawaki/projects/002/data/read/gcm/%s/%s', par.model, par.gcm.clim);
    elseif strcmp(type, 'echam'); newdir=sprintf('/project2/tas1/miyawaki/projects/002/data/read/echam/%s', par.echam.clim); end;
    if ~exist(newdir, 'dir'); mkdir(newdir); end
    filename='malr.mat';
    save(sprintf('%s/%s', newdir, filename), 'dtmdz', '-v7.3');

end
function make_malrz(type, par) % compute moist adiabatic lapse rate at every lat, time
% calculate moist adiabats
    if strcmp(type, 'era5') | strcmp(type, 'erai')
        foldername = sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/', type, par.lat_interp);
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s', type);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s', type);
    elseif strcmp(type, 'gcm')
        foldername = sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s/%s/', type, par.model, par.gcm.clim, par.lat_interp);
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s/%s', type, par.model, par.gcm.clim);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s', type, par.model, par.gcm.clim);
    elseif strcmp(type, 'echam')
        foldername = sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s/', type, par.echam.clim, par.lat_interp);
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s', type, par.echam.clim);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s', type, par.echam.clim);
    end
    load(sprintf('%s/grid.mat', prefix)); % read grid data
    load(sprintf('%s/tai.mat', prefix)); clear zgi; % read interpolated temperature

    tai = squeeze(nanmean(tai, 1)); % zonal average

    pa = repmat(par.pa', [1 length(grid.dim3.lat) 12]);
    pa = permute(pa, [2 1 3]);
    pa = pa/1e2; % convert to mb

    dtmdzz = nan([length(grid.dim3.lat) length(par.pa) 12]);

    % calculate MALR following equations 3-5 in Stone and Carlson (1979)
    dalr = 9.8; % K/km
    eps = 0.622;
    R = 0.287; % J/g/K
    cp = 1.005; % J/g/K
    es0 = 6.11; % mb
    T0 = 273; % K
    L = 2510 - 2.38*(tai-T0);
    es = es0*exp(eps*L/R.*(1/T0-1./tai));
    desdT = eps*L.*es./(R*tai.^2);

    dtmdzz = dalr * (1+eps*L.*es./(pa.*R.*tai))./(1+(eps.*L./(cp.*pa)).*(desdT));

    dtmdzz = repmat(dtmdzz, [1 1 1 length(grid.dim3.lon)]); % repeat to longitude
    dtmdzz = permute(dtmdzz, [4 1 2 3]); % bring lon to 1st

    if any(strcmp(type, {'era5', 'erai'})); newdir=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s', type);
    elseif strcmp(type, 'gcm'); newdir=sprintf('/project2/tas1/miyawaki/projects/002/data/read/gcm/%s/%s', par.model, par.gcm.clim);
    elseif strcmp(type, 'echam'); newdir=sprintf('/project2/tas1/miyawaki/projects/002/data/read/echam/%s', par.echam.clim); end;
    if ~exist(newdir, 'dir'); mkdir(newdir); end
    filename='malrz.mat';
    save(sprintf('%s/%s', newdir, filename), 'dtmdzz', '-v7.3');

end
function make_malrzsi(type, par) % compute moist adiabatic lapse rate at every lat, time
% calculate moist adiabats
    if strcmp(type, 'era5') | strcmp(type, 'erai')
        foldername = sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/', type, par.lat_interp);
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s', type);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s', type);
    elseif strcmp(type, 'gcm')
        foldername = sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s/%s/', type, par.model, par.gcm.clim, par.lat_interp);
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s/%s', type, par.model, par.gcm.clim);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s', type, par.model, par.gcm.clim);
    elseif strcmp(type, 'echam')
        foldername = sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s/', type, par.echam.clim, par.lat_interp);
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s', type, par.echam.clim);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s', type, par.echam.clim);
    end
    load(sprintf('%s/grid.mat', prefix)); % read grid data
    load(sprintf('%s/srfc.mat', prefix)); % read surface variable data
    load(sprintf('%s/tai.mat', prefix)); clear zgi; % read interpolated temperature

    tai = squeeze(nanmean(tai, 1)); % zonal average

    pa = repmat(par.pa', [1 length(grid.dim3.lat) 12]);
    pa = permute(pa, [2 1 3]);
    pa = pa/1e2; % convert to mb

    dtmdzz = nan([length(grid.dim3.lat) length(par.pa) 12]);

    % calculate MALR following equations 3-5 in Stone and Carlson (1979)
    dalr = 9.8; % K/km
    eps = 0.622;
    R = 0.287; % J/g/K
    cp = 1.005; % J/g/K
    es0 = 6.11; % mb
    T0 = 273; % K
    L = 2510 - 2.38*(tai-T0);
    es = es0*exp(eps*L/R.*(1/T0-1./tai));
    desdT = eps*L.*es./(R*tai.^2);

    dtmdzz = dalr * (1+eps*L.*es./(pa.*R.*tai))./(1+(eps.*L./(cp.*pa)).*(desdT));

    dtmdzz = repmat(dtmdzz, [1 1 1 length(grid.dim3.lon)]); % repeat to longitude
    dtmdzz = permute(dtmdzz, [4 1 2 3]); % bring lon to 1st

    if strcmp(type, 'era5') | strcmp(type, 'erai')
        ps_vert = repmat(srfc.sp, [1 1 1 size(dtmdzz, 3)]); % dims (lon x lat x time x plev)
        ps_vert = permute(ps_vert, [1 2 4 3]); % dims (lon x lat x plev x time)
        pa = double(permute(repmat(par.pa', [1 size(srfc.sp)]), [2 3 1 4]));
    elseif strcmp(type, 'gcm')
        ps_vert = repmat(srfc.ps, [1 1 1 size(dtmdzz, 3)]); % dims (lon x lat x time x plev)
        ps_vert = permute(ps_vert, [1 2 4 3]); % dims (lon x lat x plev x time)
        pa = permute(repmat(par.pa', [1 size(srfc.ps)]), [2 3 1 4]);
    elseif strcmp(type, 'echam')
        ps_vert = repmat(srfc.aps, [1 1 1 size(dtmdzz, 3)]); % dims (lon x lat x time x plev)
        ps_vert = permute(ps_vert, [1 2 4 3]); % dims (lon x lat x plev x time)
        pa = permute(repmat(par.pa', [1 size(srfc.aps)]), [2 3 1 4]);
    end

    pa = permute(pa, [3 1 2 4]); % bring height to 1st
    dtmdzz = permute(dtmdzz, [3 1 2 4]); % bring height to 1st
    pb = CmdLineProgressBar("Sorting and interpolating dtmdzz to new standard grid...");
    for lo=1:size(pa,2)
        pb.print(lo, size(pa,2));
        for la=1:size(pa,3)
            for mo=1:size(pa,4)
                dtmdzzsi(:,lo,la,mo) = interp1(pa(:,lo,la,mo)/ps_vert(lo,la,1,mo), dtmdzz(:,lo,la,mo), grid.dim3.si);
            end
        end
    end
    dtmdzzsi = permute(dtmdzzsi, [2 3 1 4]); % bring height to 3rd

    if any(strcmp(type, {'era5', 'erai'})); newdir=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s', type);
    elseif strcmp(type, 'gcm'); newdir=sprintf('/project2/tas1/miyawaki/projects/002/data/read/gcm/%s/%s', par.model, par.gcm.clim);
    elseif strcmp(type, 'echam'); newdir=sprintf('/project2/tas1/miyawaki/projects/002/data/read/echam/%s', par.echam.clim); end;
    if ~exist(newdir, 'dir'); mkdir(newdir); end
    filename='malrzsi.mat';
    save(sprintf('%s/%s', newdir, filename), 'dtmdzzsi', '-v7.3');

end
function make_dtmdz(type, par) % compute moist adiabat at every lat, lon, time
% calculate moist adiabats
    if strcmp(type, 'era5') | strcmp(type, 'erai')
        foldername = sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/', type, par.lat_interp);
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s', type);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s', type);
    elseif strcmp(type, 'gcm')
        foldername = sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s/%s/', type, par.model, par.gcm.clim, par.lat_interp);
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s/%s', type, par.model, par.gcm.clim);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s', type, par.model, par.gcm.clim);
    elseif strcmp(type, 'echam')
        foldername = sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s/', type, par.echam.clim, par.lat_interp);
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s', type, par.echam.clim);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s', type, par.echam.clim);
    end
    load(sprintf('%s/grid.mat', prefix)); % read grid data
    load(sprintf('%s/srfc.mat', prefix)); % read surface variable data

    dtmdz = nan([length(grid.dim3.lon) length(grid.dim3.lat) length(par.pa) 12]);

    if strcmp(type, 'era5') | strcmp(type, 'erai')
        pb = CmdLineProgressBar("Calculating moist adiabats...");
        for ilon = 1:length(grid.dim3.lon);
            pb.print(ilon, length(grid.dim3.lon)); % output progress of moist adiabat calculation
            for ilat = 1:length(grid.dim3.lat);
                for imon = 1:12;
                    for fn_var = fieldnames(srfc)'
                        ima.(fn_var{1}) = srfc.(fn_var{1})(ilon, ilat, imon);
                    end
                    dtmdz(ilon, ilat, :, imon) = calc_maz_dew_pa(ima, par.pa, par, type); % compute moist adiabat with RH
                end
            end
        end
    elseif strcmp(type, 'gcm')
        pb = CmdLineProgressBar("Calculating moist adiabats...");
        for ilon = 1:length(grid.dim3.lon);
            pb.print(ilon, length(grid.dim3.lon)); % output progress of moist adiabat calculation
            for ilat = 1:length(grid.dim3.lat);
                for imon = 1:12;
                    for fn_var = fieldnames(srfc)'
                        ima.(fn_var{1}) = srfc.(fn_var{1})(ilon, ilat, imon);
                    end
                    dtmdz(ilon, ilat, :, imon) = calc_maz_hurs_pa(ima, par.pa, par); % compute moist adiabat with RH
                end
            end
        end
    end

    if any(strcmp(type, {'era5', 'erai'})); newdir=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s', type);
    elseif strcmp(type, 'gcm'); newdir=sprintf('/project2/tas1/miyawaki/projects/002/data/read/gcm/%s/%s', par.model, par.gcm.clim);
    elseif strcmp(type, 'echam'); newdir=sprintf('/project2/tas1/miyawaki/projects/002/data/read/echam/%s', par.echam.clim); end;
    if ~exist(newdir, 'dir'); mkdir(newdir); end
    filename='dtmdz.mat';
    save(sprintf('%s/%s', newdir, filename), 'dtmdz', '-v7.3');

end
function make_dtmdzz(type, par) % compute moist adiabat at every lon, time
% calculate moist adiabats
    if strcmp(type, 'era5') | strcmp(type, 'erai')
        foldername = sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/', type, par.lat_interp);
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s', type);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s', type);
    elseif strcmp(type, 'gcm')
        foldername = sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s/%s', type, par.model, par.gcm.clim, par.lat_interp);
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s/%s', type, par.model, par.gcm.clim);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s', type, par.model, par.gcm.clim);
    elseif strcmp(type, 'echam')
        foldername = sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s/', type, par.echam.clim, par.lat_interp);
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s', type, par.echam.clim);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s', type, par.echam.clim);
    end
    load(sprintf('%s/grid.mat', prefix)); % read grid data
    load(sprintf('%s/srfc.mat', prefix)); % read surface variable data

    dtmdzz = nan([length(grid.dim3.lat) length(par.pa) 12]);

    for fn_var = fieldnames(srfc)'
        srfc.(fn_var{1}) = squeeze(nanmean(srfc.(fn_var{1}),1)); % zonal average
    end

    if strcmp(type, 'era5') | strcmp(type, 'erai')
        pb = CmdLineProgressBar("Calculating moist adiabats...");
        for ilat = 1:length(grid.dim3.lat);
            pb.print(ilat, length(grid.dim3.lat)); % output progress of moist adiabat calculation
            for imon = 1:12;
                for fn_var = fieldnames(srfc)'
                    ima.(fn_var{1}) = srfc.(fn_var{1})(ilat, imon);
                end
                dtmdzz(ilat, :, imon) = calc_maz_dew_pa(ima, par.pa, par, type); % compute moist adiabat with RH
            end
        end
    elseif strcmp(type, 'gcm')
        pb = CmdLineProgressBar("Calculating moist adiabats...");
        for ilat = 1:length(grid.dim3.lat);
            pb.print(ilat, length(grid.dim3.lat)); % output progress of moist adiabat calculation
            for imon = 1:12;
                for fn_var = fieldnames(srfc)'
                    ima.(fn_var{1}) = srfc.(fn_var{1})(ilat, imon);
                end
                dtmdzz(ilat, :, imon) = calc_maz_hurs_pa(ima, par.pa, par); % compute moist adiabat with RH
            end
        end
    end

    dtmdzz = repmat(dtmdzz, [1 1 1 length(grid.dim3.lon)]); % repeat in longitude
    dtmdzz = permute(dtmdzz, [4 1 2 3]); % bring lon to 1st

    if any(strcmp(type, {'era5', 'erai'})); newdir=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s', type);
    elseif strcmp(type, 'gcm'); newdir=sprintf('/project2/tas1/miyawaki/projects/002/data/read/gcm/%s/%s', par.model, par.gcm.clim);
    elseif strcmp(type, 'echam'); newdir=sprintf('/project2/tas1/miyawaki/projects/002/data/read/echam/%s', par.echam.clim); end;
    if ~exist(newdir, 'dir'); mkdir(newdir); end
    filename='dtmdzz.mat';
    save(sprintf('%s/%s', newdir, filename), 'dtmdzz', '-v7.3');

end
function make_dtmdzzsi(type, par) % compute moist adiabat at every lon, time
% calculate moist adiabats
    if strcmp(type, 'era5') | strcmp(type, 'erai')
        foldername = sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/', type, par.lat_interp);
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s', type);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s', type);
    elseif strcmp(type, 'gcm')
        foldername = sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s/%s', type, par.model, par.gcm.clim, par.lat_interp);
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s/%s', type, par.model, par.gcm.clim);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s', type, par.model, par.gcm.clim);
    elseif strcmp(type, 'echam')
        foldername = sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s/', type, par.echam.clim, par.lat_interp);
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s', type, par.echam.clim);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s', type, par.echam.clim);
    end
    load(sprintf('%s/grid.mat', prefix)); % read grid data
    load(sprintf('%s/srfc.mat', prefix)); % read surface variable data

    dtmdzzsi = nan([length(grid.dim3.lat) length(par.pa) 12]);

    for fn_var = fieldnames(srfc)'
        srfc.(fn_var{1}) = squeeze(nanmean(srfc.(fn_var{1}),1)); % zonal average
    end

    if strcmp(type, 'era5') | strcmp(type, 'erai')
        pb = CmdLineProgressBar("Calculating moist adiabats...");
        for ilat = 1:length(grid.dim3.lat);
            pb.print(ilat, length(grid.dim3.lat)); % output progress of moist adiabat calculation
            for imon = 1:12;
                for fn_var = fieldnames(srfc)'
                    ima.(fn_var{1}) = srfc.(fn_var{1})(ilat, imon);
                end
                dtmdzzsi(ilat, :, imon) = calc_maz_dew_si(ima, par.pa, par, type, grid); % compute moist adiabat with RH
            end
        end
    elseif strcmp(type, 'gcm')
        pb = CmdLineProgressBar("Calculating moist adiabats...");
        for ilat = 1:length(grid.dim3.lat);
            pb.print(ilat, length(grid.dim3.lat)); % output progress of moist adiabat calculation
            for imon = 1:12;
                for fn_var = fieldnames(srfc)'
                    ima.(fn_var{1}) = srfc.(fn_var{1})(ilat, imon);
                end
                dtmdzzsi(ilat, :, imon) = calc_maz_hurs_si(ima, par.pa, par, grid); % compute moist adiabat with RH
            end
        end
    elseif strcmp(type, 'echam')
        pb = CmdLineProgressBar("Calculating moist adiabats...");
        for ilat = 1:length(grid.dim3.lat);
            pb.print(ilat, length(grid.dim3.lat)); % output progress of moist adiabat calculation
            for imon = 1:12;
                for fn_var = fieldnames(srfc)'
                    ima.(fn_var{1}) = srfc.(fn_var{1})(ilat, imon);
                end
                dtmdzzsi(ilat, :, imon) = calc_maz_dew_si(ima, par.pa, par, type, grid); % compute moist adiabat with RH
            end
        end
    end

    dtmdzzsi = repmat(dtmdzzsi, [1 1 1 length(grid.dim3.lon)]); % repeat in longitude
    dtmdzzsi = permute(dtmdzzsi, [4 1 2 3]); % bring lon to 1st

    if any(strcmp(type, {'era5', 'erai'})); newdir=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s', type);
    elseif strcmp(type, 'gcm'); newdir=sprintf('/project2/tas1/miyawaki/projects/002/data/read/gcm/%s/%s', par.model, par.gcm.clim);
    elseif strcmp(type, 'echam'); newdir=sprintf('/project2/tas1/miyawaki/projects/002/data/read/echam/%s', par.echam.clim); end;
    if ~exist(newdir, 'dir'); mkdir(newdir); end
    filename='dtmdzzsi.mat';
    save(sprintf('%s/%s', newdir, filename), 'dtmdzzsi', '-v7.3');

end
function make_dtmdz_lat_plev(type, par) % compute moist adiabat as function of lat
% calculate moist adiabats
    if strcmp(type, 'era5') | strcmp(type, 'erai')
        foldername = sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/', type, par.lat_interp);
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s', type);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s', type);
    elseif strcmp(type, 'gcm')
        foldername = sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s/%s/', type, par.model, par.gcm.clim, par.lat_interp);
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s/%s', type, par.model, par.gcm.clim);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s', type, par.model, par.gcm.clim);
    elseif strcmp(type, 'echam')
        foldername = sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s/', type, par.echam.clim, par.lat_interp);
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s', type, par.echam.clim);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s', type, par.echam.clim);
    end
    load(sprintf('%s/grid.mat', prefix)); % read grid data
    load(sprintf('%s/srfc.mat', prefix)); % read surface variable data

    dtmdz_zt = nan([length(grid.dim3.lat) length(par.pa)]);

    for fn_var = fieldnames(srfc)'
        srfc.(fn_var{1}) = squeeze(nanmean(nanmean(srfc.(fn_var{1}),1),3)); % zonal and time mean
    end

    if strcmp(type, 'era5') | strcmp(type, 'erai')
        pb = CmdLineProgressBar("Calculating moist adiabats...");
        for ilat = 1:length(grid.dim3.lat);
            pb.print(ilat, length(grid.dim3.lat)); % output progress of moist adiabat calculation
            for fn_var = fieldnames(srfc)'
                ima.(fn_var{1}) = srfc.(fn_var{1})(ilat);
            end
            dtmdz_zt(ilat, :) = calc_maz_dew_pa(ima, par.pa, par, type); % compute moist adiabat with RH
        end
    elseif strcmp(type, 'gcm')
        pb = CmdLineProgressBar("Calculating moist adiabats...");
        for ilat = 1:length(grid.dim3.lat);
            pb.print(ilat, length(grid.dim3.lat)); % output progress of moist adiabat calculation
            for fn_var = fieldnames(srfc)'
                ima.(fn_var{1}) = srfc.(fn_var{1})(ilat);
            end
            dtmdz_zt(ilat, :) = calc_maz_hurs_pa(ima, par.pa, par); % compute moist adiabat with RH
        end
    elseif strcmp(type, 'gcm')
        pb = CmdLineProgressBar("Calculating moist adiabats...");
        for ilat = 1:length(grid.dim3.lat);
            pb.print(ilat, length(grid.dim3.lat)); % output progress of moist adiabat calculation
            for fn_var = fieldnames(srfc)'
                ima.(fn_var{1}) = srfc.(fn_var{1})(ilat);
            end
            dtmdz_zt(ilat, :) = calc_maz_dew_pa(ima, par.pa, par, type); % compute moist adiabat with RH
        end
    end

    if any(strcmp(type, {'era5', 'erai'})); newdir=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s', type);
    elseif strcmp(type, 'gcm'); newdir=sprintf('/project2/tas1/miyawaki/projects/002/data/read/gcm/%s/%s', par.model, par.gcm.clim);
    elseif strcmp(type, 'echam'); newdir=sprintf('/project2/tas1/miyawaki/projects/002/data/read/echam/%s', par.echam.clim); end;
    if ~exist(newdir, 'dir'); mkdir(newdir); end
    filename='dtmdz_zt.mat';
    save(sprintf('%s/%s', newdir, filename), 'dtmdz_zt', '-v7.3');

end
