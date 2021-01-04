function make_tempsi(type, par)
    if any(strcmp(type, {'era5', 'era5c', 'erai'}))
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s', type, par.(type).yr_span);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s', type, par.(type).yr_span);
        file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/%s/temp/%s_temp_%s.ymonmean.nc', type, type, par.(type).yr_span));
        fullpath=sprintf('%s/%s', file.folder, file.name);
        ta_orig = double(ncread(fullpath, 't'));
    elseif strcmp(type, 'merra2')
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s', type, par.(type).yr_span);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s', type, par.(type).yr_span);
        file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/%s/temp/%s_temp_%s.ymonmean.nc', type, type, par.(type).yr_span));
        fullpath=sprintf('%s/%s', file.folder, file.name);
        ta_orig = double(ncread(fullpath, 'T'));
    elseif strcmp(type, 'jra55')
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s', type, par.(type).yr_span);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s', type, par.(type).yr_span);
        file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/%s/temp/%s_tmp_%s.ymonmean.nc', type, type, par.(type).yr_span));
        fullpath=sprintf('%s/%s', file.folder, file.name);
        ta_orig = double(ncread(fullpath, 'TMP_GDS0_ISBL_S123'));
    elseif strcmp(type, 'gcm')
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s/%s', type, par.model, par.gcm.clim);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s', type, par.model, par.gcm.clim);
        var = 'ta';
        file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/gcm/%s/%s_Amon_%s_%s_r1i1p1_*.ymonmean.nc', par.model, var, par.model, par.gcm.clim));
        fullpath=sprintf('%s/%s', file.folder, file.name);
        ta_orig = double(ncread(fullpath, var));
    elseif strcmp(type, 'echam')
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
    elseif contains(type, 'echam_pl')
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s', type);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s', type);
        var = 't';
        file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/%s/ATM_*.ymonmean.nc', type));
        fullpath=sprintf('%s/%s', file.folder, file.name);
        ta_orig = double(ncread(fullpath, var));
    end

    load(sprintf('%s/grid.mat', prefix)); % read grid data
    load(sprintf('%s/srfc.mat', prefix)); % load surface data

    % create surface mask
    if strcmp(type, 'era5') | strcmp(type, 'erai') | strcmp(type, 'era5c')
        ps_vert = repmat(srfc.sp, [1 1 1 size(ta_orig, 3)]); % dims (lon x lat x time x plev)
        ps_vert = permute(ps_vert, [1 2 4 3]); % dims (lon x lat x plev x time)
        pa = double(permute(repmat(grid.dim3.plev, [1 size(srfc.sp)]), [2 3 1 4]));
    elseif strcmp(type, 'merra2')
        ps_vert = repmat(srfc.PS, [1 1 1 size(ta_orig, 3)]); % dims (lon x lat x time x plev)
        ps_vert = permute(ps_vert, [1 2 4 3]); % dims (lon x lat x plev x time)
        pa = double(permute(repmat(grid.dim3.plev, [1 size(srfc.PS)]), [2 3 1 4]));
    elseif strcmp(type, 'jra55')
        ps_vert = repmat(srfc.ps, [1 1 1 size(ta_orig, 3)]); % dims (lon x lat x time x plev)
        ps_vert = permute(ps_vert, [1 2 4 3]); % dims (lon x lat x plev x time)
        pa = double(permute(repmat(grid.dim3.plev, [1 size(srfc.ps)]), [2 3 1 4]));
    elseif strcmp(type, 'gcm')
        ps_vert = repmat(srfc.ps, [1 1 1 size(ta_orig, 3)]); % dims (lon x lat x time x plev)
        ps_vert = permute(ps_vert, [1 2 4 3]); % dims (lon x lat x plev x time)
        pa = permute(repmat(grid.dim3.plev, [1 size(srfc.ps)]), [2 3 1 4]);
    elseif strcmp(type, 'echam')
        ps_vert = repmat(srfc.aps, [1 1 1 size(ta_orig, 3)]); % dims (lon x lat x time x plev)
        ps_vert = permute(ps_vert, [1 2 4 3]); % dims (lon x lat x plev x time)
        pa = permute(repmat(grid.dim3.plev, [1 size(srfc.aps)]), [2 3 1 4]);
    elseif contains(type, 'echam_pl')
        ps_vert = repmat(srfc.aps, [1 1 1 size(ta_orig, 3)]); % dims (lon x lat x time x plev)
        ps_vert = permute(ps_vert, [1 2 4 3]); % dims (lon x lat x plev x time)
        pa = permute(repmat(grid.dim3.plev, [1 size(srfc.aps)]), [2 3 1 4]);
    end
    sm = nan(size(ta_orig));
    sm(pa < 0.9961*ps_vert) = 1;
    ta_sm = ta_orig.*sm; % filter ta with surface mask

    ps_vert = permute(ps_vert, [3 1 2 4]);
    pa = permute(pa, [3 1 2 4]);
    ta_sm = permute(ta_sm, [3 1 2 4]);

    pb = CmdLineProgressBar("Sorting and interpolating temperature to new standard grid...");
    for lo=1:size(pa,2)
        pb.print(lo, size(pa,2));
        for la=1:size(pa,3)
            for mo=1:size(pa,4)
                tmp = interp1(pa(:,lo,la,mo)./ps_vert(1,lo,la,mo), ta_sm(:,lo,la,mo), 1e-5*grid.dim3.plev);

                % add surface data
                if strcmp(type, 'era5') | strcmp(type, 'erai') | strcmp(type, 'era5c')
                    tmp(end) = srfc.t2m(lo,la,mo); % surface is the last element in era plev
                elseif strcmp(type, 'merra2')
                    tmp(1) = srfc.T2M(lo,la,mo);
                elseif strcmp(type, 'jra55')
                    tmp(1) = srfc.tas(lo,la,mo);
                elseif strcmp(type, 'gcm')
                    tmp(1) = srfc.tas(lo,la,mo);
                elseif contains(type, 'echam')
                    tmp(1) = srfc.temp2(lo,la,mo);
                end

                % only keep nonnan data and redo interpolation
                notnan = find(~isnan(squeeze(tmp)));

                % ta_si.lin(:,lo,la,mo) = interp1(1e-5*grid.dim3.plev(notnan), tmp(notnan), grid.dim3.si, 'linear', nan);
                % ta_si.cub(:,lo,la,mo) = interp1(1e-5*grid.dim3.plev(notnan), tmp(notnan), grid.dim3.si, 'pchip', nan);
                ta_si.spl(:,lo,la,mo) = interp1(1e-5*grid.dim3.plev(notnan), tmp(notnan), grid.dim3.si, 'spline', nan);
                % ta_si.mak(:,lo,la,mo) = interp1(1e-5*grid.dim3.plev(notnan), tmp(notnan), grid.dim3.si, 'makima', nan);

                clear tmp

            end
        end
    end

    % ta_si.lin = permute(ta_si.lin, [2 3 1 4]); % reorder to lon x lat x si x mon
    % ta_si.cub = permute(ta_si.cub, [2 3 1 4]); % reorder to lon x lat x si x mon
    ta_si.spl = permute(ta_si.spl, [2 3 1 4]); % reorder to lon x lat x si x mon
    % ta_si.mak = permute(ta_si.mak, [2 3 1 4]); % reorder to lon x lat x si x mon

    % ta_si = permute(ta_si, [2 3 1 4]); % reorder to lon x lat x si x mon

    if any(strcmp(type, {'era5', 'era5c', 'erai'})); newdir=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s', type, par.(type).yr_span);
    elseif strcmp(type, 'merra2'); newdir=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s', type, par.(type).yr_span);
    elseif strcmp(type, 'jra55'); newdir=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s', type, par.(type).yr_span);
    elseif strcmp(type, 'gcm'); newdir=sprintf('/project2/tas1/miyawaki/projects/002/data/read/gcm/%s/%s', par.model, par.gcm.clim);
    elseif strcmp(type, 'echam'); newdir=sprintf('/project2/tas1/miyawaki/projects/002/data/read/echam/%s', par.echam.clim);
    elseif contains(type, 'echam_pl'); newdir=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s', type); end;
    if ~exist(newdir, 'dir'); mkdir(newdir); end
    filename='ta_si.mat';
    save(sprintf('%s/%s', newdir, filename), 'ta_si', '-v7.3');
end
