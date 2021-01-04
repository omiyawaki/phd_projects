function make_psi(type, par)
    if any(strcmp(type, {'era5', 'era5c', 'erai'}))
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s', type, par.(type).yr_span);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s', type, par.(type).yr_span);
    elseif contains(type, 'merra2')
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s', type, par.(type).yr_span);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s', type, par.(type).yr_span);
    elseif contains(type, 'jra55')
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s', type, par.(type).yr_span);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s', type, par.(type).yr_span);
    elseif strcmp(type, 'gcm')
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s/%s', type, par.model, par.gcm.clim);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s', type, par.model, par.gcm.clim);
    elseif strcmp(type, 'echam')
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s', type, par.echam.clim);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s', type, par.echam.clim);
    elseif any(strcmp(type, {'echam_ml', 'echam_pl'}))
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s', type);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s', type);
        var = 'aps';
        file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/echam_ml/ATM_*.ymonmean.nc'));
        fullpath=sprintf('%s/%s', file.folder, file.name);
        ps_orig = double(ncread(fullpath, var));
    end

    load(sprintf('%s/grid.mat', prefix)); % read grid data
    load(sprintf('%s/srfc.mat', prefix)); % load surface data

    % create surface mask
    if strcmp(type, 'era5') | strcmp(type, 'erai') | strcmp(type, 'era5c')
        ps_vert = repmat(srfc.sp, [1 1 1 length(grid.dim3.plev)]); % dims (lon x lat x time x plev)
        ps_vert = permute(ps_vert, [1 2 4 3]); % dims (lon x lat x plev x time)
        pa = double(permute(repmat(grid.dim3.plev, [1 size(srfc.sp)]), [2 3 1 4]));
    elseif strcmp(type, 'merra2')
        ps_vert = repmat(srfc.PS, [1 1 1 length(grid.dim3.plev)]); % dims (lon x lat x time x plev)
        ps_vert = permute(ps_vert, [1 2 4 3]); % dims (lon x lat x plev x time)
        pa = double(permute(repmat(grid.dim3.plev, [1 size(srfc.PS)]), [2 3 1 4]));
    elseif strcmp(type, 'jra55')
        ps_vert = repmat(srfc.ps, [1 1 1 length(grid.dim3.plev)]); % dims (lon x lat x time x plev)
        ps_vert = permute(ps_vert, [1 2 4 3]); % dims (lon x lat x plev x time)
        pa = double(permute(repmat(grid.dim3.plev, [1 size(srfc.ps)]), [2 3 1 4]));
    elseif strcmp(type, 'gcm')
        ps_vert = repmat(srfc.ps, [1 1 1 length(grid.dim3.plev)]); % dims (lon x lat x time x plev)
        ps_vert = permute(ps_vert, [1 2 4 3]); % dims (lon x lat x plev x time)
        pa = permute(repmat(grid.dim3.plev, [1 size(srfc.ps)]), [2 3 1 4]);
    elseif strcmp(type, 'echam_ml')
        % compute sigma from a and b
        ps_vert = repmat(ps_orig, [1 1 1 length(grid.dim3.a)]); % dims (lon x lat x time x plev)
        ps_vert = permute(ps_vert, [1 2 4 3]); % dims (lon x lat x plev x time)
        a = permute(repmat(grid.dim3.a, [1 size(ps_orig)]), [2 3 1 4]);
        b = permute(repmat(grid.dim3.b, [1 size(ps_orig)]), [2 3 1 4]);
        pa = a + b.*ps_vert;
    elseif contains(type, 'echam')
        ps_vert = repmat(srfc.aps, [1 1 1 length(grid.dim3.plev)]); % dims (lon x lat x time x plev)
        ps_vert = permute(ps_vert, [1 2 4 3]); % dims (lon x lat x plev x time)
        pa = permute(repmat(grid.dim3.plev, [1 size(srfc.aps)]), [2 3 1 4]);
    end
    sm = nan(size(pa));
    sm(pa < ps_vert) = 1;
    pa_sm = pa.*sm; % filter pa with surface mask

    ps_vert = permute(ps_vert, [3 1 2 4]);
    pa = permute(pa, [3 1 2 4]);
    pa_sm = permute(pa_sm, [3 1 2 4]);

    pb = CmdLineProgressBar("Calculaing psi...");
    for lo=1:size(pa,2)
        pb.print(lo, size(pa,2));
        for la=1:size(pa,3)
            for mo=1:size(pa,4)
                tmp = interp1(pa(:,lo,la,mo)./ps_vert(1,lo,la,mo), pa_sm(:,lo,la,mo), 1e-5*grid.dim3.plev, 'linear');

                % add surface dapa
                if strcmp(type, 'era5') | strcmp(type, 'erai') | strcmp(type, 'era5c')
                    tmp(end) = srfc.sp(lo,la,mo);
                elseif strcmp(type, 'merra2')
                    tmp(1) = srfc.PS(lo,la,mo);
                elseif strcmp(type, 'jra55')
                    tmp(1) = srfc.ps(lo,la,mo);
                elseif strcmp(type, 'gcm')
                    tmp(1) = srfc.ps(lo,la,mo);
                elseif contains(type, 'echam')
                    tmp(1) = srfc.aps(lo,la,mo);
                end

                % only keep nonnan data and redo interpolation
                notnan = find(~isnan(squeeze(tmp)));
                pa_si(:,lo,la,mo) = interp1(1e-5*grid.dim3.plev(notnan), tmp(notnan), grid.dim3.si, 'spline', nan);

            end
        end
    end

    pa_si = permute(pa_si, [2 3 1 4]); % reorder to lon x lat x si x mon

    if any(strcmp(type, {'era5', 'era5c', 'erai'})); newdir=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s', type, par.(type).yr_span);
    elseif strcmp(type, 'merra2'); newdir=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s', type, par.(type).yr_span);
    elseif strcmp(type, 'jra55'); newdir=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s', type, par.(type).yr_span);
    elseif strcmp(type, 'gcm'); newdir=sprintf('/project2/tas1/miyawaki/projects/002/data/read/gcm/%s/%s', par.model, par.gcm.clim);
    elseif strcmp(type, 'echam'); newdir=sprintf('/project2/tas1/miyawaki/projects/002/data/read/echam/%s', par.echam.clim);
    elseif any(strcmp(type, {'echam_ml','echam_pl'})); newdir=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s', type); end;
    if ~exist(newdir, 'dir'); mkdir(newdir); end
    filename='pa_si.mat';
    save(sprintf('%s/%s', newdir, filename), 'pa_si', '-v7.3');
end
