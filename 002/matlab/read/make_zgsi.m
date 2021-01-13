function make_zgsi(type, par)
    if any(strcmp(type, {'era5', 'era5c', 'erai'}))
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s', type, par.(type).yr_span);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s', type, par.(type).yr_span);
        file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/%s/zg/%s_zg_%s.ymonmean.nc', type, type, par.(type).yr_span));
        fullpath=sprintf('%s/%s', file.folder, file.name);
        zg_orig = double(ncread(fullpath, 'z'));
    elseif contains(type, 'merra2')
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s', type, par.(type).yr_span);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s', type, par.(type).yr_span);
        file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/%s/zg/%s_zg_%s.ymonmean.nc', type, type, par.(type).yr_span));
        fullpath=sprintf('%s/%s', file.folder, file.name);
        zg_orig = double(ncread(fullpath, 'H'));
    elseif contains(type, 'jra55')
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s', type, par.(type).yr_span);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s', type, par.(type).yr_span);
        file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/%s/zg/%s_hgt_%s.ymonmean.nc', type, type, par.(type).yr_span));
        fullpath=sprintf('%s/%s', file.folder, file.name);
        zg_orig = double(ncread(fullpath, 'HGT_GDS0_ISBL_S123'));
    elseif strcmp(type, 'gcm')
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s/%s', type, par.model, par.gcm.clim);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s', type, par.model, par.gcm.clim);
        var = 'zg';
        file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/gcm/%s/%s_Amon_%s_%s_r1i1p1_*.ymonmean.nc', par.model, var, par.model, par.gcm.clim));
        fullpath=sprintf('%s/%s', file.folder, file.name);
        zg_orig = double(ncread(fullpath, var));
    elseif strcmp(type, 'echam')
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s', type, par.echam.clim);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s', type, par.echam.clim);
        var = 'geopoth';
        if contains(par.echam.clim, 'rp000')
            file=dir(sprintf('/project2/tas1/ockham/data11/tas/echam-aiv_rcc_6.1.00p1/%s/ATM_%s_0020_39.nc', par.echam.clim, par.echam.clim));
        else
            file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/echam/ATM*_%s_*.ymonmean.nc', par.echam.clim));
        end
        fullpath=sprintf('%s/%s', file.folder, file.name);
        zg_orig = double(ncread(fullpath, var));
    elseif contains(type, 'echam_pl')
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s', type);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s', type);
        var = 'geopoth';
        file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/%s/ATM_*.ymonmean.nc', type));
        fullpath=sprintf('%s/%s', file.folder, file.name);
        zg_orig = double(ncread(fullpath, var));
    end

    load(sprintf('%s/grid.mat', prefix)); % read grid data
    load(sprintf('%s/srfc.mat', prefix)); % load surface data

    if strcmp(type, 'gcm') & any(contains(par.model, {'GISS-E2-H', 'GISS-E2-R'})) % zg data in GISS-E2-H has an anomalous lat grid
        zg_lat = double(ncread(fullpath, 'lat'));
        zg_orig = permute(zg_orig, [2 1 3 4]);
        zg_orig = interp1(zg_lat, zg_orig, grid.dim3.lat);
        zg_orig = permute(zg_orig, [2 1 3 4]);
    end

    % create surface mask
    if strcmp(type, 'era5') | strcmp(type, 'erai') | strcmp(type, 'era5c')
        ps_vert = repmat(srfc.sp, [1 1 1 size(zg_orig, 3)]); % dims (lon x lat x time x plev)
        ps_vert = permute(ps_vert, [1 2 4 3]); % dims (lon x lat x plev x time)
        pa = double(permute(repmat(grid.dim3.plev, [1 size(srfc.sp)]), [2 3 1 4]));
    elseif strcmp(type, 'merra2')
        ps_vert = repmat(srfc.PS, [1 1 1 size(zg_orig, 3)]); % dims (lon x lat x time x plev)
        ps_vert = permute(ps_vert, [1 2 4 3]); % dims (lon x lat x plev x time)
        pa = double(permute(repmat(grid.dim3.plev, [1 size(srfc.PS)]), [2 3 1 4]));
    elseif strcmp(type, 'jra55')
        ps_vert = repmat(srfc.ps, [1 1 1 size(zg_orig, 3)]); % dims (lon x lat x time x plev)
        ps_vert = permute(ps_vert, [1 2 4 3]); % dims (lon x lat x plev x time)
        pa = double(permute(repmat(grid.dim3.plev, [1 size(srfc.ps)]), [2 3 1 4]));
    elseif strcmp(type, 'gcm')
        ps_vert = repmat(srfc.ps, [1 1 1 size(zg_orig, 3)]); % dims (lon x lat x time x plev)
        ps_vert = permute(ps_vert, [1 2 4 3]); % dims (lon x lat x plev x time)
        pa = permute(repmat(grid.dim3.plev, [1 size(srfc.ps)]), [2 3 1 4]);
    elseif contains(type, 'echam')
        ps_vert = repmat(srfc.aps, [1 1 1 size(zg_orig, 3)]); % dims (lon x lat x time x plev)
        ps_vert = permute(ps_vert, [1 2 4 3]); % dims (lon x lat x plev x time)
        pa = permute(repmat(grid.dim3.plev, [1 size(srfc.aps)]), [2 3 1 4]);
    end
    sm = nan(size(zg_orig));
    sm(pa < ps_vert) = 1;
    zg_sm = zg_orig.*sm; % filter zg with surface mask

    % % add tsurf dazg and interpolate to higher resolution vertical grid
    % [pa_plus zg_plus] = deal(nan([size(pa,1), size(pa,2) size(pa,3)+1 size(pa,4)])); % create empty grid with one extra vertical level
    % pa_plus(:,:,1:end-1,:) = pa; % populate with szgndard pressure grid
    % zg_plus(:,:,1:end-1,:) = zg_sm; % populate with standard geopotential
    % pa_plus(:,:,end,:) = ps_vert(:,:,1,:); % add surface pressure data into standard pressure grid
    % zg_plus(:,:,end,:) = srfc.zs(:,:,:); % add surface height data
    % pa_plus = permute(pa_plus, [3 1 2 4]); % bring plev dimension to front
    % zg_plus = permute(zg_plus, [3 1 2 4]); % bring plev dimension to front
    % [pa_plus sort_index] = sort(pa_plus, 1, 'descend'); % sort added surface pressure such that pressure decreases monotonically
    % zgi_sm = nan(length(par.si), size(pa, 1), size(pa, 2), size(pa, 4));
    % pb = CmdLineProgressBar("Sorting and interpolating zg to new standard grid...");
    % for lo=1:size(pa_plus,2)
    %     pb.print(lo, size(pa_plus,2));
    %     for la=1:size(pa_plus,3)
    %         for mo=1:size(pa_plus,4)
    %             zg_plus(:,lo,la,mo) = zg_plus(sort_index(:,lo,la,mo),lo,la,mo); % sort zg (has to be in loop because sort_index works for vector calls only)
    %             zgsi(:,lo,la,mo) = interp1(pa_plus(:,lo,la,mo)/ps_vert(lo,la,1,mo), zg_plus(:,lo,la,mo), grid.dim3.si);
    %         end
    %     end
    % end
    % clear pa_plus zg_plus; % clear unneeded variables

    ps_vert = permute(ps_vert, [3 1 2 4]);
    pa = permute(pa, [3 1 2 4]);
    zg_sm = permute(zg_sm, [3 1 2 4]);

    pb = CmdLineProgressBar("Sorting and interpolating temperature to new szgndard grid...");
    for lo=1:size(pa,2)
        pb.print(lo, size(pa,2));
        for la=1:size(pa,3)
            for mo=1:size(pa,4)
                tmp = interp1(pa(:,lo,la,mo)./ps_vert(1,lo,la,mo), zg_sm(:,lo,la,mo), 1e-5*grid.dim3.plev);

                if all(isnan(tmp))
                    zgsi(:,lo,la,mo) = nan(size(grid.dim3.si));
                else

                    if any(strcmp(type, {'era5', 'era5c', 'erai', 'jra55'}))
                        tmp(end) = srfc.zs(lo,la,mo); % add surface height (surface is last element)
                    else
                        tmp(1) = srfc.zs(lo,la,mo); % add surface height (surface is first element)
                    end

                    notnan = find(~isnan(squeeze(tmp))); % only keep nonnan data and redo interpolation

                    zgsi(:,lo,la,mo) = interp1(1e-5*grid.dim3.plev(notnan), tmp(notnan), grid.dim3.si, 'spline', nan);
                end

                clear tmp

            end
        end
    end

    zgsi = permute(zgsi, [2 3 1 4]); % reorder to lon x lat x si x mon

    if any(strcmp(type, {'era5', 'era5c', 'erai'})); newdir=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s', type, par.(type).yr_span);
    elseif strcmp(type, 'merra2'); newdir=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s', type, par.(type).yr_span);
    elseif strcmp(type, 'jra55'); newdir=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s', type, par.(type).yr_span);
    elseif strcmp(type, 'gcm'); newdir=sprintf('/project2/tas1/miyawaki/projects/002/data/read/gcm/%s/%s', par.model, par.gcm.clim);
    elseif strcmp(type, 'echam'); newdir=sprintf('/project2/tas1/miyawaki/projects/002/data/read/echam/%s', par.echam.clim);
    elseif contains(type, 'echam_pl'); newdir=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s', type); end;
    if ~exist(newdir, 'dir'); mkdir(newdir); end
    filename='zgsi.mat';
    save(sprintf('%s/%s', newdir, filename), 'zgsi', '-v7.3');
end
