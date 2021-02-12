function make_zgsi(type, par)
    zg_orig = load_zg(type, par);
    prefix = make_prefix(type, par);
    newdir = make_savedir(type, par);
    
    tmp = load(sprintf('%s/grid.mat', prefix)); grid=tmp.grid; clear tmp; % read grid data
    load(sprintf('%s/srfc.mat', prefix)); % load surface data

    if strcmp(type, 'gcm') & any(contains(par.model, {'GISS-E2-H', 'GISS-E2-R'})) % zg data in GISS-E2-H has an anomalous lat grid
        zg_lat = double(ncread(fullpath, 'lat'));
        zg_orig = permute(zg_orig, [2 1 3 4]);
        zg_orig = interp1(zg_lat, zg_orig, grid.dim3.lat);
        zg_orig = permute(zg_orig, [2 1 3 4]);
    end

    % create surface mask
    [ps_vert, pa] = make_surfmask_vars(grid, type, srfc);
    sm = nan(size(zg_orig));
    sm(pa < ps_vert) = 1;
    zg_sm = zg_orig.*sm; % filter zg with surface mask

    ps_vert = permute(ps_vert, [3 1 2 4]);
    pa = permute(pa, [3 1 2 4]);
    zg_sm = permute(zg_sm, [3 1 2 4]);

    pb = CmdLineProgressBar("Sorting and interpolating zg to new standard grid...");
    for lo=1:size(pa,2)
        pb.print(lo, size(pa,2));
        for la=1:size(pa,3)
            for mo=1:size(pa,4)
                tmp = interp1(pa(:,lo,la,mo)./ps_vert(1,lo,la,mo), zg_sm(:,lo,la,mo), 1e-5*grid.dim3.plev);

                if all(isnan(tmp))
                    zgsi(:,lo,la,mo) = nan(size(grid.dim3.si));
                else

                    if any(strcmp(type, {'era5', 'era5c', 'erai', 'jra55', 'hahn'}))
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

    filename='zgsi.mat';
    save(sprintf('%s/%s', newdir, filename), 'zgsi', '-v7.3');
end
