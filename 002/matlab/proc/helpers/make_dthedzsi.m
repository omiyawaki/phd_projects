function make_dthedzsi(type, par)
    % compute model lapse rate in lat x plev x mon (add surface data after converting to sigma)
    
    newdir = make_savedir(type, par);
    prefix = make_prefix(type, par);
    prefix_proc = make_prefix_proc(type, par);

    [thetaeq_sm, sm, pa, ps_vert, thetaeqs, srfc, grid] = load_thetaeq(type, par);
    zg = load_zg(type, par);
    
    % load(sprintf('%s/grid.mat', prefix)); % read grid data
    % load(sprintf('%s/srfc.mat', prefix)); % read surface variable data
    % load(sprintf('%s/orog.mat', prefix)); % read surface orography
    
    if strcmp(type, 'gcm') & any(contains(par.model, {'GISS-E2-H', 'GISS-E2-R'})) % zg data in GISS-E2-H has an anomalous lat grid
        var = 'zg';
        file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/gcm/%s/%s_Amon_%s_%s_r1i1p1_*.ymonmean.nc', par.model, var, par.model, par.gcm.clim));
        fullpath=sprintf('%s/%s', file.folder, file.name);
        zg_lat = double(ncread(fullpath, 'lat'));
        zg = permute(zg, [2 1 3 4]);
        zg = interp1(zg_lat, zg, grid.dim3.lat);
        zg = permute(zg, [2 1 3 4]);
    end

    coarse_si = 1e-5*grid.dim3.plev;
    [~,si0] = max(coarse_si); % surface sigma index

    zg_sm = zg .* sm;
    zg_sm = permute(zg_sm, [3 1 2 4]);
    % convert to sigma coord.
    pb = CmdLineProgressBar("Sorting and interpolating dthedz to new standard grid...");
    for lo=1:size(pa,2)
        pb.print(lo, size(pa,2));
        for la=1:size(pa,3)
            for mo=1:size(pa,4)
            
                thetaeq_tmp = interp1(pa(:,lo,la,mo)./ps_vert(1,lo,la,mo), thetaeq_sm(:,lo,la,mo), coarse_si, 'linear', nan);
                zg_tmp = interp1(pa(:,lo,la,mo)./ps_vert(1,lo,la,mo), zg_sm(:,lo,la,mo), coarse_si, 'linear', nan);

                if all(isnan(thetaeq_tmp)) | all(isnan(zg_tmp))
                    thetaeq_si(:,lo,la,mo) = nan(size(grid.dim3.si));
                    zg_si(:,lo,la,mo) = nan(size(grid.dim3.si));
                else

                % add surface data
                thetaeq_tmp = add_thetaeqs(thetaeq_tmp, lo, la, mo, type, thetaeqs);
                zg_tmp(si0) = srfc.zs(lo,la,mo);
                
                % only keep nonnan data and redo interpolation
                notnan = find(~isnan(squeeze(thetaeq_tmp)) & ~isnan(squeeze(zg_tmp)));
                thetaeq_si(:,lo,la,mo) = interp1(coarse_si(notnan), thetaeq_tmp(notnan), grid.dim3.si, 'spline', nan);
                zg_si(:,lo,la,mo) = interp1(coarse_si(notnan), zg_tmp(notnan), grid.dim3.si, 'spline', nan);
    
                end

            end
        end
    end
    clear thetaeq_sm zg_sm

    %% calculate lapse rate before taking zonal average
    %dthedzsi = -1e3*(thetaeq_si(2:end,:,:,:)-thetaeq_si(1:end-1,:,:,:))./(zg_si(2:end,:,:,:)-zg_si(1:end-1,:,:,:)); % lapse rate in K/km
    %dthedzsi = interp1(1/2*(grid.dim3.si(2:end)+grid.dim3.si(1:end-1)), dthedzsi, grid.dim3.si);
    %dthedzsi(1,:,:,:) = -1e3*(thetaeq_si(2,:,:,:)-thetaeq_si(1,:,:,:))./(zg_si(2,:,:,:)-zg_si(1,:,:,:));
    %dthedzsi = permute(dthedzsi, [2 3 1 4]);

    % calculate lapse rate after taking zonal average
    thetaeq_si = squeeze(nanmean(thetaeq_si, 2)); thetaeq_si = repmat(thetaeq_si, [1 1 1 length(grid.dim3.lon)]); % now (lev lat time lon)
    zg_si = squeeze(nanmean(zg_si, 2)); zg_si = repmat(zg_si, [1 1 1 length(grid.dim3.lon)]);
    thetaeq_si = permute(thetaeq_si, [1 4 2 3]); zg_si = permute(zg_si, [1 4 2 3]); % bring back to (lev lon lat time)
    dthedzsi = 1e3*(thetaeq_si(2:end,:,:,:)-thetaeq_si(1:end-1,:,:,:))./(zg_si(2:end,:,:,:)-zg_si(1:end-1,:,:,:)); % lapse rate in K/km
    dthedzsi = interp1(1/2*(grid.dim3.si(2:end)+grid.dim3.si(1:end-1)), dthedzsi, grid.dim3.si);
    dthedzsi(1,:,:,:) = 1e3*(thetaeq_si(2,:,:,:)-thetaeq_si(1,:,:,:))./(zg_si(2,:,:,:)-zg_si(1,:,:,:));
    dthedzsi = permute(dthedzsi, [2 3 1 4]); % (lon lat lev time)
    
    filename='dthedzsi.mat';
    save(sprintf('%s/%s', newdir, filename), 'dthedzsi', '-v7.3');

end
