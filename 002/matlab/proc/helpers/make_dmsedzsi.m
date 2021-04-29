function make_dmsedzsi(type, par)
    % compute model lapse rate in lat x plev x mon (add surface data after converting to sigma)
    
    newdir = make_savedir(type, par);
    prefix = make_prefix(type, par);
    prefix_proc = make_prefix_proc(type, par);

    [mse_sm, sm, pa, ps_vert, mses, srfc, grid] = load_mse(type, par);
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
    pb = CmdLineProgressBar("Sorting and interpolating dmsedz to new standard grid...");
    for lo=1:size(pa,2)
        pb.print(lo, size(pa,2));
        for la=1:size(pa,3)
            for mo=1:size(pa,4)
            
                mse_tmp = interp1(pa(:,lo,la,mo)./ps_vert(1,lo,la,mo), mse_sm(:,lo,la,mo), coarse_si, 'linear', nan);
                zg_tmp = interp1(pa(:,lo,la,mo)./ps_vert(1,lo,la,mo), zg_sm(:,lo,la,mo), coarse_si, 'linear', nan);

                if all(isnan(mse_tmp)) | all(isnan(zg_tmp))
                    mse_si(:,lo,la,mo) = nan(size(grid.dim3.si));
                    zg_si(:,lo,la,mo) = nan(size(grid.dim3.si));
                else

                % add surface data
                mse_tmp = add_mses(mse_tmp, lo, la, mo, type, mses);
                zg_tmp(si0) = srfc.zs(lo,la,mo);
                
                % only keep nonnan data and redo interpolation
                notnan = find(~isnan(squeeze(mse_tmp)) & ~isnan(squeeze(zg_tmp)));
                mse_si(:,lo,la,mo) = interp1(coarse_si(notnan), mse_tmp(notnan), grid.dim3.si, 'spline', nan);
                zg_si(:,lo,la,mo) = interp1(coarse_si(notnan), zg_tmp(notnan), grid.dim3.si, 'spline', nan);
    
                end

            end
        end
    end
    clear mse_sm zg_sm

    %% calculate lapse rate before taking zonal average
    %dmsedzsi = -1e3*(mse_si(2:end,:,:,:)-mse_si(1:end-1,:,:,:))./(zg_si(2:end,:,:,:)-zg_si(1:end-1,:,:,:)); % lapse rate in K/km
    %dmsedzsi = interp1(1/2*(grid.dim3.si(2:end)+grid.dim3.si(1:end-1)), dmsedzsi, grid.dim3.si);
    %dmsedzsi(1,:,:,:) = -1e3*(mse_si(2,:,:,:)-mse_si(1,:,:,:))./(zg_si(2,:,:,:)-zg_si(1,:,:,:));
    %dmsedzsi = permute(dmsedzsi, [2 3 1 4]);

    % calculate lapse rate after taking zonal average
    mse_si = squeeze(nanmean(mse_si, 2)); mse_si = repmat(mse_si, [1 1 1 length(grid.dim3.lon)]); % now (lev lat time lon)
    zg_si = squeeze(nanmean(zg_si, 2)); zg_si = repmat(zg_si, [1 1 1 length(grid.dim3.lon)]);
    mse_si = permute(mse_si, [1 4 2 3]); zg_si = permute(zg_si, [1 4 2 3]); % bring back to (lev lon lat time)
    dmsedzsi = (mse_si(2:end,:,:,:)-mse_si(1:end-1,:,:,:))./(zg_si(2:end,:,:,:)-zg_si(1:end-1,:,:,:)); % lapse rate in kJ/km
    dmsedzsi = interp1(1/2*(grid.dim3.si(2:end)+grid.dim3.si(1:end-1)), dmsedzsi, grid.dim3.si);
    dmsedzsi(1,:,:,:) = (mse_si(2,:,:,:)-mse_si(1,:,:,:))./(zg_si(2,:,:,:)-zg_si(1,:,:,:));
    dmsedzsi = permute(dmsedzsi, [2 3 1 4]); % (lon lat lev time)
    
    filename='dmsedzsi.mat';
    save(sprintf('%s/%s', newdir, filename), 'dmsedzsi', '-v7.3');

end
