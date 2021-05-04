function make_dtdzsi(type, par)
    % compute model lapse rate in lat x plev x mon (add surface data after converting to sigma)
    
    newdir = make_savedir(type, par);
    prefix = make_prefix(type, par);
    prefix_proc = make_prefix_proc(type, par);
    temp = load_temp(type, par);
    zg = load_zg(type, par);
    
    load(sprintf('%s/grid.mat', prefix)); % read grid data
    load(sprintf('%s/srfc.mat', prefix)); % read surface variable data
    % load(sprintf('%s/orog.mat', prefix)); % read surface orography
    
    if strcmp(type, 'gcm') & any(contains(par.model, {'GISS-E2-H', 'GISS-E2-R'})) % zg data in GISS-E2-H has an anomalous lat grid
        var = 'zg';
        file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/gcm/%s/%s_Amon_%s_%s_r1i1p1_%s*.ymonmean.nc', par.model, var, par.model, par.(type).clim, par.(type).yr_span));
        fullpath=sprintf('%s/%s', file.folder, file.name);
        zg_lat = double(ncread(fullpath, 'lat'));
        zg = permute(zg, [2 1 3 4]);
        zg = interp1(zg_lat, zg, grid.dim3.lat);
        zg = permute(zg, [2 1 3 4]);
    end

    coarse_si = 1e-5*grid.dim3.plev;
    [~,si0] = max(coarse_si); % surface sigma index

    [ps_vert, pa] = make_surfmask_vars(grid, type, srfc);
    surface_mask = nan(size(pa));
    surface_mask(pa < ps_vert) = 1;

    ta_sm = temp .* surface_mask;
    zg_sm = zg .* surface_mask;

    ps_vert = permute(ps_vert, [3 1 2 4]);
    pa = permute(pa, [3 1 2 4]);
    ta_sm = permute(ta_sm, [3 1 2 4]);
    zg_sm = permute(zg_sm, [3 1 2 4]);
    % convert to sigma coord.
    pb = CmdLineProgressBar("Sorting and interpolating dtdz to new standard grid...");
    for lo=1:size(pa,2)
        pb.print(lo, size(pa,2));
        for la=1:size(pa,3)
            for mo=1:size(pa,4)
            
                ta_tmp = interp1(pa(:,lo,la,mo)./ps_vert(1,lo,la,mo), ta_sm(:,lo,la,mo), coarse_si, 'linear', nan);
                zg_tmp = interp1(pa(:,lo,la,mo)./ps_vert(1,lo,la,mo), zg_sm(:,lo,la,mo), coarse_si, 'linear', nan);

                if all(isnan(ta_tmp)) | all(isnan(zg_tmp))
                    ta_si(:,lo,la,mo) = nan(size(grid.dim3.si));
                    zg_si(:,lo,la,mo) = nan(size(grid.dim3.si));
                else

                    % add surface data
                    ta_tmp = add_tas(ta_tmp, lo, la, mo, type, srfc);
                    zg_tmp(si0) = srfc.zs(lo,la,mo);

                    % compute dtdz profile at nonnan levels
                    notnan = find(~isnan(squeeze(ta_tmp)) & ~isnan(squeeze(zg_tmp)));
                    ta_tmp = ta_tmp(notnan);
                    zg_tmp = zg_tmp(notnan);
                    dtdzsi_tmp = -1e3*(ta_tmp(2:end)-ta_tmp(1:end-1))./(zg_tmp(2:end)-zg_tmp(1:end-1)); 
                    si_in = coarse_si(notnan);
                    midpoints = 1/2*(si_in(2:end) + si_in(1:end-1));

                    dtdzsi(:,lo,la,mo) = interp1(midpoints, dtdzsi_tmp, grid.dim3.si, 'spline');
                    % if abs(dtdzsi(1,lo,la,mo)) > 50
                    %     disp(dtdzsi(1,lo,la,mo))
                    % end

                    clear notnan ta_tmp zg_tmp dtdzsi_tmp si_in midpoints
    
                end

            end
        end
    end
    clear ta_sm zg_sm

    dtdzsi = permute(dtdzsi, [2 3 1 4]); % (lon lat lev time)
    
    filename='dtdzsi.mat';
    save(sprintf('%s/%s', newdir, filename), 'dtdzsi', '-v7.3');

end
