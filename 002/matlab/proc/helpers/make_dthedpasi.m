function make_dthedpasi(type, par)
    % compute model lapse rate in lat x plev x mon (add surface data after converting to sigma)
    
    newdir = make_savedir(type, par);
    prefix = make_prefix(type, par);
    prefix_proc = make_prefix_proc(type, par);

    [thetaeq_sm, sm, pa, ps_vert, thetaeqs, srfc, grid] = load_thetaeq(type, par);
    
    coarse_si = 1e-5*grid.dim3.plev;
    [~,si0] = max(coarse_si); % surface sigma index

    % convert to sigma coord.
    pb = CmdLineProgressBar("Sorting and interpolating dthedz to new standard grid...");
    for lo=1:size(pa,2)
        pb.print(lo, size(pa,2));
        for la=1:size(pa,3)
            for mo=1:size(pa,4)
            
                thetaeq_tmp = interp1(pa(:,lo,la,mo)./ps_vert(1,lo,la,mo), thetaeq_sm(:,lo,la,mo), coarse_si, 'linear', nan);
                pa_tmp = interp1(pa(:,lo,la,mo)./ps_vert(1,lo,la,mo), pa(:,lo,la,mo), coarse_si, 'linear', nan);

                if all(isnan(thetaeq_tmp)) | all(isnan(pa_tmp))
                    thetaeq_si(:,lo,la,mo) = nan(size(grid.dim3.si));
                    pa_si(:,lo,la,mo) = nan(size(grid.dim3.si));
                else

                % add surface data
                thetaeq_tmp = add_thetaeqs(thetaeq_tmp, lo, la, mo, type, thetaeqs);
                pa_tmp(si0) = ps_vert(1,lo,la,mo);
                
                % only keep nonnan data and redo interpolation
                notnan = find(~isnan(squeeze(thetaeq_tmp)) & ~isnan(squeeze(pa_tmp)));
                thetaeq_si(:,lo,la,mo) = interp1(coarse_si(notnan), thetaeq_tmp(notnan), grid.dim3.si, 'spline', nan);
                pa_si(:,lo,la,mo) = interp1(coarse_si(notnan), pa_tmp(notnan), grid.dim3.si, 'spline', nan);
    
                end

            end
        end
    end
    clear thetaeq_sm pa

    %% calculate lapse rate before taking zonal average
    %dthedpasi = -1e3*(thetaeq_si(2:end,:,:,:)-thetaeq_si(1:end-1,:,:,:))./(zg_si(2:end,:,:,:)-zg_si(1:end-1,:,:,:)); % lapse rate in K/km
    %dthedpasi = interp1(1/2*(grid.dim3.si(2:end)+grid.dim3.si(1:end-1)), dthedpasi, grid.dim3.si);
    %dthedpasi(1,:,:,:) = -1e3*(thetaeq_si(2,:,:,:)-thetaeq_si(1,:,:,:))./(zg_si(2,:,:,:)-zg_si(1,:,:,:));
    %dthedpasi = permute(dthedpasi, [2 3 1 4]);

    % calculate lapse rate after taking zonal average
    thetaeq_si = squeeze(nanmean(thetaeq_si, 2)); thetaeq_si = repmat(thetaeq_si, [1 1 1 length(grid.dim3.lon)]); % now (lev lat time lon)
    pa_si = squeeze(nanmean(pa_si, 2)); pa_si = repmat(pa_si, [1 1 1 length(grid.dim3.lon)]);
    thetaeq_si = permute(thetaeq_si, [1 4 2 3]); pa_si = permute(pa_si, [1 4 2 3]); % bring back to (lev lon lat time)
    dthedpasi = -1e2*(thetaeq_si(2:end,:,:,:)-thetaeq_si(1:end-1,:,:,:))./(pa_si(2:end,:,:,:)-pa_si(1:end-1,:,:,:)); % lapse rate in K/hPa
    dthedpasi = interp1(1/2*(grid.dim3.si(2:end)+grid.dim3.si(1:end-1)), dthedpasi, grid.dim3.si);
    dthedpasi(1,:,:,:) = -1e2*(thetaeq_si(2,:,:,:)-thetaeq_si(1,:,:,:))./(pa_si(2,:,:,:)-pa_si(1,:,:,:));
    dthedpasi = permute(dthedpasi, [2 3 1 4]); % (lon lat lev time)
    
    filename='dthedpasi.mat';
    save(sprintf('%s/%s', newdir, filename), 'dthedpasi', '-v7.3');

end
