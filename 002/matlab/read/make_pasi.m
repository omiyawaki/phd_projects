function make_pasi(type, par)
    prefix = make_prefix(type, par);
    newdir = make_savedir(type, par);
    
    tmp = load(sprintf('%s/grid.mat', prefix)); grid=tmp.grid; clear tmp; % read grid data
    load(sprintf('%s/srfc.mat', prefix)); % load surface data
    
    % create surface mask
    [ps_vert, pa] = make_surfmask_vars(grid, type, srfc);
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
                    tmp(end) = srfc.sp(lo,la,mo); % surface is last element
                elseif strcmp(type, 'merra2')
                    tmp(1) = srfc.PS(lo,la,mo);
                elseif strcmp(type, 'jra55')
                    tmp(end) = srfc.ps(lo,la,mo); % surface is last element
                elseif strcmp(type, 'gcm')
                    tmp(1) = srfc.ps(lo,la,mo);
                elseif contains(type, 'echam')
                    tmp(1) = srfc.aps(lo,la,mo);
                elseif contains(type, 'hahn')
                    tmp(end) = srfc.PS(lo,la,mo); % surface is last element
                end

                % only keep nonnan data and redo interpolation
                notnan = find(~isnan(squeeze(tmp)));
                pa_si(:,lo,la,mo) = interp1(1e-5*grid.dim3.plev(notnan), tmp(notnan), grid.dim3.si, 'spline', nan);

            end
        end
    end

    pa_si = permute(pa_si, [2 3 1 4]); % reorder to lon x lat x si x mon

    if ~exist(newdir, 'dir'); mkdir(newdir); end
    filename='pa_si.mat';
    save(sprintf('%s/%s', newdir, filename), 'pa_si', '-v7.3');
end
