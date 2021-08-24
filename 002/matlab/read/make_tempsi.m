function make_tempml(type, par)
    ta_orig = load_temp(type, par);
    prefix = make_prefix(type, par);
    newdir = make_savedir(type, par);

    tmp = load(sprintf('%s/grid.mat', prefix)); grid=tmp.grid; clear tmp; % read grid data
    load(sprintf('%s/srfc.mat', prefix)); % load surface data

    % create surface mask
    [ps_vert, pa] = make_surfmask_vars(grid, type, srfc);
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
                if strcmp(type, 'hahn')
                    plevp1 = 1e-5*[grid.dim3.plev; 1e5];
                    tmp = interp1(pa(:,lo,la,mo)./ps_vert(1,lo,la,mo), ta_sm(:,lo,la,mo), plevp1); % include 10000 Pa level to Hahn
                else
                    tmp = interp1(pa(:,lo,la,mo)./ps_vert(1,lo,la,mo), ta_sm(:,lo,la,mo), 1e-5*grid.dim3.plev);
                end

                % add surface data
                tmp = add_tas(tmp, lo, la, mo, type, srfc);

                % only keep nonnan data and redo interpolation
                notnan = find(~isnan(squeeze(tmp)));

                if strcmp(type, 'hahn')
                    ta_si.spl(:,lo,la,mo) = interp1(plevp1(notnan), tmp(notnan), grid.dim3.si, 'spline', nan);
                else
                    % ta_si.lin(:,lo,la,mo) = interp1(1e-5*grid.dim3.plev(notnan), tmp(notnan), grid.dim3.si, 'linear', nan);
                    % ta_si.cub(:,lo,la,mo) = interp1(1e-5*grid.dim3.plev(notnan), tmp(notnan), grid.dim3.si, 'pchip', nan);
                    ta_si.spl(:,lo,la,mo) = interp1(1e-5*grid.dim3.plev(notnan), tmp(notnan), grid.dim3.si, 'spline', nan);
                    % ta_si.mak(:,lo,la,mo) = interp1(1e-5*grid.dim3.plev(notnan), tmp(notnan), grid.dim3.si, 'makima', nan);
                end

                clear tmp

            end
        end
    end

    % ta_si.lin = permute(ta_si.lin, [2 3 1 4]); % reorder to lon x lat x si x mon
    % ta_si.cub = permute(ta_si.cub, [2 3 1 4]); % reorder to lon x lat x si x mon
    ta_si.spl = permute(ta_si.spl, [2 3 1 4]); % reorder to lon x lat x si x mon
    % ta_si.mak = permute(ta_si.mak, [2 3 1 4]); % reorder to lon x lat x si x mon

    % ta_si = permute(ta_si, [2 3 1 4]); % reorder to lon x lat x si x mon

    filename='ta_si.mat';
    save(sprintf('%s/%s', newdir, filename), 'ta_si', '-v7.3');
end
