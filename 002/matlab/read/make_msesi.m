function make_msesi(type, par)

    savedir = make_savedir(type, par);

    [mse_sm, sm, pa, ps_vert, mses, srfc, grid] = load_mse(type, par);

    pb = CmdLineProgressBar("Sorting and interpolating equivalent potential temperature to new standard grid...");
    for lo=1:size(pa,2)
        pb.print(lo, size(pa,2));
        for la=1:size(pa,3)
            for mo=1:size(pa,4)
                tmp = interp1(pa(:,lo,la,mo)./ps_vert(1,lo,la,mo), mse_sm(:,lo,la,mo), 1e-5*grid.dim3.plev);

                % add surface data
                % tmp(1) = mses(lo,la,mo);
                tmp = add_mses(tmp, lo, la, mo, type, mses);

                % only keep nonnan damse and redo interpolation
                notnan = find(~isnan(squeeze(tmp)));

                % mse_si.lin(:,lo,la,mo) = interp1(1e-5*grid.dim3.plev(notnan), tmp(notnan), grid.dim3.si, 'linear', nan);
                % mse_si.cub(:,lo,la,mo) = interp1(1e-5*grid.dim3.plev(notnan), tmp(notnan), grid.dim3.si, 'pchip', nan);
                mse_si.spl(:,lo,la,mo) = interp1(1e-5*grid.dim3.plev(notnan), tmp(notnan), grid.dim3.si, 'spline', nan);
                % mse_si.mak(:,lo,la,mo) = interp1(1e-5*grid.dim3.plev(notnan), tmp(notnan), grid.dim3.si, 'makima', nan);

                clear tmp

            end
        end
    end

    % mse_si.lin = permute(mse_si.lin, [2 3 1 4]); % reorder to lon x lat x si x mon
    % mse_si.cub = permute(mse_si.cub, [2 3 1 4]); % reorder to lon x lat x si x mon
    mse_si.spl = permute(mse_si.spl, [2 3 1 4]); % reorder to lon x lat x si x mon
    % mse_si.mak = permute(mse_si.mak, [2 3 1 4]); % reorder to lon x lat x si x mon

    % mse_si = permute(mse_si, [2 3 1 4]); % reorder to lon x lat x si x mon

    filename='mse_si.mat';
    save(sprintf('%s/%s', savedir, filename), 'mse_si', '-v7.3');
end
