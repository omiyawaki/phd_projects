function make_thetaeqsi(type, par)

    savedir = make_savedir(type, par);

    [thetaeq_sm, sm, pa, ps_vert, thetaeqs, srfc, grid] = load_thetaeq(type, par);

    pb = CmdLineProgressBar("Sorting and interpolating equivalent potential temperature to new standard grid...");
    for lo=1:size(pa,2)
        pb.print(lo, size(pa,2));
        for la=1:size(pa,3)
            for mo=1:size(pa,4)
                tmp = interp1(pa(:,lo,la,mo)./ps_vert(1,lo,la,mo), thetaeq_sm(:,lo,la,mo), 1e-5*grid.dim3.plev);

                % add surface data
                % tmp(1) = thetaeqs(lo,la,mo);
                tmp = add_thetaeqs(tmp, lo, la, mo, type, thetaeqs);

                % only keep nonnan dathetaeq and redo interpolation
                notnan = find(~isnan(squeeze(tmp)));

                % thetaeq_si.lin(:,lo,la,mo) = interp1(1e-5*grid.dim3.plev(notnan), tmp(notnan), grid.dim3.si, 'linear', nan);
                % thetaeq_si.cub(:,lo,la,mo) = interp1(1e-5*grid.dim3.plev(notnan), tmp(notnan), grid.dim3.si, 'pchip', nan);
                thetaeq_si.spl(:,lo,la,mo) = interp1(1e-5*grid.dim3.plev(notnan), tmp(notnan), grid.dim3.si, 'spline', nan);
                % thetaeq_si.mak(:,lo,la,mo) = interp1(1e-5*grid.dim3.plev(notnan), tmp(notnan), grid.dim3.si, 'makima', nan);

                clear tmp

            end
        end
    end

    % thetaeq_si.lin = permute(thetaeq_si.lin, [2 3 1 4]); % reorder to lon x lat x si x mon
    % thetaeq_si.cub = permute(thetaeq_si.cub, [2 3 1 4]); % reorder to lon x lat x si x mon
    thetaeq_si.spl = permute(thetaeq_si.spl, [2 3 1 4]); % reorder to lon x lat x si x mon
    % thetaeq_si.mak = permute(thetaeq_si.mak, [2 3 1 4]); % reorder to lon x lat x si x mon

    % thetaeq_si = permute(thetaeq_si, [2 3 1 4]); % reorder to lon x lat x si x mon

    filename='thetaeq_si.mat';
    save(sprintf('%s/%s', savedir, filename), 'thetaeq_si', '-v7.3');
end
