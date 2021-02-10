function make_tempsi(type, par)
    ml_info
    
    [ta, lon, lat] = load_tempml(type, par);
    prefix = make_prefix(type, par);
    newdir = make_savedir(type, par);

    tmp = load(sprintf('%s/grid.mat', prefix)); grid=tmp.grid; clear tmp; % read grid data
    load(sprintf('%s/srfc.mat', prefix)); % load surface data
    
    % make sure model level ta has same lon x lev grid as pl data
    ta = interp1(lon, ta, grid.dim3.lon);
    ta = permute(ta, [2 1 3 4]);
    ta = interp1(lat, ta, grid.dim3.lat);
    ta = permute(ta, [2 1 3 4]);
    
    ta = permute(ta, [3 1 2 4]);
    ps = rename_ps(type, srfc);
    tas = rename_tas(type, srfc);
    clear srfc

    pb = CmdLineProgressBar("Sorting and interpolating temperature to new standard grid...");
    for lo=1:size(ta,2)
        pb.print(lo, size(ta,2));
        for la=1:size(ta,3)
            for mo=1:size(ta,4)
                % add surface data
                tmp = squeeze(ta(:,lo, la, mo));
                %tmp = [nan; tmp];
                %tmp(1) = tas(lo,la,mo);
                
                % only keep nonnan data and do interpolation
                notnan = find(~isnan(squeeze(tmp)));

                ta_si.spl(:,lo,la,mo) = interp1(ml.(type).a(notnan)/squeeze(ps(lo,la,mo))+ml.(type).b(notnan), tmp(notnan), grid.dim3.si, 'spline', nan); 
                
                clear tmp

            end
        end
    end

    ta_si.spl = permute(ta_si.spl, [2 3 1 4]); % reorder to lon x lat x si x mon

    % ta_si = permute(ta_si, [2 3 1 4]); % reorder to lon x lat x si x mon

    filename='taml_si.mat';
    save(sprintf('%s/%s', newdir, filename), 'ta_si', '-v7.3');
end
