function make_thetaeqsi(type, par)
    if any(strcmp(type, {'era5', 'era5c', 'erai'}))
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s', type);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s', type);
        file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/%s/temp/%s_temp_%s.ymonmean.nc', type, type, par.(type).yr_span));
        fullpath=sprintf('%s/%s', file.folder, file.name);
        ta_orig = double(ncread(fullpath, 't'));
    elseif strcmp(type, 'merra2')
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s', type);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s', type);
        file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/%s/temp/%s_temp_%s.ymonmean.nc', type, type, par.(type).yr_span));
        fullpath=sprintf('%s/%s', file.folder, file.name);
        ta_orig = double(ncread(fullpath, 'T'));
    elseif strcmp(type, 'gcm')
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s/%s', type, par.model, par.gcm.clim);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s', type, par.model, par.gcm.clim);
        var = 'ta';
        file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/gcm/%s/%s_Amon_%s_%s_r1i1p1_*.ymonmean.nc', par.model, var, par.model, par.gcm.clim));
        fullpath=sprintf('%s/%s', file.folder, file.name);
        ta_orig = double(ncread(fullpath, var));
        var = 'hur';
        file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/gcm/%s/%s_Amon_%s_%s_r1i1p1_*.ymonmean.nc', par.model, var, par.model, par.gcm.clim));
        fullpath=sprintf('%s/%s', file.folder, file.name);
        hur_orig = double(ncread(fullpath, var));
    elseif strcmp(type, 'echam')
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s', type, par.echam.clim);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s', type, par.echam.clim);
        var = 't';
        if contains(par.echam.clim, 'rp000')
            file=dir(sprintf('/project2/tas1/ockham/data11/tas/echam-aiv_rcc_6.1.00p1/%s/ATM_%s_0020_39.nc', par.echam.clim, par.echam.clim));
        else
            file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/echam/ATM*_%s_*.ymonmean.nc', par.echam.clim));
        end
        fullpath=sprintf('%s/%s', file.folder, file.name);
        ta_orig = double(ncread(fullpath, var));
    elseif contains(type, 'echam_pl')
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s', type);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s', type);
        var = 't';
        file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/%s/ATM_*.ymonmean.nc', type));
        fullpath=sprintf('%s/%s', file.folder, file.name);
        ta_orig = double(ncread(fullpath, var));
    end

    load(sprintf('%s/grid.mat', prefix)); % read grid data
    load(sprintf('%s/srfc.mat', prefix)); % load surface data

    esat = calc_esat(ta_orig, 0); % compute saturation vapor pressure
    p = permute(repmat(grid.dim3.plev, [1 size(ta_orig,1) size(ta_orig,2) size(ta_orig,4)]), [2 3 1 4]);
    e = esat.*hur_orig/100;
    r = calc_r(p, e, par);
    pd = p - e; % partial pressure of dry air
    clear p e esat;

    % compute eq potential temperature following AMS glossary definition
    thetaeq = ta_orig .* (1e5./pd).^(par.Rd/par.cpd).*(hur_orig/100).^(-r*par.Rv/par.cpd).*exp(par.L*r./(par.cpd*ta_orig));
    clear ta_orig hur_orig;

    % surface eq pot temp
    if strcmp(type, 'era5') | strcmp(type, 'erai') | strcmp(type, 'era5c')
        tas = srfc.t2m; % surface is the last element in era plev
        ps = srfc.sp;
    elseif strcmp(type, 'merra2')
        tas = srfc.T2M;
        ps = srfc.PS;
    elseif strcmp(type, 'gcm')
        tas = srfc.tas;
        ps = srfc.ps;
        hurs = srfc.hurs;
    elseif contains(type, 'echam')
        tas = srfc.temp2;
    end

    esats = calc_esat(tas, 0); % compute saturation vapor pressure
    es = esats.*hurs/100;
    rs = calc_r(ps, es, par);
    pds = ps - es; % partial pressure of dry air
    clear ps es esats;

    thetaeqs = tas .* (1e5./pds).^(par.Rd/par.cpd).*(hurs/100).^(-rs*par.Rv/par.cpd).*exp(par.L*rs./(par.cpd*tas));
    clear tas hurs;

    % create surface mask
    if strcmp(type, 'era5') | strcmp(type, 'erai') | strcmp(type, 'era5c')
        ps_vert = repmat(srfc.sp, [1 1 1 size(thetaeq, 3)]); % dims (lon x lat x time x plev)
        ps_vert = permute(ps_vert, [1 2 4 3]); % dims (lon x lat x plev x time)
        pa = double(permute(repmat(grid.dim3.plev, [1 size(srfc.sp)]), [2 3 1 4]));
    elseif strcmp(type, 'merra2')
        ps_vert = repmat(srfc.PS, [1 1 1 size(thetaeq, 3)]); % dims (lon x lat x time x plev)
        ps_vert = permute(ps_vert, [1 2 4 3]); % dims (lon x lat x plev x time)
        pa = double(permute(repmat(grid.dim3.plev, [1 size(srfc.PS)]), [2 3 1 4]));
    elseif strcmp(type, 'gcm')
        ps_vert = repmat(srfc.ps, [1 1 1 size(thetaeq, 3)]); % dims (lon x lat x time x plev)
        ps_vert = permute(ps_vert, [1 2 4 3]); % dims (lon x lat x plev x time)
        pa = permute(repmat(grid.dim3.plev, [1 size(srfc.ps)]), [2 3 1 4]);
    elseif strcmp(type, 'echam')
        ps_vert = repmat(srfc.aps, [1 1 1 size(thetaeq, 3)]); % dims (lon x lat x time x plev)
        ps_vert = permute(ps_vert, [1 2 4 3]); % dims (lon x lat x plev x time)
        pa = permute(repmat(grid.dim3.plev, [1 size(srfc.aps)]), [2 3 1 4]);
    elseif contains(type, 'echam_pl')
        ps_vert = repmat(srfc.aps, [1 1 1 size(thetaeq, 3)]); % dims (lon x lat x time x plev)
        ps_vert = permute(ps_vert, [1 2 4 3]); % dims (lon x lat x plev x time)
        pa = permute(repmat(grid.dim3.plev, [1 size(srfc.aps)]), [2 3 1 4]);
    end
    sm = nan(size(thetaeq));
    sm(pa < 0.9961*ps_vert) = 1;
    thetaeq_sm = thetaeq.*sm; % filter thetaeq with surface mask

    ps_vert = permute(ps_vert, [3 1 2 4]);
    pa = permute(pa, [3 1 2 4]);
    thetaeq_sm = permute(thetaeq_sm, [3 1 2 4]);

    pb = CmdLineProgressBar("Sorting and interpolating equivalent potential temperature to new standard grid...");
    for lo=1:size(pa,2)
        pb.print(lo, size(pa,2));
        for la=1:size(pa,3)
            for mo=1:size(pa,4)
                tmp = interp1(pa(:,lo,la,mo)./ps_vert(1,lo,la,mo), thetaeq_sm(:,lo,la,mo), 1e-5*grid.dim3.plev);

                % add surface data
                tmp(1) = thetaeqs(lo,la,mo);

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

    if any(strcmp(type, {'era5', 'era5c', 'erai'})); newdir=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s', type);
    elseif strcmp(type, 'merra2'); newdir=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s', type);
    elseif strcmp(type, 'gcm'); newdir=sprintf('/project2/tas1/miyawaki/projects/002/data/read/gcm/%s/%s', par.model, par.gcm.clim);
    elseif strcmp(type, 'echam'); newdir=sprintf('/project2/tas1/miyawaki/projects/002/data/read/echam/%s', par.echam.clim);
    elseif conthetaeqins(type, 'echam_pl'); newdir=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s', type); end;
    if ~exist(newdir, 'dir'); mkdir(newdir); end
    filename='thetaeq_si.mat';
    save(sprintf('%s/%s', newdir, filename), 'thetaeq_si', '-v7.3');
end
