function [mse_sm, sm, pa, ps_vert, mses, srfc, grid] = load_mse(type, par)

    prefix = make_prefix(type, par);
    savedir = make_savedir(type, par);
    ta_orig = load_temp(type, par);
    hur_orig = load_rh(type, par);
    zg_orig = load_zg(type, par);

    load(sprintf('%s/grid.mat', prefix)); % read grid data
    load(sprintf('%s/srfc.mat', prefix)); % load surface data

    esat = calc_esat(ta_orig, 0); % compute saturation vapor pressure
    p = permute(repmat(grid.dim3.plev, [1 size(ta_orig,1) size(ta_orig,2) size(ta_orig,4)]), [2 3 1 4]);
    e = esat.*hur_orig/100;
    r = calc_r(p, e, par);
    q = calc_q(p, e, par);
    clear p e esat hur_orig;

    % compute eq potential temperature following AMS glossary definition
    mse = par.cpd*ta_orig + par.L*q + par.g*zg_orig;
    clear ta_orig zg_orig;

    % surface eq pot temp
    if strcmp(type, 'era5') | strcmp(type, 'erai') | strcmp(type, 'era5c')
        tas = srfc.t2m; % surface is the last element in era plev
        ps = srfc.sp;
        e2m = calc_esat(srfc.d2m, par.frz);
        esat2m = calc_esat(tas, par.frz);
        hurs = 100 * e2m./esat2m;
        clear e2m esat2m;
    elseif any(strcmp(type, {'merra2', 'merra2c'}))
        tas = srfc.T2M;
        ps = srfc.PS;
        q2m = srfc.QV2M;
        e2m = calc_e(ps, q2m, par);
        esat2m = calc_esat(tas, par.frz);
        hurs = 100 * e2m./esat2m;
        clear q2m e2m esat2m
    elseif any(strcmp(type, {'gcm', 'jra55'}))
        tas = srfc.tas;
        ps = srfc.ps;
        hurs = srfc.hurs;
    elseif contains(type, 'echam')
        tas = srfc.temp2;
    end

    esats = calc_esat(tas, 0); % compute saturation vapor pressure
    es = esats.*hurs/100;
    rs = calc_r(ps, es, par);
    qs = calc_q(ps, es, par);
    clear ps es esats;

    mses = par.cpd*tas + par.L*qs + par.g*srfc.zs;
    clear tas hurs;

    % create surface mask
    if strcmp(type, 'era5') | strcmp(type, 'erai') | strcmp(type, 'era5c')
        ps_vert = repmat(srfc.sp, [1 1 1 size(mse, 3)]); % dims (lon x lat x time x plev)
        ps_vert = permute(ps_vert, [1 2 4 3]); % dims (lon x lat x plev x time)
        pa = double(permute(repmat(grid.dim3.plev, [1 size(srfc.sp)]), [2 3 1 4]));
    elseif any(strcmp(type, {'merra2', 'merra2c'}))
        ps_vert = repmat(srfc.PS, [1 1 1 size(mse, 3)]); % dims (lon x lat x time x plev)
        ps_vert = permute(ps_vert, [1 2 4 3]); % dims (lon x lat x plev x time)
        pa = double(permute(repmat(grid.dim3.plev, [1 size(srfc.PS)]), [2 3 1 4]));
    elseif any(strcmp(type, {'gcm', 'jra55'}))
        ps_vert = repmat(srfc.ps, [1 1 1 size(mse, 3)]); % dims (lon x lat x time x plev)
        ps_vert = permute(ps_vert, [1 2 4 3]); % dims (lon x lat x plev x time)
        pa = permute(repmat(grid.dim3.plev, [1 size(srfc.ps)]), [2 3 1 4]);
    elseif strcmp(type, 'echam')
        ps_vert = repmat(srfc.aps, [1 1 1 size(mse, 3)]); % dims (lon x lat x time x plev)
        ps_vert = permute(ps_vert, [1 2 4 3]); % dims (lon x lat x plev x time)
        pa = permute(repmat(grid.dim3.plev, [1 size(srfc.aps)]), [2 3 1 4]);
    elseif contains(type, 'echam_pl')
        ps_vert = repmat(srfc.aps, [1 1 1 size(mse, 3)]); % dims (lon x lat x time x plev)
        ps_vert = permute(ps_vert, [1 2 4 3]); % dims (lon x lat x plev x time)
        pa = permute(repmat(grid.dim3.plev, [1 size(srfc.aps)]), [2 3 1 4]);
    end
    sm = nan(size(mse));
    sm(pa < 0.9961*ps_vert) = 1;
    mse_sm = mse.*sm; % filter mse with surface mask

    ps_vert = permute(ps_vert, [3 1 2 4]);
    pa = permute(pa, [3 1 2 4]);
    mse_sm = permute(mse_sm, [3 1 2 4]);

end