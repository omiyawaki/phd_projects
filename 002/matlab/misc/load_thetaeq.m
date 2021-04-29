function [thetaeq_sm, sm, pa, ps_vert, thetaeqs, srfc, grid] = load_thetaeq(type, par)

    prefix = make_prefix(type, par);
    savedir = make_savedir(type, par);
    ta_orig = load_temp(type, par);
    hur_orig = load_rh(type, par);

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
        e2m = calc_esat(srfc.d2m, par.frz);
        esat2m = calc_esat(tas, par.frz);
        hurs = 100 * e2m./esat2m;
        clear e2m esat2m;
    elseif strcmp(type, 'merra2')
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
    elseif any(strcmp(type, {'gcm', 'jra55'}))
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

end