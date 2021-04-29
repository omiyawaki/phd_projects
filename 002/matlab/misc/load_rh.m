function rh = load_rh(type, par)

    if any(strcmp(type, {'era5', 'erai', 'era5c'}))
        file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/%s/rh/%s_rh_%s.ymonmean.nc', type, type, par.(type).yr_span));
        fullpath=sprintf('%s/%s', file.folder, file.name);
        rh = double(ncread(fullpath, 'r'));
    elseif any(strcmp(type, {'merra2', 'merra2c'}))
        file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/%s/hus/%s_hus_%s.ymonmean.nc', type, type, par.(type).yr_span));
        fullpath=sprintf('%s/%s', file.folder, file.name);
        hus = double(ncread(fullpath, 'QV'));
        p = double(ncread(fullpath, 'lev'));
        p = repmat(p, [1 size(hus, 1) size(hus, 2) size(hus, 4)]);
        p = permute(p, [2 3 1 4]);
        e = calc_e(p, hus, par);

        file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/%s/temp/%s_temp_%s.ymonmean.nc', type, type, par.(type).yr_span));
        fullpath=sprintf('%s/%s', file.folder, file.name);
        temp = double(ncread(fullpath, 'T'));
        esat = calc_esat(temp, par.frz);

        rh = 100 * e./esat;
    elseif any(strcmp(type, {'jra55'}))
        file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/%s/rh/%s_rh_%s.ymonmean.nc', type, type, par.(type).yr_span));
        fullpath=sprintf('%s/%s', file.folder, file.name);
        rh = double(ncread(fullpath, 'RH_GDS0_ISBL_S123'));
    elseif strcmp(type, 'gcm')
        var = 'hur';
        file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/gcm/%s/%s_Amon_%s_%s_r1i1p1_%s*.ymonmean.nc', par.model, var, par.model, par.(type).clim, par.(type).yr_span));
        fullpath=sprintf('%s/%s', file.folder, file.name);
        rh = double(ncread(fullpath, var));
    elseif strcmp(type, 'echam')
        var = 't';
        if contains(par.echam.clim, 'echr000') | any(strcmp(par.echam.clim, par.echam.exceptions))
            file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/echam/ATM*_%s_*.ymonmean.nc', par.echam.clim));
        else
            file=dir(sprintf('/project2/tas1/ockham/data11/tas/echam-aiv_rcc_6.1.00p1/%s/ATM_%s_0020_39.nc', par.echam.clim, par.echam.clim));
        end
        fullpath=sprintf('%s/%s', file.folder, file.name);
        rh = double(ncread(fullpath, var));
    elseif strcmp(type, 'hahn')
        var = 'T';
        fprefix = make_hahn_fprefix(par);
        file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/hahn/lapserateclima/%s.%s.nc', fprefix, var));
        fullpath=sprintf('%s/%s', file.folder, file.name);
        rh = double(ncread(fullpath, 'varmo'));
    end

end
