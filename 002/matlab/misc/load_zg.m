function zg = load_zg(type, par)

    if any(strcmp(type, {'era5', 'erai', 'era5c', 'merra2', 'merra2c'}))
        file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/%s/zg/%s_zg_%s.ymonmean.nc', type, type, par.(type).yr_span));
        fullpath=sprintf('%s/%s', file.folder, file.name);
        if any(strcmp(type, {'era5', 'era5c', 'erai'}))
            zg = 1/par.g * double(ncread(fullpath, 'z')); % output is in geopotential so convert to height
        elseif strcmp(type, 'jra55')
            zg = double(ncread(fullpath, 'HGT_GDS0_ISBL_S123'));
        elseif any(strcmp(type, {'merra2', 'merra2c'}))
            zg = double(ncread(fullpath, 'H'));
        end
    elseif any(strcmp(type, {'jra55'}))
        file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/%s/zg_trop/%s_hgt_%s.ymonmean.nc', type, type, par.(type).yr_span));
        fullpath=sprintf('%s/%s', file.folder, file.name);
        zg = double(ncread(fullpath, 'HGT_GDS0_ISBL_S123'));
    elseif strcmp(type, 'gcm')
        var = 'zg';
        file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/gcm/%s/%s_Amon_%s_%s_r1i1p1_%s*.ymonmean.nc', par.model, var, par.model, par.(type).clim, par.(type).yr_span));
        fullpath=sprintf('%s/%s', file.folder, file.name);
        zg = double(ncread(fullpath, var));
    elseif strcmp(type, 'hahn')
        var = 'Z3';
        fprefix = make_hahn_fprefix(par);
        file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/hahn/lapserateclima/%s.%s.nc', fprefix, var));
        fullpath=sprintf('%s/%s', file.folder, file.name);
        zg = double(ncread(fullpath, 'varmo'));
    elseif strcmp(type, 'echam')
        var = 'geopoth';
        if contains(par.echam.clim, 'echr000') | any(strcmp(par.echam.clim, par.echam.exceptions))
            file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/echam/ATM*_%s_*.ymonmean.nc', par.echam.clim));
        else
            file=dir(sprintf('/project2/tas1/ockham/data11/tas/echam-aiv_rcc_6.1.00p1/%s/ATM_%s_0020_39.nc', par.echam.clim, par.echam.clim));
        end
        fullpath=sprintf('%s/%s', file.folder, file.name);
        zg = double(ncread(fullpath, var));
    end

end
