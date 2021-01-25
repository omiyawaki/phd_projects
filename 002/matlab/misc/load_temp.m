function temp = load_temp(type, par)

    if any(strcmp(type, {'era5', 'erai', 'era5c', 'merra2'}))
        file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/%s/temp/%s_temp_%s.ymonmean.nc', type, type, par.(type).yr_span));
        fullpath=sprintf('%s/%s', file.folder, file.name);
        if any(strcmp(type, {'era5', 'era5c', 'erai'}))
            temp = double(ncread(fullpath, 't'));
        elseif strcmp(type, 'merra2')
            temp = double(ncread(fullpath, 'T'));
        end
    elseif any(strcmp(type, {'jra55'}))
        file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/%s/temp/%s_tmp_%s.ymonmean.nc', type, type, par.(type).yr_span));
        fullpath=sprintf('%s/%s', file.folder, file.name);
        temp = double(ncread(fullpath, 'TMP_GDS0_ISBL_S123'));
    elseif strcmp(type, 'gcm')
        var = 'ta';
        file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/gcm/%s/%s_Amon_%s_%s_r1i1p1_*.ymonmean.nc', par.model, var, par.model, par.gcm.clim));
        fullpath=sprintf('%s/%s', file.folder, file.name);
        temp = double(ncread(fullpath, var));
    elseif strcmp(type, 'echam')
        var = 't';
        if contains(par.echam.clim, 'rp000')
            file=dir(sprintf('/project2/tas1/ockham/data11/tas/echam-aiv_rcc_6.1.00p1/%s/ATM_%s_0020_39.nc', par.echam.clim, par.echam.clim));
        else
            file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/echam/ATM*_%s_*.ymonmean.nc', par.echam.clim));
        end
        fullpath=sprintf('%s/%s', file.folder, file.name);
        temp = double(ncread(fullpath, var));
    elseif strcmp(type, 'hahn')
        var = 'T';
        fprefix = make_hahn_fprefix(par);
        file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/hahn/lapserateclima/%s.%s.nc', fprefix, var));
        fullpath=sprintf('%s/%s', file.folder, file.name);
        temp = double(ncread(fullpath, 'varmo'));
    end

end
