function [temp_ml, lon, lat] = load_tempml(type, par)

    if any(strcmp(type, {'era5', 'erai', 'era5c', 'merra2', 'merra2c'}))
        file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/%s/temp_ml/%s_temp_ml_%s.ymonmean.nc', type, type, par.(type).yr_span));
        fullpath=sprintf('%s/%s', file.folder, file.name);
        if any(strcmp(type, {'era5', 'era5c', 'erai'}))
            temp_ml = double(ncread(fullpath, 't'));
        elseif strcmp(type, 'merra2')
            temp_ml = double(ncread(fullpath, 'T'));
        end
        if strcmp(type, 'erai')
            lon = double(ncread(fullpath, 'longitude'));
            lat = double(ncread(fullpath, 'latitude'));
        else
            lon = double(ncread(fullpath, 'lon'));
            lat = double(ncread(fullpath, 'lat'));
        end
    elseif any(strcmp(type, {'jra55'}))
        file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/%s/temp_ml/%s_tmp_%s.ymonmean.nc', type, type, par.(type).yr_span));
        fullpath=sprintf('%s/%s', file.folder, file.name);
        temp_ml = double(ncread(fullpath, 'TMP_GDS4_HYBL_S123'));
        lon = double(ncread(fullpath, 'g4_lon_3'));
        lat = double(ncread(fullpath, 'g4_lat_2'));
    elseif strcmp(type, 'gcm')
        var = 'ta';
        file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/gcm/%s/%s_Amon_%s_%s_r1i1p1_%s*.ymonmean.nc', par.model, var, par.model, par.(type).clim, par.(type).yr_span));
        fullpath=sprintf('%s/%s', file.folder, file.name);
        temp_ml = double(ncread(fullpath, var));
    elseif strcmp(type, 'echam')
        var = 't';
        if contains(par.echam.clim, 'rp000')
            file=dir(sprintf('/project2/tas1/ockham/data11/tas/echam-aiv_rcc_6.1.00p1/%s/ATM_%s_0020_39.nc', par.echam.clim, par.echam.clim));
        else
            file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/echam/ATM*_%s_*.ymonmean.nc', par.echam.clim));
        end
        fullpath=sprintf('%s/%s', file.folder, file.name);
        temp_ml = double(ncread(fullpath, var));
    elseif strcmp(type, 'hahn')
        var = 'T';
        fprefix = make_hahn_fprefix(par);
        file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/hahn/lapserateclima/%s.%s.nc', fprefix, var));
        fullpath=sprintf('%s/%s', file.folder, file.name);
        temp_ml = double(ncread(fullpath, 'varmo'));
    end

end
