function read_lfrac(type, par)
% land fraction
    if any(strcmp(type, {'era5', 'era5c', 'erai'}))
        sftlf = double(ncread(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/%s/lfrac/%s_lfrac_%s.ymonmean.nc', type, type, par.(type).yr_span), 'lsm'));
        save(sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s/sftlf.mat', type, par.(type).yr_span), 'sftlf');
    elseif strcmp(type, 'gcm')
        if any(strcmp(par.gcm.clim, {'piControl', 'abrupt4xCO2'}))
            file=dir(sprintf('/project2/tas1/CMIP5_piControl/%s/sftlf_*.nc', par.model));
        elseif any(strcmp(par.gcm.clim, {'historical', 'rcp85'}))
            if any(strcmp(par.model, {'bcc-csm1-1-m', 'GISS-E2-H-CC', 'GISS-E2-R-CC'}))
                file=dir(sprintf('/project2/tas1/ockham/data14/tas/CMIP5_RAW/%s/piControl/atmos/fx/sftlf/r0i0p0/sftlf_*.nc', par.model));
            elseif strcmp(par.model, 'inmcm4')
                file=dir(sprintf('/project2/tas1/ockham/data9/tas/CMIP5_RAW/%s/rcp85/atmos/fx/sftlf/r0i0p0/sftlf_*.nc', par.model));
            else
                file=dir(sprintf('/project2/tas1/ockham/data9/tas/CMIP5_RAW/%s/historical/atmos/fx/sftlf/r0i0p0/sftlf_*.nc', par.model));
            end
        else
            file=dir(sprintf('/project2/tas1/ockham/data9/tas/CMIP5_RAW/%s/%s/atmos/fx/sftlf/r0i0p0/sftlf_*.nc', par.model, par.(type).clim));
        end
        fullpath=sprintf('%s/%s', file.folder, file.name);
        sftlf=double(ncread(fullpath, 'sftlf'));
        newdir=sprintf('/project2/tas1/miyawaki/projects/002/data/read/gcm/%s/%s/%s', par.model, par.(type).clim, par.(type).yr_span);
        if ~exist(newdir, 'dir'); mkdir(newdir); end
        filename='sftlf.mat';
        save(sprintf('%s/%s', newdir, filename), 'sftlf');
    elseif strcmp(type, 'merra2')
        sftlf = double(ncread(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/%s/lfrac/MERRA2_101.const_2d_asm_Nx.00000000.nc4.nc4', type), 'FRLAND'));
        save(sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s/sftlf.mat', type, par.(type).yr_span), 'sftlf');
    elseif strcmp(type, 'merra2c')
        sftlf = double(ncread(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/%s/lfrac/merra2c_lfrac.nc', type), 'FRLAND'));
        save(sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s/sftlf.mat', type, par.(type).yr_span), 'sftlf');
    elseif strcmp(type, 'jra55')
        sftlf = double(ncread(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/%s/lfrac/%s_land_%s.ymonmean.nc', type, type, par.(type).yr_span), 'LAND_GDS0_SFC'));
        save(sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s/sftlf.mat', type, par.(type).yr_span), 'sftlf');
    elseif strcmp(type, 'echam')
        file=dir(sprintf('/project2/tas1/miyawaki/models/echam/inputdata/ECHAM6/jsbach/T63/jsbach_T63GR15_11tiles_1976.nc'));
        fullpath=sprintf('%s/%s', file.folder, file.name);
        sftlf=double(ncread(fullpath, 'slm'));
        newdir=sprintf('/project2/tas1/miyawaki/projects/002/data/read/echam/%s', par.echam.clim);
        if ~exist(newdir, 'dir'); mkdir(newdir); end
        filename='sftlf.mat';
        save(sprintf('%s/%s', newdir, filename), 'sftlf');
    end
end
