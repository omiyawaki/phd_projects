function read_rad(type, ymonmean, par)

    if strcmp(ymonmean, 'ymonmean')
        ymm_in = '.ymonmean';
        ymm_out = '';
    elseif strcmp(ymonmean, 'mon')
        ymm_in = '';
        ymm_out = '_mon';
    end
    
    if strcmp(type, 'erai') | strcmp(type, 'era5') | strcmp(type, 'era5c')
        rad_vars=par.era.vars.rad;
        for i=1:length(rad_vars)
            % dimensions are (lon x lat x time)
            % time is sequenced as id(1) = jan, step 00-12, id(2) = jan, step 12-24, id(3) = feb, step 00-12, etc.
            rad.(rad_vars{i}) = double(ncread(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/%s/rad/%s_rad_%s%s.nc', type, type, par.(type).yr_span, ymm_in), rad_vars{i}));
            % the data is originally reported as J m^-2 per day, so
            % divide by 86400 s to get the conventional W m^-2 flux
            % over the full day
            rad.(rad_vars{i}) = rad.(rad_vars{i})/86400;
        end
        save(sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s/rad%s.mat', type, par.(type).yr_span, ymm_out), 'rad', 'rad_vars');
        % if contains(type, 'era5')
        %     for i=1:length(rad_vars)
        %         rad.(rad_vars{i}) = double(ncread(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/%s/rad/%s_rad_%s%s.nc', type, type, par.(type).yr_span), rad_vars{i}));
        %         rad.(rad_vars{i}) = rad.(rad_vars{i})/86400;
        %     end
        %     save(sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s/rad_2000_2012.mat', type, par.(type).yr_span), 'rad', 'rad_vars');
        % end
    elseif any(strcmp(type, {'merra2', 'merra2c'}))
        rad_vars=par.(type).vars.rad;
        for i=1:length(rad_vars)
            % dimensions are (lon x lat x time)
            rad.(rad_vars{i}) = double(ncread(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/%s/rad/%s_rad_%s%s.nc', type, type, par.(type).yr_span, ymm_in), rad_vars{i}));
        end
        save(sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s/rad%s.mat', type, par.(type).yr_span, ymm_out), 'rad', 'rad_vars');
    elseif strcmp(type, 'jra55')
        rad_vars=par.jra55.vars.rad;
        for i=1:length(rad_vars)
            if strcmp(rad_vars{i}, 'rsdt')
                rad.(rad_vars{i}) = double(ncread(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/%s/rad/%s_%s_%s%s.nc', type, type, 'dswrf', par.(type).yr_span, ymm_in), 'DSWRF_GDS0_NTAT_S130'));
            elseif strcmp(rad_vars{i}, 'rsds')
                rad.(rad_vars{i}) = double(ncread(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/%s/rad/%s_%s_%s%s.nc', type, type, 'dswrf', par.(type).yr_span, ymm_in), 'DSWRF_GDS0_SFC_S130'));
            elseif strcmp(rad_vars{i}, 'rsut')
                rad.(rad_vars{i}) = double(ncread(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/%s/rad/%s_%s_%s%s.nc', type, type, 'uswrf', par.(type).yr_span, ymm_in), 'USWRF_GDS0_NTAT_S130'));
            elseif strcmp(rad_vars{i}, 'rsus')
                rad.(rad_vars{i}) = double(ncread(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/%s/rad/%s_%s_%s%s.nc', type, type, 'uswrf', par.(type).yr_span, ymm_in), 'USWRF_GDS0_SFC_S130'));
            elseif strcmp(rad_vars{i}, 'rlds')
                rad.(rad_vars{i}) = double(ncread(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/%s/rad/%s_%s_%s%s.nc', type, type, 'dlwrf', par.(type).yr_span, ymm_in), 'DLWRF_GDS0_SFC_S130'));
            elseif strcmp(rad_vars{i}, 'rlut')
                rad.(rad_vars{i}) = double(ncread(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/%s/rad/%s_%s_%s%s.nc', type, type, 'ulwrf', par.(type).yr_span, ymm_in), 'ULWRF_GDS0_NTAT_S130'));
            elseif strcmp(rad_vars{i}, 'rlus')
                rad.(rad_vars{i}) = double(ncread(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/%s/rad/%s_%s_%s%s.nc', type, type, 'ulwrf', par.(type).yr_span, ymm_in), 'ULWRF_GDS0_SFC_S130'));
            end
        end
        save(sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s/rad%s.mat', type, par.(type).yr_span, ymm_out), 'rad', 'rad_vars');

    elseif strcmp(type, 'gcm')
        rad_vars=par.gcm.vars.rad;
        for i=1:length(par.gcm.vars.rad); var = par.gcm.vars.rad{i};
            file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/gcm/%s/%s_Amon_%s_%s_r1i1p1_%s*%s.nc', par.model, var, par.model, par.(type).clim, par.(type).yr_span, ymm_in));
            fullpath=sprintf('%s/%s', file.folder, file.name);
            rad.(var)=double(ncread(fullpath, var));
        end
        newdir=sprintf('/project2/tas1/miyawaki/projects/002/data/read/gcm/%s/%s/%s', par.model, par.(type).clim, par.(type).yr_span);
        if ~exist(newdir, 'dir'); mkdir(newdir); end
        filename=sprintf('rad%s.mat', ymm_out);
        save(sprintf('%s/%s', newdir, filename), 'rad', 'rad_vars');
    elseif strcmp(type, 'hahn')
        rad_vars=par.hahn.vars.rad;
        for i=1:length(par.hahn.vars.rad); var = par.hahn.vars.rad{i};
            fprefix = make_hahn_fprefix(par);
            file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/hahn/lapserateclima/%s.%s.nc', fprefix, var));
            fullpath=sprintf('%s/%s', file.folder, file.name);
            rad.(var)=double(ncread(fullpath, 'varmo'));
        end
        save(sprintf('/project2/tas1/miyawaki/projects/002/data/read/hahn/%s/rad%s.mat', par.hahn.clim, ymm_out), 'rad', 'rad_vars');
    elseif strcmp(type, 'echam')
        rad_vars=par.echam.vars.rad;
        for i=1:length(par.echam.vars.rad); var = par.echam.vars.rad{i};
            if contains(par.echam.clim, 'echr000') | any(strcmp(par.echam.clim, par.echam.exceptions))
                file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/echam/BOT*_%s_*%s.nc', par.echam.clim, ymm_in));
            else 
                file=dir(sprintf('/project2/tas1/ockham/data11/tas/echam-aiv_rcc_6.1.00p1/%s/BOT_%s_0020_39.nc', par.echam.clim, par.echam.clim));
            end
            fullpath=sprintf('%s/%s', file.folder, file.name);
            rad.(var)=double(ncread(fullpath, var));
        end
        save(sprintf('/project2/tas1/miyawaki/projects/002/data/read/echam/%s/rad%s.mat', par.echam.clim, ymm_out), 'rad', 'rad_vars');
    elseif strcmp(type, 'ceres')
        rad_vars = par.ceres.vars.rad;
        rad_vars_txt = par.ceres.vars.rad_txt;
        for i=1:length(rad_vars)
            rad.(rad_vars_txt{i}) = double(ncread(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/%s/CERES_EBAF_Ed4.1_Subset_%s%s.nc', type, par.(type).yr_span, ymm_in), rad_vars{i}));
        end
        save(sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/rad%s.mat', type, ymm_out), 'rad');

        for i=1:length(rad_vars)
            rad.(rad_vars_txt{i}) = double(ncread(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/%s/CERES_EBAF_Ed4.1_Subset_200101-200912%s.nc', type), rad_vars{i}, ymm_in));
        end
        save(sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/rad_2001_2009.mat', type), 'rad');
    end
end
