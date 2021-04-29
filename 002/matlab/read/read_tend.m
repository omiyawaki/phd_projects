function read_tend(type, ymonmean, par)

    if strcmp(ymonmean, 'ymonmean')
        ymm_in = '.ymonmean';
        ymm_out = '';
    elseif strcmp(ymonmean, 'mon')
        ymm_in = '';
        ymm_out = '_mon';
    end

    if strcmp(type, 'era5') | strcmp(type, 'erai') | strcmp(type, 'era5c')
        tend_vars=par.era.vars.tend;
        tend_vars_txt=par.era.vars.tend_txt;
        for i=1:length(tend_vars)
            % dimensions are (lon x lat x time)
            if strcmp(tend_vars_txt{i}, 'tend');
                tend.(tend_vars_txt{i}) = double(ncread(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/%s/tend/%s_tend_%s.ymonmean.nc', type, type, par.(type).yr_span), 'tend'));
            elseif strcmp(tend_vars_txt{i}, 'tendmon')
                tend.(tend_vars_txt{i}) = double(ncread(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/%s/tendmon/%s_tendmon_%s.ymonmean.nc', type, type, par.(type).yr_span), 'tend'));
            elseif strcmp(tend_vars_txt{i}, 'tendalt')
                tend.(tend_vars_txt{i}) = double(ncread(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/%s/tendalt/%s_tend_%s.deltat.ymonmean.nc', type, type, par.(type).yr_span), 'p62.162'));
                % the data is originally reported as J m^-2, so
                % divide by 6 hr (because data is 6 hourly) to
                % convert to W m^-2.
                tend.(tend_vars_txt{i}) = tend.(tend_vars_txt{i})/(6*3600);
            end
        end
        save(sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s/tend%s.mat', type, par.(type).yr_span, ymm_out), 'tend');
    elseif any(strcmp(type, {'merra2', 'merra2c', 'jra55'}))
        tend_vars=par.(type).vars.tend;
        tend_vars_txt=par.(type).vars.tend_txt;
        for i=1:length(tend_vars)
            % dimensions are (lon x lat x time)
            if strcmp(tend_vars_txt{i}, 'tend');
                tend.(tend_vars_txt{i}) = double(ncread(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/%s/tend/%s_tend_%s.ymonmean.nc', type, type, par.(type).yr_span), 'tend'));
            elseif strcmp(tend_vars_txt{i}, 'tendmon')
                tend.(tend_vars_txt{i}) = double(ncread(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/%s/tendmon/%s_tendmon_%s.ymonmean.nc', type, type, par.(type).yr_span), 'tend'));
            end
        end
        save(sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s/tend%s.mat', type, par.(type).yr_span, ymm_out), 'tend');
    elseif strcmp(type, 'gcm')
        tend_vars=par.gcm.vars.tend;
        tend_vars_txt=par.gcm.vars.tend_txt;
        for i=1:length(par.gcm.vars.tend); var = par.gcm.vars.tend{i}; var_txt = par.gcm.vars.tend_txt{i};
            if strcmp(var, 'tend');
                file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/gcm/%s/%s_day_%s_%s_r1i1p1_%s*%s.nc', par.model, var, par.model, par.gcm.clim, par.gcm.yr_span, ymm_in));
            elseif strcmp(var, 'tendmon')
                file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/gcm/%s/%s_Amon_%s_%s_r1i1p1_%s*%s.nc', par.model, var, par.model, par.gcm.clim, par.gcm.yr_span, ymm_in));
            end
            fullpath=sprintf('%s/%s', file.folder, file.name);
            tend.(var)=double(ncread(fullpath, 'tend'));
        end
        newdir=sprintf('/project2/tas1/miyawaki/projects/002/data/read/gcm/%s/%s/%s', par.model, par.(type).clim, par.(type).yr_span);
        if ~exist(newdir, 'dir'); mkdir(newdir); end
        filename=sprintf('tend%s.mat', ymm_out);
        save(sprintf('%s/%s', newdir, filename), 'tend', 'tend_vars');
    elseif strcmp(type, 'echam')
        tend_vars=par.echam.vars.tend;
        tend_vars_txt=par.echam.vars.tend_txt;
        for i=1:length(par.echam.vars.tend); var = par.echam.vars.tend{i}; var_txt = par.echam.vars.tend_txt{i};
            if strcmp(par.echam.clim, 'rp000172')
                file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/echam/%s/%s_%s_0009_14*%s.nc', par.echam.clim, var, par.echam.clim, ymm_in));
            else
                file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/echam/%s/%s_%s_0020_39*%s.nc', par.echam.clim, var, par.echam.clim, ymm_in));
            end
            fullpath=sprintf('%s/%s', file.folder, file.name);
            tend.(var)=double(ncread(fullpath, 'tend'));
        end
        save(sprintf('/project2/tas1/miyawaki/projects/002/data/read/echam/%s/tend%s.mat', par.echam.clim, ymm_out), 'tend');
    else
        error('MSE tendency data are not available for this data.');
    end
end
