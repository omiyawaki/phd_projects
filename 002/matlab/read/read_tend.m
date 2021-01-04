function read_tend(type, par)
    if strcmp(type, 'era5') | strcmp(type, 'erai') | strcmp(type, 'era5c')
        tend_vars=par.era.vars.tend;
        tend_vars_txt=par.era.vars.tend_txt;
        for i=1:length(tend_vars)
            % dimensions are (lon x lat x time)
            tend.(tend_vars_txt{i}) = double(ncread(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/%s/tend/%s_tend_%s.deltat.ymonmean.nc', type, type, par.(type).yr_span), tend_vars{i}));
            % the data is originally reported as J m^-2, so
            % divide by 6 hr (because data is 6 hourly) to
            % convert to W m^-2.
            tend.(tend_vars_txt{i}) = tend.(tend_vars_txt{i})/(6*3600);
        end
        save(sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s/tend.mat', type, par.(type).yr_span), 'tend', 'tend_vars_txt');

    else
        error('MSE tendency data are read only for ERA data.');
    end
end
