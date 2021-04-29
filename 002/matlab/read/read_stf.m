function read_stf(type, ymonmean, par)

    if strcmp(ymonmean, 'ymonmean')
        ymm_in = '.ymonmean';
        ymm_out = '';
    elseif strcmp(ymonmean, 'mon')
        ymm_in = '';
        ymm_out = '_mon';
    end
    
    if strcmp(type, 'era5') | strcmp(type, 'erai') | strcmp(type, 'era5c')
        stf_vars=par.era.vars.stf;
        for i=1:length(stf_vars)
            % dimensions are (lon x lat x time)
            stf.(stf_vars{i}) = double(ncread(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/%s/stf/%s_stf_%s%s.nc', type, type, par.(type).yr_span, ymm_in), stf_vars{i}));
            % the data is originally reported as J m^-2 per day, so
            % divide by 86400 s to get the conventional W m^-2 flux
            % over the full day
            stf.(stf_vars{i}) = stf.(stf_vars{i})/86400;
        end
        save(sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s/stf%s.mat', type, par.(type).yr_span, ymm_out), 'stf', 'stf_vars');

        if contains(type, 'era5')
            for i=1:length(stf_vars)
                stf.(stf_vars{i}) = double(ncread(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/%s/stf/%s_stf_%s%s.nc', type, type, par.(type).yr_span, ymm_in), stf_vars{i}));
                stf.(stf_vars{i}) = stf.(stf_vars{i})/86400;
            end
            save(sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s/stf_2000_2012.mat', type, par.(type).yr_span), 'stf', 'stf_vars');
        end
    elseif any(strcmp(type, {'merra2', 'merra2c'}))
        stf_vars=par.(type).vars.stf;
        for i=1:length(stf_vars)
            % dimensions are (lon x lat x time)
            stf.(stf_vars{i}) = double(ncread(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/%s/stf/%s_stf_%s%s.nc', type, type, par.(type).yr_span, ymm_in), stf_vars{i}));
        end
        save(sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s/stf%s.mat', type, par.(type).yr_span, ymm_out), 'stf', 'stf_vars');

    elseif strcmp(type, 'jra55')
        stf_vars=par.(type).vars.stf;
        for i=1:length(par.jra55.vars.stf); var = par.jra55.vars.stf{i};
            if strcmp(var, 'hfls')
                stf.(stf_vars{i}) = double(ncread(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/%s/stf/%s_lhtfl_%s%s.nc', type, type, par.(type).yr_span, ymm_in), 'LHTFL_GDS0_SFC_S130'));
            elseif strcmp(var, 'hfss')
                stf.(stf_vars{i}) = double(ncread(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/%s/stf/%s_shtfl_%s%s.nc', type, type, par.(type).yr_span, ymm_in), 'SHTFL_GDS0_SFC_S130'));
            end
        end
        save(sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s/stf%s.mat', type, par.(type).yr_span, ymm_out), 'stf', 'stf_vars');

    elseif strcmp(type, 'gcm')
        stf_vars=par.gcm.vars.stf;
        for i=1:length(par.gcm.vars.stf); var = par.gcm.vars.stf{i};
            file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/gcm/%s/%s_Amon_%s_%s_r1i1p1_%s*%s.nc', par.model, var, par.model, par.(type).clim, par.(type).yr_span, ymm_in));
            fullpath=sprintf('%s/%s', file.folder, file.name);
            stf.(var)=double(ncread(fullpath, var));
        end
        newdir=sprintf('/project2/tas1/miyawaki/projects/002/data/read/gcm/%s/%s/%s', par.model, par.(type).clim, par.(type).yr_span);
        if ~exist(newdir, 'dir'); mkdir(newdir); end
        filename=sprintf('stf%s.mat', ymm_out);
        save(sprintf('%s/%s', newdir, filename), 'stf', 'stf_vars');
    elseif strcmp(type, 'hahn')
        stf_vars=par.hahn.vars.stf;
        for i=1:length(par.hahn.vars.stf); var = par.hahn.vars.stf{i};
            fprefix = make_hahn_fprefix(par);
            file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/hahn/lapserateclima/%s.%s.nc', fprefix, var));
            fullpath=sprintf('%s/%s', file.folder, file.name);
            stf.(var)=double(ncread(fullpath, 'varmo'));
        end
        save(sprintf('/project2/tas1/miyawaki/projects/002/data/read/hahn/%s/stf%s.mat', par.hahn.clim, ymm_out), 'stf', 'stf_vars');
    elseif strcmp(type, 'echam')
        stf_vars=par.echam.vars.stf;
        for i=1:length(par.echam.vars.stf); var = par.echam.vars.stf{i};
            if contains(par.echam.clim, 'echr000') | any(strcmp(par.echam.clim, par.echam.exceptions))
                file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/echam/BOT*_%s_*%s.nc', par.echam.clim, ymm_in));
            else
                file=dir(sprintf('/project2/tas1/ockham/data11/tas/echam-aiv_rcc_6.1.00p1/%s/BOT_%s_0020_39.nc', par.echam.clim, par.echam.clim));
            end
            fullpath=sprintf('%s/%s', file.folder, file.name);
            stf.(var)=double(ncread(fullpath, var));
        end
        newdir=sprintf('/project2/tas1/miyawaki/projects/002/data/read/echam/%s', par.echam.clim);
        if ~exist(newdir, 'dir'); mkdir(newdir); end
        filename=sprintf('stf%s.mat', ymm_out);
        save(sprintf('%s/%s', newdir, filename), 'stf', 'stf_vars');
    end
end
