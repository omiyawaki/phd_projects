function read_hydro(type, ymonmean, par)

    if strcmp(ymonmean, 'ymonmean')
        ymm_in = '.ymonmean';
        ymm_out = '';
    elseif strcmp(ymonmean, 'mon')
        ymm_in = '';
        ymm_out = '_mon';
    end
    
    if strcmp(type, 'era5') | strcmp(type, 'erai') | strcmp(type, 'era5c')
        hydro_vars=par.era.vars.hydro;
        for i=1:length(hydro_vars)
            % dimensions are (lon x lat x time)
            hydro.(hydro_vars{i}) = double(ncread(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/%s/hydro/%s_hydro_%s%s.nc', type, type, par.(type).yr_span, ymm_in), hydro_vars{i}));
            % the data is originally reported as accumulated m (depth) over a day, so
            % divide by 86400 s and multiply by 1000 kg/m^3 to get the
            % conventional kg/m^2/s mass flux over the full day
            hydro.(hydro_vars{i}) = hydro.(hydro_vars{i})/86400*1e3;
        end
        save(sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s/hydro%s.mat', type, par.(type).yr_span, ymm_out), 'hydro', 'hydro_vars');

    elseif strcmp(type, 'merra2')
        hydro_vars=par.merra2.vars.hydro;
        for i=1:length(hydro_vars)
            % dimensions are (lon x lat x time)
            hydro.(hydro_vars{i}) = double(ncread(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/%s/hydro/%s_hydro_%s%s.nc', type, type, par.(type).yr_span, ymm_in), hydro_vars{i}));
        end
        save(sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s/hydro%s.mat', type, par.(type).yr_span, ymm_out), 'hydro', 'hydro_vars');

    elseif strcmp(type, 'jra55')
        hydro_vars=par.jra55.vars.hydro;
        for i=1:length(hydro_vars)
            if strcmp(hydro_vars{i}, 'pr')
                % dimensions are (lon x lat x time)
                hydro.(hydro_vars{i}) = double(ncread(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/%s/hydro/%s_tprat_%s%s.nc', type, type, par.(type).yr_span, ymm_in), 'TPRAT_GDS0_SFC_S130'));
                % the data is originally reported as mm per day, so
                % divide by 86400 (to convert day to s) to get the
                % conventional kg/m^2/s mass flux over the full day
                hydro.(hydro_vars{i}) = hydro.(hydro_vars{i})/86400;
            elseif strcmp(hydro_vars{i}, 'prc')
                hydro.(hydro_vars{i}) = double(ncread(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/%s/hydro/%s_cprat_%s%s.nc', type, type, par.(type).yr_span, ymm_in), 'CPRAT_GDS0_SFC_S130'));
                hydro.(hydro_vars{i}) = hydro.(hydro_vars{i})/86400;
            elseif strcmp(hydro_vars{i}, 'evspsbl')
                hydro.(hydro_vars{i}) = double(ncread(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/%s/hydro/%s_evp_%s%s.nc', type, type, par.(type).yr_span, ymm_in), 'EVP_GDS0_SFC_S130'));
                hydro.(hydro_vars{i}) = hydro.(hydro_vars{i})/86400;
            end
        end
        save(sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s/hydro%s.mat', type, par.(type).yr_span, ymm_out), 'hydro', 'hydro_vars');

    elseif strcmp(type, 'gcm')
        hydro_vars=par.gcm.vars.hydro;
        for i=1:length(par.gcm.vars.hydro); var = par.gcm.vars.hydro{i};
            file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/gcm/%s/%s_Amon_%s_%s_r1i1p1_*%s.nc', par.model, var, par.model, par.gcm.clim, ymm_in));
            fullpath=sprintf('%s/%s', file.folder, file.name);
            hydro.(var)=double(ncread(fullpath, var));
        end
        newdir=sprintf('/project2/tas1/miyawaki/projects/002/data/read/gcm/%s/%s', par.model, par.gcm.clim);
        if ~exist(newdir, 'dir'); mkdir(newdir); end
        filename=sprintf('hydro%s.mat', ymm_out);
        save(sprintf('%s/%s', newdir, filename), 'hydro', 'hydro_vars');
    elseif strcmp(type, 'hahn')
        hydro_vars=par.hahn.vars.hydro;
        for i=1:length(par.hahn.vars.hydro); var = par.hahn.vars.hydro{i};
            fprefix = make_hahn_fprefix(par);
            file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/hahn/lapserateclima/%s.%s.nc', fprefix, var));
            fullpath=sprintf('%s/%s', file.folder, file.name);
            hydro.(var)=double(ncread(fullpath, 'varmo'));
        end
        save(sprintf('/project2/tas1/miyawaki/projects/002/data/read/hahn/%s/hydro%s.mat', par.hahn.clim, ymm_out), 'hydro', 'hydro_vars');
    elseif strcmp(type, 'echam')
        hydro_vars=par.echam.vars.hydro;
        for i=1:length(par.echam.vars.hydro); var = par.echam.vars.hydro{i};
            if contains(par.echam.clim, 'echr000') | strcmp(par.echam.clim, 'rp000092')
                file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/echam/BOT*_%s_*%s.nc', par.echam.clim, ymm_in));
            else
                file=dir(sprintf('/project2/tas1/ockham/data11/tas/echam-aiv_rcc_6.1.00p1/%s/BOT_%s_0020_39.nc', par.echam.clim, par.echam.clim));
            end
            fullpath=sprintf('%s/%s', file.folder, file.name);
            hydro.(var)=double(ncread(fullpath, var));
        end
        newdir=sprintf('/project2/tas1/miyawaki/projects/002/data/read/echam/%s', par.echam.clim);
        if ~exist(newdir, 'dir'); mkdir(newdir); end
        filename=sprintf('hydro%s.mat', ymm_out);
        save(sprintf('%s/%s', newdir, filename), 'hydro', 'hydro_vars');
    end
end
