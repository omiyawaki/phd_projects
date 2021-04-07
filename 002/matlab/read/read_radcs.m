function read_radcs(type, ymonmean, par)

    if strcmp(ymonmean, 'ymonmean')
        ymm_in = '.ymonmean';
        ymm_out = '';
    elseif strcmp(ymonmean, 'mon')
        ymm_in = '';
        ymm_out = '_mon';
    end
    
    if strcmp(type, 'erai') | strcmp(type, 'era5') | strcmp(type, 'era5c')
        radcs_vars=par.era.vars.radcs;
        for i=1:length(radcs_vars)
            % dimensions are (lon x lat x time)
            % time is sequenced as id(1) = jan, step 00-12, id(2) = jan, step 12-24, id(3) = feb, step 00-12, etc.
            radcs.(radcs_vars{i}) = double(ncread(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/%s/radcs/%s_radcs_%s%s.nc', type, type, par.(type).yr_span, ymm_in), radcs_vars{i}));
            % the data is originally reported as J m^-2 per day, so
            % divide by 86400 s to get the conventional W m^-2 flux
            % over the full day
            radcs.(radcs_vars{i}) = radcs.(radcs_vars{i})/86400;
        end
        save(sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s/radcs%s.mat', type, par.(type).yr_span, ymm_out), 'radcs', 'radcs_vars');
        % if contains(type, 'era5')
        %     for i=1:length(radcs_vars)
        %         radcs.(radcs_vars{i}) = double(ncread(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/%s/radcs/%s_radcs_%s%s.nc', type, type, par.(type).yr_span), radcs_vars{i}));
        %         radcs.(radcs_vars{i}) = radcs.(radcs_vars{i})/86400;
        %     end
        %     save(sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s/radcs_2000_2012.mat', type, par.(type).yr_span), 'radcs', 'radcs_vars');
        % end
    elseif strcmp(type, 'merra2')
        radcs_vars=par.merra2.vars.radcs;
        for i=1:length(radcs_vars)
            % dimensions are (lon x lat x time)
            radcs.(radcs_vars{i}) = double(ncread(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/%s/radcs/%s_radcs_%s%s.nc', type, type, par.(type).yr_span, ymm_in), radcs_vars{i}));
        end
        save(sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s/radcs%s.mat', type, par.(type).yr_span, ymm_out), 'radcs', 'radcs_vars');
    elseif strcmp(type, 'jra55')
        radcs_vars=par.jra55.vars.radcs;
        for i=1:length(radcs_vars)
            % if strcmp(radcs_vars{i}, 'rsdt')
            %     radcs.(radcs_vars{i}) = double(ncread(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/%s/radcs/%s_%s_%s%s.nc', type, type, 'dswrf', par.(type).yr_span, ymm_in), 'DSWRF_GDS0_NTAT_S130'));
            if strcmp(radcs_vars{i}, 'rsdscs')
                radcs.(radcs_vars{i}) = double(ncread(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/%s/radcs/%s_%s_%s%s.nc', type, type, 'csdsf', par.(type).yr_span, ymm_in), 'CSDSF_GDS0_SFC_S130'));
            elseif strcmp(radcs_vars{i}, 'rsutcs')
                radcs.(radcs_vars{i}) = double(ncread(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/%s/radcs/%s_%s_%s%s.nc', type, type, 'csusf', par.(type).yr_span, ymm_in), 'CSUSF_GDS0_NTAT_S130'));
            elseif strcmp(radcs_vars{i}, 'rsuscs')
                radcs.(radcs_vars{i}) = double(ncread(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/%s/radcs/%s_%s_%s%s.nc', type, type, 'csusf', par.(type).yr_span, ymm_in), 'CSUSF_GDS0_SFC_S130'));
            elseif strcmp(radcs_vars{i}, 'rldscs')
                radcs.(radcs_vars{i}) = double(ncread(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/%s/radcs/%s_%s_%s%s.nc', type, type, 'csdlf', par.(type).yr_span, ymm_in), 'CSDLF_GDS0_SFC_S130'));
            elseif strcmp(radcs_vars{i}, 'rlutcs')
                radcs.(radcs_vars{i}) = double(ncread(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/%s/radcs/%s_%s_%s%s.nc', type, type, 'csulf', par.(type).yr_span, ymm_in), 'CSULF_GDS0_NTAT_S130'));
            % elseif strcmp(radcs_vars{i}, 'rlus')
                % radcs.(radcs_vars{i}) = double(ncread(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/%s/radcs/%s_%s_%s%s.nc', type, type, 'ulwrf', par.(type).yr_span, ymm_in), 'ULWRF_GDS0_SFC_S130'));
            end
        end
        save(sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s/radcs%s.mat', type, par.(type).yr_span, ymm_out), 'radcs', 'radcs_vars');

    elseif strcmp(type, 'gcm')
        radcs_vars=par.gcm.vars.radcs;
        for i=1:length(par.gcm.vars.radcs); var = par.gcm.vars.radcs{i};
            file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/gcm/%s/%s_Amon_%s_%s_r1i1p1_*%s.nc', par.model, var, par.model, par.gcm.clim, ymm_in));
            fullpath=sprintf('%s/%s', file.folder, file.name);
            radcs.(var)=double(ncread(fullpath, var));
        end
        newdir=sprintf('/project2/tas1/miyawaki/projects/002/data/read/gcm/%s/%s', par.model, par.gcm.clim);
        if ~exist(newdir, 'dir'); mkdir(newdir); end
        filename=sprintf('radcs%s.mat', ymm_out);
        save(sprintf('%s/%s', newdir, filename), 'radcs', 'radcs_vars');
    elseif strcmp(type, 'echam')
        radcs_vars=par.echam.vars.radcs;
        for i=1:length(par.echam.vars.radcs); var = par.echam.vars.radcs{i};
            if contains(par.echam.clim, 'echr0') | any(strcmp(par.echam.clim, par.echam.exceptions))
                file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/echam/BOT*_%s_*.ymonmean.nc', par.echam.clim));
            else 
                file=dir(sprintf('/project2/tas1/ockham/data11/tas/echam-aiv_rcc_6.1.00p1/%s/BOT_%s_0020_39.nc', par.echam.clim, par.echam.clim));
            end
            fullpath=sprintf('%s/%s', file.folder, file.name);
            radcs.(var)=double(ncread(fullpath, var));
        end
        save(sprintf('/project2/tas1/miyawaki/projects/002/data/read/echam/%s/radcs%s.mat', par.echam.clim, ymm_out), 'radcs', 'radcs_vars');
    elseif strcmp(type, 'ceres')
        radcs_vars = par.ceres.vars.radcs;
        radcs_vars_txt = par.ceres.vars.radcs_txt;
        for i=1:length(radcs_vars)
            radcs.(radcs_vars_txt{i}) = double(ncread(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/%s/CERES_EBAF_Ed4.1_Subset_%s%s.nc', type, par.(type).yr_span, ymm_in), radcs_vars{i}));
        end
        save(sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/radcs%s.mat', type, ymm_out), 'radcs');

        for i=1:length(radcs_vars)
            radcs.(radcs_vars_txt{i}) = double(ncread(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/%s/CERES_EBAF_Ed4.1_Subset_200101-200912%s.nc', type), radcs_vars{i}, ymm_in));
        end
        save(sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/radcs_2001_2009.mat', type), 'radcs');
    end
end
