function read_srfc(type, ymonmean, par)

    if strcmp(ymonmean, 'ymonmean')
        ymm_in = '.ymonmean';
        ymm_out = '';
    elseif strcmp(ymonmean, 'mon')
        ymm_in = '';
        ymm_out = '_mon';
    end
    
    if strcmp(type, 'era5') | strcmp(type, 'erai') | strcmp(type, 'era5c')
        srfc_vars=par.era.vars.srfc;
        for i=1:length(srfc_vars); var = srfc_vars{i};
            file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/%s/srfc/%s_srfc_%s%s.nc', type, type, par.(type).yr_span, ymm_in));
            fullpath=sprintf('%s/%s', file.folder, file.name);

            if strcmp(var, 'zs'); % use orography as surface geopotential if data exists
                if strcmp(type, 'erai')
                    file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/%s/orog/interim_%s.nc', type, 'orog'));
                elseif contains(type, 'era5')
                    file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/%s/orog/%s_%s_%s%s.nc', type, type, 'orog', par.(type).yr_span, ymm_in));
                end
                fullpath=sprintf('%s/%s', file.folder, file.name);
                if exist(fullpath, 'file')
                    srfc.(var)=double(ncread(fullpath, 'z')); srfc.(var) = srfc.(var)/par.g; % divide by g to get height
                    if strcmp(type, 'erai')
                        srfc.(var)=repmat(srfc.(var), [1 1 12]);
                    end
                else % create surface geopotential height using surface pressure data
                    prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s', type, par.(type).yr_span);
                    load(sprintf('%s/grid.mat', prefix)); % read grid data
                    file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/%s/zg/%s_zg_%s%s.nc', type, type, par.(type).yr_span, ymm_in));
                    fullpath=sprintf('%s/%s', file.folder, file.name);
                    zg = double(ncread(fullpath, 'z'));
                    zg = permute(zg, [3 1 2 4]);
                    pb=CmdLineProgressBar("Calculating zs..."); % track progress of this loop
                    for lo = 1:length(grid.dim2.lon)
                        pb.print(lo, length(grid.dim2.lon));
                        for la = 1:length(grid.dim2.lat)
                            for mo = 1:12
                                srfc.zs(lo,la,mo) = interp1(grid.dim3.plev, zg(:,lo,la,mo), srfc.sp(lo,la,mo), 'linear', 'extrap');
                            end
                        end
                    end
                end
            else
                srfc.(var) = double(squeeze(ncread(fullpath, srfc_vars{i})));
            end

        end
        save(sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s/srfc%s.mat', type, par.(type).yr_span, ymm_out), 'srfc', 'srfc_vars');

    elseif strcmp(type, 'merra2')
        srfc_vars=par.merra2.vars.srfc;
        for i=1:length(srfc_vars); var = srfc_vars{i};
            file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/%s/srfc/%s_srfc_%s%s.nc', type, type, par.(type).yr_span, ymm_in));
            fullpath=sprintf('%s/%s', file.folder, file.name);

            if strcmp(var, 'zs'); % use orography as surface geopotential if data exists
                file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/%s/orog/MERRA2_*.nc4', type));
                fullpath=sprintf('%s/%s', file.folder, file.name);
                if exist(fullpath, 'file')
                    srfc.(var) = double(ncread(fullpath, 'PHIS')); srfc.(var)=srfc.(var)/par.g; % convert geopotential to height
                    srfc.(var) = repmat(srfc.(var), [1 1 12]);
                else % create surface geopotential height using surface pressure data
                    prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s', type, par.(type).yr_span);
                    load(sprintf('%s/grid.mat', prefix)); % read grid data
                    file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/%s/zg/%s_zg_%s%s.nc', type, type, par.(type).yr_span, ymm_in));
                    fullpath=sprintf('%s/%s', file.folder, file.name);
                    zg = double(ncread(fullpath, 'H'));
                    zg = permute(zg, [3 1 2 4]);
                    pb=CmdLineProgressBar("Calculating zs..."); % track progress of this loop
                    for lo = 1:length(grid.dim2.lon)
                        pb.print(lo, length(grid.dim2.lon));
                        for la = 1:length(grid.dim2.lat)
                            for mo = 1:12
                                srfc.zs(lo,la,mo) = interp1(grid.dim3.plev, zg(:,lo,la,mo), srfc.PS(lo,la,mo), 'linear', 'extrap');
                            end
                        end
                    end
                end
            else
                srfc.(var) = double(squeeze(ncread(fullpath, srfc_vars{i})));
            end

        end
        save(sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s/srfc%s.mat', type, par.(type).yr_span, ymm_out), 'srfc', 'srfc_vars');

    elseif strcmp(type, 'jra55')
        srfc_vars=par.jra55.vars.srfc;
        for i=1:length(srfc_vars); var = srfc_vars{i};
            if strcmp(var, 'ps')
                file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/%s/srfc/%s_pres_%s%s.nc', type, type, par.(type).yr_span, ymm_in));
                fullpath=sprintf('%s/%s', file.folder, file.name);
                srfc.(var) = double(squeeze(ncread(fullpath, 'PRES_GDS0_SFC_S123')));
            elseif strcmp(var, 'tas')
                file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/%s/srfc/%s_tmp_%s%s.nc', type, type, par.(type).yr_span, ymm_in));
                fullpath=sprintf('%s/%s', file.folder, file.name);
                srfc.(var) = double(squeeze(ncread(fullpath, 'TMP_GDS0_HTGL_S123')));
            elseif strcmp(var, 'hurs')
                file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/%s/srfc/%s_rh_%s%s.nc', type, type, par.(type).yr_span, ymm_in));
                fullpath=sprintf('%s/%s', file.folder, file.name);
                srfc.(var) = double(squeeze(ncread(fullpath, 'RH_GDS0_HTGL_S123')));
            elseif strcmp(var, 'zs'); % use orography as surface geopotential if data exists
                file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/%s/orog/jra55_orog_*%s.nc', type, ymm_in));
                fullpath=sprintf('%s/%s', file.folder, file.name);
                srfc.(var) = double(ncread(fullpath, 'GP_GDS0_SFC')); srfc.(var)=srfc.(var)/par.g; % convert geopotential to height
                srfc.(var) = repmat(srfc.(var), [1 1 12]);
            end

        end
        save(sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s/srfc%s.mat', type, par.(type).yr_span, ymm_out), 'srfc', 'srfc_vars');

    elseif strcmp(type, 'gcm')
        load(sprintf('/project2/tas1/miyawaki/projects/002/data/read/gcm/%s/%s/grid.mat', par.model, par.gcm.clim));

        srfc_vars=par.gcm.vars.srfc;
        for i=1:length(par.gcm.vars.srfc); var = par.gcm.vars.srfc{i};
            file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/gcm/%s/%s_Amon_%s_%s_r1i1p1_*%s.nc', par.model, var, par.model, par.gcm.clim, ymm_in));
            fullpath=sprintf('%s/%s', file.folder, file.name);
            if ~exist(fullpath)
                if strcmp(var, 'hurs')
                    file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/gcm/%s/%s_Amon_%s_%s_r1i1p1_*.nc', par.model, 'hur', par.model, par.gcm.clim));
                    fullpath=sprintf('%s/%s', file.folder, file.name);
                    hur=double(ncread(fullpath, 'hur'));
                    % if ~isequal(grid.dim2.lon, grid.dim3.lon); hur=interp1(grid.dim3.lon, hur, grid.dim2.lon); end; % interpolate to 2D lon if different from 3D
                    % hur=permute(hur, [2 1 3 4]); % bring lat to first dim
                    % if ~isequal(grid.dim2.lat, grid.dim3.lat); hur=interp1(grid.dim3.lat, hur, grid.dim2.lat); end;
                    hur=permute(hur, [3 1 2 4]); % bring plev to first dim
                    pb=CmdLineProgressBar("Calculating hurs..."); % track progress of this loop
                    for id_lon=1:length(grid.dim2.lon)
                        pb.print(id_lon, length(grid.dim2.lon));
                        for id_lat=1:length(grid.dim2.lat)
                            for id_time=1:size(srfc.ps, 3)
                                srfc.hurs(id_lon, id_lat, id_time)=interp1(grid.dim3.plev, hur(:,id_lon,id_lat,id_time), srfc.ps(id_lon, id_lat, id_time), 'linear', 'extrap');
                            end
                        end
                    end
                srfc.hurs = permute(srfc.hurs, [2 1 3]);
                srfc.hurs = interp1(grid.dim2.lat, srfc.hurs, grid.dim3.lat);
                srfc.hurs = permute(srfc.hurs, [2 1 3]);

                elseif strcmp(var, 'zs')
                    file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/%s_raw/%s/%s_*.nc', par.gcm.clim, par.model, 'orog'));
                    fullpath=sprintf('%s/%s', file.folder, file.name);
                    if exist(fullpath, 'file')
                        srfc.(var)=double(ncread(fullpath, 'orog'));
                        srfc.(var)=repmat(srfc.(var), [1 1 12]);
                    else
                        % create surface geopotential height using surface pressure data
                        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s/%s', type, par.model, par.gcm.clim);
                        load(sprintf('%s/grid.mat', prefix)); % read grid data
                        file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/gcm/%s/%s_Amon_%s_%s_r1i1p1_*%s.nc', par.model, 'zg', par.model, par.gcm.clim, ymm_in));
                        fullpath=sprintf('%s/%s', file.folder, file.name);
                        zg = double(ncread(fullpath, 'zg'));
                        if any(contains(par.model, {'GISS-E2-H', 'GISS-E2-R'})) % zg data in GISS-E2-H has an anomalous lat grid
                            zg_lat = double(ncread(fullpath, 'lat'));
                            zg = permute(zg, [2 1 3 4]);
                            zg = interp1(zg_lat, zg, grid.dim3.lat);
                            zg = permute(zg, [2 1 3 4]);
                        end
                        zg = permute(zg, [3 1 2 4]);
                        pb=CmdLineProgressBar("Calculating zs..."); % track progress of this loop
                        for lo = 1:length(grid.dim3.lon)
                            pb.print(lo, length(grid.dim2.lon));
                            for la = 1:length(grid.dim3.lat)
                                for mo = 1:12
                                    srfc.zs(lo,la,mo) = interp1(grid.dim3.plev, zg(:,lo,la,mo), srfc.ps(lo,la,mo), 'linear', 'extrap');
                                end
                            end
                        end
                    end

                else
                    error(sprintf('The file for variable %s does not exist. Check in the raw data folder to see if you forgot to download the file.'))
                end
            else

            srfc.(var)=double(ncread(fullpath, var));

            if ~isequal(grid.dim2.lat, grid.dim3.lat)
                srfc.(var) = permute(srfc.(var), [2 1 3]);
                srfc.(var) = interp1(grid.dim2.lat, srfc.(var), grid.dim3.lat);
                srfc.(var) = permute(srfc.(var), [2 1 3]);
            end


            end
        end
        newdir=sprintf('/project2/tas1/miyawaki/projects/002/data/read/gcm/%s/%s', par.model, par.gcm.clim);
        if ~exist(newdir, 'dir'); mkdir(newdir); end
        filename=sprintf('srfc%s.mat', ymm_out);
        save(sprintf('%s/%s', newdir, filename), 'srfc', 'srfc_vars');
    elseif strcmp(type, 'hahn')
        srfc_vars=par.hahn.vars.srfc;
        for i=1:length(par.hahn.vars.srfc); var = par.hahn.vars.srfc{i};
            if strcmp(var, 'zs'); % create surface geopotential height using surface pressure data
                prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s', type, par.(type).clim);
                load(sprintf('%s/grid.mat', prefix)); % read grid data
                zg = load_zg(type, par);
                zg = permute(zg, [3 1 2 4]);
                pb=CmdLineProgressBar("Calculating zs..."); % track progress of this loop
                for lo = 1:length(grid.dim2.lon)
                    pb.print(lo, length(grid.dim2.lon));
                    for la = 1:length(grid.dim2.lat)
                        for mo = 1:12
                            srfc.zs(lo,la,mo) = interp1(grid.dim3.plev, zg(:,lo,la,mo), srfc.PS(lo,la,mo), 'linear', 'extrap');
                        end
                    end
                end
            else
                fprefix = make_hahn_fprefix(par);
                file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/hahn/lapserateclima/%s.%s.nc', fprefix, var));
                fullpath=sprintf('%s/%s', file.folder, file.name);
                srfc.(var)=double(ncread(fullpath, 'varmo'));
            end

        end
        save(sprintf('/project2/tas1/miyawaki/projects/002/data/read/hahn/%s/srfc%s.mat', par.hahn.clim, ymm_out), 'srfc', 'srfc_vars');
    elseif strcmp(type, 'echam')
        srfc_vars=par.echam.vars.srfc;
        for i=1:length(par.echam.vars.srfc); var = par.echam.vars.srfc{i};
            if contains(par.echam.clim, 'rp000')
                file=dir(sprintf('/project2/tas1/ockham/data11/tas/echam-aiv_rcc_6.1.00p1/%s/BOT_%s_0020_39.nc', par.echam.clim, par.echam.clim));
            else
                file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/echam/BOT*_%s_*%s.nc', par.echam.clim, ymm_in));
            end
            fullpath=sprintf('%s/%s', file.folder, file.name);
            srfc.(var)=double(ncread(fullpath, var));

            if strcmp(var, 'aps'); % create surface geopotential height using surface pressure data
                prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s', type, par.echam.clim);
                load(sprintf('%s/grid.mat', prefix)); % read grid data
                if contains(par.echam.clim, 'rp000')
                    file=dir(sprintf('/project2/tas1/ockham/data11/tas/echam-aiv_rcc_6.1.00p1/%s/ATM_%s_0020_39.nc', par.echam.clim, par.echam.clim));
                else
                    file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/echam/ATM*_%s_*%s.nc', par.echam.clim, ymm_in));
                end
                fullpath=sprintf('%s/%s', file.folder, file.name);
                zg = double(ncread(fullpath, 'geopoth'));
                zg = permute(zg, [3 1 2 4]);
                pb=CmdLineProgressBar("Calculating zs..."); % track progress of this loop
                for lo = 1:length(grid.dim2.lon)
                    pb.print(lo, length(grid.dim2.lon));
                    for la = 1:length(grid.dim2.lat)
                        for mo = 1:12
                            srfc.zs(lo,la,mo) = interp1(grid.dim3.plev, zg(:,lo,la,mo), srfc.aps(lo,la,mo), 'linear', 'extrap');
                        end
                    end
                end
            end

        end
        newdir=sprintf('/project2/tas1/miyawaki/projects/002/data/read/echam/%s', par.echam.clim);
        if ~exist(newdir, 'dir'); mkdir(newdir); end
        filename=sprintf('srfc%s.mat', ymm_out);
        save(sprintf('%s/%s', newdir, filename), 'srfc', 'srfc_vars');
    elseif strcmp(type, 'echam_ml')
        srfc_vars=par.echam.vars.srfc;
        for i=1:length(par.echam.vars.srfc); var = par.echam.vars.srfc{i};
            file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/echam_ml/BOT*%s.nc', ymm_in));
            fullpath=sprintf('%s/%s', file.folder, file.name);
            srfc.(var)=double(ncread(fullpath, var));

            if strcmp(var, 'aps'); % create surface geopotential height using surface pressure data
                prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s', type);
                load(sprintf('%s/grid.mat', prefix)); % read grid data
                file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/echam_ml/ATM*%s.nc', ymm_in));
                fullpath=sprintf('%s/%s', file.folder, file.name);
                zg = double(ncread(fullpath, 'geopoth'));
                srfc.zs(:,:,:) = squeeze(zg(:,:,1,:));
            end

        end
        newdir=sprintf('/project2/tas1/miyawaki/projects/002/data/read/echam_ml');
        if ~exist(newdir, 'dir'); mkdir(newdir); end
        filename=sprintf('srfc%s.mat', ymm_out);
        save(sprintf('%s/%s', newdir, filename), 'srfc', 'srfc_vars');
    elseif strcmp(type, 'echam_pl')
        srfc_vars=par.echam.vars.srfc;
        for i=1:length(par.echam.vars.srfc); var = par.echam.vars.srfc{i};
            file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/echam_pl/BOT*%s.nc', ymm_in));
            fullpath=sprintf('%s/%s', file.folder, file.name);
            srfc.(var)=double(ncread(fullpath, var));

            if strcmp(var, 'aps'); % create surface geopotential height using surface pressure data
                prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s', type);
                load(sprintf('%s/grid.mat', prefix)); % read grid data
                file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/echam_pl/ATM*%s.nc', ymm_in));
                fullpath=sprintf('%s/%s', file.folder, file.name);
                zg = double(ncread(fullpath, 'geopoth'));
                srfc.zs(:,:,:) = squeeze(zg(:,:,1,:));
            end

        end
        newdir=sprintf('/project2/tas1/miyawaki/projects/002/data/read/echam_pl');
        if ~exist(newdir, 'dir'); mkdir(newdir); end
        filename=sprintf('srfc%s.mat', ymm_out);
        save(sprintf('%s/%s', newdir, filename), 'srfc', 'srfc_vars');
    end
end
