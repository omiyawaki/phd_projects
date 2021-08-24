function read_albedo(type, par)
% read sea ice depth
    if strcmp(type, 'echam')

        if contains(par.echam.clim, 'rp000')
            file=dir(sprintf('/project2/tas1/ockham/data11/tas/echam-aiv_rcc_6.1.00p1/%s/BOT_%s_0020_39.nc', par.echam.clim, par.echam.clim));
        else
            file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/echam/BOT*_%s_*.ymonmean.nc', par.echam.clim));
        end
        fullpath=sprintf('%s/%s', file.folder, file.name);
        albedo=double(ncread(fullpath, 'albedo'));
        newdir=sprintf('/project2/tas1/miyawaki/projects/002/data/read/echam/%s', par.echam.clim);
        if ~exist(newdir, 'dir'); mkdir(newdir); end
        filename='albedo.mat';
        save(sprintf('%s/%s', newdir, filename), 'albedo');

    elseif any(strcmp(type, {'era5', 'era5c', 'erai'}))

        albedo = double(ncread(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/%s/albedo/%s_albedo_%s.ymonmean.nc', type, type, par.(type).yr_span), 'fal'));
        save(sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s/albedo.mat', type, par.(type).yr_span), 'albedo');

    elseif strcmp(type, 'hahn')
        fprefix = make_hahn_fprefix(par);

        % net SW surface
        file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/hahn/lapserateclima/%s.%s.nc', fprefix, 'FSNS')); 
        fullpath=sprintf('%s/%s', file.folder, file.name);
        fsns=double(ncread(fullpath, 'varmo'));
        
        % downward SW surface
        file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/hahn/lapserateclima/%s.%s.nc', fprefix, 'FSDS')); 
        fullpath=sprintf('%s/%s', file.folder, file.name);
        fsds=double(ncread(fullpath, 'varmo'));
        
        % compute upward SW surface
        fsus = fsds - fsns;

        % compute surface albedo
        albedo = fsus./fsds;

        save(sprintf('/project2/tas1/miyawaki/projects/002/data/read/hahn/%s/albedo.mat', par.hahn.clim), 'albedo');

    % else
    %     if any(strcmp(type, {'era5', 'era5c', 'erai'}))
    %         prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s', type);
    %         prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s', type);
    %     elseif strcmp(type, 'gcm')
    %         prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s/%s', type, par.model, par.gcm.clim);
    %         prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s', type, par.model, par.gcm.clim);
    %     end
    %     load(sprintf('%s/grid.mat', prefix)); % read grid data
    %     load(sprintf('%s/rad.mat', prefix)); % read radiation data

    %     if strcmp(type, 'erai') | strcmp(type, 'era5') | strcmp(type, 'era5c')
    %         alb = rad.rsus./rad.rsds; % TODO download shortwave radiation up and down separately for ERA
    %     elseif strcmp(type, 'gcm')
    %         alb = rad.rsus./rad.rsds;
    %     end

    %     if any(strcmp(type, {'era5', 'era5c', 'erai'})); newdir=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s', type);
    %     elseif strcmp(type, 'gcm'); newdir=sprintf('/project2/tas1/miyawaki/projects/002/data/read/gcm/%s/%s', par.model, par.gcm.clim); end;
    %     if ~exist(newdir, 'dir'); mkdir(newdir); end
    %     filename='alb.mat';
    %     save(sprintf('%s/%s', newdir, filename), 'alb', '-v7.3');

    end
end
