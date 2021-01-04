function proc_rcae(type, par)
    for f = {'flux', 'flux_t', 'flux_z'}; ftype = f{1};
        if strcmp(type, 'era5') | strcmp(type, 'erai') | strcmp(type, 'era5c')
            filename = sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s.mat', type, par.lat_interp, ftype); % read ERA5 zonally averaged flux
            foldername = sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/eps_%g_ga_%g/', type, par.lat_interp, par.ep, par.ga);
            prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s', type, par.(type).yr_span);
            prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s', type, par.lat_interp);
        elseif strcmp(type, 'gcm')
            filename = sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s/%s/%s.mat', type, par.model, par.gcm.clim, par.lat_interp, ftype); % read gcm zonally averaged flux
            foldername = sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s/%s/eps_%g_ga_%g/', type, par.model, par.gcm.clim, par.lat_interp, par.ep, par.ga);
            prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s/%s', type, par.model, par.gcm.clim);
            prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s/%s', type, par.model, par.gcm.clim, par.lat_interp);
        elseif strcmp(type, 'echam')
            filename = sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s/%s.mat', type, par.echam.clim, par.lat_interp, ftype); % read echam zonally averaged flux
            foldername = sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s/eps_%g_ga_%g/', type, par.echam.clim, par.lat_interp, par.ep, par.ga);
            prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s', type, par.echam.clim);
            prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s', type, par.echam.clim, par.lat_interp);
        end
        if ~exist(filename); error(sprintf('Data does not exist. Please run proc_%s.m first.', ftype)); else
            load(filename);
        end

        % load(sprintf('%s/masks.mat', prefix_proc)); % load land and ocean masks
        load(sprintf('%s/vh.mat', prefix_proc)); % load atmospheric heat transport
        load(sprintf('%s/vh_mon.mat', prefix_proc)); % load atmospheric heat transport

        % identify locations of RCE and RAE
        if strcmp(ftype, 'flux')
            rcae = def_rcae(type, flux, vh_mon, par); % lon x lat x time structure
            rcae_rc = def_rcae_recomp_r1(type, flux, vh_mon, par); % lon x lat x time structure
            printname = [foldername 'rcae.mat'];
            printname_rc = [foldername 'rcae_rc.mat'];
        elseif strcmp(ftype, 'flux_t')
            % for l = {'lo', 'l', 'o'}; land = l{1};
            for l = {'lo'}; land = l{1};
                % lon x lat structure over various time averages
                for t = {'ann', 'djf', 'jja', 'mam', 'son'}; time = t{1};
                    rcae_t.(land).(time) = def_rcae(type, flux_t.(land).(time), vh.(land).(time), par);
                    rcae_rc_t.(land).(time) = def_rcae_recomp_r1(type, flux_t.(land).(time), vh.(land).(time), par);
                end % time
            end % land

            printname = [foldername 'rcae_t.mat'];
            printname_rc = [foldername 'rcae_rc_t.mat'];
        elseif strcmp(ftype, 'flux_z') % show lat x mon structure of RCAE
            % for l = {'lo', 'l', 'o'}; land = l{1};
            for l = {'lo'}; land = l{1};
                rcae_z.(land) = def_rcae(type, flux_z.(land), vh_mon.(land), par);
                rcae_rc_z.(land) = def_rcae_recomp_r1(type, flux_z.(land), vh_mon.(land), par);
                printname = [foldername 'rcae_z.mat'];
                printname_rc = [foldername 'rcae_rc_z.mat'];
            end
        end

        % save rcae data
        if ~exist(foldername, 'dir')
            mkdir(foldername)
        end
        if strcmp(ftype, 'flux'); save(printname, 'rcae', 'lat', '-v7.3'); save(printname_rc, 'rcae_rc', 'lat', '-v7.3');
        elseif strcmp(ftype, 'flux_t'); save(printname, 'rcae_t', 'lat', '-v7.3'); save(printname_rc, 'rcae_rc_t', 'lat', '-v7.3');
        elseif strcmp(ftype, 'flux_z'); save(printname, 'rcae_z', 'lat', '-v7.3'); save(printname_rc, 'rcae_rc_z', 'lat', '-v7.3'); end
    end

end % process RCE/RAE regimes
