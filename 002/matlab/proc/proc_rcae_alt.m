function proc_rcae_alt(type, par)
    foldername = make_savedir_proc_ep(type, par);
    prefix = make_prefix(type, par);
    prefix_proc = make_prefix_proc(type, par);

    load(sprintf('%s/vh.mat', prefix_proc)); % load atmospheric heat transport
    load(sprintf('%s/vh_mon.mat', prefix_proc)); % load atmospheric heat transport

    for f = {'flux', 'flux_t', 'flux_z'}; ftype = f{1};
        % if strcmp(type, 'era5') | strcmp(type, 'erai') | strcmp(type, 'era5c')
        %     filename = sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s.mat', type, par.lat_interp, ftype); % read ERA5 zonally averaged flux
        %     foldername = sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/eps_%g_ga_%g/', type, par.lat_interp, par.ep, par.ga);
        %     prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s', type, par.(type).yr_span);
        %     prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s', type, par.lat_interp);
        % elseif strcmp(type, 'gcm')
        %     filename = sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s/%s/%s.mat', type, par.model, par.gcm.clim, par.lat_interp, ftype); % read gcm zonally averaged flux
        %     foldername = sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s/%s/eps_%g_ga_%g/', type, par.model, par.gcm.clim, par.lat_interp, par.ep, par.ga);
        %     prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s/%s', type, par.model, par.gcm.clim);
        %     prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s/%s', type, par.model, par.gcm.clim, par.lat_interp);
        % elseif strcmp(type, 'echam')
        %     filename = sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s/%s.mat', type, par.echam.clim, par.lat_interp, ftype); % read echam zonally averaged flux
        %     foldername = sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s/eps_%g_ga_%g/', type, par.echam.clim, par.lat_interp, par.ep, par.ga);
        %     prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s', type, par.echam.clim);
        %     prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s', type, par.echam.clim, par.lat_interp);
        % end
        % if ~exist(filename); error(sprintf('Data does not exist. Please run proc_%s.m first.', ftype)); else
        %     load(filename);
        % end


        % load(sprintf('%s/masks.mat', prefix_proc)); % load land and ocean masks
        load(sprintf('%s/%s.mat', prefix_proc, ftype)); % load land and ocean masks

        % identify locations of RCE and RAE
        if strcmp(ftype, 'flux')
            rcae_alt = def_rcae_alt(type, flux, vh_mon, par); % lon x lat x time structure
            rcae_alt_rc = def_rcae_alt_recomp_r1(type, flux, vh_mon, par); % lon x lat x time structure
            printname = [foldername 'rcae_alt.mat'];
            printname_rc = [foldername 'rcae_alt_rc.mat'];
        elseif strcmp(ftype, 'flux_t')
            % for l = {'lo', 'l', 'o'}; land = l{1};
            for l = {'lo'}; land = l{1};
                % lon x lat structure over various time averages
                % for t = {'ann', 'djf', 'jja', 'mam', 'son'}; time = t{1};
                for t = {'ann'}; time = t{1};
                    rcae_alt_t.(land).(time) = def_rcae_alt(type, flux_t.(land).(time), vh.(land).(time), par);
                    rcae_alt_rc_t.(land).(time) = def_rcae_alt_recomp_r1(type, flux_t.(land).(time), vh.(land).(time), par);
                end % time
            end % land

            printname = [foldername 'rcae_alt_t.mat'];
            printname_rc = [foldername 'rcae_alt_rc_t.mat'];
        elseif strcmp(ftype, 'flux_z') % show lat x mon structure of RCAE_ALT
            % for l = {'lo', 'l', 'o'}; land = l{1};
            for l = {'lo'}; land = l{1};
                rcae_alt_z.(land) = def_rcae_alt(type, flux_z.(land), vh_mon.(land), par);
                rcae_alt_rc_z.(land) = def_rcae_alt_recomp_r1(type, flux_z.(land), vh_mon.(land), par);
                printname = [foldername 'rcae_alt_z.mat'];
                printname_rc = [foldername 'rcae_alt_rc_z.mat'];
            end
        end

        % save rcae_alt data
        if ~exist(foldername, 'dir')
            mkdir(foldername)
        end
        if strcmp(ftype, 'flux'); save(printname, 'rcae_alt', 'lat', '-v7.3'); save(printname_rc, 'rcae_alt_rc', 'lat', '-v7.3');
        elseif strcmp(ftype, 'flux_t'); save(printname, 'rcae_alt_t', 'lat', '-v7.3'); save(printname_rc, 'rcae_alt_rc_t', 'lat', '-v7.3');
        elseif strcmp(ftype, 'flux_z'); save(printname, 'rcae_alt_z', 'lat', '-v7.3'); save(printname_rc, 'rcae_alt_rc_z', 'lat', '-v7.3'); end
    end

end % process RCE/RAE regimes (all divergence allowed for RCE)
