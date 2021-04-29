function proc_thetaeq_mon_lat(type, par)
    
    % if strcmp(type, 'era5') | strcmp(type, 'erai') | strcmp(type, 'era5c')
    %     foldername = sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s', type, par.(type).yr_span, par.lat_interp);
    %     prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s', type, par.(type).yr_span);
    %     prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s', type, par.(type).yr_span);
    % elseif strcmp(type, 'gcm')
    %     foldername = sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s/%s/', type, par.model, par.gcm.clim, par.lat_interp);
    %     prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s/%s', type, par.model, par.gcm.clim);
    %     prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s', type, par.model, par.gcm.clim);
    % elseif strcmp(type, 'echam')
    %     foldername = sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s/', type, par.echam.clim, par.lat_interp);
    %     prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s', type, par.echam.clim);
    %     prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s', type, par.echam.clim);
    % elseif strcmp(type, 'echam_ml')
    %     foldername = sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s', type, par.(type).yr_span, par.lat_interp);
    %     prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s', type, par.(type).yr_span);
    %     prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s', type, par.(type).yr_span);
    % elseif strcmp(type, 'echam_pl')
    %     foldername = sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s', type, par.(type).yr_span, par.lat_interp);
    %     prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s', type, par.(type).yr_span);
    %     prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s', type, par.(type).yr_span);
    % end

    prefix = make_prefix(type, par);
    foldername = make_savedir_proc(type, par);

    load(sprintf('%s/grid.mat', prefix)); % read grid data
    load(sprintf('%s/thetaeq_si.mat', prefix)); thetaeqsi_orig = thetaeq_si.spl; clear thetaeq_si; % read temp in si coordinates
    load(sprintf('%s/pa_si.mat', prefix)); pasi_orig = pa_si; clear pa_si; % read temp in si coordinates
    load(sprintf('%s/srfc.mat', prefix)); % load surface data
    % load(sprintf('%s/%s/masks.mat', prefix_proc, par.lat_interp)); % load land and ocean masks

    if strcmp(par.lat_interp, 'std')
        lat = par.lat_std;
    else
        lat = grid.dim3.lat;
    end

    % interpolate thetaeq to sthetaeqndard lat grid
    thetaeqsi_orig = permute(thetaeqsi_orig, [2 1 3 4]);
    thetaeqsi_orig = interp1(grid.dim3.lat, thetaeqsi_orig, lat);
    thetaeqsi_orig = permute(thetaeqsi_orig, [2 1 3 4]);

    pasi_orig = permute(pasi_orig, [2 1 3 4]);
    pasi_orig = interp1(grid.dim3.lat, pasi_orig, lat);
    pasi_orig = permute(pasi_orig, [2 1 3 4]);

    thetaeqsi_sm.lo = thetaeqsi_orig; % surface is already masked in sthetaeqndard sigma coordinates
    pasi_sm.lo = pasi_orig; % surface is already masked in spndard sigma coordinates

    thetaeqsi_sm.lo = permute(thetaeqsi_sm.lo, [1 2 4 3]); % bring plev to last dimension

    pasi_sm.lo = permute(pasi_sm.lo, [1 2 4 3]); % bring plev to last dimension

    % mask_t.land = nanmean(mask.land, 3);
    % mask_t.ocean = nanmean(mask.ocean, 3);

    for l = {'lo'}; land = l{1}; % over land, over ocean, or both
        thetaeqsi.(land)= squeeze(mean(thetaeqsi_sm.(land), 1)); % zonal average
        pasi.(land)= squeeze(mean(pasi_sm.(land), 1)); % zonal average
    end

    for l = {'lo'}; land = l{1}; % over land, over ocean, or both
        % thetaeqke time averages
        for t = {'ann', 'djf', 'jja', 'mam', 'son'}; time = t{1};
            if strcmp(time, 'ann')
                thetaeqsi_t.(land).(time) = squeeze(mean(thetaeqsi_sm.(land), 3));
                pasi_t.(land).(time) = squeeze(mean(pasi_sm.(land), 3));
            elseif strcmp(time, 'djf')
                thetaeqsi_shift.(land) = circshift(thetaeqsi_sm.(land), 1, 3);
                pasi_shift.(land) = circshift(pasi_sm.(land), 1, 3);
                thetaeqsi_t.(land).(time) = squeeze(mean(thetaeqsi_shift.(land)(:,:,1:3,:), 3));
                pasi_t.(land).(time) = squeeze(mean(pasi_shift.(land)(:,:,1:3,:), 3));
            elseif strcmp(time, 'jja')
                thetaeqsi_t.(land).(time) = squeeze(mean(thetaeqsi_sm.(land)(:,:,6:8,:), 3));
                pasi_t.(land).(time) = squeeze(mean(pasi_sm.(land)(:,:,6:8,:), 3));
            elseif strcmp(time, 'mam')
                thetaeqsi_t.(land).(time) = squeeze(mean(thetaeqsi_sm.(land)(:,:,3:5,:), 3));
                pasi_t.(land).(time) = squeeze(mean(pasi_sm.(land)(:,:,3:5,:), 3));
            elseif strcmp(time, 'son')
                thetaeqsi_t.(land).(time) = squeeze(mean(thetaeqsi_sm.(land)(:,:,9:11,:), 3));
                pasi_t.(land).(time) = squeeze(mean(pasi_sm.(land)(:,:,9:11,:), 3));
            end
        end
    end

    % save filtered dathetaeq
    printname = [foldername 'thetaeq_mon_lat'];
    if ~exist(foldername, 'dir')
        mkdir(foldername)
    end
    if par.do_surf; save(printname, 'thetaeqsi', 'lat');
    else save(printname, 'thetaeqsi', 'pasi', 'lat', '-v7.3'); end

    printname = [foldername 'thetaeq_lon_lat'];
    if ~exist(foldername, 'dir')
        mkdir(foldername)
    end
    if par.do_surf; save(printname, 'thetaeqsi_t', 'lat');
    else save(printname, 'thetaeqsi_t', 'pasi_t', 'lat', '-v7.3'); end
end % compute mon x lat temperature field
