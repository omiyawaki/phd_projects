function proc_dtdz_mon_lat(type, par)
% calculate moist adiabats
    if strcmp(type, 'era5') | strcmp(type, 'erai') | strcmp(type, 'era5c')
        foldername = sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s/', type, par.(type).yr_span, par.lat_interp);
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s', type, par.(type).yr_span);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s', type, par.(type).yr_span);
    elseif strcmp(type, 'gcm')
        foldername = sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s/%s/', type, par.model, par.gcm.clim, par.lat_interp);
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s/%s', type, par.model, par.gcm.clim);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s', type, par.model, par.gcm.clim);
    elseif strcmp(type, 'echam')
        foldername = sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s/', type, par.echam.clim, par.lat_interp);
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s', type, par.echam.clim);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s', type, par.echam.clim);
    end
    load(sprintf('%s/grid.mat', prefix)); % read grid data
    load(sprintf('%s/srfc.mat', prefix)); % read surface variable data
    load(sprintf('%s/dtdzsi.mat', prefix));
    % load(sprintf('%s/%s/masks.mat', prefix_proc, par.lat_interp)); % load land and ocean masks

    % % Land/ocean filter 3D variables
    % mask.land_vert = repmat(mask.land, [1 1 1 size(dtdzsi, 3)]); % expand land mask to vertical dim
    % mask.land_vert = permute(mask.land, [1 2 4 3]); % place vertical dim where it belongs
    % mask.ocean_vert = repmat(mask.ocean, [1 1 1 size(dtdzsi, 3)]); % expand ocean mask to vertical dim
    % mask.ocean_vert = permute(mask.ocean, [1 2 4 3]); % place vertical dim where it belongs

    dtdz_sm.lo = dtdzsi;
    % dtdz_sm.l = dtdz_sm.lo.*mask.ocean_vert; % filter dtdz with surface mask
    % dtdz_sm.o = dtdz_sm.lo.*mask.land_vert; % filter dtdz with surface mask
    dtdz_sm.lo = permute(dtdz_sm.lo, [1 2 4 3]); % bring plev to last dimension
    % dtdz_sm.l = permute(dtdz_sm.l, [1 2 4 3]); % bring plev to last dimension
    % dtdz_sm.o = permute(dtdz_sm.o, [1 2 4 3]); % bring plev to last dimension

    if strcmp(par.lat_interp, 'std')
        lat = par.lat_std;
    else
        lat = grid.dim3.lat;
    end

    pb = CmdLineProgressBar("Calculating lapse rates...");
    % for l = {'lo', 'l', 'o'}; land = l{1};
    for l = {'lo'}; land = l{1};
        for ilat = 1:length(lat);
            pb.print(ilat, length(lat));
            for imon = 1:12;
                dtdz.(land) = squeeze(nanmean(dtdz_sm.(land))); % zonal mean
            end

            for ilon = 1:length(grid.dim3.lon);
                for t = {'ann', 'djf', 'jja', 'mam', 'son'}; time = t{1};
                    if strcmp(time, 'ann')
                        dtdz_t.(land).(time) = squeeze(nanmean(dtdz_sm.(land), 3));
                    elseif strcmp(time, 'djf')
                        dtdz_shift.(land) = circshift(dtdz_sm.(land), 1, 3);
                        dtdz_t.(land).(time) = squeeze(nanmean(dtdz_shift.(land)(:,:,1:3,:), 3));
                    elseif strcmp(time, 'jja')
                        dtdz_t.(land).(time) = squeeze(nanmean(dtdz_sm.(land)(:,:,6:8,:), 3));
                    elseif strcmp(time, 'mam')
                        dtdz_t.(land).(time) = squeeze(nanmean(dtdz_sm.(land)(:,:,3:5,:), 3));
                    elseif strcmp(time, 'son')
                        dtdz_t.(land).(time) = squeeze(nanmean(dtdz_sm.(land)(:,:,9:11,:), 3));
                    end
                end
            end
        end
    end

    % save data into mat file
    printname = [foldername 'dtdz_mon_lat.mat'];
    if ~exist(foldername, 'dir')
        mkdir(foldername)
    end
    save(printname, 'dtdz');

    printname = [foldername 'dtdz_lon_lat.mat'];
    if ~exist(foldername, 'dir')
        mkdir(foldername)
    end
    save(printname, 'dtdz_t');

end % compute mon x lat moist adiabat field
