function save_lfrac(type, par)

    foldername = make_savedir_proc(type, par);
    prefix = make_prefix(type, par);

    load(sprintf('%s/grid.mat', prefix)); % read grid data
    if strcmp(par.lat_interp, 'std')
        lat = par.lat_std;
    else
        lat = grid.dim3.lat;
    end

    if strcmp(type, 'echam')
        if any(contains(par.echam.clim, {'rp000'}))
            ;
        else
            load(sprintf('%s/sftlf.mat', prefix)); % load land fraction data
        end
    else
        load(sprintf('%s/sftlf.mat', prefix)); % load land fraction data
    end

    if any(strcmp(type, {'era5c', 'erai', 'era5'}))
        lfrac = sftlf;
    elseif strcmp(type, 'echam') & any(contains(par.echam.clim, {'rp000'}))
        lfrac = zeros([length(grid.dim2.lon) length(grid.dim2.lat)]);
    elseif any(strcmp(type, {'gcm', 'merra2c', 'merra2', 'jra55', 'echam'})) % repeat to monthly dimension
        lfrac=sftlf; lfrac=repmat(lfrac,[1 1 12]);
    else
        lfrac = remove_ocean(lat, grid.dim3.lon, 12);
    end

    % save lfracs
    printname = [foldername 'lfrac'];
    if ~exist(foldername, 'dir')
        mkdir(foldername)
    end
    save(printname, 'lfrac');

end % compute lfracs once and save it for later use
