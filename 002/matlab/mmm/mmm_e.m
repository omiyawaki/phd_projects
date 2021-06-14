function mmm_e(type, par)
    lat = par.lat;

    % output info
    foldername = make_savedir(type, par);
    filename = 'e';
    load(sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s/%s/grid.mat', type, par.outname, par.clim));

    par.lat_interp = 'native'; % input files will be in native grid

    e_vec = {'e'};
    for s = 1:length(e_vec); svn = e_vec{s};
        e_list.(svn) = nan(length(par.model_list), length(grid.dim2.lon), length(grid.dim2.lat), 12);
        e_mmm.(svn) = nan(length(grid.dim2.lon), length(grid.dim2.lat), 12);
        e_std.(svn) = nan(length(grid.dim2.lon), length(grid.dim2.lat), 12);
        e_min.(svn) = nan(length(grid.dim2.lon), length(grid.dim2.lat), 12);
        e_max.(svn) = nan(length(grid.dim2.lon), length(grid.dim2.lat), 12);
        e_25.(svn) = nan(length(grid.dim2.lon), length(grid.dim2.lat), 12);
        e_75.(svn) = nan(length(grid.dim2.lon), length(grid.dim2.lat), 12);
    end


    pb = CmdLineProgressBar("Creating the multi-model mean...");
    for k=1:length(par.model_list); par.model = par.model_list{k};
        pb.print(k, length(par.model_list)); % output progress of moist adiabat calculonion

        if strcmp(type, 'gcm')
            type_in = type;
        else
            type_in = par.model;
        end

        % input info
        prefix = make_prefix(type_in, par);
        grid0 = load(sprintf('%s/grid.mat', prefix));
        srfc_0 = load(sprintf('%s/srfc.mat', prefix));

        if any(strcmp(type_in, {'era5c', 'era5', 'erai'})); % convert dew point temp to humidity
            dew = srfc_0.srfc.d2m;
            ps = srfc_0.srfc.sp;
            tas = srfc_0.srfc.t2m;
            srfc_0.esat = 1e2*(1.0007 + 3.46e-6*ps/100).*6.1121.*exp(17.966*(tas-273.15)./(247.15+(tas-273.15)));
            srfc_0.e = calc_esat(dew, par.frz);
        elseif any(strcmp(type_in, {'merra2', 'merra2c'})); % convert specific humidity to relative humidity
            q = srfc_0.srfc.QV2M;
            ps = srfc_0.srfc.PS;
            tas = srfc_0.srfc.T2M;
            srfc_0.esat = 1e2*(1.0007 + 3.46e-6*ps/100).*6.1121.*exp(17.966*(tas-273.15)./(247.15+(tas-273.15)));
            srfc_0.e = calc_e(ps, q, par);
        elseif any(strcmp(type_in, {'jra55', 'gcm'})); % convert specific humidity to relative humidity
            hurs = srfc_0.srfc.hurs;
            ps = srfc_0.srfc.ps;
            tas = srfc_0.srfc.tas;
            srfc_0.esat = 1e2*(1.0007 + 3.46e-6*ps/100).*6.1121.*exp(17.966*(tas-273.15)./(247.15+(tas-273.15)));
            srfc_0.e = srfc_0.esat.*hurs/100;
        end

        for s = 1:length(e_vec); svn = e_vec{s};
            tmp = srfc_0.(svn);
            tmp = interp1(grid0.grid.dim2.lon, tmp, grid.dim2.lon);
            tmp = permute(tmp, [2 1 3]);
            if size(grid0.grid.dim2.lat,1) ~= size(tmp,1)
                tmp = interp1(grid0.grid.dim3.lat, tmp, grid.dim2.lat);
            else
                tmp = interp1(grid0.grid.dim2.lat, tmp, grid.dim2.lat);
            end
            tmp = permute(tmp, [2 1 3]);
            e_list.(svn)(k,:,:,:) = tmp;
        end

    end % models

    for s = 1:length(e_vec); svn = e_vec{s};
        e_mmm.(svn) = squeeze(nanmean(e_list.(svn), 1));
        e_std.(svn) = squeeze(nanstd(e_list.(svn), 1));
        e_min.(svn) = squeeze(min(e_list.(svn), [], 1));
        e_max.(svn) = squeeze(max(e_list.(svn), [], 1));
        e_25.(svn) = squeeze(prctile(e_list.(svn), [25], 1));
        e_75.(svn) = squeeze(prctile(e_list.(svn), [75], 1));
    end

    e = e_mmm;

    % save data into mat file
    if ~exist(foldername, 'dir')
        mkdir(foldername)
    end
    save(sprintf('%s/%s', foldername, filename), 'e', 'e_std', 'e_min', 'e_max', 'e_25', 'e_75', 'lat', '-v7.3');

end
 
