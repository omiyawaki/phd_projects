function mmm_srfc(type, par)
    lat = par.lat;

    % output info
    foldername = make_savedir(type, par);
    filename = 'srfc.mat';
    load(sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s/%s/%s/grid.mat', type, par.outname, par.clim, par.(type).yr_span));

    par.lat_interp = 'native'; % input files will be in native grid

    srfc_vec = par.rea.vars.srfc;
    for s = 1:length(srfc_vec); svn = srfc_vec{s};
        if ~strcmp(svn, 'hurs') % omit hurs because not all reanalyses have this
            srfc_list.(svn) = nan(length(par.model_list), length(grid.dim2.lon), length(grid.dim2.lat), 12);
            srfc_mmm.(svn) = nan(length(grid.dim2.lon), length(grid.dim2.lat), 12);
            srfc_std.(svn) = nan(length(grid.dim2.lon), length(grid.dim2.lat), 12);
            srfc_min.(svn) = nan(length(grid.dim2.lon), length(grid.dim2.lat), 12);
            srfc_max.(svn) = nan(length(grid.dim2.lon), length(grid.dim2.lat), 12);
            srfc_25.(svn) = nan(length(grid.dim2.lon), length(grid.dim2.lat), 12);
            srfc_75.(svn) = nan(length(grid.dim2.lon), length(grid.dim2.lat), 12);
        end
    end


    pb = CmdLineProgressBar("Creating the multi-model mean...");
    for k=1:length(par.model_list); par.model = par.model_list{k};
        pb.print(k, length(par.model_list)); % output progress of moist adiabat calculonion

        if strcmp(type, 'gcm')
            type_in = type;
        else
            type_in = par.model
        end

        % input info
        prefix = make_prefix(type_in, par);
        grid0 = load(sprintf('%s/grid.mat', prefix));
        srfc_0 = load(sprintf('%s/%s', prefix, filename));

        srfc_vec0 = par.(type_in).vars.srfc;
        for s = 1:length(srfc_vec0);
            svn_in = srfc_vec0{s}
            svn = srfc_vec{s}
            if any(strcmp(svn_in, {'d2m'})); % convert dew point temp to humidity
                dew = srfc_0.srfc.(svn_in);
                tas = srfc_0.srfc.t2m;
                e = calc_esat(dew, par.frz);
                esat = calc_esat(tas, par.frz);
                tmp = 100 * e./esat;
            elseif any(strcmp(svn_in, {'QV2M'})); % convert specific humidity to relative humidity
                q = srfc_0.srfc.(svn_in);
                ps = srfc_0.srfc.PS;
                tas = srfc_0.srfc.T2M;
                e = calc_e(ps, q, par);
                esat = calc_esat(tas, par.frz);
                tmp = 100 * e./esat;
            else
                tmp = srfc_0.srfc.(svn_in);
            end
            tmp = interp1(grid0.grid.dim2.lon, tmp, grid.dim2.lon);
            tmp = permute(tmp, [2 1 3]);
            if size(grid0.grid.dim2.lat,1) ~= size(tmp,1)
                tmp = interp1(grid0.grid.dim3.lat, tmp, grid.dim2.lat);
            else
                tmp = interp1(grid0.grid.dim2.lat, tmp, grid.dim2.lat);
            end
            tmp = permute(tmp, [2 1 3]);
            srfc_list.(svn)(k,:,:,:) = tmp;
        end

    end % models

    srfc_vec = par.rea.vars.srfc;
    for s = 1:length(srfc_vec); svn = srfc_vec{s};
        srfc_mmm.(svn) = squeeze(nanmean(srfc_list.(svn), 1));
        srfc_std.(svn) = squeeze(nanstd(srfc_list.(svn), 1));
        srfc_min.(svn) = squeeze(min(srfc_list.(svn), [], 1));
        srfc_max.(svn) = squeeze(max(srfc_list.(svn), [], 1));
        srfc_25.(svn) = squeeze(prctile(srfc_list.(svn), [25], 1));
        srfc_75.(svn) = squeeze(prctile(srfc_list.(svn), [75], 1));
    end

    srfc = srfc_mmm;

    % save energy flux data into mat file
    if ~exist(foldername, 'dir')
        mkdir(foldername)
    end
    save(sprintf('%s/%s', foldername, filename), 'srfc', 'srfc_std', 'srfc_min', 'srfc_max', 'srfc_25', 'srfc_75', 'lat', '-v7.3');

end
 
