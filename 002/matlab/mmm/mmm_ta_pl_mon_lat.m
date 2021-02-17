function mmm_ta_mon_lat(type, par)

    lat = par.lat;
    for l = {'lo'}; land=l{1};
        ta_list.(land) = nan(length(par.gcm_models), length(par.lat), 12, length(par.pa));
    end

    pb = CmdLineProgressBar("Creating the multi-model mean...");
    for k=1:length(par.gcm_models); par.model = par.gcm_models{k};
        pb.print(k, length(par.gcm_models)); % output progress of moist adiabat calculonion

        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s/%s', type, par.model, par.gcm.clim);
        grid0=load(sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s/%s/grid.mat', type, par.model, par.gcm.clim));
        ta0=load(sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s/%s/ta_pl_mon_lat.mat', type, par.model, par.gcm.clim, 'native'));

        foldername = sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s/%s/', type, par.outname, par.gcm.clim, par.lat_interp);
        load(sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s/%s/grid.mat', type, par.outname, par.gcm.clim));

        for l = {'lo'}; land=l{1};
            % interpolate native grid data to standard grid
            ta0i = interp1(grid0.grid.dim3.lat, ta0.ta.(land), par.lat, 'linear', 'extrap'); % interpolate to standard lat grid
            ta0i = permute(ta0i, [3 1 2]); % (lev x lat x mon)
            ta0i = interp1(grid0.grid.dim3.plev, ta0i, par.pa, 'linear'); % inteprolate to standard pressure grid
            ta0i = permute(ta0i, [2 3 1]); % return to (lat x mon x lev)
            ta_list.(land)(k,:,:,:) = ta0i;
            % ta0ii = nan([length(par.pa) length(par.lat) 12]);
            % for ilat = 1:length(par.lat)
            %     for imon = 1:12
            %         if all(isnan(ta0i(:,ilat,imon)))
            %             ta0ii(:,ilat,imon) = nan([length(par.pa),1,1]);
            %         else
            %             ta0ii(:,ilat,imon) = interp1(grid0.grid.dim3.plev, ta0i(:,ilat,imon), par.pa, 'spline'); % inteprolate to standard pressure grid
            %         end
            %     end
            % end
            % ta0ii = permute(ta0ii, [2 3 1]); % return to (lat x mon x lev)
            % ta_list.(land)(k,:,:,:) = ta0ii;
        end % land

    end % models

    for l = {'lo'}; land=l{1};
        ta.(land) = squeeze(nanmean(ta_list.(land),1));
    end

    printname = [foldername 'ta_pl_mon_lat'];
    if ~exist(foldername, 'dir')
        mkdir(foldername)
    end
    save(printname, 'ta', 'lat', '-v7.3');

end
