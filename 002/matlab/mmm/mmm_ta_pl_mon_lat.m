function mmm_ta_mon_lat(type, par)

    % output info
    foldername = make_savedir_proc(type, par);
    load(sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s/%s/grid.mat', type, par.outname, par.clim));

    par.lat_interp = 'native';

    lat = par.lat;
    for l = {'lo'}; land=l{1};
        ta_list.(land) = nan(length(par.model_list), length(par.lat), 12, length(par.pa));
        ta_mmm.(land) = nan(length(par.lat), 12, length(par.pa));
        ta_std.(land) = nan(length(par.lat), 12, length(par.pa));
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
        prefix_proc = make_prefix_proc(type_in, par);
        grid0 = load(sprintf('%s/grid.mat', prefix));
        ta0 = load(sprintf('%s/ta_pl_mon_lat.mat', prefix_proc));

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
        ta_mmm.(land) = squeeze(nanmean(ta_list.(land),1));
        if strcmp(type, 'rea')
            ta_std.(land) = squeeze(range(ta_list.(land),1));
        else
            ta_std.(land) = squeeze(nanstd(ta_list.(land),1));
        end
    end

    ta = ta_mmm;

    printname = [foldername 'ta_pl_mon_lat'];
    if ~exist(foldername, 'dir')
        mkdir(foldername)
    end
    save(printname, 'ta', 'ta_std', 'lat', '-v7.3');

end
