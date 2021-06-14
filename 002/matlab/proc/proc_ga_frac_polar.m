function proc_ga_frac_polar(type, par)
    make_dirs(type, par)

    prefix = make_prefix(type, par);
    prefix_proc = make_prefix_proc(type, par);
    plotdir = make_plotdir(type, par);

    load(sprintf('%s/grid.mat', prefix)); % read grid data
    % load(sprintf('%s/sftlf.mat', prefix)); % read land fraction data
    load(sprintf('%s/ga_frac_mon_lat.mat', prefix_proc)); % load lat x mon RCAE data
    % load(sprintf('%s/%s/flux_t.mat', prefix_proc, par.lat_interp)); % load lat x lon RCAE data
    % landdata = load('/project2/tas1/miyawaki/matlab/landmask/land_mask.mat');
    % par.land = landdata.land_mask; par.landlat = landdata.landlat; par.landlon = landdata.landlon;

    % sftlf = nanmean(sftlf, 1); % zonal average
    % sftlf = repmat(sftlf', [1 12]); % expand land fraction data to time
    
    lat_bound_list = [-80 80];

    for lb = 1:length(lat_bound_list); par.lat_bound = lat_bound_list(lb);
    
        [lat, clat, clat_mon, par] = make_polar_lat(par);

        savename = sprintf('%s/ga_frac_poleward_of_lat_%g', prefix_proc, par.lat_bound);
                
        for l = par.land_list; land = l{1};
            if strcmp(land, 'lo'); land_text = 'L+O';
            elseif strcmp(land, 'l'); land_text = 'L';
            elseif strcmp(land, 'o'); land_text = 'O';
            end
            if strcmp(type, 'echam');
                if strcmp(par.echam.clim, '20170908')
                    land_text = 'Snowball';
                else
                    land_text = par.echam.(par.echam.clim);
                end
            end;

            ga_frac_lat.(land) = interp1(grid.dim3.lat, ga_frac.(land), lat);
            clat_mon_lev = repmat(clat_mon, [1, 1, length(grid.dim3.si)]);
            ga_frac_lat.(land) = squeeze(nansum(ga_frac_lat.(land) .* clat_mon_lev, 1)/nansum(clat));

        end % for land

        save(savename, 'ga_frac_lat', '-v7.3')

    end % lat bound

end % function
