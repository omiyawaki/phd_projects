function proc_ga_malr_diff_midlatitude_line(type, par)
    make_dirs(type, par)

    prefix = make_prefix(type, par);
    prefix_proc = make_prefix_proc(type, par);
    plotdir = make_plotdir(type, par);
    foldername = make_savedir_si_bl(type, par);

    load(sprintf('%s/grid.mat', prefix)); % read grid data
    % load(sprintf('%s/sftlf.mat', prefix)); % read land fraction data
    load(sprintf('%sga_malr_diff_si_mon_lat_%g.mat', foldername, par.si_up)); % load lat x mon RCAE data
    % landdata = load('/project2/tas1/miyawaki/matlab/landmask/land_mask.mat');
    % par.land = landdata.land_mask; par.landlat = landdata.landlat; par.landlon = landdata.landlon;

    % sftlf = nanmean(sftlf, 1); % zonal average
    % sftlf = repmat(sftlf', [1 12]); % expand land fraction data to time

    lat_bound_list = [-10 10];
    lat_center = 50;

    for lb = 1:length(lat_bound_list); par.lat_bound = lat_bound_list(lb);
    
        [lat, clat, clat_mon, par] = make_midlatitude_lat(lat_center, par);

        savename = sprintf('%sga_malr_diff_midlatitude_lat_%g_to_%g_%g.mat', foldername, par.lat_center-par.lat_bound, par.lat_center+par.lat_bound, par.si_up);
                
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
            [mesh_lat, mesh_mon] = meshgrid(1:12, lat);

            ga_frac = ga_malr_diff.(land);
            ga_frac = interp1(grid.dim3.lat, ga_frac, lat);
            ga_frac_lat.(land) = nansum(ga_frac.*clat_mon)/nansum(clat);

            clear ga_frac

        end % for land

        save(savename, 'ga_frac_lat', '-v7.3')

    end % lat bound

end % function
