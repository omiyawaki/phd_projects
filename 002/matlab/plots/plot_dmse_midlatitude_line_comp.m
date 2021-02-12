function plot_dmse_midlatitude_line_comp(type, par)
    make_dirs(type, par)

    prefix = make_prefix(type, par);
    prefix_proc = make_prefix_proc(type, par);
    plotdir = make_plotdir(type, par);

    load(sprintf('%s/grid.mat', prefix)); % read grid data

    prefix_era5c = make_prefix('era5c', par);
    prefix_merra = make_prefix('merra2', par);
    prefix_jra55 = make_prefix('jra55', par);

    tmp = load(sprintf('%s/grid.mat', prefix_era5c)); grid_era5c = tmp.grid; clear tmp; % read grid data
    tmp = load(sprintf('%s/grid.mat', prefix_merra)); grid_merra = tmp.grid; clear tmp; % read grid data
    tmp = load(sprintf('%s/grid.mat', prefix_jra55)); grid_jra55 = tmp.grid; clear tmp; % read grid data

    par.lat_interp = 'native';
    prefix_proc_era5c = make_prefix_proc('era5c', par);
    prefix_proc_merra = make_prefix_proc('merra2', par);
    prefix_proc_jra55 = make_prefix_proc('jra55', par);

    lat_bound_list = [-10 10];
    lat_center = 50;

    for lb = 1:length(lat_bound_list); par.lat_bound = lat_bound_list(lb);

        [lat, clat, clat_mon, par] = make_midlatitude_lat(lat_center, par);

        tmp = load(sprintf('%s/dmse_midlatitude_lat_%g_to_%g', prefix_proc, par.lat_center-par.lat_bound, par.lat_center+par.lat_bound)); dmse = tmp.dmse; dmse_std = tmp.dmse_std; clear tmp;
        tmp = load(sprintf('%s/dmse_midlatitude_lat_%g_to_%g', prefix_proc_era5c, par.lat_center-par.lat_bound, par.lat_center+par.lat_bound)); dmse_era5c = tmp.dmse; clear tmp;
        tmp = load(sprintf('%s/dmse_midlatitude_lat_%g_to_%g', prefix_proc_merra, par.lat_center-par.lat_bound, par.lat_center+par.lat_bound)); dmse_merra = tmp.dmse; clear tmp;
        tmp = load(sprintf('%s/dmse_midlatitude_lat_%g_to_%g', prefix_proc_jra55, par.lat_center-par.lat_bound, par.lat_center+par.lat_bound)); dmse_jra55 = tmp.dmse; clear tmp;

        % for l = {'lo', 'l', 'o'}; land = l{1};
        for l = {'lo'}; land = l{1};
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
            [mesh_lat, mesh_mon] = meshgrid(1:12, grid.dim2.lat);

            f_vec = assign_fw(type, par);
            for f = f_vec; fw = f{1};
            
                par.folder = sprintf('%s/dmse-comp/%s/%s/0_midlatitude_lat_%g_to_%g', plotdir, fw, land, par.lat_center-par.lat_bound, par.lat_center+par.lat_bound);
                if ~exist(par.folder, 'dir'); mkdir(par.folder); end;

                plot_dmse_comp(dmse.ra_lat.(land).(fw), dmse_std.ra_lat.(land).(fw), dmse_era5c.ra_lat.(land).(fw), dmse_merra.ra_lat.(land).(fw), dmse_jra55.ra_lat.(land).(fw), '$R_a$ (Wm$^{-2}$)', 'ra', "", type, fw, par);

                plot_dmse_comp(dmse.lh_lat.(land).(fw), dmse_std.lh_lat.(land).(fw), dmse_era5c.lh_lat.(land).(fw), dmse_merra.lh_lat.(land).(fw), dmse_jra55.lh_lat.(land).(fw), 'LH (Wm$^{-2}$)', 'lh', "", type, fw, par);

                plot_dmse_comp(dmse.sh_lat.(land).(fw), dmse_std.sh_lat.(land).(fw), dmse_era5c.sh_lat.(land).(fw), dmse_merra.sh_lat.(land).(fw), dmse_jra55.sh_lat.(land).(fw), 'SH (Wm$^{-2}$)', 'sh', "", type, fw, par);

            end % for mse dse
        end % for land
    end % lat bound

end % for function
