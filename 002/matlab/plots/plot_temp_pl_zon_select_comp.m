function plot_temp_zon_select_comp(type, par)
% temperature profiles at selected locations and month
    make_dirs(type, par)

    prefix = make_prefix(type, par);
    prefix_proc = make_prefix_proc(type, par);
    plotdir = make_plotdir(type, par);
    
    load(sprintf('%s/grid.mat', prefix));
    load(sprintf('%s/ta_pl_mon_lat.mat', prefix_proc));
    if par.ma
        load(sprintf('%s/ma_mon_lat.mat', prefix_proc));
    end

    par.lat_interp = 'native';
    prefix_era5c = make_prefix('era5c', par);
    prefix_merra = make_prefix('merra2', par);
    prefix_jra55 = make_prefix('jra55', par);

    tmp = load(sprintf('%s/grid.mat', prefix_era5c)); grid_era5c = tmp.grid; clear tmp;
    tmp = load(sprintf('%s/grid.mat', prefix_merra)); grid_merra = tmp.grid; clear tmp;
    tmp = load(sprintf('%s/grid.mat', prefix_jra55)); grid_jra55 = tmp.grid; clear tmp;

    prefix_proc_era5c = make_prefix_proc('era5c', par);
    prefix_proc_merra = make_prefix_proc('merra2', par);
    prefix_proc_jra55 = make_prefix_proc('jra55', par);

    tmp = load(sprintf('%s/ta_pl_mon_lat.mat', prefix_proc_era5c)); ta_era5c = tmp.ta; lat_era5c = tmp.lat; clear tmp;
    tmp = load(sprintf('%s/ta_pl_mon_lat.mat', prefix_proc_merra)); ta_merra = tmp.ta; lat_merra = tmp.lat; clear tmp;
    tmp = load(sprintf('%s/ta_pl_mon_lat.mat', prefix_proc_jra55)); ta_jra55 = tmp.ta; lat_jra55 = tmp.lat; clear tmp;

    % load HARA data
    load('/project2/tas1/miyawaki/projects/002/data/raw/hara/HARA_avgsnd_80N-90N.mat')

    mon_list = [1 4 6 7 10];

    % area averaged profiles
    lat_bound_list = [80 -80];
    pole_lat = 83;

    for lb = 1:length(lat_bound_list); par.lat_bound = lat_bound_list(lb);

        [alat, clat, clat_mon, par] =  make_polar_lat(par, pole_lat);

        % for l = {'lo', 'l', 'o'}; land = l{1};
        for l = {'lo'}; land = l{1};
            if strcmp(land, 'lo'); land_text = 'L+O';
            elseif strcmp(land, 'l'); land_text = 'L';
            elseif strcmp(land, 'o'); land_text = 'O';
            end

        for m = mon_list; month = m(1);
            if month==1; mon_str = 'January';
            elseif month==4; mon_str = 'April';
            elseif month==6; mon_str = 'June';
            elseif month==7; mon_str = 'July';
            elseif month==10; mon_str = 'October'; end;

                ta_mon.(land) = squeeze(ta.(land)(:,month,:));
                ta_mon_era5c.(land) = squeeze(ta_era5c.(land)(:,month,:));
                ta_mon_merra.(land) = squeeze(ta_merra.(land)(:,month,:));
                ta_mon_jra55.(land) = squeeze(ta_jra55.(land)(:,month,:));

                clat_lev = repmat(clat', [1 length(grid.dim3.plev)]);
                clat_lev_era5c = repmat(clat', [1 length(grid_era5c.dim3.plev)]);
                clat_lev_merra = repmat(clat', [1 length(grid_merra.dim3.plev)]);
                clat_lev_jra55 = repmat(clat', [1 length(grid_jra55.dim3.plev)]);

                ta_area(:,m) = nansum(clat_lev.* interp1(grid.dim3.lat, ta_mon.(land), alat) )./nansum(clat); % area averaged sounding
                ta_area_era5c(:,m) = nansum(clat_lev_era5c.* interp1(grid_era5c.dim3.lat, ta_mon_era5c.(land), alat) )./nansum(clat); % area averaged sounding
                ta_area_merra(:,m) = nansum(clat_lev_merra.* interp1(grid_merra.dim3.lat, ta_mon_merra.(land), alat) )./nansum(clat); % area averaged sounding
                ta_area_jra55(:,m) = nansum(clat_lev_jra55.* interp1(grid_jra55.dim3.lat, ta_mon_jra55.(land), alat) )./nansum(clat); % area averaged sounding

                % NH HIGH ONLY
                figure(); clf; hold all; box on;
                cmip5 = plot(ta_area(:,m), 1e-2*par.pa, 'color', 'k');
                era5c = plot(ta_area_era5c(:,m), 1e-2*grid_era5c.dim3.plev, 'color', par.blue);
                merra = plot(ta_area_merra(:,m), 1e-2*grid_merra.dim3.plev, 'color', par.orange);
                jra55 = plot(ta_area_jra55(:,m), 1e-2*grid_jra55.dim3.plev, 'color', par.green);
                hara = plot(hdat.tmon(:,m), hdat.pstd, '--', 'color', par.maroon);
                xlabel('T (K)'); ylabel('$p$ (hPa)');
                title(sprintf('$\\phi=%g^\\circ$ to $%g^\\circ$, %s', par.lat_bound, par.lat_pole, mon_str));
                legend([cmip5 era5c merra jra55 hara], 'CMIP5', 'ERA5', 'MERRA2', 'JRA55', 'HARA', 'location', 'eastoutside');
                axis('tight');
                set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_wide)
                set(gca, 'fontsize', par.fs, 'xlim', [205 inf], 'ydir', 'reverse', 'ytick', [0:100:1000], 'ylim', [200 1000], 'xminortick', 'on')
                print(sprintf('%s/temp_zon_sel_comp/%s/%g/area_pl_%g', plotdir, land, month, par.lat_bound), '-dpng', '-r300');
                close;
                end % month

        end % land

    end % lat bound
    
    % latitudinal band profiles
    lat_pole = 88;
    lat_mid = 45;

    % for l = {'lo', 'l', 'o'}; land = l{1};
    for l = {'lo'}; land = l{1};
        if strcmp(land, 'lo'); land_text = 'Land + Ocean';
        elseif strcmp(land, 'l'); land_text = 'Land';
        elseif strcmp(land, 'o'); land_text = 'Ocean';
        end
        for m = mon_list; month = m(1);
            if month==1; mon_str = 'January';
            elseif month==4; mon_str = 'April';
            elseif month==6; mon_str = 'June';
            elseif month==7; mon_str = 'July';
            elseif month==10; mon_str = 'October'; end;

            ta_mon.(land) = squeeze(ta.(land)(:,month,:));
            ta_mon_era5c.(land) = squeeze(ta_era5c.(land)(:,month,:));
            ta_mon_merra.(land) = squeeze(ta_merra.(land)(:,month,:));
            ta_mon_jra55.(land) = squeeze(ta_jra55.(land)(:,month,:));

            if par.ma
                masi_mon.(land) = squeeze(masi.(land)(:,month,:));

                % remove moist adiabat data below initialization level
                if ~strcmp(par.ma_init, 'surf')
                    masi_mon.(land)(:, grid.dim3.si>par.ma_init) = nan;
                end
            end

            ta_sp(:,m) = interp1(lat, ta_mon.(land), -lat_pole); % sounding at -lat_pole S
            ta_np(:,m) = interp1(lat, ta_mon.(land), lat_pole); % sounding at lat_pole N
            ta_smid(:,m) = interp1(lat, ta_mon.(land), -lat_mid); % sounding at -lat_mid S
            ta_nmid(:,m) = interp1(lat, ta_mon.(land), lat_mid); % sounding at lat_mid N
            ta_eq(:,m) = interp1(lat, ta_mon.(land), 0); % sounding at equator

            ta_np_era5c(:,m) = interp1(lat_era5c, ta_mon_era5c.(land), lat_pole); % sounding at lat_pole N
            ta_np_merra(:,m) = interp1(lat_merra, ta_mon_merra.(land), lat_pole); % sounding at lat_pole N
            ta_np_jra55(:,m) = interp1(lat_jra55, ta_mon_jra55.(land), lat_pole); % sounding at lat_pole N

            ta_sp_era5c(:,m) = interp1(lat_era5c, ta_mon_era5c.(land), -lat_pole); % sounding at lat_pole N
            ta_sp_merra(:,m) = interp1(lat_merra, ta_mon_merra.(land), -lat_pole); % sounding at lat_pole N
            ta_sp_jra55(:,m) = interp1(lat_jra55, ta_mon_jra55.(land), -lat_pole); % sounding at lat_pole N

            if par.ma
                masi_sp(:,m) = interp1(lat, masi_mon.(land), -lat_pole); % sounding at -lat_pole S
                masi_np(:,m) = interp1(lat, masi_mon.(land), lat_pole); % sounding at lat_pole N
                masi_smid(:,m) = interp1(lat, masi_mon.(land), -lat_mid); % sounding at -lat_mid S
                masi_nmid(:,m) = interp1(lat, masi_mon.(land), lat_mid); % sounding at lat_mid N
                masi_eq(:,m) = interp1(lat, masi_mon.(land), 0); % sounding at equator
            end

            % NH HIGH ONLY
            figure(); clf; hold all; box on;
            cmip5 = plot(ta_np(:,m), 1e-2*par.pa, 'color', 'k');
            era5c = plot(ta_np_era5c(:,m), 1e-2*grid_era5c.dim3.plev, 'color', par.blue);
            merra = plot(ta_np_merra(:,m), 1e-2*grid_merra.dim3.plev, 'color', par.orange);
            jra55 = plot(ta_np_jra55(:,m), 1e-2*grid_jra55.dim3.plev, 'color', par.green);
            hara = plot(hdat.tmon(:,m), hdat.pstd, '--', 'color', par.maroon);
            xlabel('T (K)'); ylabel('$p$ (hPa)');
            title(sprintf('$\\phi=%g^\\circ$, %s', lat_pole, mon_str));
            legend([cmip5 era5c merra jra55 hara], 'CMIP5', 'ERA5', 'MERRA2', 'JRA55', 'HARA', 'location', 'eastoutside');
            axis('tight');
            set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_wide)
            set(gca, 'fontsize', par.fs, 'xlim', [205 inf], 'ydir', 'reverse', 'ytick', [0:100:1000], 'ylim', [200 1000], 'xminortick', 'on')
            print(sprintf('%s/temp_zon_sel_comp/%s/%g/np_pl', plotdir, land, month), '-dpng', '-r300');
            close;

            % SH HIGH ONLY
            figure(); clf; hold all; box on;
            cmip5 = plot(ta_sp(:,m), 1e-2*par.pa, 'color', 'k');
            era5c = plot(ta_sp_era5c(:,m), 1e-2*grid_era5c.dim3.plev, 'color', par.blue);
            merra = plot(ta_sp_merra(:,m), 1e-2*grid_merra.dim3.plev, 'color', par.orange);
            jra55 = plot(ta_sp_jra55(:,m), 1e-2*grid_jra55.dim3.plev, 'color', par.green);
            hara = plot(hdat.tmon(:,m), hdat.pstd, '--', 'color', par.maroon);
            xlabel('T (K)'); ylabel('$p$ (hPa)');
            title(sprintf('$\\phi=%g^\\circ$, %s', lat_pole, mon_str));
            legend([cmip5 era5c merra jra55 hara], 'CMIP5', 'ERA5', 'MERRA2', 'JRA55', 'HARA', 'location', 'eastoutside');
            axis('tight');
            set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_wide)
            set(gca, 'fontsize', par.fs, 'xlim', [205 inf], 'ydir', 'reverse', 'ytick', [0:100:1000], 'ylim', [200 1000], 'xminortick', 'on')
            print(sprintf('%s/temp_zon_sel_comp/%s/%g/sp_pl', plotdir, land, month), '-dpng', '-r300');
            close;

        end % months

    end % land

end
