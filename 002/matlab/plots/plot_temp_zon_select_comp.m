function plot_temp_zon_select_comp(type, par)
% temperature profiles at selected locations and month
    make_dirs(type, par)

    prefix = make_prefix(type, par);
    prefix_proc = make_prefix_proc(type, par);
    plotdir = make_plotdir(type, par);
    
    load(sprintf('%s/grid.mat', prefix));
    load(sprintf('%s/ta_mon_lat.mat', prefix_proc));
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

    tmp = load(sprintf('%s/ta_mon_lat.mat', prefix_proc_era5c)); tasi_era5c = tmp.tasi; lat_era5c = tmp.lat; clear tmp;
    tmp = load(sprintf('%s/ta_mon_lat.mat', prefix_proc_merra)); tasi_merra = tmp.tasi; lat_merra = tmp.lat; clear tmp;
    tmp = load(sprintf('%s/ta_mon_lat.mat', prefix_proc_jra55)); tasi_jra55 = tmp.tasi; lat_jra55 = tmp.lat; clear tmp;

    lat_pole = 85;
    lat_mid = 45;

    % for l = {'lo', 'l', 'o'}; land = l{1};
    for l = {'lo'}; land = l{1};
        if strcmp(land, 'lo'); land_text = 'Land + Ocean';
        elseif strcmp(land, 'l'); land_text = 'Land';
        elseif strcmp(land, 'o'); land_text = 'Ocean';
        end
        for m = [1 4 7 10]; month = m(1);
            if month==1; mon_str = 'January';
            elseif month==4; mon_str = 'April';
            elseif month==6; mon_str = 'June';
            elseif month==7; mon_str = 'July';
            elseif month==10; mon_str = 'October'; end;

            tasi_mon.(land) = squeeze(tasi.(land)(:,month,:));
            tasi_mon_era5c.(land) = squeeze(tasi_era5c.(land)(:,month,:));
            tasi_mon_merra.(land) = squeeze(tasi_merra.(land)(:,month,:));
            tasi_mon_jra55.(land) = squeeze(tasi_jra55.(land)(:,month,:));

            if par.ma
                masi_mon.(land) = squeeze(masi.(land)(:,month,:));

                % remove moist adiabat data below initialization level
                if ~strcmp(par.ma_init, 'surf')
                    masi_mon.(land)(:, grid.dim3.si>par.ma_init) = nan;
                end
            end

            tasi_sp(:,m) = interp1(lat, tasi_mon.(land), -lat_pole); % sounding at -lat_pole S
            tasi_np(:,m) = interp1(lat, tasi_mon.(land), lat_pole); % sounding at lat_pole N
            tasi_smid(:,m) = interp1(lat, tasi_mon.(land), -lat_mid); % sounding at -lat_mid S
            tasi_nmid(:,m) = interp1(lat, tasi_mon.(land), lat_mid); % sounding at lat_mid N
            tasi_eq(:,m) = interp1(lat, tasi_mon.(land), 0); % sounding at equator

            tasi_np_era5c(:,m) = interp1(lat_era5c, tasi_mon_era5c.(land), lat_pole); % sounding at lat_pole N
            tasi_np_merra(:,m) = interp1(lat_merra, tasi_mon_merra.(land), lat_pole); % sounding at lat_pole N
            tasi_np_jra55(:,m) = interp1(lat_jra55, tasi_mon_jra55.(land), lat_pole); % sounding at lat_pole N

            if par.ma
                masi_sp(:,m) = interp1(lat, masi_mon.(land), -lat_pole); % sounding at -lat_pole S
                masi_np(:,m) = interp1(lat, masi_mon.(land), lat_pole); % sounding at lat_pole N
                masi_smid(:,m) = interp1(lat, masi_mon.(land), -lat_mid); % sounding at -lat_mid S
                masi_nmid(:,m) = interp1(lat, masi_mon.(land), lat_mid); % sounding at lat_mid N
                masi_eq(:,m) = interp1(lat, masi_mon.(land), 0); % sounding at equator
            end

            % NH HIGH ONLY
            figure(); clf; hold all; box on;
            cmip5 = plot(tasi_np(:,m), grid.dim3.si, 'color', 'k');
            era5c = plot(tasi_np_era5c(:,m), grid_era5c.dim3.si, 'color', par.blue);
            merra = plot(tasi_np_merra(:,m), grid_merra.dim3.si, 'color', par.orange);
            jra55 = plot(tasi_np_jra55(:,m), grid_jra55.dim3.si, 'color', par.green);
            xlabel('T (K)'); ylabel('$\sigma$ (unitless)');
            title(sprintf('$\\phi=%g^\\circ$, %s', lat_pole, mon_str));
            legend([cmip5 era5c merra jra55], 'CMIP5', 'ERA5', 'MERRA2', 'JRA55', 'location', 'eastoutside');
            axis('tight');
            set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_wide)
            set(gca, 'fontsize', par.fs, 'xlim', [205 inf], 'ydir', 'reverse', 'ytick', 1e-3*[0:100:1000], 'ylim', 1e-3*[200 1000], 'xminortick', 'on')
            print(sprintf('%s/temp_zon_sel_comp/%s/%g/np', plotdir, land, month), '-dpng', '-r300');
            close;

        end

    end
end
