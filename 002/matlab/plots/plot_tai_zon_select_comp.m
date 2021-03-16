function plot_tai_zon_select_comp(type, par)
% temperature profiles at selected locations and month
    make_dirs(type, par)

    prefix = make_prefix(type, par);
    prefix_proc = make_prefix_proc(type, par);
    plotdir = make_plotdir(type, par);
    
    load(sprintf('%s/grid.mat', prefix));
    load(sprintf('%s/tai_mon_lat.mat', prefix_proc));
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

    tmp = load(sprintf('%s/tai_mon_lat.mat', prefix_proc_era5c)); tai_era5c = tmp.tai; lat_era5c = tmp.lat; clear tmp;
    tmp = load(sprintf('%s/tai_mon_lat.mat', prefix_proc_merra)); tai_merra = tmp.tai; lat_merra = tmp.lat; clear tmp;
    tmp = load(sprintf('%s/tai_mon_lat.mat', prefix_proc_jra55)); tai_jra55 = tmp.tai; lat_jra55 = tmp.lat; clear tmp;

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

                tai_mon.(land) = squeeze(tai.(land)(:,month,:));
                tai_mon_era5c.(land) = squeeze(tai_era5c.(land)(:,month,:));
                tai_mon_merra.(land) = squeeze(tai_merra.(land)(:,month,:));
                tai_mon_jra55.(land) = squeeze(tai_jra55.(land)(:,month,:));

                % clat_lev = repmat(clat', [1 length(par.pa)]);
                % clat_lev_era5c = repmat(clat', [1 length(par.pa)]);
                % clat_lev_merra = repmat(clat', [1 length(par.pa)]);
                % clat_lev_jra55 = repmat(clat', [1 length(par.pa)]);

                taii = interp1(grid.dim3.lat, tai_mon.(land), alat);
                taii_era5c = interp1(grid_era5c.dim3.lat, tai_mon_era5c.(land), alat);
                taii_merra = interp1(grid_merra.dim3.lat, tai_mon_merra.(land), alat);
                taii_jra55 = interp1(grid_jra55.dim3.lat, tai_mon_jra55.(land), alat);

                nanf = ones(size(taii));
                nanf_era5c = ones(size(taii_era5c));
                nanf_merra = ones(size(taii_merra));
                nanf_jra55 = ones(size(taii_jra55));

                nanf(isnan(taii)) = nan;
                nanf_era5c(isnan(taii_era5c)) = nan;
                nanf_merra(isnan(taii_merra)) = nan;
                nanf_jra55(isnan(taii_jra55)) = nan;

                tai_area = nan([length(par.pa) 12]);

                for ilev = 1:length(par.pa); 
                    tai_area(ilev,m) = nansum(clat'.* taii(:,ilev))/nansum(clat'.*nanf(:,ilev)); % area averaged sounding
                    tai_area_era5c(ilev,m) = nansum(clat'.* taii_era5c(:,ilev))/nansum(clat'.*nanf_era5c(:,ilev)); % area averaged sounding
                    tai_area_merra(ilev,m) = nansum(clat'.* taii_merra(:,ilev))/nansum(clat'.*nanf_merra(:,ilev)); % area averaged sounding
                    tai_area_jra55(ilev,m) = nansum(clat'.* taii_jra55(:,ilev))/nansum(clat'.*nanf_jra55(:,ilev)); % area averaged sounding
                end

                % NH HIGH ONLY
                figure(); clf; hold all; box on;
                cmip5 = plot(tai_area(:,m), 1e-2*par.pa, 'color', 'k');
                era5c = plot(tai_area_era5c(:,m), 1e-2*par.pa, 'color', par.blue);
                merra = plot(tai_area_merra(:,m), 1e-2*par.pa, 'color', par.orange);
                jra55 = plot(tai_area_jra55(:,m), 1e-2*par.pa, 'color', par.green);
                hara = plot(hdat.tmon(:,m), hdat.pstd, '--', 'color', par.maroon);
                xlabel('T (K)'); ylabel('$p$ (hPa)');
                title(sprintf('$\\phi=%g^\\circ$ to $%g^\\circ$, %s', par.lat_bound, par.lat_pole, mon_str));
                legend([cmip5 era5c merra jra55 hara], 'CMIP5', 'ERA5', 'MERRA2', 'JRA55', 'HARA', 'location', 'eastoutside');
                axis('tight');
                set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_wide)
                set(gca, 'fontsize', par.fs, 'xlim', [205 inf], 'ydir', 'reverse', 'ytick', [0:100:1000], 'ylim', [200 1000], 'xminortick', 'on')
                print(sprintf('%s/temp_zon_sel_comp/%s/%g/area_tai_%g', plotdir, land, month, par.lat_bound), '-dpng', '-r300');
                close;
                end % month

        end % land

    end % lat bound
    
    % latitudinal band profiles
    lat_pole = 82;
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

            tai_mon.(land) = squeeze(tai.(land)(:,month,:));
            tai_mon_era5c.(land) = squeeze(tai_era5c.(land)(:,month,:));
            tai_mon_merra.(land) = squeeze(tai_merra.(land)(:,month,:));
            tai_mon_jra55.(land) = squeeze(tai_jra55.(land)(:,month,:));

            if par.ma
                masi_mon.(land) = squeeze(masi.(land)(:,month,:));

                % remove moist adiabat data below initialization level
                if ~strcmp(par.ma_init, 'surf')
                    masi_mon.(land)(:, grid.dim3.si>par.ma_init) = nan;
                end
            end

            tai_sp(:,m) = interp1(lat, tai_mon.(land), -lat_pole); % sounding at -lat_pole S
            tai_np(:,m) = interp1(lat, tai_mon.(land), lat_pole); % sounding at lat_pole N
            tai_smid(:,m) = interp1(lat, tai_mon.(land), -lat_mid); % sounding at -lat_mid S
            tai_nmid(:,m) = interp1(lat, tai_mon.(land), lat_mid); % sounding at lat_mid N
            tai_eq(:,m) = interp1(lat, tai_mon.(land), 0); % sounding at equator

            tai_np_era5c(:,m) = interp1(lat_era5c, tai_mon_era5c.(land), lat_pole); % sounding at lat_pole N
            tai_np_merra(:,m) = interp1(lat_merra, tai_mon_merra.(land), lat_pole); % sounding at lat_pole N
            tai_np_jra55(:,m) = interp1(lat_jra55, tai_mon_jra55.(land), lat_pole); % sounding at lat_pole N

            tai_sp_era5c(:,m) = interp1(lat_era5c, tai_mon_era5c.(land), -lat_pole); % sounding at lat_pole N
            tai_sp_merra(:,m) = interp1(lat_merra, tai_mon_merra.(land), -lat_pole); % sounding at lat_pole N
            tai_sp_jra55(:,m) = interp1(lat_jra55, tai_mon_jra55.(land), -lat_pole); % sounding at lat_pole N

            if par.ma
                masi_sp(:,m) = interp1(lat, masi_mon.(land), -lat_pole); % sounding at -lat_pole S
                masi_np(:,m) = interp1(lat, masi_mon.(land), lat_pole); % sounding at lat_pole N
                masi_smid(:,m) = interp1(lat, masi_mon.(land), -lat_mid); % sounding at -lat_mid S
                masi_nmid(:,m) = interp1(lat, masi_mon.(land), lat_mid); % sounding at lat_mid N
                masi_eq(:,m) = interp1(lat, masi_mon.(land), 0); % sounding at equator
            end

            % NH HIGH ONLY
            figure(); clf; hold all; box on;
            cmip5 = plot(tai_np(:,m), 1e-2*par.pa, 'color', 'k');
            era5c = plot(tai_np_era5c(:,m), 1e-2*par.pa, 'color', par.blue);
            merra = plot(tai_np_merra(:,m), 1e-2*par.pa, 'color', par.orange);
            jra55 = plot(tai_np_jra55(:,m), 1e-2*par.pa, 'color', par.green);
            hara = plot(hdat.tmon(:,m), hdat.pstd, '--', 'color', par.maroon);
            xlabel('T (K)'); ylabel('$p$ (hPa)');
            title(sprintf('$\\phi=%g^\\circ$, %s', lat_pole, mon_str));
            legend([cmip5 era5c merra jra55 hara], 'CMIP5', 'ERA5', 'MERRA2', 'JRA55', 'HARA', 'location', 'eastoutside');
            axis('tight');
            set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_wide)
            set(gca, 'fontsize', par.fs, 'xlim', [205 inf], 'ydir', 'reverse', 'ytick', [0:100:1000], 'ylim', [200 1000], 'xminortick', 'on')
            print(sprintf('%s/temp_zon_sel_comp/%s/%g/np_tai', plotdir, land, month), '-dpng', '-r300');
            close;

            % SH HIGH ONLY
            figure(); clf; hold all; box on;
            cmip5 = plot(tai_sp(:,m), 1e-2*par.pa, 'color', 'k');
            era5c = plot(tai_sp_era5c(:,m), 1e-2*par.pa, 'color', par.blue);
            merra = plot(tai_sp_merra(:,m), 1e-2*par.pa, 'color', par.orange);
            jra55 = plot(tai_sp_jra55(:,m), 1e-2*par.pa, 'color', par.green);
            hara = plot(hdat.tmon(:,m), hdat.pstd, '--', 'color', par.maroon);
            xlabel('T (K)'); ylabel('$p$ (hPa)');
            title(sprintf('$\\phi=%g^\\circ$, %s', lat_pole, mon_str));
            legend([cmip5 era5c merra jra55 hara], 'CMIP5', 'ERA5', 'MERRA2', 'JRA55', 'HARA', 'location', 'eastoutside');
            axis('tight');
            set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_wide)
            set(gca, 'fontsize', par.fs, 'xlim', [205 inf], 'ydir', 'reverse', 'ytick', [0:100:1000], 'ylim', [200 1000], 'xminortick', 'on')
            print(sprintf('%s/temp_zon_sel_comp/%s/%g/sp_tai', plotdir, land, month), '-dpng', '-r300');
            close;

        end % months

    end % land

end
