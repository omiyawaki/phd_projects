function plot_temp_zon_select_comp_ml(type, par)
% temperature profiles at selected locations and month
    make_dirs(type, par)

    prefix = make_prefix(type, par);
    prefix_proc = make_prefix_proc(type, par);
    plotdir = make_plotdir(type, par);
    
    load(sprintf('%s/grid.mat', prefix));
    load(sprintf('%s/ta_mon_lat.mat', prefix_proc));
    tmp = load(sprintf('%s/taml_mon_lat.mat', prefix_proc)); tamlsi = tmp.tasi; clear tmp;
    if par.ma
        load(sprintf('%s/ma_mon_lat.mat', prefix_proc));
    end

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
            tamlsi_mon.(land) = squeeze(tamlsi.(land)(:,month,:));

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

            tamlsi_sp(:,m) = interp1(lat, tamlsi_mon.(land), -lat_pole); % sounding at -lat_pole S
            tamlsi_np(:,m) = interp1(lat, tamlsi_mon.(land), lat_pole); % sounding at lat_pole N
            tamlsi_smid(:,m) = interp1(lat, tamlsi_mon.(land), -lat_mid); % sounding at -lat_mid S
            tamlsi_nmid(:,m) = interp1(lat, tamlsi_mon.(land), lat_mid); % sounding at lat_mid N
            tamlsi_eq(:,m) = interp1(lat, tamlsi_mon.(land), 0); % sounding at equator

            if par.ma
                masi_sp(:,m) = interp1(lat, masi_mon.(land), -lat_pole); % sounding at -lat_pole S
                masi_np(:,m) = interp1(lat, masi_mon.(land), lat_pole); % sounding at lat_pole N
                masi_smid(:,m) = interp1(lat, masi_mon.(land), -lat_mid); % sounding at -lat_mid S
                masi_nmid(:,m) = interp1(lat, masi_mon.(land), lat_mid); % sounding at lat_mid N
                masi_eq(:,m) = interp1(lat, masi_mon.(land), 0); % sounding at equator
            end

            % NH HIGH ONLY
            figure(); clf; hold all; box on;
            pl = plot(tasi_np(:,m), grid.dim3.si, 'color', 'k');
            ml = plot(tamlsi_np(:,m), grid.dim3.si, 'color', par.blue);
            xlabel('T (K)'); ylabel('$\sigma$ (unitless)');
            make_title_type_mon_lat(type, mon_str, lat_pole, par);
            legend([pl ml], 'pl', 'ml', 'location', 'eastoutside');
            axis('tight');
            set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_wide)
            set(gca, 'fontsize', par.fs, 'xlim', [205 inf], 'ydir', 'reverse', 'ytick', 1e-3*[0:100:1000], 'ylim', 1e-3*[200 1000], 'xminortick', 'on')
            print(sprintf('%s/temp_zon_sel_comp/%s/%g/np_ml', plotdir, land, month), '-dpng', '-r300');
            close;

            % SH HIGH ONLY
            figure(); clf; hold all; box on;
            pl = plot(tasi_sp(:,m), grid.dim3.si, 'color', 'k');
            ml = plot(tamlsi_sp(:,m), grid.dim3.si, 'color', par.blue);
            xlabel('T (K)'); ylabel('$\sigma$ (unitless)');
            make_title_type_mon_lat(type, mon_str, -lat_pole, par);
            legend([pl ml], 'pl', 'ml', 'location', 'eastoutside');
            axis('tight');
            set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_wide)
            set(gca, 'fontsize', par.fs, 'xlim', [205 inf], 'ydir', 'reverse', 'ytick', 1e-3*[0:100:1000], 'ylim', 1e-3*[200 1000], 'xminortick', 'on')
            print(sprintf('%s/temp_zon_sel_comp/%s/%g/sp_ml', plotdir, land, month), '-dpng', '-r300');
            close;

            % EQ ONLY
            figure(); clf; hold all; box on;
            pl = plot(tasi_eq(:,m), grid.dim3.si, 'color', 'k');
            ml = plot(tamlsi_eq(:,m), grid.dim3.si, 'color', par.blue);
            xlabel('T (K)'); ylabel('$\sigma$ (unitless)');
            make_title_type_mon_lat(type, mon_str, 0, par);
            legend([pl ml], 'pl', 'ml', 'location', 'eastoutside');
            axis('tight');
            set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_wide)
            set(gca, 'fontsize', par.fs, 'xlim', [205 inf], 'ydir', 'reverse', 'ytick', 1e-3*[0:100:1000], 'ylim', 1e-3*[200 1000], 'xminortick', 'on')
            print(sprintf('%s/temp_zon_sel_comp/%s/%g/eq_ml', plotdir, land, month), '-dpng', '-r300');
            close;

            % NH MID ONLY
            figure(); clf; hold all; box on;
            pl = plot(tasi_nmid(:,m), grid.dim3.si, 'color', 'k');
            ml = plot(tamlsi_nmid(:,m), grid.dim3.si, 'color', par.blue);
            xlabel('T (K)'); ylabel('$\sigma$ (unitless)');
            make_title_type_mon_lat(type, mon_str, lat_mid, par);
            legend([pl ml], 'pl', 'ml', 'location', 'eastoutside');
            axis('tight');
            set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_wide)
            set(gca, 'fontsize', par.fs, 'xlim', [205 inf], 'ydir', 'reverse', 'ytick', 1e-3*[0:100:1000], 'ylim', 1e-3*[200 1000], 'xminortick', 'on')
            print(sprintf('%s/temp_zon_sel_comp/%s/%g/nmid_ml', plotdir, land, month), '-dpng', '-r300');
            close;

            % NH MID ONLY
            figure(); clf; hold all; box on;
            pl = plot(tasi_smid(:,m), grid.dim3.si, 'color', 'k');
            ml = plot(tamlsi_smid(:,m), grid.dim3.si, 'color', par.blue);
            xlabel('T (K)'); ylabel('$\sigma$ (unitless)');
            make_title_type_mon_lat(type, mon_str, -lat_mid, par);
            legend([pl ml], 'pl', 'ml', 'location', 'eastoutside');
            axis('tight');
            set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_wide)
            set(gca, 'fontsize', par.fs, 'xlim', [205 inf], 'ydir', 'reverse', 'ytick', 1e-3*[0:100:1000], 'ylim', 1e-3*[200 1000], 'xminortick', 'on')
            print(sprintf('%s/temp_zon_sel_comp/%s/%g/smid_ml', plotdir, land, month), '-dpng', '-r300');
            close;

        end

    end
end
