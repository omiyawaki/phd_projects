function plot_dthedpa_zon_select(type, par)
% temperature profiles at selected locations and month
    make_dirs(type, par)

    prefix = make_prefix(type, par);
    prefix_proc = make_prefix_proc(type, par);
    plotdir = make_plotdir(type, par);
    load(sprintf('%s/grid.mat', prefix));
    load(sprintf('%s/dthedpa_mon_lat.mat', prefix_proc));

    lat_pole = 85;
    lat_mid = 45;

    % for l = {'lo', 'l', 'o'}; land = l{1};
    for l = {'lo'}; land = l{1};
        if strcmp(land, 'lo'); land_text = 'Land + Ocean';
        elseif strcmp(land, 'l'); land_text = 'Land';
        elseif strcmp(land, 'o'); land_text = 'Ocean';
        end
        for m = [1 7]; month = m(1);
            if month==1; mon_str = 'January';
            elseif month==7; mon_str = 'July'; end;

            dthedpa_mon.(land) = squeeze(dthedpa.(land)(:,month,:));

            dthedpa_sp = interp1(lat, dthedpa_mon.(land), -lat_pole); % sounding at -lat_pole S
            dthedpa_np = interp1(lat, dthedpa_mon.(land), lat_pole); % sounding at lat_pole N
            dthedpa_smid = interp1(lat, dthedpa_mon.(land), -lat_mid); % sounding at -lat_mid S
            dthedpa_nmid = interp1(lat, dthedpa_mon.(land), lat_mid); % sounding at lat_mid N
            dthedpa_eq = interp1(lat, dthedpa_mon.(land), 0); % sounding at equator

            % NH MID and HIGH ONLY
            figure(); clf; hold all; box on;
            if m == 1
                h_np = plot(dthedpa_np, grid.dim3.si, 'color', par.blue);
                h_nmid = plot(dthedpa_nmid, grid.dim3.si, 'color', 0.25*[1 1 1]);
            elseif m == 7
                h_np = plot(dthedpa_np, grid.dim3.si, 'color', 0.25*[1 1 1]);
                h_nmid = plot(dthedpa_nmid, grid.dim3.si, 'color', par.orange);
            end
            xlabel('$\frac{\partial \theta_e}{\partial p}$ (K hPa$^{-1}$)'); ylabel('$\sigma$ (unitless)');
            make_title_type_mon(type, mon_str, par);
            legend([h_np h_nmid], sprintf('%g N', lat_pole), sprintf('%g N', lat_mid), 'location', 'northeast');
            axis('tight');
            set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_sq)
            set(gca, 'fontsize', par.fs, 'xlim', 0.2*[-1 1], 'ydir', 'reverse', 'ytick', 1e-3*[0:100:1000], 'ylim', 1e-3*[0 1000], 'xminortick', 'on')
            print(sprintf('%s/thetaeq_zon_sel/%s/%g/dthedpa_nh_only', plotdir, land, month), '-dpng', '-r300');
            close;

            % SH MID and HIGH ONLY
            figure(); clf; hold all; box on;
            h_sp = plot(dthedpa_sp, grid.dim3.si, 'color', par.blue);
            h_smid = plot(dthedpa_smid, grid.dim3.si, 'color', 0.25*[1 1 1]);
            xlabel('$\frac{\partial \theta_e}{\partial p}$ (K hPa$^{-1}$)'); ylabel('$\sigma$ (unitless)');
            make_title_type_mon(type, mon_str, par);
            legend([h_sp h_smid], sprintf('%g S', lat_pole), sprintf('%g S', lat_mid), 'location', 'northeast');
            axis('tight');
            set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_sq)
            set(gca, 'fontsize', par.fs, 'xlim', 0.2*[-1 1], 'ydir', 'reverse', 'ytick', 1e-3*[0:100:1000], 'ylim', 1e-3*[0 1000], 'xminortick', 'on')
            print(sprintf('%s/thetaeq_zon_sel/%s/%g/dthedpa_sh_only', plotdir, land, month), '-dpng', '-r300');
            close;

        end
    end
end
