function plot_temp_zon_select(type, par)
% temperature profiles at selected locations and month
    make_dirs(type, par)

    prefix = make_prefix(type, par);
    prefix_proc = make_prefix_proc(type, par);
    plotdir = make_plotdir(type, par);
    subplotdir = 'ga_zon_sel';
    
    load(sprintf('%s/grid.mat', prefix));
    load(sprintf('%s/dtdzsi.mat', prefix));
    load(sprintf('%s/malrsi.mat', prefix));

    dtdzsi = squeeze(nanmean(dtdzsi, 1));
    dtmdzsi = squeeze(nanmean(dtmdzsi, 1));

    lat_pole = 85;
    lat_mid = 45;

    for m = [1 4 6 7 8 9 10]; month = m(1);
        if month==1; mon_str = 'January';
        elseif month==4; mon_str = 'April';
        elseif month==6; mon_str = 'June';
        elseif month==7; mon_str = 'July';
        elseif month==8; mon_str = 'August';
        elseif month==9; mon_str = 'September';
        elseif month==10; mon_str = 'October'; end;

        dtdzsi_mon = squeeze(dtdzsi(:,:,month));
        dtmdzsi_mon = squeeze(dtmdzsi(:,:, month));

        dtdzsi_sp(:,m) = interp1(grid.dim3.lat, dtdzsi_mon, -lat_pole); % sounding at -lat_pole S
        dtdzsi_np(:,m) = interp1(grid.dim3.lat, dtdzsi_mon, lat_pole); % sounding at lat_pole N
        dtdzsi_smid(:,m) = interp1(grid.dim3.lat, dtdzsi_mon, -lat_mid); % sounding at -lat_mid S
        dtdzsi_nmid(:,m) = interp1(grid.dim3.lat, dtdzsi_mon, lat_mid); % sounding at lat_mid N
        dtdzsi_eq(:,m) = interp1(grid.dim3.lat, dtdzsi_mon, 0); % sounding at equator

        dtmdzsi_sp(:,m) = interp1(grid.dim3.lat, dtmdzsi_mon, -lat_pole); % sounding at -lat_pole S
        dtmdzsi_np(:,m) = interp1(grid.dim3.lat, dtmdzsi_mon, lat_pole); % sounding at lat_pole N
        dtmdzsi_smid(:,m) = interp1(grid.dim3.lat, dtmdzsi_mon, -lat_mid); % sounding at -lat_mid S
        dtmdzsi_nmid(:,m) = interp1(grid.dim3.lat, dtmdzsi_mon, lat_mid); % sounding at lat_mid N
        dtmdzsi_eq(:,m) = interp1(grid.dim3.lat, dtmdzsi_mon, 0); % sounding at equator

        % NH HIGH ONLY
        figure(); clf; hold all; box on;
        if m == 1 | m == 8 | m== 9
            h_np = plot(dtdzsi_np(:,m), grid.dim3.si, 'color', par.blue);
            plot(dtmdzsi_np(:,m), grid.dim3.si, ':', 'color', par.blue);
        elseif m==6 | m == 7
            h_np = plot(dtdzsi_np(:,m), grid.dim3.si, 'color', 0.25*[1 1 1]);
            plot(dtmdzsi_np(:,m), grid.dim3.si, ':', 'color', 0.25*[1 1 1]);
        end
        xlabel('$\Gamma$ (K/km)'); ylabel('$\sigma$ (unitless)');
        make_title_type_mon(type, mon_str, par);
        axis('tight');
        set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_sq)
        set(gca, 'fontsize', par.fs, 'xlim', [-10 50], 'ydir', 'reverse', 'ytick', 1e-3*[0:100:1000], 'ylim', 1e-3*[200 1000], 'xminortick', 'on')
        folder = sprintf('%s/%s/%02d/', plotdir, subplotdir, month);
        if ~exist(folder,'dir');
            mkdir(folder)
        end
        print(sprintf('%snp', folder), '-dpng', '-r300');
        close;

        % NH MID ONLY
        figure(); clf; hold all; box on;
        if m == 1 | m==8 | m==9
            h_nmid = plot(dtdzsi_nmid(:,m), grid.dim3.si, 'color', 0.25*[1 1 1]);
            plot(dtmdzsi_nmid(:,m), grid.dim3.si, ':', 'color', 0.25*[1 1 1]);
        elseif m==6 | m == 7
            h_nmid = plot(dtdzsi_nmid(:,m), grid.dim3.si, 'color', par.orange);
            plot(dtmdzsi_nmid(:,m), grid.dim3.si, ':', 'color', par.orange);
        end
        xlabel('$\Gamma$ (K/km)'); ylabel('$\sigma$ (unitless)');
        make_title_type_mon_lat(type, mon_str, lat_mid, par);
        axis('tight');
        set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_sq)
        set(gca, 'fontsize', par.fs, 'xlim', [-10 50], 'ydir', 'reverse', 'ytick', 1e-3*[0:100:1000], 'ylim', 1e-3*[200 1000], 'xminortick', 'on')
        print(sprintf('%snhmid', folder), '-dpng', '-r300');
        close;

        % SH HIGH ONLY
        figure(); clf; hold all; box on;
        h_sp = plot(dtdzsi_sp(:,m), grid.dim3.si, 'color', par.blue);
        plot(dtmdzsi_sp(:,m), grid.dim3.si, ':', 'color', par.blue);
        xlabel('$\Gamma$ (K/km)'); ylabel('$\sigma$ (unitless)');
        make_title_type_mon(type, mon_str, par);
        axis('tight');
        set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_sq)
        set(gca, 'fontsize', par.fs, 'xlim', [-10 50], 'ydir', 'reverse', 'ytick', 1e-3*[0:100:1000], 'ylim', 1e-3*[200 1000], 'xminortick', 'on')
        print(sprintf('%ssp', folder), '-dpng', '-r300');
        close;

        % SH MID ONLY
        figure(); clf; hold all; box on;
        h_smid = plot(dtdzsi_smid(:,m), grid.dim3.si, 'color', 0.25*[1 1 1]);
        plot(dtmdzsi_smid(:,m), grid.dim3.si, ':', 'color', 0.25*[1 1 1]);
        xlabel('$\Gamma$ (K/km)'); ylabel('$\sigma$ (unitless)');
        make_title_type_mon_lat(type, mon_str, -lat_mid, par);
        axis('tight');
        set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_sq)
        set(gca, 'fontsize', par.fs, 'xlim', [-10 50], 'ydir', 'reverse', 'ytick', 1e-3*[0:100:1000], 'ylim', 1e-3*[200 1000], 'xminortick', 'on')
        print(sprintf('%sshmid', folder), '-dpng', '-r300');
        close;

        % EQUATOR ONLY
        figure(); clf; hold all; box on;
        h_eq = plot(dtdzsi_eq(:,m), grid.dim3.si, 'color', par.orange);
        plot(dtmdzsi_eq(:,m), grid.dim3.si, ':', 'color', par.orange);
        xlabel('$\Gamma$ (K/km)'); ylabel('$\sigma$ (unitless)');
        make_title_type_mon_lat(type, mon_str, 0, par);
        axis('tight');
        set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_sq)
        set(gca, 'fontsize', par.fs, 'xlim', [-10 50], 'ydir', 'reverse', 'ytick', 1e-3*[0:100:1000], 'ylim', 1e-3*[200 1000], 'xminortick', 'on')
        print(sprintf('%seq', folder), '-dpng', '-r300');
        close;

    end

end
