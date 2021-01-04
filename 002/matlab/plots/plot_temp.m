function plot_temp(type, par)

    % load data
    [~, ~, ~, lat, par] = load_flux(type, par);
    if strcmp(type, 'era5') | strcmp(type, 'erai') | strcmp(type, 'era5c')
        load(sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/grid.mat', type));
        load(sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/eps_%g_ga_%g/ta_si.mat', type, par.lat_interp, par.ep, par.ga));
        load(sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/eps_%g_ga_%g/ma_si.mat', type, par.lat_interp, par.ep, par.ga));
        plev = grid.dim3.plev/100; % convert Pa to hPa
    elseif strcmp(type, 'gcm')
        load(sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s/%s/grid.mat', type, par.model, par.gcm.clim));
        load(sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s/%s/eps_%g_ga_%g/ta_si.mat', type, par.model, par.gcm.clim, par.lat_interp, par.ep, par.ga));
        load(sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s/%s/eps_%g_ga_%g/ma_si.mat', type, par.model, par.gcm.clim, par.lat_interp, par.ep, par.ga));
        plev = grid.dim3.plev/100; % convert Pa to hPa
    elseif strcmp(type, 'echam')
        load(sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s/grid.mat', type, par.echam.clim));
        load(sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s/eps_%g_ga_%g/ta_si.mat', type, par.echam.clim, par.lat_interp, par.ep, par.ga));
        load(sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s/eps_%g_ga_%g/ma_si.mat', type, par.echam.clim, par.lat_interp, par.ep, par.ga));
        plev = grid.dim3.plev/100; % convert Pa to hPa
    end

    make_dirs_ep(type, par)

    for f = {'mse'}; fw = f{1};
        for c = fieldnames(ta_si.rce.tp.(fw))'; crit = c{1};
            % for t = {'ann', 'djf', 'jja', 'mam', 'son'}; time = t{1};
            for t = {'ann'}; time = t{1};
                % for l = {'lo', 'l', 'o'}; land = l{1};
                for l = {'lo'}; land = l{1};
                    if strcmp(land, 'lo'); land_text = 'Land + Ocean';
                    elseif strcmp(land, 'l'); land_text = 'Land';
                    elseif strcmp(land, 'o'); land_text = 'Ocean';
                    end

                % RCE and RAE separated into NH and SH
                    figure(); clf; hold all;
                    h_rce_tp = plot(ta_si.rce.tp.(fw).(crit).(land).(time), grid.dim3.si, 'color', par.maroon);
                    h_rce_nh = plot(ta_si.rce.nh.(fw).(crit).(land).(time), grid.dim3.si, '-', 'color', par.orange);
                    h_rce_sh = plot(ta_si.rce.sh.(fw).(crit).(land).(time), grid.dim3.si, '--', 'color', par.orange);
                    h_rae_nh = plot(ta_si.rae.nh.(fw).(crit).(land).(time), grid.dim3.si, '-', 'color', par.blue);
                    h_rae_sh = plot(ta_si.rae.sh.(fw).(crit).(land).(time), grid.dim3.si, '--', 'color', par.blue);
                    xlabel('T (K)'); ylabel('$\sigma$ (unitless)');
                    title(sprintf('%s, %s', upper(time), land_text));
                    legend([h_rce_tp h_rce_nh h_rce_sh h_rae_nh h_rae_sh], 'Tropical RCE', 'NH ML RCE', 'SH ML RCE', 'NH RAE', 'SH RAE', 'location', 'northeast');
                    axis('tight');
                    set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_sq)
                    set(gca, 'fontsize', par.fs, 'ydir', 'reverse', 'yscale', 'linear', 'ytick', 1e-3*[10 20 50 100 200 300 400:200:1000], 'ylim', 1e-3*[50 1000], 'xminortick', 'on')
                    hline(0, '-k');
                    print(sprintf('%s/eps_%g_ga_%g/%s/%s/%s/%s/temp/rcae_all', par.plotdir, par.ep, par.ga, fw, crit, land, time), '-dpng', '-r300');
                    close;

                % Difference between RCE temperature profile and moist adiabat
                    figure(); clf; hold all;
                    line([0 0], [100 1000], 'linewidth', 0.5, 'color', 'k');
                    h_rce_all = plot(ta_si.rce.all.(fw).(crit).(land).(time)-ma_si.rce.all.(fw).(crit).(land).(time).ta, grid.dim3.si, 'color', par.orange);
                    xlabel('$T-T_m$ (K)'); ylabel('$\sigma$ (unitless)');
                    title(sprintf('RCE, %s, %s', upper(time), land_text));
                    % legend([h_rce_tp, h_rce_nh, h_rce_sh], 'Tropics', 'NH ML', 'SH ML', 'location', 'southeast');
                    axis('tight');
                    set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_sq)
                    set(gca, 'fontsize', par.fs, 'ydir', 'reverse', 'yscale', 'linear', 'ytick', 1e-3*[10 20 50 100 200 300 400:200:1000], 'ylim', 1e-3*[200 1000], 'xtick', [-5:5:40], 'xminortick', 'on')
                    hline(0, '-k');
                    print(sprintf('%s/eps_%g_ga_%g/%s/%s/%s/%s/temp/rce_all_diff', par.plotdir, par.ep, par.ga, fw, crit, land, time), '-dpng', '-r300');
                    close;

                % Difference between RCE temperature profile and moist adiabat
                    figure(); clf; hold all;
                    line([0 0], [100 1000], 'linewidth', 0.5, 'color', 'k');
                    h_rce_tp = plot(ta_si.rce.tp.(fw).(crit).(land).(time)-ma_si.rce.tp.(fw).(crit).(land).(time).ta, grid.dim3.si, 'color', par.maroon);
                    h_rce_nh = plot(ta_si.rce.nh.(fw).(crit).(land).(time)-ma_si.rce.nh.(fw).(crit).(land).(time).ta, grid.dim3.si, 'color', par.orange);
                    h_rce_sh = plot(ta_si.rce.sh.(fw).(crit).(land).(time)-ma_si.rce.sh.(fw).(crit).(land).(time).ta, grid.dim3.si, '--', 'color', par.orange);
                    xlabel('$T-T_m$ (K)'); ylabel('$\sigma$ (unitless)');
                    title(sprintf('Tropical RCE, %s, %s', upper(time), land_text));
                    % legend([h_rce_tp, h_rce_nh, h_rce_sh], 'Tropics', 'NH ML', 'SH ML', 'location', 'southeast');
                    axis('tight');
                    set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_sq)
                    set(gca, 'fontsize', par.fs, 'ydir', 'reverse', 'yscale', 'linear', 'ytick', 1e-3*[10 20 50 100 200 300 400:200:1000], 'ylim', 1e-3*[200 1000], 'xtick', [-5:5:40], 'xminortick', 'on')
                    hline(0, '-k');
                    print(sprintf('%s/eps_%g_ga_%g/%s/%s/%s/%s/temp/rce_diff', par.plotdir, par.ep, par.ga, fw, crit, land, time), '-dpng', '-r300');
                    close;

                % All RAE compared with moist adiabat
                    figure(); clf; hold all;
                    h_rae = plot(ta_si.rae.all.(fw).(crit).(land).(time), grid.dim3.si, 'color', par.blue);
                    % h_rae_ma_si = plot(ma_si.rae.all.(fw).(crit).(land).(time).ta, grid.dim3.si, ':', 'color', par.orange);
                    xlabel('T (K)'); ylabel('$\sigma$ (unitless)');
                    title(sprintf('RAE, %s, %s', upper(time), land_text));
                    % if strcmp(type, 'era5') | strcmp(type, 'erai')
                    %     legend([h_rae, h_rae_ma_si], upper(type), 'Moist adiabat', 'location', 'southwest');
                    % elseif strcmp(type, 'gcm')
                    %     legend([h_rae, h_rae_ma_si], par.model, 'Moist adiabat', 'location', 'southwest');
                    % end
                    axis('tight');
                    set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_sq)
                    set(gca, 'fontsize', par.fs, 'ydir', 'reverse', 'yscale', 'linear', 'ytick', 1e-3*[10 20 50 100 200 300 400:200:1000], 'ylim', 1e-3*[50 1000], 'xminortick', 'on')
                    hline(0, '-k');
                    print(sprintf('%s/eps_%g_ga_%g/%s/%s/%s/%s/temp/rae_all', par.plotdir, par.ep, par.ga, fw, crit, land, time), '-dpng', '-r300');
                    close;

                % All RCE compared with moist adiabat
                    figure(); clf; hold all;
                    h_rce = plot(ta_si.rce.all.(fw).(crit).(land).(time), grid.dim3.si, 'color', par.orange);
                    h_rce_ma_si = plot(ma_si.rce.all.(fw).(crit).(land).(time).ta, grid.dim3.si, ':', 'color', par.orange);
                    xlabel('T (K)'); ylabel('$\sigma$ (unitless)');
                    title(sprintf('RCE, %s, %s', upper(time), land_text));
                    if strcmp(type, 'era5') | strcmp(type, 'erai') | strcmp(type, 'era5c')
                        legend([h_rce, h_rce_ma_si], upper(type), 'Moist adiabat', 'location', 'southwest');
                    elseif strcmp(type, 'gcm')
                        legend([h_rce, h_rce_ma_si], par.model, 'Moist adiabat', 'location', 'southwest');
                    end
                    axis('tight');
                    set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_sq)
                    set(gca, 'fontsize', par.fs, 'xlim', [200 inf], 'ydir', 'reverse', 'yscale', 'linear', 'ytick', 1e-3*[10 20 50 100 200 300 400:200:1000], 'ylim', 1e-3*[50 1000], 'xminortick', 'on')
                    hline(0, '-k');
                    print(sprintf('%s/eps_%g_ga_%g/%s/%s/%s/%s/temp/rce_all', par.plotdir, par.ep, par.ga, fw, crit, land, time), '-dpng', '-r300');
                    close;

                % Tropical RCE compared with moist adiabat
                    figure(); clf; hold all;
                    h_rce = plot(ta_si.rce.tp.(fw).(crit).(land).(time), grid.dim3.si, 'color', par.maroon);
                    h_rce_ma_si = plot(ma_si.rce.tp.(fw).(crit).(land).(time).ta, grid.dim3.si, ':', 'color', par.maroon);
                    xlabel('T (K)'); ylabel('$\sigma$ (unitless)');
                    title(sprintf('Tropical RCE, %s, %s', upper(time), land_text));
                    if strcmp(type, 'era5') | strcmp(type, 'erai') | strcmp(type, 'merra2') | strcmp(type, 'era5c')
                        legend([h_rce, h_rce_ma_si], upper(type), 'Moist adiabat', 'location', 'southwest');
                    elseif strcmp(type, 'gcm')
                        legend([h_rce, h_rce_ma_si], par.model, 'Moist adiabat', 'location', 'southwest');
                    end
                    axis('tight');
                    set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_sq)
                    set(gca, 'fontsize', par.fs, 'ydir', 'reverse', 'yscale', 'linear', 'ytick', 1e-3*[10 20 50 100 200 300 400:200:1000], 'ylim', 1e-3*[50 1000], 'xminortick', 'on')
                    hline(0, '-k');
                    print(sprintf('%s/eps_%g_ga_%g/%s/%s/%s/%s/temp/rce_tp', par.plotdir, par.ep, par.ga, fw, crit, land, time), '-dpng', '-r300');
                    close;

                % NH RCE compared with moist adiabat
                    figure(); clf; hold all;
                    h_rce = plot(ta_si.rce.nh.(fw).(crit).(land).(time), grid.dim3.si, 'color', par.orange);
                    h_rce_ma_si = plot(ma_si.rce.nh.(fw).(crit).(land).(time).ta, grid.dim3.si, ':', 'color', par.orange);
                    xlabel('T (K)'); ylabel('$\sigma$ (unitless)');
                    title(sprintf('NH RCE, %s, %s', upper(time), land_text));
                    if strcmp(type, 'era5') | strcmp(type, 'erai') | strcmp(type, 'merra2') | strcmp(type, 'era5c')
                        legend([h_rce, h_rce_ma_si], upper(type), 'Moist adiabat', 'location', 'southwest');
                    elseif strcmp(type, 'gcm')
                        legend([h_rce, h_rce_ma_si], par.model, 'Moist adiabat', 'location', 'southwest');
                    end
                    axis('tight');
                    set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_sq)
                    set(gca, 'fontsize', par.fs, 'ydir', 'reverse', 'yscale', 'linear', 'ytick', 1e-3*[10 20 50 100 200 300 400:200:1000], 'ylim', 1e-3*[50 1000], 'xminortick', 'on')
                    hline(0, '-k');
                    print(sprintf('%s/eps_%g_ga_%g/%s/%s/%s/%s/temp/rce_nh', par.plotdir, par.ep, par.ga, fw, crit, land, time), '-dpng', '-r300');
                    close;

                % SH RCE compared with moist adiabat
                    figure(); clf; hold all;
                    h_rce = plot(ta_si.rce.sh.(fw).(crit).(land).(time), grid.dim3.si, 'color', par.orange);
                    h_rce_ma_si = plot(ma_si.rce.sh.(fw).(crit).(land).(time).ta, grid.dim3.si, ':', 'color', par.orange);
                    xlabel('T (K)'); ylabel('$\sigma$ (unitless)');
                    title(sprintf('SH RCE, %s, %s', upper(time), land_text));
                    if strcmp(type, 'era5') | strcmp(type, 'erai') | strcmp(type, 'merra2') | strcmp(type, 'era5c')
                        legend([h_rce, h_rce_ma_si], upper(type), 'Moist adiabat', 'location', 'southwest');
                    elseif strcmp(type, 'gcm')
                        legend([h_rce, h_rce_ma_si], par.model, 'Moist adiabat', 'location', 'southwest');
                    end
                    axis('tight');
                    set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_sq)
                    set(gca, 'fontsize', par.fs, 'ydir', 'reverse', 'yscale', 'linear', 'ytick', 1e-3*[10 20 50 100 200 300 400:200:1000], 'ylim', 1e-3*[50 1000], 'xminortick', 'on')
                    hline(0, '-k');
                    print(sprintf('%s/eps_%g_ga_%g/%s/%s/%s/%s/temp/rce_sh', par.plotdir, par.ep, par.ga, fw, crit, land, time), '-dpng', '-r300');
                    close;

                % NH RCAE compared with moist adiabat
                    figure(); clf; hold all;
                    h_rcae = plot(ta_si.rcae.nh.(fw).(crit).(land).(time), grid.dim3.si, 'color', 0.25*[1 1 1]);
                    h_rcae_ma_si = plot(ma_si.rcae.nh.(fw).(crit).(land).(time).ta, grid.dim3.si, ':', 'color', 0.25*[1 1 1]);
                    xlabel('T (K)'); ylabel('$\sigma$ (unitless)');
                    title(sprintf('NH RCAE, %s, %s', upper(time), land_text));
                    if strcmp(type, 'era5') | strcmp(type, 'erai') | strcmp(type, 'merra2') | strcmp(type, 'era5c')
                        legend([h_rcae, h_rcae_ma_si], upper(type), 'Moist adiabat', 'location', 'southwest');
                    elseif strcmp(type, 'gcm')
                        legend([h_rcae, h_rcae_ma_si], par.model, 'Moist adiabat', 'location', 'southwest');
                    end
                    axis('tight');
                    set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_sq)
                    set(gca, 'fontsize', par.fs, 'ydir', 'reverse', 'yscale', 'linear', 'ytick', 1e-3*[10 20 50 100 200 300 400:200:1000], 'ylim', 1e-3*[50 1000], 'xminortick', 'on')
                    hline(0, '-k');
                    print(sprintf('%s/eps_%g_ga_%g/%s/%s/%s/%s/temp/rcae_nh', par.plotdir, par.ep, par.ga, fw, crit, land, time), '-dpng', '-r300');
                    close;

                % SH RCAE compared with moist adiabat
                    figure(); clf; hold all;
                    h_rcae = plot(ta_si.rcae.sh.(fw).(crit).(land).(time), grid.dim3.si, 'color', 0.25*[1 1 1]);
                    h_rcae_ma_si = plot(ma_si.rcae.sh.(fw).(crit).(land).(time).ta, grid.dim3.si, ':', 'color', 0.25*[1 1 1]);
                    xlabel('T (K)'); ylabel('$\sigma$ (unitless)');
                    title(sprintf('SH RCAE, %s, %s', upper(time), land_text));
                    if strcmp(type, 'era5') | strcmp(type, 'erai') | strcmp(type, 'merra2') | strcmp(type, 'era5c')
                        legend([h_rcae, h_rcae_ma_si], upper(type), 'Moist adiabat', 'location', 'southwest');
                    elseif strcmp(type, 'gcm')
                        legend([h_rcae, h_rcae_ma_si], par.model, 'Moist adiabat', 'location', 'southwest');
                    end
                    axis('tight');
                    set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_sq)
                    set(gca, 'fontsize', par.fs, 'ydir', 'reverse', 'yscale', 'linear', 'ytick', 1e-3*[10 20 50 100 200 300 400:200:1000], 'ylim', 1e-3*[50 1000], 'xminortick', 'on')
                    hline(0, '-k');
                    print(sprintf('%s/eps_%g_ga_%g/%s/%s/%s/%s/temp/rcae_sh', par.plotdir, par.ep, par.ga, fw, crit, land, time), '-dpng', '-r300');
                    close;

                % NH RAE compared with moist adiabat
                    figure(); clf; hold all;
                    h_rae = plot(ta_si.rae.nh.(fw).(crit).(land).(time), grid.dim3.si, 'color', par.blue);
                    xlabel('T (K)'); ylabel('$\sigma$ (unitless)');
                    title(sprintf('NH RAE, %s, %s', upper(time), land_text));
                    axis('tight');
                    set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_sq)
                    set(gca, 'fontsize', par.fs, 'ydir', 'reverse', 'yscale', 'linear', 'ytick', 1e-3*[10 20 50 100 200 300 400:200:1000], 'ylim', 1e-3*[50 1000], 'xminortick', 'on')
                    hline(0, '-k');
                    print(sprintf('%s/eps_%g_ga_%g/%s/%s/%s/%s/temp/rae_nh', par.plotdir, par.ep, par.ga, fw, crit, land, time), '-dpng', '-r300');
                    close;

                % SH RAE compared with moist adiabat
                    figure(); clf; hold all;
                    h_rae = plot(ta_si.rae.sh.(fw).(crit).(land).(time), grid.dim3.si, 'color', par.blue);
                    xlabel('T (K)'); ylabel('$\sigma$ (unitless)');
                    title(sprintf('SH RAE, %s, %s', upper(time), land_text));
                    axis('tight');
                    set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_sq)
                    set(gca, 'fontsize', par.fs, 'ydir', 'reverse', 'yscale', 'linear', 'ytick', 1e-3*[10 20 50 100 200 300 400:200:1000], 'ylim', 1e-3*[50 1000], 'xminortick', 'on')
                    hline(0, '-k');
                    print(sprintf('%s/eps_%g_ga_%g/%s/%s/%s/%s/temp/rae_sh', par.plotdir, par.ep, par.ga, fw, crit, land, time), '-dpng', '-r300');
                    close;

                % ALL NH compared with moist adiabat
                    figure(); clf; hold all;
                    h_rce = plot(ta_si.rce.nh.(fw).(crit).(land).(time), grid.dim3.si, 'color', par.orange);
                    h_rce_ma_si = plot(ma_si.rce.nh.(fw).(crit).(land).(time).ta, grid.dim3.si, ':', 'color', par.orange);
                    h_rcae = plot(ta_si.rcae.sh.(fw).(crit).(land).(time), grid.dim3.si, 'color', 0.25*[1 1 1]);
                    h_rcae_ma_si = plot(ma_si.rcae.sh.(fw).(crit).(land).(time).ta, grid.dim3.si, ':', 'color', 0.25*[1 1 1]);
                    h_rae = plot(ta_si.rae.nh.(fw).(crit).(land).(time), grid.dim3.si, 'color', par.blue);
                    xlabel('T (K)'); ylabel('$\sigma$ (unitless)');
                    title(sprintf('NH, %s', upper(time)));
                    legend([h_rce, h_rcae, h_rae], 'RCE', 'RCAE', 'RAE', 'location', 'southwest');
                    axis('tight');
                    set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_sq)
                    set(gca, 'fontsize', par.fs, 'xlim', [200 inf], 'ydir', 'reverse', 'yscale', 'linear', 'ytick', 1e-3*[10 20 50 100 200 300 400:200:1000], 'ylim', 1e-3*[50 1000], 'xminortick', 'on')
                    hline(0, '-k');
                    print(sprintf('%s/eps_%g_ga_%g/%s/%s/%s/%s/temp/all_nh', par.plotdir, par.ep, par.ga, fw, crit, land, time), '-dpng', '-r300');
                    close;

                % ALL SH compared with moist adiabat
                    figure(); clf; hold all;
                    h_rce = plot(ta_si.rce.sh.(fw).(crit).(land).(time), grid.dim3.si, 'color', par.orange);
                    h_rce_ma_si = plot(ma_si.rce.sh.(fw).(crit).(land).(time).ta, grid.dim3.si, ':', 'color', par.orange);
                    h_rcae = plot(ta_si.rcae.sh.(fw).(crit).(land).(time), grid.dim3.si, 'color', 0.25*[1 1 1]);
                    h_rcae_ma_si = plot(ma_si.rcae.sh.(fw).(crit).(land).(time).ta, grid.dim3.si, ':', 'color', 0.25*[1 1 1]);
                    h_rae = plot(ta_si.rae.sh.(fw).(crit).(land).(time), grid.dim3.si, 'color', par.blue);
                    xlabel('T (K)'); ylabel('$\sigma$ (unitless)');
                    title(sprintf('SH, %s', upper(time)));
                    legend([h_rce, h_rcae, h_rae], 'RCE', 'RCAE', 'RAE', 'location', 'southwest');
                    axis('tight');
                    set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_sq)
                    set(gca, 'fontsize', par.fs, 'xlim', [200 inf], 'ydir', 'reverse', 'yscale', 'linear', 'ytick', 1e-3*[10 20 50 100 200 300 400:200:1000], 'ylim', 1e-3*[50 1000], 'xminortick', 'on')
                    hline(0, '-k');
                    print(sprintf('%s/eps_%g_ga_%g/%s/%s/%s/%s/temp/all_sh', par.plotdir, par.ep, par.ga, fw, crit, land, time), '-dpng', '-r300');
                    close;

                end % land
                % land and ocean comp
                    % figure(); clf; hold all;
                    % h_rce_l = plot(ta_si.rce.nh.(fw).(crit).l.(time)-ma_si.rce.nh.(fw).(crit).l.(time).ta, grid.dim3.si, 'color', par.orange);
                    % h_rce_o = plot(ta_si.rce.nh.(fw).(crit).o.(time)-ma_si.rce.nh.(fw).(crit).o.(time).ta, grid.dim3.si, ':', 'color', par.orange);
                    % xlabel('$T-T_m$ (K)'); ylabel('$\sigma$ (unitless)');
                    % title(sprintf('Tropical RCE, %s, %s', upper(time), land_text));
                    % legend([h_rce_l, h_rce_o], 'Land', 'Ocean', 'location', 'southeast');
                    % axis('tight');
                    % set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_sq)
                    % set(gca, 'fontsize', par.fs, 'ydir', 'reverse', 'yscale', 'linear', 'ytick', 1e-3*[10 20 50 100 200 300 400:200:1000], 'ylim', 1e-3*[100 1000], 'xminortick', 'on')
                    % hline(0, '-k');
                    % print(sprintf('%s/eps_%g_ga_%g/%s/%s/%s/%s/temp/rce_lo_diff', par.plotdir, par.ep, par.ga, fw, crit, 'lo', time), '-dpng', '-r300');
                    % close;
            end % time avg
        end % RCE/RAE definition
    end % MSE/DSE framework
    % Legend for RCE and RAE separated into NH and SH
    figure(); clf; hold all;
    h_rce_tp = plot(ta_si.rce.tp.(fw).(crit).(land).(time), plev, 'color', par.maroon);
    h_rce_nh = plot(ta_si.rce.nh.(fw).(crit).(land).(time), plev, '-', 'color', par.orange);
    h_rce_sh = plot(ta_si.rce.sh.(fw).(crit).(land).(time), plev, '--', 'color', par.orange);
    h_rae_nh = plot(ta_si.rae.nh.(fw).(crit).(land).(time), plev, '-', 'color', par.blue);
    h_rae_sh = plot(ta_si.rae.sh.(fw).(crit).(land).(time), plev, '--', 'color', par.blue);
    xlabel('T (K)'); ylabel('$\sigma$ (unitless)');
    legend([h_rce_tp h_rce_nh h_rce_sh h_rae_nh h_rae_sh], 'Tropical RCE ($\pm 30^\circ$)', 'NH ML RCE ($>+30^\circ$)', 'SH ML RCE ($<-30^\circ$)', 'NH RAE', 'SH RAE', 'location', 'eastoutside');
    axis('tight');
    set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_wide)
    set(gca, 'fontsize', par.fs, 'ydir', 'reverse', 'yscale', 'linear', 'ytick', [10 20 50 100 200 300 400:200:1000], 'ylim', [50 1000], 'xminortick', 'on')
    hline(0, '-k');
    print(sprintf('%s/legends/temp', par.plotdir), '-dpng', '-r300');
    close;
    % Legend for moist adiabat comparisons
    figure(); clf; hold all;
    h_rce_tp = plot(ta_si.rce.tp.(fw).(crit).(land).(time), plev, '-k');
    h_rce_tp_ma_si = plot(ma_si.rce.tp.(fw).(crit).(land).(time).ta, plev, ':k');
    xlabel('T (K)'); ylabel('$\sigma$ (unitless)');
    if strcmp(type, 'era5') | strcmp(type, 'erai') | strcmp(type, 'era5c')
        legend([h_rce_tp h_rce_tp_ma_si], upper(type), 'Moist adiabat', 'location', 'eastoutside');
    elseif strcmp(type, 'gcm')
        legend([h_rce_tp h_rce_tp_ma_si], par.model, 'Moist adiabat', 'location', 'eastoutside');
    end
    title(upper(sprintf('%s', time)));
    axis('tight');
    set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_sq)
    set(gca, 'fontsize', par.fs, 'ydir', 'reverse', 'yscale', 'linear', 'ytick', [10 20 50 100 200 300 400:200:1000], 'ylim', [50 1000], 'xminortick', 'on')
    hline(0, '-k');
    print(sprintf('%s/legends/temp_ma_si', par.plotdir), '-dpng', '-r300');
    close;
end
