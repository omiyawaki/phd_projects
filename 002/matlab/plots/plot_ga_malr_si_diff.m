function plot_ga_malr_si_diff(type, par)
% sigma coordinates
    make_dirs(type, par)

    % load data
    % [~, ~, ~, lat, par] = load_flux(type, par);
    if strcmp(type, 'era5') | strcmp(type, 'erai') | strcmp(type, 'era5c')
        par.plotdir = sprintf('./figures/%s/%s/%s', type, par.(type).yr_span, par.lat_interp);
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s', type, par.(type).yr_span);
    elseif strcmp(type, 'gcm')
        par.plotdir = sprintf('./figures/%s/%s/%s/%s', type, par.model, par.gcm.clim, par.lat_interp);
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s/%s', type, par.model, par.gcm.clim);
    elseif strcmp(type, 'echam')
        par.plotdir = sprintf('./figures/%s/%s/%s', type, par.(type).yr_span, par.lat_interp);
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s', type, par.echam.clim);
    elseif any(strcmp(type, {'echam_ml', 'echam_pl'}))
        par.plotdir = sprintf('./figures/%s/%s/%s', type, par.(type).yr_span, par.lat_interp);
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s', type, par.(type).yr_span);
    end
    load(sprintf('%s/grid.mat', prefix)); % read grid data
    load(sprintf('%s/dtdzsi.mat', prefix)); dtdzzsi = dtdzsi; clear dtdzsi; % read model lapse rate
    load(sprintf('%s/malrsi.mat', prefix)); dtmdzzsi = dtmdzsi; clear dtmdzsi; % read moist adiabatic lapse rate

    dtdz_zt = squeeze(nanmean(nanmean(dtdzzsi, 1), 4));
    dtmdz_zt = squeeze(nanmean(nanmean(dtmdzzsi, 1), 4));
    dtdz_jan = squeeze(nanmean(dtdzzsi(:,:,:,1), 1));
    dtmdz_jan = squeeze(nanmean(dtmdzzsi(:,:,:,1), 1));
    dtdz_jul = squeeze(nanmean(dtdzzsi(:,:,:,7), 1));
    dtmdz_jul = squeeze(nanmean(dtmdzzsi(:,:,:,7), 1));

    diff = (dtmdzzsi - dtdzzsi)./dtmdzzsi * 1e2; % percentage difference of moist adiabatic vs GCM lapse rates

    % LAT x PLEV
    diff_zt = squeeze(nanmean(nanmean(diff, 1), 4)); % take zonal and annual mean
    diff_jan = squeeze(nanmean(diff(:,:,:,1), 1)); % take zonal and annual mean
    diff_jul = squeeze(nanmean(diff(:,:,:,7), 1)); % take zonal and annual mean

    % % vertically average
    % dtdz_ztv = permute(dtdz_zt, [2 1]); % bring height front
    % dtdz_ztv = interp1(par.pa, dtdz_ztv, 1e2*linspace(1000,200,100)); % prepare to average between 1000-200 hPa
    % dtdz_ztv = squeeze(nanmean(dtdz_ztv,1)); % take vertical average
    % dtmdz_ztv = permute(dtmdz_zt, [2 1]); % bring height front
    % dtmdz_ztv = interp1(par.pa, dtmdz_ztv, 1e2*linspace(1000,200,100)); % prepare to average between 1000-200 hPa
    % dtmdz_ztv = squeeze(nanmean(dtmdz_ztv,1)); % take vertical average
    % dtdz_janv = permute(dtdz_jan, [2 1]); % bring height front
    % dtdz_janv = interp1(par.pa, dtdz_janv, 1e2*linspace(1000,200,100)); % prepare to average between 1000-200 hPa
    % dtdz_janv = squeeze(nanmean(dtdz_janv,1)); % take vertical average
    % dtmdz_janv = permute(dtmdz_jan, [2 1]); % bring height front
    % dtmdz_janv = interp1(par.pa, dtmdz_janv, 1e2*linspace(1000,200,100)); % prepare to average between 1000-200 hPa
    % dtmdz_janv = squeeze(nanmean(dtmdz_janv,1)); % take vertical average
    % dtdz_julv = permute(dtdz_jul, [2 1]); % bring height front
    % dtdz_julv = interp1(par.pa, dtdz_julv, 1e2*linspace(1000,200,100)); % prepare to average between 1000-200 hPa
    % dtdz_julv = squeeze(nanmean(dtdz_julv,1)); % take vertical average
    % dtmdz_julv = permute(dtmdz_jul, [2 1]); % bring height front
    % dtmdz_julv = interp1(par.pa, dtmdz_julv, 1e2*linspace(1000,200,100)); % prepare to average between 1000-200 hPa
    % dtmdz_julv = squeeze(nanmean(dtmdz_julv,1)); % take vertical average

    % % lat of climatological and moist adiabatic lapse rate
    % figure(); clf; hold all;
    % plot(grid.dim3.lat(1:6:end), dtdz_ztv(1:6:end), 'xk', 'markersize', 4);
    % plot(grid.dim3.lat(1:6:end), dtmdz_ztv(1:6:end), '^k', 'markersize', 4);
    % xlabel('Latitude (deg)'); ylabel('Lapse rate (K/km)')
    % legend('$\Gamma$', '$\Gamma_m$', 'location', 'southwest');
    % if strcmp(type, 'era5') | strcmp(type, 'erai')
    %     title(sprintf('%s', upper(type)));
    % elseif strcmp(type, 'gcm')
    %     title(sprintf('%s', par.model));
    % end
    % set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos)
    % set(gca, 'xlim', [-10 80], 'xtick', [0:15:75], 'ylim', [0 9], 'ytick', [0:2:8]);
    % print(sprintf('%s/ga_diff/ga_comp_lat', par.plotdir), '-dpng', '-r300');
    % close;

    % figure(); clf; hold all;
    % plot(grid.dim3.lat(1:6:end), dtdz_janv(1:6:end), 'xk', 'markersize', 4);
    % plot(grid.dim3.lat(1:6:end), dtmdz_janv(1:6:end), '^k', 'markersize', 4);
    % xlabel('Latitude (deg)'); ylabel('Lapse rate (K/km)')
    % legend('$\Gamma$', '$\Gamma_m$', 'location', 'southwest');
    % if strcmp(type, 'era5') | strcmp(type, 'erai')
    %     title(sprintf('%s, January', upper(type)));
    % elseif strcmp(type, 'gcm')
    %     title(sprintf('%s, January', par.model));
    % end
    % set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos)
    % set(gca, 'xlim', [-10 80], 'xtick', [0:15:75], 'ylim', [0 9], 'ytick', [0:2:8]);
    % print(sprintf('%s/ga_diff/ga_comp_lat_jan', par.plotdir), '-dpng', '-r300');
    % close;

    % figure(); clf; hold all;
    % plot(grid.dim3.lat(1:6:end), dtdz_julv(1:6:end), 'xk', 'markersize', 4);
    % plot(grid.dim3.lat(1:6:end), dtmdz_julv(1:6:end), '^k', 'markersize', 4);
    % xlabel('Latitude (deg)'); ylabel('Lapse rate (K/km)')
    % legend('$\Gamma$', '$\Gamma_m$', 'location', 'southwest');
    % if strcmp(type, 'era5') | strcmp(type, 'erai')
    %     title(sprintf('%s, July', upper(type)));
    % elseif strcmp(type, 'gcm')
    %     title(sprintf('%s, July', par.model));
    % end
    % set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos)
    % set(gca, 'xlim', [-10 80], 'xtick', [0:15:75], 'ylim', [0 9], 'ytick', [0:2:8]);
    % print(sprintf('%s/ga_diff/ga_comp_lat_jul', par.plotdir), '-dpng', '-r300');
    % close;

    [mesh_si, mesh_lat] = meshgrid(grid.dim3.lat, grid.dim3.si);

    % lat x plev of lapse rate percentage diff
    figure(); clf; hold all;
    cmp = colCog(40);
    colormap(cmp);
    [C, h] = contourf(mesh_si, mesh_lat, diff_zt', [-200 -100 -50 -20 -10 0 10 20 50 100 200]);
    clabel(C, h, [-100 -50 -20 -10 0 10 20 50 100], 'fontsize', 6, 'interpreter', 'latex');
    caxis([-100 100]);
    xlabel('Latitude (deg)'); ylabel('$\sigma$ (unitless)');
    if strcmp(type, 'era5') | strcmp(type, 'erai') | strcmp(type, 'era5c')
        title(sprintf('%s, Annual', upper(type)));
    elseif strcmp(type, 'gcm')
        title(sprintf('%s, Annual', par.model));
    end
    cb = colorbar('limits', [-100 100], 'ytick', [-100:20:100], 'location', 'eastoutside');
    cb.TickLabelInterpreter = 'latex'; cb.Label.Interpreter = 'latex';
    ylabel(cb, '$(\Gamma_m - \Gamma)/\Gamma_m$ (\%)');
    set(gca, 'xlim', [-90 90], 'xtick', [-90:30:90], 'ylim', [0.2 1], 'ytick', [0.2:0.1:1], 'ydir', 'reverse', 'yscale', 'linear', 'yminortick', 'on', 'tickdir', 'out');
    print(sprintf('%s/ga_malr_diff/ga_malr_diff_lat_sigma', par.plotdir), '-dpng', '-r300');
    close;

    % lat x plev of lapse rate percentage diff
    figure(); clf; hold all;
    cmp = colCog(40);
    colormap(cmp);
    [C, h] = contourf(mesh_si, mesh_lat, diff_jul', [-200 -100 -50 -20 -10 0 10 20 50 100 200]);
    clabel(C, h, [-100 -50 -20 -10 0 10 20 50 100], 'fontsize', 6, 'interpreter', 'latex');
    caxis([-100 100]);
    xlabel('Latitude (deg)'); ylabel('$\sigma$ (unitless)');
    if strcmp(type, 'era5') | strcmp(type, 'erai') | strcmp(type, 'era5c')
        title(sprintf('%s, July', upper(type)));
    elseif strcmp(type, 'gcm')
        title(sprintf('%s, July', par.model));
    end
    cb = colorbar('limits', [-100 100], 'ytick', [-100:20:100], 'location', 'eastoutside');
    cb.TickLabelInterpreter = 'latex'; cb.Label.Interpreter = 'latex';
    ylabel(cb, '$(\Gamma_m - \Gamma)/\Gamma_m$ (\%)');
    set(gca, 'xlim', [-90 90], 'xtick', [-90:30:90], 'ylim', [0.2 1], 'ytick', [0.2:0.1:1], 'ydir', 'reverse', 'yscale', 'linear', 'yminortick', 'on', 'tickdir', 'out');
    print(sprintf('%s/ga_malr_diff/ga_malr_diff_lat_sigma_jul', par.plotdir), '-dpng', '-r300');
    close;

    % lat x plev of lapse rate percentage diff
    figure(); clf; hold all;
    cmp = colCog(40);
    colormap(cmp);
    [C, h] = contourf(mesh_si, mesh_lat, diff_jan', [-200 -100 -50 -20 -10 0 10 20 50 100 200]);
    clabel(C, h, [-100 -50 -20 -10 0 10 20 50 100], 'fontsize', 6, 'interpreter', 'latex');
    caxis([-100 100]);
    xlabel('Latitude (deg)'); ylabel('$\sigma$ (unitless)');
    if strcmp(type, 'era5') | strcmp(type, 'erai') | strcmp(type, 'era5c')
        title(sprintf('%s, January', upper(type)));
    elseif strcmp(type, 'gcm')
        title(sprintf('%s, January', par.model));
    end
    cb = colorbar('limits', [-100 100], 'ytick', [-100:20:100], 'location', 'eastoutside');
    cb.TickLabelInterpreter = 'latex'; cb.Label.Interpreter = 'latex';
    ylabel(cb, '$(\Gamma_m - \Gamma)/\Gamma_m$ (\%)');
    set(gca, 'xlim', [-90 90], 'xtick', [-90:30:90], 'ylim', [0.2 1], 'ytick', [0.2:0.1:1], 'ydir', 'reverse', 'yscale', 'linear', 'yminortick', 'on', 'tickdir', 'out');
    print(sprintf('%s/ga_malr_diff/ga_malr_diff_lat_sigma_jan', par.plotdir), '-dpng', '-r300');
    close;

    % lat x sigma of lapse rate percentage diff
    figure(); clf; hold all;
    cmp = colCog(40);
    colormap(cmp);
    [C, h] = contour(mesh_si, mesh_lat, diff_zt', [-200 -100 -50 -20 -10 0 10 20 50 100 200], 'color', 'k');
    clabel(C, h, [-100 -50 -20 -10 0 10 20 50 100], 'fontsize', 6, 'interpreter', 'latex');
    caxis([-100 100]);
    xlabel('Latitude (deg)'); ylabel('$\sigma$ (unitless)');
    if strcmp(type, 'era5') | strcmp(type, 'erai') | strcmp(type, 'era5c')
        title(sprintf('%s, Annual', upper(type)));
    elseif strcmp(type, 'gcm')
        title(sprintf('%s, Annual', par.model));
    end
    % cb = colorbar('limits', [-100 100], 'ytick', [-100:20:100], 'location', 'eastoutside');
    % cb.TickLabelInterpreter = 'latex'; cb.Label.Interpreter = 'latex';
    % ylabel(cb, '$(\Gamma_m - \Gamma)/\Gamma_m$ (\%)');
    set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_sq)
    set(gca, 'xlim', [-10 75], 'xtick', [-10:10:70 75], 'ylim', [0.2 1], 'ytick', [0.2:0.1:1], 'ydir', 'reverse', 'yminortick', 'on', 'tickdir', 'out');
    print(sprintf('%s/ga_malr_diff/ga_malr_diff_lat_sigma_sc', par.plotdir), '-dpng', '-r300');
    close;

    % lat x sigma of lapse rate percentage diff
    figure(); clf; hold all;
    cmp = colCog(40);
    colormap(cmp);
    [C, h] = contour(mesh_si, mesh_lat, diff_jan', [-200 -100 -50 -20 -10 0 10 20 50 100 200], 'color', 'k');
    clabel(C, h, [-100 -50 -20 -10 0 10 20 50 100], 'fontsize', 6, 'interpreter', 'latex');
    caxis([-100 100]);
    xlabel('Latitude (deg)'); ylabel('$\sigma$ (unitless)');
    if strcmp(type, 'era5') | strcmp(type, 'erai') | strcmp(type, 'era5c')
        title(sprintf('%s, January', upper(type)));
    elseif strcmp(type, 'gcm')
        title(sprintf('%s, January', par.model));
    end
    % cb = colorbar('limits', [-100 100], 'ytick', [-100:20:100], 'location', 'eastoutside');
    % cb.TickLabelInterpreter = 'latex'; cb.Label.Interpreter = 'latex';
    % ylabel(cb, '$(\Gamma_m - \Gamma)/\Gamma_m$ (\%)');
    set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_sq)
    set(gca, 'xlim', [-10 75], 'xtick', [-10:10:70 75], 'ylim', [0.2 1], 'ytick', [0.2:0.1:1], 'ydir', 'reverse', 'yminortick', 'on', 'tickdir', 'out');
    print(sprintf('%s/ga_malr_diff/ga_malr_diff_lat_sigma_sc_jan', par.plotdir), '-dpng', '-r300');
    close;

    % lat x sigma of lapse rate percentage diff
    figure(); clf; hold all;
    cmp = colCog(40);
    colormap(cmp);
    [C, h] = contour(mesh_si, mesh_lat, diff_jul', [-100 -50 -20 -10 0 10 20 50 100 200], 'color', 'k');
    clabel(C, h, [-50 -20 -10 0 10 20 50 100], 'fontsize', 6, 'interpreter', 'latex');
    caxis([-100 100]);
    xlabel('Latitude (deg)'); ylabel('$\sigma$ (unitless)');
    if strcmp(type, 'era5') | strcmp(type, 'erai') | strcmp(type, 'era5c')
        title(sprintf('%s, July', upper(type)));
    elseif strcmp(type, 'gcm')
        title(sprintf('%s, July', par.model));
    end
    % cb = colorbar('limits', [-100 100], 'ytick', [-100:20:100], 'location', 'eastoutside');
    % cb.TickLabelInterpreter = 'latex'; cb.Label.Interpreter = 'latex';
    % ylabel(cb, '$(\Gamma_m - \Gamma)/\Gamma_m$ (\%)');
    set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_sq)
    set(gca, 'xlim', [-10 75], 'xtick', [-10:10:70 75], 'ylim', [0.2 1], 'ytick', [0.2:0.1:1], 'ydir', 'reverse', 'yminortick', 'on', 'tickdir', 'out');
    print(sprintf('%s/ga_malr_diff/ga_malr_diff_lat_sigma_sc_jul', par.plotdir), '-dpng', '-r300');
    close;

    % lat x sigma of model lapse rate
    figure(); clf; hold all;
    cmp = colCog(10);
    colormap(cmp);
    [C, h] = contourf(mesh_si, mesh_lat, dtdz_zt', [-20 -10:2:10 20], 'linecolor', 'w');
    clabel(C, h, [-10:2:10], 'fontsize', 6, 'interpreter', 'latex', 'color', 'w');
    caxis([-10 10]);
    xlabel('Latitude (deg)'); ylabel('$\sigma$ (unitless)');
    if strcmp(type, 'era5') | strcmp(type, 'erai') | strcmp(type, 'era5c')
        title(sprintf('%s, Annual', upper(type)));
    elseif strcmp(type, 'gcm')
        title(sprintf('%s, Annual', par.model));
    end
    cb = colorbar('limits', [-10 10], 'ytick', [-10:2:10], 'location', 'eastoutside');
    cb.TickLabelInterpreter = 'latex'; cb.Label.Interpreter = 'latex';
    ylabel(cb, '$\Gamma$ (K/km)');
    set(gca, 'xlim', [-90 90], 'xtick', [-90:30:90], 'ylim', [0.2 1], 'ytick', [0.2:0.1:1], 'ydir', 'reverse', 'yscale', 'linear', 'yminortick', 'on', 'tickdir', 'out');
    print(sprintf('%s/ga_malr_diff/ga_lat_sigma', par.plotdir), '-dpng', '-r300');
    close;

    % lat x sigma of model lapse rate JANUARY
    figure(); clf; hold all;
    cmp = colCog(10);
    colormap(cmp);
    [C, h] = contourf(mesh_si, mesh_lat, dtdz_jan', [-20 -10:2:10 20], 'linecolor', 'w');
    clabel(C, h, [-10:2:10], 'fontsize', 6, 'interpreter', 'latex', 'color', 'w');
    caxis([-10 10]);
    xlabel('Latitude (deg)'); ylabel('$\sigma$ (unitless)');
    if strcmp(type, 'era5') | strcmp(type, 'erai') | strcmp(type, 'era5c')
        title(sprintf('%s, January', upper(type)));
    elseif strcmp(type, 'gcm')
        title(sprintf('%s, January', par.model));
    end
    cb = colorbar('limits', [-10 10], 'ytick', [-10:2:10], 'location', 'eastoutside');
    cb.TickLabelInterpreter = 'latex'; cb.Label.Interpreter = 'latex';
    ylabel(cb, '$\Gamma$ (K/km)');
    set(gca, 'xlim', [-90 90], 'xtick', [-90:30:90], 'ylim', [0.2 1], 'ytick', [0.2:0.1:1], 'ydir', 'reverse', 'yscale', 'linear', 'yminortick', 'on', 'tickdir', 'out');
    print(sprintf('%s/ga_malr_diff/ga_lat_sigma_jan', par.plotdir), '-dpng', '-r300');
    close;

    % lat x sigma of model lapse rate JULY
    figure(); clf; hold all;
    cmp = colCog(10);
    colormap(cmp);
    [C, h] = contourf(mesh_si, mesh_lat, dtdz_jul', [-20 -10:2:10 20], 'linecolor', 'w');
    clabel(C, h, [-10:2:10], 'fontsize', 6, 'interpreter', 'latex', 'color', 'w');
    caxis([-10 10]);
    xlabel('Latitude (deg)'); ylabel('$\sigma$ (unitless)');
    if strcmp(type, 'era5') | strcmp(type, 'erai') | strcmp(type, 'era5c')
        title(sprintf('%s, July', upper(type)));
    elseif strcmp(type, 'gcm')
        title(sprintf('%s, July', par.model));
    end
    cb = colorbar('limits', [-10 10], 'ytick', [-10:2:10], 'location', 'eastoutside');
    cb.TickLabelInterpreter = 'latex'; cb.Label.Interpreter = 'latex';
    ylabel(cb, '$\Gamma$ (K/km)');
    set(gca, 'xlim', [-90 90], 'xtick', [-90:30:90], 'ylim', [0.2 1], 'ytick', [0.2:0.1:1], 'ydir', 'reverse', 'yscale', 'linear', 'yminortick', 'on', 'tickdir', 'out');
    print(sprintf('%s/ga_malr_diff/ga_lat_sigma_jul', par.plotdir), '-dpng', '-r300');
    close;

    % lat x sigma of model lapse rate 10 S to 75 N
    figure(); clf; hold all; box on;
    cmp = colCog(40);
    colormap(cmp);
    [C, h] = contour(mesh_si, mesh_lat, dtdz_zt', [0:9], 'color', 'k');
    clabel(C, h, [0 2 4 5 6 7 8], 'fontsize', 6, 'interpreter', 'latex');
    caxis([-10 10]);
    xlabel('Latitude (deg)'); ylabel('$\sigma$ (unitless)');
    if strcmp(type, 'era5') | strcmp(type, 'erai') | strcmp(type, 'era5c')
        title(sprintf('%s, Annual', upper(type)));
    elseif strcmp(type, 'gcm')
        title(sprintf('%s, Annual', par.model));
    end
    % cb = colorbar('limits', [-10 10], 'ytick', [-10:2:10], 'location', 'eastoutside');
    % cb.TickLabelInterpreter = 'latex'; cb.Label.Interpreter = 'latex';
    % ylabel(cb, '$\Gamma$ (K/km)');
    set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_sq)
    set(gca, 'xlim', [-10 75], 'xtick', [-10:10:70 75], 'ylim', [0.2 1], 'ytick', [0.2:0.1:1], 'ydir', 'reverse', 'yminortick', 'on', 'tickdir', 'out');
    print(sprintf('%s/ga_malr_diff/ga_lat_sigma_sc', par.plotdir), '-dpng', '-r300');
    close;

    % lat x sigma of model lapse rate 10 S to 75 N JANUARY
    figure(); clf; hold all; box on;
    cmp = colCog(40);
    colormap(cmp);
    [C, h] = contour(mesh_si, mesh_lat, dtdz_jan', [-5:9], 'color', 'k');
    clabel(C, h, [-3 -1 2 4 5 6 7 8], 'fontsize', 6, 'interpreter', 'latex');
    caxis([-10 10]);
    xlabel('Latitude (deg)'); ylabel('$\sigma$ (unitless)');
    if strcmp(type, 'era5') | strcmp(type, 'erai') | strcmp(type, 'era5c')
        title(sprintf('%s, January', upper(type)));
    elseif strcmp(type, 'gcm')
        title(sprintf('%s, January', par.model));
    end
    % cb = colorbar('limits', [-10 10], 'ytick', [-10:2:10], 'location', 'eastoutside');
    % cb.TickLabelInterpreter = 'latex'; cb.Label.Interpreter = 'latex';
    % ylabel(cb, '$\Gamma$ (K/km)');
    set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_sq)
    set(gca, 'xlim', [-10 75], 'xtick', [-10:10:70 75], 'ylim', [0.2 1], 'ytick', [0.2:0.1:1], 'ydir', 'reverse', 'yminortick', 'on', 'tickdir', 'out');
    print(sprintf('%s/ga_malr_diff/ga_lat_sigma_sc_jan', par.plotdir), '-dpng', '-r300');
    close;

    % lat x sigma of model lapse rate 10 S to 75 N JULY
    figure(); clf; hold all; box on;
    cmp = colCog(40);
    colormap(cmp);
    [C, h] = contour(mesh_si, mesh_lat, dtdz_jul', [0:9], 'color', 'k');
    clabel(C, h, [0 1 2 3 4 5 6 7 8], 'fontsize', 6, 'interpreter', 'latex');
    caxis([-10 10]);
    xlabel('Latitude (deg)'); ylabel('$\sigma$ (unitless)');
    if strcmp(type, 'era5') | strcmp(type, 'erai') | strcmp(type, 'era5c')
        title(sprintf('%s, July', upper(type)));
    elseif strcmp(type, 'gcm')
        title(sprintf('%s, July', par.model));
    end
    % cb = colorbar('limits', [-10 10], 'ytick', [-10:2:10], 'location', 'eastoutside');
    % cb.TickLabelInterpreter = 'latex'; cb.Label.Interpreter = 'latex';
    % ylabel(cb, '$\Gamma$ (K/km)');
    set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_sq)
    set(gca, 'xlim', [-10 75], 'xtick', [-10:10:70 75], 'ylim', [0.2 1], 'ytick', [0.2:0.1:1], 'ydir', 'reverse', 'yminortick', 'on', 'tickdir', 'out');
    print(sprintf('%s/ga_malr_diff/ga_lat_sigma_sc_jul', par.plotdir), '-dpng', '-r300');
    close;

    % lat x sigma of moist adiabatic lapse rate
    figure(); clf; hold all;
    cmp = colCog(10);
    colormap(cmp);
    [C, h] = contourf(mesh_si, mesh_lat, dtmdz_zt', -10:2:10, 'w');
    clabel(C, h, [-10:2:10], 'fontsize', 6, 'interpreter', 'latex', 'color', 'w');
    caxis([-10 10]);
    xlabel('Latitude (deg)'); ylabel('$\sigma$ (unitless)');
    if strcmp(type, 'era5') | strcmp(type, 'erai') | strcmp(type, 'era5c')
        title(sprintf('%s, Annual', upper(type)));
    elseif strcmp(type, 'gcm')
        title(sprintf('%s, Annual', par.model));
    end
    cb = colorbar('limits', [0 10], 'ytick', [-10:2:10], 'location', 'eastoutside');
    cb.TickLabelInterpreter = 'latex'; cb.Label.Interpreter = 'latex';
    ylabel(cb, '$\Gamma_m$ (K/km)');
    set(gca, 'xlim', [-90 90], 'xtick', [-90:30:90], 'ylim', [0.2 1], 'ytick', [0.2:0.1:1], 'ydir', 'reverse', 'yscale', 'linear', 'yminortick', 'on', 'tickdir', 'out');
    print(sprintf('%s/ga_malr_diff/ga_m_lat_sigma', par.plotdir), '-dpng', '-r300');
    close;

    % lat x sigma of moist adiabatic lapse rate 10 S to 75 N
    figure(); clf; hold all; box on;
    cmp = colCog(10);
    colormap(cmp);
    [C, h] = contour(mesh_si, mesh_lat, dtmdz_zt', [5 6 7 8 9], 'color', 'k');
    clabel(C, h, [5 6 7 8 9], 'fontsize', 6, 'interpreter', 'latex');
    caxis([-10 10]);
    xlabel('Latitude (deg)'); ylabel('$\sigma$ (unitless)');
    if strcmp(type, 'era5') | strcmp(type, 'erai') | strcmp(type, 'era5c')
        title(sprintf('%s, Annual', upper(type)));
    elseif strcmp(type, 'gcm')
        title(sprintf('%s, Annual', par.model));
    end
    % cb = colorbar('limits', [-10 10], 'ytick', [-10:2:10], 'location', 'eastoutside');
    % cb.TickLabelInterpreter = 'latex'; cb.Label.Interpreter = 'latex';
    % ylabel(cb, '$\Gamma_m$ (K/km)');
    set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_sq)
    set(gca, 'xlim', [-10 75], 'xtick', [-10:10:70 75], 'ylim', [0.2 1], 'ytick', [0.2:0.1:1], 'ydir', 'reverse', 'yminortick', 'on', 'tickdir', 'out');
    print(sprintf('%s/ga_malr_diff/ga_m_lat_sigma_sc', par.plotdir), '-dpng', '-r300');
    close;

end
