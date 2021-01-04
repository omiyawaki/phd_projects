function plot_r1z_lat(type, par)
% Compare GCM R1 with ERA5

    grid_era5 = load(sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/grid.mat', 'era5'));
    grid_gcm = load(sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s/%s/grid.mat', type, par.model, par.gcm.clim));
    flux_era5 = load(sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/flux_zt.mat', 'era5', par.lat_interp));
    flux_gcm = load(sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s/%s/flux_zt.mat', type, par.model, par.gcm.clim, par.lat_interp));
    par.plotdir = sprintf('./figures/%s/%s/%s/%s', type, par.model, par.gcm.clim, par.lat_interp);

    % R1Z
    figure(); clf; hold all; box on;
    r1z_era5=flux_era5.flux_zt.lo.ann.res.mse./flux_era5.flux_zt.lo.ann.ra.mse;
    r1z_gcm=flux_gcm.flux_zt.lo.ann.res.mse./flux_gcm.flux_zt.lo.ann.ra.mse;
    ylim_lo = min([r1z_era5, r1z_gcm]); if isnan(ylim_lo)|ylim_lo==0; ylim_lo = -1; end;
    ylim_up = max([r1z_era5, r1z_gcm]); if isnan(ylim_up)|ylim_up==0; ylim_up = 1; end
    raemin = par.ga;
    vertices = [-90 raemin; 90 raemin; 90 ylim_up; -90 ylim_up];
    patch(vertices(:,1), vertices(:,2), par.blue, 'edgecolor', 'none', 'facealpha', 0.5);
    rcemax = par.ep;
    vertices = [-90 rcemax; 90 rcemax; 90 raemin; -90 raemin];
    patch(vertices(:,1), vertices(:,2), 0.75*[1 1 1], 'edgecolor', 'none', 'facealpha', 0.5);
    vertices = [-90 ylim_lo; 90 ylim_lo; 90 rcemax; -90 rcemax];
    patch(vertices(:,1), vertices(:,2), par.orange, 'edgecolor', 'none', 'facealpha', 0.5);
    % plot(grid_gcm.grid.dim2.lat,r1z_gcm, '-k');
    plot(grid_era5.grid.dim2.lat,r1z_era5, '-k');
    % legend('RAE', 'RCAE', 'RCE', sprintf('%s', par.model), 'ERA5', 'location', 'eastoutside')
    xlabel('latitude (deg)');
    ylabel('$R_1$ (unitless)')
    axis('tight');
    set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_wide)
    set(gca, 'fontsize', par.fs, 'xtick', [-90:30:90], 'xminortick', 'on', 'yminortick', 'on')
    print(sprintf('%s/energy-flux/%s/%s/%s-r1z-comp-era5', par.plotdir, 'lo', 'ann', 'mse'), '-dpng', '-r300');
    close;

end
