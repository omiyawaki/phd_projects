function plot_tas(type, par)
% plot 2 m temperature
    make_dirs(type, par)

    % load data
    prefix = make_prefix(type, par);
    plotdir = make_plotdir(type, par);

    load(sprintf('%s/grid.mat', prefix)); % read grid data
    load(sprintf('%s/srfc.mat', prefix));

    tas = rename_tas(type, srfc);

    [mesh_lat, mesh_mon] = meshgrid([1:12], grid.dim2.lat);

    tas = squeeze(nanmean(tas,1)); % zonal mean
    tas_ann = squeeze(nanmean(tas,2)); % annual mean

    % mon x lat of tasow depth
    figure(); clf; hold all; box on;
    % cmp = colCog(20);
    % colormap(cmp);
    [C,h] = contour(mesh_lat, mesh_mon, tas, [100:10:400]);
    % [C,h] = contour(mesh_lat, mesh_mon, tas, [270 270], 'linecolor', 'w');
    clabel(C, h, [100 150 200 250 270 300], 'fontsize', 6, 'interpreter', 'latex');
    ylabel('Latitude (deg)');
    make_title_type(type, par);
    if strcmp(type, 'echam') & ~contains(par.echam.clim, 'echr')
        caxis([150 270]);
        cb = colorbar('limits', [150 270], 'ytick', [150:10:270], 'location', 'eastoutside');
    else
        cb = colorbar('location', 'eastoutside');
    end
    cb.TickLabelInterpreter = 'latex'; cb.Label.Interpreter = 'latex';
    ylabel(cb, sprintf('$T_{2\\,\\mathrm{m}}$ (K)'));
    set(gca, 'xlim', [1 12], 'xtick', [1:12], 'xticklabels', par.monlabelnh, 'ylim', [-90 90], 'ytick', [-90:30:90], 'yminortick', 'on', 'tickdir', 'out');
    print(sprintf('%s/tas/tas_mon_lat', plotdir), '-dpng', '-r300');
    close;

    % mon x lat of surface pressure
    figure(); clf; hold all; box on;
    plot(grid.dim2.lat, tas_ann, 'k');
    xlabel('Latitude (deg)');
    ylabel('$T_{\mathrm{2\, m}}$ (hPa)');
    make_title_type(type, par);
    set(gca, 'xlim', [-90 90], 'xtick', [-90:30:90], 'ylim', [230 300], 'xminortick', 'on', 'yminortick', 'on', 'tickdir', 'out');
    print(sprintf('%s/tas/tas_lat', plotdir), '-dpng', '-r300');
    close;

end
