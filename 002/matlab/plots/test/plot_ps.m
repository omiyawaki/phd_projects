function plot_ps(type, par)
% plot surface pressure
    make_dirs(type, par)

    % load data
    prefix = make_prefix(type, par);
    plotdir = make_plotdir(type, par);

    load(sprintf('%s/grid.mat', prefix)); % read grid data
    load(sprintf('%s/srfc.mat', prefix));

    ps = rename_ps(type, srfc);

    ps = squeeze(nanmean(ps,1)); % zonal mean
    ps = squeeze(nanmean(ps,2)); % annual mean

    % mon x lat of surface pressure
    figure(); clf; hold all; box on;
    plot(grid.dim2.lat, ps/100, 'k');
    xlabel('Latitude (deg)');
    ylabel('$p_{s}$ (hPa)');
    make_title_type(type, par);
    set(gca, 'xlim', [-90 90], 'xtick', [-90:30:90], 'ylim', [600 1100], 'xminortick', 'on', 'yminortick', 'on', 'tickdir', 'out');
    print(sprintf('%s/ps/ps_lat', plotdir), '-dpng', '-r300');
    close;

end
