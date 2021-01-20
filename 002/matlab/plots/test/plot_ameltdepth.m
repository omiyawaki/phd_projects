function plot_ameltdepth(type, par)
% plot flwacow depth
    make_dirs(type, par)

    prefix = make_prefix(type, par);
    prefix_proc = make_prefix_proc(type, par);
    plotdir = make_plotdir(type, par);

    load(sprintf('%s/grid.mat', prefix)); % read grid data
    load(sprintf('%s/ameltdepth.mat', prefix)); ameltdepth = echamvar; % read clear sky albedo data

    [mesh_lat, mesh_mon] = meshgrid([1:12], grid.dim2.lat);

    ameltdepth = squeeze(nanmean(ameltdepth,1)); % zonal mean

    % mon x lat of ameltdepthow depth
    figure(); clf; hold all;
    cmp = colCog(50);
    colormap(cmp);
    contourf(mesh_lat, mesh_mon, ameltdepth, [0:0.01:0.25], 'linecolor', 'none');
    % contourf(mesh_lat, mesh_mon, ameltdepth, 'linecolor', 'none');
    % [C,h] = contour(mesh_lat, mesh_mon, ameltdepth, [1,1], 'linecolor', 'k');
    % clabel(C, h, 1, 'fontsize', 6, 'interpreter', 'latex', 'color', 'k');
    caxis([-0.25 0.25]);
    ylabel('Latitude (deg)');
    if strcmp(type, 'era5') | strcmp(type, 'erai') | strcmp(type, 'era5c')
        title(sprintf('%s', upper(type)));
    elseif strcmp(type, 'echam')
        title(sprintf('%s, %s', upper(type), par.echam.(par.echam.clim)));
    elseif strcmp(type, 'gcm')
        title(sprintf('%s', par.model));
    end
    cb = colorbar('limits', [0 0.25], 'ytick', [0:0.05:0.25], 'location', 'eastoutside');
    %cb = colorbar('location', 'eastoutside');
    cb.TickLabelInterpreter = 'latex'; cb.Label.Interpreter = 'latex';
    ylabel(cb, sprintf('Meltpond depth (m)'));
    set(gca, 'xlim', [1 12], 'xtick', [1:12], 'xticklabels', par.monlabel, 'ylim', [-90 90], 'ytick', [-90:30:90], 'yminortick', 'on', 'tickdir', 'out');
    print(sprintf('%s/ameltdepth/ameltdepth_mon_lat', plotdir), '-dpng', '-r300');
    close;

end
