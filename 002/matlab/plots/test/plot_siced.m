function plot_siced(type, par)
% plot sicedow depth
    make_dirs(type, par)

    % load data
    [~, ~, ~, lat, par] = load_flux(type, par);
    if strcmp(type, 'era5') | strcmp(type, 'erai') | strcmp(type, 'era5c')
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s', type, par.(type).yr_span);
    elseif strcmp(type, 'gcm')
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s/%s', type, par.model, par.gcm.clim);
    elseif strcmp(type, 'echam')
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s', type, par.echam.clim);
    end
    load(sprintf('%s/grid.mat', prefix)); % read grid data
    load(sprintf('%s/siced.mat', prefix)); % read clear sky albedo data

    [mesh_lat, mesh_mon] = meshgrid([1:12], lat);

    siced = squeeze(nanmean(siced,1)); % zonal mean

    % mon x lat of sicedow depth
    figure(); clf; hold all;
    cmp = colCog(20);
    colormap(flip(cmp));
    [C,h] = contourf(mesh_lat, mesh_mon, siced, [0:0.5:10 30 40 100], 'linecolor', 'none');
    % clabel(C, h, 0:0.5:10, 'fontsize', 6, 'interpreter', 'latex', 'color', 0.75*[1 1 1]);
    caxis([-5 5]);
    ylabel('Latitude (deg)');
    if strcmp(type, 'era5') | strcmp(type, 'erai') | strcmp(type, 'era5c')
        title(sprintf('%s', upper(type)));
    elseif strcmp(type, 'echam')
        title(sprintf('%s, %s', upper(type), par.echam.(par.echam.clim)));
    elseif strcmp(type, 'gcm')
        title(sprintf('%s', par.model));
    end
    cb = colorbar('limits', [0 5], 'ytick', [0:0.5:5], 'location', 'eastoutside');
    cb.TickLabelInterpreter = 'latex'; cb.Label.Interpreter = 'latex';
    ylabel(cb, sprintf('Sea ice depth (m)'));
    set(gca, 'xlim', [1 12], 'xtick', [1:12], 'xticklabels', par.monlabel, 'ylim', [-90 90], 'ytick', [-90:30:90], 'yminortick', 'on', 'tickdir', 'out');
    print(sprintf('%s/siced/siced_mon_lat', par.plotdir), '-dpng', '-r300');
    close;

end
