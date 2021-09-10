function plot_dtemp_diff_mon_lat(type, par)

    if ~strcmp(type, 'gcm')
        error('This analysis only works for CMIP5 historical')
    end

    partitle = par;
    partitle.gcm.clim = 'RCP8.5$-$historical';

    % plot temperature response ratio

    make_dirs(type, par)

    prefix = make_prefix(type, par);
    prefix_proc = make_prefix_proc(type, par);
    plotdir = make_plotdir(type, par);

    load(sprintf('%s/grid.mat', prefix)); % read grid data
    load(sprintf('%s/dtempsi_mon_lat.mat', prefix_proc)); % load lat x mon RCAE_ALT data

    dtemp = permute(dtempsi.lo, [3 1 2]);
    dtemp_u = squeeze(interp1(grid.dim3.si, dtemp, 0.3)); % upper trop dtemp
    dtemp_l = squeeze(interp1(grid.dim3.si, dtemp, 1)); % lower trop dtemp
    dtemp_r = dtemp_u./dtemp_l; % upper to lower temperature response ratio

    [mesh_lat, mesh_mon] = meshgrid([1:12], lat);

    for l = par.land_list; land = l{1};
        if strcmp(land, 'lo'); land_text = 'Land + Ocean';
        elseif strcmp(land, 'l'); land_text = 'Land';
        elseif strcmp(land, 'o'); land_text = 'Ocean';
        end

        % dtemp_r mon x lat
        figure(); clf; hold all; box on;
        cmp = colCog(40);
        colormap(flip(cmp));
        contourf(mesh_lat, mesh_mon, dtemp_r, -0.5:0.1:4, 'linecolor', 'none');
        [C,h]=contour(mesh_lat, mesh_mon, dtemp_r, -0.5:0.5:4, '-w', 'linewidth', 0.1);
        % contour(mesh_lat, mesh_mon, dtemp_r, [0 0], 'color', 0.75*[1 1 1]);
        clabel(C, h, 0:0.5:4, 'fontsize', 6, 'interpreter', 'latex', 'color', 0.75*[1 1 1]);
        caxis([-1 3]);
        ylabel('Latitude (deg)');
        make_title_type(type, partitle);
        cb = colorbar('limits', [0 3], 'ytick', [0:0.2:3], 'location', 'eastoutside');
        cb.TickLabelInterpreter = 'latex'; cb.Label.Interpreter = 'latex';
        ylabel(cb, sprintf('$\\Delta T_{0.3}/\\Delta T_{1.0}$ (unitless)'));
        set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_verywide);
        set(gca, 'xlim', [1 12], 'xtick', [1:12], 'xticklabels', par.monlabelnh, 'ylim', [-90 90], 'ytick', [-90:30:90], 'yminortick', 'on', 'tickdir', 'out');
        folder = sprintf('%s/dtempr/%s/', plotdir, land);
        if ~exist(folder, 'dir')
            mkdir(folder)
        end
        print(sprintf('%sdtempr_mon_lat', folder), '-dpng', '-r300');
        close;

        % dtemp_u mon x lat
        figure(); clf; hold all; box on;
        cmp = colCog(40);
        colormap(cmp);
        contourf(mesh_lat, mesh_mon, dtemp_u, -5:0.25:10, 'linecolor', 'none');
        [C,h]=contour(mesh_lat, mesh_mon, dtemp_u, -5:0.5:10, '-w', 'linewidth', 0.1);
        % contour(mesh_lat, mesh_mon, dtemp_r, [0 0], 'color', 0.75*[1 1 1]);
        clabel(C, h, 0:2:10, 'fontsize', 6, 'interpreter', 'latex', 'color', 0.75*[1 1 1]);
        caxis([-10 10]);
        ylabel('Latitude (deg)');
        make_title_type(type, partitle);
        cb = colorbar('limits', [0 10], 'ytick', [0:2:10], 'location', 'eastoutside');
        cb.TickLabelInterpreter = 'latex'; cb.Label.Interpreter = 'latex';
        ylabel(cb, sprintf('$\\Delta T_{0.3}$ (K)'));
        set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_verywide);
        set(gca, 'xlim', [1 12], 'xtick', [1:12], 'xticklabels', par.monlabelnh, 'ylim', [-90 90], 'ytick', [-90:30:90], 'yminortick', 'on', 'tickdir', 'out');
        folder = sprintf('%s/dtempr/%s/', plotdir, land);
        if ~exist(folder, 'dir')
            mkdir(folder)
        end
        print(sprintf('%sdtempu_mon_lat', folder), '-dpng', '-r300');
        close;

        % dtemp_l mon x lat
        figure(); clf; hold all; box on;
        cmp = colCog(40);
        colormap(cmp);
        contourf(mesh_lat, mesh_mon, dtemp_l, 0:0.5:20, 'linecolor', 'none');
        [C,h]=contour(mesh_lat, mesh_mon, dtemp_l, 0:1:20, '-w', 'linewidth', 0.1);
        % contour(mesh_lat, mesh_mon, dtemp_r, [0 0], 'color', 0.75*[1 1 1]);
        clabel(C, h, 0:2:20, 'fontsize', 6, 'interpreter', 'latex', 'color', 0.75*[1 1 1]);
        caxis([-20 20]);
        ylabel('Latitude (deg)');
        make_title_type(type, partitle);
        cb = colorbar('limits', [0 20], 'ytick', [0:5:20], 'location', 'eastoutside');
        cb.TickLabelInterpreter = 'latex'; cb.Label.Interpreter = 'latex';
        ylabel(cb, sprintf('$\\Delta T_{1.0}$ (K)'));
        set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_verywide);
        set(gca, 'xlim', [1 12], 'xtick', [1:12], 'xticklabels', par.monlabelnh, 'ylim', [-90 90], 'ytick', [-90:30:90], 'yminortick', 'on', 'tickdir', 'out');
        folder = sprintf('%s/dtempr/%s/', plotdir, land);
        if ~exist(folder, 'dir')
            mkdir(folder)
        end
        print(sprintf('%sdtempl_mon_lat', folder), '-dpng', '-r300');
        close;

    end
    
end
