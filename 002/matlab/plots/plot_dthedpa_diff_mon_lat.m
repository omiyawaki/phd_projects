function plot_dthedpa_mon_lat(type, par)

    % plot inversion strength
    make_dirs_si_bl(type, par)

    prefix = make_prefix(type, par);
    prefix_proc = make_prefix_proc(type, par);
    plotdir = make_plotdir(type, par);

    load(sprintf('%s/si_bl_%g/dthedpa_si_mon_lat_%g.mat', prefix_proc, par.si_bl, par.si_up));

    [mesh_lat, mesh_mon] = meshgrid([1:12], lat);

    for l = par.land_list; land = l{1};
        if strcmp(land, 'lo'); land_text = 'Land + Ocean';
        elseif strcmp(land, 'l'); land_text = 'Land';
        elseif strcmp(land, 'o'); land_text = 'Ocean';
        end

        % mon x lat of diff
        figure(); clf; hold all; box on;
        cmp = colCog(40);
        colormap(flipud(cmp));
        contourf(mesh_lat, mesh_mon, dthedpa.(land), -0.2:0.005:0.2, 'linecolor', 'none');
        [C,h]=contour(mesh_lat, mesh_mon, dthedpa.(land), -0.2:0.005:0.2, '-w', 'linewidth', 0.1);
        contour(mesh_lat, mesh_mon, dthedpa.(land), [0 0], 'color', 0.75*[1 1 1]);
        clabel(C, h, -100:10:100, 'fontsize', 6, 'interpreter', 'latex', 'color', 0.75*[1 1 1]);
        caxis([-0.1 0.1]);
        ylabel('Latitude (deg)');
        make_title_type(type, par);
        cb = colorbar('limits', [-0.1 0.1], 'ytick', [-0.1:0.02:0.1], 'location', 'eastoutside');
        cb.TickLabelInterpreter = 'latex'; cb.Label.Interpreter = 'latex';
        ylabel(cb, sprintf('$\\left[-\\frac{\\partial \\theta_e}{\\partial p}\\right]_{%g-%g}$ (\\%%)', par.si_bl, par.si_up));
        set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_verywide);
        set(gca, 'xlim', [1 12], 'xtick', [1:12], 'xticklabels', par.monlabelnh, 'ylim', [-90 90], 'ytick', [-90:30:90], 'yminortick', 'on', 'tickdir', 'out');
        print(sprintf('%s/dthedpa/si_bl_%g/%s/dthedpa_mon_lat_%g.png', plotdir, par.si_bl, land, par.si_up), '-dpng', '-r300');
        close;

        % ANOT mon x lat of diff
        figure(); clf; hold all; box on;
        cmp = colCog(40);
        colormap(flipud(cmp));
        contourf(mesh_lat, mesh_mon, dthedpa.(land), -0.2:0.005:0.2, 'linecolor', 'none');
        [C,h]=contour(mesh_lat, mesh_mon, dthedpa.(land), -0.2:0.005:0.2, '-w', 'linewidth', 0.1);
        contour(mesh_lat, mesh_mon, dthedpa.(land), par.ga_thresh*[1 1], 'color', par.orange, 'linewidth', 2);
        contour(mesh_lat, mesh_mon, dthedpa.(land), [0 0], 'color', 0.75*[1 1 1]);
        clabel(C, h, -100:10:100, 'fontsize', 6, 'interpreter', 'latex', 'color', 0.75*[1 1 1]);
        caxis([-0.1 0.1]);
        ylabel('Latitude (deg)');
        make_title_type(type, par);
        cb = colorbar('limits', [-50 50], 'ytick', [-50:10:50], 'location', 'eastoutside');
        cb.TickLabelInterpreter = 'latex'; cb.Label.Interpreter = 'latex';
        ylabel(cb, sprintf('$\\left[-\\frac{\\partial \\theta_e}{\\partial p}\\right]_{%g-%g}$ (\\%%)', par.si_bl, par.si_up));
        set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_verywide);
        set(gca, 'xlim', [1 12], 'xtick', [1:12], 'xticklabels', par.monlabelnh, 'ylim', [-90 90], 'ytick', [-90:30:90], 'yminortick', 'on', 'tickdir', 'out');
        print(sprintf('%s/dthedpa/si_bl_%g/%s/dthedpa_mon_lat_%g_anot.png', plotdir, par.si_bl, land, par.si_up), '-dpng', '-r300');
        close;

    end
    
end
