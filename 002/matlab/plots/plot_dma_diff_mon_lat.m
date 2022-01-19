function plot_dt_norm_mon_lat(type, par)

    % plot inversion strength
    make_dirs_si_bl(type, par)

    prefix = make_prefix(type, par);
    prefix_proc = make_prefix_proc(type, par);
    plotdir = make_plotdir(type, par);

    load(sprintf('%s/si_bl_%g/dt_norm_si_mon_lat_%g.mat', prefix_proc, par.si_bl, par.si_up));

    [mesh_lat, mesh_mon] = meshgrid([1:12], lat);

    for l = par.land_list; land = l{1};
        if strcmp(land, 'lo'); land_text = 'Land + Ocean';
        elseif strcmp(land, 'l'); land_text = 'Land';
        elseif strcmp(land, 'o'); land_text = 'Ocean';
        end

        % mon x lat of diff
        figure(); clf; hold all; box on;
        cmp = colCog(20);
        colormap(flipud(cmp));
        contourf(mesh_lat, mesh_mon, dt_norm.(land), -1000:10:1000, 'linecolor', 'none');
        [C,h]=contour(mesh_lat, mesh_mon, dt_norm.(land), -1000:10:1000, '-w', 'linewidth', 0.1);
        contour(mesh_lat, mesh_mon, dt_norm.(land), [0 0], 'color', 0.75*[1 1 1]);
        clabel(C, h, -100:10:100, 'fontsize', 6, 'interpreter', 'latex', 'color', 0.75*[1 1 1]);
        caxis([-100 100]);
        ylabel('Latitude (deg)');
        make_title_type(type, par);
        cb = colorbar('limits', [-50 50], 'ytick', [-50:10:50], 'location', 'eastoutside');
        cb.TickLabelInterpreter = 'latex'; cb.Label.Interpreter = 'latex';
        ylabel(cb, sprintf('$\\left\\langle(\\Delta T_m - \\Delta T)/\\Delta T_m\\right\\rangle_{%g-%g}$ (\\%%)', par.si_bl, par.si_up));
        set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_verywide);
        set(gca, 'xlim', [1 12], 'xtick', [1:12], 'xticklabels', par.monlabelnh, 'ylim', [-90 90], 'ytick', [-90:30:90], 'yminortick', 'on', 'tickdir', 'out');
        % print(sprintf('%s/dt_norm/si_bl_%g/%s/dt_norm_mon_lat_%g.png', plotdir, par.si_bl, land, par.si_up), '-dpng', '-r300');
        folder = sprintf('%s/dt_norm/si_bl_%g/%s/', plotdir, par.si_bl, land);
        if ~exist(folder, 'dir')
            mkdir(folder)
        end
        print(sprintf('%sdt_norm_mon_lat_%g.png', folder, par.si_up), '-dpng', '-r300');
        if par.make_tikz & par.si_bl == 0.7
            data = [mesh_lat(:) mesh_mon(:) dt_norm.(land)(:) ];
            matlab2tikz(sprintf('%sdt_norm_mon_lat_%g.tex', folder, par.si_up));
            save(sprintf('%sdt_norm_mon_lat_%g.dat', folder, par.si_up), 'data', '-ASCII')
            clear data
        end
        close;

        % ANOT mon x lat of diff
        figure(); clf; hold all; box on;
        cmp = colCog(20);
        colormap(flipud(cmp));
        contourf(mesh_lat, mesh_mon, dt_norm.(land), -1000:5:1000, 'linecolor', 'none');
        [C,h]=contour(mesh_lat, mesh_mon, dt_norm.(land), -1000:5:1000, '-w', 'linewidth', 0.1);
        contour(mesh_lat, mesh_mon, dt_norm.(land), par.ga_thresh*[1 1], 'color', par.orange, 'linewidth', 2);
        contour(mesh_lat, mesh_mon, dt_norm.(land), [0 0], 'color', 0.75*[1 1 1]);
        clabel(C, h, -100:10:100, 'fontsize', 6, 'interpreter', 'latex', 'color', 0.75*[1 1 1]);
        caxis([-50 50]);
        ylabel('Latitude (deg)');
        make_title_type(type, par);
        cb = colorbar('limits', [-50 50], 'ytick', [-50:10:50], 'location', 'eastoutside');
        cb.TickLabelInterpreter = 'latex'; cb.Label.Interpreter = 'latex';
        ylabel(cb, sprintf('$\\left\\langle(\\Delta T_m - \\Gamma)/\\Delta T_m\\right\\rangle_{%g-%g}$ (\\%%)', par.si_bl, par.si_up));
        set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_verywide);
        set(gca, 'xlim', [1 12], 'xtick', [1:12], 'xticklabels', par.monlabelnh, 'ylim', [-90 90], 'ytick', [-90:30:90], 'yminortick', 'on', 'tickdir', 'out');
        % print(sprintf('%s/dt_norm/si_bl_%g/%s/dt_norm_mon_lat_%g_anot.png', plotdir, par.si_bl, land, par.si_up), '-dpng', '-r300');
        print(sprintf('%sdt_norm_mon_lat_%g_anot.png', folder, par.si_up), '-dpng', '-r300');
        close;

    end
    
end
