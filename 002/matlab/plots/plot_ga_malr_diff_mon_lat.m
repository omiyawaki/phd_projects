function plot_ga_malr_diff_mon_lat(type, par)

    % plot inversion strength
    make_dirs_si_bl(type, par)

    prefix = make_prefix(type, par);
    prefix_proc = make_prefix_proc(type, par);
    plotdir = make_plotdir(type, par);

    load(sprintf('%s/si_bl_%g/ga_malr_diff_si_mon_lat_%g.mat', prefix_proc, par.si_bl, par.si_up));
    load(sprintf('%s/si_bl_%g/ga_malr_bl_diff_si_mon_lat.mat', prefix_proc, par.si_bl));
    load(sprintf('%s/si_bl_%g/ga_dalr_bl_diff_si_mon_lat.mat', prefix_proc, par.si_bl));

    [mesh_lat, mesh_mon] = meshgrid([1:12], lat);

    for l = par.land_list; land = l{1};
        if strcmp(land, 'lo'); land_text = 'Land + Ocean';
        elseif strcmp(land, 'l'); land_text = 'Land';
        elseif strcmp(land, 'o'); land_text = 'Ocean';
        end

        % DALR BL mon x lat of diff
        figure(); clf; hold all; box on;
        cmp = colCog(20);
        colormap(flipud(cmp));
        contourf(mesh_lat, mesh_mon, ga_dalr_bl_diff.(land), -1000:10:2000, 'linecolor', 'none');
        [C,h]=contour(mesh_lat, mesh_mon, ga_dalr_bl_diff.(land), -1000:10:2000, '-w', 'linewidth', 0.1);
        contour(mesh_lat, mesh_mon, ga_dalr_bl_diff.(land), [0 0], 'color', 0.75*[1 1 1]);
        clabel(C, h, -100:20:100, 'fontsize', 6, 'interpreter', 'latex', 'color', 0.75*[1 1 1]);
        caxis([-100 100]);
        ylabel('Latitude (deg)');
        make_title_type(type, par);
        cb = colorbar('limits', [0 100], 'ytick', [0:10:100], 'location', 'eastoutside');
        cb.TickLabelInterpreter = 'latex'; cb.Label.Interpreter = 'latex';
        ylabel(cb, sprintf('$\\left\\langle(\\Gamma_d - \\Gamma)/\\Gamma_d\\right\\rangle_{1.0-%g}$ (\\%%)', par.si_bl));
        set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_verywide);
        set(gca, 'xlim', [1 12], 'xtick', [1:12], 'xticklabels', par.monlabelnh, 'ylim', [-90 90], 'ytick', [-90:30:90], 'yminortick', 'on', 'tickdir', 'out');
        folder = sprintf('%s/ga_malr_diff/si_bl_%g/%s/', plotdir, par.si_bl, land);
        if ~exist(folder, 'dir')
            mkdir(folder)
        end
        print(sprintf('%sga_dalr_bl_diff_mon_lat', folder), '-dpng', '-r300');
        close;

        % ANOT DALR BL mon x lat of diff
        figure(); clf; hold all; box on;
        cmp = colCog(20);
        colormap(flipud(cmp));
        contourf(mesh_lat, mesh_mon, ga_dalr_bl_diff.(land), -1000:10:2000, 'linecolor', 'none');
        [C,h]=contour(mesh_lat, mesh_mon, ga_dalr_bl_diff.(land), -1000:10:2000, '-w', 'linewidth', 0.1);
        contour(mesh_lat, mesh_mon, ga_dalr_bl_diff.(land), par.ga_bl_thresh*[1 1], 'color', par.cyan, 'linewidth', 2);
        contour(mesh_lat, mesh_mon, ga_dalr_bl_diff.(land), [0 0], 'color', 0.75*[1 1 1]);
        clabel(C, h, -100:20:100, 'fontsize', 6, 'interpreter', 'latex', 'color', 0.75*[1 1 1]);
        caxis([-100 100]);
        ylabel('Latitude (deg)');
        make_title_type(type, par);
        cb = colorbar('limits', [0 100], 'ytick', [0:10:100], 'location', 'eastoutside');
        cb.TickLabelInterpreter = 'latex'; cb.Label.Interpreter = 'latex';
        ylabel(cb, sprintf('$\\left\\langle(\\Gamma_d - \\Gamma)/\\Gamma_d\\right\\rangle_{1.0-%g}$ (\\%%)', par.si_bl));
        set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_verywide);
        set(gca, 'xlim', [1 12], 'xtick', [1:12], 'xticklabel', par.monlabelnh, 'ylim', [-90 90], 'ytick', [-90:30:90], 'yminortick', 'on', 'tickdir', 'out');
        % print(sprintf('%s/ga_malr_diff/si_bl_%g/%s/ga_dalr_bl_diff_mon_lat_anot', plotdir, par.si_bl, land), '-dpng', '-r300');
        print(sprintf('%sga_dalr_bl_diff_mon_lat_anot', folder), '-dpng', '-r300');
        close;

        % MALR BL mon x lat of diff
        figure(); clf; hold all; box on;
        cmp = colCog(20);
        colormap(flipud(cmp));
        contourf(mesh_lat, mesh_mon, ga_malr_bl_diff.(land), -1000:10:2000, 'linecolor', 'none');
        [C,h]=contour(mesh_lat, mesh_mon, ga_malr_bl_diff.(land), -1000:10:2000, '-w', 'linewidth', 0.1);
        contour(mesh_lat, mesh_mon, ga_malr_bl_diff.(land), [0 0], 'color', 0.75*[1 1 1]);
        clabel(C, h, -100:20:100, 'fontsize', 6, 'interpreter', 'latex', 'color', 0.75*[1 1 1]);
        caxis([-100 100]);
        ylabel('Latitude (deg)');
        make_title_type(type, par);
        cb = colorbar('limits', [-100 100], 'ytick', [-100:20:100], 'location', 'eastoutside');
        cb.TickLabelInterpreter = 'latex'; cb.Label.Interpreter = 'latex';
        ylabel(cb, sprintf('$\\left\\langle(\\Gamma_m - \\Gamma)/\\Gamma_m\\right\\rangle_{1.0-%g}$ (\\%%)', par.si_bl));
        set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_verywide);
        set(gca, 'xlim', [1 12], 'xtick', [1:12], 'xticklabels', par.monlabelnh, 'ylim', [-90 90], 'ytick', [-90:30:90], 'yminortick', 'on', 'tickdir', 'out');
        % print(sprintf('%s/ga_malr_diff/si_bl_%g/%s/ga_malr_bl_diff_mon_lat', plotdir, par.si_bl, land), '-dpng', '-r300');
        print(sprintf('%sga_malr_bl_diff_mon_lat', folder), '-dpng', '-r300');
        if par.make_tikz & par.si_bl == 0.9
            data = [mesh_lat(:) mesh_mon(:) ga_malr_bl_diff.(land)(:) ];
            matlab2tikz(sprintf('%sga_malr_bl_diff_mon_lat.tex', folder));
            save(sprintf('%sga_malr_bl_diff_mon_lat.dat', folder), 'data', '-ASCII')
            clear data
        end
        close;

        % ANOT MALR BL mon x lat of diff
        figure(); clf; hold all; box on;
        cmp = colCog(20);
        colormap(flipud(cmp));
        contourf(mesh_lat, mesh_mon, ga_malr_bl_diff.(land), -1000:10:2000, 'linecolor', 'none');
        [C,h]=contour(mesh_lat, mesh_mon, ga_malr_bl_diff.(land), -1000:10:2000, '-w', 'linewidth', 0.1);
        contour(mesh_lat, mesh_mon, ga_malr_bl_diff.(land), par.ga_bl_thresh*[1 1], 'color', par.cyan, 'linewidth', 2);
        contour(mesh_lat, mesh_mon, ga_malr_bl_diff.(land), [0 0], 'color', 0.75*[1 1 1]);
        clabel(C, h, -100:20:100, 'fontsize', 6, 'interpreter', 'latex', 'color', 0.75*[1 1 1]);
        caxis([-100 100]);
        ylabel('Latitude (deg)');
        make_title_type(type, par);
        cb = colorbar('limits', [-100 100], 'ytick', [-100:20:100], 'location', 'eastoutside');
        cb.TickLabelInterpreter = 'latex'; cb.Label.Interpreter = 'latex';
        ylabel(cb, sprintf('$\\left\\langle(\\Gamma_m - \\Gamma)/\\Gamma_m\\right\\rangle_{1.0-%g}$ (\\%%)', par.si_bl));
        set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_verywide);
        set(gca, 'xlim', [1 12], 'xtick', [1:12], 'xticklabel', par.monlabelnh, 'ylim', [-90 90], 'ytick', [-90:30:90], 'yminortick', 'on', 'tickdir', 'out');
        % print(sprintf('%s/ga_malr_diff/si_bl_%g/%s/ga_malr_bl_diff_mon_lat_anot', plotdir, par.si_bl, land), '-dpng', '-r300');
        print(sprintf('%sga_malr_bl_diff_mon_lat_anot', folder), '-dpng', '-r300');
        close;

        % mon x lat of diff
        figure(); clf; hold all; box on;
        cmp = colCog(20);
        colormap(flipud(cmp));
        contourf(mesh_lat, mesh_mon, ga_malr_diff.(land), -1000:5:1000, 'linecolor', 'none');
        [C,h]=contour(mesh_lat, mesh_mon, ga_malr_diff.(land), -1000:5:1000, '-w', 'linewidth', 0.1);
        contour(mesh_lat, mesh_mon, ga_malr_diff.(land), [0 0], 'color', 0.75*[1 1 1]);
        clabel(C, h, -100:10:100, 'fontsize', 6, 'interpreter', 'latex', 'color', 0.75*[1 1 1]);
        caxis([-50 50]);
        ylabel('Latitude (deg)');
        make_title_type(type, par);
        cb = colorbar('limits', [-50 50], 'ytick', [-50:10:50], 'location', 'eastoutside');
        cb.TickLabelInterpreter = 'latex'; cb.Label.Interpreter = 'latex';
        ylabel(cb, sprintf('$\\left\\langle(\\Gamma_m - \\Gamma)/\\Gamma_m\\right\\rangle_{%g-%g}$ (\\%%)', par.si_bl, par.si_up));
        set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_verywide);
        set(gca, 'xlim', [1 12], 'xtick', [1:12], 'xticklabels', par.monlabelnh, 'ylim', [-90 90], 'ytick', [-90:30:90], 'yminortick', 'on', 'tickdir', 'out');
        % print(sprintf('%s/ga_malr_diff/si_bl_%g/%s/ga_malr_diff_mon_lat_%g.png', plotdir, par.si_bl, land, par.si_up), '-dpng', '-r300');
        print(sprintf('%sga_malr_diff_mon_lat_%g.png', folder, par.si_up), '-dpng', '-r300');
        if par.make_tikz & par.si_bl == 0.7
            data = [mesh_lat(:) mesh_mon(:) ga_malr_diff.(land)(:) ];
            matlab2tikz(sprintf('%sga_malr_diff_mon_lat_%g.tex', folder, par.si_up));
            save(sprintf('%sga_malr_diff_mon_lat_%g.dat', folder, par.si_up), 'data', '-ASCII')
            clear data
        end
        close;

        % ANOT mon x lat of diff
        figure(); clf; hold all; box on;
        cmp = colCog(20);
        colormap(flipud(cmp));
        contourf(mesh_lat, mesh_mon, ga_malr_diff.(land), -1000:5:1000, 'linecolor', 'none');
        [C,h]=contour(mesh_lat, mesh_mon, ga_malr_diff.(land), -1000:5:1000, '-w', 'linewidth', 0.1);
        contour(mesh_lat, mesh_mon, ga_malr_diff.(land), par.ga_thresh*[1 1], 'color', par.orange, 'linewidth', 2);
        contour(mesh_lat, mesh_mon, ga_malr_diff.(land), [0 0], 'color', 0.75*[1 1 1]);
        clabel(C, h, -100:10:100, 'fontsize', 6, 'interpreter', 'latex', 'color', 0.75*[1 1 1]);
        caxis([-50 50]);
        ylabel('Latitude (deg)');
        make_title_type(type, par);
        cb = colorbar('limits', [-50 50], 'ytick', [-50:10:50], 'location', 'eastoutside');
        cb.TickLabelInterpreter = 'latex'; cb.Label.Interpreter = 'latex';
        ylabel(cb, sprintf('$\\left\\langle(\\Gamma_m - \\Gamma)/\\Gamma_m\\right\\rangle_{%g-%g}$ (\\%%)', par.si_bl, par.si_up));
        set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_verywide);
        set(gca, 'xlim', [1 12], 'xtick', [1:12], 'xticklabels', par.monlabelnh, 'ylim', [-90 90], 'ytick', [-90:30:90], 'yminortick', 'on', 'tickdir', 'out');
        % print(sprintf('%s/ga_malr_diff/si_bl_%g/%s/ga_malr_diff_mon_lat_%g_anot.png', plotdir, par.si_bl, land, par.si_up), '-dpng', '-r300');
        print(sprintf('%sga_malr_diff_mon_lat_%g_anot.png', folder, par.si_up), '-dpng', '-r300');
        close;

    end
    
end
