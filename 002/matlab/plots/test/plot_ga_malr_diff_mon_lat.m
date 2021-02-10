function plot_ga_malr_diff_mon_lat(type, par)

    % plot inversion strength
    make_dirs_si_bl(type, par)

    prefix = make_prefix(type, par);
    prefix_proc = make_prefix_proc(type, par);
    plotdir = make_plotdir(type, par);

    load(sprintf('%s/si_bl_%g/ga_malr_diff_si_mon_lat.mat', prefix_proc, par.si_bl));
    load(sprintf('%s/si_bl_%g/ga_dalr_bl_diff_si_mon_lat.mat', prefix_proc, par.si_bl));

    [mesh_lat, mesh_mon] = meshgrid([1:12], lat);

%for l = {'lo', 'l', 'o'}; land = l{1};
    for l = {'lo'}; land = l{1};
        if strcmp(land, 'lo'); land_text = 'Land + Ocean';
        elseif strcmp(land, 'l'); land_text = 'Land';
        elseif strcmp(land, 'o'); land_text = 'Ocean';
        end

        % DALR BL mon x lat of diff
        figure(); clf; hold all; box on;
        cmp = colCog(20);
        colormap(flipud(cmp));
        contourf(mesh_lat, mesh_mon, ga_dalr_bl_diff.(land), -100:10:100, 'linecolor', 'none');
        [C,h]=contour(mesh_lat, mesh_mon, ga_dalr_bl_diff.(land), -100:10:100, '-w');
        contour(mesh_lat, mesh_mon, ga_dalr_bl_diff.(land), [0 0], 'color', 0.75*[1 1 1]);
        clabel(C, h, -100:20:100, 'fontsize', 6, 'interpreter', 'latex', 'color', 0.75*[1 1 1]);
        caxis([-100 100]);
        xlabel('Month'); ylabel('Latitude (deg)');
        make_title_type(type, par);
        cb = colorbar('limits', [0 100], 'ytick', [0:10:100], 'location', 'eastoutside');
        cb.TickLabelInterpreter = 'latex'; cb.Label.Interpreter = 'latex';
        ylabel(cb, sprintf('$\\left[(\\Gamma_d - \\Gamma)/\\Gamma_d\\right]_{1.0-%g}$ (\\%%)', par.si_bl));
        set(gca, 'xlim', [1 12], 'xtick', [1:12], 'xticklabels', par.monlabelnh, 'ylim', [-90 90], 'ytick', [-90:30:90], 'yminortick', 'on', 'tickdir', 'out');
        print(sprintf('%s/ga_malr_diff/si_bl_%g/%s/ga_dalr_bl_diff_mon_lat', plotdir, par.si_bl, land), '-dpng', '-r300');
        close;

        % ANOT DALR BL mon x lat of diff
        figure(); clf; hold all; box on;
        cmp = colCog(20);
        colormap(flipud(cmp));
        contourf(mesh_lat, mesh_mon, ga_dalr_bl_diff.(land), -100:10:100, 'linecolor', 'none');
        contour(mesh_lat, mesh_mon, ga_dalr_bl_diff.(land), par.ga_bl_thresh*[1 1], 'color', 'r', 'linewidth', 2);
        [C,h]=contour(mesh_lat, mesh_mon, ga_dalr_bl_diff.(land), -100:10:100, '-w');
        contour(mesh_lat, mesh_mon, ga_dalr_bl_diff.(land), [0 0], 'color', 0.75*[1 1 1]);
        clabel(C, h, -100:20:100, 'fontsize', 6, 'interpreter', 'latex', 'color', 0.75*[1 1 1]);
        caxis([-100 100]);
        xlabel('Month'); ylabel('Latitude (deg)');
        make_title_type(type, par);
        cb = colorbar('limits', [0 100], 'ytick', [0:10:100], 'location', 'eastoutside');
        cb.TickLabelInterpreter = 'latex'; cb.Label.Interpreter = 'latex';
        ylabel(cb, sprintf('$\\left[(\\Gamma_d - \\Gamma)/\\Gamma_d\\right]_{1.0-%g}$ (\\%%)', par.si_bl));
        set(gca, 'xlim', [1 12], 'xtick', [1:12], 'xticklabel', par.monlabelnh, 'ylim', [-90 90], 'ytick', [-90:30:90], 'yminortick', 'on', 'tickdir', 'out');
        print(sprintf('%s/ga_malr_diff/si_bl_%g/%s/ga_dalr_bl_diff_mon_lat_anot', plotdir, par.si_bl, land), '-dpng', '-r300');
        close;
        
        % mon x lat of diff
        figure(); clf; hold all; box on;
        cmp = colCog(12);
        colormap(flipud(cmp));
        contourf(mesh_lat, mesh_mon, ga_malr_diff.(land), -100:10:100, 'linecolor', 'none');
        [C,h]=contour(mesh_lat, mesh_mon, ga_malr_diff.(land), -100:10:100, '-w');
        contour(mesh_lat, mesh_mon, ga_malr_diff.(land), [0 0], 'color', 0.75*[1 1 1]);
        clabel(C, h, -100:10:100, 'fontsize', 6, 'interpreter', 'latex', 'color', 0.75*[1 1 1]);
        caxis([-60 60]);
        xlabel('Month'); ylabel('Latitude (deg)');
        make_title_type(type, par);
        cb = colorbar('limits', [-40 60], 'ytick', [-40:10:60], 'location', 'eastoutside');
        cb.TickLabelInterpreter = 'latex'; cb.Label.Interpreter = 'latex';
        ylabel(cb, sprintf('$\\left[(\\Gamma_m - \\Gamma)/\\Gamma_m\\right]_{%g-%g}$ (\\%%)', par.si_bl, par.si_up));
        set(gca, 'xlim', [1 12], 'xtick', [1:12], 'xticklabels', par.monlabelnh, 'ylim', [-90 90], 'ytick', [-90:30:90], 'yminortick', 'on', 'tickdir', 'out');
        print(sprintf('%s/ga_malr_diff/si_bl_%g/%s/ga_malr_diff_mon_lat', plotdir, par.si_bl, land), '-dpng', '-r300');
        close;

        % ANOT mon x lat of diff
        figure(); clf; hold all; box on;
        cmp = colCog(12);
        colormap(flipud(cmp));
        contourf(mesh_lat, mesh_mon, ga_malr_diff.(land), -100:10:100, 'linecolor', 'none');
        contour(mesh_lat, mesh_mon, ga_malr_diff.(land), par.ga_thresh*[1 1], 'color', 'r', 'linewidth', 2);
        [C,h]=contour(mesh_lat, mesh_mon, ga_malr_diff.(land), -100:10:100, '-w');
        contour(mesh_lat, mesh_mon, ga_malr_diff.(land), [0 0], 'color', 0.75*[1 1 1]);
        clabel(C, h, -100:10:100, 'fontsize', 6, 'interpreter', 'latex', 'color', 0.75*[1 1 1]);
        caxis([-60 60]);
        xlabel('Month'); ylabel('Latitude (deg)');
        make_title_type(type, par);
        cb = colorbar('limits', [-40 60], 'ytick', [-40:10:60], 'location', 'eastoutside');
        cb.TickLabelInterpreter = 'latex'; cb.Label.Interpreter = 'latex';
        ylabel(cb, sprintf('$\\left[(\\Gamma_m - \\Gamma)/\\Gamma_m\\right]_{%g-%g}$ (\\%%)', par.si_bl, par.si_up));
        set(gca, 'xlim', [1 12], 'xtick', [1:12], 'xticklabels', par.monlabelnh, 'ylim', [-90 90], 'ytick', [-90:30:90], 'yminortick', 'on', 'tickdir', 'out');
        print(sprintf('%s/ga_malr_diff/si_bl_%g/%s/ga_malr_diff_mon_lat_anot', plotdir, par.si_bl, land), '-dpng', '-r300');
        close;

    end
    
end
