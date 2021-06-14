function proc_ga_frac_midlatitude(type, par)
    make_dirs(type, par)

    prefix = make_prefix(type, par);
    prefix_proc = make_prefix_proc(type, par);
    plotdir = make_plotdir(type, par);

    load(sprintf('%s/grid.mat', prefix)); % read grid data
    
    lat_bound_list = [-10 10];
    lat_center = 50;

    [mesh_si, mesh_mon] = meshgrid(1:12, grid.dim3.si);

    for lb = 1:length(lat_bound_list); par.lat_bound = lat_bound_list(lb);
    
        [lat, clat, clat_mon, par] = make_midlatitude_lat(lat_center, par);

        load(sprintf('%s/ga_frac_midlatitude_lat_%g_to_%g', prefix_proc, par.lat_center-par.lat_bound, par.lat_center+par.lat_bound));
        foldername = sprintf('0_midlatitude_lat_%g_to_%g', par.lat_center-par.lat_bound, par.lat_center+par.lat_bound);
                

        for l = par.land_list; land = l{1};
            if strcmp(land, 'lo'); land_text = 'L+O';
            elseif strcmp(land, 'l'); land_text = 'L';
            elseif strcmp(land, 'o'); land_text = 'O';
            end
            if strcmp(type, 'echam');
                if strcmp(par.echam.clim, '20170908')
                    land_text = 'Snowball';
                else
                    land_text = par.echam.(par.echam.clim);
                end
            end;

            % min(min(ga_frac_lat.(land)))
            % max(max(ga_frac_lat.(land)))
            figure(); clf; hold all; box on;
            cmp = colCog(40);
            colormap(cmp);
            contourf(mesh_si, mesh_mon, 100 * circshift(ga_frac_lat.(land), par.shiftby, 1)', [-300:5:200], 'linecolor', 'w', 'linewidth', 0.1);
            contour(mesh_si, mesh_mon, 100 * circshift(ga_frac_lat.(land), par.shiftby, 1)', [0 0], 'linecolor', 0.75*[1 1 1], 'linewidth', 0.1);
            [C,h] = contour(mesh_si, mesh_mon, 100 * circshift(ga_frac_lat.(land), par.shiftby, 1)', [-300:10:200], 'linecolor', 'w', 'linewidth', 0.1);
            clabel(C, h, 'fontsize', 6, 'interpreter', 'latex', 'color', 0.75*[1 1 1]);
            make_title_type_lat(type, par.lat_center-par.lat_bound, par.lat_center+par.lat_bound, par);
            caxis([-100 100]);
            cb = colorbar('limits', [-30 100], 'ytick', [-30:10:100], 'location', 'eastoutside');
            cb.TickLabelInterpreter = 'latex'; cb.Label.Interpreter = 'latex';
            ylabel(cb, '$(\Gamma_m - \Gamma)/\Gamma_m$ (\%)');
            ylabel('$\sigma$ (unitless)');
            folder = sprintf('%s/ga_frac/%s/%s', plotdir, land, foldername);
            set(gca, 'ydir', 'reverse', 'ylim', [0.3 1], 'ytick', [0.3:0.1:1], 'xlim', [1 12], 'xtick', [1:12], 'xticklabels', par.monlabel, 'yminortick', 'on', 'tickdir', 'out');
            if ~exist(folder, 'dir')
                mkdir(folder);
            end
            print(sprintf('%s/ga_frac_mon_lev', folder), '-dpng', '-r300');
            close;

            figure(); clf; hold all; box on;
            cmp = colCog(40);
            colormap(cmp);
            contourf(mesh_si, mesh_mon, 100 * circshift(ga_frac_lat.(land), par.shiftby, 1)', [-300:5:200], 'linecolor', 'w', 'linewidth', 0.1);
            contour(mesh_si, mesh_mon, 100 * circshift(ga_frac_lat.(land), par.shiftby, 1)', [0 0], 'linecolor', 0.75*[1 1 1], 'linewidth', 0.1);
            [C,h] = contour(mesh_si, mesh_mon, 100 * circshift(ga_frac_lat.(land), par.shiftby, 1)', [-300:10:200], 'linecolor', 'w', 'linewidth', 0.1);
            clabel(C, h, 'fontsize', 6, 'interpreter', 'latex', 'color', 0.75*[1 1 1]);
            make_title_type_lat(type, par.lat_center-par.lat_bound, par.lat_center+par.lat_bound, par);
            caxis([-100 100]);
            ylabel('$\sigma$ (unitless)');
            folder = sprintf('%s/ga_frac/%s/%s', plotdir, land, foldername);
            set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos)
            set(gca, 'ydir', 'reverse', 'ylim', [0.3 1], 'ytick', [0.3:0.1:1], 'xlim', [1 12], 'xtick', [1:12], 'xticklabels', par.monlabel, 'yminortick', 'on', 'tickdir', 'out');
            if ~exist(folder, 'dir')
                mkdir(folder);
            end
            print(sprintf('%s/ga_frac_mon_lev_no_cb', folder), '-dpng', '-r300');
            close;

        end % for land

    end % lat bound

end % function
