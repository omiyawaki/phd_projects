function plot_flux(type, par)
    make_dirs(type, par)

    prefix = make_prefix(type, par);
    prefix_proc = make_prefix_proc(type, par);
    plotdir = make_plotdir(type, par);

    load(sprintf('%s/grid.mat', prefix)); % read grid data
    load(sprintf('%s/flux_z.mat', prefix_proc)); % load lat x mon RCAE data
    load(sprintf('%s/flux_t.mat', prefix_proc)); % load lat x lon RCAE data
    % load(sprintf('%s/%s/ga_malr_diff_si_mon_lat.mat', prefix_proc, par.lat_interp)); ga_malr_diff_orig = ga_malr_diff; clear ga_diff;
    % load(sprintf('%s/%s/ga_dalr_bl_diff_si_mon_lat.mat', prefix_proc, par.lat_interp)); ga_dalr_bl_diff_orig = ga_dalr_bl_diff; clear ga_bl_diff;
    % if strcmp(type, 'gcm') & ~any(strcmp(par.gcm.clim, {'hist-pi'})) & ~contains(par.model, 'mmm')
    %     load(sprintf('%s/sftlf.mat', prefix)); % read land fraction data
    % else
        landdata = load('/project2/tas1/miyawaki/matlab/landmask/land_mask.mat');
        par.land = landdata.land_mask; par.landlat = landdata.landlat; par.landlon = landdata.landlon;
    % end

    for l = par.land_list; land = l{1};
        if strcmp(land, 'lo'); land_text = 'Land + Ocean';
        elseif strcmp(land, 'l'); land_text = 'Land';
        elseif strcmp(land, 'o'); land_text = 'Ocean';
        end

        % ga_malr_diff = ga_malr_diff_orig.(land);
        % ga_dalr_bl_diff = ga_dalr_bl_diff_orig.(land);

        [mesh_lat, mesh_mon] = meshgrid(1:12, lat);
        f_vec = assign_fw(type, par);
        for f = f_vec; fw = f{1};

            % lat x mon dependence of RCE and RAE
            if strcmp(fw, 'dse'); var_text = '$\partial_t s + \nabla \cdot F_s$';
            elseif strcmp(fw, 'mse'); var_text = '$\nabla \cdot F_s$';
            else var_text = '$\partial_t h + \nabla \cdot F_m$'; end
            figure(); clf; hold all; box on;
            cmp = colCog(30);
            colormap(cmp);
            contourf(mesh_lat, mesh_mon, flux_z.(land).res.(fw), [-300 -150:10:150 300], 'linecolor', 'w', 'linewidth', 0.1);
            contour(mesh_lat, mesh_mon, flux_z.(land).res.(fw), [0 0], 'color', 0.75*[1 1 1]);
            make_title_type(type, par);
            caxis([-150 150]);
            cb = colorbar('limits', [-150 150], 'ytick', [-140:20:140], 'location', 'eastoutside');
            cb.TickLabelInterpreter = 'latex'; cb.Label.Interpreter = 'latex';
            ylabel(cb, sprintf('%s (Wm$^{-2}$)', var_text));
            ylabel('Latitude (deg)');
            set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_verywide);
            set(gca, 'xlim', [1 12], 'xtick', [1:12], 'xticklabels', par.monlabelnh, 'ylim', [-90 90], 'ytick', [-90:30:90], 'yminortick', 'on', 'tickdir', 'out');
            print(sprintf('%s/flux/%s/%s/0_div_mon_lat', plotdir, fw, land), '-dpng', '-r300');
            close;

            if any(strcmp(type, {'era5', 'era5c', 'erai', 'merra2c', 'jra55', 'gcm'}))

                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % UNCOMMENT FOR TEND
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                % var_text = '$\partial_t h$';
                % figure(); clf; hold all; box on;
                % cmp = colCog(100);
                % colormap(cmp);
                % contourf(mesh_lat, mesh_mon, flux_z.(land).tend, [-300 -50:1:50 300], 'linecolor', 'w', 'linewidth', 0.1);
                % contour(mesh_lat, mesh_mon, flux_z.(land).tend, [0 0], 'color', 0.75*[1 1 1]);
                % make_title_type(type, par);
                % caxis([-50 50]);
                % cb = colorbar('limits', [-50 50], 'ytick', [-50:10:50], 'location', 'eastoutside');
                % cb.TickLabelInterpreter = 'latex'; cb.Label.Interpreter = 'latex';
                % ylabel(cb, sprintf('%s (Wm$^{-2}$)', var_text));
                % ylabel('Latitude (deg)');
                % set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_verywide);
                % set(gca, 'xlim', [1 12], 'xtick', [1:12], 'xticklabels', par.monlabelnh, 'ylim', [-90 90], 'ytick', [-90:30:90], 'yminortick', 'on', 'tickdir', 'out');
                % print(sprintf('%s/flux/%s/%s/0_tend_mon_lat', plotdir, fw, land), '-dpng', '-r300');
                % close;

                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % UNCOMMENT FOR TEND
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                % var_text = '$\nabla \cdot F_m$';
                % figure(); clf; hold all; box on;
                % cmp = colCog(60);
                % colormap(cmp);
                % contourf(mesh_lat, mesh_mon, flux_z.(land).divfm, [-300 -150:5:150 300], 'linecolor', 'w', 'linewidth', 0.1);
                % contour(mesh_lat, mesh_mon, flux_z.(land).divfm, [0 0], 'color', 0.75*[1 1 1]);
                % make_title_type(type, par);
                % caxis([-150 150]);
                % cb = colorbar('limits', [-150 150], 'ytick', [-140:20:140], 'location', 'eastoutside');
                % cb.TickLabelInterpreter = 'latex'; cb.Label.Interpreter = 'latex';
                % ylabel(cb, sprintf('%s (Wm$^{-2}$)', var_text));
                % ylabel('Latitude (deg)');
                % set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_verywide);
                % set(gca, 'xlim', [1 12], 'xtick', [1:12], 'xticklabels', par.monlabelnh, 'ylim', [-90 90], 'ytick', [-90:30:90], 'yminortick', 'on', 'tickdir', 'out');
                % print(sprintf('%s/flux/%s/%s/0_divfm_mon_lat', plotdir, fw, land), '-dpng', '-r300');
                % close;
            end

            % lat x mon dependence of RCE and RAE
            var_text = '$R_a$';
            figure(); clf; hold all; box on;
            if strcmp(type, 'echam')
                cmp = colCog(60);
            else
                cmp = colCog(30);
            end
            colormap(cmp);
            if strcmp(type, 'echam')
                contourf(mesh_lat, mesh_mon, flux_z.(land).ra.(fw), [-300 -150:5:150 300], 'linecolor', 'none');
            else
                contourf(mesh_lat, mesh_mon, flux_z.(land).ra.(fw), [-300 -150:10:150 300], 'linecolor', 'w', 'linewidth', 0.1);
            end
            contour(mesh_lat, mesh_mon, flux_z.(land).ra.(fw), [0 0], 'color', 0.75*[1 1 1]);
            make_title_type(type, par);
            caxis([-150 150]);
            if strcmp(type, 'echam')
                cb = colorbar('limits', [-50 50], 'ytick', [-50:10:50], 'location', 'eastoutside');
            else
                cb = colorbar('limits', [-150 150], 'ytick', [-150:25:150], 'location', 'eastoutside');
            end
            cb.TickLabelInterpreter = 'latex'; cb.Label.Interpreter = 'latex';
            ylabel(cb, sprintf('%s (Wm$^{-2}$)', var_text));
            ylabel('Latitude (deg)');
            set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_verywide);
            set(gca, 'xlim', [1 12], 'xtick', [1:12], 'xticklabels', par.monlabelnh, 'ylim', [-90 90], 'ytick', [-90:30:90], 'yminortick', 'on', 'tickdir', 'out');
            print(sprintf('%s/flux/%s/%s/0_ra_mon_lat', plotdir, fw, land), '-dpng', '-r300');
            close;

            % lat x mon dependence of RCE and RAE
            [lh, sh] = rename_stf(type, flux_z, land);
            var_text = '$\mathrm{LH}$';
            figure(); clf; hold all; box on;
            cmp = colCog(60);
            colormap(cmp);
            contourf(mesh_lat, mesh_mon, lh, [-300 -150:5:150 300], 'linecolor', 'w', 'linewidth', 0.1);
            % contour(mesh_lat, mesh_mon, lh, [0 0], 'color', 0.75*[1 1 1]);
            make_title_type(type, par);
            caxis([-150 150]);
            cb = colorbar('limits', [-150 150], 'ytick', [-140:20:140], 'location', 'eastoutside');
            cb.TickLabelInterpreter = 'latex'; cb.Label.Interpreter = 'latex';
            ylabel(cb, sprintf('%s (Wm$^{-2}$)', var_text));
            ylabel('Latitude (deg)');
            set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_verywide);
            set(gca, 'xlim', [1 12], 'xtick', [1:12], 'xticklabels', par.monlabelnh, 'ylim', [-90 90], 'ytick', [-90:30:90], 'yminortick', 'on', 'tickdir', 'out');
            print(sprintf('%s/flux/%s/%s/0_lh_mon_lat',plotdir, fw, land), '-dpng', '-r300');
            close;

            % lat x mon dependence of RCE and RAE
            var_text = '$\mathrm{SH}$';
            figure(); clf; hold all; box on;
            cmp = colCog(40);
            colormap(cmp);
            contourf(mesh_lat, mesh_mon, sh, [-300 -50:2.5:50 300], 'linecolor', 'w', 'linewidth', 0.1);
            contour(mesh_lat, mesh_mon, sh, [0 0], 'color', 0.75*[1 1 1]);
            make_title_type(type, par);
            caxis([-50 50]);
            cb = colorbar('limits', [-50 50], 'ytick', [-50:10:50], 'location', 'eastoutside');
            cb.TickLabelInterpreter = 'latex'; cb.Label.Interpreter = 'latex';
            ylabel(cb, sprintf('%s (Wm$^{-2}$)', var_text));
            ylabel('Latitude (deg)');
            set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_verywide);
            set(gca, 'xlim', [1 12], 'xtick', [1:12], 'xticklabels', par.monlabelnh, 'ylim', [-90 90], 'ytick', [-90:30:90], 'yminortick', 'on', 'tickdir', 'out');
            print(sprintf('%s/flux/%s/%s/0_sh_mon_lat',plotdir, fw, land), '-dpng', '-r300');
            close;

            if strcmp(type, 'gcm') & strcmp(par.gcm.clim, 'hist-pi')
                % R1 lat x mon dependence of RCE and RAE
                var_text = '$R_1^{\mathrm{historical}}-R_1^{\mathrm{piControl}}$';
                figure(); clf; hold all; box on;
                cmp = colCog(10);
                colormap(flipud(cmp));
                contourf(mesh_lat, mesh_mon, flux_z.(land).r1.(fw), [-0.1 -0.05:0.01:0.05 0.1], 'linecolor', 'w');
                [C, h] = contour(mesh_lat, mesh_mon, flux_z.(land).r1.(fw), [par.ep par.ga], 'linecolor', 'w', 'linewidth', 0.25);
                clabel(C, h, 'fontsize', 6, 'interpreter', 'latex', 'color', 'k');
                if contains(par.model, 'mmm')
                    title(sprintf('CMIP5 %s', par.gcm.clim));
                else
                    title(sprintf('%s', par.model));
                end
                caxis([-0.05 0.05]);
                cb = colorbar('limits', [-0.05 0.05], 'ytick', [-0.05:0.01:0.05], 'location', 'eastoutside');
                cb.TickLabelInterpreter = 'latex'; cb.Label.Interpreter = 'latex';
                ylabel(cb, sprintf('%s (unitless)', var_text));
                ylabel('Latitude (deg)');
                set(gca, 'xlim', [1 12], 'xtick', [1:12], 'xticklabels', par.monlabelnh, 'ylim', [-90 90], 'ytick', [-90:30:90], 'yminortick', 'on', 'tickdir', 'out');
                print(sprintf('%s/flux/%s/%s/0_r1_mon_lat', plotdir, fw, land), '-dpng', '-r300');
                close;

            else
                % R1 lat x mon dependence of RCE and RAE
                var_text = '$R_1$';
                figure(); clf; hold all; box on;
                cmp = colCog(10);
                colormap(flipud(cmp));
                contourf(mesh_lat, mesh_mon, flux_z.(land).r1.(fw), [-16 -8 -4 -2:0.2:2 4 8 16], 'linecolor', 'w');
                contour(mesh_lat, mesh_mon, flux_z.(land).r1.(fw), [0 0], 'color', 0.75*[1 1 1]);
                contour(mesh_lat, mesh_mon, flux_z.(land).r1.(fw), par.ep*[1 1], 'linecolor', par.orange, 'linewidth', 1.5);
                contour(mesh_lat, mesh_mon, flux_z.(land).r1.(fw), par.ga*[1 1], 'linecolor', par.cyan, 'linewidth', 1.5);
                [C, h] = contour(mesh_lat, mesh_mon, flux_z.(land).r1.(fw), [par.ep par.ga], 'linecolor', 'w', 'linewidth', 0.25);
                clabel(C, h, 'fontsize', 6, 'interpreter', 'latex', 'color', 'k');
                make_title_type(type, par);
                caxis([-1 1]);
                cb = colorbar('limits', [-1 1], 'ytick', [-1:0.2:1], 'location', 'eastoutside');
                cb.TickLabelInterpreter = 'latex'; cb.Label.Interpreter = 'latex';
                ylabel(cb, sprintf('%s (unitless)', var_text));
                ylabel('Latitude (deg)');
                set(gca, 'xlim', [1 12], 'xtick', [1:12], 'xticklabels', par.monlabelnh, 'ylim', [-90 90], 'ytick', [-90:30:90], 'yminortick', 'on', 'tickdir', 'out');
                print(sprintf('%s/flux/%s/%s/0_r1_mon_lat', plotdir, fw, land), '-dpng', '-r300');
                close;

                % R1z lat x mon dependence of RCE and RAE
                if ~any(strcmp(fw, {'mse_old', 'dse_old'}))
                    var_text = '$R_1^*$';
                else
                    var_text = '$R_1$';
                end
                figure(); clf; hold all; box on;
                cmp = colCog(20);
                colormap(flipud(cmp));
                r1z = flux_z.(land).res.(fw)./flux_z.(land).ra.(fw);
                contourf(mesh_lat, mesh_mon, r1z, [-16 -8 -4 -2:0.1:2 4 8 16], 'linecolor', 'w', 'linewidth', 0.1);
                contour(mesh_lat, mesh_mon,  flux_z.(land).res.(fw)./flux_z.(land).ra.(fw), [0 0], 'color', 0.75*[1 1 1]);
                [C, h] = contour(mesh_lat, mesh_mon, flux_z.(land).res.(fw)./flux_z.(land).ra.(fw), [-16 -8 -4 -2:0.2:2 4 8 16], 'linecolor', 'w', 'linewidth', 0.1);
                contour(mesh_lat, mesh_mon,  flux_z.(land).res.(fw)./flux_z.(land).ra.(fw), [-16 -8 -4 -2:0.1:2 4 8 16], 'color', 'w', 'linewidth', 0.1);
                % [C, h] = contour(mesh_lat, mesh_mon, flux_z.(land).res.(fw)./flux_z.(land).ra.(fw), [par.ep par.ga], 'linecolor', 'w', 'linewidth', 0.25);
                contour(mesh_lat, mesh_mon,  flux_z.(land).res.(fw)./flux_z.(land).ra.(fw), par.ep*[1 1], 'linecolor', par.orange, 'linewidth', 2);
                contour(mesh_lat, mesh_mon,  flux_z.(land).res.(fw)./flux_z.(land).ra.(fw), par.ga*[1 1], 'linecolor', par.cyan, 'linewidth', 2);
                clabel(C, h, 'fontsize', 6, 'interpreter', 'latex', 'color', 0.75*[1 1 1]);
                make_title_type(type, par);
                caxis([-1 1]);
                cb = colorbar('limits', [-1 1], 'ytick', [-1:0.2:1], 'location', 'eastoutside');
                cb.TickLabelInterpreter = 'latex'; cb.Label.Interpreter = 'latex';
                ylabel(cb, sprintf('%s (unitless)', var_text));
                ylabel('Latitude (deg)');
                set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_verywide);
                set(gca, 'xlim', [1 12], 'xtick', [1:12], 'xticklabels', par.monlabelnh, 'ylim', [-90 90], 'ytick', [-90:30:90], 'yminortick', 'on', 'tickdir', 'out');
                print(sprintf('%s/flux/%s/%s/0_r1z_mon_lat', plotdir, fw, land), '-dpng', '-r300');
                data = [mesh_lat(:) mesh_mon(:) r1z(:) ];
                save(sprintf('%s/flux/%s/%s/0_r1z_mon_lat.dat', plotdir, fw, land), 'data', '-ASCII')
                % if par.make_tikz
                %     matlab2tikz(sprintf('%s/flux/%s/%s/0_r1z_mon_lat.tex', plotdir, fw, land));
                % end
                close;

                % R1z lat x mon dependence of RCE and RAE
                if ~strcmp(fw, 'mse_old')
                    var_text = '$R_1^*$';
                else
                    var_text = '$R_1$';
                end
                figure(); clf; hold all; box on;
                cmp = colCog(10);
                colormap(flipud(cmp));
                contourf(mesh_lat, mesh_mon, flux_z.(land).res.(fw)./flux_z.(land).ra.(fw), par.ga*[1 1], 'edgecolor', 'none');
                % contourf(mesh_lat, mesh_mon, flux_z.(land).res.(fw)./flux_z.(land).ra.(fw), par.ep*[1 1], 'edgecolor', 'none');
                contour(mesh_lat, mesh_mon, flux_z.(land).res.(fw)./flux_z.(land).ra.(fw), [-16 -8 -4 -2:0.2:2 4 8 16], 'linecolor', 'k');
                % contour(mesh_lat, mesh_mon,  flux_z.(land).res.(fw)./flux_z.(land).ra.(fw), [0 0], 'color', 0.75*[1 1 1]);
                % contour(mesh_lat, mesh_mon,  flux_z.(land).res.(fw)./flux_z.(land).ra.(fw), par.ep*[1 1], 'linecolor', par.orange, 'linewidth', 1.5);
                % contour(mesh_lat, mesh_mon,  flux_z.(land).res.(fw)./flux_z.(land).ra.(fw), par.ga*[1 1], 'linecolor', par.cyan, 'linewidth', 1.5);
                % [C, h] = contour(mesh_lat, mesh_mon, flux_z.(land).res.(fw)./flux_z.(land).ra.(fw), [par.ep par.ga], 'linecolor', 'w', 'linewidth', 0.25);
                % clabel(C, h, 'fontsize', 6, 'interpreter', 'latex', 'color', 'k');
                make_title_type(type, par);
                % caxis([-1 1]);
                % cb = colorbar('limits', [-1 1], 'ytick', [-1:0.2:1], 'location', 'eastoutside');
                % cb.TickLabelInterpreter = 'latex'; cb.Label.Interpreter = 'latex';
                % ylabel(cb, sprintf('%s (unitless)', var_text));
                ylabel('Latitude (deg)');
                set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos);
                set(gca, 'xlim', [1 12], 'xtick', [1:12], 'xticklabels', par.monlabelnh, 'ylim', [-90 90], 'ytick', [-90:30:90], 'yminortick', 'on', 'tickdir', 'out');
                print(sprintf('%s/flux/%s/%s/0_r1z_mon_lat_alt', plotdir, fw, land), '-dpng', '-r300');
                close;

                % R2z lat x mon dependence of RCE and RAE
                if ~strcmp(fw, 'mse_old')
                    var_text = '$R_2^*$';
                else
                    var_text = '$R_2$';
                end
                figure(); clf; hold all; box on;
                cmp = colCog(10);
                colormap(flipud(cmp));
                contourf(mesh_lat, mesh_mon, flux_z.(land).stf.(fw)./flux_z.(land).ra.(fw), [-16 -8 -4 -2:0.2:2 4 8 16], 'linecolor', 'w');
                contour(mesh_lat, mesh_mon,  flux_z.(land).stf.(fw)./flux_z.(land).ra.(fw), [0 0], 'color', 0.75*[1 1 1]);
                contour(mesh_lat, mesh_mon,  flux_z.(land).stf.(fw)./flux_z.(land).ra.(fw), (par.ep-1)*[1 1], 'linecolor', par.orange, 'linewidth', 1.5);
                contour(mesh_lat, mesh_mon,  flux_z.(land).stf.(fw)./flux_z.(land).ra.(fw), (par.ga-1)*[1 1], 'linecolor', par.cyan, 'linewidth', 1.5);
                [C, h] = contour(mesh_lat, mesh_mon, flux_z.(land).stf.(fw)./flux_z.(land).ra.(fw), [par.ep par.ga], 'linecolor', 'w', 'linewidth', 0.25);
                clabel(C, h, 'fontsize', 6, 'interpreter', 'latex', 'color', 'k');
                make_title_type(type, par);
                caxis([-2 2]);
                cb = colorbar('limits', [-1-1 1-1], 'ytick', [-1-1:0.2:1-1], 'location', 'eastoutside');
                cb.TickLabelInterpreter = 'latex'; cb.Label.Interpreter = 'latex';
                ylabel(cb, sprintf('%s (unitless)', var_text));
                ylabel('Latitude (deg)');
                set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_verywide);
                set(gca, 'xlim', [1 12], 'xtick', [1:12], 'xticklabels', par.monlabelnh, 'ylim', [-90 90], 'ytick', [-90:30:90], 'yminortick', 'on', 'tickdir', 'out');
                print(sprintf('%s/flux/%s/%s/0_r2z_mon_lat', plotdir, fw, land), '-dpng', '-r300');
                close;

                % DEVIATION R1z lat x mon dependence of RCE and RAE
                var_text = '$\Delta R_1$';
                figure(); clf; hold all; box on;
                cmp = colCog(10);
                colormap(flipud(cmp));
                dev = flux_z.(land).res.(fw)./flux_z.(land).ra.(fw) - repmat(nanmean(flux_z.(land).res.(fw)./flux_z.(land).ra.(fw), 2), [1 12]);
                contourf(mesh_lat, mesh_mon, dev, 1/2*[-16 -8 -4 -2 -1:0.2:1 2 4 8 16], 'linecolor', 'w');
                contour(mesh_lat, mesh_mon,  dev, [0 0], 'color', 0.75*[1 1 1]);
                contour(mesh_lat, mesh_mon,  flux_z.(land).res.(fw)./flux_z.(land).ra.(fw), par.ep*[1 1], 'linecolor', par.orange, 'linewidth', 1.5);
                contour(mesh_lat, mesh_mon,  flux_z.(land).res.(fw)./flux_z.(land).ra.(fw), par.ga*[1 1], 'linecolor', par.cyan, 'linewidth', 1.5);
                [C, h] = contour(mesh_lat, mesh_mon, dev, [par.ep par.ga], 'linecolor', 'w', 'linewidth', 0.25);
                % clabel(C, h, 'fontsize', 6, 'interpreter', 'latex', 'color', 'k');
                make_title_type(type, par);
                caxis([-0.5 0.5]);
                cb = colorbar('limits', [-0.5 0.5], 'ytick', [-0.5:0.1:0.5], 'location', 'eastoutside');
                cb.TickLabelInterpreter = 'latex'; cb.Label.Interpreter = 'latex';
                ylabel(cb, sprintf('%s (unitless)', var_text));
                ylabel('Latitude (deg)');
                set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos);
                set(gca, 'xlim', [1 12], 'xtick', [1:12], 'xticklabels', par.monlabelnh, 'ylim', [-90 90], 'ytick', [-90:30:90], 'yminortick', 'on', 'tickdir', 'out');
                print(sprintf('%s/flux/%s/%s/0_dr1z_mon_lat', plotdir, fw, land), '-dpng', '-r300');
                close;

                % % GA MALR and GA_BL DALR contour lat x mon dependence of RCE and RAE
                % figure(); clf; hold all;
                % cmp = colCog(10);
                % colormap(flipud(cmp));
                % contourf(mesh_lat, mesh_mon, flux_z.(land).r1.(fw), [-16 -8 -4 -2 -1:0.2:1 2 4 8 16], 'linecolor', 'w');
                % contour(mesh_lat, mesh_mon, flux_z.(land).r1.(fw), [0 0], 'color', 0.75*[1 1 1]);
                % [C, h] = contour(mesh_lat, mesh_mon, flux_z.(land).r1.(fw), [par.ep par.ga], 'linecolor', 0.5*[1 1 1], 'linewidth', 1);
                % clabel(C, h, 'fontsize', 6, 'interpreter', 'latex', 'color', 0.5*[1 1 1]);
                % contour(mesh_lat, mesh_mon, ga_malr_diff, par.ga_thresh*[1 1], 'color', par.orange, 'linewidth', 2);
                % contour(mesh_lat, mesh_mon, ga_dalr_bl_diff, par.ga_bl_thresh*[1 1], 'color', par.blue, 'linewidth', 2);
                % if any(strcmp(type, {'era5', 'erai'})); title(sprintf('%s, $\\Gamma_m - \\Gamma / \\Gamma_m < %g \\%%$, $(\\Gamma_d - \\Gamma )/ \\Gamma_d > %g \\%%$', upper(type), par.ga_thresh, par.ga_bl_thresh));
                % elseif strcmp(type, 'merra2'); title(sprintf('%s, $\\Gamma_m - \\Gamma / \\Gamma_m < %g \\%%$, $(\\Gamma_d - \\Gamma )/ \\Gamma_d > %g \\%%$', upper(type), par.ga_thresh, par.ga_bl_thresh));
                % elseif strcmp(type, 'gcm'); title(sprintf('%s, $\\Gamma_m - \\Gamma / \\Gamma_m < %g \\%%$, $(\\Gamma_d - \\Gamma )/ \\Gamma_d > %g \\%%$', par.model, par.ga_thresh, par.ga_bl_thresh));
                % elseif any(strcmp(type, {'echam'})); title(sprintf('%s, $\\Gamma_m - \\Gamma / \\Gamma_m < %g \\%%$, $(\\Gamma_d - \\Gamma )/ \\Gamma_d > %g \\%%$', upper(type), par.ga_thresh, par.ga_bl_thresh)); end;
                % caxis([-1 1]);
                % cb = colorbar('limits', [-0.4 1], 'ytick', [-1:0.2:1], 'location', 'eastoutside');
                % cb.TickLabelInterpreter = 'latex'; cb.Label.Interpreter = 'latex';
                % ylabel(cb, sprintf('%s (unitless)', var_text));
                % xlabel('Month'); ylabel('Latitude (deg)');
                % set(gca, 'xlim', [1 12], 'xtick', [1:12], 'xticklabels', par.monlabel, 'ylim', [-90 90], 'ytick', [-90:30:90], 'yminortick', 'on', 'tickdir', 'out');
                % print(sprintf('%s/flux/%s/%s/0_r1_mon_lat_ga_malr_ga_dalr_bl_overlay', plotdir, fw, land), '-dpng', '-r300');
                % close;


                [mesh_ll_lat_2d, mesh_ll_lon_2d] = meshgrid(grid.dim2.lon, grid.dim2.lat);
                [mesh_ll_lat, mesh_ll_lon] = meshgrid(grid.dim3.lon, lat);
                for t = {'ann', 'djf', 'jja', 'mam', 'son'}; time = t{1};
                % for t = {'ann'}; time = t{1};
                    % lat x lon of RCE and RAE
                    if strcmp(fw, 'dse'); var_text = '$\partial_t s + \nabla \cdot F_s$ (Wm$^{-2}$)';
                    elseif strcmp(fw, 'mse'); var_text = '$\nabla \cdot F_s$ (Wm$^{-2}$)';
                    else var_text = '$\partial_t h + \nabla \cdot F_m$ (Wm$^{-2}$)'; end
                    figure(); clf; hold all;
                    cmp = colCog(40);
                    colormap(cmp);
                    contourf(mesh_ll_lat, mesh_ll_lon, flux_t.(land).(time).res.(fw)', [-200:10:200], 'linecolor', 'none');
                    caxis([-200 200]);
                    cb = colorbar('limits', [-200 60], 'ytick', [-200:20:60], 'location', 'eastoutside');
                    cb.TickLabelInterpreter = 'latex'; cb.Label.Interpreter = 'latex';
                    ylabel(cb, var_text);
                    % if strcmp(type, 'gcm') & ~any(strcmp(par.gcm.clim, {'hist-pi'})) & ~contains(par.model, 'mmm')
                    %     contour(mesh_ll_lat_2d, mesh_ll_lon_2d, sftlf', 0.5*[1 1], 'w', 'linewidth', 1);
                    %     contour(mesh_ll_lat_2d, mesh_ll_lon_2d, sftlf', 0.5*[1 1], 'k');
                    % else
                        contour(par.landlon+180, par.landlat, par.land, [1 1], 'k');
                    % end
                    make_title_type_time(type, time, par);
                    xlabel('Longitude (deg)'); ylabel('Latitude (deg)');
                    set(gca, 'xlim', [0 360], 'xtick', [0:60:360], 'ylim', [-90 90], 'ytick', [-90:30:90], 'yminortick', 'on', 'tickdir', 'out');
                    print(sprintf('%s/flux/%s/%s/%s/div_lat_lon', plotdir, fw, land, time), '-dpng', '-r300');
                    close;

                    if any(strcmp(type, {'era5c', 'era5', 'erai', 'merra2c', 'jra55', 'gcm'}))

                        % var_text = '$\nabla \cdot F_m$ (Wm$^{-2}$)';
                        % figure(); clf; hold all;
                        % cmp = colCog(40);
                        % colormap(cmp);
                        % contourf(mesh_ll_lat, mesh_ll_lon, flux_t.(land).(time).divfm', [-200:10:200], 'linecolor', 'none');
                        % caxis([-200 200]);
                        % cb = colorbar('limits', [-200 60], 'ytick', [-200:20:60], 'location', 'eastoutside');
                        % cb.TickLabelInterpreter = 'latex'; cb.Label.Interpreter = 'latex';
                        % ylabel(cb, var_text);
                        % if strcmp(type, 'gcm') & ~any(strcmp(par.gcm.clim, {'hist-pi'})) & ~contains(par.model, 'mmm')
                        %     contour(mesh_ll_lat, mesh_ll_lon, sftlf', 0.5*[1 1], 'w', 'linewidth', 1);
                        %     contour(mesh_ll_lat, mesh_ll_lon, sftlf', 0.5*[1 1], 'k');
                        % else
                        %     contour(par.landlon+180, par.landlat, par.land, [1 1], 'k');
                        % end
                        % make_title_type_time(type, time, par);
                        % xlabel('Longitude (deg)'); ylabel('Latitude (deg)');
                        % set(gca, 'xlim', [0 360], 'xtick', [0:60:360], 'ylim', [-90 90], 'ytick', [-90:30:90], 'yminortick', 'on', 'tickdir', 'out');
                        % print(sprintf('%s/flux/%s/%s/%s/divfm_lat_lon', plotdir, fw, land, time), '-dpng', '-r300');
                        % close;

                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        % UNCOMMENT FOR TEND
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                        % var_text = '$\partial_t h$ (Wm$^{-2}$)';
                        % figure(); clf; hold all;
                        % cmp = colCog(20);
                        % colormap(cmp);
                        % contourf(mesh_ll_lat, mesh_ll_lon, flux_t.(land).(time).tend', [-50:5:50], 'linecolor', 'none');
                        % caxis([-50 50]);
                        % cb = colorbar('limits', [-50 50], 'ytick', [-50:10:50], 'location', 'eastoutside');
                        % cb.TickLabelInterpreter = 'latex'; cb.Label.Interpreter = 'latex';
                        % ylabel(cb, var_text);
                        % % if strcmp(type, 'gcm') & ~any(strcmp(par.gcm.clim, {'hist-pi'})) & ~contains(par.model, 'mmm')
                        % %     contour(mesh_ll_lat_2d, mesh_ll_lon_2d, sftlf', 0.5*[1 1], 'w', 'linewidth', 1);
                        % %     contour(mesh_ll_lat_2d, mesh_ll_lon_2d, sftlf', 0.5*[1 1], 'k');
                        % % else
                        %     contour(par.landlon+180, par.landlat, par.land, [1 1], 'k');
                        % % end
                        % make_title_type_time(type, time, par);
                        % xlabel('Longitude (deg)'); ylabel('Latitude (deg)');
                        % set(gca, 'xlim', [0 360], 'xtick', [0:60:360], 'ylim', [-90 90], 'ytick', [-90:30:90], 'yminortick', 'on', 'tickdir', 'out');
                        % print(sprintf('%s/flux/%s/%s/%s/tend_lat_lon', plotdir, fw, land, time), '-dpng', '-r300');
                        % close;

                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        % UNCOMMENT FOR TEND
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                    end

                    var_text = '$R_1$ (unitless)';
                    figure(); clf; hold all;
                    cmp = colCog(20);
                    colormap(flipud(cmp));
                    contourf(mesh_ll_lat, mesh_ll_lon, flux_t.(land).(time).r1.(fw)', [-16 -8 -4 -2 -1:0.1:1 2 4 8 16], 'linecolor', 'none');
                    % if strcmp(type, 'gcm') & ~any(strcmp(par.gcm.clim, {'hist-pi'})) & ~contains(par.model, 'mmm')
                    %     contour(mesh_ll_lat_2d, mesh_ll_lon_2d, sftlf', 0.5*[1 1], 'w', 'linewidth', 1);
                    %     contour(mesh_ll_lat_2d, mesh_ll_lon_2d, sftlf', 0.5*[1 1], 'k');
                    % else
                        contour(par.landlon+180, par.landlat, par.land, [1 1], 'k');
                    % end
                    caxis([-1 1]);
                    cb = colorbar('limits', [-1 1], 'ytick', [-1:0.2:1], 'location', 'eastoutside');
                    cb.TickLabelInterpreter = 'latex'; cb.Label.Interpreter = 'latex';
                    ylabel(cb, var_text);
                    make_title_type_time(type, time, par);
                    xlabel('Longitude (deg)'); ylabel('Latitude (deg)');
                    set(gca, 'xlim', [0 360], 'xtick', [0:60:360], 'ylim', [-90 90], 'ytick', [-90:30:90], 'yminortick', 'on', 'tickdir', 'out');
                    print(sprintf('%s/flux/%s/%s/%s/r1_lat_lon', plotdir, fw, land, time), '-dpng', '-r300');
                    close;

                    var_text = '$R_a$ (Wm$^{-2}$)';
                    figure(); clf; hold all;
                    cmp = colCog(40);
                    colormap(cmp);
                    contourf(mesh_ll_lat, mesh_ll_lon, flux_t.(land).(time).ra.(fw)', [-200:10:0], 'linecolor', 'none');
                    % if strcmp(type, 'gcm') & ~any(strcmp(par.gcm.clim, {'hist-pi'})) & ~contains(par.model, 'mmm')
                    %     contour(mesh_ll_lat_2d, mesh_ll_lon_2d, sftlf', 0.5*[1 1], 'w', 'linewidth', 1);
                    %     contour(mesh_ll_lat_2d, mesh_ll_lon_2d, sftlf', 0.5*[1 1], 'k');
                    % else
                        contour(par.landlon+180, par.landlat, par.land, [1 1], 'k');
                    % end
                    caxis([-200 200]);
                    cb = colorbar('limits', [-160 0], 'ytick', [-160:20:0], 'location', 'eastoutside');
                    cb.TickLabelInterpreter = 'latex'; cb.Label.Interpreter = 'latex';
                    ylabel(cb, var_text);
                    make_title_type_time(type, time, par);
                    xlabel('Longitude (deg)'); ylabel('Latitude (deg)');
                    set(gca, 'xlim', [0 360], 'xtick', [0:60:360], 'ylim', [-90 90], 'ytick', [-90:30:90], 'yminortick', 'on', 'tickdir', 'out');
                    print(sprintf('%s/flux/%s/%s/%s/ra_lat_lon', plotdir, fw, land, time), '-dpng', '-r300');
                    close;

                    [lh, sh] = rename_stf(type, flux_t, land, time);
                    var_text = 'LH (Wm$^{-2}$)';
                    figure(); clf; hold all;
                    cmp = colCog(80);
                    colormap(cmp);
                    contourf(mesh_ll_lat, mesh_ll_lon, lh', [-200:5:200], 'linecolor', 'none');
                    % if strcmp(type, 'gcm') & ~any(strcmp(par.gcm.clim, {'hist-pi'})) & ~contains(par.model, 'mmm')
                    %     contour(mesh_ll_lat_2d, mesh_ll_lon_2d, sftlf', 0.5*[1 1], 'w', 'linewidth', 1);
                    %     contour(mesh_ll_lat_2d, mesh_ll_lon_2d, sftlf', 0.5*[1 1], 'k');
                    % else
                        contour(par.landlon+180, par.landlat, par.land, [1 1], 'k');
                    % end
                    caxis([-200 200]);
                    cb = colorbar('limits', [0 160], 'ytick', [0:20:160], 'location', 'eastoutside');
                    cb.TickLabelInterpreter = 'latex'; cb.Label.Interpreter = 'latex';
                    ylabel(cb, var_text);
                    make_title_type_time(type, time, par);
                    xlabel('Longitude (deg)'); ylabel('Latitude (deg)');
                    set(gca, 'xlim', [0 360], 'xtick', [0:60:360], 'ylim', [-90 90], 'ytick', [-90:30:90], 'yminortick', 'on', 'tickdir', 'out');
                    print(sprintf('%s/flux/%s/%s/%s/lh_lat_lon', plotdir, fw, land, time), '-dpng', '-r300');
                    close;

                    var_text = 'SH (Wm$^{-2}$)';
                    figure(); clf; hold all;
                    cmp = colCog(40);
                    colormap(cmp);
                    contourf(mesh_ll_lat, mesh_ll_lon, sh', [-100:5:100], 'linecolor', 'none');
                    % if strcmp(type, 'gcm') & ~any(strcmp(par.gcm.clim, {'hist-pi'})) & ~contains(par.model, 'mmm')
                    %     contour(mesh_ll_lat_2d, mesh_ll_lon_2d, sftlf', 0.5*[1 1], 'w', 'linewidth', 1);
                    %     contour(mesh_ll_lat_2d, mesh_ll_lon_2d, sftlf', 0.5*[1 1], 'k');
                    % else
                        contour(par.landlon+180, par.landlat, par.land, [1 1], 'k');
                    % end
                    caxis([-100 100]);
                    cb = colorbar('limits', [-100 100], 'ytick', [-100:20:100], 'location', 'eastoutside');
                    cb.TickLabelInterpreter = 'latex'; cb.Label.Interpreter = 'latex';
                    ylabel(cb, var_text);
                    make_title_type_time(type, time, par);
                    xlabel('Longitude (deg)'); ylabel('Latitude (deg)');
                    set(gca, 'xlim', [0 360], 'xtick', [0:60:360], 'ylim', [-90 90], 'ytick', [-90:30:90], 'yminortick', 'on', 'tickdir', 'out');
                    print(sprintf('%s/flux/%s/%s/%s/sh_lat_lon', plotdir, fw, land, time), '-dpng', '-r300');
                    close;

                end % for time

            end % if hist-pi


        end % for mse dse
    end % for land
end
