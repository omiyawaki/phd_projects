function plot_dflux(type, par)

    if ~strcmp(type, 'gcm')
        error('This analysis only applies for CMIP5 historical')
    end

    partitle = par;
    partitle.gcm.clim = 'RCP8.5$-$historical';

    make_dirs(type, par)

    prefix = make_prefix(type, par);
    prefix_proc = make_prefix_proc(type, par);
    plotdir = make_plotdir(type, par);

    load(sprintf('%s/grid.mat', prefix)); % read grid data

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Load temperature response amplification
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    load(sprintf('%s/dtempsi_mon_lat.mat', prefix_proc)); % load lat x mon RCAE_ALT data

    dtemp = permute(dtempsi.lo, [3 1 2]);
    dtemp_u = squeeze(interp1(grid.dim3.si, dtemp, 0.3)); % upper trop dtemp
    dtemp_l = squeeze(interp1(grid.dim3.si, dtemp, 1)); % lower trop dtemp
    dtemp_r = dtemp_u./dtemp_l; % upper to lower temperature response ratio

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Load R1 info
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    tmp=load(sprintf('%s/flux_z.mat', prefix_proc)); flux_z1 = tmp.flux_z; lat = tmp.lat; clear tmp; % load lat x mon RCAE data
    tmp=load(sprintf('%s/flux_t.mat', prefix_proc)); flux_t1 = tmp.flux_t; clear tmp; % load lat x lon RCAE data
    % load(sprintf('%s/%s/ga_malr_diff_si_mon_lat.mat', prefix_proc, par.lat_interp)); ga_malr_diff_orig = ga_malr_diff; clear ga_diff;
    % load(sprintf('%s/%s/ga_dalr_bl_diff_si_mon_lat.mat', prefix_proc, par.lat_interp)); ga_dalr_bl_diff_orig = ga_dalr_bl_diff; clear ga_bl_diff;
    % if strcmp(type, 'gcm') & ~any(strcmp(par.gcm.clim, {'hist-pi'})) & ~contains(par.model, 'mmm')
    %     load(sprintf('%s/sftlf.mat', prefix)); % read land fraction data
    % else

    par2 = par;
    par2.gcm.clim = 'rcp85'; % choose either piControl, historical, or abrupt4xCO2
    par2.gcm.yr_span = '207001-209912'; % number of years that I am considering in the GCM climatology
    prefix2 = make_prefix(type, par2);
    prefix_proc2 = make_prefix_proc(type, par2);

    tmp=load(sprintf('%s/flux_z.mat', prefix_proc2)); flux_z2 = tmp.flux_z; clear tmp; % load lat x mon RCAE data
    % tmp=load(sprintf('%s/flux_t.mat', prefix_proc2)); flux_t2 = tmp.flux_t; clear tmp; % load lat x lon RCAE data

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
            flux_z.(land).r1.(fw) = flux_z2.(land).r1.(fw) - flux_z1.(land).r1.(fw);
            flux_z.(land).ra.(fw) = flux_z2.(land).ra.(fw) - flux_z1.(land).ra.(fw);
            flux_z.(land).res.(fw) = flux_z2.(land).res.(fw) - flux_z1.(land).res.(fw);
            flux_z.(land).stf.(fw) = flux_z2.(land).stf.(fw) - flux_z1.(land).stf.(fw);

            % lat x mon dependence of RCE and RAE
            var_text = '$\Delta (\partial_t h + \nabla \cdot F_m)$';
            figure(); clf; hold all; box on;
            cmp = colCog(30);
            colormap(cmp);
            contourf(mesh_lat, mesh_mon, flux_z.(land).res.(fw), 0.1*[-300 -150:10:150 300], 'linecolor', 'w', 'linewidth', 0.1);
            contour(mesh_lat, mesh_mon, flux_z.(land).res.(fw), [0 0], 'color', 0.75*[1 1 1]);
            make_title_type(type, partitle);
            caxis(0.1*[-150 150]);
            cb = colorbar('limits', 0.1*[-150 150], 'ytick', 0.1*[-140:20:140], 'location', 'eastoutside');
            cb.TickLabelInterpreter = 'latex'; cb.Label.Interpreter = 'latex';
            ylabel(cb, sprintf('%s (Wm$^{-2}$)', var_text));
            ylabel('Latitude (deg)');
            set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_verywide);
            set(gca, 'xlim', [1 12], 'xtick', [1:12], 'xticklabels', par.monlabelnh, 'ylim', [-90 90], 'ytick', [-90:30:90], 'yminortick', 'on', 'tickdir', 'out');
            print(sprintf('%s/dflux/%s/%s/0_div_mon_lat', plotdir, fw, land), '-dpng', '-r300');
            close;

            % lat x mon dependence of STF
            var_text = '$\Delta(\mathrm{LH+SH})$';
            figure(); clf; hold all; box on;
            cmp = colCog(30);
            colormap(cmp);
            flux_z.(land)
            contourf(mesh_lat, mesh_mon, flux_z.(land).stf.(fw), 0.1*[-300:10:300], 'linecolor', 'w', 'linewidth', 0.1);
            contour(mesh_lat, mesh_mon, flux_z.(land).stf.(fw), [0 0], 'color', 0.75*[1 1 1]);
            make_title_type(type, partitle);
            caxis(0.1*[-150 150]);
            cb = colorbar('limits', 0.1*[-150 150], 'ytick', 0.1*[-140:20:140], 'location', 'eastoutside');
            cb.TickLabelInterpreter = 'latex'; cb.Label.Interpreter = 'latex';
            ylabel(cb, sprintf('%s (Wm$^{-2}$)', var_text));
            ylabel('Latitude (deg)');
            set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_verywide);
            set(gca, 'xlim', [1 12], 'xtick', [1:12], 'xticklabels', par.monlabelnh, 'ylim', [-90 90], 'ytick', [-90:30:90], 'yminortick', 'on', 'tickdir', 'out');
            print(sprintf('%s/dflux/%s/%s/0_stf_mon_lat', plotdir, fw, land), '-dpng', '-r300');
            close;

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
                contourf(mesh_lat, mesh_mon, flux_z.(land).ra.(fw), [-30 -15:0.5:15 30], 'linecolor', 'none');
            else
                contourf(mesh_lat, mesh_mon, flux_z.(land).ra.(fw), [-30 -15:1:15 30], 'linecolor', 'w', 'linewidth', 0.1);
            end
            contour(mesh_lat, mesh_mon, flux_z.(land).ra.(fw), [0 0], 'color', 0.75*[1 1 1]);
            make_title_type(type, partitle);
            caxis([-15 15]);
            if strcmp(type, 'echam')
                cb = colorbar('limits', [-50 50], 'ytick', [-50:10:50], 'location', 'eastoutside');
            else
                cb = colorbar('limits', [-15 15], 'ytick', [-15:2.5:15], 'location', 'eastoutside');
            end
            cb.TickLabelInterpreter = 'latex'; cb.Label.Interpreter = 'latex';
            ylabel(cb, sprintf('%s (Wm$^{-2}$)', var_text));
            ylabel('Latitude (deg)');
            set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_verywide);
            set(gca, 'xlim', [1 12], 'xtick', [1:12], 'xticklabels', par.monlabelnh, 'ylim', [-90 90], 'ytick', [-90:30:90], 'yminortick', 'on', 'tickdir', 'out');
            print(sprintf('%s/dflux/%s/%s/0_ra_mon_lat', plotdir, fw, land), '-dpng', '-r300');
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
                    title(sprintf('CMIP5 %s', partitle.gcm.clim));
                else
                    title(sprintf('%s', partitle.model));
                end
                caxis([-0.05 0.05]);
                cb = colorbar('limits', [-0.05 0.05], 'ytick', [-0.05:0.01:0.05], 'location', 'eastoutside');
                cb.TickLabelInterpreter = 'latex'; cb.Label.Interpreter = 'latex';
                ylabel(cb, sprintf('%s (unitless)', var_text));
                ylabel('Latitude (deg)');
                set(gca, 'xlim', [1 12], 'xtick', [1:12], 'xticklabels', par.monlabelnh, 'ylim', [-90 90], 'ytick', [-90:30:90], 'yminortick', 'on', 'tickdir', 'out');
                print(sprintf('%s/dflux/%s/%s/0_r1_mon_lat', plotdir, fw, land), '-dpng', '-r300');
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
                make_title_type(type, partitle);
                caxis([-1 1]);
                cb = colorbar('limits', [-1 1], 'ytick', [-1:0.2:1], 'location', 'eastoutside');
                cb.TickLabelInterpreter = 'latex'; cb.Label.Interpreter = 'latex';
                ylabel(cb, sprintf('%s (unitless)', var_text));
                ylabel('Latitude (deg)');
                set(gca, 'xlim', [1 12], 'xtick', [1:12], 'xticklabels', par.monlabelnh, 'ylim', [-90 90], 'ytick', [-90:30:90], 'yminortick', 'on', 'tickdir', 'out');
                print(sprintf('%s/dflux/%s/%s/0_r1_mon_lat', plotdir, fw, land), '-dpng', '-r300');
                close;

                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % R1z lat x mon dependence of RCE and RAE
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                if ~any(strcmp(fw, {'mse_old', 'dse_old'}))
                    var_text = '$\Delta R_1^*$';
                else
                    var_text = '$\Delta R_1$';
                end
                figure(); clf; hold all; box on;
                cmp = colCog(20);
                colormap(flipud(cmp));
                r1z = flux_z2.(land).res.(fw)./flux_z2.(land).ra.(fw) - flux_z1.(land).res.(fw)./flux_z1.(land).ra.(fw);
                contourf(mesh_lat, mesh_mon, r1z, [-0.5:0.01:0.5], 'linecolor', 'w', 'linewidth', 0.1);
                contour(mesh_lat, mesh_mon, flux_z2.(land).res.(fw)./flux_z2.(land).ra.(fw) - flux_z1.(land).res.(fw)./flux_z1.(land).ra.(fw), [0 0], 'color', 0.75*[1 1 1]);
                [C, h] = contour(mesh_lat, mesh_mon, flux_z2.(land).res.(fw)./flux_z2.(land).ra.(fw) - flux_z1.(land).res.(fw)./flux_z1.(land).ra.(fw), [-0.1:0.05:0.1], 'linecolor', 'w', 'linewidth', 0.1);
                contour(mesh_lat, mesh_mon, flux_z2.(land).res.(fw)./flux_z2.(land).ra.(fw) - flux_z1.(land).res.(fw)./flux_z1.(land).ra.(fw), [-0.1:0.01:0.1], 'color', 'w', 'linewidth', 0.1);
                % [C, h] = contour(mesh_lat, mesh_mon, flux_z.(land).res.(fw)./flux_z.(land).ra.(fw), [par.ep par.ga], 'linecolor', 'w', 'linewidth', 0.25);
                contour(mesh_lat, mesh_mon,  flux_z1.(land).res.(fw)./flux_z1.(land).ra.(fw), par.ep*[1 1], 'linecolor', par.orange, 'linewidth', 2);
                contour(mesh_lat, mesh_mon,  flux_z1.(land).res.(fw)./flux_z1.(land).ra.(fw), par.ga*[1 1], 'linecolor', par.cyan, 'linewidth', 2);
                clabel(C, h, 'fontsize', 6, 'interpreter', 'latex', 'color', 0.75*[1 1 1]);
                % make_title_type(type, partitle);
                title('CMIP5 RCP8.5$-$historical');
                caxis([-0.1 0.1]);
                cb = colorbar('limits', [-0.1 0.1], 'ytick', [-0.1:0.02:0.1], 'location', 'eastoutside');
                cb.TickLabelInterpreter = 'latex'; cb.Label.Interpreter = 'latex';
                ylabel(cb, sprintf('%s (unitless)', var_text));
                ylabel('Latitude (deg)');
                set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_verywide);
                set(gca, 'xlim', [1 12], 'xtick', [1:12], 'xticklabels', par.monlabelnh, 'ylim', [-90 90], 'ytick', [-90:30:90], 'yminortick', 'on', 'tickdir', 'out');
                print(sprintf('%s/dflux/%s/%s/0_r1z_mon_lat', plotdir, fw, land), '-dpng', '-r300');
                data = [mesh_lat(:) mesh_mon(:) r1z(:) ];
                save(sprintf('%s/dflux/%s/%s/0_r1z_mon_lat.dat', plotdir, fw, land), 'data', '-ASCII')
                % if par.make_tikz
                %     matlab2tikz(sprintf('%s/flux/%s/%s/0_r1z_mon_lat.tex', plotdir, fw, land));
                % end
                close;

                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % (TEMPR OVERLAY) R1z lat x mon dependence of RCE and RAE
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                figure(); clf; hold all; box on;
                cmp = colCog(40);
                colormap(cmp);
                contourf(mesh_lat, mesh_mon, dtemp_r, -0.5:0.1:4, 'linecolor', 'none');
                [C,h]=contour(mesh_lat, mesh_mon, dtemp_r, -0.5:0.5:4, '-w', 'linewidth', 0.1);
                % contour(mesh_lat, mesh_mon, dtemp_r, [0 0], 'color', 0.75*[1 1 1]);

                contour(mesh_lat, mesh_mon, flux_z2.(land).res.(fw)./flux_z2.(land).ra.(fw) - flux_z1.(land).res.(fw)./flux_z1.(land).ra.(fw), [0:0.025:0.1], 'linecolor', 'k', 'linewidth', 0.3);
                contour(mesh_lat, mesh_mon, flux_z2.(land).res.(fw)./flux_z2.(land).ra.(fw) - flux_z1.(land).res.(fw)./flux_z1.(land).ra.(fw), [-0.1:0.025:0], 'linestyle', '--', 'linecolor', 'k', 'linewidth', 0.3);
                contour(mesh_lat, mesh_mon,  flux_z1.(land).res.(fw)./flux_z1.(land).ra.(fw), par.ep*[1 1], 'linecolor', par.orange, 'linewidth', 2);
                contour(mesh_lat, mesh_mon,  flux_z1.(land).res.(fw)./flux_z1.(land).ra.(fw), par.ga*[1 1], 'linecolor', par.blue, 'linewidth', 2);

                contour(mesh_lat, mesh_mon,  flux_z2.(land).res.(fw)./flux_z2.(land).ra.(fw), par.ep*[1 1], 'linestyle', '--', 'linecolor', par.orange, 'linewidth', 2);
                contour(mesh_lat, mesh_mon,  flux_z2.(land).res.(fw)./flux_z2.(land).ra.(fw), par.ga*[1 1], 'linestyle', '--', 'linecolor', par.blue, 'linewidth', 2);

                clabel(C, h, 0:0.5:4, 'fontsize', 6, 'interpreter', 'latex', 'color', 0.75*[1 1 1]);
                caxis([-1 3]);
                ylabel('Latitude (deg)');
                make_title_type(type, partitle);
                cb = colorbar('limits', [0 3], 'ytick', [0:0.2:3], 'location', 'eastoutside');
                cb.TickLabelInterpreter = 'latex'; cb.Label.Interpreter = 'latex';
                ylabel(cb, sprintf('$\\Delta T_{0.3}/\\Delta T_{1.0}$ (unitless)'));
                set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_wide);
                set(gca, 'xlim', [1 12], 'xtick', [1:12], 'xticklabels', par.monlabelnh, 'ylim', [-90 90], 'ytick', [-90:30:90], 'yminortick', 'on', 'tickdir', 'out');
                folder = sprintf('%s/dtempr/%s/', plotdir, land);
                if ~exist(folder, 'dir')
                    mkdir(folder)
                end
                print(sprintf('%s/dflux/%s/%s/0_r1z_mon_lat_overlay', plotdir, fw, land), '-dpng', '-r300');
                data = [mesh_lat(:) mesh_mon(:) r1z(:) ];
                save(sprintf('%s/dflux/%s/%s/0_r1z_mon_lat_overlay.dat', plotdir, fw, land), 'data', '-ASCII')
                close;

                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % (TEMPR OVERLAY ALT) R1z lat x mon dependence of RCE and RAE
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                figure(); clf; hold all; box on;
                cmp = colCog(40);
                colormap(cmp);
                contourf(mesh_lat, mesh_mon, dtemp_r, -0.5:0.1:4, 'linecolor', 'none');
                [C,h]=contour(mesh_lat, mesh_mon, dtemp_r, -0.5:0.5:4, '-w', 'linewidth', 0.1);
                % contour(mesh_lat, mesh_mon, dtemp_r, [0 0], 'color', 0.75*[1 1 1]);

                contour(mesh_lat, mesh_mon,  flux_z1.(land).res.(fw)./flux_z1.(land).ra.(fw), par.ep*[1 1], 'linecolor', par.orange, 'linewidth', 2);
                contour(mesh_lat, mesh_mon,  flux_z1.(land).res.(fw)./flux_z1.(land).ra.(fw), par.ga*[1 1], 'linecolor', par.blue, 'linewidth', 2);

                contour(mesh_lat, mesh_mon,  flux_z2.(land).res.(fw)./flux_z2.(land).ra.(fw), par.ep*[1 1], 'linestyle', '--', 'linecolor', par.orange, 'linewidth', 2);
                contour(mesh_lat, mesh_mon,  flux_z2.(land).res.(fw)./flux_z2.(land).ra.(fw), par.ga*[1 1], 'linestyle', '--', 'linecolor', par.blue, 'linewidth', 2);

                clabel(C, h, 0:0.5:4, 'fontsize', 6, 'interpreter', 'latex', 'color', 0.75*[1 1 1]);
                caxis([-1 3]);
                ylabel('Latitude (deg)');
                make_title_type(type, partitle);
                cb = colorbar('limits', [0 3], 'ytick', [0:0.2:3], 'location', 'eastoutside');
                cb.TickLabelInterpreter = 'latex'; cb.Label.Interpreter = 'latex';
                ylabel(cb, sprintf('$\\Delta T_{0.3}/\\Delta T_{1.0}$ (unitless)'));
                set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_wide);
                set(gca, 'xlim', [1 12], 'xtick', [1:12], 'xticklabels', par.monlabelnh, 'ylim', [-90 90], 'ytick', [-90:30:90], 'yminortick', 'on', 'tickdir', 'out');
                folder = sprintf('%s/dtempr/%s/', plotdir, land);
                if ~exist(folder, 'dir')
                    mkdir(folder)
                end
                print(sprintf('%s/dflux/%s/%s/0_r1z_mon_lat_overlay_alt', plotdir, fw, land), '-dpng', '-r300');
                data = [mesh_lat(:) mesh_mon(:) r1z(:) ];
                save(sprintf('%s/dflux/%s/%s/0_r1z_mon_lat_overlay_alt.dat', plotdir, fw, land), 'data', '-ASCII')
                close;

                % if ~any(strcmp(fw, {'mse_old', 'dse_old'}))
                %     var_text = '$\Delta R_1^*$';
                % else
                %     var_text = '$\Delta R_1$';
                % end
                % figure(); clf; hold all; box on;
                % cmp = colCog(20);
                % colormap(flipud(cmp));
                % r1z = flux_z2.(land).res.(fw)./flux_z2.(land).ra.(fw) - flux_z1.(land).res.(fw)./flux_z1.(land).ra.(fw);
                % contourf(mesh_lat, mesh_mon, r1z, [-0.5:0.01:0.5], 'linecolor', 'w', 'linewidth', 0.1);
                % contour(mesh_lat, mesh_mon, flux_z2.(land).res.(fw)./flux_z2.(land).ra.(fw) - flux_z1.(land).res.(fw)./flux_z1.(land).ra.(fw), [0 0], 'color', 0.75*[1 1 1]);
                % [C, h] = contour(mesh_lat, mesh_mon, flux_z2.(land).res.(fw)./flux_z2.(land).ra.(fw) - flux_z1.(land).res.(fw)./flux_z1.(land).ra.(fw), [-0.1:0.05:0.1], 'linecolor', 'w', 'linewidth', 0.1);
                % contour(mesh_lat, mesh_mon, flux_z2.(land).res.(fw)./flux_z2.(land).ra.(fw) - flux_z1.(land).res.(fw)./flux_z1.(land).ra.(fw), [-0.1:0.01:0.1], 'color', 'w', 'linewidth', 0.1);
                % % [C, h] = contour(mesh_lat, mesh_mon, flux_z.(land).res.(fw)./flux_z.(land).ra.(fw), [par.ep par.ga], 'linecolor', 'w', 'linewidth', 0.25);
                % clabel(C, h, 'fontsize', 6, 'interpreter', 'latex', 'color', 0.75*[1 1 1]);
                % % make_title_type(type, partitle);
                % title('CMIP5 RCP8.5$-$historical');
                % caxis([-0.1 0.1]);
                % cb = colorbar('limits', [-0.1 0.1], 'ytick', [-0.1:0.02:0.1], 'location', 'eastoutside');
                % cb.TickLabelInterpreter = 'latex'; cb.Label.Interpreter = 'latex';
                % ylabel(cb, sprintf('%s (unitless)', var_text));
                % ylabel('Latitude (deg)');
                % set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_verywide);
                % set(gca, 'xlim', [1 12], 'xtick', [1:12], 'xticklabels', par.monlabelnh, 'ylim', [-90 90], 'ytick', [-90:30:90], 'yminortick', 'on', 'tickdir', 'out');
                % % if par.make_tikz
                % %     matlab2tikz(sprintf('%s/flux/%s/%s/0_r1z_mon_lat.tex', plotdir, fw, land));
                % % end
                % close;

            end % if hist-pi


        end % for mse dse
    end % for land
end
