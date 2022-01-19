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
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Load LR deviation info
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    tmp=load(sprintf('%s/si_bl_%g/ga_malr_diff_si_mon_lat_%g.mat', prefix_proc, par.si_bl, par.si_up)); lr_ft1 = tmp.ga_malr_diff; clear tmp; % load lat x mon free tropospheric LR deviation data
    tmp=load(sprintf('%s/si_bl_%g/ga_malr_bl_diff_si_mon_lat.mat', prefix_proc, par.si_bl)); lr_bl1 = tmp.ga_malr_bl_diff; clear tmp; % load lat x mon boundary layer LR deviation data
    
    tmp=load(sprintf('%s/si_bl_%g/ga_malr_diff_si_mon_lat_%g.mat', prefix_proc2, par.si_bl, par.si_up)); lr_ft2 = tmp.ga_malr_diff; clear tmp; % load lat x mon free tropospheric LR deviation data
    tmp=load(sprintf('%s/si_bl_%g/ga_malr_bl_diff_si_mon_lat.mat', prefix_proc2, par.si_bl)); lr_bl2 = tmp.ga_malr_bl_diff; clear tmp; % load lat x mon boundary layer LR deviation data

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
            % if alternative R1 definition, use different epsilon
            if strcmp(fw, 'mse')
                disp('test')
                par.ep = 0.3;
            end

            flux_z.(land).r1.(fw) = flux_z2.(land).r1.(fw) - flux_z1.(land).r1.(fw);
            flux_z.(land).ra.(fw) = flux_z2.(land).ra.(fw) - flux_z1.(land).ra.(fw);
            flux_z.(land).res.(fw) = flux_z2.(land).res.(fw) - flux_z1.(land).res.(fw);
            flux_z.(land).stf.(fw) = flux_z2.(land).stf.(fw) - flux_z1.(land).stf.(fw);
            
            ft1 = lr_ft1.(land);
            bl1 = lr_bl1.(land);
            ft2 = lr_ft2.(land);
            bl2 = lr_bl2.(land);

            if strcmp(type, 'gcm') & strcmp(par.gcm.clim, 'hist-pi')
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % R1 lat x mon dependence of RCE and RAE
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
                r1z = flux_z2.(land).res.(fw)./flux_z2.(land).ra.(fw) - flux_z1.(land).res.(fw)./flux_z1.(land).ra.(fw);

                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % (TEMPR OVERLAY -- STIPPLE SMALL R1) R1z lat x mon dependence of RCE and RAE
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                figure(); clf; hold all; box on;
                cmp = colCog(40);
                colormap(cmp);
                contourf(mesh_lat, mesh_mon, dtemp_r, -0.5:0.1:4, 'linecolor', 'none');
                [C,h]=contour(mesh_lat, mesh_mon, dtemp_r, -0.5:0.5:4, '-w', 'linewidth', 0.1);
                contour(mesh_lat, mesh_mon, dtemp_r, [0 0], 'color', 0.75*[1 1 1]);

                clabel(C, h, 0:0.5:4, 'fontsize', 6, 'interpreter', 'latex', 'color', 0.75*[1 1 1]);
                caxis([-1 3]);
                ylabel('Latitude (deg)');
                make_title_type(type, partitle);
                cb = colorbar('limits', [0 3], 'ytick', [0:0.2:3], 'location', 'eastoutside');
                cb.TickLabelInterpreter = 'latex'; cb.Label.Interpreter = 'latex';

                contour(mesh_lat, mesh_mon,  flux_z1.(land).res.(fw)./flux_z1.(land).ra.(fw), par.ep*[1 1], 'linecolor', par.orange, 'linewidth', 2);
                contour(mesh_lat, mesh_mon,  flux_z1.(land).res.(fw)./flux_z1.(land).ra.(fw), par.ga*[1 1], 'linecolor', par.blue, 'linewidth', 2);

                contour(mesh_lat, mesh_mon,  flux_z2.(land).res.(fw)./flux_z2.(land).ra.(fw), par.ep*[1 1], 'linestyle', '--', 'linecolor', par.orange, 'linewidth', 2);
                contour(mesh_lat, mesh_mon,  flux_z2.(land).res.(fw)./flux_z2.(land).ra.(fw), par.ga*[1 1], 'linestyle', '--', 'linecolor', par.blue, 'linewidth', 2);

                ylabel(cb, sprintf('$\\Delta T_{0.3}/\\Delta T_{1.0}$ (unitless)'));
                set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_wide);
                set(gca, 'xlim', [1 12], 'xtick', [1:12], 'xticklabels', par.monlabelnh, 'ylim', [-90 90], 'ytick', [-90:30:90], 'yminortick', 'on', 'tickdir', 'out');
                folder = sprintf('%s/dtempr/%s/', plotdir, land);
                if ~exist(folder, 'dir')
                    mkdir(folder)
                end

                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % CHANGE IN R1
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                tol = 0.03;
                dr1 = flux_z2.(land).res.(fw)./flux_z2.(land).ra.(fw) - flux_z1.(land).res.(fw)./flux_z1.(land).ra.(fw);
                % [c2,h2] = contourf(mesh_lat, mesh_mon, -abs(dr1), -tol*[1, 1]);

                tolp = 10;
                dr1_f = 100*dr1 ./ (flux_z1.(land).res.(fw)./flux_z1.(land).ra.(fw));
                [c2,h2] = contourf(mesh_lat, mesh_mon, -abs(dr1_f), -tolp*[1, 1]);
                set(h2,'linestyle','none','Tag','HatchingRegion');
                ax1 = gca;
                ax2 = copyobj(ax1,figure);
                % Example 1: Default hatching
                hp = findobj(ax1,'Tag','HatchingRegion');
                hh = hatchfill2(hp,'single','LineWidth',1,'Fill','off', 'HatchDensity',100, 'HatchColor', par.gray);

                % filename = sprintf('%s/dflux/%s/%s/0_r1_mon_lat', plotdir, fw, land);
                % export_fig filename -pdf;

                % print(sprintf('%s/dflux/%s/%s/0_r1z_mon_lat_overlay', plotdir, fw, land), '-dpng', '-r300');
                data = [mesh_lat(:) mesh_mon(:) r1z(:) ];
                save(sprintf('%s/dflux/%s/%s/0_r1z_mon_lat_overlay.dat', plotdir, fw, land), 'data', '-ASCII')

                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % (TEMPR OVERLAY -- STIPPLE SMALL DYNAMIC R1 CONTRIB) R1z lat x mon dependence of RCE and RAE
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                figure(); clf; hold all; box on;
                cmp = colCog(40);
                colormap(cmp);
                contourf(mesh_lat, mesh_mon, dtemp_r, -0.5:0.1:4, 'linecolor', 'none');
                [C,h]=contour(mesh_lat, mesh_mon, dtemp_r, -0.5:0.5:4, '-w', 'linewidth', 0.1);
                contour(mesh_lat, mesh_mon, dtemp_r, [0 0], 'color', 0.75*[1 1 1]);

                clabel(C, h, 0:0.5:4, 'fontsize', 6, 'interpreter', 'latex', 'color', 0.75*[1 1 1]);
                caxis([-1 3]);
                ylabel('Latitude (deg)');
                make_title_type(type, partitle);
                cb = colorbar('limits', [0 3], 'ytick', [0:0.2:3], 'location', 'eastoutside');
                cb.TickLabelInterpreter = 'latex'; cb.Label.Interpreter = 'latex';

                contour(mesh_lat, mesh_mon,  flux_z1.(land).res.(fw)./flux_z1.(land).ra.(fw), par.ep*[1 1], 'linecolor', par.orange, 'linewidth', 2);
                contour(mesh_lat, mesh_mon,  flux_z1.(land).res.(fw)./flux_z1.(land).ra.(fw), par.ga*[1 1], 'linecolor', par.blue, 'linewidth', 2);

                contour(mesh_lat, mesh_mon,  flux_z2.(land).res.(fw)./flux_z2.(land).ra.(fw), par.ep*[1 1], 'linestyle', '--', 'linecolor', par.orange, 'linewidth', 2);
                contour(mesh_lat, mesh_mon,  flux_z2.(land).res.(fw)./flux_z2.(land).ra.(fw), par.ga*[1 1], 'linestyle', '--', 'linecolor', par.blue, 'linewidth', 2);

                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % DYNAMIC COMPONENT OF CHANGE IN R1
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                tol = 0.1;
                r1_z1 = flux_z1.(land).res.(fw)./flux_z1.(land).ra.(fw);
                r1_z2 = flux_z2.(land).res.(fw)./flux_z2.(land).ra.(fw);
                dr1 = r1_z2 - r1_z1;
                dra = flux_z2.(land).ra.(fw) - flux_z1.(land).ra.(fw);
                dres = flux_z2.(land).res.(fw) - flux_z1.(land).res.(fw);
                rad = -r1_z1 .* (dra./flux_z1.(land).ra.(fw));
                dyn = r1_z1 .* (dres./flux_z1.(land).res.(fw));
                resi = dr1 - rad - dyn;
                [c2,h2] = contourf(mesh_lat, mesh_mon, dyn, tol*[1, 1]);

                % tolp = 3;
                % dr1_f = 100*dr1 ./ (flux_z1.(land).res.(fw)./flux_z1.(land).ra.(fw));
                % [c2,h2] = contourf(mesh_lat, mesh_mon, -abs(dr1_f), -tolp*[1, 1]);
                set(h2,'linestyle','none','Tag','HatchingRegion');
                ax1 = gca;
                ax2 = copyobj(ax1,figure);
                % Example 1: Default hatching
                hp = findobj(ax1,'Tag','HatchingRegion');
                hh = hatchfill2(hp,'single','LineWidth',1,'Fill','off', 'HatchDensity',100, 'HatchColor', par.gray);

                ylabel(cb, sprintf('$\\Delta T_{0.3}/\\Delta T_{1.0}$ (unitless)'));
                set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_wide);
                set(gca, 'xlim', [1 12], 'xtick', [1:12], 'xticklabels', par.monlabelnh, 'ylim', [-90 90], 'ytick', [-90:30:90], 'yminortick', 'on', 'tickdir', 'out');
                folder = sprintf('%s/dtempr/%s/', plotdir, land);
                if ~exist(folder, 'dir')
                    mkdir(folder)
                end

                % filename = sprintf('%s/dflux/%s/%s/0_r1_mon_lat', plotdir, fw, land);
                % export_fig filename -pdf;

                % print(sprintf('%s/dflux/%s/%s/0_r1z_mon_lat_overlay', plotdir, fw, land), '-dpng', '-r300');
                data = [mesh_lat(:) mesh_mon(:) r1z(:) ];
                save(sprintf('%s/dflux/%s/%s/0_r1z_mon_lat_dynr1_overlay.dat', plotdir, fw, land), 'data', '-ASCII')

                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % plot the dynamic component of r1 change
                % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                [c2,h2] = contourf(mesh_lat, mesh_mon, dyn, tol*[1, 1]);
                var_text = '$\frac{\Delta (\partial_t m+ \partial_y(vm))}{(\partial_t m+ \partial_y(vm))_{\mathrm{control}}}$';
                figure(); clf; hold all; box on;
                cmp = colCog(30);
                colormap(flipud(cmp));
                contourf(mesh_lat, mesh_mon, dyn, [-0.2:0.01:0.2], 'linecolor', 'w', 'linewidth', 0.25);
                contour(mesh_lat, mesh_mon, dyn, [0 0], 'linecolor', par.gray, 'linewidth', 0.25);
                clabel(C, h, 'fontsize', 6, 'interpreter', 'latex', 'color', 'k');
                make_title_type(type, partitle);
                caxis([-0.15 0.15]);
                cb = colorbar('limits', [-0.1 0.1], 'ytick', [-0.1:0.02:0.1], 'location', 'eastoutside');
                cb.TickLabelInterpreter = 'latex'; cb.Label.Interpreter = 'latex';
                ylabel(cb, sprintf('%s (unitless)', var_text));
                ylabel('Latitude (deg)');
                set(gca, 'xlim', [1 12], 'xtick', [1:12], 'xticklabels', par.monlabelnh, 'ylim', [-90 90], 'ytick', [-90:30:90], 'yminortick', 'on', 'tickdir', 'out');
                print(sprintf('%s/dflux/%s/%s/0_dyn_mon_lat', plotdir, fw, land), '-dpng', '-r300');
                close;


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
                % (TEMPR OVERLAY ALT) LR DEV lat x mon dependence of MA and SI
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                figure(); clf; hold all; box on;
                cmp = colCog(40);
                colormap(cmp);
                contourf(mesh_lat, mesh_mon, dtemp_r, -0.5:0.1:4, 'linecolor', 'none');
                [C,h]=contour(mesh_lat, mesh_mon, dtemp_r, -0.5:0.5:4, '-w', 'linewidth', 0.1);
                % contour(mesh_lat, mesh_mon, dtemp_r, [0 0], 'color', 0.75*[1 1 1]);

                contour(mesh_lat, mesh_mon,  ft1, 13*[1 1], 'linecolor', par.orange, 'linewidth', 2);
                contour(mesh_lat, mesh_mon,  bl1, 100*[1 1], 'linecolor', par.blue, 'linewidth', 2);

                contour(mesh_lat, mesh_mon,  ft2, 13*[1 1], 'linestyle', '--', 'linecolor', par.orange, 'linewidth', 2);
                contour(mesh_lat, mesh_mon,  bl2, 100*[1 1], 'linestyle', '--', 'linecolor', par.blue, 'linewidth', 2);

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
                print(sprintf('%s/dflux/%s/%s/0_ga_malr_diff_mon_lat_overlay_alt', plotdir, fw, land), '-dpng', '-r300');
                data = [mesh_lat(:) mesh_mon(:) r1z(:) ];
                save(sprintf('%s/dflux/%s/%s/0_ga_malr_diff_mon_lat_overlay_alt.dat', plotdir, fw, land), 'data', '-ASCII')
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
               
                r1z1 = flux_z1.(land).res.(fw)./flux_z1.(land).ra.(fw);
                r1z2 = flux_z2.(land).res.(fw)./flux_z2.(land).ra.(fw);

                contour(mesh_lat, mesh_mon,  r1z1, par.ep*[1 1], 'linecolor', par.orange, 'linewidth', 2);
                contour(mesh_lat, mesh_mon,  r1z1, par.ga*[1 1], 'linecolor', par.blue, 'linewidth', 2);

                contour(mesh_lat, mesh_mon,  r1z2, par.ep*[1 1], 'linestyle', '--', 'linecolor', par.orange, 'linewidth', 2);
                contour(mesh_lat, mesh_mon,  r1z2, par.ga*[1 1], 'linestyle', '--', 'linecolor', par.blue, 'linewidth', 2);

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
                data = [mesh_lat(:) mesh_mon(:) dtemp_r(:)];
                save(sprintf('%s/dflux/%s/%s/0_dtempr_mon_lat_final.dat', plotdir, fw, land), 'data', '-ASCII')
                clear data;
                data = [mesh_lat(:) mesh_mon(:) r1z1(:)];
                save(sprintf('%s/dflux/%s/%s/0_r1z1_mon_lat_final.dat', plotdir, fw, land), 'data', '-ASCII')
                clear data;
                data = [mesh_lat(:) mesh_mon(:) r1z2(:)];
                save(sprintf('%s/dflux/%s/%s/0_r1z2_mon_lat_final.dat', plotdir, fw, land), 'data', '-ASCII')
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
