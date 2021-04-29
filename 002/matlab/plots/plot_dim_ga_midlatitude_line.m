function plot_dim_ga_midlatitude_line(type, par)
    make_dirs_si_bl(type, par)

    prefix = make_prefix(type, par);
    prefix_proc = make_prefix_proc(type, par);
    plotdir = make_plotdir(type, par);

    load(sprintf('%s/grid.mat', prefix)); % read grid data
    load(sprintf('%s/flux_z.mat', prefix_proc)); % load lat x mon RCAE_ALT data
    load(sprintf('%s/ga_frac_mon_lat.mat', prefix_proc));

    %% SPATIAL COMPARISON
    mon_list = [1 7];

    %% TEMPORAL COMPARISON
    lat_list = [10 -10 80 -80];
    lat_center = 50;
    
    for ep = 1:length(par.ep_swp); par.ep = par.ep_swp(ep);
        for ga = 1:length(par.ga_swp); par.ga = par.ga_swp(ga);

            for f = par.(type).fw; fw = f{1};
                if strcmp(fw, 'mse_old'); r1col = 'k';
                elseif strcmp(fw, 'mse_lat'); r1col = par.blue;
                elseif strcmp(fw, 'mse'); r1col = par.orange;
                end
                for l = par.land_list; land = l{1};
                % for l = {'lo'}; land = l{1};
                    if strcmp(land, 'lo'); land_text = 'Land + Ocean';
                    elseif strcmp(land, 'l'); land_text = 'Land';
                    elseif strcmp(land, 'o'); land_text = 'Ocean';
                    end

                    % for mn = 1:length(mon_list); mon_eval = mon_list(mn);
                    %     mon_str = make_mon_str(mon_eval);

                    %     % select month
                    %     r1 = flux_z.(land).res.(fw)(:,mon_eval)./flux_z.(land).ra.(fw)(:,mon_eval);
                    %     ga = 100 * ga_frac.(land)(:,mon_eval,:); % multiply by 100 to express in percentage
                    %     ga = permute(ga, [3 1 2]);
                    %     ga_v = interp1(grid.dim3.si, ga, linspace(1, par.si_up, 101));
                    %     ga_v = squeeze(nanmean(ga_v, 1))';
                    %     ga_ft = interp1(grid.dim3.si, ga, linspace(par.si_bl, par.si_up, 101));
                    %     ga_ft = squeeze(nanmean(ga_ft, 1))';
                    %     ga_bl = interp1(grid.dim3.si, ga, linspace(1, par.si_bl, 101));
                    %     % ga_bl = interp1(grid.dim3.si, ga, linspace(1, 0.9, 101));
                    %     ga_bl = squeeze(nanmean(ga_bl, 1))';

                    %     figure(); clf; hold all; box on;
                    %     yyaxis left;
                    %     ylim_lo = -1;
                    %     ylim_up = 1.6;
                    %     rcemax = par.ep;
                    %     vertices = [-90 ylim_lo; 90 ylim_lo; 90 rcemax; -90 rcemax];
                    %     patch(vertices(:,1), vertices(:,2), par.orange, 'edgecolor', 'none', 'facealpha', 0.5);
                    %     raemin = par.ga;
                    %     vertices = [-90 raemin; 90 raemin; 90 ylim_up; -90 ylim_up];
                    %     patch(vertices(:,1), vertices(:,2), par.blue, 'edgecolor', 'none', 'facealpha', 0.5);
                    %     plot(grid.dim2.lat, r1, '-', 'color', r1col);
                    %     ylabel('$R_1$ (unitless)');
                    %     set(gca, 'yminortick', 'on', 'ylim', [-1 1.6]);
                    %     yyaxis right;
                    %     plot(grid.dim3.lat, ga_v, '-', 'color', par.green);
                    %     ylabel(sprintf('$\\left\\langle(\\Gamma_m - \\Gamma)/\\Gamma_m\\right\\rangle_{%g}^{%g}$ (\\%%)', 1, par.si_up));
                    %     xlabel('Latitude (deg)');
                    %     make_title_type_mon(type, mon_str, par);
                    %     set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos)
                    %     set(gca, 'fontsize', par.fs, 'xlim', [-90 90], 'xtick', [-90:30:90], 'xminortick', 'on', 'yminortick', 'on', 'ylim', [-40 80])
                    %     ax = gca;
                    %     ax.YAxis(1).Color = 'k';
                    %     ax.YAxis(2).Color = par.green;
                    %     print(sprintf('%s/ga_malr_diff/si_bl_%g/%s/%s/r1_ga_lat_%g_%02d.png', plotdir, par.si_bl, fw, land, par.si_up, mon_eval), '-dpng', '-r300');
                    %     close;

                    %     figure(); clf; hold all; box on; grid on;
                    %     yyaxis left;
                    %     ylim_lo = -0.9;
                    %     ylim_up = 1.6;
                    %     rcemax = par.ep;
                    %     vertices = [-90 ylim_lo; 90 ylim_lo; 90 rcemax; -90 rcemax];
                    %     patch(vertices(:,1), vertices(:,2), par.orange, 'edgecolor', 'none', 'facealpha', 0.5);
                    %     raemin = par.ga;
                    %     vertices = [-90 raemin; 90 raemin; 90 ylim_up; -90 ylim_up];
                    %     patch(vertices(:,1), vertices(:,2), par.blue, 'edgecolor', 'none', 'facealpha', 0.5);
                    %     plot(grid.dim2.lat, r1, '-', 'color', r1col);
                    %     ylabel('$R_1$ (unitless)');
                    %     set(gca, 'yminortick', 'on', 'ylim', [-0.9 0.7]);
                    %     yyaxis right;
                    %     plot(grid.dim3.lat, ga_ft, '-', 'color', par.maroon);
                    %     ylabel(sprintf('$\\left\\langle(\\Gamma_m - \\Gamma)/\\Gamma_m\\right\\rangle_{%g}^{%g}$ (\\%%)', par.si_bl, par.si_up));
                    %     xlabel('Latitude (deg)');
                    %     make_title_type_mon(type, mon_str, par);
                    %     set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos)
                    %     set(gca, 'fontsize', par.fs, 'xlim', [-75 75], 'xtick', [-75:15:75], 'xminortick', 'on', 'yminortick', 'on', 'ylim', [-13 31])
                    %     ax = gca;
                    %     ax.YAxis(1).Color = 'k';
                    %     ax.YAxis(2).Color = par.maroon;
                    %     print(sprintf('%s/ga_malr_diff/si_bl_%g/%s/%s/r1_gaft_lat_%g_%02d.png', plotdir, par.si_bl, fw, land, par.si_up, mon_eval), '-dpng', '-r300');
                    %     close;

                    %     figure(); clf; hold all; box on; grid on;
                    %     yyaxis left;
                    %     ylim_lo = -1;
                    %     ylim_up = 1.6;
                    %     rcemax = par.ep;
                    %     vertices = [-90 ylim_lo; 90 ylim_lo; 90 rcemax; -90 rcemax];
                    %     patch(vertices(:,1), vertices(:,2), par.orange, 'edgecolor', 'none', 'facealpha', 0.5);
                    %     raemin = par.ga;
                    %     vertices = [-90 raemin; 90 raemin; 90 ylim_up; -90 ylim_up];
                    %     patch(vertices(:,1), vertices(:,2), par.blue, 'edgecolor', 'none', 'facealpha', 0.5);
                    %     plot(grid.dim2.lat, r1, '-', 'color', r1col);
                    %     ylabel('$R_1$ (unitless)');
                    %     set(gca, 'yminortick', 'on', 'ylim', [0.55 1.6]);
                    %     yyaxis right;
                    %     plot(grid.dim3.lat, ga_bl, '-', 'color', par.darkblue);
                    %     ylabel(sprintf('$\\left\\langle(\\Gamma_m - \\Gamma)/\\Gamma_m\\right\\rangle_{%g}^{%g}$ (\\%%)', 1, par.si_bl));
                    %     xlabel('Latitude (deg)');
                    %     make_title_type_mon(type, mon_str, par);
                    %     set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_short)
                    %     set(gca, 'fontsize', par.fs, 'xlim', [-90 90], 'xtick', [-90:30:90], 'xminortick', 'on', 'yminortick', 'on', 'ylim', [0 300])
                    %     ax = gca;
                    %     ax.YAxis(1).Color = 'k';
                    %     ax.YAxis(2).Color = par.darkblue;
                    %     print(sprintf('%s/ga_malr_diff/si_bl_%g/%s/%s/r1_gabl_lat_%g_%02d.png', plotdir, par.si_bl, fw, land, par.si_up, mon_eval), '-dpng', '-r300');
                    %     close;

                    %     clear r1 ga ga_ft

                    % end

                    for lb = 1:length(lat_list); par.lat_bound = lat_list(lb);

                        if par.lat_bound > 0
                            hemi = 'nh';
                            shiftby = 0;
                            monlabel = par.monlabelnh;
                        else
                            hemi = 'sh';
                            shiftby = 6;
                            monlabel = par.monlabelsh;
                        end

                        if abs(par.lat_bound) > 45  % polar
                            lat_print = 'polar';
                            [lat, clat, clat_mon, par] = make_polar_lat(par);
                            lat1 = par.lat_bound;
                            lat2 = par.lat_pole;
                        else % midlats
                            lat_print = 'mid';
                            [lat, clat, clat_mon, par] = make_midlatitude_lat(lat_center, par);
                            lat1 = par.lat_center-par.lat_bound;
                            lat2 = par.lat_center+par.lat_bound;
                        end

                        % interpolate at select latitude
                        rce_crit = interp1(grid.dim2.lat, flux_z.(land).res.(fw)./flux_z.(land).ra.(fw), lat);
                        [lh, sh] = rename_stf(type, flux_z, land);
                        % rae_crit = interp1(grid.dim2.lat, flux_z.(land).lw+lh+sh, lat);
                        rae_crit = interp1(grid.dim2.lat, flux_z.(land).res.(fw)./flux_z.(land).ra.(fw), lat);
                        ga = 100 * interp1(grid.dim3.lat, ga_frac.(land), lat); % multiply by 100 to express in percentage
                        rce_crit = nansum(clat_mon .* rce_crit)/nansum(clat); % area averaged weighting
                        rae_crit = nansum(clat_mon .* rae_crit)/nansum(clat); % area averaged weighting
                        ga = nansum(clat_mon .* ga)/nansum(clat); % area averaged weighting

                        ga = permute(ga, [3 1 2]);
                        ga_v = interp1(grid.dim3.si, ga, linspace(1, par.si_up, 101));
                        ga_v = squeeze(nanmean(ga_v, 1))';
                        ga_ft = interp1(grid.dim3.si, ga, linspace(par.si_bl, par.si_up, 101));
                        ga_ft = squeeze(nanmean(ga_ft, 1))';
                        ga_bl = interp1(grid.dim3.si, ga, linspace(1, par.si_bl, 101));
                        ga_bl = squeeze(nanmean(ga_bl, 1))';

                        figure(); clf; hold all; box on;

                        yyaxis left;
                        % ylim_lo = -1;
                        % ylim_up = 1.6;
                        % rcemax = par.ep;
                        % vertices = [1 ylim_lo; 12 ylim_lo; 12 rcemax; 1 rcemax];
                        % patch(vertices(:,1), vertices(:,2), par.orange, 'edgecolor', 'none', 'facealpha', 0.5);
                        % raemin = par.ga;
                        % vertices = [1 raemin; 12 raemin; 12 ylim_up; 1 ylim_up];
                        % patch(vertices(:,1), vertices(:,2), par.blue, 'edgecolor', 'none', 'facealpha', 0.5);
                        if abs(par.lat_bound) > 45
                            plot(1:12, circshift(rae_crit, shiftby), '-', 'color', r1col);
                            ylabel('LH$+$SH (W m$^{-2}$)');
                            set(gca, 'yminortick', 'on', 'ylim', [-inf inf]);
                        else
                            plot(1:12, circshift(rce_crit, shiftby), '-', 'color', r1col);
                            ylabel('$\partial_t h + \nabla\cdot F_m$ (W m$^{-2}$)');
                            set(gca, 'yminortick', 'on', 'ylim', [-inf inf]);
                        end

                        yyaxis right;
                            plot(1:12, circshift(ga_v, shiftby), '-', 'color', par.green);
                        ylabel(sprintf('$\\left\\langle(\\Gamma_m - \\Gamma)/\\Gamma_m\\right\\rangle_{%g}^{%g}$ (\\%%)', 1, par.si_up));

                        make_title_type_lat(type, lat1, lat2, par);
                        set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos)
                        set(gca, 'fontsize', par.fs, 'xlim', [1 12], 'xtick', [1:12], 'xticklabel', monlabel, 'yminortick', 'on', 'ylim', [-inf inf]) 
                        ax = gca;
                        ax.YAxis(1).Color = 'k';
                        ax.YAxis(2).Color = par.green;
                        print(sprintf('%s/ga_malr_diff/si_bl_%g/%s/%s/dim_ga_mon_%g_%s_%s.png', plotdir, par.si_bl, fw, land, par.si_up, hemi, lat_print), '-dpng', '-r300');
                        close;

                        figure(); clf; hold all; box on;

                        yyaxis left;
                        % ylim_lo = -1;
                        % ylim_up = 1.6;
                        % rcemax = par.ep;
                        % vertices = [1 ylim_lo; 12 ylim_lo; 12 rcemax; 1 rcemax];
                        % patch(vertices(:,1), vertices(:,2), par.orange, 'edgecolor', 'none', 'facealpha', 0.5);
                        % raemin = par.ga;
                        % vertices = [1 raemin; 12 raemin; 12 ylim_up; 1 ylim_up];
                        % patch(vertices(:,1), vertices(:,2), par.blue, 'edgecolor', 'none', 'facealpha', 0.5);
                        if abs(par.lat_bound) > 45
                            plot(1:12, circshift(rae_crit, shiftby), '-', 'color', r1col);
                            ylabel('LH$+$SH (W m$^{-2}$)');
                            set(gca, 'yminortick', 'on', 'ylim', [-inf inf]);
                        else
                            plot(1:12, circshift(rce_crit, shiftby), '-', 'color', r1col);
                            ylabel('$\partial_t h + \nabla\cdot F_m$ (W m$^{-2}$)');
                            set(gca, 'yminortick', 'on', 'ylim', [-inf inf]);
                        end

                        yyaxis right;
                        ax = gca;
                        ax.YAxis(1).Color = 'k';
                        if abs(par.lat_bound) > 45
                            plot(1:12, circshift(ga_bl, shiftby), '-', 'color', par.darkblue);
                            ylabel(sprintf('$\\left\\langle(\\Gamma_m - \\Gamma)/\\Gamma_m\\right\\rangle_{%g}^{%g}$ (\\%%)', 1, par.si_bl));
                            if strcmp(type, 'echam')
                                set(gca, 'fontsize', par.fs, 'xlim', [1 12], 'xtick', [1:12], 'xticklabel', monlabel, 'yminortick', 'on', 'ylim', [-50 450]) 
                            else
                                set(gca, 'fontsize', par.fs, 'xlim', [1 12], 'xtick', [1:12], 'xticklabel', monlabel, 'yminortick', 'on', 'ylim', [0 300]) 
                            end
                            ax.YAxis(2).Color = par.darkblue;
                        else
                            plot(1:12, circshift(ga_ft, shiftby), '-', 'color', par.maroon);
                            ylabel(sprintf('$\\left\\langle(\\Gamma_m - \\Gamma)/\\Gamma_m\\right\\rangle_{%g}^{%g}$ (\\%%)', par.si_bl, par.si_up));
                            set(gca, 'fontsize', par.fs, 'xlim', [1 12], 'xtick', [1:12], 'xticklabel', monlabel, 'yminortick', 'on', 'ylim', [-13 31]) 
                            ax.YAxis(2).Color = par.maroon;
                        end

                        make_title_type_lat(type, lat1, lat2, par);
                        set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos)
                        print(sprintf('%s/ga_malr_diff/si_bl_%g/%s/%s/dim_gablft_mon_%g_%s_%s.png', plotdir, par.si_bl, fw, land, par.si_up, hemi, lat_print), '-dpng', '-r300');
                        close;

                        clear r1 ga ga_ft

                    end % lat bounds

                end % land/ocean
            end % framework

        end
    end

end
