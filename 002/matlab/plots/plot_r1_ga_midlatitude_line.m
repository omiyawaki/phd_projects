function plot_r1_ga_midlatitude_line(type, par)
    make_dirs_si_bl(type, par)

    prefix = make_prefix(type, par);
    prefix_proc = make_prefix_proc(type, par);
    plotdir = make_plotdir(type, par);

    load(sprintf('%s/grid.mat', prefix)); % read grid data
    % load(sprintf('%s/flux_z.mat', prefix_proc)); % load lat x mon RCAE_ALT data
    % load(sprintf('%s/ga_frac_mon_lat.mat', prefix_proc));

    lat_list = [10 -10];
    lat_center = 50;
    
    for ep = 1:length(par.ep_swp); par.ep = par.ep_swp(ep);
        for ga = 1:length(par.ga_swp); par.ga = par.ga_swp(ga);

            for f = par.(type).fw; fw = f{1};
                if strcmp(fw, 'mse_old'); r1col = 'k';
                elseif strcmp(fw, 'mse_lat'); r1col = par.blue;
                elseif strcmp(fw, 'mse'); r1col = par.orange;
                end
                for l = par.land_list; land = l{1};
                    if strcmp(land, 'lo'); land_text = 'Land + Ocean';
                    elseif strcmp(land, 'l'); land_text = 'Land';
                    elseif strcmp(land, 'o'); land_text = 'Ocean';
                    end

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

                        lat_print = 'mid';
                        [lat, clat, clat_mon, par] = make_midlatitude_lat(lat_center, par);
                        lat1 = par.lat_center-par.lat_bound;
                        lat2 = par.lat_center+par.lat_bound;

                        load(sprintf('%s/dr1_midlatitude_lat_%g_to_%g.mat', prefix_proc, par.lat_center-par.lat_bound, par.lat_center+par.lat_bound));
                        load(sprintf('%s/si_bl_%g/ga_malr_diff_midlatitude_lat_%g_to_%g_%g.mat', prefix_proc, par.si_bl, par.lat_center-par.lat_bound, par.lat_center+par.lat_bound, par.si_up));

                        figure(); clf; hold all; box on;

                        yyaxis left;
                        ylim_lo = -1;
                        ylim_up = 2.1;
                        rcemax = par.ep;
                        vertices = [1 ylim_lo; 12 ylim_lo; 12 rcemax; 1 rcemax];
                        patch(vertices(:,1), vertices(:,2), par.orange, 'edgecolor', 'none', 'facealpha', 0.5);
                        raemin = par.ga;
                        vertices = [1 raemin; 12 raemin; 12 ylim_up; 1 ylim_up];
                        patch(vertices(:,1), vertices(:,2), par.blue, 'edgecolor', 'none', 'facealpha', 0.5);
                        if strcmp(type, 'rea')
                            plot(1:12, circshift(dr1.r1z_lat.(land).(fw), shiftby), '-', 'color', r1col);
                            mon2 = [1:12, fliplr(1:12)];
                            r1z_spr = nanmean(dr1.r1z_lat.(land).(fw)) + [circshift(dr1_max.dr1z_lat.(land).(fw), par.shiftby, 2), fliplr(circshift(dr1_min.dr1z_lat.(land).(fw), par.shiftby, 2))];
                            fill(mon2, r1z_spr, 0.5*[1 1 1], 'facealpha', 0.2, 'edgealpha', 0.2);
                            if strcmp(fw, 'mse_old')
                                % plot(1:12, circshift(dr1.r1z_lat.(land).mse, shiftby), '--', 'color', r1col);
                                % r1z_spr = nanmean(dr1.r1z_lat.(land).mse) + [circshift(dr1_max.dr1z_lat.(land).mse, par.shiftby, 2), fliplr(circshift(dr1_min.dr1z_lat.(land).mse, par.shiftby, 2))];
                                % fill(mon2, r1z_spr, 0.5*[1 1 1], 'facealpha', 0.2, 'edgealpha', 0.2);
                            end
                        elseif strcmp(type, 'gcm') & contains(par.model, 'mmm')
                            plot(1:12, circshift(dr1.r1z_lat.(land).(fw), shiftby), '-', 'color', r1col);
                            mon2 = [1:12, fliplr(1:12)];
                            r1z_spr = nanmean(dr1.r1z_lat.(land).(fw)) + [circshift(dr1_75.dr1z_lat.(land).(fw), par.shiftby, 2), fliplr(circshift(dr1_25.dr1z_lat.(land).(fw), par.shiftby, 2))];
                            fill(mon2, r1z_spr, 0.5*[1 1 1], 'facealpha', 0.2, 'edgealpha', 0.2);
                            if strcmp(fw, 'mse_old')
                                % plot(1:12, circshift(dr1.r1z_lat.(land).mse, shiftby), '--', 'color', r1col);
                                % r1z_spr = nanmean(dr1.r1z_lat.(land).mse) +  [circshift(dr1_75.dr1z_lat.(land).mse, par.shiftby, 2), fliplr(circshift(dr1_25.dr1z_lat.(land).mse, par.shiftby, 2))];
                                % fill(mon2, r1z_spr, 0.5*[1 1 1], 'facealpha', 0.2, 'edgealpha', 0.2);
                            end
                        else
                            plot(1:12, circshift(dr1.r1z_lat.(land).(fw), shiftby), '-', 'color', r1col);
                            if strcmp(fw, 'mse_old') & ~strcmp(type, 'hahn')
                                % plot(1:12, circshift(dr1.r1z_lat.(land).mse, shiftby), '--', 'color', r1col);
                            end
                        end
                        ylabel('$R_1$ (unitless)');
                        if abs(par.lat_bound) > 45
                            set(gca, 'yminortick', 'on', 'ylim', [0.5 2.1]);
                        else
                            set(gca, 'yminortick', 'on', 'ylim', [-0.4 0.8]);
                            % set(gca, 'yminortick', 'on', 'ylim', [-0.7 0.7]);
                        end

                        yyaxis right;
                        ax = gca;
                        ax.YAxis(1).Color = 'k';
                        if strcmp(type, 'rea')
                            plot(1:12, circshift(ga_frac_lat_mmm.(land), shiftby), '-', 'color', par.orange);
                            mon2 = [1:12, fliplr(1:12)];
                            ga_frac_spr = [circshift(ga_frac_lat_max.(land), par.shiftby, 2), fliplr(circshift(ga_frac_lat_min.(land), par.shiftby, 2))];
                            fill(mon2, ga_frac_spr, par.orange, 'facealpha', 0.2, 'edgealpha', 0.2);
                        elseif strcmp(type, 'gcm') & contains(par.model, 'mmm')
                            plot(1:12, circshift(ga_frac_lat_mmm.(land), shiftby), '-', 'color', par.orange);
                            mon2 = [1:12, fliplr(1:12)];
                            ga_frac_spr = [circshift(ga_frac_lat_75.(land), par.shiftby, 2), fliplr(circshift(ga_frac_lat_25.(land), par.shiftby, 2))];
                            fill(mon2, ga_frac_spr, par.orange, 'facealpha', 0.2, 'edgealpha', 0.2);
                        else
                            plot(1:12, circshift(ga_frac_lat.(land), shiftby), '-', 'color', par.orange);
                        end
                        ylabel(sprintf('$\\left\\langle(\\Gamma_m - \\Gamma)/\\Gamma_m\\right\\rangle_{%g}^{%g}$ (\\%%)', par.si_bl, par.si_up));
                        set(gca, 'fontsize', par.fs, 'xlim', [1 12], 'xtick', [1:12], 'xticklabel', monlabel, 'yminortick', 'on', 'ylim', [-4.86 38]) 
                        % set(gca, 'fontsize', par.fs, 'xlim', [1 12], 'xtick', [1:12], 'xticklabel', monlabel, 'yminortick', 'on', 'ylim', [-10 30]) 
                        ax.YAxis(2).Color = par.orange;

                        make_title_type_lat(type, lat1, lat2, par);
                        set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos)
                        print(sprintf('%s/ga_malr_diff/si_bl_%g/%s/%s/r1_gablft_mon_%g_%s_%s.png', plotdir, par.si_bl, fw, land, par.si_up, hemi, lat_print), '-dpng', '-r300');
                        if par.make_tikz & par.lat_bound > 0 & par.si_bl == 0.7 & par.si_up == 0.3 & strcmp(fw, 'mse_old')
                            matlab2tikz(sprintf('%s/ga_malr_diff/si_bl_%g/%s/%s/r1_gablft_mon_%g_%s_%s.tex', plotdir, par.si_bl, fw, land, par.si_up, hemi, lat_print));
                        end
                        close;

                        if strcmp(fw, 'mse_old') & ~strcmp(type, 'hahn')
                            figure(); clf; box on; hold all;
                            axis off;
                            axis([10,11,10,11])
                            plot(1:12, circshift(dr1.r1z_lat.(land).(fw), shiftby), '-', 'color', r1col);
                            plot(1:12, circshift(dr1.r1z_lat.(land).mse, shiftby), '--', 'color', r1col);
                            legend('$\frac{\partial_t m + \nabla\cdot F_m}{R_a}$', '$\frac{\nabla\cdot F_m}{R_a}$', 'location', 'northwest', 'orientation', 'horizontal');
                            set(gcf, 'paperunits', 'inches', 'paperposition', [0 0 2.7 0.5])
                            print(sprintf('%s/ga_malr_diff/si_bl_%g/%s/%s/legend.png', plotdir, par.si_bl, fw, land), '-dpng', '-r300');
                        end

                        clear r1 ga ga_ft

                    end % lat bounds

                end % land/ocean
            end % framework

        end
    end

end
