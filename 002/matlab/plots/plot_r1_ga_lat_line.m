function plot_r1_ga_midlatitude_line(type, par)
    make_dirs_si_bl(type, par)

    prefix = make_prefix(type, par);
    prefix_proc = make_prefix_proc(type, par);
    plotdir = make_plotdir(type, par);

    load(sprintf('%s/grid.mat', prefix)); % read grid data
    load(sprintf('%s/flux_zt.mat', prefix_proc)); % load lat x mon RCAE_ALT data
    load(sprintf('%s/si_bl_%g/ga_malr_diff_si_lat_%g.mat', prefix_proc, par.si_bl, par.si_up));
    load(sprintf('%s/si_bl_%g/ga_malr_bl_diff_si_lat.mat', prefix_proc, par.si_bl));
    load(sprintf('%s/si_bl_%g/ga_dalr_bl_diff_si_lat.mat', prefix_proc, par.si_bl));

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

                    for t = {'ann'}; time = t{1};
                    % for t = {'ann', 'djf', 'jja', 'mam', 'son'}; time = t{1};

                        figure(); clf; hold all; box on;
                        if strcmp(type, 'rea')
                            r1z_min = flux_zt_min.(land).(time).r1z.(fw);
                            r1z_max = flux_zt_max.(land).(time).r1z.(fw);
                            r1z_l = r1z_min;
                            r1z_u = r1z_max;

                            ga_ft_min = ga_malr_diff_zt_min.(land).(time);
                            ga_ft_max = ga_malr_diff_zt_max.(land).(time);
                            ga_ft_l = ga_ft_min;
                            ga_ft_u = ga_ft_max;
                        elseif (any(strcmp(type, {'gcm'})) & strcmp(par.model, 'mmm'))
                            r1z_25 = flux_zt_25.(land).(time).r1z.(fw);
                            r1z_75 = flux_zt_75.(land).(time).r1z.(fw);
                            r1z_l = r1z_25;
                            r1z_u = r1z_75;

                            ga_ft_25 = ga_malr_diff_zt_25.(land).(time);
                            ga_ft_75 = ga_malr_diff_zt_75.(land).(time);
                            ga_ft_l = ga_ft_25;
                            ga_ft_u = ga_ft_75;
                        end
                        yyaxis left;
                        ylim_lo = -0.6;
                        ylim_up = 1.8;
                        rcemax = par.ep;
                        vertices = [-90 ylim_lo; 90 ylim_lo; 90 rcemax; -90 rcemax];
                        patch(vertices(:,1), vertices(:,2), par.orange, 'edgecolor', 'none', 'facealpha', 0.5);
                        raemin = par.ga;
                        vertices = [-90 raemin; 90 raemin; 90 ylim_up; -90 ylim_up];
                        patch(vertices(:,1), vertices(:,2), par.blue, 'edgecolor', 'none', 'facealpha', 0.5);
                        if strcmp(type, 'rea') | (any(strcmp(type, {'gcm'})) & strcmp(par.model, 'mmm'))
                            lat2 = [grid.dim3.lat', fliplr(grid.dim3.lat')];
                            inbtw = [r1z_l, fliplr(r1z_u)];
                            fill(lat2, inbtw, 'k', 'facealpha', 0.2, 'edgealpha', 0.2);
                        end
                        plot(grid.dim2.lat, flux_zt.(land).(time).r1z.(fw), '-', 'color', r1col);
                        ylabel('$R_1$ (unitless)');
                        set(gca, 'yminortick', 'on', 'ylim', [ylim_lo ylim_up]);
                        yyaxis right;
                        if strcmp(type, 'rea') | (any(strcmp(type, {'gcm'})) & strcmp(par.model, 'mmm'))
                            lat2 = [grid.dim3.lat', fliplr(grid.dim3.lat')];
                            inbtw = [ga_ft_l, fliplr(ga_ft_u)];
                            fill(lat2, inbtw, par.orange, 'facealpha', 0.2, 'edgealpha', 0.2);
                        end
                        disp(sprintf('Minimum deviation for fw=%s, time=%s, lower=%g, and upper=%g is %g%%.', fw, upper(time), par.si_bl, par.si_up, nanmin(nanmin(ga_malr_diff_zt.(land).(time)))))
                        plot(grid.dim3.lat, ga_malr_diff_zt.(land).(time), '-', 'color', par.orange);
                        ylabel(sprintf('$\\left\\langle(\\Gamma_m - \\Gamma)/\\Gamma_m\\right\\rangle_{%g}^{%g}$ (\\%%)', par.si_bl, par.si_up));
                        xlabel('Latitude (deg)');
                        make_title_type_time(type, time, par);
                        set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_wide)
                        set(gca, 'fontsize', par.fs, 'xlim', [-90 90], 'xtick', [-90:30:90], 'xminortick', 'on', 'yminortick', 'on', 'ylim', [-7 61.57])
                        ax = gca;
                        ax.YAxis(1).Color = 'k';
                        ax.YAxis(2).Color = par.orange;
                        folder = sprintf('%s/ga_malr_diff/si_bl_%g/%s/%s/%s', plotdir, par.si_bl, fw, land, time);
                        if ~exist(folder, 'dir')
                            mkdir(folder);
                        end
                        print(sprintf('%s/r1_gaft_lat_%g.png', folder, par.si_up), '-dpng', '-r300');
                        if par.make_tikz & par.si_bl == 0.7 & strcmp(fw, 'mse_old')
                            matlab2tikz(sprintf('%s/r1_gaft_lat_%g.tex', folder, par.si_up));
                        end
                        close;

                        figure(); clf; hold all; box on;
                        if strcmp(type, 'rea')
                            r1z_min = flux_zt_min.(land).(time).r1z.(fw);
                            r1z_max = flux_zt_max.(land).(time).r1z.(fw);
                            r1z_l = r1z_min;
                            r1z_u = r1z_max;

                            ga_bl_min = ga_malr_bl_diff_zt_min.(land).(time);
                            ga_bl_max = ga_malr_bl_diff_zt_max.(land).(time);
                            ga_bl_l = ga_bl_min;
                            ga_bl_u = ga_bl_max;
                        elseif (any(strcmp(type, {'gcm'})) & strcmp(par.model, 'mmm'))
                            r1z_25 = flux_zt_25.(land).(time).r1z.(fw);
                            r1z_75 = flux_zt_75.(land).(time).r1z.(fw);
                            r1z_l = r1z_25;
                            r1z_u = r1z_75;

                            ga_bl_25 = ga_malr_bl_diff_zt_25.(land).(time);
                            ga_bl_75 = ga_malr_bl_diff_zt_75.(land).(time);
                            ga_bl_l = ga_bl_25;
                            ga_bl_u = ga_bl_75;
                        end
                        yyaxis left;
                        ylim_lo = -0.6;
                        ylim_up = 1.8;
                        rcemax = par.ep;
                        vertices = [-90 ylim_lo; 90 ylim_lo; 90 rcemax; -90 rcemax];
                        patch(vertices(:,1), vertices(:,2), par.orange, 'edgecolor', 'none', 'facealpha', 0.5);
                        raemin = par.ga;
                        vertices = [-90 raemin; 90 raemin; 90 ylim_up; -90 ylim_up];
                        patch(vertices(:,1), vertices(:,2), par.blue, 'edgecolor', 'none', 'facealpha', 0.5);
                        if strcmp(type, 'rea') | (any(strcmp(type, {'gcm'})) & strcmp(par.model, 'mmm'))
                            lat2 = [grid.dim3.lat', fliplr(grid.dim3.lat')];
                            inbtw = [r1z_l, fliplr(r1z_u)];
                            fill(lat2, inbtw, 'k', 'facealpha', 0.2, 'edgealpha', 0.2);
                        end
                        plot(grid.dim2.lat, flux_zt.(land).(time).r1z.(fw), '-', 'color', r1col);
                        ylabel('$R_1$ (unitless)');
                        set(gca, 'yminortick', 'on', 'ylim', [ylim_lo ylim_up]);
                        yyaxis right;
                        if strcmp(type, 'rea') | (any(strcmp(type, {'gcm'})) & strcmp(par.model, 'mmm'))
                            lat2 = [grid.dim3.lat', fliplr(grid.dim3.lat')];
                            inbtw = [ga_bl_l, fliplr(ga_bl_u)];
                            fill(lat2, inbtw, par.blue, 'facealpha', 0.2, 'edgealpha', 0.2);
                        end
                        plot(grid.dim3.lat, ga_malr_bl_diff_zt.(land).(time), '-', 'color', par.blue);
                        ylabel(sprintf('$\\left\\langle(\\Gamma_m - \\Gamma)/\\Gamma_m\\right\\rangle_{%g}^{%g}$ (\\%%)', 1, par.si_bl));
                        xlabel('Latitude (deg)');
                        make_title_type_time(type, time, par);
                        set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_wide)
                        set(gca, 'fontsize', par.fs, 'xlim', [-90 90], 'xtick', [-90:30:90], 'xminortick', 'on', 'yminortick', 'on', 'ylim', [-316.67 350])
                        ax = gca;
                        ax.YAxis(1).Color = 'k';
                        ax.YAxis(2).Color = par.blue;
                        folder = sprintf('%s/ga_malr_diff/si_bl_%g/%s/%s/%s', plotdir, par.si_bl, fw, land, time);
                        if ~exist(folder, 'dir')
                            mkdir(folder);
                        end
                        print(sprintf('%s/r1_gabl_lat.png', folder), '-dpng', '-r300');
                        if par.make_tikz & par.si_bl == 0.9 & strcmp(fw, 'mse_old')
                            matlab2tikz(sprintf('%s/r1_gabl_lat.tex', folder));
                        end
                        close;

                    end

                end % land/ocean
            end % framework

        end
    end

end
