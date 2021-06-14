function plot_ga_frac_ann(type, par)

    prefix = make_prefix(type, par);
    prefix_proc = make_prefix_proc(type, par);
    plotdir = make_plotdir(type, par);

    load(sprintf('%s/grid.mat', prefix)); % read grid data
    load(sprintf('%s/flux_z.mat', prefix_proc)); % load lat x mon RCAE_ALT data
    load(sprintf('%s/ga_frac_mon_lat.mat', prefix_proc));

    % if strcmp(type, 'gcm') & contains(par.model, 'mmm')
    %     grid.dim3.lat = grid.dim3.lat';
    % end

    make_dirs_ep(type, par)

    time = 'ann';
    crit = 'def';

    [mesh_si, mesh_lat] = meshgrid(grid.dim3.lat, grid.dim3.si);

    for f = {'mse_old'}; fw = f{1};
        for l = par.land_list; land = l{1};
            if strcmp(land, 'lo'); land_text = 'Land + Ocean';
            elseif strcmp(land, 'l'); land_text = 'Land';
            elseif strcmp(land, 'o'); land_text = 'Ocean';
            end

            % take annual mean
            r1_ann = nanmean(flux_z.(land).res.(fw)./flux_z.(land).ra.(fw), 2);
            ga_fr_ann = 100 * squeeze(nanmean(ga_frac.(land), 2)); % multiply by 100 to express as percentage
            if strcmp(type, 'rea') | ( strcmp(type, 'gcm') & strcmp(par.model, 'mmm') )
                ga_fr_std_ann = 100 * squeeze(nanmean(ga_frac_std.(land), 2));
                ga_fr_min_ann = 100 * squeeze(nanmean(ga_frac_min.(land), 2));
                ga_fr_max_ann = 100 * squeeze(nanmean(ga_frac_max.(land), 2));
                ga_fr_25_ann = 100 * squeeze(nanmean(ga_frac_25.(land), 2));
                ga_fr_75_ann = 100 * squeeze(nanmean(ga_frac_75.(land), 2));
            end

            % lat x lev contour of lapse rate deviation
            figure(); clf; hold all; box on;
            cmp = colCog(40);
            colormap(cmp);
            contourf(mesh_si, mesh_lat, ga_fr_ann', [-300:5:100 150:50:1000], 'linecolor', 'w', 'linewidth', 0.1);
            contour(mesh_si, mesh_lat, ga_fr_ann', [0 0], 'linecolor', 0.75*[1 1 1], 'linewidth', 0.1);
            [C,h] = contour(mesh_si, mesh_lat, ga_fr_ann', [-300:10:100 200:100:1000], 'linecolor', 'w', 'linewidth', 0.1);
            clabel(C, h, 'fontsize', 6, 'interpreter', 'latex', 'color', 0.75*[1 1 1]);
            make_title_type_time(type, time, par);
            caxis([-100 100]);
            cb = colorbar('limits', [-100 100], 'ytick', [-100:20:100], 'location', 'eastoutside');
            cb.TickLabelInterpreter = 'latex'; cb.Label.Interpreter = 'latex';
            ylabel(cb, '$(\Gamma_m - \Gamma)/\Gamma_m$ (\%)');
            xlabel('Latitude (deg)');
            ylabel('$\sigma$ (unitless)');
            set(gca, 'ydir', 'reverse', 'ylim', [0.3 1], 'ytick', [0.3:0.1:1], 'xlim', [-90 90], 'xtick', [-90:30:90], 'yminortick', 'on', 'tickdir', 'out', 'xminortick', 'on');
            folder = sprintf('%s/ga_frac/%s/%s', plotdir, land, time);
            if ~exist(folder, 'dir')
                mkdir(folder);
            end
            print(sprintf('%s/ga_frac_lat_lev', folder), '-dpng', '-r300');
            close;

            if strcmp(type, 'rea') | (strcmp(type, 'gcm') & strcmp(par.model, 'mmm'))
                % locate rce, rcae, and rae for NH and SH
                idx_rce = find(r1_ann<=par.ep);
                idx_rcae = find(r1_ann>par.ep);
                idx_rae = find(r1_ann>=par.ga);

                idx_rce_nh = find(r1_ann<=par.ep & grid.dim3.lat>0);
                idx_rcae_nh = find(r1_ann>par.ep & r1_ann<par.ga & grid.dim3.lat>0);
                idx_rae_nh = find(r1_ann>=par.ga & grid.dim3.lat>0);

                idx_rce_sh = find(r1_ann<=par.ep & grid.dim3.lat<0);
                idx_rcae_sh = find(r1_ann>par.ep & r1_ann<par.ga & grid.dim3.lat<0);
                idx_rae_sh = find(r1_ann>=par.ga & grid.dim3.lat<0);

                % take area averaged temperature profile
                ga_fr_rce = nansum( ga_fr_ann(idx_rce, :) .* repmat(cosd(grid.dim3.lat(idx_rce)), [1 size(ga_fr_ann,2)]) )/nansum(cosd(grid.dim3.lat(idx_rce)));
                ga_fr_rcae = nansum( ga_fr_ann(idx_rcae, :) .* repmat(cosd(grid.dim3.lat(idx_rcae)), [1 size(ga_fr_ann,2)]) )/nansum(cosd(grid.dim3.lat(idx_rcae)));
                ga_fr_rae = nansum( ga_fr_ann(idx_rae, :) .* repmat(cosd(grid.dim3.lat(idx_rae)), [1 size(ga_fr_ann,2)]) )/nansum(cosd(grid.dim3.lat(idx_rae)));

                ga_fr_std_rce = nansum( ga_fr_std_ann(idx_rce, :) .* repmat(cosd(grid.dim3.lat(idx_rce)), [1 size(ga_fr_std_ann,2)]) )/nansum(cosd(grid.dim3.lat(idx_rce)));
                ga_fr_std_rcae = nansum( ga_fr_std_ann(idx_rcae, :) .* repmat(cosd(grid.dim3.lat(idx_rcae)), [1 size(ga_fr_std_ann,2)]) )/nansum(cosd(grid.dim3.lat(idx_rcae)));
                ga_fr_std_rae = nansum( ga_fr_std_ann(idx_rae, :) .* repmat(cosd(grid.dim3.lat(idx_rae)), [1 size(ga_fr_std_ann,2)]) )/nansum(cosd(grid.dim3.lat(idx_rae)));

                ga_fr_min_rce = nansum( ga_fr_min_ann(idx_rce, :) .* repmat(cosd(grid.dim3.lat(idx_rce)), [1 size(ga_fr_min_ann,2)]) )/nansum(cosd(grid.dim3.lat(idx_rce)));
                ga_fr_min_rcae = nansum( ga_fr_min_ann(idx_rcae, :) .* repmat(cosd(grid.dim3.lat(idx_rcae)), [1 size(ga_fr_min_ann,2)]) )/nansum(cosd(grid.dim3.lat(idx_rcae)));
                ga_fr_min_rae = nansum( ga_fr_min_ann(idx_rae, :) .* repmat(cosd(grid.dim3.lat(idx_rae)), [1 size(ga_fr_min_ann,2)]) )/nansum(cosd(grid.dim3.lat(idx_rae)));

                ga_fr_max_rce = nansum( ga_fr_max_ann(idx_rce, :) .* repmat(cosd(grid.dim3.lat(idx_rce)), [1 size(ga_fr_max_ann,2)]) )/nansum(cosd(grid.dim3.lat(idx_rce)));
                ga_fr_max_rcae = nansum( ga_fr_max_ann(idx_rcae, :) .* repmat(cosd(grid.dim3.lat(idx_rcae)), [1 size(ga_fr_max_ann,2)]) )/nansum(cosd(grid.dim3.lat(idx_rcae)));
                ga_fr_max_rae = nansum( ga_fr_max_ann(idx_rae, :) .* repmat(cosd(grid.dim3.lat(idx_rae)), [1 size(ga_fr_max_ann,2)]) )/nansum(cosd(grid.dim3.lat(idx_rae)));

                ga_fr_25_rce = nansum( ga_fr_25_ann(idx_rce, :) .* repmat(cosd(grid.dim3.lat(idx_rce)), [1 size(ga_fr_25_ann,2)]) )/nansum(cosd(grid.dim3.lat(idx_rce)));
                ga_fr_25_rcae = nansum( ga_fr_25_ann(idx_rcae, :) .* repmat(cosd(grid.dim3.lat(idx_rcae)), [1 size(ga_fr_25_ann,2)]) )/nansum(cosd(grid.dim3.lat(idx_rcae)));
                ga_fr_25_rae = nansum( ga_fr_25_ann(idx_rae, :) .* repmat(cosd(grid.dim3.lat(idx_rae)), [1 size(ga_fr_25_ann,2)]) )/nansum(cosd(grid.dim3.lat(idx_rae)));

                ga_fr_75_rce = nansum( ga_fr_75_ann(idx_rce, :) .* repmat(cosd(grid.dim3.lat(idx_rce)), [1 size(ga_fr_75_ann,2)]) )/nansum(cosd(grid.dim3.lat(idx_rce)));
                ga_fr_75_rcae = nansum( ga_fr_75_ann(idx_rcae, :) .* repmat(cosd(grid.dim3.lat(idx_rcae)), [1 size(ga_fr_75_ann,2)]) )/nansum(cosd(grid.dim3.lat(idx_rcae)));
                ga_fr_75_rae = nansum( ga_fr_75_ann(idx_rae, :) .* repmat(cosd(grid.dim3.lat(idx_rae)), [1 size(ga_fr_75_ann,2)]) )/nansum(cosd(grid.dim3.lat(idx_rae)));

                ga_fr_rce_nh = nansum( ga_fr_ann(idx_rce_nh, :) .* repmat(cosd(grid.dim3.lat(idx_rce_nh)), [1 size(ga_fr_ann,2)]) )/nansum(cosd(grid.dim3.lat(idx_rce_nh)));
                ga_fr_rcae_nh = nansum( ga_fr_ann(idx_rcae_nh, :) .* repmat(cosd(grid.dim3.lat(idx_rcae_nh)), [1 size(ga_fr_ann,2)]) )/nansum(cosd(grid.dim3.lat(idx_rcae_nh)));
                ga_fr_rae_nh = nansum( ga_fr_ann(idx_rae_nh, :) .* repmat(cosd(grid.dim3.lat(idx_rae_nh)), [1 size(ga_fr_ann,2)]) )/nansum(cosd(grid.dim3.lat(idx_rae_nh)));

                ga_fr_std_rce_nh = nansum( ga_fr_std_ann(idx_rce_nh, :) .* repmat(cosd(grid.dim3.lat(idx_rce_nh)), [1 size(ga_fr_std_ann,2)]) )/nansum(cosd(grid.dim3.lat(idx_rce_nh)));
                ga_fr_std_rcae_nh = nansum( ga_fr_std_ann(idx_rcae_nh, :) .* repmat(cosd(grid.dim3.lat(idx_rcae_nh)), [1 size(ga_fr_std_ann,2)]) )/nansum(cosd(grid.dim3.lat(idx_rcae_nh)));
                ga_fr_std_rae_nh = nansum( ga_fr_std_ann(idx_rae_nh, :) .* repmat(cosd(grid.dim3.lat(idx_rae_nh)), [1 size(ga_fr_std_ann,2)]) )/nansum(cosd(grid.dim3.lat(idx_rae_nh)));

                ga_fr_min_rce_nh = nansum( ga_fr_min_ann(idx_rce_nh, :) .* repmat(cosd(grid.dim3.lat(idx_rce_nh)), [1 size(ga_fr_min_ann,2)]) )/nansum(cosd(grid.dim3.lat(idx_rce_nh)));
                ga_fr_min_rcae_nh = nansum( ga_fr_min_ann(idx_rcae_nh, :) .* repmat(cosd(grid.dim3.lat(idx_rcae_nh)), [1 size(ga_fr_min_ann,2)]) )/nansum(cosd(grid.dim3.lat(idx_rcae_nh)));
                ga_fr_min_rae_nh = nansum( ga_fr_min_ann(idx_rae_nh, :) .* repmat(cosd(grid.dim3.lat(idx_rae_nh)), [1 size(ga_fr_min_ann,2)]) )/nansum(cosd(grid.dim3.lat(idx_rae_nh)));

                ga_fr_max_rce_nh = nansum( ga_fr_max_ann(idx_rce_nh, :) .* repmat(cosd(grid.dim3.lat(idx_rce_nh)), [1 size(ga_fr_max_ann,2)]) )/nansum(cosd(grid.dim3.lat(idx_rce_nh)));
                ga_fr_max_rcae_nh = nansum( ga_fr_max_ann(idx_rcae_nh, :) .* repmat(cosd(grid.dim3.lat(idx_rcae_nh)), [1 size(ga_fr_max_ann,2)]) )/nansum(cosd(grid.dim3.lat(idx_rcae_nh)));
                ga_fr_max_rae_nh = nansum( ga_fr_max_ann(idx_rae_nh, :) .* repmat(cosd(grid.dim3.lat(idx_rae_nh)), [1 size(ga_fr_max_ann,2)]) )/nansum(cosd(grid.dim3.lat(idx_rae_nh)));

                ga_fr_25_rce_nh = nansum( ga_fr_25_ann(idx_rce_nh, :) .* repmat(cosd(grid.dim3.lat(idx_rce_nh)), [1 size(ga_fr_25_ann,2)]) )/nansum(cosd(grid.dim3.lat(idx_rce_nh)));
                ga_fr_25_rcae_nh = nansum( ga_fr_25_ann(idx_rcae_nh, :) .* repmat(cosd(grid.dim3.lat(idx_rcae_nh)), [1 size(ga_fr_25_ann,2)]) )/nansum(cosd(grid.dim3.lat(idx_rcae_nh)));
                ga_fr_25_rae_nh = nansum( ga_fr_25_ann(idx_rae_nh, :) .* repmat(cosd(grid.dim3.lat(idx_rae_nh)), [1 size(ga_fr_25_ann,2)]) )/nansum(cosd(grid.dim3.lat(idx_rae_nh)));

                ga_fr_75_rce_nh = nansum( ga_fr_75_ann(idx_rce_nh, :) .* repmat(cosd(grid.dim3.lat(idx_rce_nh)), [1 size(ga_fr_75_ann,2)]) )/nansum(cosd(grid.dim3.lat(idx_rce_nh)));
                ga_fr_75_rcae_nh = nansum( ga_fr_75_ann(idx_rcae_nh, :) .* repmat(cosd(grid.dim3.lat(idx_rcae_nh)), [1 size(ga_fr_75_ann,2)]) )/nansum(cosd(grid.dim3.lat(idx_rcae_nh)));
                ga_fr_75_rae_nh = nansum( ga_fr_75_ann(idx_rae_nh, :) .* repmat(cosd(grid.dim3.lat(idx_rae_nh)), [1 size(ga_fr_75_ann,2)]) )/nansum(cosd(grid.dim3.lat(idx_rae_nh)));

                ga_fr_rce_sh = nansum( ga_fr_ann(idx_rce_sh, :) .* repmat(cosd(grid.dim3.lat(idx_rce_sh)), [1 size(ga_fr_ann,2)]) )/nansum(cosd(grid.dim3.lat(idx_rce_sh)));
                ga_fr_rcae_sh = nansum( ga_fr_ann(idx_rcae_sh, :) .* repmat(cosd(grid.dim3.lat(idx_rcae_sh)), [1 size(ga_fr_ann,2)]) )/nansum(cosd(grid.dim3.lat(idx_rcae_sh)));
                ga_fr_rae_sh = nansum( ga_fr_ann(idx_rae_sh, :) .* repmat(cosd(grid.dim3.lat(idx_rae_sh)), [1 size(ga_fr_ann,2)]) )/nansum(cosd(grid.dim3.lat(idx_rae_sh)));

                ga_fr_std_rce_sh = nansum( ga_fr_std_ann(idx_rce_sh, :) .* repmat(cosd(grid.dim3.lat(idx_rce_sh)), [1 size(ga_fr_std_ann,2)]) )/nansum(cosd(grid.dim3.lat(idx_rce_sh)));
                ga_fr_std_rcae_sh = nansum( ga_fr_std_ann(idx_rcae_sh, :) .* repmat(cosd(grid.dim3.lat(idx_rcae_sh)), [1 size(ga_fr_std_ann,2)]) )/nansum(cosd(grid.dim3.lat(idx_rcae_sh)));
                ga_fr_std_rae_sh = nansum( ga_fr_std_ann(idx_rae_sh, :) .* repmat(cosd(grid.dim3.lat(idx_rae_sh)), [1 size(ga_fr_std_ann,2)]) )/nansum(cosd(grid.dim3.lat(idx_rae_sh)));

                ga_fr_min_rce_sh = nansum( ga_fr_min_ann(idx_rce_sh, :) .* repmat(cosd(grid.dim3.lat(idx_rce_sh)), [1 size(ga_fr_min_ann,2)]) )/nansum(cosd(grid.dim3.lat(idx_rce_sh)));
                ga_fr_min_rcae_sh = nansum( ga_fr_min_ann(idx_rcae_sh, :) .* repmat(cosd(grid.dim3.lat(idx_rcae_sh)), [1 size(ga_fr_min_ann,2)]) )/nansum(cosd(grid.dim3.lat(idx_rcae_sh)));
                ga_fr_min_rae_sh = nansum( ga_fr_min_ann(idx_rae_sh, :) .* repmat(cosd(grid.dim3.lat(idx_rae_sh)), [1 size(ga_fr_min_ann,2)]) )/nansum(cosd(grid.dim3.lat(idx_rae_sh)));

                ga_fr_max_rce_sh = nansum( ga_fr_max_ann(idx_rce_sh, :) .* repmat(cosd(grid.dim3.lat(idx_rce_sh)), [1 size(ga_fr_max_ann,2)]) )/nansum(cosd(grid.dim3.lat(idx_rce_sh)));
                ga_fr_max_rcae_sh = nansum( ga_fr_max_ann(idx_rcae_sh, :) .* repmat(cosd(grid.dim3.lat(idx_rcae_sh)), [1 size(ga_fr_max_ann,2)]) )/nansum(cosd(grid.dim3.lat(idx_rcae_sh)));
                ga_fr_max_rae_sh = nansum( ga_fr_max_ann(idx_rae_sh, :) .* repmat(cosd(grid.dim3.lat(idx_rae_sh)), [1 size(ga_fr_max_ann,2)]) )/nansum(cosd(grid.dim3.lat(idx_rae_sh)));

                ga_fr_25_rce_sh = nansum( ga_fr_25_ann(idx_rce_sh, :) .* repmat(cosd(grid.dim3.lat(idx_rce_sh)), [1 size(ga_fr_25_ann,2)]) )/nansum(cosd(grid.dim3.lat(idx_rce_sh)));
                ga_fr_25_rcae_sh = nansum( ga_fr_25_ann(idx_rcae_sh, :) .* repmat(cosd(grid.dim3.lat(idx_rcae_sh)), [1 size(ga_fr_25_ann,2)]) )/nansum(cosd(grid.dim3.lat(idx_rcae_sh)));
                ga_fr_25_rae_sh = nansum( ga_fr_25_ann(idx_rae_sh, :) .* repmat(cosd(grid.dim3.lat(idx_rae_sh)), [1 size(ga_fr_25_ann,2)]) )/nansum(cosd(grid.dim3.lat(idx_rae_sh)));

                ga_fr_75_rce_sh = nansum( ga_fr_75_ann(idx_rce_sh, :) .* repmat(cosd(grid.dim3.lat(idx_rce_sh)), [1 size(ga_fr_75_ann,2)]) )/nansum(cosd(grid.dim3.lat(idx_rce_sh)));
                ga_fr_75_rcae_sh = nansum( ga_fr_75_ann(idx_rcae_sh, :) .* repmat(cosd(grid.dim3.lat(idx_rcae_sh)), [1 size(ga_fr_75_ann,2)]) )/nansum(cosd(grid.dim3.lat(idx_rcae_sh)));
                ga_fr_75_rae_sh = nansum( ga_fr_75_ann(idx_rae_sh, :) .* repmat(cosd(grid.dim3.lat(idx_rae_sh)), [1 size(ga_fr_75_ann,2)]) )/nansum(cosd(grid.dim3.lat(idx_rae_sh)));

                % ALL NH compared with moist adiabat
                figure(); clf; hold all; box on;
                line([0 0], [0 1], 'color', 'k', 'linewidth', 0.5);
                line(100*[1 1], [0 1], 'linestyle', '--', 'color', 'k', 'linewidth', 0.5);
                if strcmp(type, 'rea')
                    ga_fr_rce_nh_u = ga_fr_max_rce_nh;
                    ga_fr_rce_nh_l = ga_fr_min_rce_nh;
                    ga_fr_rcae_nh_u = ga_fr_max_rcae_nh;
                    ga_fr_rcae_nh_l = ga_fr_min_rcae_nh;
                    ga_fr_rae_nh_u = ga_fr_max_rae_nh;
                    ga_fr_rae_nh_l = ga_fr_min_rae_nh;
                elseif (strcmp(type, 'gcm') & strcmp(par.model, 'mmm'))
                    ga_fr_rce_nh_u = ga_fr_75_rce_nh;
                    ga_fr_rce_nh_l = ga_fr_25_rce_nh;
                    ga_fr_rcae_nh_u = ga_fr_75_rcae_nh;
                    ga_fr_rcae_nh_l = ga_fr_25_rcae_nh;
                    ga_fr_rae_nh_u = ga_fr_75_rae_nh;
                    ga_fr_rae_nh_l = ga_fr_25_rae_nh;
                end

                if strcmp(type, 'rea') | (strcmp(type, 'gcm') & strcmp(par.model, 'mmm'))
                    ga_fr_rce_nh_2 = [ga_fr_rce_nh_l'; flipud(ga_fr_rce_nh_u')];
                    ga_fr_rcae_nh_2 = [ga_fr_rcae_nh_l'; flipud(ga_fr_rcae_nh_u')];
                    ga_fr_rae_nh_2 = [ga_fr_rae_nh_l'; flipud(ga_fr_rae_nh_u')];
                    si2 = [grid.dim3.si'; flipud(grid.dim3.si')];

                    fill(ga_fr_rae_nh_2, si2, par.blue, 'facealpha', 0.2, 'edgealpha', 0.2);
                    fill(ga_fr_rcae_nh_2, si2, 0.25*[1 1 1], 'facealpha', 0.2, 'edgealpha', 0.2);
                    fill(ga_fr_rce_nh_2, si2, par.orange, 'facealpha', 0.2, 'edgealpha', 0.2);
                end
                h_rce = plot(ga_fr_rce_nh, grid.dim3.si, 'color', par.orange);
                h_rcae = plot(ga_fr_rcae_nh, grid.dim3.si, 'color', 0.25*[1 1 1]);
                h_rae = plot(ga_fr_rae_nh, grid.dim3.si, 'color', par.blue);
                text(ga_fr_rce_nh(60)-30, grid.dim3.si(60), '\textbf{RCE}', 'color', par.orange, 'fontsize', 6);
                text(ga_fr_rcae_nh(45), grid.dim3.si(45), '\textbf{RCAE}', 'color', 0.25*[1 1 1], 'fontsize', 6);
                text(ga_fr_rae_nh(30)+5, grid.dim3.si(30), '\textbf{RAE}', 'color', par.blue, 'fontsize', 6);
                % draw arrow and label as inversion 
                arrows(105, 0.85, 40, 0, 'Cartesian', [6e-4, 0.1, 0.05, 1e-5]);
                text(105, 0.8, 'Inversion', 'fontsize', 7);
                xlabel('$\frac{\Gamma_m-\Gamma}{\Gamma_m}$ (\%)'); ylabel('$\sigma$ (unitless)');
                title(sprintf('NH, %s', upper(time)));
                % legend([h_rce, h_rcae, h_rae], 'RCE', 'RCAE', 'RAE', 'location', 'southwest');
                axis('tight');
                set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_sq)
                set(gca, 'fontsize', par.fs, 'xlim', [-50 150], 'ydir', 'reverse', 'yscale', 'linear', 'ytick', [0:0.1:1], 'ylim', [0.3 1], 'xminortick', 'on')
                folder = sprintf('%s/eps_%g_ga_%g/%s/%s/%s/%s/ga_fr/', plotdir, par.ep, par.ga, fw, crit, land, time);
                if ~exist(folder, 'dir')
                    mkdir(folder);
                end
                print(sprintf('%sall_nh_ann', folder), '-dpng', '-r300');
                if par.make_tikz
                    matlab2tikz(sprintf('%sall_nh_ann.tex', folder));
                end
                close;

                % vertically averaged deviation
                vavg_dev_rce = nanmean(interp1(grid.dim3.si, ga_fr_rce, linspace(0.9, 0.3, 101)));
                disp(sprintf('The vertically integrated lapse rate deviation from MALR over all RCE is %g%%.', vavg_dev_rce))
                vavg_dev_rcae = nanmean(interp1(grid.dim3.si, ga_fr_rcae, linspace(0.9, 0.3, 101)));
                disp(sprintf('The vertically integrated lapse rate deviation from MALR over all RCAE is %g%%.', vavg_dev_rcae))

                % vertically averaged deviation
                vavg_dev_rce = nanmean(interp1(grid.dim3.si, ga_fr_rce_nh, linspace(0.9, 0.3, 101)));
                disp(sprintf('The vertically integrated lapse rate deviation from MALR over NH RCE is %g%%.', vavg_dev_rce))
                vavg_dev_rcae = nanmean(interp1(grid.dim3.si, ga_fr_rcae_nh, linspace(0.9, 0.3, 101)));
                disp(sprintf('The vertically integrated lapse rate deviation from MALR over NH RCAE is %g%%.', vavg_dev_rcae))

                % ALL SH compared with moist adiabat
                figure(); clf; hold all; box on;
                line([0 0], [0 1], 'color', 'k', 'linewidth', 0.5);
                line(100*[1 1], [0 1], 'linestyle', '--', 'color', 'k', 'linewidth', 0.5);
                if strcmp(type, 'rea')
                    ga_fr_rce_sh_u = ga_fr_max_rce_sh;
                    ga_fr_rce_sh_l = ga_fr_min_rce_sh;
                    ga_fr_rcae_sh_u = ga_fr_max_rcae_sh;
                    ga_fr_rcae_sh_l = ga_fr_min_rcae_sh;
                    ga_fr_rae_sh_u = ga_fr_max_rae_sh;
                    ga_fr_rae_sh_l = ga_fr_min_rae_sh;
                elseif (strcmp(type, 'gcm') & strcmp(par.model, 'mmm'))
                    ga_fr_rce_sh_u = ga_fr_75_rce_sh;
                    ga_fr_rce_sh_l = ga_fr_25_rce_sh;
                    ga_fr_rcae_sh_u = ga_fr_75_rcae_sh;
                    ga_fr_rcae_sh_l = ga_fr_25_rcae_sh;
                    ga_fr_rae_sh_u = ga_fr_75_rae_sh;
                    ga_fr_rae_sh_l = ga_fr_25_rae_sh;
                end

                if strcmp(type, 'rea') | (strcmp(type, 'gcm') & strcmp(par.model, 'mmm'))
                    ga_fr_rce_sh_2 = [ga_fr_rce_sh_l'; flipud(ga_fr_rce_sh_u')];
                    ga_fr_rcae_sh_2 = [ga_fr_rcae_sh_l'; flipud(ga_fr_rcae_sh_u')];
                    ga_fr_rae_sh_2 = [ga_fr_rae_sh_l'; flipud(ga_fr_rae_sh_u')];
                    si2 = [grid.dim3.si'; flipud(grid.dim3.si')];

                    fill(ga_fr_rae_sh_2, si2, par.blue, 'facealpha', 0.2, 'edgealpha', 0.2);
                    fill(ga_fr_rcae_sh_2, si2, 0.25*[1 1 1], 'facealpha', 0.2, 'edgealpha', 0.2);
                    fill(ga_fr_rce_sh_2, si2, par.orange, 'facealpha', 0.2, 'edgealpha', 0.2);
                end
                h_rce = plot(ga_fr_rce_sh, grid.dim3.si, 'color', par.orange);
                h_rcae = plot(ga_fr_rcae_sh, grid.dim3.si, 'color', 0.25*[1 1 1]);
                h_rae = plot(ga_fr_rae_sh, grid.dim3.si, 'color', par.blue);
                text(ga_fr_rce_sh(60)-30, grid.dim3.si(60), '\textbf{RCE}', 'color', par.orange, 'fontsize', 6);
                text(ga_fr_rcae_sh(45), grid.dim3.si(45), '\textbf{RCAE}', 'color', 0.25*[1 1 1], 'fontsize', 6);
                text(ga_fr_rae_sh(30)+5, grid.dim3.si(30), '\textbf{RAE}', 'color', par.blue, 'fontsize', 6);
                % draw arrow and label as inversion 
                arrows(105, 0.85, 40, 0, 'Cartesian', [6e-4, 0.1, 0.05, 1e-5]);
                text(105, 0.8, 'Inversion', 'fontsize', 7);
                xlabel('$\frac{\Gamma_m-\Gamma}{\Gamma_m}$ (\%)'); ylabel('$\sigma$ (unitless)');
                title(sprintf('SH, %s', upper(time)));
                % legend([h_rce, h_rcae, h_rae], 'RCE', 'RCAE', 'RAE', 'location', 'southwest');
                axis('tight');
                set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_sq)
                set(gca, 'fontsize', par.fs, 'xlim', [-50 150], 'ydir', 'reverse', 'yscale', 'linear', 'ytick', [0:0.1:1], 'ylim', [0.3 1], 'xminortick', 'on')
                print(sprintf('%sall_sh_ann', folder), '-dpng', '-r300');
                if par.make_tikz
                    matlab2tikz(sprintf('%sall_sh_ann.tex', folder));
                end
                close;

                % vertically averaged deviation
                vavg_dev_rce = nanmean(interp1(grid.dim3.si, ga_fr_rce_sh, linspace(0.9, 0.3, 101)));
                disp(sprintf('The vertically integrated lapse rate deviation from MALR over SH RCE is %g%%.', vavg_dev_rce))
                vavg_dev_rcae = nanmean(interp1(grid.dim3.si, ga_fr_rcae_sh, linspace(0.9, 0.3, 101)));
                disp(sprintf('The vertically integrated lapse rate deviation from MALR over SH RCAE is %g%%.', vavg_dev_rcae))
            end

        end % land
    end % MSE/DSE framework
    % Legend for RCE and RAE separated into NH and SH
end
